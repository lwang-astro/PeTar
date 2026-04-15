#pragma once
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "Common/Float.h"
#include "Common/binary_tree.h"
#include "changeover.hpp"
#include "AR/force.h"
#include "hard_ptcl.hpp"
#include "Hermite/hermite_particle.h"
#include "ar_perturber.hpp"
#ifdef BSE_BASE
#include "bse_interface.h"
#include "two_body_tide.hpp"
#endif
#ifdef DISK_STAR_MERGER
#include "disk_star_merger.hpp"
#endif
#include "external_hard.hpp"

#ifdef SDAR_PN
#include "pn.hpp"
#endif
//#define GR_PRECESSION

//! AR interaction clas
class ARInteraction{
public:
    typedef H4::ParticleH4<PtclHard> H4Ptcl;
    Float eps_sq; ///> softening parameter
    Float gravitational_constant;
    int interrupt_detection_option;    // 0: no interruption; 1: merge when the pair distance is less than the sum of two members' radii; 2: record binary status instead of merger
#ifdef STELLAR_EVOLUTION
    Float time_interrupt_max;
#ifdef BSE_BASE
    int stellar_evolution_option;
    bool stellar_evolution_write_flag;
    BSEManager bse_manager;
    TwoBodyTide tide;
    std::ofstream fout_sse; ///> log file for SSE event
    std::ofstream fout_bse; ///> log file for BSE event
#else
#ifdef DISK_STAR_MERGER
    DiskStarMergerManager disk_star_merger_manager;
#endif
    std::ofstream fout_interrupt; ///> log file for interrupted binary
#endif
#endif
#ifdef EXTERNAL_HARD
    ExternalHardForce *ext_force; // external hard to calculate perturbation
#endif
#ifdef SDAR_PN
    PostNewtonian pn; // PN force for AR (compile-time enabled by SDAR_PN)
#endif



    ARInteraction(): eps_sq(Float(-1.0)), gravitational_constant(Float(-1.0)), interrupt_detection_option(0)
#ifdef STELLAR_EVOLUTION
                   , time_interrupt_max(NUMERIC_FLOAT_MAX) 
#ifdef BSE_BASE
                   , stellar_evolution_option(0), stellar_evolution_write_flag(false), bse_manager(), tide(), fout_sse(), fout_bse()
#else
#ifdef DISK_STAR_MERGER
                   , disk_star_merger_manager()
#endif
                   , fout_interrupt()  
#endif
#endif
#ifdef EXTERNAL_HARD
                   , ext_force(NULL)
#endif
                   
    {}

    //! (Necessary) check whether publicly initialized parameters are correctly set
    /*! \return true: all parmeters are correct. In this case no parameters, return true;
     */
    bool checkParams() {
        ASSERT(eps_sq>=0.0);
        ASSERT(gravitational_constant>0.0);
        ASSERT(interrupt_detection_option>=0 && interrupt_detection_option<=2);
#ifdef SDAR_PN
        ASSERT(pn.checkParams());
#endif
#ifdef STELLAR_EVOLUTION
        ASSERT(time_interrupt_max>=0.0);
#ifdef BSE_BASE
        ASSERT(stellar_evolution_option==0 || (stellar_evolution_option==1 && bse_manager.checkParams()) || (stellar_evolution_option==2 && bse_manager.checkParams() && tide.checkParams()));
        ASSERT(!stellar_evolution_write_flag||(stellar_evolution_write_flag&&fout_sse.is_open()));
        ASSERT(!stellar_evolution_write_flag||(stellar_evolution_write_flag&&fout_bse.is_open()));
#else
#ifdef DISK_STAR_MERGER
        ASSERT(disk_star_merger_manager.checkParams());
#endif
        ASSERT(interrupt_detection_option==0||(interrupt_detection_option>0&&fout_interrupt.is_open()));
#endif
#endif
        return true;
    }

    //! print parameters
    void print(std::ostream & _fout) const{
        _fout<<"eps_sq : "<<eps_sq<<std::endl
             <<"G      : "<<gravitational_constant<<std::endl
             <<"Interrupt_opt: "<<interrupt_detection_option<<std::endl;
#ifdef SDAR_PN
        pn.print(_fout);
#endif
#ifdef STELLAR_EVOLUTION
#ifdef BSE_BASE
        _fout<<"SE_opt : "<<stellar_evolution_option<<std::endl;
#endif
#ifdef DISK_STAR_MERGER
        disk_star_merger_manager.print(_fout);
#endif
#endif
    }    

    //! (Necessary) calculate inner member acceleration, potential and inverse time transformation function gradient and factor for kick (two-body case)
    /*!
      @param[out] _f1: force for particle 1 to store the calculation results (in acc_in[3] for acceleration and gtgrad[3] for gradient, notice acc/gtgard are overwritten, not accummulating old values)
      @param[out] _f2: force for particle 2
      @param[out] _epot: total inner potential energy
      @param[in] _p1: particle 1
      @param[in] _p2: particle 2
      @param[in] _pos_offset: position offset need to be added to calculate dr
      \return the inverse time transformation factor (gt_kick_inv) for kick step
    */
    inline Float calcInnerAccPotAndGTKickInvTwo(AR::Force& _f1, AR::Force& _f2, Float& _epot, const PtclHard& _p1, const PtclHard& _p2, const Float* _pos_offset) {
        // acceleration
        const Float mass1 = _p1.mass;
        const auto& pos1 = _p1.pos;

        const Float mass2 = _p2.mass;
        const auto& pos2 = _p2.pos;

        Float gm1 = gravitational_constant*mass1;
        Float gm2 = gravitational_constant*mass2;
        Float gm1m2 = gm1*mass2;

#ifdef USE_CM_FRAME
        Float dr[3] = {pos2[0] -pos1[0] + _pos_offset[0],
                       pos2[1] -pos1[1] + _pos_offset[1],
                       pos2[2] -pos1[2] + _pos_offset[2]};
#else
        Float dr[3] = {pos2[0] -pos1[0],
                       pos2[1] -pos1[1],
                       pos2[2] -pos1[2]};
#endif
        Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
        Float inv_r = 1.0/sqrt(r2);
        Float inv_r3 = inv_r*inv_r*inv_r;

        Float* acc1 = _f1.acc_in;
        Float* acc2 = _f2.acc_in;

#ifdef AR_CHANGEOVER
        auto& ch1 = _p1.changeover;
        auto& ch2 = _p2.changeover;

        Float r = r2*invr;
        Float k = ChangeOver::calcAcc0WTwo(ch1,ch2,r);
        Float kpot = ChangeOver::calcPotWTwo(ch1,ch2,r);

        Float gmor3_1 = gm2*inv_r3*k;
        Float gmor3_2 = gm1*inv_r3*k;

        Float inv_rk = inv_r*kpot;
        Float gm1or =  gm1*inv_rk;
        Float gm2or =  gm2*inv_rk;
        Float gm1m2or = gm1m2*inv_rk;
#else
        Float gmor3_1 = gm2*inv_r3;
        Float gmor3_2 = gm1*inv_r3;

        Float gm1or =  gm1*inv_r;
        Float gm2or =  gm2*inv_r;
        Float gm1m2or = gm1m2*inv_r;
#endif

        acc1[0] = gmor3_1 * dr[0];
        acc1[1] = gmor3_1 * dr[1];
        acc1[2] = gmor3_1 * dr[2];

        _f1.pot_in = -gm2or;

        acc2[0] = - gmor3_2 * dr[0];
        acc2[1] = - gmor3_2 * dr[1];
        acc2[2] = - gmor3_2 * dr[2];

        _f2.pot_in = -gm1or;


#ifdef AR_TTL
        // trans formation function gradient
#ifdef AR_CHANGEOVER
        Float gm1m2or3 = gm1m2*inv_r3*k;
#else
        Float gm1m2or3 = gm1m2*inv_r3;
#endif
        Float* gtgrad1 = _f1.gtgrad;
        Float* gtgrad2 = _f2.gtgrad;
        gtgrad1[0] = gm1m2or3 * dr[0];
        gtgrad1[1] = gm1m2or3 * dr[1];
        gtgrad1[2] = gm1m2or3 * dr[2];

        gtgrad2[0] = - gtgrad1[0];
        gtgrad2[1] = - gtgrad1[1];
        gtgrad2[2] = - gtgrad1[2];
#endif

        // potential energy
        _epot = - gm1m2or;

        // transformation factor for kick
        Float gt_kick_inv = gm1m2or;

        return gt_kick_inv;
    }

    //! calculate inner member acceleration, potential and inverse time transformation function gradient and factor for kick
    /*!
      @param[out] _force: force array to store the calculation results (in acc_in[3] for acceleration and gtgrad[3] for gradient, notice acc/gtgard may need to reset zero to avoid accummulating old values)
      @param[out] _epot: total inner potential energy
      @param[in] _particles: member particle array
      @param[in] _n_particle: number of member particles
      \return the inverse time transformation factor (gt_kick_inv) for kick step
    */
    inline Float calcInnerAccPotAndGTKickInv(AR::Force* _force, Float& _epot, const PtclHard* _particles, const int _n_particle) {
        _epot = Float(0.0);
        Float gt_kick_inv = Float(0.0);

        for (int i=0; i<_n_particle; i++) {
            const Float massi = _particles[i].mass;
            const Float* posi = &_particles[i].pos.x;
            Float* acci = _force[i].acc_in;
            acci[0] = acci[1] = acci[2] = Float(0.0);

#ifdef AR_TTL
            Float* gtgradi = _force[i].gtgrad;
            gtgradi[0] = gtgradi[1] = gtgradi[2] = Float(0.0);
#endif

            Float poti = Float(0.0);
            Float gtki = Float(0.0);

            for (int j=0; j<_n_particle; j++) {
                if (i==j) continue;
                const Float massj = _particles[j].mass;
                const Float* posj = &_particles[j].pos.x;
                Float dr[3] = {posj[0] -posi[0],
                               posj[1] -posi[1],
                               posj[2] -posi[2]};
                Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
                Float inv_r = 1.0/sqrt(r2);
                Float inv_r3 = inv_r*inv_r*inv_r;
#ifdef AR_CHANGEOVER
                Float r = r2*inv_r;
                const Float kpot  = ChangeOver::calcPotWTwo(_particles[i].changeover, _particles[j].changeover, r);
                const Float k     = ChangeOver::calcAcc0WTwo(_particles[i].changeover, _particles[j].changeover, r);

                Float gmor3 = gravitational_constant*massj*inv_r3*k;
                Float gmor = gravitational_constant*massj*inv_r*kpot;
#else
                Float gmor3 = gravitational_constant*massj*inv_r3;
                Float gmor = gravitational_constant*massj*inv_r;
#endif
                acci[0] += gmor3 * dr[0];
                acci[1] += gmor3 * dr[1];
                acci[2] += gmor3 * dr[2];

#ifdef AR_TTL
                Float gmimjor3 = massi*gmor3;
                gtgradi[0] += gmimjor3 * dr[0];
                gtgradi[1] += gmimjor3 * dr[1];
                gtgradi[2] += gmimjor3 * dr[2];
#endif

                poti -= gmor;
                gtki += gmor;

            }
            _epot += poti * massi;
            gt_kick_inv += gtki * massi;
        }
        _epot   *= 0.5;
        gt_kick_inv *= 0.5;

        return gt_kick_inv;
    }

    //! (Necessary) calculate acceleration from perturber and the perturbation factor for slowdown calculation
    /*!@param[out] _force: force array to store the calculation results (in acc_pert[3], notice acc_pert may need to reset zero to avoid accummulating old values)
      @param[in] _particles: member particle array
      @param[in] _n_particle: number of member particles
      @param[in] _particle_cm: center-of-mass particle
      @param[in] _perturber: pertuber container
      @param[in] _time: current time
    */
    void calcAccPert(AR::Force* _force, const PtclHard* _particles, const int _n_particle, const H4Ptcl& _particle_cm, const ARPerturber& _perturber, const Float _time) {
        static const Float inv3 = 1.0 / 3.0;

        // perturber force
        const int n_pert = _perturber.neighbor_address.getSize();
        const int n_pert_single = _perturber.n_neighbor_single;
        const int n_pert_group = _perturber.n_neighbor_group;

        if (n_pert>0) {

            Float time = _time;

            auto* pert_adr = _perturber.neighbor_address.getDataAddress();

            Float xp[n_pert][3], xcm[3], m[n_pert];
            ChangeOver* changeover[n_pert_single+1];
            H4::NBAdr<PtclHard>::Group* ptclgroup[n_pert_group+1];

            int n_single_count=0;
            int n_group_count=0;
            for (int j=0; j<n_pert; j++) {
                H4::NBAdr<PtclHard>::Single* pertj;
                int k; // index of predicted data
                if (pert_adr[j].type==H4::NBType::group) {
                    pertj = &(((H4::NBAdr<PtclHard>::Group*)pert_adr[j].adr)->cm);
                    k = n_group_count + n_pert_single;
                    ptclgroup[n_group_count] = (H4::NBAdr<PtclHard>::Group*)pert_adr[j].adr;
                    n_group_count++;
                }
                else {
                    pertj = (H4::NBAdr<PtclHard>::Single*)pert_adr[j].adr;
                    k = n_single_count;
                    changeover[n_single_count] = &pertj->changeover;
                    n_single_count++;
                }

                Float dt = time - pertj->time;
                //ASSERT(dt>=-1e-7);
                xp[k][0] = pertj->pos[0] + dt*(pertj->vel[0] + 0.5*dt*(pertj->acc0[0] + inv3*dt*pertj->acc1[0]));
                xp[k][1] = pertj->pos[1] + dt*(pertj->vel[1] + 0.5*dt*(pertj->acc0[1] + inv3*dt*pertj->acc1[1]));
                xp[k][2] = pertj->pos[2] + dt*(pertj->vel[2] + 0.5*dt*(pertj->acc0[2] + inv3*dt*pertj->acc1[2]));


                m[k] = pertj->mass;
            }
            ASSERT(n_single_count == n_pert_single);
            ASSERT(n_group_count == n_pert_group);

            Float dt = time - _particle_cm.time;
            //ASSERT(dt>=0.0);
            xcm[0] = _particle_cm.pos[0] + dt*(_particle_cm.vel[0] + 0.5*dt*(_particle_cm.acc0[0] + inv3*dt*_particle_cm.acc1[0]));
            xcm[1] = _particle_cm.pos[1] + dt*(_particle_cm.vel[1] + 0.5*dt*(_particle_cm.acc0[1] + inv3*dt*_particle_cm.acc1[1]));
            xcm[2] = _particle_cm.pos[2] + dt*(_particle_cm.vel[2] + 0.5*dt*(_particle_cm.acc0[2] + inv3*dt*_particle_cm.acc1[2]));


            Float acc_pert_cm[3]={0.0, 0.0, 0.0};
            Float mcm = 0.0;
            // if (_perturber.need_resolve_flag) {
            // calculate component perturbation
            for (int i=0; i<_n_particle; i++) {
                Float* acc_pert = _force[i].acc_pert;
                Float& pot_pert = _force[i].pot_pert;
                auto& pi = _particles[i];
                auto& chi = pi.changeover;
                acc_pert[0] = acc_pert[1] = acc_pert[2] = Float(0.0);
                pot_pert = 0.0;

                Float xi[3];
                xi[0] = pi.pos[0] + xcm[0];
                xi[1] = pi.pos[1] + xcm[1];
                xi[2] = pi.pos[2] + xcm[2];

                // single perturber
                for (int j=0; j<n_pert_single; j++) {
                    Float dr[3] = {xp[j][0] - xi[0],
                                   xp[j][1] - xi[1],
                                   xp[j][2] - xi[2]};
                    Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] + eps_sq;
                    Float r  = sqrt(r2);
                    Float k  = ChangeOver::calcAcc0WTwo(chi, *changeover[j], r);
                    Float r3 = r*r2;
                    Float gm = gravitational_constant*m[j];
                    Float gmor3 = gm/r3 * k;

                    acc_pert[0] += gmor3 * dr[0];
                    acc_pert[1] += gmor3 * dr[1];
                    acc_pert[2] += gmor3 * dr[2];

                }
                // group perturber
                for (int j=n_pert_single; j<n_pert; j++) {
                    Float dr[3] = {xp[j][0] - xi[0],
                                   xp[j][1] - xi[1],
                                   xp[j][2] - xi[2]};
                    Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] + eps_sq;
                    Float r  = sqrt(r2);
                    const int jk = j-n_pert_single;
                    auto* ptcl_mem = ptclgroup[jk]->getDataAddress();
                    Float mk = 0.0;
                    for (int k=0; k<ptclgroup[jk]->getSize(); k++) {
                        mk += ptcl_mem[k].mass * ChangeOver::calcAcc0WTwo(chi, ptcl_mem[k].changeover, r);
                    }
                    Float r3 = r*r2;
                    Float gmor3 = gravitational_constant*mk/r3;

                    acc_pert[0] += gmor3 * dr[0];
                    acc_pert[1] += gmor3 * dr[1];
                    acc_pert[2] += gmor3 * dr[2];

                }

#ifdef EXTERNAL_HARD
                if (ext_force->isEnabled()) {
            
                    auto pgi = pi; 
                    auto gcm = _perturber.global_cm;

                    pgi.pos[0] += xcm[0] + gcm->pos[0];
                    pgi.pos[1] += xcm[1] + gcm->pos[1];
                    pgi.pos[2] += xcm[2] + gcm->pos[2];
                    // Here pi.vel has half kick step delay
                    pgi.vel[0] += pi.vel[0] + _particle_cm.vel[0] + dt*(_particle_cm.acc0[0] + 0.5*dt*_particle_cm.acc1[0]) + gcm->vel[0];
                    pgi.vel[1] += pi.vel[1] + _particle_cm.vel[1] + dt*(_particle_cm.acc0[1] + 0.5*dt*_particle_cm.acc1[1]) + gcm->vel[1];
                    pgi.vel[2] += pi.vel[2] + _particle_cm.vel[2] + dt*(_particle_cm.acc0[2] + 0.5*dt*_particle_cm.acc1[2]) + gcm->vel[2];

                    ext_force->calcAccJerkExternal(acc_pert, NULL, pgi, false);
                }
#endif

                acc_pert_cm[0] += pi.mass *acc_pert[0];
                acc_pert_cm[1] += pi.mass *acc_pert[1];
                acc_pert_cm[2] += pi.mass *acc_pert[2];

                mcm += pi.mass;

            }
//#ifdef AR_DEBUG
//            ASSERT(abs(mcm-_particle_cm.mass)<1e-10);
//#endif
                
            // （n_pert>0 can add corss term in PN with oher particle) (treat PN at same status as external hard)

            // get cm perturbation (exclude soft pert)
            acc_pert_cm[0] /= mcm;
            acc_pert_cm[1] /= mcm;
            acc_pert_cm[2] /= mcm;

            // remove cm. perturbation
            for (int i=0; i<_n_particle; i++) {
                Float* acc_pert = _force[i].acc_pert;
                Float& pot_pert = _force[i].pot_pert;
                const auto& pi = _particles[i];
                acc_pert[0] -= acc_pert_cm[0];
                acc_pert[1] -= acc_pert_cm[1];
                acc_pert[2] -= acc_pert_cm[2];
                pot_pert -= acc_pert[0]*pi.pos[0] + acc_pert[1]*pi.pos[1] + acc_pert[2]*pi.pos[2];

#ifdef SOFT_PERT
                if(_perturber.soft_pert!=NULL) {
                    // avoid too large perturbation force if system is disruptted
                    if (pi.pos*pi.pos<pi.changeover.getRout()*pi.changeover.getRout()) {
                        _perturber.soft_pert->eval(acc_pert, pi.pos);
                        pot_pert += _perturber.soft_pert->evalPot(pi.pos);
                    }
                }
#endif
            }

        }
        else {

            for(int i=0; i<_n_particle; i++) {
                Float* acc_pert = _force[i].acc_pert;
                Float& pot_pert = _force[i].pot_pert;
                acc_pert[0] = acc_pert[1] = acc_pert[2] = pot_pert = Float(0.0);
            }

#ifdef EXTERNAL_HARD
            if (ext_force->isEnabled()) {
            
                Float dt = _time - _particle_cm.time;
                Float xcm[3],vcm[3];
                //ASSERT(dt>=0.0);
                xcm[0] = _particle_cm.pos[0] + dt*(_particle_cm.vel[0] + 0.5*dt*(_particle_cm.acc0[0] + inv3*dt*_particle_cm.acc1[0]));
                xcm[1] = _particle_cm.pos[1] + dt*(_particle_cm.vel[1] + 0.5*dt*(_particle_cm.acc0[1] + inv3*dt*_particle_cm.acc1[1]));
                xcm[2] = _particle_cm.pos[2] + dt*(_particle_cm.vel[2] + 0.5*dt*(_particle_cm.acc0[2] + inv3*dt*_particle_cm.acc1[2]));
                vcm[0] = _particle_cm.vel[0] + dt*(_particle_cm.acc0[0] + 0.5*dt*_particle_cm.acc1[0]);
                vcm[1] = _particle_cm.vel[1] + dt*(_particle_cm.acc0[1] + 0.5*dt*_particle_cm.acc1[1]);
                vcm[2] = _particle_cm.vel[2] + dt*(_particle_cm.acc0[2] + 0.5*dt*_particle_cm.acc1[2]);

                Float acc_pert_cm[3]={0.0, 0.0, 0.0};
                Float mcm = 0.0;
                for (int i=0; i<_n_particle; i++) {
                    
                    const auto& pi = _particles[i];
                    auto pgi = pi; 
                    auto gcm = _perturber.global_cm;

                    pgi.pos[0] += xcm[0] + gcm->pos[0];
                    pgi.pos[1] += xcm[1] + gcm->pos[1];
                    pgi.pos[2] += xcm[2] + gcm->pos[2];
                    // Here pi.vel has half kick step delay
                    pgi.vel[0] += vcm[0] + gcm->vel[0];
                    pgi.vel[1] += vcm[1] + gcm->vel[1];
                    pgi.vel[2] += vcm[2] + gcm->vel[2];

                    Float* acc_pert = _force[i].acc_pert;
                    ext_force->calcAccJerkExternal(acc_pert, NULL, pgi, false);

                    acc_pert_cm[0] += pi.mass *acc_pert[0];
                    acc_pert_cm[1] += pi.mass *acc_pert[1];
                    acc_pert_cm[2] += pi.mass *acc_pert[2];

                    mcm += pi.mass;
                }
                acc_pert_cm[0] /= mcm;
                acc_pert_cm[1] /= mcm;
                acc_pert_cm[2] /= mcm;

                // remove cm. perturbation
                for (int i=0; i<_n_particle; i++) {
                    Float* acc_pert = _force[i].acc_pert;
                    Float& pot_pert = _force[i].pot_pert;
                    const auto& pi = _particles[i];
                    acc_pert[0] -= acc_pert_cm[0]; 
                    acc_pert[1] -= acc_pert_cm[1];        
                    acc_pert[2] -= acc_pert_cm[2]; 
                
                    pot_pert -= acc_pert[0]*pi.pos[0] + acc_pert[1]*pi.pos[1] + acc_pert[2]*pi.pos[2];
                }
            }
#endif

#ifdef SOFT_PERT
            if(_perturber.soft_pert!=NULL) {
                for(int i=0; i<_n_particle; i++) {
                    Float* acc_pert = _force[i].acc_pert;
                    Float& pot_pert = _force[i].pot_pert;
                    const auto& pi = _particles[i];
                    //acc_pert[0] = acc_pert[1] = acc_pert[2] = pot_pert = Float(0.0);
                    // avoid too large perturbation force if system is disruptted
                    if (pi.pos*pi.pos<pi.changeover.getRout()*pi.changeover.getRout()) {
                        _perturber.soft_pert->eval(acc_pert, pi.pos);
                        pot_pert += _perturber.soft_pert->evalPot(pi.pos);
                    }
                }
            }
#endif
        }

#ifdef SDAR_PN
        for (int i=0; i<_n_particle; i++) {
            const auto& pi = _particles[i];
            for (int j=i+1; j<_n_particle; j++) {
                const auto& pj = _particles[j];
                Float dr[3] = {pj.pos[0]-pi.pos[0],
                               pj.pos[1]-pi.pos[1],
                               pj.pos[2]-pi.pos[2]};
                Float dv[3] = {pj.vel[0]-pi.vel[0],
                               pj.vel[1]-pi.vel[1],
                               pj.vel[2]-pi.vel[2]};
                Float v2 = dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2];
                bool used_pn_orders[6] = {false,false,false,false,false,false};
                if (pn.setUsedPNOrders(used_pn_orders, v2, pi.mass, pj.mass)) {
                    Float ai[6][3], aj[6][3];
                    pn.calcAccJerkPN(ai, aj, NULL, NULL, NULL, NULL, pi.mass, pj.mass, dr, dv, NULL, NULL, used_pn_orders, false);
                    pn.sumAccJerkPN(&_force[i].acc_pert[0], NULL, ai, NULL, 1);
                    pn.sumAccJerkPN(&_force[j].acc_pert[0], NULL, aj, NULL, 1);
                }
            }
        }
#endif

    }
    

    //! calculate perturbation from c.m. acceleration
    Float calcPertFromForcePot(const Float* _force, const Float& _pot) {
        Float force2 = _force[0]*_force[0]+_force[1]*_force[1]+_force[2]*_force[2];
#ifdef AR_SLOWDOWN_PERT_R4
        Float inv_r = -force2/_pot;
        return sqrt(force2)*inv_r*inv_r*inv_r/gravitational_constant;
#else
        return -force2/(_pot*gravitational_constant);
#endif
    }

    //! calculate perturbation from binary tree
    static Float calcPertFromBinary(const COMM::Binary& _bin) {
        Float apo = _bin.semi*(1.0+_bin.ecc);
        Float apo2 = apo*apo;
#ifdef AR_SLOWDOWN_PERT_R4
        return (_bin.m1*_bin.m2)/(apo2*apo2);
#else
        return (_bin.m1*_bin.m2)/(apo2*apo);
#endif
    }

    //! calculate perturbation from distance to perturber and masses of particle and perturber
    static Float calcPertFromMR(const Float _r, const Float _mp, const Float _mpert) {
        Float r2 = _r*_r;
#ifdef AR_SLOWDOWN_PERT_R4
        return _mp*_mpert/(r2*r2);
#else
        return (_mp*_mpert)/(r2*_r);
#endif
    }

#if (defined AR_SLOWDOWN_ARRAY) || (defined AR_SLOWDOWN_TREE)

    //! calculate slowdown timescale
    void calcSlowDownTimeScale(Float& _t_min_sq, const Float dv[3], const Float dr[3], const Float& r, const Float& gm) {

        Float r2 = r*r;
        Float v2 = dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2];
        Float drdv = dr[0]*dv[0] + dr[1]*dv[1] + dr[2]*dv[2];

        Float semi = 1.0/(2.0/r - v2/gm);
        //hyperbolic, directly use velocity v
        if (semi<0)
            _t_min_sq = std::min(_t_min_sq, r2/v2);
        else {
            Float ra_fact = (1 - r/semi);
            Float e2 = drdv*drdv/(gm*semi) + ra_fact*ra_fact; // ecc^2
            Float r_vrmax = semi*(1-e2);
            if (r<r_vrmax) {
                // avoid decrese of vr once the orbit pass, calculate vr max at cos(E)=e (r==semi*(1-e^2))
                // vr_max = sqrt(er*(drdv^2*er + r*vcr2^2))/(G(m1+m2)r)
                //        = e*sqrt[G(m1+m2)/(a*(1-e^2)]
                Float vrmax_sq = e2*gm/r_vrmax;
                //Float rv2 = r*v2;
                //Float er = 2*gm - rv2;
                //Float vcr2 = gm - rv2;
                //Float vrmax_sq = er*(drdv*drdv*er + r*vcr2*vcr2)/(gm*gm*r2);
                _t_min_sq = std::min(_t_min_sq, semi*semi/vrmax_sq);
            }
            else {
                // r/vr
                Float rovr = r2/abs(drdv);
                _t_min_sq = std::min(_t_min_sq, rovr*rovr);
            }
        }
    }

    //! calculate slowdown perturbation and timescale from particle j to particle i
    /*!
      @param[out] _pert_out: perturbation from particle j
      @param[out] _t_min_sq: timescale limit from particle j
      @param[in] _pi: particle i (cm of binary)
      @param[in] _pj: particle j
     */
    void calcSlowDownPertOne(Float& _pert_out, Float& _t_min_sq, const PtclHard& pi, const PtclHard& pj) {
        Float dr[3] = {pj.pos[0] - pi.pos[0],
                       pj.pos[1] - pi.pos[1],
                       pj.pos[2] - pi.pos[2]};
        Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
        Float r = sqrt(r2);
        _pert_out += calcPertFromMR(r, pi.mass, pj.mass);

#ifdef AR_SLOWDOWN_TIMESCALE
        Float dv[3] = {pj.vel[0] - pi.vel[0],
                       pj.vel[1] - pi.vel[1],
                       pj.vel[2] - pi.vel[2]};

        // identify whether hyperbolic or closed orbit
        Float gm = gravitational_constant*(pi.mass+pj.mass);

        calcSlowDownTimeScale(_t_min_sq, dv, dr, r, gm);
        // force dependent method
        // min sqrt(r^3/(G m))
        //Float gmor3 = (mp+mcm)*r*r2/(sdt->G*mp*mcm);
        //sdt->trf2_min =  std::min(sdt->trf2_min, gmor3);
#endif
    }

    //! (Necessary) calculate slowdown perturbation from external effect on inner binary
    /*!
      @param[in] _bin: binary tree of member particles
    */
    void calcSlowDownPertExt(Float& _pert_out, const AR::BinaryTree<PtclHard>& _bin) {
        // perturbation from binary tree
#ifdef SDAR_PN
        auto* p1 = _bin.getMember(0);
        auto& vel1 = p1->vel;
        auto& pos1 = p1->pos;
        auto* p2 = _bin.getMember(1);
        auto& vel2 = p2->vel;
        auto& pos2 = p2->pos;

        Float dv[3] = {vel2[0]-vel1[0],
                       vel2[1]-vel1[1],
                       vel2[2]-vel1[2]};
        Float dr[3] = {pos2[0]-pos1[0],
                       pos2[1]-pos1[1],
                       pos2[2]-pos1[2]};
        Float v2 = dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2];
        Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
        Float ratio = pn.calcPN1OverNewton(v2);
        Float apo = _bin.semi*(1.0+_bin.ecc);
        _pert_out += ratio*p1->mass*p2->mass/(r2*apo);
#endif
    }
    

    //! (Necessary) calculate slowdown perturbation and timescale
    /*!
      @param[out] _pert_out: perturbation
      @param[out] _t_min_sq: timescale limit
      @param[in] _time: physical time for prediction
      @param[in] _bin: binary tree of member particles
      @param[in] _particle_cm: center-of-mass particle
      @param[in] _perturber: pertuber container
    */
    void calcSlowDownPert(Float& _pert_out, Float& _t_min_sq, const Float& _time, const AR::BinaryTree<PtclHard>& _bin, const H4Ptcl& _particle_cm, const ARPerturber& _perturber) {
        static const Float inv3 = 1.0 / 3.0;

        const int n_pert = _perturber.neighbor_address.getSize();

        if (n_pert>0) {

            auto* pert_adr = _perturber.neighbor_address.getDataAddress();

            Float xp[3], xcm[3];
            Float dt = _time - _particle_cm.time;
            //ASSERT(dt>=0.0);
            xcm[0] = _particle_cm.pos[0] + dt*(_particle_cm.vel[0] + 0.5*dt*(_particle_cm.acc0[0] + inv3*dt*_particle_cm.acc1[0]));
            xcm[1] = _particle_cm.pos[1] + dt*(_particle_cm.vel[1] + 0.5*dt*(_particle_cm.acc0[1] + inv3*dt*_particle_cm.acc1[1]));
            xcm[2] = _particle_cm.pos[2] + dt*(_particle_cm.vel[2] + 0.5*dt*(_particle_cm.acc0[2] + inv3*dt*_particle_cm.acc1[2]));

            Float mcm = _particle_cm.mass;
            auto& chi = _particle_cm.changeover;

#ifdef AR_SLOWDOWN_TIMESCALE
            // velocity dependent method
            Float vp[3], vcm[3];

            vcm[0] = _particle_cm.vel[0] + dt*(_particle_cm.acc0[0] + 0.5*dt*_particle_cm.acc1[0]);
            vcm[1] = _particle_cm.vel[1] + dt*(_particle_cm.acc0[1] + 0.5*dt*_particle_cm.acc1[1]);
            vcm[2] = _particle_cm.vel[2] + dt*(_particle_cm.acc0[2] + 0.5*dt*_particle_cm.acc1[2]);
#endif

            for (int j=0; j<n_pert; j++) {
                H4::NBAdr<PtclHard>::Single* pertj;
                if (pert_adr[j].type==H4::NBType::group) pertj = &(((H4::NBAdr<PtclHard>::Group*)pert_adr[j].adr)->cm);
                else pertj = (H4::NBAdr<PtclHard>::Single*)pert_adr[j].adr;

                Float dt = _time - pertj->time;
                //ASSERT(dt>=0.0);
                xp[0] = pertj->pos[0] + dt*(pertj->vel[0] + 0.5*dt*(pertj->acc0[0] + inv3*dt*pertj->acc1[0]));
                xp[1] = pertj->pos[1] + dt*(pertj->vel[1] + 0.5*dt*(pertj->acc0[1] + inv3*dt*pertj->acc1[1]));
                xp[2] = pertj->pos[2] + dt*(pertj->vel[2] + 0.5*dt*(pertj->acc0[2] + inv3*dt*pertj->acc1[2]));

                Float mj = pertj->mass;

                auto& chj = pertj->changeover;

                Float dr[3] = {xp[0] - xcm[0],
                               xp[1] - xcm[1],
                               xp[2] - xcm[2]};

                Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] + eps_sq;
                Float r = sqrt(r2);
                Float k  = ChangeOver::calcAcc0WTwo(chi, chj, r);
                _pert_out += calcPertFromMR(r, mcm, k*mj);

#ifdef AR_SLOWDOWN_TIMESCALE
                // velocity dependent method
                vp[0] = pertj->vel[0] + dt*(pertj->acc0[0] + 0.5*dt*pertj->acc1[0]);
                vp[1] = pertj->vel[1] + dt*(pertj->acc0[1] + 0.5*dt*pertj->acc1[1]);
                vp[2] = pertj->vel[2] + dt*(pertj->acc0[2] + 0.5*dt*pertj->acc1[2]);

                Float dv[3] = {vp[0] - vcm[0],
                               vp[1] - vcm[1],
                               vp[2] - vcm[2]};

                // identify whether hyperbolic or closed orbit
                Float gm = gravitational_constant*(mcm+mj);

                calcSlowDownTimeScale(_t_min_sq, dv, dr, r, gm);
#endif
            }
        }

        // add soft perturbation
        _pert_out += _perturber.soft_pert_min;

        // external perturbation on inner binary
        calcSlowDownPertExt(_pert_out, _bin);
    }
#endif

    //! (Necessary) modify one particle function
    /*!
      @param[in] _p: particle
      @param[in] _time_now: current time (not physical time, NB unit)
      @param[in] _time_end: required evolved time (not physical time, do not directly use, NB unit)
      \return 0: no modification; 1: modify mass; 2: modify mass and velocity; 3: mass become zero
     */
    template <class Tparticle>
    int modifyOneParticle(Tparticle& _p, const Float& _time_now, const Float& _time_end) {
#ifdef STELLAR_EVOLUTION
        // sample of mass loss
        //if (_p.time_interrupt<_time_end) {
        //    _p.dm = -_p.mass*1e-4;
        //    _p.mass +=_p.dm;
        //    _p.time_interrupt = _time_end+1e-5;
        //    return true;
        //}
#ifdef BSE_BASE
        // SSE/BSE stellar evolution
        if (_p.time_interrupt<=_time_end&&stellar_evolution_option>0) {
            ASSERT(bse_manager.checkParams());

            int modify_flag = 1;

            // time_record and time_interrupt have offsets, thus use difference to obtain true dt
            Float dt = _time_end - _p.time_record;

            // evolve star
            StarParameterOut output;
            StarParameter star_bk = _p.star;
            int event_flag = bse_manager.evolveStar(_p.star, output, dt);

            // error
            if (event_flag<0) {
                std::cerr<<"SSE Error: ID= "<<_p.id;
                _p.star.print(std::cerr);
                output.print(std::cerr);
                std::cerr<<std::endl;
                DATADUMP("dump_sse_error");
                std::cout<<std::flush;
                std::cerr<<std::flush;
                abort();
            }

            // if expected time not reach, record actually evolved time
            double dt_miss = bse_manager.getDTMiss(output);
            _p.time_record += dt-dt_miss;

            // estimate next time to check
            _p.time_interrupt = std::min(_p.time_record + bse_manager.getTimeStepStar(_p.star), time_interrupt_max);

            // record mass change (if loss, negative)
            double dm = bse_manager.getMassLoss(output);
            _p.dm += dm;
            if (dm==0.0) modify_flag = 0;

            // change mass in main data
            _p.mass = bse_manager.getMass(_p.star);

            // set merger check radius
            _p.radius = bse_manager.getMergerRadius(_p.star);

            // type change
            if (stellar_evolution_write_flag&&event_flag>=1) {
#pragma omp critical
                {
                    fout_sse<<"Type_change ";
                    //bse_manager.printTypeChange(fout_sse, _p.star, output);
                    fout_sse<<std::setw(WRITE_WIDTH)<<_p.id;
                    star_bk.printColumnAscii(fout_sse, WRITE_WIDTH);
                    _p.star.printColumnAscii(fout_sse, WRITE_WIDTH);
                    //output.printColumnAscii(fout_sse, WRITE_WIDTH);
                    fout_sse<<std::endl;
                }
            }

            // add velocity change if exist
            if (event_flag==2) {
                double dv[3];
                double dvabs=bse_manager.getVelocityChange(dv, output);
                assert(dvabs>0);
                for (int k=0; k<3; k++) _p.vel[k] += dv[k];
                modify_flag = 2;
                if (stellar_evolution_write_flag) {
#pragma omp critical
                    {
                        fout_sse<<"SN_kick "
                                <<std::setw(WRITE_WIDTH)<<_p.id
                                <<std::setw(WRITE_WIDTH)<<dvabs*bse_manager.vscale;
                        _p.star.printColumnAscii(fout_sse, WRITE_WIDTH);
                        fout_sse<<std::endl;
                    }
                }
            }
            // if mass become zero, set to unused for removing
            if (_p.mass==0.0) {
                _p.group_data.artificial.setParticleTypeToUnused(); // necessary to identify particle to remove
                modify_flag = 3;
            }

            return modify_flag;
        }
        
#endif // BSE_BASE
#ifdef DISK_STAR_MERGER
        if (_p.getBinaryInterruptState() != BinaryInterruptState::delaycollision) {
            // call mass change function
            if (_p.time_interrupt<=_time_end) {
                int modify_flag = disk_star_merger_manager.calcMassChange(&_p, _time_end, time_interrupt_max);
                return modify_flag;
            }
            else return 0;
        }
#endif

#endif // STELLAR_EVOLUTION
        return 0;
    }

    //! (Necessary) modify the orbits and interrupt check
    /*! check the inner left binary whether their separation is smaller than particle radius sum and become close, if true, set one component stauts to merger with cm mass and the other unused with zero mass. Return the binary tree address
      @param[in] _bin_interrupt: interrupt binary information: adr: binary tree address; time_now: current physical time; time_end: integration finishing time; status: interrupt status: change, merge,none
      @param[in] _bin: binarytree to check iteratively
      \return 0: no modification; 2: modified; 3: destroyed
    */
    int modifyAndInterruptIter(AR::InterruptBinary<PtclHard>& _bin_interrupt, AR::BinaryTree<PtclHard>& _bin) {
        int modify_return = 0;
#ifdef STELLAR_EVOLUTION
        int modify_branch[2];
        if (_bin.getMemberN()>2) {
            for (int k=0; k<2; k++) {
                if (_bin.isMemberTree(k)) {
                    modify_branch[k] = modifyAndInterruptIter(_bin_interrupt, *_bin.getMemberAsTree(k));
                    modify_return = std::max(modify_return, modify_branch[k]);
                }
                else {
                    // if member is star, evolve single star 
                    bool evolve_single_flag = false;
#ifdef BSE_BASE
                    evolve_single_flag = (stellar_evolution_option>0);
#endif
#ifdef DISK_STAR_MERGER
                    evolve_single_flag = true;
#endif
                    if (evolve_single_flag) {
                        modify_branch[k] = modifyOneParticle(*_bin.getMember(k), _bin.getMember(k)->time_record, _bin_interrupt.time_now);
                        modify_return = std::max(modify_return, modify_branch[k]);
                        // if status not set, set to change
                        if (modify_branch[k]>0&&_bin_interrupt.status == AR::InterruptStatus::none) {
                            _bin_interrupt.status = AR::InterruptStatus::change;
                            _bin_interrupt.setBinaryTreeAddress(&_bin);
                        }
                    }
                }
            }
            // ensure to record the root binary tree to include all changed members, if only record binary information (interrupt_detection_option == 2), should not do this
            if (modify_branch[0]>0&&modify_branch[1]>0 && interrupt_detection_option!=2) {
                _bin_interrupt.setBinaryTreeAddress(&_bin);
            }
            if (_bin_interrupt.status == AR::InterruptStatus::destroy) {
                // if both branch has destroyed, set destroy status, otherwise set merge status
                if (!(modify_branch[0]==3&&modify_branch[1]==3))
                    _bin_interrupt.status = AR::InterruptStatus::merge;
            }
        }
        else {

#ifdef DISK_STAR_MERGER
            for (int k=0; k<2; k++) {
                modify_branch[k] = modifyOneParticle(*_bin.getMember(k), _bin.getMember(k)->time_record, _bin_interrupt.time_now);
                modify_return = std::max(modify_return, modify_branch[k]);
                // if status not set, set to change
                if (modify_branch[k]>0&&_bin_interrupt.status == AR::InterruptStatus::none) {
                    _bin_interrupt.status = AR::InterruptStatus::change;
                    _bin_interrupt.setBinaryTreeAddress(&_bin);
                }
            }
#endif  

            auto* p1 = _bin.getLeftMember();
            auto* p2 = _bin.getRightMember();

#ifdef BSE_BASE
            auto postProcess =[&](StarParameterOut* out, Float* pos_cm, Float*vel_cm, Float& semi, Float& ecc, int binary_type_final=0, double vkick[][4]=NULL) {
                // if status not set, set to change
                if (_bin_interrupt.status == AR::InterruptStatus::none)
                    _bin_interrupt.status = AR::InterruptStatus::change;

                // set return flag >0
                modify_return = 2;
                p1->time_record = _bin_interrupt.time_now - bse_manager.getDTMiss(out[0]);
                p2->time_record = _bin_interrupt.time_now - bse_manager.getDTMiss(out[1]);
                // estimate next time to check
                p1->time_interrupt = std::min(p1->time_record + bse_manager.getTimeStepBinary(p1->star, p2->star, semi, ecc, binary_type_final), time_interrupt_max);
                p2->time_interrupt = p1->time_interrupt;

#ifdef GR_PRECESSION
    // calculate GR precession angular frequency 
                auto calcOmegaGR = [&](Float _m1, Float _m2, Float _semi, Float _ecc)->Float {
                    const Float c = bse_manager.getSpeedOfLight();
                    const Float c2 = c*c;
                    const Float mean_motion  = sqrt(gravitational_constant*(_m1+_m2)/fabs(_semi*_semi*_semi));
                    return 3.0 * gravitational_constant * (_m1 + _m2)*mean_motion / (c2 * _semi * (1.0 - _ecc * _ecc));
                };



                // calculate GR precession timescale and choose appropriate interrupt timestep
                Float Omega_GR = calcOmegaGR(p1->mass, p2->mass, _bin.semi, _bin.ecc);
                std::cout<<"_bin.semi: "<<_bin.semi<<std::endl;
                std::cout<<"_bin.ecc: "<<_bin.ecc<<std::endl;
                std::cout<<"bse_manager.getSpeedOfLight(): "<<bse_manager.getSpeedOfLight()<<std::endl;
                Float time_precession = 6.28318530717958647692 / Omega_GR; // 2*pi / Omega_GR
                Float period = _bin.period;
                Float N_precession = time_precession / period; // how many orbits to rotate 2pi

                if (N_precession > 18.0) {
                    // 
                    p1->time_interrupt = p1->time_record + 0.1f/Omega_GR;
                    //std::min({p1->time_record + bse_manager.getTimeStepBinary(p1->star, p2->star, semi, ecc, binary_type_final), p1->time_record + 0.1f/Omega_GR, time_interrupt_max})
                    std::cout<<"bse_manager.getTimeStepBinary: "<<bse_manager.getTimeStepBinary(p1->star, p2->star, semi, ecc, binary_type_final)<<std::endl;
                    std::cout<< 'time_interrupt(0.1f/Omega_GR):'<<0.1f/Omega_GR<<std::endl;
                }
                else {
                    // precession timescale is short: use precession timescale but enforce a minimum
                    time_precession = std::max(time_precession, 1e-2f*period);
                    p1->time_interrupt = std::min({p1->time_record + bse_manager.getTimeStepBinary(p1->star, p2->star, semi, ecc, binary_type_final), p1->time_record + time_precession, time_interrupt_max});
                }
                p2->time_interrupt = p1->time_interrupt;

#endif

                // reset collision state since binary orbit changes
                if (p1->getBinaryInterruptState()== BinaryInterruptState::collision)
                    p1->setBinaryInterruptState(BinaryInterruptState::none);
                if (p2->getBinaryInterruptState()== BinaryInterruptState::collision)
                    p2->setBinaryInterruptState(BinaryInterruptState::none);

                // set binary status (this is done in new/end group in Hermite group info printing, should not be used here)
                //p1->setBinaryPairID(p2->id);
                //p2->setBinaryPairID(p1->id);
                p1->setBinaryInterruptState(static_cast<BinaryInterruptState>(binary_type_final));
                p2->setBinaryInterruptState(static_cast<BinaryInterruptState>(binary_type_final));

                // record mass change (if loss, negative)
                // dm is used to correct energy, thus must be correctly set, use += since it may change mass before merge
                p1->dm += bse_manager.getMassLoss(out[0]);
                p2->dm += bse_manager.getMassLoss(out[1]);

                // update masses
                p1->mass = bse_manager.getMass(p1->star);
                p2->mass = bse_manager.getMass(p2->star);

                // set merger check radius
                p1->radius = bse_manager.getMergerRadius(p1->star);
                p2->radius = bse_manager.getMergerRadius(p2->star);


                bool mass_zero_flag = false;

                // if both mass becomes zero, set destroy state
                if (p1->mass==0.0&&p2->mass==0.0) {
                    p1->group_data.artificial.setParticleTypeToUnused(); // necessary to identify particle to remove
                    p2->group_data.artificial.setParticleTypeToUnused(); // necessary to identify particle to remove
                    _bin_interrupt.status = AR::InterruptStatus::destroy;
                    modify_return = 3;
                    mass_zero_flag = true;
                }
                else {
                    // if mass become zero, set to unused for removing and merger status
                    if (p1->mass==0.0) {
                        p1->group_data.artificial.setParticleTypeToUnused(); // necessary to identify particle to remove
                        _bin_interrupt.status = AR::InterruptStatus::merge;
                        p1->setBinaryInterruptState(BinaryInterruptState::none);
                        p2->setBinaryInterruptState(BinaryInterruptState::none);
                        // set new particle position and velocity to be the original cm
                        for (int k=0; k<3; k++) {
                            p2->pos[k] = pos_cm[k];
                            p2->vel[k] = vel_cm[k];
                        }
                        modifyOneParticle(*p2, p2->time_record, _bin_interrupt.time_now);
                        mass_zero_flag = true;
                    }

                    if (p2->mass==0.0) {
                        p2->group_data.artificial.setParticleTypeToUnused(); // necessary to identify particle to remove
                        _bin_interrupt.status = AR::InterruptStatus::merge;
                        p1->setBinaryInterruptState(BinaryInterruptState::none);
                        p2->setBinaryInterruptState(BinaryInterruptState::none);
                        // set new particle position and velocity to be the original cm
                        for (int k=0; k<3; k++) {
                            p1->pos[k] = pos_cm[k];
                            p1->vel[k] = vel_cm[k];
                        }
                        modifyOneParticle(*p1, p1->time_record, _bin_interrupt.time_now);
                        mass_zero_flag = true;
                    }

                    // case when velocity kick appears
                    bool kick_flag = false;
                    if (vkick!=NULL) {
                        for (int k=0; k<2; k++) {
                            if (vkick[k][3]>0.0) {
                                auto* pk = _bin.getMember(k);
                                kick_flag = true;    
                                for (int i=0; i<3; i++) pk->vel[i] += vkick[k][i];
                            }
                        }
                    }

                    if (!kick_flag && !mass_zero_flag) {
                        // case for elliptic case
                        if (ecc>=0.0&&ecc<=1.0) {
                            // obtain full orbital parameters
                            //_bin.calcOrbit(gravitational_constant);
                            // update new period, ecc
//#pragma omp critical
//                            std::cerr<<"Event: "<<event_flag<<" "<<_bin.period<<" "<<period<<" "<<_bin.ecc<<" "<<ecc<<std::endl;
                            _bin.semi = semi;
                            ASSERT(_bin.semi>0);
                            _bin.ecc = ecc;
                            _bin.m1 = p1->mass;
                            _bin.m2 = p2->mass;
                            //if (((ecc-ecc_bk)/(1-ecc)>0.01||(period-period_bk)/period>1e-2)) {
                            // kepler orbit to particles using the same ecc anomaly
                            
                            _bin.calcParticles(gravitational_constant);
                            p1->pos += _bin.pos;
                            p2->pos += _bin.pos;
                            p1->vel += _bin.vel;
                            p2->vel += _bin.vel;
                            //}
                        }
                        // in case of disruption but no kick
                        else {
                            // obtain full orbital parameters
                            // _bin.calcOrbit(gravitational_constant);
                            // assume energy no change
                            if(_bin.semi>0) _bin.semi = -_bin.semi;
                            _bin.ecc = ecc;
                            _bin.m1 = p1->mass;
                            _bin.m2 = p2->mass;
                            ASSERT(ecc>=1.0);
                            // kepler orbit to particles using the same ecc anomaly
                            //if ((ecc-ecc_bk)/ecc>1e-6) {
                            _bin.calcParticles(gravitational_constant);
                            p1->pos += _bin.pos;
                            p2->pos += _bin.pos;
                            p1->vel += _bin.vel;
                            p2->vel += _bin.vel;
                            //}
                        }
                    }
                }
                if (mass_zero_flag) {
                    DATADUMP("dump_binary_merger");
                }

            };

            COMM::Vector3<Float> pos_red(p2->pos[0] - p1->pos[0], p2->pos[1] - p1->pos[1], p2->pos[2] - p1->pos[2]);
            COMM::Vector3<Float> vel_red(p2->vel[0] - p1->vel[0], p2->vel[1] - p1->vel[1], p2->vel[2] - p1->vel[2]);
            Float drdv = pos_red * vel_red;
            if (stellar_evolution_option>0) {
                int binary_type_p1 = static_cast<int>(p1->getBinaryInterruptState());
                int binary_type_p2 = static_cast<int>(p2->getBinaryInterruptState());
                int binary_type_init = 0;
                if (binary_type_p1==binary_type_p2) binary_type_init = binary_type_p1;

                // check whether need to update based on time step
                double time_check = std::min(p1->time_interrupt, p2->time_interrupt);
                Float dt1 = _bin_interrupt.time_now - p1->time_record;
                Float dt2 = _bin_interrupt.time_now - p2->time_record;
                // if next time to check > time_now, do not evolve by setting dt = 0;
                if (time_check>_bin_interrupt.time_now) dt1 = dt2 = 0.0;

                // check whether bse is needed
                Float dr[3] = {p1->pos[0] - p2->pos[0],
                               p1->pos[1] - p2->pos[1],
                               p1->pos[2] - p2->pos[2]};
                Float dr2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];

                bool check_flag = bse_manager.isCallBSENeeded(p1->star, p2->star, dr2, _bin.semi, _bin.ecc, dt1, dt2, binary_type_init);

	        	if (check_flag) {
                    ASSERT(bse_manager.checkParams());
                    // record address of modified binary
                    _bin_interrupt.setBinaryTreeAddress(&_bin);

                    // first evolve two components to the same starting time
                    if (p1->time_record!=p2->time_record) {
                        if (p1->time_record<p2->time_record) {
                            p1->time_interrupt = p1->time_record;
                            modifyOneParticle(*p1, p1->time_record, p2->time_record);
                        }
                        else {
                            p2->time_interrupt = p2->time_record;
                            modifyOneParticle(*p2, p2->time_record, p1->time_record);
                        }
                    }
                    Float dt = _bin_interrupt.time_now - std::max(p1->time_record,p2->time_record);
                    ASSERT(dt>0);

                    StarParameterOut out[2];
                    _bin.calcOrbit(gravitational_constant);
                    Float ecc = _bin.ecc;
                    //Float ecc_bk = ecc;
                    Float semi = _bin.semi;
                    //Float semi_bk =semi;
                    Float mtot = p1->mass+p2->mass;
                    Float period = _bin.period;
                    //Float period_bk = period;

                    // backup c.m. information
                    Float pos_cm[3], vel_cm[3];
                    for (int k=0; k<3; k++) {
                        pos_cm[k] = (p1->mass*p1->pos[k] + p2->mass*p2->pos[k])/mtot;
                        vel_cm[k] = (p1->mass*p1->vel[k] + p2->mass*p2->vel[k])/mtot;
                    }
                    // backup star
                    StarParameter p1_star_bk = p1->star;
                    StarParameter p2_star_bk = p2->star;

                    BinaryEvent bin_event;
                    // loop until the time_end reaches
                    int event_flag = bse_manager.evolveBinary(p1->star, p2->star, out[0], out[1], semi, period, ecc, &_bin.am.x, &pos_red.x, bin_event, binary_type_init, dt);

                    // error
                    if (event_flag<0) {
                        std::cerr<<"BSE Error! ";
                        std::cerr<<" ID="<<p1->id<<" "<<p2->id<<" ";
                        std::cerr<<" semi[R*]: "
                                 <<_bin.semi*bse_manager.rscale
                                 <<" ecc: "<<_bin.ecc
                                 <<" period[days]: "<<_bin.period*bse_manager.tscale*bse_manager.year_to_day;
                        std::cerr<<" Init: Star1: ";
                        p1_star_bk.print(std::cerr);
                        std::cerr<<" Star2: ";
                        p2_star_bk.print(std::cerr);
                        std::cerr<<" final: Star1: ";
                        p1->star.print(std::cerr);
                        out[0].print(std::cerr);
                        std::cerr<<" Star2: ";
                        p2->star.print(std::cerr);
                        out[1].print(std::cerr);
                        std::cerr<<std::endl;
                        DATADUMP("dump_bse_error");
                        std::cout<<std::flush;
                        std::cerr<<std::flush;
                        abort();
                    }

                    // check binary type and print event information
                    int binary_type_final=0;
                    int nmax = bin_event.getEventNMax();
                    int binary_type_init = bin_event.getType(bin_event.getEventIndexInit());
                    for (int i=0; i<nmax; i++) {
                        int binary_type = bin_event.getType(i);
                        if (binary_type>0) {
                            bool first_event = (i==0);
                            if (stellar_evolution_write_flag) {
                                if ((first_event&&binary_type_init!=binary_type)||!first_event) {
                                    //if (!(binary_type_init==11&&(binary_type==3||binary_type==11))) {// avoid repeating printing Start Roche and BSS
#pragma omp critical
                                    {
                                        bse_manager.printBinaryEventColumnOne(fout_bse, bin_event, i, WRITE_WIDTH);
                                        fout_bse<<std::setw(WRITE_WIDTH)<<p1->id
                                                <<std::setw(WRITE_WIDTH)<<p2->id
                                                <<std::setw(WRITE_WIDTH)<<drdv*bse_manager.rscale*bse_manager.vscale
                                                <<std::setw(WRITE_WIDTH)<<_bin.r*bse_manager.rscale;
                                        fout_bse<<std::endl;
                                    }
                                }
                            }
                            //if (binary_type==10) {
                            //    DATADUMP("co_dump");
                            //    abort();
                            //}
                            //if (binary_type==3&&p1->time_record>=2499.0114383817) {
                            //    DATADUMP("re_dump");
                            //    abort();
                            //}

                            //if (vkick[3]>0||vkick[7]>0) event_flag = 3; // kick
                            if (binary_type>0) event_flag = std::max(event_flag, 1); // type change
                            if (bse_manager.isMassTransfer(binary_type)) event_flag = std::max(event_flag, 2); // orbit change
                            else if (bse_manager.isDisrupt(binary_type)) event_flag = std::max(event_flag, 3); // disrupt
                            else if (bse_manager.isMerger(binary_type)) {
                                event_flag = std::max(event_flag, 4); // Merger
                                if (bse_manager.isGWMerger(binary_type)) event_flag = std::max(event_flag, 6); // GW Merger
                            }
                            else if (bse_manager.isNoRemnant(binary_type)) event_flag = std::max(event_flag, 5); // no Remnant
                            binary_type_final = binary_type;

                        }
                        else if(binary_type<0) break;
                    }

                    // check event_flag
                    if (event_flag<=2) {
                        ASSERT(bse_manager.getMass(p1->star)>0 && bse_manager.getMass(p2->star)>0); 
                    }

                    // check SN event
                    Float vkick[2][4];
                    for (int k=0; k<2; k++) {
                        auto* pk = _bin.getMember(k);
                        vkick[k][3] = bse_manager.getVelocityChange(vkick[k], out[k]);    
                        if (vkick[k][3]>0) {
#pragma omp critical 
                            {
                                if (event_flag==6) fout_bse<<"GW_kick ";
                                else fout_bse<<"SN_kick ";
                                fout_bse<<std::setw(WRITE_WIDTH)<<p1->id
                                        <<std::setw(WRITE_WIDTH)<<p2->id
                                        <<std::setw(WRITE_WIDTH)<<k+1
                                        <<std::setw(WRITE_WIDTH)<<vkick[k][3]*bse_manager.vscale;
                                pk->star.printColumnAscii(fout_bse, WRITE_WIDTH);
                                fout_bse<<std::endl;
                            }
                        }
                    }

        

                    // update semi
                    mtot = bse_manager.getMass(p1->star) + bse_manager.getMass(p2->star);
                    semi = COMM::Binary::periodToSemi(period, mtot, gravitational_constant);

#ifdef GR_PRECESSION
                    auto calcOmegaGR = [&](Float _m1, Float _m2, Float _semi, Float _ecc)->Float {
                    const Float c = bse_manager.getSpeedOfLight();
                    const Float c2 = c*c;
                    const Float mean_motion  = sqrt(gravitational_constant*(_m1+_m2)/fabs(_semi*_semi*_semi));
                    return 3.0 * gravitational_constant * (_m1 + _m2)*mean_motion / (c2 * _semi * (1.0 - _ecc * _ecc));
                    };
                    Float Omega_GR = calcOmegaGR(p1->mass, p2->mass, semi, ecc);
                    _bin.rot_self += Omega_GR*dt;


                    std::cout<<"dt"<<dt<<std::endl;
                    std::cout<<"rot_self: "<<_bin.rot_self<<std::endl;
                    std::cout<<"_bin_interrupt.time_now: "<<_bin_interrupt.time_now<<std::endl;
                    std::cout<<"p1->time_record "<<p1->time_record<<std::endl;
                    std::cout<<"p1->time_interrupt "<<p1->time_interrupt<<std::endl;




#endif

                    // change p1 and p2 due to output from stellar evolution
                    postProcess(out, pos_cm, vel_cm, semi, ecc, binary_type_final, vkick);
                }
            }
#endif // BSE_BASE

            // dynamical merger and tide check
            if (_bin_interrupt.status!=AR::InterruptStatus::merge&&_bin_interrupt.status!=AR::InterruptStatus::destroy) {

                auto merge = [&](const Float& dr, const Float& t_peri, const Float& sd_factor, std::string logmessage = "Dynamic_merge" ) {
                    _bin_interrupt.setBinaryTreeAddress(&_bin);

#ifdef BSE_BASE
                    //Float m1_bk = p1->mass;
                    //Float m2_bk = p2->mass;
                    // backup original data for print

                    // first evolve two components to the current time
                    if (stellar_evolution_option>0) {
                        if (p1->time_record<_bin_interrupt.time_now) {
                            p1->time_interrupt = p1->time_record; // force SSE evolution
                            modifyOneParticle(*p1, p1->time_record, _bin_interrupt.time_now);
                        }
                        if (p2->time_record<_bin_interrupt.time_now) {
                            p2->time_interrupt = p2->time_record; // force SSE evolution
                            modifyOneParticle(*p2, p2->time_record, _bin_interrupt.time_now);
                        }

                        ASSERT(p1->star.tphys==p2->star.tphys);
                        //ASSERT(p1->star.mass>0&&p2->star.mass>0); // one may have SNe before merge

                        StarParameter p1_star_bk, p2_star_bk;
                        StarParameterOut out[2];
                        Float pos_cm[3], vel_cm[3];
                        Float mtot = p1->mass+p2->mass;

                        for (int k=0; k<3; k++) {
                            pos_cm[k] = (p1->mass*p1->pos[k] + p2->mass*p2->pos[k])/mtot;
                            vel_cm[k] = (p1->mass*p1->vel[k] + p2->mass*p2->vel[k])/mtot;
                        }

                        p1_star_bk = p1->star;
                        p2_star_bk = p2->star;
                        // call BSE function to merge two stars
                        Float semi = _bin.semi;
                        Float ecc = _bin.ecc;
                        bse_manager.merge(p1->star, p2->star, out[0], out[1], semi, ecc);

                        postProcess(out, pos_cm, vel_cm, semi, ecc);
                        if (stellar_evolution_write_flag&&(p1->mass==0.0||p2->mass==0.0)) {
#pragma omp critical
                            {
                                fout_bse<<logmessage<<" "
                     			        <<std::setw(WRITE_WIDTH)<<p1->id
                                        <<std::setw(WRITE_WIDTH)<<p2->id
                                        <<std::setw(WRITE_WIDTH)<<_bin.period*bse_manager.tscale*bse_manager.year_to_day
                                        <<std::setw(WRITE_WIDTH)<<_bin.semi*bse_manager.rscale
                                        <<std::setw(WRITE_WIDTH)<<_bin.ecc;
#ifndef DYNAMIC_MERGER_LESS_OUTPUT
                                fout_bse<<std::setw(WRITE_WIDTH)<<dr*bse_manager.rscale
                                        <<std::setw(WRITE_WIDTH)<<t_peri*bse_manager.tscale*bse_manager.year_to_day
                                        <<std::setw(WRITE_WIDTH)<<sd_factor;
#endif
                                // before
                                p1_star_bk.printColumnAscii(fout_bse, WRITE_WIDTH);
                                p2_star_bk.printColumnAscii(fout_bse, WRITE_WIDTH);
                                // after
                                p1->star.printColumnAscii(fout_bse, WRITE_WIDTH);
                                p2->star.printColumnAscii(fout_bse, WRITE_WIDTH);
                                fout_bse<<std::endl;

                                std::string dump_name = "dump_" + logmessage;
                                DATADUMP(dump_name.c_str());
                            }
                        }
                    }
#else //not BSE_BASE
                    // print data
#pragma omp critical
                    {
                        _bin_interrupt.printColumnAscii(fout_interrupt, WRITE_WIDTH, true);
                        fout_interrupt<<std::endl;

                        DATADUMP("dump_interrupt"); 
                    }

                    // set return flag >0
                    modify_return = 2;

                    // merge two particles
                    if (interrupt_detection_option == 1) {
                        p1->time_record = _bin_interrupt.time_now;
                        p2->time_record = _bin_interrupt.time_now;

#ifdef DISK_STAR_MERGER
                        // use disk star merger
                        disk_star_merger_manager.calcMergerProperties(p1, p2, _bin_interrupt.time_now);

                        if (p1->mass ==0.0) p1->group_data.artificial.setParticleTypeToUnused(); // necessary to identify particle to remove
                        if (p2->mass ==0.0) p2->group_data.artificial.setParticleTypeToUnused(); // necessary to identify particle to remove
#else
                        // merge two particles
                        Float mcm = p1->mass + p2->mass;
                        for (int k=0; k<3; k++) {
                            p1->pos[k] = (p1->mass*p1->pos[k] + p2->mass*p2->pos[k])/mcm;
                            p1->vel[k] = (p1->mass*p1->vel[k] + p2->mass*p2->vel[k])/mcm;
                        }
                        p1->dm += p2->mass;
                        p2->dm -= p2->mass;

                        p1->mass = mcm;
                        p2->mass = 0.0;

                        p2->radius = 0.0;
                        p1->mass += p2->mass;

                        p2->group_data.artificial.setParticleTypeToUnused(); // necessary to identify particle to remove
#endif

                        // reset collision state since binary orbit changes
                        p1->setBinaryInterruptState(BinaryInterruptState::none);
                        p2->setBinaryInterruptState(BinaryInterruptState::none);

                        if (_bin_interrupt.status != AR::InterruptStatus::destroy) 
                            _bin_interrupt.status = AR::InterruptStatus::merge;
                    }
                    // record particle information, only set status
                    else if (interrupt_detection_option == 2) {
                        p1->setBinaryInterruptState(BinaryInterruptState::collision);
                        p2->setBinaryInterruptState(BinaryInterruptState::collision);
                    }

#endif // end BSE_BASE
                    //p1->setBinaryPairID(0);
                    //p2->setBinaryPairID(0);
                };

#ifndef BSE_BASE
                // delayed merger
                if (p1->getBinaryInterruptState()== BinaryInterruptState::delaycollision && 
                    p2->getBinaryInterruptState()== BinaryInterruptState::delaycollision &&
                    (p1->time_interrupt<_bin_interrupt.time_now && p2->time_interrupt<_bin_interrupt.time_now) &&
                    (p1->getBinaryPairID()==p2->id||p2->getBinaryPairID()==p1->id)) {
                    Float dr[3] = {p1->pos[0] - p2->pos[0], 
                                   p1->pos[1] - p2->pos[1], 
                                   p1->pos[2] - p2->pos[2]};
                    Float dr2  = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
                    merge(std::sqrt(dr2), 0.0, 1.0);
                }
                else if (p1->getBinaryInterruptState() != BinaryInterruptState::collision && p2->getBinaryInterruptState() != BinaryInterruptState::collision) {
                    // check merger
                    Float radius = p1->radius + p2->radius;
                    // slowdown case
                    if (_bin.slowdown.getSlowDownFactor()>1.0) {
                        ASSERT(_bin.semi>0.0);
                        Float drdv;
                        _bin.particleToSemiEcc(_bin.semi, _bin.ecc, _bin.r, drdv, *_bin.getLeftMember(), *_bin.getRightMember(), gravitational_constant);
                        Float peri = _bin.semi*(1 - _bin.ecc);
                        if (peri<radius) {
                            Float ecc_anomaly  = _bin.calcEccAnomaly(_bin.r);
                            Float mean_anomaly = _bin.calcMeanAnomaly(ecc_anomaly, _bin.ecc);
                            Float mean_motion  = sqrt(gravitational_constant*_bin.mass/(fabs(_bin.semi*_bin.semi*_bin.semi)));
                            Float t_peri = mean_anomaly/mean_motion;
                            if (drdv<0 && t_peri<_bin_interrupt.time_end-_bin_interrupt.time_now) {
                                Float dr[3] = {p1->pos[0] - p2->pos[0],
                                               p1->pos[1] - p2->pos[1],
                                               p1->pos[2] - p2->pos[2]};
                                Float dr2  = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
                                _bin_interrupt.time_now += t_peri;
                                merge(std::sqrt(dr2), t_peri, _bin.slowdown.getSlowDownFactor());
                            }
                            else if (_bin.semi>0||(_bin.semi<0&&drdv<0)) {
                                // ensure to set pair id for delayed collision
                                //p1->setBinaryPairID(p2->id);
                                //p2->setBinaryPairID(p1->id);
                                p1->setBinaryInterruptState(BinaryInterruptState::delaycollision);
                                p2->setBinaryInterruptState(BinaryInterruptState::delaycollision);
                                p1->time_interrupt = std::min(_bin_interrupt.time_now + drdv<0 ? t_peri : (_bin.period - t_peri), time_interrupt_max);
                                p2->time_interrupt = p1->time_interrupt;

                            }
                        }
                    }
                    else { // no slowdown case, check separation directly
                        Float dr[3] = {p1->pos[0] - p2->pos[0],
                                       p1->pos[1] - p2->pos[1],
                                       p1->pos[2] - p2->pos[2]};
                        Float dr2  = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
                        if (dr2<radius*radius) merge(std::sqrt(dr2), 0.0, 1.0);
                    }
                }

#else // BSE_BASE
                // Check TDE 
                //This my revised implementation (GIU), we do not distinguish a priori about the
                //even type since the class _bin can handle both close and open orbit.
                //Notice that this class is defined in SDAR::binary_trww
                //One thing I don't uderstand is if we have to check for the slowdown as in the previous section

                if (p1->getBinaryInterruptState()== BinaryInterruptState::tde && 
                    p2->getBinaryInterruptState()== BinaryInterruptState::tde &&
                    (p1->time_interrupt<_bin_interrupt.time_now && p2->time_interrupt<_bin_interrupt.time_now) &&
                    (p1->getBinaryPairID()==p2->id||p2->getBinaryPairID()==p1->id)) {
                    Float dr[3] = {p1->pos[0] - p2->pos[0], 
                                   p1->pos[1] - p2->pos[1], 
                                   p1->pos[2] - p2->pos[2]};
                    Float dr2  = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
                    merge(std::sqrt(dr2), 0.0, 1.0, "Binary_TDE");
                }
                else {
                    Float r_tde = bse_manager.getTDERadius(p1->star, p2->star);
                    if (r_tde>0) {
                        /********* Hyperbolic and binary TDE ****************/
                        //r_tde>r_peri, at certain point along the orbit the distance of the two objects is whitin the tidal radius
                        Float r_peri = _bin.semi * (1-_bin.ecc);
                        if (r_tde > r_peri){
                            //Estimate distance at this point //Petar units
                            Float dr[3] = {p1->pos[0] - p2->pos[0],
                                        p1->pos[1] - p2->pos[1],
                                        p1->pos[2] - p2->pos[2]};
                            Float dr2  = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];

                            //check if the two stars are already within the tidal radius
                            if (dr2 < r_tde*r_tde){
                                //Check if in hyperbolic or close orbit
                                if (_bin.semi<0){merge(std::sqrt(dr2), 0 , _bin.slowdown.getSlowDownFactor(), "Hyperbolic_TDE");}
                                else if (_bin.semi>0){merge(std::sqrt(dr2), 0,  _bin.slowdown.getSlowDownFactor(), "Binary_TDE");}
                            }
                        
                            //check if the star is approaching the BH/NS (drdv<0)
                            if (drdv<0) {
                                //Here we have three choice:
                                //A: if this conditiond is satisfied always merge (but this can create situation in which
                                //we create an instantenous TDE of very distant objects, very poor approximation)
                                //if  (_bin.semi<0) merge(std::sqrt(dr2), 0.0, _bin.slowdown.getSlowDownFactor(), "Hyperbolic_TDE_A: ");
                                //else if (_bin.semi>0) merge(std::sqrt(dr2), 0.0, _bin.slowdown.getSlowDownFactor(), "Binary_TDE_A: ");

                                //B: A simple improvement merge only if the current distance is withinn 3 times the r_tde distance
                                //if (dr2<9*r_tde*r_tde){
                                //merge(std::sqrt(dr2), 0.0, _bin.slowdown.getSlowDownFactor(), "Hyperbolic_TDE: ");
                                //}

                                //C: More realistic option, estimate the time needed to reach the pericentre and check
                                //if this is within the current check timestep
                                Float ecc_anomaly  = _bin.calcEccAnomaly(_bin.r);
                                Float mean_anomaly = _bin.calcMeanAnomaly(ecc_anomaly, _bin.ecc);
                                Float mean_motion  = sqrt(gravitational_constant*_bin.mass/(fabs(_bin.semi*_bin.semi*_bin.semi)));
                                Float t_peri = mean_anomaly/mean_motion;

                                if (t_peri<_bin_interrupt.time_end-_bin_interrupt.time_now) {
                                    //Check if in hyperbolic or close orbit
                                    if (_bin.semi<0){merge(std::sqrt(dr2), t_peri, _bin.slowdown.getSlowDownFactor(), "Hyperbolic_TDE");}
                                    else if (_bin.semi>0){merge(std::sqrt(dr2), t_peri, _bin.slowdown.getSlowDownFactor(), "Binary_TDE");}
                                }
                                else if(_bin.slowdown.getSlowDownFactor()>1.0){
                                    p1->setBinaryPairID(p2->id);
                                    p2->setBinaryPairID(p1->id);
                                    p1->setBinaryInterruptState(BinaryInterruptState::tde);
                                    p2->setBinaryInterruptState(BinaryInterruptState::tde);
                                    p1->time_interrupt = std::min(_bin_interrupt.time_now + drdv<0 ? t_peri : (_bin.period - t_peri), time_interrupt_max);
                                    p2->time_interrupt = p1->time_interrupt;
                                }
                                //D: In the most realistic we shoud check the time needed to reach the position r=r_tde, not r=rp
                            }
                        }
                    }
                }

                // Check hyperbolic merger
                // in bse case, handle binary merger in bse, only check hyperbolic merger
                if (_bin.semi<0.0) {
                    // check merger
                    Float radius = p1->radius + p2->radius;
                    
                    Float dr[3] = {p1->pos[0] - p2->pos[0], 
                        p1->pos[1] - p2->pos[1], 
                        p1->pos[2] - p2->pos[2]};
                    Float dr2  = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
                    if (dr2<radius*radius) merge(std::sqrt(dr2), 0.0, 1.0);
                }

        //             // LSO merger check
        //                 if (_bin.semi > 0.0 && _bin.ecc < 1.0 && p1->mass > 0 && p2->mass > 0) {
        //                     Float lso = 6.0 * (p1->radius + p2->radius);
        //                     if (_bin.semi * (1.0 - _bin.ecc) < lso) {
        //                         Float dr[3] = {p1->pos[0] - p2->pos[0], 
        //                                     p1->pos[1] - p2->pos[1], 
        //                                     p1->pos[2] - p2->pos[2]};
        //                         Float dr2  = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
        //                         merge(std::sqrt(dr2), 0.0, 1.0);
                                
        // #pragma omp critical
        //                         {
        //                             fout_bse<<"LSO_merger "
        //                                     <<std::setw(WRITE_WIDTH)<<_bin_interrupt.time_now
        //                                     <<std::setw(WRITE_WIDTH)<<p1->id
        //                                     <<std::setw(WRITE_WIDTH)<<p2->id
        //                                     <<std::setw(WRITE_WIDTH)<<_bin.semi
        //                                     <<std::setw(WRITE_WIDTH)<<_bin.ecc
        //                                     <<std::setw(WRITE_WIDTH)<<_bin.semi * (1.0 - _bin.ecc)
        //                                     <<std::setw(WRITE_WIDTH)<<lso
        //                                     <<std::endl;
        //                         }
        //                     // After merge, status could be merge or destroy, we should skip the tide part
        //                     if (_bin_interrupt.status == AR::InterruptStatus::merge || _bin_interrupt.status == AR::InterruptStatus::destroy) return modify_return;
        //                 }
        //             }

                // tide energy loss
                if (stellar_evolution_option==2 && p1->mass>0 && p2->mass>0) {
                    if (drdv<0) { // when two star approach each other; reset tide status
                        if (p1->getBinaryInterruptState() == BinaryInterruptState::tide) {
                            p1->setBinaryInterruptState(BinaryInterruptState::none);
                        }
                        if (p2->getBinaryInterruptState() == BinaryInterruptState::tide) {
                            p2->setBinaryInterruptState(BinaryInterruptState::none);
                        }
                    }
                    else { // modify orbit based on energy loss
                        int binary_type_p1 = static_cast<int>(p1->getBinaryInterruptState());
                        int binary_type_p2 = static_cast<int>(p2->getBinaryInterruptState());
                        long long int pair_id1 = p1->getBinaryPairID();
                        long long int pair_id2 = p2->getBinaryPairID();
                        bool tide_flag = true;
                        if ((binary_type_p1 != binary_type_p2) || (pair_id1 != p2->id) || (pair_id2 != p1->id)) tide_flag = false;
                        else if (bse_manager.isMassTransfer(binary_type_p1) 
                                 || bse_manager.isMerger(binary_type_p1) 
                                 || bse_manager.isNoRemnant(binary_type_p1) 
                                 || bse_manager.isDisrupt(binary_type_p1))
                            tide_flag = false;

                        bool change_flag=false;
                        if (tide_flag) {
                            Float poly_type1=0, poly_type2=0;
                            Float Etid=0, Ltid=0;
                            Float semi = _bin.semi;
                            Float ecc  = _bin.ecc;
                            Float rad1 = bse_manager.getStellarRadius(p1->star);
                            Float rad2 = bse_manager.getStellarRadius(p2->star);
                            if (p1->star.isCompactObject() && p2->star.isCompactObject()) {
                                if (_bin.semi<0) {
                                    bool merge_flag = tide.evolveOrbitHyperbolicGW(_bin, Etid, Ltid);
                                    if (merge_flag) {
                                        Float dr[3] = {p1->pos[0] - p2->pos[0],
                                                       p1->pos[1] - p2->pos[1],
                                                       p1->pos[2] - p2->pos[2]};
                                        Float dr2  = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
                                        merge(std::sqrt(dr2), 0.0, 1.0, "GW_tide_merge");
                                    }
                                    else change_flag = true;
                                }
                            }
                            else if (p1->star.isNotNSBH() || p2->star.isNotNSBH()) {
                                poly_type1 = (p1->star.isBeforeGiant()) ? 3.0 : 1.5;
                                poly_type2 = (p2->star.isBeforeGiant()) ? 3.0 : 1.5;
                                Etid = tide.evolveOrbitDynamicalTide(_bin, rad1, rad2, poly_type1, poly_type2);
                                change_flag = (Etid>0);
                                // for slowdown case, repeating tide effect based on slowdown factor
                                Float sd_factor_ext = _bin.slowdown.getSlowDownFactor() - 1.5;
                                if (change_flag && sd_factor_ext>0) {
                                    for (Float k=0; k<sd_factor_ext; k=k+1.0) {
                                        Float etid_k = tide.evolveOrbitDynamicalTide(_bin, rad1, rad2, poly_type1, poly_type2);
                                        if (etid_k==0) break;
                                        Etid += etid_k;
                                    }
                                }
                            }

                            if (change_flag) {

                                // record address of modified binary
                                _bin_interrupt.setBinaryTreeAddress(&_bin);

                                // if status not set, set to change
                                if (_bin_interrupt.status == AR::InterruptStatus::none)
                                    _bin_interrupt.status = AR::InterruptStatus::change;
                                _bin.calcParticles(gravitational_constant);
                                p1->pos += _bin.pos;
                                p2->pos += _bin.pos;
                                p1->vel += _bin.vel;
                                p2->vel += _bin.vel;

                                //p1->setBinaryPairID(p2->id);
                                //p2->setBinaryPairID(p1->id);
                                p1->setBinaryInterruptState(BinaryInterruptState::tide);
                                p2->setBinaryInterruptState(BinaryInterruptState::tide);

                                modify_return = 2;

#pragma omp critical
                                {
                                    fout_bse<<"Tide "
                                            <<std::setw(WRITE_WIDTH)<<_bin_interrupt.time_now
                                            <<std::setw(WRITE_WIDTH)<<p1->id
                                            <<std::setw(WRITE_WIDTH)<<p2->id
                                            <<std::setw(WRITE_WIDTH)<<pair_id1
                                            <<std::setw(WRITE_WIDTH)<<pair_id2
                                            <<std::setw(WRITE_WIDTH)<<binary_type_p1
                                            <<std::setw(WRITE_WIDTH)<<binary_type_p2
                                            <<std::setw(WRITE_WIDTH)<<poly_type1
                                            <<std::setw(WRITE_WIDTH)<<poly_type2
                                            <<std::setw(WRITE_WIDTH)<<drdv
                                            <<std::setw(WRITE_WIDTH)<<semi //old
                                            <<std::setw(WRITE_WIDTH)<<ecc  //old
                                            <<std::setw(WRITE_WIDTH)<<Etid
                                            <<std::setw(WRITE_WIDTH)<<Ltid;
                                    _bin.BinarySlowDown::printColumnAscii(fout_bse, WRITE_WIDTH);
                                            p1->star.printColumnAscii(fout_bse, WRITE_WIDTH);
                                            p2->star.printColumnAscii(fout_bse, WRITE_WIDTH);
                                    fout_bse<<std::endl;
                                }

                            }
                        }
                    }
                }
#endif // BSE_BASE
            }
        }
#endif

#ifdef SDAR_PN


#endif
        return modify_return;
    }

#ifndef AR_TTL
    //! (Necessary) calcualte the inverse time transformation factor for drift
    /*! The time transformation factor for drift only depends on (kinetic energy - total energy)
      @param[in] _ekin_minus_etot: ekin - etot
    */
    Float calcGTDriftInv(Float _ekin_minus_etot) {
        return _ekin_minus_etot;
    }
#endif   

    //! (Necessary) calculate the time transformed Hamiltonian
    /*! calculate the time transformed Hamiltonian
      @param[in] _ekin_minus_etot: ekin - etot
    */
    Float calcH(Float _ekin_minus_etot, Float _epot) {
        return log(_ekin_minus_etot) - log(-_epot);
    }

    //! write class data to file with binary format
    /*! @param[in] _fp: FILE type file for output
     */
    void writeBinary(FILE *_fp) const {
        fwrite(&eps_sq, sizeof(Float),1,_fp);
        fwrite(&gravitational_constant, sizeof(Float),1,_fp);
        fwrite(&interrupt_detection_option, sizeof(int),1,_fp);
#ifdef STELLAR_EVOLUTION
        fwrite(&time_interrupt_max, sizeof(Float),1,_fp);
#ifdef BSE_BASE
        fwrite(&stellar_evolution_option, sizeof(int),1,_fp);
        fwrite(&stellar_evolution_write_flag, sizeof(bool),1,_fp);
        fwrite(&tide, sizeof(TwoBodyTide),1,_fp);
#endif
#endif
#ifdef SDAR_PN
        pn.writeBinary(_fp);
#endif
    }

    void printColumnBinary(std::ostream& _fout) const {
        _fout.write(reinterpret_cast<const char*>(&eps_sq), sizeof(Float));
        _fout.write(reinterpret_cast<const char*>(&gravitational_constant), sizeof(Float));
        _fout.write(reinterpret_cast<const char*>(&interrupt_detection_option), sizeof(int));
#ifdef STELLAR_EVOLUTION
        _fout.write(reinterpret_cast<const char*>(&time_interrupt_max), sizeof(Float));
#ifdef BSE_BASE
        _fout.write(reinterpret_cast<const char*>(&stellar_evolution_option), sizeof(int));
        _fout.write(reinterpret_cast<const char*>(&stellar_evolution_write_flag), sizeof(bool));
        _fout.write(reinterpret_cast<const char*>(&tide), sizeof(TwoBodyTide));
#endif
#endif
#ifdef SDAR_PN
    pn.printColumnBinary(_fout);
#endif
    }

    //! read class data to file with binary format
    /*! @param[in] _fp: FILE type file for reading
     */
    void readBinary(FILE *_fin) {
        size_t rcount = fread(&eps_sq, sizeof(Float),1,_fin);
        if (rcount<1) {
            std::cerr<<"Error: ARInteraction:readBinary: get eps_sq fails!\n";
            abort();
        }
        rcount = fread(&gravitational_constant, sizeof(Float),1,_fin);
        if (rcount<1) {
            std::cerr<<"Error: ARInteraction:readBinary: get gravitational_constant fails!\n";
            abort();
        }
        rcount = fread(&interrupt_detection_option, sizeof(int),1,_fin);
        if (rcount<1) {
            std::cerr<<"Error: ARInteraction:readBinary: get interrupt_detection_option fails!\n";
            abort();
        }
#ifdef STELLAR_EVOLUTION
        rcount = fread(&time_interrupt_max, sizeof(Float),1,_fin);
        if (rcount<1) {
            std::cerr<<"Error: ARInteraction:readBinary: get time_interrupt_max fails!\n";
            abort();
        }
#ifdef BSE_BASE
        rcount = fread(&stellar_evolution_option, sizeof(int),1,_fin);
        if (rcount<1) {
            std::cerr<<"Error: ARInteraction:readBinary: get stellar_evolution_option fails!\n";
            abort();
        }
        rcount = fread(&stellar_evolution_write_flag, sizeof(bool),1,_fin);
        if (rcount<1) {
            std::cerr<<"Error: ARInteraction:readBinary: get stellar_evolution_write_flag fails!\n";
            abort();
        }
        rcount = fread(&tide, sizeof(TwoBodyTide),1,_fin);
        if (rcount<1) {
            std::cerr<<"Error: ARInteraction:readBinary: get tide fails!\n";
            abort();
        }
#endif
#endif
#ifdef SDAR_PN
        pn.readBinary(_fin);
#endif
    }    

    void readBinary(std::istream& _fin) {
        _fin.read(reinterpret_cast<char*>(&eps_sq), sizeof(Float));
        if (!_fin) {
            std::cerr<<"Error: ARInteraction:readBinary: get eps_sq fails!\n";
            abort();
        }
        _fin.read(reinterpret_cast<char*>(&gravitational_constant), sizeof(Float));
        if (!_fin) {
            std::cerr<<"Error: ARInteraction:readBinary: get gravitational_constant fails!\n";
            abort();
        }
        _fin.read(reinterpret_cast<char*>(&interrupt_detection_option), sizeof(int));
        if (!_fin) {
            std::cerr<<"Error: ARInteraction:readBinary: get interrupt_detection_option fails!\n";
            abort();
        }
#ifdef STELLAR_EVOLUTION
        _fin.read(reinterpret_cast<char*>(&time_interrupt_max), sizeof(Float));
        if (!_fin) {
            std::cerr<<"Error: ARInteraction:readBinary: get time_interrupt_max fails!\n";
            abort();
        }
#ifdef BSE_BASE
        _fin.read(reinterpret_cast<char*>(&stellar_evolution_option), sizeof(int));
        if (!_fin) {
            std::cerr<<"Error: ARInteraction:readBinary: get stellar_evolution_option fails!\n";
            abort();
        }
        _fin.read(reinterpret_cast<char*>(&stellar_evolution_write_flag), sizeof(bool));
        if (!_fin) {
            std::cerr<<"Error: ARInteraction:readBinary: get stellar_evolution_write_flag fails!\n";
            abort();
        }
        _fin.read(reinterpret_cast<char*>(&tide), sizeof(TwoBodyTide));
        if (!_fin) {
            std::cerr<<"Error: ARInteraction:readBinary: get tide fails!\n";
            abort();
        }
#endif
#endif
#ifdef SDAR_PN
        pn.readBinary(_fin);
#endif
    }
};
