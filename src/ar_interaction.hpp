#pragma once
#include <cmath>
#include "Common/Float.h"
#include "Common/binary_tree.h"
#include "changeover.hpp"
#include "AR/force.h"
#include "hard_ptcl.hpp"
#include "Hermite/hermite_particle.h"
#include "ar_perturber.hpp"
#include "two_body_tide.hpp"
#ifdef BSE_BASE
#include "bse_interface.h"
#endif

//! AR interaction clas
class ARInteraction{
public:
    typedef H4::ParticleH4<PtclHard> H4Ptcl;
    Float eps_sq; ///> softening parameter
    Float gravitational_constant;
#ifdef STELLAR_EVOLUTION
    int stellar_evolution_option;
    bool stellar_evolution_write_flag;
#ifdef BSE_BASE
    BSEManager bse_manager;
    TwoBodyTide tide;
    std::ofstream fout_sse; ///> log file for SSE event
    std::ofstream fout_bse; ///> log file for BSE event

    ARInteraction(): eps_sq(Float(-1.0)), gravitational_constant(Float(-1.0)), stellar_evolution_option(1), stellar_evolution_write_flag(true), bse_manager(), fout_sse(), fout_bse() {}
#else
    ARInteraction(): eps_sq(Float(-1.0)), gravitational_constant(Float(-1.0)), stellar_evolution_option(0), stellar_evolution_write_flag(true) {}
#endif
#else
    ARInteraction(): eps_sq(Float(-1.0)), gravitational_constant(Float(-1.0)) {}
#endif

    //! (Necessary) check whether publicly initialized parameters are correctly set
    /*! \return true: all parmeters are correct. In this case no parameters, return true;
     */
    bool checkParams() {
        ASSERT(eps_sq>=0.0);
        ASSERT(gravitational_constant>0.0);
#ifdef BSE_BASE
        ASSERT(stellar_evolution_option==0 || (stellar_evolution_option==1 && bse_manager.checkParams()) || (stellar_evolution_option==2 && bse_manager.checkParams() && tide.checkParams()));
        ASSERT(!stellar_evolution_write_flag||(stellar_evolution_write_flag&&fout_sse.is_open()));
        ASSERT(!stellar_evolution_write_flag||(stellar_evolution_write_flag&&fout_bse.is_open()));
#endif
        return true;
    }        

    //! print parameters
    void print(std::ostream & _fout) const{
        _fout<<"eps_sq : "<<eps_sq<<std::endl
             <<"G      : "<<gravitational_constant<<std::endl;
#ifdef STELLAR_EVOLUTION
        _fout<<"SE_opt : "<<stellar_evolution_option<<std::endl;
#endif
    }    

    //! (Necessary) calculate inner member acceleration, potential and inverse time transformation function gradient and factor for kick (two-body case)
    /*!
      @param[out] _f1: force for particle 1 to store the calculation results (in acc_in[3] for acceleration and gtgrad[3] for gradient, notice acc/gtgard are overwritten, not accummulating old values)
      @param[out] _f2: force for particle 2
      @param[out] _epot: total inner potential energy
      @param[in] _p1: particle 1
      @param[in] _p2: particle 2
      \return the inverse time transformation factor (gt_kick_inv) for kick step
    */
    inline Float calcInnerAccPotAndGTKickInvTwo(AR::Force& _f1, AR::Force& _f2, Float& _epot, const PtclHard& _p1, const PtclHard& _p2) {
        // acceleration
        const Float mass1 = _p1.mass;
        const Float* pos1 = &_p1.pos.x;

        const Float mass2 = _p2.mass;
        const Float* pos2 = &_p2.pos.x;

        Float gm1 = gravitational_constant*mass1;
        Float gm2 = gravitational_constant*mass2;
        Float gm1m2 = gm1*mass2;
        
        Float dr[3] = {pos2[0] -pos1[0],
                       pos2[1] -pos1[1],
                       pos2[2] -pos1[2]};
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
            ChangeOver* changeover[n_pert_single];
            H4::NBAdr<PtclHard>::Group* ptclgroup[n_pert_group];

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
                const auto& pi = _particles[i];
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

                acc_pert_cm[0] += pi.mass *acc_pert[0];
                acc_pert_cm[1] += pi.mass *acc_pert[1];
                acc_pert_cm[2] += pi.mass *acc_pert[2];

                mcm += pi.mass;

            }
//#ifdef AR_DEBUG
//            ASSERT(abs(mcm-_particle_cm.mass)<1e-10);
//#endif
                
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
#ifdef SOFT_PERT
            if(_perturber.soft_pert!=NULL) {
                for(int i=0; i<_n_particle; i++) {
                    Float* acc_pert = _force[i].acc_pert;
                    Float& pot_pert = _force[i].pot_pert;
                    const auto& pi = _particles[i];
                    acc_pert[0] = acc_pert[1] = acc_pert[2] = pot_pert = Float(0.0);
                    // avoid too large perturbation force if system is disruptted
                    if (pi.pos*pi.pos<pi.changeover.getRout()*pi.changeover.getRout()) {
                        _perturber.soft_pert->eval(acc_pert, pi.pos);
                        pot_pert += _perturber.soft_pert->evalPot(pi.pos);
                    }
                }
            }
#endif
        }
    }
    
    //! (Necessary) calculate acceleration from internal members and perturbers
    /*! The Force class acc_pert should be updated
      @param[out] _force: force array to store the calculation results (in acc_pert[3], notice acc_pert may need to reset zero to avoid accummulating old values)
      @param[out] _epot: potential 
      @param[in] _particles: member particle array
      @param[in] _n_particle: number of member particles
      @param[in] _particle_cm: center-of-mass particle
      @param[in] _perturber: pertuber container
      @param[in] _time: current time
      \return perturbation energy to calculate slowdown factor
    */
    Float calcAccPotAndGTKickInv(AR::Force* _force, Float& _epot, const PtclHard* _particles, const int _n_particle, const H4Ptcl& _particle_cm, const ARPerturber& _perturber, const Float _time) {
        // inner force
        Float gt_kick_inv;
        if (_n_particle==2) gt_kick_inv = calcInnerAccPotAndGTKickInvTwo(_force[0], _force[1], _epot, _particles[0], _particles[1]);
        else gt_kick_inv = calcInnerAccPotAndGTKickInv(_force, _epot, _particles, _n_particle);

        calcAccPert(_force, _particles, _n_particle, _particle_cm, _perturber, _time);

        return gt_kick_inv;
    }


    //! calculate perturbation from c.m. acceleration
    Float calcPertFromForce(const Float* _force, const Float _mp, const Float _mpert) {
        Float force2 = _force[0]*_force[0]+_force[1]*_force[1]+_force[2]*_force[2];
#ifdef AR_SLOWDOWN_PERT_R4
        return force2/(gravitational_constant*_mp*_mpert);
#else
        Float force = sqrt(force2)/gravitational_constant;
        return sqrt(force/(_mp*_mpert))*force;
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

    //! (Necessary) calculate slowdown perturbation and timescale
    /*!
      @param[out] _pert_out: perturbation 
      @param[out] _t_min_sq: timescale limit 
      @param[in] _time: physical time for prediction
      @param[in] _particle_cm: center-of-mass particle
      @param[in] _perturber: pertuber container
    */
    void calcSlowDownPert(Float& _pert_out, Float& _t_min_sq, const Float& _time, const H4Ptcl& _particle_cm, const ARPerturber& _perturber) {
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
            _p.time_interrupt = _p.time_record + bse_manager.getTimeStepStar(_p.star);

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
                    star_bk.printColumn(fout_sse, WRITE_WIDTH);
                    _p.star.printColumn(fout_sse, WRITE_WIDTH);
                    //output.printColumn(fout_sse, WRITE_WIDTH);
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
                        _p.star.printColumn(fout_sse, WRITE_WIDTH);
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
#ifdef BSE_BASE
                // if member is star, evolve single star using SSE
                else if (stellar_evolution_option>0) {
                    ASSERT(bse_manager.checkParams());
                    modify_branch[k] = modifyOneParticle(*_bin.getMember(k), _bin.getMember(k)->time_record, _bin_interrupt.time_now);
                    modify_return = std::max(modify_return, modify_branch[k]);
                    // if status not set, set to change
                    if (modify_branch[k]>0&&_bin_interrupt.status == AR::InterruptStatus::none) {
                        _bin_interrupt.status = AR::InterruptStatus::change;
                        _bin_interrupt.adr = &_bin;
                    }
                }
#endif
            }
            // ensure to record the root binary tree to include all changed members
            if (modify_branch[0]>0&&modify_branch[1]>0) {
                _bin_interrupt.adr = &_bin;
            }
            if (_bin_interrupt.status == AR::InterruptStatus::destroy) {
                // if both branch has destroyed, set destroy status, otherwise set merge status
                if (!(modify_branch[0]==3&&modify_branch[1]==3))
                    _bin_interrupt.status = AR::InterruptStatus::merge;
            }
        }
        else {
            auto* p1 = _bin.getLeftMember();
            auto* p2 = _bin.getRightMember();

            COMM::Vector3<Float> pos_red(p2->pos[0] - p1->pos[0], p2->pos[1] - p1->pos[1], p2->pos[2] - p1->pos[2]);
            COMM::Vector3<Float> vel_red(p2->vel[0] - p1->vel[0], p2->vel[1] - p1->vel[1], p2->vel[2] - p1->vel[2]);
            Float drdv = pos_red * vel_red;

#ifdef BSE_BASE
            auto postProcess =[&](StarParameterOut* out, Float* pos_cm, Float*vel_cm, Float& semi, Float& ecc, int binary_type_final) {
                // if status not set, set to change
                if (_bin_interrupt.status == AR::InterruptStatus::none) 
                    _bin_interrupt.status = AR::InterruptStatus::change;
                    
                // set return flag >0
                modify_return = 2;

                p1->time_record = _bin_interrupt.time_now - bse_manager.getDTMiss(out[0]);
                p2->time_record = _bin_interrupt.time_now - bse_manager.getDTMiss(out[1]);

                // estimate next time to check 
                p1->time_interrupt = p1->time_record + bse_manager.getTimeStepBinary(p1->star, p2->star, semi, ecc, binary_type_final);
                p2->time_interrupt = p1->time_interrupt;

                // reset collision state since binary orbit changes
                if (p1->getBinaryInterruptState()== BinaryInterruptState::collision) 
                    p1->setBinaryInterruptState(BinaryInterruptState::none);
                if (p2->getBinaryInterruptState()== BinaryInterruptState::collision)
                    p2->setBinaryInterruptState(BinaryInterruptState::none);

                // set binary status
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

                    // case when SN kick appears
                    bool kick_flag=false;
                    for (int k=0; k<2; k++) {
                        double dv[4];
                        auto* pk = _bin.getMember(k);
                        dv[3] = bse_manager.getVelocityChange(dv,out[k]);
                        if (dv[3]>0) {
                            kick_flag=true;
#pragma omp critical 
                            {
                                fout_bse<<"SN_kick "
                                        <<std::setw(WRITE_WIDTH)<<p1->id
                                        <<std::setw(WRITE_WIDTH)<<p2->id
                                        <<std::setw(WRITE_WIDTH)<<k+1
                                        <<std::setw(WRITE_WIDTH)<<dv[3]*bse_manager.vscale;
                                pk->star.printColumn(fout_bse, WRITE_WIDTH);
                                fout_bse<<std::endl;
                            }
                            for (int k=0; k<3; k++) pk->vel[k] += dv[k];
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
                    DATADUMP("dump_merger");
                }

            };

            bool check_flag = false;
            if (stellar_evolution_option>0) {
                int binary_type_p1 = static_cast<int>(p1->getBinaryInterruptState());
                int binary_type_p2 = static_cast<int>(p2->getBinaryInterruptState());
                int binary_type_init = 0;
                if (binary_type_p1==binary_type_p2) binary_type_init = binary_type_p1;

                // check whether need to update
                double time_check = std::min(p1->time_interrupt, p2->time_interrupt);
                Float dt = _bin_interrupt.time_now - std::max(p1->time_record,p2->time_record);
                if (time_check<=_bin_interrupt.time_now&&dt>0) check_flag = true;

                if (!check_flag) {
                    // pre simple check whether calling BSE is needed
                    Float t_record_min = std::min(p1->time_record,p2->time_record);
                    

                    if (t_record_min<_bin_interrupt.time_now&&_bin.semi>0 && (!bse_manager.isMassTransfer(binary_type_init)) && (!bse_manager.isDisrupt(binary_type_init))) {
                        check_flag=bse_manager.isCallBSENeeded(p1->star, p2->star, _bin.semi, _bin.ecc, dt);
                    }
                }

                if (check_flag) {
                    ASSERT(bse_manager.checkParams());
                    // record address of modified binary
                    _bin_interrupt.adr = &_bin;

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
                    int event_flag = bse_manager.evolveBinary(p1->star, p2->star, out[0], out[1], semi, period, ecc, bin_event, binary_type_init, dt);

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
                        bse_manager.dumpRandConstant("bse.rand.par");
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
                            else if (bse_manager.isMassTransfer(binary_type)) event_flag = std::max(event_flag, 2); // orbit change
                            else if (bse_manager.isDisrupt(binary_type)) event_flag = std::max(event_flag, 3); // disrupt
                            else if (bse_manager.isMerger(binary_type)) event_flag = std::max(event_flag, 4); // Merger
                            binary_type_final = binary_type;
                        }
                        else if(binary_type<0) break;
                    }

                    // update semi
                    mtot = bse_manager.getMass(p1->star) + bse_manager.getMass(p2->star);
                    semi = COMM::Binary::periodToSemi(period, mtot, gravitational_constant);

                    postProcess(out, pos_cm, vel_cm, semi, ecc, binary_type_final);
                }
            }
#endif // BSE_BASE

            // dynamical merger and tide check
            if (_bin_interrupt.status!=AR::InterruptStatus::merge&&_bin_interrupt.status!=AR::InterruptStatus::destroy) {

                auto merge = [&](const Float& dr, const Float& t_peri, const Float& sd_factor) {
                    _bin_interrupt.adr = &_bin;
                
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

                        postProcess(out, pos_cm, vel_cm, semi, ecc, 0);
                        if (stellar_evolution_write_flag&&(p1->mass==0.0||p2->mass==0.0)) {
#pragma omp critical
                            {
                                fout_bse<<"Dynamic_merge: "
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
                                p1_star_bk.printColumn(fout_bse, WRITE_WIDTH);
                                p2_star_bk.printColumn(fout_bse, WRITE_WIDTH);
                                // after
                                p1->star.printColumn(fout_bse, WRITE_WIDTH);
                                p2->star.printColumn(fout_bse, WRITE_WIDTH);
                                fout_bse<<std::endl;

                                DATADUMP("dump_merger");
                                
                            }
                        }
                    }
#else //not BSE_BASE
                    // print data
                    std::cerr<<"Binary Merge: time: "<<_bin_interrupt.time_now<<std::endl;
                    _bin.Binary::printColumnTitle(std::cerr);
                    //PtclHard::printColumnTitle(std::cerr);
                    //PtclHard::printColumnTitle(std::cerr);
                    std::cerr<<std::endl;
                    _bin.Binary::printColumn(std::cerr);
                    //p1->printColumn(std::cerr);
                    //p2->printColumn(std::cerr);
                    std::cerr<<std::endl;

                    // set return flag >0
                    modify_return = 2;

                    p1->time_record = _bin_interrupt.time_now;
                    p2->time_record = _bin_interrupt.time_now;
            
                    // new particle data
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

                    if (_bin_interrupt.status == AR::InterruptStatus::none) 
                        _bin_interrupt.status = AR::InterruptStatus::merge;

                    // reset collision state since binary orbit changes
                    p1->setBinaryInterruptState(BinaryInterruptState::none);
                    p2->setBinaryInterruptState(BinaryInterruptState::none);

                    p2->group_data.artificial.setParticleTypeToUnused(); // necessary to identify particle to remove
#endif
                    //p1->setBinaryPairID(0);
                    //p2->setBinaryPairID(0);
                };
                
                // delayed merger
                if (p1->getBinaryInterruptState()== BinaryInterruptState::collision && 
                    p2->getBinaryInterruptState()== BinaryInterruptState::collision &&
                    (p1->time_interrupt<_bin_interrupt.time_now && p2->time_interrupt<_bin_interrupt.time_now) &&
                    (p1->getBinaryPairID()==p2->id||p2->getBinaryPairID()==p1->id)) {
                    Float dr[3] = {p1->pos[0] - p2->pos[0], 
                                   p1->pos[1] - p2->pos[1], 
                                   p1->pos[2] - p2->pos[2]};
                    Float dr2  = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
                    merge(std::sqrt(dr2), 0.0, 1.0);
                }
                else {
                    // check merger
                    Float radius = p1->radius + p2->radius;
#ifndef BSE_BASE
                    // slowdown case
                    if (_bin.slowdown.getSlowDownFactor()>1.0) {
                        ASSERT(_bin.semi>0.0);
                        Float drdv;
                        _bin.particleToSemiEcc(_bin.semi, _bin.ecc, _bin.r, drdv, *_bin.getLeftMember(), *_bin.getRightMember(), gravitational_constant);
                        Float peri = _bin.semi*(1 - _bin.ecc);
                        if (peri<radius && p1->getBinaryPairID()!=p2->id&&p2->getBinaryPairID()!=p1->id) {
                            Float ecc_anomaly  = _bin.calcEccAnomaly(_bin.r);
                            Float mean_anomaly = _bin.calcMeanAnomaly(ecc_anomaly, _bin.ecc);
                            Float mean_motion  = sqrt(gravitational_constant*_bin.mass/(fabs(_bin.semi*_bin.semi*_bin.semi))); 
                            Float t_peri = mean_anomaly/mean_motion;
                            if (drdv<0 && t_peri<_bin_interrupt.time_end-_bin_interrupt.time_now) {
                                Float dr[3] = {p1->pos[0] - p2->pos[0], 
                                               p1->pos[1] - p2->pos[1], 
                                               p1->pos[2] - p2->pos[2]};
                                Float dr2  = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
                                merge(std::sqrt(dr2), t_peri, _bin.slowdown.getSlowDownFactor());
                            }
                            else if (_bin.semi>0||(_bin.semi<0&&drdv<0)) {
                                p1->setBinaryPairID(p2->id);
                                p2->setBinaryPairID(p1->id);
                                p1->setBinaryInterruptState(BinaryInterruptState::collision);
                                p2->setBinaryInterruptState(BinaryInterruptState::collision);
                                p1->time_interrupt = _bin_interrupt.time_now + drdv<0 ? t_peri : (_bin.period - t_peri);
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
#else
                    // in bse case, handle binary merger in bse, only check hyperbolic merger
                    if (_bin.semi<0.0) {
                        Float dr[3] = {p1->pos[0] - p2->pos[0], 
                                       p1->pos[1] - p2->pos[1], 
                                       p1->pos[2] - p2->pos[2]};
                        Float dr2  = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
                        if (dr2<radius*radius) merge(std::sqrt(dr2), 0.0, 1.0);
                    }
#endif
                }


#ifdef BSE_BASE
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
                                 || bse_manager.isDisrupt(binary_type_p1)
                                 || binary_type_p1 == 14)
                            tide_flag = false;

                        bool change_flag=false;
                        if (tide_flag) {
                            Float poly_type1=0, poly_type2=0;
                            Float Etid=0, Ltid=0;
                            Float semi = _bin.semi;
                            Float ecc  = _bin.ecc;
                            Float rad1 = bse_manager.getStellarRadius(p1->star);
                            Float rad2 = bse_manager.getStellarRadius(p2->star);
                            if (p1->star.kw>=10&&p1->star.kw<15&&p2->star.kw>=10&&p2->star.kw<15) {
                                if (_bin.semi<0) {
                                    bool merge_flag = tide.evolveOrbitHyperbolicGW(_bin, Etid, Ltid);
                                    if (merge_flag) {
                                        Float dr[3] = {p1->pos[0] - p2->pos[0], 
                                                       p1->pos[1] - p2->pos[1], 
                                                       p1->pos[2] - p2->pos[2]};
                                        Float dr2  = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
                                        merge(std::sqrt(dr2), 0.0, 1.0);
                                    }
                                    else change_flag = true;
                                }
                            }
                            else if (std::min(p1->star.kw, p2->star.kw)<13) {
                                poly_type1 = (p1->star.kw<=2) ? 3.0 : 1.5;
                                poly_type2 = (p2->star.kw<=2) ? 3.0 : 1.5;
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

                                _bin_interrupt.adr = &_bin;

                                // if status not set, set to change
                                if (_bin_interrupt.status == AR::InterruptStatus::none) 
                                    _bin_interrupt.status = AR::InterruptStatus::change;
                                _bin.calcParticles(gravitational_constant);
                                p1->pos += _bin.pos;
                                p2->pos += _bin.pos;
                                p1->vel += _bin.vel;
                                p2->vel += _bin.vel;

                                p1->setBinaryPairID(p2->id);
                                p2->setBinaryPairID(p1->id);
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
                                    _bin.BinarySlowDown::printColumn(fout_bse, WRITE_WIDTH);
                                    p1->star.printColumn(fout_bse, WRITE_WIDTH);
                                    p2->star.printColumn(fout_bse, WRITE_WIDTH);
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

    //! (Necessary) calculate the time transformed Hamiltonian
    /*! calculate the time transformed Hamiltonian
      @param[in] _ekin_minus_etot: ekin - etot
    */
    Float calcH(Float _ekin_minus_etot, Float _epot) {
        return log(_ekin_minus_etot) - log(-_epot);
    }
#endif   

    //! write class data to file with binary format
    /*! @param[in] _fp: FILE type file for output
     */
    void writeBinary(FILE *_fp) const {
        fwrite(&eps_sq, sizeof(Float),1,_fp);
        fwrite(&gravitational_constant, sizeof(Float),1,_fp);
#ifdef STELLAR_EVOLUTION
        fwrite(&stellar_evolution_option, sizeof(int),1,_fp);
        fwrite(&stellar_evolution_write_flag, sizeof(bool),1,_fp);
#endif
    }

    //! read class data to file with binary format
    /*! @param[in] _fp: FILE type file for reading
     */
    void readBinary(FILE *_fin) {
        size_t rcount = fread(&eps_sq, sizeof(Float),1,_fin);
        rcount += fread(&gravitational_constant, sizeof(Float),1,_fin);
        if (rcount<2) {
            std::cerr<<"Error: Data reading fails! requiring data number is 2, only obtain "<<rcount<<".\n";
            abort();
        }
#ifdef STELLAR_EVOLUTION
        rcount += fread(&stellar_evolution_option, sizeof(int),1,_fin);
        rcount += fread(&stellar_evolution_write_flag, sizeof(bool),1,_fin);
        if (rcount<4) {
            std::cerr<<"Error: Data reading fails! requiring data number is 4, only obtain "<<rcount<<".\n";
            abort();
        }
#endif
    }    
};

