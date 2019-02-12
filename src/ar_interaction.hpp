#pragma once

#include <cmath>
#include "AR/Float.h"
#include "ptclh4.hpp"
#include "AR/force.h"
#include "ar_perturber.hpp"

//! AR interaction clas
class ARInteraction{
public:
    Float eps_sq; ///> softening parameter
    ChangeOver changeover; ///> changover control

    //! (Necessary) calculate acceleration from perturber and the perturbation factor for slowdown calculation
    /*! The Force class acc_pert should be updated
      @param[out] _force: force array to store the calculation results (in acc_pert[3], notice acc_pert may need to reset zero to avoid accummulating old values)
      @param[in] _particles: member particle array
      @param[in] _n_particle: number of member particles
      @param[in] _particle_cm: center-of-mass particle
      @param[in] _perturber: pertuber container
      \return perturbation force to calculate slowdown factor
    */
    Float calcAccAndSlowDownPert(Force* _force, const Ptcl* _particles, const int _n_particle, const PtclH4& _particle_cm, const Perturber& _perturber) {
        for (int i=0; i<_n_particle; i++) {
            Float* acc_pert = _force[i].acc_pert;
            acc_pert[0] = acc_pert[1] = acc_pert[2] = Float(0.0);
        }            
        return 0.0;
    }

    //! (Necessary) calculate inner member acceleration, potential and time transformation function gradient and factor for kick (two-body case)
    /*!
      @param[out] _force: force array to store the calculation results (in acc_in[3] for acceleration and gtgrad[3] for gradient, notice acc/gtgard may need to reset zero to avoid accummulating old values)
      @param[out] _epot: total inner potential energy
      @param[in] _particles: member particle array
      @param[in] _n_particle: number of member particles
      \return the time transformation factor (gt_kick) for kick step
    */
    Float calcAccPotAndGTKickTwo(Force* _force, Float& _epot, const Particle* _particles, const int _n_particle) {
        assert(_n_particle==2);

        // acceleration
        const Float mass1 = _particles[0].mass;
        const Float* pos1 = _particles[0].pos;

        const Float mass2 = _particles[1].mass;
        const Float* pos2 = _particles[1].pos;

        Float m1m2 = mass1*mass2;
        
        Float dr[3] = {pos2[0] -pos1[0],
                       pos2[1] -pos1[1],
                       pos2[2] -pos1[2]};
        Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
        Float inv_r = 1.0/sqrt(r2);
        Float inv_r3 = inv_r*inv_r*inv_r;

        Float* acc1 = _force[0].acc_in;
        Float* acc2 = _force[1].acc_in;

#ifdef AR_CHANGEOVER
        const Float kpot  = pars->changeover.calcPotW(rij);
        const Float k     = pars->changeover.calcAcc0W(rij);
        Float mor3_1 = mass2*inv_r3*k;
        Float mor3_2 = mass1*inv_r3*k;

        Float m1m2or = m1m2*inv_r*kpot;
#else
        Float mor3_1 = mass2*inv_r3;
        Float mor3_2 = mass1*inv_r3;

        Float m1m2or = m1m2*inv_r;
#endif

        acc1[0] = mor3_1 * dr[0];
        acc1[1] = mor3_1 * dr[1];
        acc1[2] = mor3_1 * dr[2];

        acc2[0] = - mor3_2 * dr[0];
        acc2[1] = - mor3_2 * dr[1];
        acc2[2] = - mor3_2 * dr[2];

#ifdef AR_TTL 
        // trans formation function gradient
#ifdef AR_CHANGEOVER
        Float m1m2or3 = m1m2*inv_r3*k;
#else
        Float m1m2or3 = m1m2*inv_r3;
#endif
        Float* gtgrad1 = _force[0].gtgrad;
        Float* gtgrad2 = _force[1].gtgrad;

        gtgrad1[0] = m1m2or3 * dr[0];
        gtgrad1[1] = m1m2or3 * dr[1];
        gtgrad1[2] = m1m2or3 * dr[2];

        gtgrad2[0] = - gtgrad1[0];
        gtgrad2[1] = - gtgrad1[1];
        gtgrad2[2] = - gtgrad1[2];
#endif

        // potential energy
        _epot = - m1m2or;

        // transformation factor for kick
        Float gt_kick = 1.0/m1m2or;

        return gt_kick;
    }

    //! (Necessary) calculate inner member acceleration, potential and time transformation function gradient and factor for kick
    /*!
      @param[out] _force: force array to store the calculation results (in acc_in[3] for acceleration and gtgrad[3] for gradient, notice acc/gtgard may need to reset zero to avoid accummulating old values)
      @param[out] _epot: total inner potential energy
      @param[in] _particles: member particle array
      @param[in] _n_particle: number of member particles
      \return the time transformation factor (gt_kick) for kick step
    */
    Float calcAccPotAndGTKick(Force* _force, Float& _epot, const Particle* _particles, const int _n_particle) {
        _epot = Float(0.0);
        Float gt_kick = Float(0.0);

#ifdef AR_CHANGEOVER
        const Float kpot  = pars->changeover.calcPotW(rij);
        const Float k     = pars->changeover.calcAcc0W(rij);
#endif

        for (int i=0; i<_n_particle; i++) {
            const Float massi = _particles[i].mass;
            const Float* posi = _particles[i].pos;
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
                const Float* posj = _particles[j].pos; 
                Float dr[3] = {posj[0] -posi[0],
                               posj[1] -posi[1],
                               posj[2] -posi[2]};
                Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
                Float inv_r = 1.0/sqrt(r2);
                Float inv_r3 = inv_r*inv_r*inv_r;
#ifdef AR_CHANGEOVER
                Float mor3 = massj*inv_r3*k;
                Float mor = massj*inv_r*kpot;
#else
                Float mor3 = massj*inv_r3;
                Float mor = massj*inv_r;
#endif
                acci[0] += mor3 * dr[0];
                acci[1] += mor3 * dr[1];
                acci[2] += mor3 * dr[2];

#ifdef AR_TTL                     
                Float mimjor3 = massi*mor3;
                gtgradi[0] += mimjor3 * dr[0];
                gtgradi[1] += mimjor3 * dr[1];
                gtgradi[2] += mimjor3 * dr[2];
#endif

                poti -= mor;
                gtki += mor;
                    
            }
            _epot += poti * massi;
            gt_kick += gtki * massi;
        }
        _epot   *= 0.5;
        gt_kick = 2.0/gt_kick;

        return gt_kick;
    }

#ifndef AR_TTL
    //! (Necessary) calcualte the time transformation factor for drift
    /*! The time transformation factor for drift only depends on (kinetic energy - total energy)
      @param[in] _ekin_minus_etot: ekin - etot
    */
    Float calcGTDrift(Float _ekin_minus_etot) {
        return 1.0/_ekin_minus_etot;
    }
#endif   
};

//! Newtonian acceleration from particle p to particle i (function type of ::ARC::pair_Ap)
/*! 
  @param[out] Ai: acceleration vector. \f$Aij[1:3] = m_i m_p (xp[1:3]-xi[1:3]) / |xp-xi|^3 \f$.
  @param[in]  time: next time
  @param[in]  p:  chain particle list
  @param[in]  np: chain particle number
  @param[in]  pert: perturber address list (first should be center-of-mass of AR)
  @param[in]  pf:  perturber force array
  @param[in]  pars: ARC pars including rin, rout and perturbation kepler spline interpolation class
*/
template<class Tptcl, class Tpert, class Tforce, class extpar>
PS::S32 Newtonian_extA_pert (double3* acc, const PS::F64 time, Tptcl* p, const PS::S32 np, Tpert* pert, Tforce* pf, const PS::S32 npert, extpar* pars){
#ifdef HARD_DEBUG
    if(npert<1) {
        std::cerr<<"Error: perturber number not enough, at least c.m. is needed!"<<std::endl;
        abort();
    }
#endif
    static const PS::F64 inv3 = 1.0 / 3.0;
    PS::F64 xp[npert][3];

    for(int i=0; i<npert; i++) {
        PS::F64 dt = time - pert[i]->time;
        xp[i][0] = pert[i]->pos.x + dt*(pert[i]->vel.x + 0.5*dt*(pert[i]->acc0.x + inv3*dt*pert[i]->acc1.x));
        xp[i][1] = pert[i]->pos.y + dt*(pert[i]->vel.y + 0.5*dt*(pert[i]->acc0.y + inv3*dt*pert[i]->acc1.y));
        xp[i][2] = pert[i]->pos.z + dt*(pert[i]->vel.z + 0.5*dt*(pert[i]->acc0.z + inv3*dt*pert[i]->acc1.z));

        //xp[i][0] = pert[i]->pos.x + dt*(pert[i]->vel.x + 0.5*dt*(pert[i]->acc0.x + inv3*dt*(pert[i]->acc1.x + 0.25*dt*(pert[i]->acc2.x + 0.2*dt*pert[i]->acc3.x))));
        //xp[i][1] = pert[i]->pos.y + dt*(pert[i]->vel.y + 0.5*dt*(pert[i]->acc0.y + inv3*dt*(pert[i]->acc1.y + 0.25*dt*(pert[i]->acc2.x + 0.2*dt*pert[i]->acc3.x))));
        //xp[i][2] = pert[i]->pos.z + dt*(pert[i]->vel.z + 0.5*dt*(pert[i]->acc0.z + inv3*dt*(pert[i]->acc1.z + 0.25*dt*(pert[i]->acc2.x + 0.2*dt*pert[i]->acc3.x))));
        //xp[i] = pert[i]->pos + dt*
        //    (pert[i]->vel* + 0.5*dt*(
        //        pert[i]->acc0 + inv3*dt*(
        //            pert[i]->acc1 + 0.25*dt*(
        //                pert[i]->acc2 + 0.2*dt*pert[i]->acc3))));
        //xp[i] = pert[i]->pos;
    }

#ifdef HARD_DEBUG
    PS::F64 mt = 0.0;
    for(int i=0; i<np; i++) mt += p[i].mass;
#ifdef HARD_DEBUG_DUMP
    if (abs(mt-pert[0]->mass)>=1e-10) return 6;
#else
    assert(abs(mt-pert[0]->mass)<1e-10);
#endif
#endif

    for(int i=0; i<np; i++) {
        PS::F64 xi[3];
        xi[0] = p[i].pos.x + xp[0][0];
        xi[1] = p[i].pos.y + xp[0][1];
        xi[2] = p[i].pos.z + xp[0][2];

        acc[i][0] = -pf[0]->acc0.x; 
        acc[i][1] = -pf[0]->acc0.y;        
        acc[i][2] = -pf[0]->acc0.z; 
//            acc[i][0] = acc[i][1] = acc[i][2] = 0.0;
        for(int j=1; j<npert; j++) {

#ifdef HARD_DEBUG
#ifdef HARD_DEBUG_DUMP
            if (p[i].id==pert[j]->id) return 7;
#else
            assert(p[i].id!=pert[j]->id);
#endif
#endif            
            PS::F64 dx[3];
            dx[0] = xp[j][0] - xi[0];
            dx[1] = xp[j][1] - xi[1];
            dx[2] = xp[j][2] - xi[2];
            
            //std::cerr<<"i = "<<i<<" j = "<<j<<" dx = "<<dx<<std::endl;
            //PS::F64 mi = p[i].mass;
            PS::F64 mp = pert[j]->mass;
            PS::F64 dr2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2] + pars->eps2;
            PS::F64 dr  = std::sqrt(dr2);
            PS::F64 dr3 = dr*dr2;
            PS::F64 mor3 = mp/dr3;

#ifdef HARD_DEBUG
            if (std::isnan(mor3)) {
                std::cerr<<"Error: Nan detected in ARC pert force, xi["<<i<<"]="<<xi<<" xj["<<j<<"]="<<xp[j]<<" pi.id="<<p[i].id<<" pj.id="<<pert[j]->id<<" cm.id="<<pert[0]->id<<" eps2="<<pars->eps2<<std::endl;
                abort();
            }
#endif
            // smpars[2:3]: rcut_out, rcut_in
            //  const double k   = cutoff_poly_3rd(dr, smpars[0], smpars[1]);
            //  const double kdx = cutoff_poly_3rd_dr(dr, dx, smpars[0], smpars[1]);
            //  const double kdy = cutoff_poly_3rd_dr(dr, dy, smpars[0], smpars[1]);
            //  const double kdz = cutoff_poly_3rd_dr(dr, dz, smpars[0], smpars[1]);  
            // const PS::F64 r_out = pars->rout;
            //  const PS::F64 r_out = std::max(pi.r_out, pp.r_out);
            const PS::F64 r_in  = pars->rin;
            const PS::F64 r_oi_inv = pars->r_oi_inv;
            const PS::F64 r_A   = pars->r_A;
            //const PS::F64 k     = CalcW(dr/r_out, r_in/r_out);
            const PS::F64 k  = cutoff_poly_3rd(dr, r_oi_inv, r_A, r_in);

            //Pij = - mi*mp / dr * (1-k);

            // Aij[0] = mp * dx / dr3 * (1-k) + Pij * (1-kdx);
            // Aij[1] = mp * dy / dr3 * (1-k) + Pij * (1-kdy);
            // Aij[2] = mp * dz / dr3 * (1-k) + Pij * (1-kdz);
            acc[i][0] += mor3 * dx[0] * k;
            acc[i][1] += mor3 * dx[1] * k;
            acc[i][2] += mor3 * dx[2] * k;
        }

#ifdef SOFT_PERT
        // soft perturbation
#ifdef TIDAL_TENSOR
        pars->eval(acc[i], p[i].pos);
#else
        if(p[i].status==0) pars->eval(acc[i], time);
#endif
#endif
    }

    return 0;
}
/// end Newtonian cut force (L.Wang)

template<class Tptcl, class Tpert, class Tforce, class extpar>
PS::S32 Newtonian_extA_soft (double3* acc, const PS::F64 time, Tptcl* p, const PS::S32 np, Tpert* pert, Tforce* pf, const PS::S32 npert, extpar* pars){
#ifdef HARD_DEBUG
    if(npert>1) {
        std::cerr<<"Error: perturboer number should be zero !"<<std::endl;
        abort();
    }
#endif
    for(int i=0; i<np; i++) {
        acc[i][0] = acc[i][1] = acc[i][2] = 0.0;
#ifdef SOFT_PERT
#ifdef TIDAL_TENSOR
        pars->eval(acc[i], p[i].pos);
#else
        if(p[i].status==0) pars->eval(acc[i], time);
#endif
#endif
    }

    return 0;
}

