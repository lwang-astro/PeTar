#pragma once
#include "transit.h"

#ifdef HARD_DEBUG_DEEP_CHECK
#include "AR.h"
#endif

typedef double double3[3];

inline void CalcAccPotShortWithLinearCutoff(const PS::F64vec & posi,
					    PS::F64vec & acci_pla,
					    PS::F64 & poti_tot,
					    const PS::F64vec & posj,
					    const PS::F64 massj,
                        const PS::F64 mass_bkj,
                        const PS::S32 pot_control_flag,
					    const PS::F64 eps2,
					    const PS::F64 rcut_oi_inv,
                        const PS::F64 rcut_A,
                        const PS::F64 rcut_out,
					    const PS::F64 rcut_in){
    const PS::F64vec rij = posi - posj;
    const PS::F64 r2 = rij * rij;
    const PS::F64 r2_eps = r2 + eps2;
    const PS::F64 rcut2_out = rcut_out * rcut_out;
    const PS::F64 R = 1.0/sqrt(r2_eps);
    const PS::F64 Rm = massj * R;
    const PS::F64 R2 = R * R;
    const PS::F64 Rm3 = Rm * R2;
    const PS::F64 r_eps = R * r2_eps;
    const PS::F64 k = cutoff_poly_3rd(r_eps, rcut_oi_inv, rcut_A, rcut_in);
    const PS::F64 r2_max = (r2_eps > rcut2_out) ? r2_eps : rcut2_out;
    const PS::F64 R_max = 1.0/sqrt(r2_max);
    const PS::F64 Rm_max = massj * R_max;
    const PS::F64 R2_max = R_max * R_max;
    const PS::F64 Rm3_max = Rm_max * R2_max;
    if(pot_control_flag==0) poti_tot -= (Rm - Rm_max);   // single, remove cutoff, obtain total potential
    else if(pot_control_flag==1) poti_tot -= (mass_bkj * R - Rm_max); // member mass is zero, use backup value
    else poti_tot += Rm_max; // artifical, should be exclude for orbital artifical particles, since it is inside neighbor, Rm_max cancel it out, the contribution become zero
    acci_pla -= (Rm3*(1-k) - Rm3_max)*rij;
}

// Pure Newtonian pair force function
template <class Tptcl, class extpar>
int Newtonian_AW (double Aij[3], double &Pij, double pWij[3], double &Wij, const double xij[3], const Tptcl &pi, const Tptcl &pj, extpar* pars) {
    const double invrij = 1.0/std::sqrt(xij[0]*xij[0]+xij[1]*xij[1]+xij[2]*xij[2]+pars->eps2);
    const double mi = pi.mass;
    const double mj = pj.mass;
    
    const double mimj = mi*mj; // m_i*m_i
    const double mmij = mimj;

    Pij = - mimj*invrij;
    Wij =   mmij*invrij;   // Transformation coefficient
    
    const double invrij3=invrij*invrij*invrij;
    double mor3 = mj * invrij3;
    Aij[0] = mor3 * xij[0];   
    Aij[1] = mor3 * xij[1]; 
    Aij[2] = mor3 * xij[2]; 

    mor3 = mmij * invrij3;
    pWij[0] = mor3 * xij[0];
    pWij[1] = mor3 * xij[1];
    pWij[2] = mor3 * xij[2];

  return 0;
}

/// start Newtonian cut force (L.Wang)
//! Newtonian acceleration with cutoff function and \f$\partial W_{ij}/\partial \mathbf{x}_i\f$ from particle j to particle i (function type of \link #ARC::pair_AW \endlink) 
/*!          @param[out] Aij: Newtonian acceleration vector for i particle from particle j. \f$Aij[1:3] = m_i m_j xij[1:3] / |xij|^3 (1-k) + Pij (1-kdot[1:3])\f$.
             @param[out] Pij: Newtonian potential of i from j. \f$ Pij = -m_i m_j /|xij| (1-k)\f$
             @param[out] pWij: TTL time transformation function partial derivates (component from j to i) \f$\partial W_{ij}/\partial \mathbf{x}_i\f$ (used for TTL method). \f$pWij[1:3] = mm_{ij} xij[1:3] /|xij|^3 \f$. (Total value is \f$\frac{\partial W}{\partial \mathbf{x}_i} = \sum_{j} mm_{ij} \mathbf{x}_{ij}/|\mathbf{x}_{ij}|^3\f$)
             @param[out] Wij: TTL time transformation function component with i,j (used for TTL method) \f$Wij = mm_{ij} /|xij|^3\f$ total value is \f$ W = \sum_{i<j} mm_{ij} /|xij| \f$
             @param[in] xij: relative position vector [1:3] \f$ \mathbf{x_j} - \mathbf{x_i} \f$
             @param[in] pi: particle i.
             @param[in] pj: particle j.
             @param[in] smpars: array of double[2]; First element is rcut_out, second element is rcut_in.
             \return status 0
*/
template <class Tptcl, class extpar>
int Newtonian_cut_AW (double Aij[3], double &Pij, double pWij[3], double &Wij, const double xij[3], const Tptcl &pi, const Tptcl &pj, extpar* pars) {
  // distance
  const double rij = std::sqrt(xij[0]*xij[0]+xij[1]*xij[1]+xij[2]*xij[2]+pars->eps2);
  const double mi = pi.mass;
  const double mj = pj.mass;

  // smpars[2:3]: rcut_out, rcut_in
//  const double k   = cutoff_poly_3rd(rij, smpars[0], smpars[1]);
//  const double kdx = cutoff_poly_3rd_dr(rij, xij[0], smpars[0], smpars[1]);
//  const double kdy = cutoff_poly_3rd_dr(rij, xij[1], smpars[0], smpars[1]);
//  const double kdz = cutoff_poly_3rd_dr(rij, xij[2], smpars[0], smpars[1]);
  //const double r_out = std::max(pi.r_out,pj.r_out);
  const double r_out = pars->rout;
  const double r_in  = pars->rin;
  const double r_oi_inv = pars->r_oi_inv;
  const double r_A   = pars->r_A;

#ifdef INTEGRATED_CUTOFF_FUNCTION
  const double kpot  = 1.0-CalcW(rij/r_out, r_in/r_out);
#else
  const double pot_off = pars->pot_off;
  const double kpot  = cutoff_pot(rij, r_oi_inv, r_A, r_in);
#endif
  const double k     = cutoff_poly_3rd(rij, r_oi_inv, r_A, r_in);

  // smooth coefficients
//  const double mm2=smpars[0];
//  const double epi=smpars[1];

  // mass parameters
  const double mimj = mi*mj; // m_i*m_i
  const double mmij = mimj;
//#ifdef HARD_DEBUG
//  if(mmij==0.5) mmij = 0.0;
//#endif
//  if (mm2>0 && epi>0) {
//    // Wij = mm2 if m_i*m_i < epi*m'^2; 0 otherwise;
//    if (mimj<epi*mm2) mmij = mm2;
//    else mmij = 0;
//  }
//  else {
//    mmij = mimj;    // Wij = m_i*m_i
//  }
  
#ifdef INTEGRATED_CUTOFF_FUNCTION
  Pij = - mimj / rij * kpot;  // Potential energy
#else
  Pij = - mimj*(1.0/ rij * kpot - pot_off);
  if (rij>r_out) Pij = 0.0;
#endif
//#ifdef DEBUG
//  std::cerr<<"Pij = "<<Pij<<" rij = "<<rij<<" k="<<k<<" r_in="<<r_in<<" r_out="<<r_out<<std::endl;
//#endif
  Wij = mmij / rij;   // Transformation coefficient
        
  // Acceleration
  const double rij3 = rij*rij*rij;
  double mor3 = mj / rij3;
//  Aij[0] = mor3 * xij[0] * (1-k) + Pij*(1-kdx);
//  Aij[1] = mor3 * xij[1] * (1-k) + Pij*(1-kdy);
//  Aij[2] = mor3 * xij[2] * (1-k) + Pij*(1-kdz);
  Aij[0] = mor3 * xij[0] * k;   
  Aij[1] = mor3 * xij[1] * k; 
  Aij[2] = mor3 * xij[2] * k; 

  // dW/dr
  mor3 = mmij / rij3;
  pWij[0] = mor3 * xij[0];
  pWij[1] = mor3 * xij[1];
  pWij[2] = mor3 * xij[2];

  return 0;
}

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
void Newtonian_extA_pert (double3* acc, const PS::F64 time, Tptcl* p, const PS::S32 np, Tpert* pert, Tforce* pf, const PS::S32 npert, extpar* pars){
#ifdef HARD_DEBUG
    if(npert<=1) {
        std::cerr<<"Error: perturboer number not enough !"<<std::endl;
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
    assert(mt==pert[0]->mass);
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
}
/// end Newtonian cut force (L.Wang)

template<class Tptcl, class Tpert, class Tforce, class extpar>
void Newtonian_extA_soft (double3* acc, const PS::F64 time, Tptcl* p, const PS::S32 np, Tpert* pert, Tforce* pf, const PS::S32 npert, extpar* pars){
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
}

// period of two-body motion
template<class extpar>
double Newtonian_timescale(const double m1, const double m2, const double dx[], const double dv[], extpar* par) {
    const double dr2  = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
    const double dv2  = dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2];
    const double dr   = std::sqrt(dr2);
    const double m    = m1+m2;

    const double semi = 1.0/(2.0/dr - dv2/m);

    if (semi<0) {
        const double peri = 0.1*std::sqrt(dr2*dr/(2.0*m));
        //      std::cout<<"dr="<<dr<<" semi="<<semi<<" peri="<<peri<<std::endl;
        return peri;
    }
    else {
        const double rdot = dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2];
        const double dr_semi = 1.0 - dr/semi;
        const double ecc  = std::sqrt(dr_semi*dr_semi + rdot*rdot/(m*semi));

        const double twopi= 6.28;
        const double peri = twopi*std::abs(semi)*std::sqrt(std::abs(semi)/m);

        //      std::cout<<"dr="<<dr<<" semi="<<semi<<" ecc="<<ecc<<" peri="<<peri<<std::endl;
        return std::max(std::sqrt(std::abs(1.0-ecc)),0.01)*peri;
    }
}

inline void CalcAcc0Acc1R2Cutoff(const PS::F64vec posi,
                                 const PS::F64vec veli,
                                 PS::F64vec & acci,
                                 PS::F64vec & jrki,
                                 PS::F64 & r2,
                                 const PS::F64vec posj, 
                                 const PS::F64vec velj,
                                 const PS::F64 massj,
                                 const PS::F64 eps2,
                                 const PS::F64 rcut_out,
                                 const PS::F64 rcut_in,
                                 const PS::F64 rcut_oi_inv,
                                 const PS::F64 rcut_A){
    const PS::F64vec rij = posi - posj;
    r2 = rij*rij;
    const PS::F64 r2_eps = r2 + eps2;
    if(r2_eps <= rcut_out*rcut_out){
        const PS::F64vec vij = veli - velj;
        const PS::F64 rijvij = rij * vij;
        const PS::F64 r_eps = sqrt(r2_eps);
        const PS::F64 R = 1.0/r_eps;
        //const PS::F64 R = 1.0 / sqrt(r2_eps);
        const PS::F64 R2 = R*R;
        const PS::F64 R3 = R2*R;
        const PS::F64 A = (rijvij)*R2;
        const PS::F64 K = cutoff_poly_3rd(r_eps, rcut_oi_inv, rcut_A, rcut_in);
        const PS::F64 Kdot = cutoff_poly_3rd_dot(r_eps, rijvij, rcut_oi_inv, rcut_A, rcut_in);
        const PS::F64 MR3 = massj*R3;
        const PS::F64vec F0 = -MR3*rij*K;
        const PS::F64vec F1 = -MR3*vij*K - 3.0*A*F0 - MR3*rij*Kdot;
        acci += F0;
        jrki += F1;
    }
}

#ifdef HARD_DEBUG_DEEP_CHECK
template<class Tptcl, class Tpert, class Tforce, class extpar>
void Newtonian_extA_test (double3* acc, const PS::F64 time, Tptcl* p, const PS::S32 np, Tpert* pert, Tforce* pf, const PS::S32 npert, extpar* pars){
    if(npert>1) {
        static const PS::F64 inv3 = 1.0 / 3.0;
        PS::F64vec xp[npert];

        ARC::chainpars ARC_control;
        ARC_control.setA(Newtonian_cut_AW<Tptcl,extpar>,Newtonian_extA<Tptcl,Tpert,Tforce,extpar>,Newtonian_timescale<extpar>);
        ARC_control.setabg(0,1,0);
        ARC_control.setErr(1e-10,1e-24,1e-6);
        ARC_control.setIterSeq(20,3,20);
        ARC_control.setIntp(1);
        ARC_control.setIterConst(0);
        ARC_control.setAutoStep(3);

        extpar Int_pars;
        Int_pars.rin = pars->rin;
        Int_pars.rout = pars->rout;
        Int_pars.eps2 = pars->eps2;
        
        ARC::chain<Tptcl> c(3);
        Tptcl pc[3];
        pc[0] = p[0];
        pc[1] = p[1];
        pc[2] = Tptcl(*pert[1]);
        c.linkP(3,pc);

        pc[0].pos += pert[0]->pos;
        pc[1].pos += pert[0]->pos;

        c.init(pert[0]->time,ARC_control,&Int_pars);
        PS::F64 ds_up_limit = 0.25*(time-pert[0]->time)/c.calc_dt_X(1.0,ARC_control);
        PS::F64 ds = std::min(ds_up_limit,c.calc_next_step_custom(ARC_control, &Int_pars));
        
        Tpert* tp;
        Tforce* tf;
        PS::F64 dscoff=1.0;
        PS::S32 converge_count=0;
        PS::S32 error_count=0;
        bool modify_step_flag=false;
        bool final_flag=false;

        while(time-c.getTime()>ARC_control.dterr*c.getTime()) {
            PS::F64 dsf = c.extrapolation_integration(ds, ARC_control, time, &Int_pars, tp, tf, 0);
            if (dsf<0) {
                final_flag=true;
                converge_count++;
                if (converge_count>5&&time-c.getTime()>ARC_control.dterr*100) {
                    std::cerr<<"Error: Time synchronization fails!\nStep size ds: "<<ds<<"\nEnding physical time: "<<time<<"\nTime difference: "<<time-c.getTime()<<"\nR_in: "<<Int_pars.rin<<"\nR_out: "<<Int_pars.rout<<"\n";
                    c.dump("ARC_dump.dat");
                    ARC_control.dump("ARC_dump.par");
                    c.print(std::cerr);
                    abort();
                }
                else ds *= -dsf;
            }
            else if (dsf==0) {
                c.info->ErrMessage(std::cerr);
                error_count++;
                if(error_count>4) {
                    std::cerr<<"Error: Too much error appear!\nStep size ds: "<<ds<<"\nEnding physical time: "<<time<<"\nTime difference: "<<time-c.getTime()<<"\nR_in: "<<Int_pars.rin<<"\nR_out: "<<Int_pars.rout<<"\n";
                    c.dump("ARC_dump.dat");
                    ARC_control.dump("ARC_dump.par");
                    c.print(std::cerr);
                    abort();
                }
                if (c.info->status==5) {
                    dscoff = 0.25;
                    ds *= dscoff;
                }
                else if (c.info->status==4) ds = std::min(dscoff*c.calc_next_step_custom(ARC_control, &Int_pars),ds_up_limit);
                else ds *= 0.1;
                modify_step_flag=true;
            }
            else  {
                if (final_flag) {
                    if (converge_count>6&&time-c.getTime()>ARC_control.dterr*100) {
                        std::cerr<<"Error: Time synchronization fails!\nStep size ds: "<<ds<<"\nEnding physical time: "<<time<<"\nTime difference: "<<time-c.getTime()<<"\nR_in: "<<Int_pars.rin<<"\nR_out: "<<Int_pars.rout<<"\n";
                        c.dump("ARC_dump.dat");
                        ARC_control.dump("ARC_dump.par");
                        c.print(std::cerr);
                        abort();
                    }
                    converge_count++;
                }
                else if (modify_step_flag&&error_count==0) {
                    ds = std::min(dscoff*c.calc_next_step_custom(ARC_control, &Int_pars),ds_up_limit);
                    modify_step_flag=false;
                }
                // reducing error counter if integration success, this is to avoid the significant change of step may cause some issue
                if(error_count>0) error_count--;
            }
        }

            
        c.resolve();
        xp[0] = (pc[0].pos*pc[0].mass + pc[1].pos*pc[1].mass)/(pc[0].mass + pc[1].mass);
        xp[1] = pc[2].pos;

        for(int i=0; i<npert; i++) {
            PS::F64 dt = time - pert[i]->time;
            xp[i] = pert[i]->pos + dt*(pert[i]->vel + 0.5*dt*(pert[i]->acc0 + inv3*dt*pert[i]->acc1));
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
        assert(mt==pert[0]->mass);
#endif

        for(int i=0; i<np; i++) {
            PS::F64vec xi = p[i].pos + xp[0];
            acc[i][0] = -pf[0]->acc0[0]; 
            acc[i][1] = -pf[0]->acc0[1];        
            acc[i][2] = -pf[0]->acc0[2]; 
//            acc[i][0] = acc[i][1] = acc[i][2] = 0.0;
            for(int j=1; j<npert; j++) {
            
                PS::F64vec dx = xp[j] - xi;
                //std::cerr<<"i = "<<i<<" j = "<<j<<" dx = "<<dx<<std::endl;
                //PS::F64 mi = p[i].mass;
                PS::F64 mp = pert[j]->mass;
                PS::F64 dr2 = dx*dx + pars->eps2;
                PS::F64 dr  = std::sqrt(dr2);
                PS::F64 dr3 = dr*dr2;
                PS::F64 mor3 = mp/dr3;

                // smpars[2:3]: rcut_out, rcut_in
                //  const double k   = cutoff_poly_3rd(dr, smpars[0], smpars[1]);
                //  const double kdx = cutoff_poly_3rd_dr(dr, dx, smpars[0], smpars[1]);
                //  const double kdy = cutoff_poly_3rd_dr(dr, dy, smpars[0], smpars[1]);
                //  const double kdz = cutoff_poly_3rd_dr(dr, dz, smpars[0], smpars[1]);  
                const PS::F64 r_out = pars->rout;
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
        }
    }
    else {
        for(int i=0; i<np; i++) acc[i][0] = acc[i][1] = acc[i][2] = 0.0;
    }

    // soft perturbation
}
#endif

