#pragma once
class HardForce{
    PS::F64 rin_;
    PS::F64 rout_;
    PS::F64 eps2_;
public:
    void setParam(const PS::F64 _rin,
                  const PS::F64 _rout,
                  const PS::F64 _eps2 = 0.0){
        rin_ = _rin;
        rout_ = _rout;
        eps2_ = _eps2;
    }
};
inline PS::F64 cutoff_poly_3rd(const PS::F64 rij,
                               const PS::F64 rout,
                               const PS::F64 rin){
    PS::F64 inv_dr = 1.0 / (rout-rin);
    PS::F64 x = (rij - rin)*inv_dr;
    x = (x < 1.0) ? x : 1.0;
    x = (x > 0.0) ? x : 0.0;
    PS::F64 x2 = x*x;
    PS::F64 x4 = x2*x2;
    PS::F64 k = (((-20.0*x+70.0)*x-84.0)*x+35.0)*x4;
    return k;
}


inline PS::F64 cutoff_poly_3rd(const PS::F64 rij,
                               const PS::F64 rout,
                               const PS::F64 rin,
                               const PS::F64 inv_dr){
    PS::F64 x = (rij - rin)*inv_dr;
    x = (x < 1.0) ? x : 1.0;
    x = (x > 0.0) ? x : 0.0;
    PS::F64 x2 = x*x;
    PS::F64 x4 = x2*x2;
    PS::F64 k = (((-20.0*x+70.0)*x-84.0)*x+35.0)*x4;
    return k;
}

/// start cutoff dr (L.Wang)
inline PS::F64 cutoff_poly_3rd_dr(const PS::F64 rij,
                                  const PS::F64 rx,
                                  const PS::F64 rout,
                                  const PS::F64 rin){
    PS::F64 inv_dr = 1.0 / (rout-rin);
    PS::F64 x = (rij - rin)*inv_dr;
    x = (x < 1.0) ? x : 1.0;
    x = (x > 0.0) ? x : 0.0;
    PS::F64 x2 = x*x;
    PS::F64 x3 = x2*x;
    PS::F64 k = (((-140.0*x+420.0)*x-420.0)*x+140.0)*x3*inv_dr*rx/rij;
    return k;
}
/// end cutoff dr (L.Wang)

inline void CalcPotShort(const PS::F64vec & posi,
			 PS::F64 & poti_tot,
			 const PS::F64vec & posj,
			 const PS::F64 massj,
			 const PS::F64 eps2,
			 const PS::F64 rcut_out,
			 const PS::F64 rcut_in){
    const PS::F64vec rij = posi - posj;
    const PS::F64 r2 = rij * rij;
    const PS::F64 r2_eps = r2 + eps2;
    //const PS::F64 inv_dr = 1.0 / (rcut_out-rcut_in);
    PS::F64 R = 1.0/sqrt(r2_eps);
    PS::F64 r_eps = R * r2_eps;
    //PS::F64 k = cutoff_poly_3rd(r_eps, rcut_out, rcut_in, inv_dr);
    PS::F64 k = cutoff_poly_3rd(r_eps, rcut_out, rcut_in);
    poti_tot -= massj * R * (1.0-k);
}

inline void CalcAccPotShortWithLinearCutoff(const PS::F64vec & posi,
					    PS::F64vec & acci_pla,
					    PS::F64 & poti_tot,
					    const PS::F64vec & posj,
					    const PS::F64 massj,
					    const PS::F64 eps2,
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
    const PS::F64 k = cutoff_poly_3rd(r_eps, rcut_out, rcut_in);
    const PS::F64 r2_max = (r2_eps > rcut2_out) ? r2_eps : rcut2_out;
    const PS::F64 R_max = 1.0/sqrt(r2_max);
    const PS::F64 Rm_max = massj * R_max;
    const PS::F64 R2_max = R_max * R_max;
    const PS::F64 Rm3_max = Rm_max * R2_max;
    poti_tot -= (Rm - Rm_max);
    acci_pla -= (Rm3*k - Rm3_max)*rij;
}


inline PS::F64 cutoff_poly_3rd_dot(const PS::F64 &rij,
                                   const PS::F64 &rijvij,
                                   const PS::F64 &_rout,
                                   const PS::F64 &_rin){
    PS::F64 rout = _rout;
    PS::F64 rin = _rin;
    PS::F64 inv_dr = 1.0/(rout-rin);
    PS::F64 x = (rij - rin)*inv_dr;
    PS::F64 xdot = rijvij/rij*inv_dr;
    PS::F64 Kdot = 0.0;
    if(x <= 0.0)
        Kdot = 0.0;
    else if(1.0 <= x)
        Kdot = 0.0;
    else{
        PS::F64 x2 = x*x;
        PS::F64 x3 = x2*x;
        PS::F64 x4 = x2*x2;
        PS::F64 x5 = x4*x;
        PS::F64 x6 = x4*x2;
        Kdot = (-140.0*x6 + 420.0*x5 - 420.0*x4 + 140.0*x3) * xdot;
    }
    return Kdot;
}

#if 1
inline void CalcAcc0AndAcc1Cutoff(const PS::F64vec posi,
                                  const PS::F64vec veli,
                                  PS::F64vec & acci,
                                  PS::F64vec & jrki,
                                  const PS::F64vec posj, 
                                  const PS::F64vec velj, 
                                  const PS::F64 massj, 
                                  const PS::F64 eps2, 
                                  const PS::F64 rcut_out,
                                  const PS::F64 rcut_in){
    const PS::F64vec rij = posi - posj;
    const PS::F64 r2_eps = rij*rij + eps2;
    if(r2_eps <= rcut_out*rcut_out){
        const PS::F64vec vij = veli - velj;
        const PS::F64 rijvij = rij * vij;
        const PS::F64 r_eps = sqrt(r2_eps);
        const PS::F64 R = 1.0/r_eps;
	//const PS::F64 R = 1.0 / sqrt(r2_eps);
        const PS::F64 R2 = R*R;
        const PS::F64 R3 = R2*R;
        const PS::F64 A = (rijvij)*R2;
#ifdef FORDEBUG
        PS::F64 K = 0.0; // for debug
        PS::F64 Kdot = 0.0; // for debug
#else
        const PS::F64 K = cutoff_poly_3rd(r_eps, rcut_out, rcut_in);
        const PS::F64 Kdot = cutoff_poly_3rd_dot(r_eps, rijvij, rcut_out, rcut_in);
#endif
        const PS::F64vec F0 = -massj*R3*rij*(1.0-K);
        const PS::F64vec F1 = -massj*R3*vij*(1.0-K) - 3.0*A*F0 + massj*R3*rij*Kdot;
        acci += F0;
        jrki += F1;
    }
}
#endif

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
                                 const PS::F64 rcut_in){
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
#ifdef FORDEBUG
        PS::F64 K = 0.0; // for debug
        PS::F64 Kdot = 0.0; // for debug
#else
        const PS::F64 K = cutoff_poly_3rd(r_eps, rcut_out, rcut_in);
        const PS::F64 Kdot = cutoff_poly_3rd_dot(r_eps, rijvij, rcut_out, rcut_in);
#endif
        const PS::F64vec F0 = -massj*R3*rij*(1.0-K);
        const PS::F64vec F1 = -massj*R3*vij*(1.0-K) - 3.0*A*F0 + massj*R3*rij*Kdot;
        acci += F0;
        jrki += F1;
    }
}

inline void CalcAcc0Acc1AndR2CutoffPair(const PS::F64 massi,
					const PS::F64vec posi,
					const PS::F64vec veli,
					PS::F64vec & acci,
					PS::F64vec & jrki,
					const PS::F64 massj,
					const PS::F64vec posj,
					const PS::F64vec velj,
					PS::F64vec & accj,
					PS::F64vec & jrkj,
					PS::F64 & r2,
					const PS::F64 eps2,
					const PS::F64 rcut_out,
					const PS::F64 rcut_in){
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
	//#ifdef FORDEBUG
#if 0
        PS::F64 K = 0.0; // for debug
        PS::F64 Kdot = 0.0; // for debug
#else
        const PS::F64 K = cutoff_poly_3rd(r_eps, rcut_out, rcut_in);
        const PS::F64 Kdot = cutoff_poly_3rd_dot(r_eps, rijvij, rcut_out, rcut_in);
#endif
        const PS::F64vec F0 = -R3*rij*(1.0-K);
        const PS::F64vec F1 = -R3*vij*(1.0-K) - 3.0*A*F0 + R3*rij*Kdot;
        acci += massj*F0;
        jrki += massj*F1;
	accj -= massi*F0;
        jrkj -= massi*F1;
    }
}

inline void CalcAcc0Acc1AndPotCutoff(const PS::F64vec posi,
                                  const PS::F64vec veli,
                                  PS::F64vec & acci,
                                  PS::F64vec & jrki,
                                  PS::F64 & poti,
                                  const PS::F64vec posj, 
                                  const PS::F64vec velj, 
                                  const PS::F64 massj,
                                  const PS::F64 eps2,
                                  const PS::F64 rcut_out,
                                  const PS::F64 rcut_in){
    const PS::F64vec rij = posi - posj;
    const PS::F64 r2_eps = rij*rij + eps2;
#ifdef FORDEBUG
    if(1){
#else
    if(r2_eps <= rcut_out*rcut_out){
#endif
        const PS::F64vec vij = veli - velj;
        const PS::F64 rijvij = rij * vij;
        const PS::F64 r_eps = sqrt(r2_eps);
        const PS::F64 R = 1.0/r_eps;
	//const PS::F64 R = 1.0 / sqrt(r2_eps);
        const PS::F64 R2 = R * R;
        const PS::F64 R3 = R2 * R;
        const PS::F64 A = rijvij * R2;
#ifdef FORDEBUG
        PS::F64 K = 0.0; // for debug
        PS::F64 Kdot = 0.0; // for debug
        poti -= massj * R;
#else
        const PS::F64 K = cutoff_poly_3rd(r_eps, rcut_out, rcut_in);
        const PS::F64 Kdot = cutoff_poly_3rd_dot(r_eps, rijvij, rcut_out, rcut_in);
#endif
        const PS::F64vec F0 = -massj*R3*rij*(1.0-K);
        const PS::F64vec F1 = -massj*R3*vij*(1.0-K) - 3.0*A*F0 + massj*R3*rij*Kdot;
        acci += F0;
        jrki += F1;
    }
}

inline void CalcAcc0Acc1(const PS::F64vec posi,
                         const PS::F64vec veli,
                         PS::F64vec & acci,
                         PS::F64vec & jrki,
                         const PS::F64vec posj,
                         const PS::F64vec velj,
                         const PS::F64 massj,
                         const PS::F64 eps2=0.0){
    const PS::F64vec rij = posi - posj;
    const PS::F64vec vij = veli - velj;
    const PS::F64 r2 = rij*rij + eps2;
    const PS::F64 r_inv = 1.0 / sqrt(r2);
    const PS::F64 r2_inv = r_inv * r_inv;
    const PS::F64 r3_inv = r2_inv * r_inv;
    const PS::F64 m_r3 = massj * r3_inv;
    const PS::F64vec F0 = -m_r3*rij;
    const PS::F64vec F1 = -m_r3*vij - 3.0*rij*vij*r2_inv*F0;
    acci += F0;
    jrki += F1;
}

inline void CalcAcc0Acc1AndR2(const PS::F64vec posi,
                              const PS::F64vec veli,
                              PS::F64vec & acci,
                              PS::F64vec & jrki,
                              const PS::F64vec posj,
                              const PS::F64vec velj,
                              const PS::F64 massj,
                              const PS::F64 eps2){
    const PS::F64vec rij = posi - posj;
    const PS::F64vec vij = veli - velj;
    //const PS::F64 r2_eps = rij*rij + eps2;
    const PS::F64 r_inv = 1.0 / sqrt(rij * rij);
    const PS::F64 r2_inv = r_inv * r_inv;
    const PS::F64 r3_inv = r2_inv * r_inv;
    const PS::F64 m_r3 = massj * r3_inv;
    const PS::F64vec F0 = -m_r3*rij;
    const PS::F64vec F1 = -m_r3*vij - 3.0*rij*vij*r2_inv*F0;
    acci += F0;
    jrki += F1;
}

/* 
class CalcW{
private:
    PS::F64 A7, A6, A5, A4, A3, A2, A1, A0, B1, A1_dash;
    PS::F64 q; // rin/rout
public:
    void setParam(const PS::F64 r_in, const PS::F64 r_out){
    }
};
*/

//#ifdef CALC_HARD_ENERGY
// y: reps/rout [reps: sqrt(r^2+eps^2)], q: rin/rout
inline PS::F64 CalcW(const PS::F64 y, const PS::F64 q=0.1){
     PS::F64 q2 = q*q;
     PS::F64 q3 = q2*q;
     PS::F64 q4 = q2*q2;
     PS::F64 q5 = q3*q2;
     PS::F64 q6 = q3*q3;
     PS::F64 q7 = q4*q3;
     PS::F64 denominator = (q-1.0)*(q-1.0)*(q-1.0)*(q-1.0)*(q-1.0)*(q-1.0)*(q-1.0);
     PS::F64 A7 = 20.0/denominator/-6;
     PS::F64 A6 = (-70.0*q - 70.0)/denominator/-5;
     PS::F64 A5 = (84.0*q2 + 252.0*q + 84.0)/denominator/-4;
     PS::F64 A4 = (-35.0*q3 - 315.0*q2 - 315.0*q - 35.0)/denominator/-3;
     PS::F64 A3 = (140.0*q3 + 420.0*q2 + 140.0*q)/denominator/-2;
     PS::F64 A2 = (-210*q3 - 210.0*q2)/denominator/-1;
     PS::F64 A1 = (140*q3)/denominator*-1;
     PS::F64 A0 = (-35.0*q4 + 21.0*q5 - 7.0*q6 + q7)/denominator;
     PS::F64 x = 1.0; // x=rout/rout
     PS::F64 B1 = 1.0 - ( (((((((A7*x + A6)*x + A5)*x + A4)*x + A3)*x + A2)*x + A1*log(x))*x) + A0 ); // to W(r>rout) = 1.0
     PS::F64 A1_dash = -7*(60*q3*log(q) - q6 + 9.0*q5 - 45.0*q4 + 45.0*q2 - 9.0*q + 1.0)/(3.0*denominator);
     if(y <= q) return A1_dash*y;
     else if(y >= 1.0) return 1.0;
     else return (((((((A7*y + A6)*y + A5)*y + A4)*y + A3)*y + A2)*y + A1*log(y) + B1)*y) + A0;
}
//#endif

template<class Tptcl, class Teng>
void CalcEnergyHard(const Tptcl ptcl[], const PS::S32 n_tot, Teng & eng, 
		    const PS::F64 r_out, const PS::F64 r_in,
		    const PS::F64 eps_sq = 0.0){
    eng.kin = eng.pot = eng.tot = 0.0;
    for(PS::S32 i=0; i<n_tot; i++){
	eng.kin += 0.5 * ptcl[i].mass * ptcl[i].vel * ptcl[i].vel;
	for(PS::S32 j=i+1; j<n_tot; j++){
	    PS::F64vec rij = ptcl[i].pos - ptcl[j].pos;
	    PS::F64 dr = sqrt(rij*rij + eps_sq);
	    eng.pot -= ptcl[i].mass/dr*(1.0 - CalcW(dr/r_out, r_in/r_out));
	}
    }
    eng.tot = eng.kin + eng.pot;
}

template<class Tptcl, class Teng>
void CalcEnergyHard(const Tptcl ptcl[], const PS::S32 n_tot, Teng & eng, 
		    const PS::F64 r_out, const PS::F64 r_in,
		    const PS::F64 mass_sun, const PS::F64vec pos_sun, const PS::F64vec vel_sun,
		    const PS::F64 eps_sq = 0.0){
    eng.kin = eng.pot = eng.tot = 0.0;
    for(PS::S32 i=0; i<n_tot; i++){
	eng.kin += 0.5 * ptcl[i].mass * ptcl[i].vel * ptcl[i].vel;
	for(PS::S32 j=i+1; j<n_tot; j++){
	    PS::F64vec rij = ptcl[i].pos - ptcl[j].pos;
	    PS::F64 dr = sqrt(rij*rij + eps_sq);
	    eng.pot -= ptcl[i].mass*ptcl[j].mass/dr*(1.0 - CalcW(dr/r_out, r_in/r_out));
	}
	PS::F64vec ri0 = ptcl[i].pos - pos_sun;
	eng.pot -= ptcl[i].mass*mass_sun / sqrt(ri0*ri0);
    }
    eng.tot = eng.kin + eng.pot;
}


template<class Tptcl, class Teng>
void CalcEnergyHardOneParticle(const Tptcl ptcl[], const PS::S32 n_tot, Teng & eng, 
			       const PS::F64 r_out, const PS::F64 r_in,
			       const PS::F64 mass_sun, const PS::F64vec pos_sun, const PS::F64vec vel_sun,
			       const PS::S32 adr_i,
			       const PS::F64 eps_sq = 0.0){
    eng.kin = eng.pot = eng.tot = 0.0;
    eng.kin += 0.5 * ptcl[adr_i].mass * ptcl[adr_i].vel * ptcl[adr_i].vel;
    for(PS::S32 j=0; j<n_tot; j++){
	if(adr_i == j) continue;
	PS::F64vec rij = ptcl[adr_i].pos - ptcl[j].pos;
	PS::F64 dr = sqrt(rij*rij + eps_sq);
	eng.pot -= ptcl[adr_i].mass*ptcl[j].mass/dr*(1.0 - CalcW(dr/r_out, r_in/r_out));
    }
    PS::F64vec ri0 = ptcl[adr_i].pos - pos_sun;
    eng.pot -= ptcl[adr_i].mass*mass_sun / sqrt(ri0*ri0);
    eng.tot = eng.kin + eng.pot;
}


/// start Newtonian cut force (L.Wang)
//! Newtonian acceleration with cutoff function and \f$\partial W_{ij}/\partial \mathbf{x}_i\f$ from particle j to particle i (function type of \link #ARC::pair_AW \endlink) 
/*!          @param[out] Aij: Newtonian acceleration vector for i particle from particle j. \f$Aij[1:3] = m_i m_j xij[1:3] / |xij|^3 (1-k) + Pij (1-kdot[1:3])\f$.
             @param[out] Pij: Newtonian potential of i from j. \f$ Pij = -m_i m_j /|xij| (1-k)\f$
             @param[out] pWij: TTL time transformation function partial derivates (component from j to i) \f$\partial W_{ij}/\partial \mathbf{x}_i\f$ (used for TTL method). \f$pWij[1:3] = mm_{ij} xij[1:3] /|xij|^3 \f$. (Total value is \f$\frac{\partial W}{\partial \mathbf{x}_i} = \sum_{j} mm_{ij} \mathbf{x}_{ij}/|\mathbf{x}_{ij}|^3\f$)
             @param[out] Wij: TTL time transformation function component with i,j (used for TTL method) \f$Wij = mm_{ij} /|xij|^3\f$ total value is \f$ W = \sum_{i<j} mm_{ij} /|xij| \f$
             @param[in] xij: relative position vector [1:3] \f$ \mathbf{x_j} - \mathbf{x_i} \f$
             @param[in] mi: particle i mass.
             @param[in] mj: particle j mass.
             @param[in] smpars: array of double[2]; First element is rcut_out, second element is rcut_in.
*/
void Newtonian_cut_AW (double Aij[3], double &Pij, double pWij[3], double &Wij, const double xij[3], const double &mi, const double &mj, const double* smpars) {
  // distance
  const double rij = std::sqrt(xij[0]*xij[0]+xij[1]*xij[1]+xij[2]*xij[2]);  

  // smpars[2:3]: rcut_out, rcut_in
//  const double k   = cutoff_poly_3rd(rij, smpars[0], smpars[1]);
//  const double kdx = cutoff_poly_3rd_dr(rij, xij[0], smpars[0], smpars[1]);
//  const double kdy = cutoff_poly_3rd_dr(rij, xij[1], smpars[0], smpars[1]);
//  const double kdz = cutoff_poly_3rd_dr(rij, xij[2], smpars[0], smpars[1]);
  const double r_out = smpars[0];
  const double r_in  = smpars[1];
  const double k = CalcW(rij/r_out, r_in/r_out);
  const double kdot = cutoff_poly_3rd(rij, r_out, r_in);

  // smooth coefficients
//  const double mm2=smpars[0];
//  const double epi=smpars[1];

  // mass parameters
  const double mimj = mi*mj; // m_i*m_i
  const double mmij = mimj;
//  if (mm2>0 && epi>0) {
//    // Wij = mm2 if m_i*m_i < epi*m'^2; 0 otherwise;
//    if (mimj<epi*mm2) mmij = mm2;
//    else mmij = 0;
//  }
//  else {
//    mmij = mimj;    // Wij = m_i*m_i
//  }
  
  Pij = - mimj / rij * (1-k);  // Potential energy
  Wij = mmij / rij;   // Transformation coefficient
        
  // Acceleration
  const double rij3 = rij*rij*rij;
  double mor3 = mj / rij3;
//  Aij[0] = mor3 * xij[0] * (1-k) + Pij*(1-kdx);
//  Aij[1] = mor3 * xij[1] * (1-k) + Pij*(1-kdy);
//  Aij[2] = mor3 * xij[2] * (1-k) + Pij*(1-kdz);
  Aij[0] = mor3 * xij[0] * (1-kdot);   
  Aij[1] = mor3 * xij[1] * (1-kdot); 
  Aij[2] = mor3 * xij[2] * (1-kdot); 

  // dW/dr
  mor3 = mmij / rij3;
  pWij[0] = mor3 * xij[0];
  pWij[1] = mor3 * xij[1];
  pWij[2] = mor3 * xij[2];
  
}

//! Newtonian acceleration from particle p to particle i (function type of ::ARC::pair_Ap)
/*! 
  @param[out]  Aij: acceleration vector. \f$Aij[1:3] = m_i m_p (xp[1:3]-xi[1:3]) / |xp-xi|^3 \f$.
  @param[out]  Pij: potential. \f$ Pij = - m_i m_p /|xp-xi|^3\f$
  @param[in]  xi: position vector i.
  @param[in]  xp: position vector p.
  @param[in]  mi: particle mass i.
  @param[in]  mp: particle mass p.
  @param[in] smpars: array of double[2]; First element is rcut_out, second element is rcut_in.
 */
void Newtonian_cut_Ap (double Aij[3], double &Pij, const double xi[3], const double xp[3], const double &mi, const double &mp, const double* smpars){
  double dx = xp[0] - xi[0];
  double dy = xp[1] - xi[1];
  double dz = xp[2] - xi[2];

  double dr2 = dx*dx + dy*dy + dz*dz;
  double dr  = std::sqrt(dr2);
  double dr3 = dr*dr2;

  // smpars[2:3]: rcut_out, rcut_in
//  const double k   = cutoff_poly_3rd(dr, smpars[0], smpars[1]);
//  const double kdx = cutoff_poly_3rd_dr(dr, dx, smpars[0], smpars[1]);
//  const double kdy = cutoff_poly_3rd_dr(dr, dy, smpars[0], smpars[1]);
//  const double kdz = cutoff_poly_3rd_dr(dr, dz, smpars[0], smpars[1]);  
  const double r_out = smpars[0];
  const double r_in  = smpars[1];
  const double k = CalcW(dr/r_out, r_in/r_out);
  const double kdot = cutoff_poly_3rd(dr, r_out, r_in);

  Pij = - mi*mp / dr * (1-k);

//  Aij[0] = mp * dx / dr3 * (1-k) + Pij * (1-kdx);
//  Aij[1] = mp * dy / dr3 * (1-k) + Pij * (1-kdy);
//  Aij[2] = mp * dz / dr3 * (1-k) + Pij * (1-kdz);
  Aij[0] = mp * dx / dr3 * (1-kdot);
  Aij[1] = mp * dy / dr3 * (1-kdot);
  Aij[2] = mp * dz / dr3 * (1-kdot);

}
/// end Newtonian cut force (L.Wang)
