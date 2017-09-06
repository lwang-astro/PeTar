#pragma once

#ifdef INTEGRATED_CUTOFF_FUNCTION
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

inline PS::F64 cutoff_poly_3rd(const PS::F64 rij,
                               const PS::F64 r_oi_inv,   // r_oi_inv = 1/(r_out-r_in)
                               const PS::F64 r_A,        // r_A = (r_out-r_in)/(r_out-r_in)
                               const PS::F64 rin){
    //PS::F64 inv_dr = 1.0 / (rout-rin);
    PS::F64 x = (rij - rin)*r_oi_inv;
    x = (x < 1.0) ? x : 1.0;
    x = (x > 0.0) ? x : 0.0;
    PS::F64 x2 = x*x;
    PS::F64 x4 = x2*x2;
    PS::F64 k = 1-(((-20.0*x+70.0)*x-84.0)*x+35.0)*x4;
    return k;
}

inline PS::F64 cutoff_poly_3rd_dot(const PS::F64 rij,
                                   const PS::F64 rijvij,
                                   const PS::F64 r_oi_inv,
                                   const PS::F64 r_A,        // r_A = (r_out-r_in)/(r_out-r_in)
                                   const PS::F64 rin){
    //PS::F64 rout = _rout;
    //PS::F64 rin = _rin;
    //PS::F64 inv_dr = 1.0/(rout-rin);
    PS::F64 x = (rij - rin)*r_oi_inv;
    PS::F64 xdot = rijvij/rij*r_oi_inv;
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
        Kdot = -(-140.0*x6 + 420.0*x5 - 420.0*x4 + 140.0*x3) * xdot;
    }
    return Kdot;
}

#else

inline PS::F64 cutoff_poly_3rd(const PS::F64 rij,
                               const PS::F64 r_oi_inv,   // r_oi_inv = 1/(r_out-r_in)
                               const PS::F64 r_A,        // r_A = (r_out-r_in)/(r_out-r_in)
                               const PS::F64 rin){
    //PS::F64 inv_dr = 1.0 / (rout-rin);
    PS::F64 x = (rij - rin)*r_oi_inv;
    x = (x < 1.0) ? x : 1.0;
    x = (x > 0.0) ? x : 0.0;
    PS::F64 x_1 = x - 1;
    PS::F64 x_2 = x_1*x_1;
    PS::F64 x_4 = x_2*x_2;
    PS::F64 x2 = x*x;
    PS::F64 x3 = x2*x;
    PS::F64 x4 = x2*x2;
    PS::F64 k = x_4*(1.0 + 4.0*x + 10.0*x2 + 20.0*x3 + 35.0*r_A*x4);
    return k;
}

inline PS::F64 cutoff_poly_3rd_dot(const PS::F64 rij,
                                   const PS::F64 rijvij,
                                   const PS::F64 r_oi_inv,
                                   const PS::F64 r_A,        // r_A = (r_out-r_in)/(r_out-r_in)
                                   const PS::F64 rin){
    //PS::F64 inv_dr = 1.0/(rout-rin);
    PS::F64 x = (rij - rin)*r_oi_inv;
    PS::F64 xdot = rijvij/rij*r_oi_inv;
    PS::F64 Kdot = 0.0;
    if(x > 0.0 && x < 1.0) {
        PS::F64 x3 = x*x*x;
        PS::F64 x_1 = x - 1;
        PS::F64 x_3 = x_1*x_1*x_1;
        Kdot = r_A*280.0*x3*(rin*r_oi_inv + x)*x_3*xdot;
    }
    return Kdot;
}

inline PS::F64 cutoff_pot(const PS::F64 rij, 
                          const PS::F64 r_oi_inv, 
                          const PS::F64 r_A, 
                          const PS::F64 r_in){    // rA = (r_out-r_in)/(r_out+r_in)
    PS::F64 x = (rij - r_in)*r_oi_inv;
    PS::F64 k = 1.0;
    if(x >= 1.0 ) k += r_A;
    else if(x > 0.0) {
        PS::F64 x2 = x*x;
        PS::F64 x3 = x2*x;
        PS::F64 x5 = x2*x3;
        k -= r_A*x5*(5.0*x3 - 20.0*x2 + 28.0*x - 14.0);
    }
    return k;
}

#endif
