#pragma once

//inline PS::F64 CalcK(const PS::F64 rij,
//                     const PS::F64 rout,
//                     const PS::F64 rin){
//    PS::F64 inv_dr = 1.0 / (rout-rin);
//    PS::F64 x = (rij - rin)*inv_dr;
//    x = (x < 1.0) ? x : 1.0;
//    x = (x > 0.0) ? x : 0.0;
//    PS::F64 x2 = x*x;
//    PS::F64 x4 = x2*x2;
//    PS::F64 k = (((-20.0*x+70.0)*x-84.0)*x+35.0)*x4;
//    std::max( std::min(k, 1.0), 0.0);
//    return k;
//}


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


//inline PS::F64 cutoff_poly_3rd(const PS::F64 rij,
//                               const PS::F64 rout,
//                               const PS::F64 rin,
//                               const PS::F64 inv_dr){
//    PS::F64 x = (rij - rin)*inv_dr;
//    x = (x < 1.0) ? x : 1.0;
//    x = (x > 0.0) ? x : 0.0;
//    PS::F64 x2 = x*x;
//    PS::F64 x4 = x2*x2;
//    PS::F64 k = (((-20.0*x+70.0)*x-84.0)*x+35.0)*x4;
//    return k;
//}

/// start cutoff dr (L.Wang)
//inline PS::F64 cutoff_poly_3rd_dr(const PS::F64 rij,
//                                  const PS::F64 rx,
//                                  const PS::F64 rout,
//                                  const PS::F64 rin){
//    PS::F64 inv_dr = 1.0 / (rout-rin);
//    PS::F64 x = (rij - rin)*inv_dr;
//    x = (x < 1.0) ? x : 1.0;
//    x = (x > 0.0) ? x : 0.0;
//    PS::F64 x2 = x*x;
//    PS::F64 x3 = x2*x;
//    PS::F64 k = (((-140.0*x+420.0)*x-420.0)*x+140.0)*x3*inv_dr*rx/rij;
//    return k;
//}
/// end cutoff dr (L.Wang)


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
