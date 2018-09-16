#pragma once

#include"matrix3.hpp"

const PS::F64 PI = 4.0*atan(1.0);

// a: semi-major axis
// l: mean anomaly
// e: eccentricity
// u: eccentric anomaly
// n: mean mortion

class Binary{
public:
    // semi: semi-major axis
    // ecc: eccentricity
    // inc: inclination
    // OMG: rotational angle 1
    // omg: rotational angle 2
    // tperi: time to peri-center
    // peri: period
    // ecca: eccentricty anomaly
    // m1: mass 1
    // m2: mass 2
    // tstep: integration step estimation
    // stable_factor: indicate whether the system is stable, >0: switch on slowdown, <=0: no slowdown
    // fpert: perturbation from neighbors (acceleration)
    // am: r \cross v, rotational angular momentum vector (without mass)
    PS::F64 semi, ecc, inc, OMG, omg, tperi, peri, ecca, m1, m2, tstep, stable_factor, fpert;
    PS::F64vec am;

    void dump(FILE *fp) {
        fwrite(this, sizeof(*this),1,fp);
    }

    void read(FILE *fp) {
        size_t rcount = fread(this, sizeof(*this),1,fp);
        if (rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            abort();
        }
    }

    void print(std::ostream& os) {
        os<<"semi: semi-major axis               : "<<semi<<std::endl
          <<"ecc: eccentricity                 : "<<ecc<<std::endl
          <<"inc: inclination                  : "<<inc<<std::endl
          <<"OMG: rotational angle 1           : "<<OMG<<std::endl
          <<"omg: rotational angle 2           : "<<omg<<std::endl
          <<"tperi: time to peri-center        : "<<tperi<<std::endl
          <<"peri: period                      : "<<peri<<std::endl
          <<"ecca: eccentricty anomaly         : "<<ecca<<std::endl
          <<"m1: mass 1                        : "<<m1<<std::endl
          <<"m2: mass 2                        : "<<m2<<std::endl
          <<"tstep: step estimation            : "<<tstep<<std::endl
          <<"stable_factor: stable factor      : "<<stable_factor<<std::endl;
    }
    void print(std::ostream& os, const PS::F64 width=14, const bool title_flag=false) {
        if(title_flag)
            os<<std::setw(width)<<"semi"
              <<std::setw(width)<<"ecc"
              <<std::setw(width)<<"inc"
              <<std::setw(width)<<"plane_rot"
              <<std::setw(width)<<"self_rot"
              <<std::setw(width)<<"tperi"
              <<std::setw(width)<<"period"
              <<std::setw(width)<<"ecca"
              <<std::setw(width)<<"m1"
              <<std::setw(width)<<"m2"
              <<std::setw(width)<<"tstep"
              <<std::setw(width)<<"stable"
              <<std::endl;
        os<<std::setw(width)<<semi
          <<std::setw(width)<<ecc
          <<std::setw(width)<<inc
          <<std::setw(width)<<OMG
          <<std::setw(width)<<omg
          <<std::setw(width)<<tperi
          <<std::setw(width)<<peri
          <<std::setw(width)<<ecca
          <<std::setw(width)<<m1
          <<std::setw(width)<<m2
          <<std::setw(width)<<tstep
          <<std::setw(width)<<stable_factor
          <<std::endl;
    }
};

template <class Tptcl>
class PtclTree: public Binary, public Tptcl {
public:
    Tptcl* member[2];
};

double keplereq(const double l,  const double e, const double u){
    return (u - e*sin(u) - l);
}

double keplereq_dot(const double e, const double u){
    return ( 1.0 - e*cos(u) );
}

// return eccentric anomaly
double solve_keplereq(const double l,
                      const double e){
    double u0 = l;
    double u1;
    int loop = 0;
    while(1){
        loop++;
        double su0 = sin(u0);
        double cu0 = sqrt(1.0 - su0*su0);
        //u1 = u0 - keplereq(l, e, u0)/keplereq_dot(e, u0);
        u1 = u0 - ((u0-e*su0-l)/(1.0 - e*cu0));
        //if( fabs(u1-u0) < 1e-13 ){ return u1; }
        if( fabs(u1-u0) < 1e-15 ){ return u1; }
        else{ u0 = u1; }
        if (loop>1e5) {
            std::cerr<<"Error: kepler solver cannot converge to find correct eccentricity anomaly!\n";
            abort();
        }
    }
}
//double solve_keplereq(const double l,
//                      const double e){
//    double twpi=atan(1.0)*8;
//    double u0 = fmod(l,twpi);
//    //    if(e>0.8) u0 = 3.14159265;
//    double u1;
//    int loop = 0;
//    //    std::cerr<<"l="<<u0<<" e="<<e<<std::endl;
//    while(1){
//        loop++;
//        double su0 = sin(u0);
//        double cu0 = sqrt(1.0 - su0*su0);
//        //u1 = u0 - keplereq(l, e, u0)/keplereq_dot(e, u0);
//        u1 = u0 - ((u0-e*su0-l)/(1.0 - e*cu0));
//        //if( fabs(u1-u0) < 1e-13 ){ return u1; }
//        //        std::cerr<<"u1="<<u1<<"; u0="<<u0<<std::endl;
//        //if( fabs(u1-u0) < 1e-15 ){ return u1; }
//        //if( fabs(fmod(u1,twpi)-fmod(u0,twpi)) < 1e-15 ){ return u1; }
//        double utmp = u1-u0;
//        utmp /= twpi;
//        utmp -= round(utmp);
//        utmp *= twpi;
//        if(utmp < 1e-13){ return u1; }
//        else{ u0 = u1; }
//    }
//}

void OrbParam2PosVel(PS::F64vec & pos0,       PS::F64vec & pos1,
                     PS::F64vec & vel0,       PS::F64vec & vel1,
                     const double mass0, const double mass1,
                     const double semi,    const double ecc,
                     const double inc,   const double OMG,
                     const double omg,   const double u){
    double m_tot = mass0 + mass1;
    double n = sqrt( m_tot / (semi*semi*semi) ); // mean mortion
    double cosu = cos(u);
    double sinu = sin(u);
    double c0 = sqrt(1.0 - ecc*ecc);
    PS::F64vec pos_star(semi*(cosu - ecc), semi*c0*sinu, 0.0);
    PS::F64vec vel_star(-semi*n*sinu/(1.0-ecc*cosu), semi*n*c0*cosu/(1.0-ecc*cosu), 0.0);
    Matrix3<PS::F64> rot;
    rot.rotation(inc, OMG, omg);
    PS::F64vec pos_red = rot*pos_star;
    PS::F64vec vel_red = rot*vel_star;
    pos0 = - mass1 / m_tot * pos_red;
    pos1 =  mass0 / m_tot * pos_red;
    vel0 = - mass1 / m_tot * vel_red;
    vel1 =  mass0 / m_tot * vel_red;

}

//! Orbit to position and velocity
/* @param[out]: _p1: particle 1
   @param[out]: _p2: particle 2
   @param[in]: _bin: binary parameter
 */
template <class Tptcl>
void OrbParam2PosVel(Tptcl& _p1, Tptcl& _p2, const Binary& _bin) {
    double m_tot = _p1.mass + _p2.mass;
    double n = sqrt( m_tot / (_bin.semi*_bin.semi*_bin.semi) );
    double cosu = cos(_bin.ecca);
    double sinu = sin(_bin.ecca);
    double c0 = sqrt(1.0 - _bin.ecc*_bin.ecc);
    PS::F64vec pos_star(_bin.semi*(cosu - _bin.ecc), _bin.semi*c0*sinu, 0.0);
    PS::F64vec vel_star(-_bin.semi*n*sinu/(1.0-_bin.ecc*cosu), _bin.semi*n*c0*cosu/(1.0-_bin.ecc*cosu), 0.0);
    Matrix3<PS::F64> rot;
    rot.rotation(_bin.inc, _bin.OMG, _bin.omg);
    PS::F64vec pos_red = rot*pos_star;
    PS::F64vec vel_red = rot*vel_star;
    _p1.pos = -_p2.mass / m_tot * pos_red;
    _p2.pos =  _p1.mass / m_tot * pos_red;
    _p1.vel = -_p2.mass / m_tot * vel_red;
    _p2.vel =  _p1.mass / m_tot * vel_red;

}

// return eccentric anomayl
// tperi is time to reach pericenter
double PosVel2OrbParam(double & semi,    double & ecc,
                       double & inc,   double & OMG,
                       double & omg,   double & tperi,
                       double & peri,  PS::F64vec & am,
                       const PS::F64vec & pos0, const PS::F64vec & pos1,
                       const PS::F64vec & vel0, const PS::F64vec & vel1,
                       const double mass0, const double mass1){
    double m_tot = mass0 + mass1;
    PS::F64vec pos_red = pos1 - pos0;
    PS::F64vec vel_red = vel1 - vel0;
    double r_sq = pos_red * pos_red;
    double r = sqrt(r_sq);
    double inv_dr = 1.0 / r;
    double v_sq = vel_red * vel_red;
    semi = 1.0 / (2.0*inv_dr - v_sq / m_tot);
    //    assert(semi > 0.0);
    am = pos_red ^ vel_red;
    inc = atan2( sqrt(am.x*am.x+am.y*am.y), am.z);
    OMG = atan2(am.x, -am.y);

    PS::F64vec pos_bar, vel_bar;
    double cosOMG = cos(OMG);
    double sinOMG = sin(OMG);
    double cosinc = cos(inc);
    double sininc = sin(inc);
    pos_bar.x =   pos_red.x*cosOMG + pos_red.y*sinOMG;
    pos_bar.y = (-pos_red.x*sinOMG + pos_red.y*cosOMG)*cosinc + pos_red.z*sininc;
    pos_bar.z = 0.0;
    vel_bar.x =   vel_red.x*cosOMG + vel_red.y*sinOMG;
    vel_bar.y = (-vel_red.x*sinOMG + vel_red.y*cosOMG)*cosinc + vel_red.z*sininc;
    vel_bar.z = 0.0;
    double h = sqrt(am*am);
    double ecccosomg =  h/m_tot*vel_bar.y - pos_bar.x*inv_dr;
    double eccsinomg = -h/m_tot*vel_bar.x - pos_bar.y*inv_dr;
    ecc = sqrt( ecccosomg*ecccosomg + eccsinomg*eccsinomg );
    omg = atan2(eccsinomg, ecccosomg);
    double phi = atan2(pos_bar.y, pos_bar.x); // f + omg (f: true anomaly)
    double f = phi - omg;
    double sinu = r*sin(f) / (semi*sqrt(1.0 - ecc*ecc));
    double cosu = (r*cos(f) / semi) + ecc;
    double u = atan2(sinu, cosu); // eccentric anomaly
    double n = sqrt(m_tot/(semi*semi*semi)); // mean mortion
    peri = 8.0*std::atan(1.0)/n;
    double l = u - ecc*sin(u);  // mean anomaly
    tperi = l / n; 
    return u;
}

//! position velocity to orbit
/* @param[out]: _bin: binary parameter
   @param[in]:  _p1: particle 1
   @param[in]:  _p2: particle 2
 */
template <class Tptcl>
void PosVel2OrbParam(Binary& _bin, const Tptcl& _p1, const Tptcl& _p2){
    double m_tot = _p1.mass + _p2.mass;
    PS::F64vec pos_red = _p2.pos - _p1.pos;
    PS::F64vec vel_red = _p2.vel - _p1.vel;
    double r_sq = pos_red * pos_red;
    double r = sqrt(r_sq);
    double inv_dr = 1.0 / r;
    double v_sq = vel_red * vel_red;
    _bin.semi = 1.0 / (2.0*inv_dr - v_sq / m_tot);
    //    assert(semi > 0.0);
    _bin.am = pos_red ^ vel_red;
    _bin.inc = atan2( sqrt(_bin.am.x*_bin.am.x+_bin.am.y*_bin.am.y), _bin.am.z);
    _bin.OMG = atan2(_bin.am.x, -_bin.am.y);

    PS::F64vec pos_bar, vel_bar;
    double cosOMG = cos(_bin.OMG);
    double sinOMG = sin(_bin.OMG);
    double cosinc = cos(_bin.inc);
    double sininc = sin(_bin.inc);
    pos_bar.x =   pos_red.x*cosOMG + pos_red.y*sinOMG;
    pos_bar.y = (-pos_red.x*sinOMG + pos_red.y*cosOMG)*cosinc + pos_red.z*sininc;
    pos_bar.z = 0.0;
    vel_bar.x =   vel_red.x*cosOMG + vel_red.y*sinOMG;
    vel_bar.y = (-vel_red.x*sinOMG + vel_red.y*cosOMG)*cosinc + vel_red.z*sininc;
    vel_bar.z = 0.0;
    double h = sqrt(_bin.am*_bin.am);
    double ecccosomg =  h/m_tot*vel_bar.y - pos_bar.x*inv_dr;
    double eccsinomg = -h/m_tot*vel_bar.x - pos_bar.y*inv_dr;
    _bin.ecc = sqrt( ecccosomg*ecccosomg + eccsinomg*eccsinomg );
    _bin.omg = atan2(eccsinomg, ecccosomg);
    double phi = atan2(pos_bar.y, pos_bar.x); // f + omg (f: true anomaly)
    double f = phi - _bin.omg;
    double sinu = r*sin(f) / (_bin.semi*sqrt(1.0 - _bin.ecc*_bin.ecc));
    double cosu = (r*cos(f) / _bin.semi) + _bin.ecc;
    _bin.ecca = atan2(sinu, cosu); // eccentric anomaly
    double n = sqrt(m_tot/(_bin.semi*_bin.semi*_bin.semi)); // mean mortion
    _bin.peri = 8.0*std::atan(1.0)/n;
    double l = _bin.ecca - _bin.ecc*sin(_bin.ecca);  // mean anomaly
    _bin.tperi = l / n; 
}



//void PosVel2SemiEccInc(double & semi,    double & ecc,
//                     double & inc,   
//                     const PS::F64vec & pos0, const PS::F64vec & pos1,
//                     const PS::F64vec & vel0, const PS::F64vec & vel1,
//                     const double mass0, const double mass1){
//    double m_tot = mass0 + mass1;
//    PS::F64vec pos_red = pos1 - pos0;
//    PS::F64vec vel_red = vel1 - vel0;
//    double r_sq = pos_red * pos_red;
//    double r = sqrt(r_sq);
//    double inv_dr = 1.0 / r;
//    double v_sq = vel_red * vel_red;
//    semi = 1.0 / (2.0*inv_dr - v_sq / m_tot);
//    PS::F64vec AM = pos_red ^ vel_red;
//    inc = atan2( sqrt(AM.x*AM.x+AM.y*AM.y), AM.z);
//    PS::F64 OMG = atan2(AM.x, -AM.y);
//    PS::F64vec pos_bar, vel_bar;
//    double cosOMG = cos(OMG);
//    double sinOMG = sin(OMG);
//    double cosinc = cos(inc);
//    double sininc = sin(inc);
//    pos_bar.x =   pos_red.x*cosOMG + pos_red.y*sinOMG;
//    pos_bar.y = (-pos_red.x*sinOMG + pos_red.y*cosOMG)*cosinc + pos_red.z*sininc;
//    pos_bar.z = 0.0;
//    vel_bar.x =   vel_red.x*cosOMG + vel_red.y*sinOMG;
//    vel_bar.y = (-vel_red.x*sinOMG + vel_red.y*cosOMG)*cosinc + vel_red.z*sininc;
//    vel_bar.z = 0.0;
//    double h = sqrt(AM*AM);
//    double ecccosomg =  h/m_tot*vel_bar.y - pos_bar.x*inv_dr;
//    double eccsinomg = -h/m_tot*vel_bar.x - pos_bar.y*inv_dr;
//    ecc = sqrt( ecccosomg*ecccosomg + eccsinomg*eccsinomg );
//}
// 
//// return semi major axis
//double PosVel2Ax(const PS::F64vec & pos0, const PS::F64vec & pos1,
//                 const PS::F64vec & vel0, const PS::F64vec & vel1,
//                 const double mass0, const double mass1){
//    //static const PS::F64 PI = 4.0 * atan(1.0);
//    double m_tot = mass0 + mass1;
//    //double m_red = (mass0*mass1)/m_tot;
//    PS::F64vec pos_red = pos1 - pos0;
//    PS::F64vec vel_red = vel1 - vel0;
//    double r_sq = pos_red * pos_red;
//    double r = sqrt(r_sq);
//    double inv_dr = 1.0 / r;
//    double v_sq = vel_red * vel_red;
//    double semi = 1.0 / (2.0*inv_dr - v_sq / m_tot);
//    return semi;
//}

void PosVel2SemiEcc(double & semi,
                  double & ecc,
                  const PS::F64vec & pos0, const PS::F64vec & pos1,
                  const PS::F64vec & vel0, const PS::F64vec & vel1,
                  const double mass0, const double mass1) {
    double m_tot = mass0 + mass1;
    PS::F64vec pos_red = pos1 - pos0;
    PS::F64vec vel_red = vel1 - vel0;
    double r_sq = pos_red * pos_red;
    double r = sqrt(r_sq);
    double inv_dr = 1.0 / r;
    double v_sq = vel_red * vel_red;
    semi = 1.0 / (2.0*inv_dr - v_sq / m_tot);

    double rs = 1.0 - r/semi;
    double rv = pos_red * vel_red;
    ecc = std::sqrt(rs*rs + rv*rv/(m_tot*semi));
}

void DriveKeplerOrbParam(PS::F64vec & pos0,
                 PS::F64vec & pos1,
                 PS::F64vec & vel0,
                 PS::F64vec & vel1,
                 const PS::F64 mass0,
                 const PS::F64 mass1,
                 const PS::F64 dt,
                 const PS::F64 semi,
                 const PS::F64 ecc,
                 const PS::F64 inc,
                 const PS::F64 OMG,
                 const PS::F64 omg,
                 const PS::F64 peri,
                 const PS::F64 ecc_anomaly) {
  PS::F64 freq = sqrt( (mass0+mass1) / (semi*semi*semi) );
  PS::F64 mean_anomaly_old = ecc_anomaly - ecc*sin(ecc_anomaly);
  PS::F64 dt_tmp = dt - (int)(dt/peri)*peri;
  PS::F64 mean_anomaly_new = freq * dt_tmp + mean_anomaly_old; // mean anomaly
  PS::F64 ecc_anomaly_new =  solve_keplereq(mean_anomaly_new, ecc); // eccentric anomaly
  OrbParam2PosVel(pos0, pos1, vel0, vel1, mass0, mass1,
                    semi, ecc, inc, OMG, omg, ecc_anomaly_new);
}

void DriveKepler(const PS::F64 mass0,
                 const PS::F64 mass1,
                 PS::F64vec & pos0,
                 PS::F64vec & pos1,
                 PS::F64vec & vel0,
                 PS::F64vec & vel1,
                 const PS::F64 dt){
    PS::F64 semi, ecc, inc, OMG, omg, tperi, peri, E;
    PS::F64vec am;
    E = PosVel2OrbParam(semi, ecc, inc, OMG, omg, tperi, peri, am,
                            pos0, pos1, vel0, vel1, mass0, mass1);
    DriveKeplerOrbParam(pos0, pos1, vel0, vel1, mass0, mass1, dt, semi, ecc, inc, OMG, omg, peri, E);
}


void DriveKeplerRestricted(PS::F64 mass0,
                           PS::F64vec & pos0,
                           PS::F64vec & pos1,
                           PS::F64vec & vel0,
                           PS::F64vec & vel1,
                           PS::F64 dt){
    DriveKepler(mass0, (PS::F64)0.0, pos0, pos1, vel0, vel1, dt);
}

  //! Calculate the center-of-mass for a group of particles
  /*! Calculate the center-of-mass information for a group of particles. Also the particles can be shifted to their center-of-mass frame.
    The particle class should contain member functions setPos(), setVel(), setMass(), getPos(), getVel(), getMass()
    @param[in] cm: particle type data for storing the center-of-mass information
    @param[in] p: particle list
    @param[in] num: number of particles
    @param[in] fshift: shifting flag to indicate whether particle \a p should be shifted to their center-of-mass. If false (defaulted), no shifting.
   */
template <class particle>
void calc_center_of_mass(particle &cm, particle p[], const int num, const bool fshift=false) {
  PS::F64vec cmr=0;
  PS::F64vec cmv=0;
  PS::F64 cmm = 0;
  for (int i=0;i<num;i++) {
    cmr += p[i].pos * p[i].mass;
    cmv += p[i].vel * p[i].mass;
    cmm += p[i].mass;
  }
  cmr /= cmm; 
  cmv /= cmm; 
      
  cm.mass = cmm;
  cm.pos = cmr;
  cm.vel = cmv;

  // shifting
  if (fshift) {
    for (int i=0;i<num;i++) {
      p[i].pos = p[i].pos - cmr;
      p[i].vel = p[i].vel - cmv;
    }
  }
}

template <class particle>
void calc_center_of_mass(particle &cm, particle* p[], const int num, const bool fshift=false) {
  PS::F64vec cmr=0;
  PS::F64vec cmv=0;
  PS::F64 cmm = 0;
  for (int i=0;i<num;i++) {
    cmr += p[i]->pos * p[i]->mass;
    cmv += p[i]->vel * p[i]->mass;
    cmm += p[i]->mass;
  }
  cmr /= cmm; 
  cmv /= cmm; 
      
  cm.mass = cmm;
  cm.pos = cmr;
  cm.vel = cmv;

  // shifting
  if (fshift) {
    for (int i=0;i<num;i++) {
      p[i]->pos = p[i]->pos - cmr;
      p[i]->vel = p[i]->vel - cmv;
    }
  }
}

// center of mass shift
template <class particle>
void center_of_mass_shift(particle &cm, particle p[], const int num) {
  for (int i=0;i<num;i++) {
    p[i].pos -= cm.pos;
    p[i].vel -= cm.vel;
  }
}

// center of mass shift
template <class particle>
void center_of_mass_shift(particle &cm, particle* p[], const int num) {
  for (int i=0;i<num;i++) {
    p[i]->pos -= cm.pos;
    p[i]->vel -= cm.vel;
  }
}

// correct the particle p position and velocity by adding center-of-mass information
template <class particle>
void center_of_mass_correction(particle &cm, particle p[], const int num) {
  for (int i=0;i<num;i++) {
    p[i].pos += cm.pos;
    p[i].vel += cm.vel;
  }
}

template <class particle>
void center_of_mass_correction(particle &cm, particle* p[], const int num) {
  for (int i=0;i<num;i++) {
    p[i]->pos += cm.pos;
    p[i]->vel += cm.vel;
  }
}

//void getCenterOfMass(Tptcl &cm, Tptcl* p, const PS::S32 adr[], const PS::S32 num) {
//    const PS::S32 ist=adr[0];
//    PS::F64vec cmr=0;
//    PS::F64vec cmv=0;
//    PS::F64 cmm = 0;
//    cm.r_out = p[ist].r_out;
//    cm.id = -p[ist].id;
//    // cm.id_cluster = p[ist].id_cluster;
//    // cm.adr_org = p[ist].adr_org;
//    cm.id_group = -p[ist].id;
//        
//    for (int k=0;k<num;k++) {
//        const PS::S32 i = adr[k];
//        cmr += p[i].pos * p[i].mass;
//        cmv += p[i].vel * p[i].mass;
//        cmm += p[i].mass;
//    }
//    cmr /= cmm; 
//    cmv /= cmm; 
//      
//    cm.mass = cmm;
//    cm.pos = cmr;
//    cm.vel = cmv;
// 
//}

//! Find the minimum distant particle and same the sqrt distance at first and index at second of r2_list
/* @param[in,out] _member_list: group particle member index, will be reordered by the minimum distance chain.
   @param[in,out] _r2_list: a pair. first is the square distance of closes neighbor (next particle i+1 in _member_list); second is equal to particle current index i
   @param[in] _n: number of members
   @param[in] _ptcl_org: original particle data
 */
template <class Tptcl>
void calcMinDisList(PS::S32 _member_list[], 
                    std::pair<PS::F32, PS::S32> _r2_list[],
                    const PS::S32 _n,
                    Tptcl* _ptcl_org) {
    for (int i=0; i<_n-1; i++) {
        PS::S32 k = _member_list[i];
        PS::S32 jc=-1;
        PS::F64 r2min = PS::LARGE_FLOAT;
        for(int j=i+1; j<_n; j++) {
            PS::F64vec dr = _ptcl_org[k].pos - _ptcl_org[_member_list[j]].pos;
            PS::F64 r2 = dr*dr;
            if(r2<r2min) {
                r2min = r2;
                jc = j;
            }
        }
#ifdef HARD_DEBUG
        assert(jc>=0);
#endif
        if (jc!=i+1) {
            PS::S32 jtmp = _member_list[i+1];
            _member_list[i+1] = _member_list[jc];
            _member_list[jc]  = jtmp;
        }
        _r2_list[i].first = r2min;
        _r2_list[i].second = i;
    }
}

bool pairLess(const std::pair<PS::F32,PS::S32> & a, const std::pair<PS::F32,PS::S32> & b)  {
    return a.first < b.first;
}

//! Find systems and generate artificial particles
/*  @param[in,out] _bins: binary tree, size of n_members
    @param[in,out] _member_list: group particle member index, will be reordered by the minimum distance chain.
    @param[in,out] _ptcl_org: particle data, the status of members are modified to 1
    @param[in] _dt_tree: tree time step, used for calculating r_search for bin
 */
template<class Tptcl>
void keplerTreeGenerator(PtclTree<Tptcl> _bins[],   // make sure bins.size = n_members!
                         PS::S32 _member_list[], // make sure list.size = n_members!
                         const PS::S32 _n_members,
                         Tptcl* _ptcl_org,
                         const PS::F64 _dt_tree){

    std::pair<PS::F32,PS::S32> r2_list[_n_members];
    // reorder _member_list by minimum distance of each particles, and save square minimum distance r2min and index i in _member_list (not particle index in _ptcl_org) in r2_list 
    calcMinDisList(_member_list, r2_list, _n_members, _ptcl_org);
    // sort r2_list by r2min 
    if(_n_members>2) std::sort(r2_list, r2_list+_n_members-1, pairLess);

    // tree root for each binary pair
    PtclTree<Tptcl>* bin_host[_n_members]; 
    for(auto &p : bin_host) p=NULL;
    
    // check binary from the closest pair to longest pair
    for(int i=0; i<_n_members-1; i++) {
        // get the pair index k in _member_list with sorted r2_list, i represent the sorted r2min order.
        PS::S32 k = r2_list[i].second;
        Tptcl* p[2];
        // if no tree root assign, set member 1 to particle and their host to current bins i
        if(bin_host[k]==NULL) {
            p[0] = &_ptcl_org[_member_list[k]];
            bin_host[k] = &_bins[i];
            p[0]->status=1; // set status to 1 for counting members
        }
        // if tree root already exist, member 1 assigned to tree root bin_host[k], tree members' bin_host -> current bins i
        else {
            p[0] = bin_host[k];
            PS::S32 ki = k;
            while(bin_host[ki]==p[0]&&ki>=0) bin_host[ki--] = &_bins[i];
        }

        // if no tree root assign, set member 2 to particle and their host to current bins i
        if(bin_host[k+1]==NULL) {
            p[1] = &_ptcl_org[_member_list[k+1]];
            bin_host[k+1] = &_bins[i];
            p[1]->status=1; // set status to 1 for counting members
        }
        // if tree root already exist, member 2 assigned to tree root bin_host[k], tree members' bin_host -> current bins i
        else {
            p[1] = bin_host[k+1];
            PS::S32 ki = k+1;
            while(bin_host[ki]==p[1]&&ki>=0) bin_host[ki++] = &_bins[i];
        }

        // calculate binary parameter
        //_bins[i].ecca=PosVel2OrbParam(_bins[i].semi, _bins[i].ecc, _bins[i].inc, _bins[i].OMG, _bins[i].omg, _bins[i].tperi, _bins[i].peri, _bins[i].axis, p[0]->pos, p[1]->pos, p[0]->vel, p[1]->vel, p[0]->mass, p[1]->mass);
        PosVel2OrbParam(_bins[i], *p[0], *p[1]);
        // calculate center-of-mass particle
        calc_center_of_mass(*(Tptcl*)&_bins[i], p, 2);
        _bins[i].member[0] = p[0];
        _bins[i].member[1] = p[1];
        _bins[i].m1 = p[0]->mass;
        _bins[i].m2 = p[1]->mass;
        PS::F64 dr = 1.0-_bins[i].ecc*cos(_bins[i].ecca);
        PS::F64 apo= 1.0+_bins[i].ecc;
        _bins[i].fpert = fabs(p[0]->mass_bk - p[1]->mass_bk)/dr*apo; // for perturbation, use acceleration difference to estimate
        _bins[i].mass_bk = (p[0]->mass*p[0]->mass_bk + p[1]->mass*p[1]->mass_bk)/_bins[i].mass; // for c.m. perturbation, use m1 a1 + m2 a2 = mcm acm
        _bins[i].id = p[0]->id;
        _bins[i].tstep = -1.0;
        _bins[i].stable_factor = 0.0;
        //_bins[i].status = _bins[i].id;
        _bins[i].status = _bins[i].member[0]->status + _bins[i].member[1]->status;  // counting total number of members in the leafs
        //bins[i].r_search = std::max(p[0]->r_search,p[1]->r_search);
        _bins[i].calcRSearch(_dt_tree);
    }

#ifdef HARD_DEBUG
    for(int i=0; i<_n_members; i++) assert(bin_host[i]==&_bins[_n_members-2]); // check whether all bin_host point to the last of bins
#endif
}

//! Three-body stability function 
/* Use Myllaeri et al. (2018, MNRAS, 476, 830) stability criterion to check whether the system is stable.
   @param[in] _m1: inner binary mass 1
   @param[in] _m1: inner binary mass 2
   @param[in] _semi_in: inner semi-major axis
   @param[in] _ecc_in: inner eccentricity
   @param[in] _mout: outer binary mass 
   @param[in] _ecc_out: outer eccentricity
   @param[in] _pec_out: peri center distance of outer oribit
   @param[in] _period_out: outer binary period
   @param[in] _incline: inclination angle between inner and outer orbit (radians)
   @param[in] _dt: time interval for stable check
   \return stability factor <1 stable; >1 unstable
 */
PS::F64 stab3body(const PS::F64 _m1, 
                  const PS::F64 _m2, 
                  const PS::F64 _semi_in,
                  const PS::F64 _ecc_in, 
                  const PS::F64 _mout, 
                  const PS::F64 _ecc_out, 
                  const PS::F64 _pec_out,
                  const PS::F64 _period_out,
                  const PS::F64 _incline,
                  const PS::F64 _dt) {

    //Adopt 10,000 outer orbits for random walk time-scale.

    PS::F64 fac = 1.0 - 2.0*_ecc_in/3.0*(1.0 - 0.5*std::pow(_ecc_in,2)) 
        - 0.3*std::cos(_incline)*(1.0 - 0.5*_ecc_in + 2.0*std::cos(_incline)*(1.0 - 2.5*std::pow(_ecc_in,1.5) - std::cos(_incline)));

    PS::F64 g = std::sqrt(std::max(_m1,_m2) /(_m1 + _m2))*(1.0 + _mout/(_m1 + _m2));
    
    PS::F64 q = 1.52*std::pow(std::sqrt(std::min(_dt,10000.0)/_period_out)/(1.0 - _ecc_out),1.0/6.0)*std::pow(fac*g,1.0/3.0);

    PS::F64 rp = _pec_out/_semi_in;
    
    PS::F64 stab = q/rp;
    return stab;
    
}


//! check two-body stability
/*            return     tstep                          stable_factor          case
   Unstable:    false     -1                                -1             hyperbolic orbit
   stable:      true      pi/4*sqrt(semi/(m1+m2))*m1*m2   pert/(inner max)   apo-center < r_bin or (apo-center<r_crit & ecc>0.9 & peri-center <r_bin)
                false     -1                                -1             others
                
   @param[in,out] _bin: binary information
   @param[in] _rbin: binary distance criterion
   @param[in] _rcrit: distance criterion
   @param[in] _dt_tree: tree time step
 */
template<class Tptcl>
bool stab2check(PtclTree<Tptcl> &_bin, const PS::F64 _rbin, const PS::F64 _rcrit, const PS::F64 _dt_tree) {

#ifdef STABLE_CHECK_DEBUG
    std::cerr<<"STAB2 semi="<<_bin.semi
             <<" ecc="<<_bin.ecc
             <<" m1="<<_bin.m1
             <<" m2="<<_bin.m2
             <<" period="<<_bin.peri
             <<" apo="<<_bin.semi*(1.0+_bin.ecc)
             <<" pec="<<_bin.semi*(1.0-_bin.ecc)
             <<std::endl;
#endif
    // hyperbolic case
    PS::F64 semi = _bin.semi;
    if(semi<0) {
        _bin.tstep = -1.0;
        _bin.stable_factor = -1;
#ifdef STABLE_CHECK_DEBUG
        std::cerr<<"STAB2 reject: Hyperbolic unstable _bin.semi<0"<<std::endl;
#endif        
        return false;
    }

    // binary case
    PS::F64 apo = _bin.semi*(1.0+_bin.ecc);
    PS::F64 pec = _bin.semi*(1.0-_bin.ecc);

    if (apo<_rbin||(apo<_rcrit&&_bin.ecc>0.9&&pec<_rbin)) {
        PS::F64 m1 = _bin.m1;
        PS::F64 m2 = _bin.m2;
        PS::F64 mcm = m1+m2;
    
//    if(apo>rbin) {
//        PS::F64 pec=bin.semi*(1.0-bin.ecc);
//        if(pec>0.01*rbin||(apo>rmax&&bin.ecc<0.99)) {
//            bin.tstep = 0.78539816339*std::sqrt(bin.semi/(bin.m1+bin.m2))*bin.m1*bin.m2;
//            bin.stable_factor = 1;
//#ifdef STABLE_CHECK_DEBUG
//            std::cerr<<"res=0 apo>rbin("<<rbin<<"), pec>0.01*rbin||apo>rmax("<<rmax<<")"<<std::endl;
//#endif        
//            return 0;
//        }
//    }
    
        //ARC step estimation: pi/4*sqrt(semi/(m1+m2))*m1*m2
        _bin.tstep = 0.78539816339*std::sqrt(semi/mcm)*m1*m2;  

        // perturbation/inner acceleration 
        PS::F64 fpert_ratio = _bin.fpert*(apo*apo)/mcm;
        _bin.stable_factor = fpert_ratio;

        // for almost no perturbation and large period case, stable_factor=-1 to switch off slowdown
        if(fpert_ratio<1e-6) {
            if(_bin.peri>0.125*_dt_tree) _bin.stable_factor = -1;
        }
    
#ifdef STABLE_CHECK_DEBUG
        std::cerr<<"STAB2 accept: Stable pert factor: "<<_bin.stable_factor<<std::endl;
#endif
        
        return true;
    }
    else{
        _bin.tstep = -1.0;
        _bin.stable_factor = -1;
#ifdef STABLE_CHECK_DEBUG
        std::cerr<<"STAB2 reject: Stable but too large orbit, apo: "<<apo<<" ecc: "<<_bin.ecc<<" pec: "<<pec<<std::endl;
#endif 
        return false;
    }
}

//! Three-body (B-S) stability check
/*             return     tstep   stable_factor            case
   Unstable:    false     -1         -1                  hyperbolic outer orbit
                false     -1         -1                  peri-center outer > r_out
                true      inner      -stab3              stab3 >1
    stable:     true      inner      -1                  stab3 <1 & unpert & outer period > 1/8 dt_tree
                false     -1         -1                   --      & acceleration ratio(out/in) <1e-6 and outer period >1e-4 dt_tree
                true      inner      fpert/(outer max)    --      & other cases
   @param[in] _bout: outer orbit parameter
   @param[in] _bin: inner orbit parameter
   @param[in] _rbin: binary radius criterion from input
   @param[in] _rin: inner radius of soft-hard changeover function
   @param[in] _rout: outer radius of soft-hard changeover function
   @param[in] _dt_tree: tree time step
 */
template<class Tptcl>
bool stab3check(PtclTree<Tptcl> &_bout, PtclTree<Tptcl> &_bin, const PS::F64 _rbin, const PS::F64 _rin, const PS::F64 _rout, const PS::F64 _dt_tree) {
#ifdef STABLE_CHECK_DEBUG
    std::cerr<<"STAB3 bout semi="<<_bout.semi
             <<" ecc="<<_bout.ecc
             <<" m1="<<_bout.m1
             <<" m2="<<_bout.m2
             <<" period="<<_bout.peri
             <<" apo="<<_bout.semi*(1.0+_bout.ecc)
             <<" pec="<<_bout.semi*(1.0-_bout.ecc)
             <<" bin semi="<<_bin.semi
             <<" ecc="<<_bin.ecc
             <<" m1="<<_bin.m1
             <<" m2="<<_bin.m2
             <<" period="<<_bin.peri
             <<" apo="<<_bin.semi*(1.0+_bin.ecc)
             <<" pec="<<_bin.semi*(1.0-_bin.ecc)
             <<std::endl;
#endif
    // hyperbolic okuter orbit
    if(_bout.semi<0) {
        _bout.tstep=-1.0;
        _bout.stable_factor=-1;
#ifdef STABLE_CHECK_DEBUG
        std::cerr<<"STAB3 reject: Unstable, Outer body hyperbolic, semi_out: "<<_bout.semi<<std::endl;
#endif
        return false;
    }
    PS::F64 apo_out=_bout.semi*(1.0+_bout.ecc);
    PS::F64 pec_out=_bout.semi*(1.0-_bout.ecc);
    PS::F64 apo_in=_bin.semi*(1.0+_bin.ecc);
    //_bout.tstep = _bin.tstep;

    // too large orbit
    if(pec_out>_rout) {
        _bout.tstep=-1.0;
        _bout.stable_factor=-1;
#ifdef STABLE_CHECK_DEBUG
        std::cerr<<"STAB3 reject: Too large outer orbit, pec_out: "<<pec_out<<std::endl;
#endif
        return false;
    } 

    // too large period
    if(_bout.peri>0.25*_dt_tree) {
        _bout.tstep=-1.0;
        _bout.stable_factor=-1;
#ifdef STABLE_CHECK_DEBUG
        std::cerr<<"STAB3 reject: Too large outer period, period_out: "<<_bout.peri<<std::endl;
#endif
        return false;
    }
        
    // stability check
    // inclination between inner and outer orbit
    PS::F64 incline=std::acos(std::min(1.0, _bout.am*_bin.am/std::sqrt((_bout.am*_bout.am)*(_bin.am*_bin.am))));
    PS::F64 stab3 = stab3body(_bin.m1, _bin.m2, _bin.semi, _bin.ecc, _bout.mass, _bout.ecc, pec_out, _bout.peri, incline, _dt_tree);
    if(stab3>1) {
        // Unstable case
        _bout.stable_factor = -stab3;
#ifdef STABLE_CHECK_DEBUG
        std::cerr<<"STAB3 accept: Unstable, stab3: "<<stab3<<std::endl;
#endif
    }
    else {
        // for large period ratio (acceleration ratio < 1.0e-6), avoid triple system in ARC
        PS::F64 acc_out = _bout.mass/(pec_out*pec_out);
        PS::F64 acc_in  = _bin.mass/(apo_in*apo_in);
        if (acc_out<1.0e-6* acc_in && _bout.peri >1.0e-4*_dt_tree) {
            _bout.tstep=-1.0;
            _bout.stable_factor=-1;
#ifdef STABLE_CHECK_DEBUG
            std::cerr<<"STAB3 reject: Too large period ratio, acc_out/acc_in: "<<acc_out/acc_in<<" period_out: "<<_bout.peri<<std::endl;
#endif
            return false;
        }

        // stable case
        PS::F64 fpert_ratio = _bout.fpert*(apo_out*apo_out)/_bout.mass;
        _bout.stable_factor = fpert_ratio;
        // for unperturbed and large period case, stable_factor=-1 to switch off slowdown
        if(fpert_ratio<1e-6) {
            if(_bout.peri>0.125*_dt_tree) {
                _bout.stable_factor = -1.0;
#ifdef STABLE_CHECK_DEBUG
                std::cerr<<"STAB3 accept: Unstable, large period, period_out: "<<_bout.peri<<std::endl;
#endif
            }
        }
#ifdef STABLE_CHECK_DEBUG
        if(_bout.stable_factor>0) 
            std::cerr<<"STAB3 accept: Stable, fpert_ratio: "<<fpert_ratio<<" stab3: "<<stab3<<" acc_out/acc_in: "<<acc_out/acc_in<<std::endl;
#endif
    }
    _bout.tstep = _bin.tstep;
    return true;

//    if(peri_bout>apo_bin) {
//        PS::F64 dr_p_a = peri_bout*peri_bout-apo_bin*apo_bin;
//        PS::F64 r_pert = 4.0*peri_bout*apo_bin/(dr_p_a*dr_p_a);
//        PS::F64 pert_ratio = _bout.mass/_bin.mass*r_pert*apo_bin*apo_bin;
//        if(pert_ratio>0.01) {
//            _bout.stable_factor=0;
//#ifdef STABLE_CHECK_DEBUG
//            std::cerr<<"True: pert_ratio>0.01: "<<pert_ratio<<std::endl;
//#endif
//            return true;
//        }
//        else _bout.stable_factor=1;
//    }
//    else {
//        _bout.stable_factor=0;
//#ifdef STABLE_CHECK_DEBUG
//            std::cerr<<"True: peri_bout<apo_bin: "<<peri_bout;
//#endif
//        return true;
//    }
//    if(apo_bout>_rin&&(peri_bout>0.01*_rbin||apo_bout>_rout)) {
//        _bout.stable_factor = 1;
//#ifdef STABLE_CHECK_DEBUG
//        std::cerr<<"False: apo_bout>_rin&&(peri_bout>0.01*_rbin||apo_bout>_rout): apo_bout="<<apo_bout<<" peri_bout=<<"<<std::endl;
//#endif
//        return false;
//    }
//#ifdef STABLE_CHECK_DEBUG
//    std::cerr<<"True"<<std::endl;
//#endif
//    return true;
}

//! Four-body (B-B) stability check
/*             return     tstep   stable_factor      case
   Unstable:    false     -1         -1          hyperbolic outer orbit
                false     -1         -1          peri-center outer > r_out
                false     -1         -1          period outer > 0.25 * dt_tree
                true      inner      -stab3_max  stab3_1 >0.8 || stab3_2 > 0.8 & apo_out <= r_out
                false     -1         -1                                        & apo_out >  r_out
    stable:     true      inner      -1          stab3_1 <=0.8 & stab3_2 <=0.8 & unpert & outer period > 1/8 dt_tree
                false     -1         -1                                    --  & acceleration ratio(out/in) <1e-6 and outer period >1e-4 dt_tree
                true      inner      fpert/(m_out/apo_out^2)               --  & other cases
   @param[in] _bout: outer orbit parameter
   @param[in] _bin1: first inner orbit parameter
   @param[in] _bin2: second inner orbit parameter
   @param[in] _rbin: binary radius criterion from input
   @param[in] _rin: inner radius of soft-hard changeover function
   @param[in] _rout: outer radius of soft-hard changeover function
   @param[in] _dt_tree: tree time step
 */
template<class Tptcl>
bool stab4check(PtclTree<Tptcl> &_bout, PtclTree<Tptcl> &_bin1, PtclTree<Tptcl> &_bin2, const PS::F64 _rbin, const PS::F64 _rin, const PS::F64 _rout, const PS::F64 _dt_tree) {
#ifdef STABLE_CHECK_DEBUG
    std::cerr<<"STAB4 bout semi="<<_bout.semi
             <<" ecc="<<_bout.ecc
             <<" m1="<<_bout.m1
             <<" m2="<<_bout.m2
             <<" period="<<_bout.peri
             <<" apo="<<_bout.semi*(1.0+_bout.ecc)
             <<" pec="<<_bout.semi*(1.0-_bout.ecc)
             <<" bin1 semi="<<_bin1.semi
             <<" ecc="<<_bin1.ecc
             <<" m1="<<_bin1.m1
             <<" m2="<<_bin1.m2
             <<" period="<<_bin1.peri
             <<" apo="<<_bin1.semi*(1.0+_bin1.ecc)
             <<" pec="<<_bin1.semi*(1.0-_bin1.ecc)
             <<" bin2 semi="<<_bin2.semi
             <<" ecc="<<_bin2.ecc
             <<" m1="<<_bin2.m1
             <<" m2="<<_bin2.m2
             <<" period="<<_bin2.peri
             <<" apo="<<_bin2.semi*(1.0+_bin2.ecc)
             <<" pec="<<_bin2.semi*(1.0-_bin2.ecc)
             <<std::endl;
#endif
    // hyperbolic outer orbit
    if(_bout.semi<0) {
        _bout.tstep=-1.0;
        _bout.stable_factor=-1;
#ifdef STABLE_CHECK_DEBUG
        std::cerr<<"STAB4 reject: Unstable, Outer body hyperbolic, semi_out: "<<_bout.semi<<std::endl;
#endif
        return false;
    }

    PS::F64 apo_out=_bout.semi*(1.0+_bout.ecc);
    PS::F64 pec_out=_bout.semi*(1.0-_bout.ecc);
    PS::F64 apo_in1=_bin1.semi*(1.0+_bin1.ecc);
    PS::F64 apo_in2=_bin2.semi*(1.0+_bin2.ecc);
    // the seperation is large for perturbation estimation
    //_bout.tstep=std::min(_bin1.tstep,_bin2.tstep);
    
    // too large orbit
    if(pec_out>_rout) {
        _bout.tstep=-1.0;
        _bout.stable_factor=-1;
#ifdef STABLE_CHECK_DEBUG
        std::cerr<<"STAB4 reject: Unstable, Too large outer orbit, pec_out: "<<pec_out<<std::endl;
#endif
        return false;
    } 

    // too large period
    if(_bout.peri>0.25*_dt_tree) {
        _bout.tstep=-1.0;
        _bout.stable_factor=-1;
#ifdef STABLE_CHECK_DEBUG
        std::cerr<<"STAB4 reject: Too large outer period, period_out: "<<_bout.peri<<std::endl;
#endif
        return false;
    }
    
    // stability check
    // inclination between inner and outer orbit
    PS::F64 incl1=std::acos(std::min(1.0, _bout.am*_bin1.am/std::sqrt((_bout.am*_bout.am)*(_bin1.am*_bin1.am))));
    PS::F64 incl2=std::acos(std::min(1.0, _bout.am*_bin2.am/std::sqrt((_bout.am*_bout.am)*(_bin2.am*_bin2.am))));
    PS::F64 stab3_1 = stab3body(_bin1.m1, _bin1.m2, _bin1.semi, _bin1.ecc, _bout.mass, _bout.ecc, pec_out, _bout.peri, incl1, _dt_tree);
    PS::F64 stab3_2 = stab3body(_bin2.m1, _bin2.m2, _bin2.semi, _bin2.ecc, _bout.mass, _bout.ecc, pec_out, _bout.peri, incl2, _dt_tree);
    if(stab3_1>0.8||stab3_2>0.8) {
        // Unstable case
        if (apo_out>_rout) {
            _bout.tstep=-1.0;
            _bout.stable_factor=-1;
#ifdef STABLE_CHECK_DEBUG
            std::cerr<<"STAB4 reject: Unstable, too large outer orbit, stab3_1: "<<stab3_1<<" stab3_2: "<<stab3_2<<std::endl;
#endif
            return false;
        }
        _bout.stable_factor = -std::max(stab3_1,stab3_2);
#ifdef STABLE_CHECK_DEBUG
        std::cerr<<"STAB4 accept: Unstable, stab3_1: "<<stab3_1<<" stab3_2: "<<stab3_2<<std::endl;
#endif
    }
    else {
        // for large period ratio (acceleration ratio < 1.0e-6), avoid quad system in ARC
        PS::F64 acc_out = _bout.mass/(pec_out*pec_out);
        PS::F64 acc_in1  = _bin1.mass/(apo_in1*apo_in1);
        PS::F64 acc_in2  = _bin2.mass/(apo_in2*apo_in2);
        PS::F64 acc_in_max = std::max(acc_in1, acc_in2);
        if (acc_out<1.0e-6* acc_in_max && _bout.peri >1.0e-4*_dt_tree) {
            _bout.tstep=-1.0;
            _bout.stable_factor=-1;
#ifdef STABLE_CHECK_DEBUG
            std::cerr<<"STAB4 reject: Too large period ratio, acc_out/acc_in: "<<acc_out/acc_in_max<<" period_out: "<<_bout.peri<<std::endl;
#endif
            return false;
        }

        // stable case
        PS::F64 fpert_ratio = _bout.fpert*(apo_out*apo_out)/_bout.mass;
        _bout.stable_factor = fpert_ratio;
        // for unperturbed and large period case, stable_factor=-1 to switch off slowdown
        if(fpert_ratio<1e-6) {
            if(_bout.peri>0.125*_dt_tree) {
                _bout.stable_factor = -1.0;
#ifdef STABLE_CHECK_DEBUG
                std::cerr<<"STAB4 accept: Unstable, large period, period_out: "<<_bout.peri<<std::endl;
#endif
            }
        }
#ifdef STABLE_CHECK_DEBUG
        if(_bout.stable_factor>0) 
            std::cerr<<"STAB4 accept: Stable, fpert_ratio: "<<fpert_ratio<<" stab3_1: "<<stab3_1<<" stab3_2: "<<stab3_2<<" acc_out/acc_in: "<<acc_out/acc_in_max<<std::endl;
#endif
    }
    _bout.tstep = std::min(_bin1.tstep,_bin2.tstep);
    return true;

//    if(peri_bout>apo_bin1+apo_bin2) {
//        /* Assume all particles are along one line
//           peri_bout -> p3
//           apo_bin1/2 -> a1/2
//           Then strongest perurbation distance is p3-(a1+a2)
//           Weakest p3+(a1+a2)
//           The perturbation is estimated as
//           M[1/2] * 4 p3*(a1+a2)/(p3^2-(a1+a2)^2)
//           The inner force
//           M1/a1^2 | M2/a2^2
//         */
//        PS::F64 apo12 = apo_bin1 + apo_bin2;
//        PS::F64 dr_p_a = peri_bout*peri_bout-apo12*apo12;
//        PS::F64 r_pert = 4.0*peri_bout*apo12/(dr_p_a*dr_p_a);
//        PS::F64 pert_ratio = std::max(bin1.mass/bin2.mass*r_pert*apo_bin2*apo_bin2,
//                                      bin2.mass/bin1.mass*r_pert*apo_bin1*apo_bin1);
//        if(pert_ratio>0.01) {
//            bout.stable_factor=0;
//            return true;
//        }
//        else bout.stable_factor=1;
//    }
//    else{
//        // Unstable
//        bout.stable_factor=0;
//        return true;
//    }
//    if(apo_bout>rin&&(peri_bout>0.01*rbin||apo_bout>rout)) {
//        bout.stable_factor = 1;
//        return false;
//    }
//    return true;
}

//! Stability check for hierarchtical tree
/* Check stability of each level in _bins, save stable system in _stab_bins. The tstep and stable_factor of _stab_bins are updated
   @param[out] _stab_bins: stable system array
   @param[in]  _bins: hierarchtical tree
   @param[in]  _rbin: binary detection criterion radius
   @param[in]  _rin: inner radius of soft-hard changeover function
   @param[in]  _rout: outer radius of soft-hard changeover function
   @param[in]  _dt_tree: tree time step for calculating r_search
 */
template<class Tptcl>
bool stabilityCheck(PS::ReallocatableArray<PtclTree<Tptcl>*> &_stab_bins, 
                    PtclTree<Tptcl> &_bins, 
                    const PS::F64 _rbin, 
                    const PS::F64 _rin, 
                    const PS::F64 _rout, 
                    const PS::F64 _dt_tree) {
    if(_bins.member[0]->status>1) {
        // B-B system
        if(_bins.member[1]->status>1) {
            //PS::F64 fs0 = stab4check(_bins, *(PtclTree<Tptcl>*)_bins.member[0], *(PtclTree<Tptcl>*)_bins.member[1], rbin, rin);
            bool fs1 = stabilityCheck<Tptcl>(_stab_bins, *(PtclTree<Tptcl>*)_bins.member[0], _rbin, _rin, _rout, _dt_tree);
            bool fs2 = stabilityCheck<Tptcl>(_stab_bins, *(PtclTree<Tptcl>*)_bins.member[1], _rbin, _rin, _rout, _dt_tree);
            if(fs1&&fs2) {
                bool fs0 = stab4check(_bins, *(PtclTree<Tptcl>*)_bins.member[0], *(PtclTree<Tptcl>*)_bins.member[1], _rbin, _rin, _rout, _dt_tree);
                if(fs0) return true;
            }
            if(fs1) _stab_bins.push_back((PtclTree<Tptcl>*)_bins.member[0]);
            if(fs2) _stab_bins.push_back((PtclTree<Tptcl>*)_bins.member[1]);
            return false;
        }
        else { // B-S system
            //PS::F64 fs0 = stab3check(_bins, *(PtclTree<Tptcl>*)_bins.member[0], _rbin, _rin);
            bool fs1 = stabilityCheck<Tptcl>(_stab_bins, *(PtclTree<Tptcl>*)_bins.member[0], _rbin, _rin, _rout, _dt_tree);
            if(fs1) {
                bool fs0 = stab3check(_bins, *(PtclTree<Tptcl>*)_bins.member[0], _rbin, _rin, _rout, _dt_tree);
                if(fs0) return true;
            }
            if(fs1) _stab_bins.push_back((PtclTree<Tptcl>*)_bins.member[0]);
            return false;
        }
    }
    else {
        if(_bins.member[1]->status>1) { // S-B system
            //PS::F64 fs0 = stab3check(_bins, *(PtclTree<Tptcl>*)_bins.member[1], _rbin, _rin);
            bool fs1 = stabilityCheck<Tptcl>(_stab_bins, *(PtclTree<Tptcl>*)_bins.member[1], _rbin, _rin, _rout, _dt_tree);
            if(fs1) {
                bool fs0 = stab3check(_bins, *(PtclTree<Tptcl>*)_bins.member[1], _rbin, _rin, _rout, _dt_tree);
                if(fs0) return true;
            }
            if(fs1) _stab_bins.push_back((PtclTree<Tptcl>*)_bins.member[1]);
            return false;
        }
        else { // Binary
            bool fs = stab2check(_bins, _rbin, _rout, _dt_tree);
            return fs;
        }
    }
}

