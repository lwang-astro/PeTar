#pragma once

#include"matrix3.hpp"

// a: semi-major axis
// l: mean anomaly
// e: eccentricity
// u: eccentric anomaly
// n: mean mortion

double keplereq(const double l,  const double e, const double u){
    return (u - e*sin(u) - l);
}

double keplereq_dot(const double e, const double u){
    return ( 1.0 - e*cos(u) );
}

// return eccentric anomaly
double solve_keplereq(const double l,
                      const double e){
    double twpi=atan(1.0)*8;
    double u0 = fmod(l,twpi);
    //    if(e>0.8) u0 = 3.14159265;
    double u1;
    int loop = 0;
    //    std::cerr<<"l="<<u0<<" e="<<e<<std::endl;
    while(1){
        loop++;
        double su0 = sin(u0);
        double cu0 = sqrt(1.0 - su0*su0);
        //u1 = u0 - keplereq(l, e, u0)/keplereq_dot(e, u0);
        u1 = u0 - ((u0-e*su0-l)/(1.0 - e*cu0));
        //if( fabs(u1-u0) < 1e-13 ){ return u1; }
        //        std::cerr<<"u1="<<u1<<"; u0="<<u0<<std::endl;
        //if( fabs(u1-u0) < 1e-15 ){ return u1; }
        //if( fabs(fmod(u1,twpi)-fmod(u0,twpi)) < 1e-15 ){ return u1; }
        double utmp = u1-u0;
        utmp /= twpi;
        utmp -= round(utmp);
        utmp *= twpi;
        if(utmp < 1e-13){ return u1; }
        else{ u0 = u1; }
    }
}

void OrbParam2PosVel(PS::F64vec & pos0,       PS::F64vec & pos1,
                     PS::F64vec & vel0,       PS::F64vec & vel1,
                     const double mass0, const double mass1,
                     const double ax,    const double ecc,
                     const double inc,   const double OMG,
                     const double omg,   const double u){
    double m_tot = mass0 + mass1;
    double n = sqrt( m_tot / (ax*ax*ax) );
    double cosu = cos(u);
    double sinu = sin(u);
    double c0 = sqrt(1.0 - ecc*ecc);
    PS::F64vec pos_star(ax*(cosu - ecc), ax*c0*sinu, 0.0);
    PS::F64vec vel_star(-ax*n*sinu/(1.0-ecc*cosu), ax*n*c0*cosu/(1.0-ecc*cosu), 0.0);
    Matrix3<PS::F64> rot;
    rot.rotation(inc, OMG, omg);
    PS::F64vec pos_red = rot*pos_star;
    PS::F64vec vel_red = rot*vel_star;
    pos0 = - mass1 / m_tot * pos_red;
    pos1 =  mass0 / m_tot * pos_red;
    vel0 = - mass1 / m_tot * vel_red;
    vel1 =  mass0 / m_tot * vel_red;

}

// return eccentric anomayl
// tperi is not just tperi
double PosVel2OrbParam(double & ax,    double & ecc,
                       double & inc,   double & OMG,
                       double & omg,   double & tperi,
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
    ax = 1.0 / (2.0*inv_dr - v_sq / m_tot);
    //    assert(ax > 0.0);
    PS::F64vec AM = pos_red ^ vel_red;
    inc = atan2( sqrt(AM.x*AM.x+AM.y*AM.y), AM.z);
    OMG = atan2(AM.x, -AM.y);

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
    double h = sqrt(AM*AM);
    double ecccosomg =  h/m_tot*vel_bar.y - pos_bar.x*inv_dr;
    double eccsinomg = -h/m_tot*vel_bar.x - pos_bar.y*inv_dr;
    ecc = sqrt( ecccosomg*ecccosomg + eccsinomg*eccsinomg );
    omg = atan2(eccsinomg, ecccosomg);
    double phi = atan2(pos_bar.y, pos_bar.x); // f + omg (f: true anomaly)
    double f = phi - omg;
    double sinu = r*sin(f) / (ax*sqrt(1.0 - ecc*ecc));
    double cosu = (r*cos(f) / ax) + ecc;
    double u = atan2(sinu, cosu); // eccentric anomaly
    double n = sqrt(m_tot/ax*ax*ax); // mean mortion
    double l = u - ecc*sin(u);
    tperi = l / n;
    return u;
}


double PosVel2OrbParamDebug(double & ax,    double & ecc,
			    double & inc,   double & OMG,
			    double & omg,   double & tperi,
			    const PS::F64vec & pos0, const PS::F64vec & pos1,
			    const PS::F64vec & vel0, const PS::F64vec & vel1,
			    const double mass0, const double mass1,
			    std::ofstream & fout){
    double m_tot = mass0 + mass1;
    PS::F64vec pos_red = pos1 - pos0;
    PS::F64vec vel_red = vel1 - vel0;
    double r_sq = pos_red * pos_red;
    double r = sqrt(r_sq);
    double inv_dr = 1.0 / r;
    double v_sq = vel_red * vel_red;
    fout<<"2.0*inv_dr= "<<2.0*inv_dr<<" v_sq / m_tot= "<<v_sq / m_tot<<std::endl;
    ax = 1.0 / (2.0*inv_dr - v_sq / m_tot);
    fout<<"ax= "<<ax<<std::endl;
    PS::F64vec AM = pos_red ^ vel_red;
    fout<<"AM= "<<AM<<std::endl;
    inc = atan2( sqrt(AM.x*AM.x+AM.y*AM.y), AM.z);
    fout<<"inc= "<<inc<<std::endl;
    OMG = atan2(AM.x, -AM.y);
    fout<<"OMG= "<<OMG<<std::endl;

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
    double h = sqrt(AM*AM);
    double ecccosomg =  h/m_tot*vel_bar.y - pos_bar.x*inv_dr;
    double eccsinomg = -h/m_tot*vel_bar.x - pos_bar.y*inv_dr;
    ecc = sqrt( ecccosomg*ecccosomg + eccsinomg*eccsinomg );
    fout<<"ecc= "<<ecc<<std::endl;
    omg = atan2(eccsinomg, ecccosomg);
    fout<<"omg= "<<omg<<std::endl;
    double phi = atan2(pos_bar.y, pos_bar.x); // f + omg (f: true anomaly)
    double f = phi - omg;
    double sinu = r*sin(f) / (ax*sqrt(1.0 - ecc*ecc));
    double cosu = (r*cos(f) / ax) + ecc;
    double u = atan2(sinu, cosu); // eccentric anomaly
    double n = sqrt(m_tot/ax*ax*ax); // mean mortion
    double l = u - ecc*sin(u);
    tperi = l / n;
    fout<<"tper= "<<tperi<<std::endl;
    return u;
}




void PosVel2AxEccInc(double & ax,    double & ecc,
                     double & inc,   
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
    ax = 1.0 / (2.0*inv_dr - v_sq / m_tot);
    PS::F64vec AM = pos_red ^ vel_red;
    inc = atan2( sqrt(AM.x*AM.x+AM.y*AM.y), AM.z);
    PS::F64 OMG = atan2(AM.x, -AM.y);
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
    double h = sqrt(AM*AM);
    double ecccosomg =  h/m_tot*vel_bar.y - pos_bar.x*inv_dr;
    double eccsinomg = -h/m_tot*vel_bar.x - pos_bar.y*inv_dr;
    ecc = sqrt( ecccosomg*ecccosomg + eccsinomg*eccsinomg );
}

// return semi major axis
double PosVel2Ax(const PS::F64vec & pos0, const PS::F64vec & pos1,
                 const PS::F64vec & vel0, const PS::F64vec & vel1,
                 const double mass0, const double mass1){
    //static const PS::F64 PI = 4.0 * atan(1.0);
    double m_tot = mass0 + mass1;
    //double m_red = (mass0*mass1)/m_tot;
    PS::F64vec pos_red = pos1 - pos0;
    PS::F64vec vel_red = vel1 - vel0;
    double r_sq = pos_red * pos_red;
    double r = sqrt(r_sq);
    double inv_dr = 1.0 / r;
    double v_sq = vel_red * vel_red;
    double ax = 1.0 / (2.0*inv_dr - v_sq / m_tot);
    return ax;
}


/*
void DriveKepler(const PS::F64 mass0,
                 const PS::F64 mass1,
                 PS::F64vec & pos0,
                 PS::F64vec & pos1,
                 PS::F64vec & vel0,
                 PS::F64vec & vel1,
                 const PS::F64 dt){
    //static __thread PS::F64 ax, ecc, inc, OMG, omg, tperi, u_old, n, l_old, l_new, u_new;
    PS::F64 ax, ecc, inc, OMG, omg, tperi, u_old, n, l_old, l_new, u_new;
    u_old = PosVel2OrbParam(ax, ecc, inc, OMG, omg, tperi,
				    pos0, pos1, vel0, vel1, mass0, mass1);
    n = sqrt( (mass0+mass1) / (ax*ax*ax) );
    l_old = u_old - ecc*sin(u_old);
    l_new = n * dt + l_old; // mean anomaly
    u_new = solve_keplereq(l_new, ecc); // eccentric anomaly

    OrbParam2PosVel(pos0, pos1, vel0, vel1, mass0, mass1,
		    ax, ecc, inc, OMG, omg, u_new);
}
*/

void DriveKepler(const PS::F64 mass0,
                 const PS::F64 mass1,
                 PS::F64vec & pos0,
                 PS::F64vec & pos1,
                 PS::F64vec & vel0,
                 PS::F64vec & vel1,
                 const PS::F64 dt){
    PS::F64 ax, ecc, inc, OMG, omg, tperi, u_old, n, l_old, l_new, u_new;
    u_old = PosVel2OrbParam(ax, ecc, inc, OMG, omg, tperi,
                            pos0, pos1, vel0, vel1, mass0, mass1);
    n = sqrt( (mass0+mass1) / (ax*ax*ax) );
    l_old = u_old - ecc*sin(u_old);
    l_new = n * dt + l_old; // mean anomaly
    u_new = solve_keplereq(l_new, ecc); // eccentric anomaly
    OrbParam2PosVel(pos0, pos1, vel0, vel1, mass0, mass1,
                    ax, ecc, inc, OMG, omg, u_new);
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
  double cmr[3]={};
  double cmv[3]={};
  double cmm = 0;
  for (std::size_t i=0;i<num;i++) {
    const double *ri = p[i].getPos();
    const double *vi = p[i].getVel();
    const double mi = p[i].getMass();
    cmr[0] += ri[0] * mi;
    cmr[1] += ri[1] * mi;
    cmr[2] += ri[2] * mi;

    cmv[0] += vi[0] * mi;
    cmv[1] += vi[1] * mi;
    cmv[2] += vi[2] * mi;

    cmm += mi;
  }
  cmr[0] /= cmm; 
  cmr[1] /= cmm; 
  cmr[2] /= cmm; 

  cmv[0] /= cmm; 
  cmv[1] /= cmm; 
  cmv[2] /= cmm;
      
  cm.setMass(cmm);
  cm.setPos(cmr[0],cmr[1],cmr[2]);
  cm.setVel(cmv[0],cmv[1],cmv[2]);

  // shifting
  if (fshift) {
    for (std::size_t i=0;i<num;i++) {
      const double *ri = p[i].getPos();
      const double *vi = p[i].getVel();
      p[i].setPos(ri[0] - cmr[0],
                  ri[1] - cmr[1],
                  ri[2] - cmr[2]);
      p[i].setVel(vi[0] - cmv[0],
                  vi[1] - cmv[1],
                  vi[2] - cmv[2]);
    }
  }
}

// correct the particle p position and velocity by adding center-of-mass information
template <class particle>
void center_of_mass_correction(particle &cm, particle p[], const int num) {
  const double *rc = cm.getPos();
  const double *vc = cm.getVel();
  for (std::size_t i=0;i<num;i++) {
    const double *ri = p[i].getPos();
    p[i].setPos(ri[0] + rc[0],
                ri[1] + rc[1],
                ri[2] + rc[2]);
    const double *vi = p[i].getVel();
    p[i].setVel(vi[0] + vc[0],
                vi[1] + vc[1],
                vi[2] + vc[2]);
  }
}
