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
    assert(ax > 0.0);
#ifdef DEBUG_PRINT_PLANET
    if(ax <= 0.0){
        std::cerr<<"ax="<<ax<<" pos1="<<pos1<<" vel1="<<vel1<<" mass1"<<mass1<<std::endl;
    }
#endif
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
