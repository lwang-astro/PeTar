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
    // ax: semi-major axis
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
    // stable_factor: indicate whether the system is stable
    // fpert: perturbation from neighbors
    PS::F64 ax, ecc, inc, OMG, omg, tperi, peri, ecca, m1, m2, tstep, stable_factor, fpert;

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

    void print(std::ostream& os, const PS::F64 width) {
        os<<"ax: semi-major axis               : "<<ax<<std::endl
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
          <<"stable_factor: step estimation    : "<<stable_factor<<std::endl;
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
// tperi is time to reach pericenter
double PosVel2OrbParam(double & ax,    double & ecc,
                       double & inc,   double & OMG,
                       double & omg,   double & tperi,
                       double & peri,
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
    double n = sqrt(m_tot/(ax*ax*ax)); // mean mortion
    peri = 8.0*std::atan(1.0)/n;
    double l = u - ecc*sin(u);  // mean anomaly
    tperi = l / n; 
    return u;
}



//void PosVel2AxEccInc(double & ax,    double & ecc,
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
//    ax = 1.0 / (2.0*inv_dr - v_sq / m_tot);
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
//    double ax = 1.0 / (2.0*inv_dr - v_sq / m_tot);
//    return ax;
//}

void PosVel2AxEcc(double & ax,
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
    ax = 1.0 / (2.0*inv_dr - v_sq / m_tot);

    double rs = 1.0 - r/ax;
    double rv = pos_red * vel_red;
    ecc = std::sqrt(rs*rs + rv*rv/(m_tot*ax));
}

void DriveKeplerOrbParam(PS::F64vec & pos0,
                 PS::F64vec & pos1,
                 PS::F64vec & vel0,
                 PS::F64vec & vel1,
                 const PS::F64 mass0,
                 const PS::F64 mass1,
                 const PS::F64 dt,
                 const PS::F64 ax,
                 const PS::F64 ecc,
                 const PS::F64 inc,
                 const PS::F64 OMG,
                 const PS::F64 omg,
                 const PS::F64 peri,
                 const PS::F64 ecc_anomaly) {
  PS::F64 freq = sqrt( (mass0+mass1) / (ax*ax*ax) );
  PS::F64 mean_anomaly_old = ecc_anomaly - ecc*sin(ecc_anomaly);
  PS::F64 dt_tmp = dt - (int)(dt/peri)*peri;
  PS::F64 mean_anomaly_new = freq * dt_tmp + mean_anomaly_old; // mean anomaly
  PS::F64 ecc_anomaly_new =  solve_keplereq(mean_anomaly_new, ecc); // eccentric anomaly
  OrbParam2PosVel(pos0, pos1, vel0, vel1, mass0, mass1,
                    ax, ecc, inc, OMG, omg, ecc_anomaly_new);
}

void DriveKepler(const PS::F64 mass0,
                 const PS::F64 mass1,
                 PS::F64vec & pos0,
                 PS::F64vec & pos1,
                 PS::F64vec & vel0,
                 PS::F64vec & vel1,
                 const PS::F64 dt){
    PS::F64 ax, ecc, inc, OMG, omg, tperi, peri, E;
    E = PosVel2OrbParam(ax, ecc, inc, OMG, omg, tperi, peri,
                            pos0, pos1, vel0, vel1, mass0, mass1);
    DriveKeplerOrbParam(pos0, pos1, vel0, vel1, mass0, mass1, dt, ax, ecc, inc, OMG, omg, peri, E);
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

template <class Tptcl>
void calcMinDisList(PS::S32 member_list[], 
                    std::pair<PS::F32, PS::S32> r2_list[],
                    const PS::S32 n,
                    Tptcl* ptcl_org) {
    for (int i=0; i<n-1; i++) {
        PS::S32 k = member_list[i];
        PS::S32 jc=-1;
        PS::F64 r2min = PS::LARGE_FLOAT;
        for(int j=i+1; j<n; j++) {
            PS::F64vec dr = ptcl_org[k].pos - ptcl_org[member_list[j]].pos;
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
            PS::S32 jtmp = member_list[i+1];
            member_list[i+1] = member_list[jc];
            member_list[jc]  = jtmp;
        }
        r2_list[i].first = r2min;
        r2_list[i].second = i;
    }
}

bool pairLess(const std::pair<PS::F32,PS::S32> & a, const std::pair<PS::F32,PS::S32> & b)  {
    return a.first < b.first;
}

template<class Tptcl>
void keplerTreeGenerator(PtclTree<Tptcl> bins[],   // make sure bins.size = n_members!
                         PS::S32 member_list[], // make sure list.size = n_members!
                         const PS::S32 n_members,
                         Tptcl* ptcl_org,
                         const PS::F64 dt_tree){

    std::pair<PS::F32,PS::S32> r2_list[n_members];
    calcMinDisList(member_list,r2_list, n_members, ptcl_org);
    if(n_members>2) std::sort(r2_list, r2_list+n_members-1, pairLess);

    PtclTree<Tptcl>* bin_host[n_members];
    for(auto &p : bin_host) p=NULL;
    
    for(int i=0; i<n_members-1; i++) {
        PS::S32 k = r2_list[i].second;
        Tptcl* p[2];
        if(bin_host[k]==NULL) {
            p[0] = &ptcl_org[member_list[k]];
            bin_host[k] = &bins[i];
        }
        else {
            p[0] = bin_host[k];
            PS::S32 ki = k;
            while(bin_host[ki]==p[0]&&ki>=0) bin_host[ki--] = &bins[i];
        }

        if(bin_host[k+1]==NULL) {
            p[1] = &ptcl_org[member_list[k+1]];
            bin_host[k+1] = &bins[i];
        }
        else {
            p[1] = bin_host[k+1];
            PS::S32 ki = k+1;
            while(bin_host[ki]==p[1]&&ki>=0) bin_host[ki++] = &bins[i];
        }
        
        bins[i].ecca=PosVel2OrbParam(bins[i].ax, bins[i].ecc, bins[i].inc, bins[i].OMG, bins[i].omg, bins[i].tperi, bins[i].peri, p[0]->pos, p[1]->pos, p[0]->vel, p[1]->vel, p[0]->mass, p[1]->mass);
        calc_center_of_mass(*(Tptcl*)&bins[i], p, 2);
        bins[i].member[0] = p[0];
        bins[i].member[1] = p[1];
        bins[i].m1 = p[0]->mass;
        bins[i].m2 = p[1]->mass;
        PS::F64 dr = 1.0-bins[i].ecc*cos(bins[i].ecca);
        PS::F64 apo= 1.0+bins[i].ecc;
        bins[i].fpert = fabs(p[0]->mass_bk - p[1]->mass_bk)/dr*apo;
        bins[i].mass_bk = p[0]->mass_bk + p[1]->mass_bk;
        bins[i].id = p[0]->id;
        bins[i].tstep = -1.0;
        bins[i].stable_factor = 0.0;
        bins[i].status = bins[i].id;
        //bins[i].r_search = std::max(p[0]->r_search,p[1]->r_search);
        bins[i].calcRSearch(dt_tree);
    }

#ifdef HARD_DEBUG
    for(int i=0; i<n_members; i++) assert(bin_host[i]==&bins[n_members-2]);
#endif
}

//! return integration step estimation if good
/* bin.tstep is //pi/4*sqrt(ax/(m1+m2))*m1*m2;
   bin.stable_factor = 1
   if ax<0, bin.tstep = -1.0, bin.stable_factor = 0
   return 0: unstable; 1: closed orbit; 2: closed orbit but strong perturbation
 */
template<class Tptcl>
PS::S32 stab2check(PtclTree<Tptcl> &bin, const PS::F64 rbin, const PS::F64 rmax) {
#ifdef STABLE_CHECK_DEBUG
    std::cerr<<"STAB2 ax="<<bin.ax<<" ecc="<<bin.ecc<<" m1="<<bin.m1<<" m2="<<bin.m2<<" period="<<bin.peri
             <<" apo="<<bin.ax*(1.0+bin.ecc)<<" peri="<<bin.ax*(1.0-bin.ecc)
             <<std::endl;
#endif
    if(bin.ax<0) {
        bin.tstep = -1.0;
        bin.stable_factor = 0;
#ifdef STABLE_CHECK_DEBUG
        std::cerr<<"res=0 bin.ax<0"<<std::endl;
#endif        
        return 0;
    }
    PS::F64 apo = bin.ax*(1.0+bin.ecc);
    if(apo>rbin) {
        PS::F64 pec=bin.ax*(1.0-bin.ecc);
        if(pec>0.01*rbin||apo>rmax) {
            bin.tstep = 0.78539816339*std::sqrt(bin.ax/(bin.m1+bin.m2))*bin.m1*bin.m2;
            bin.stable_factor = 1;
            return 0;
#ifdef STABLE_CHECK_DEBUG
            std::cerr<<"res=0 apo>rbin<0"<<std::endl;
#endif        
        }
    }
    
    //PS::F64 mrate = (bin.m1>bin.m2)?bin.m2/bin.m1:bin.m1/bin.m2;
    //return bin.peri*std::sqrt(mrate)*pow(1.0-bin.ecc,0.41666666);
    //  return std::max(std::sqrt(std::abs(1.0-bin.ecc)),0.01)*peri;
    bin.tstep=0.78539816339*std::sqrt(bin.ax/(bin.m1+bin.m2))*bin.m1*bin.m2;  //pi/4*sqrt(ax/(m1+m2))*m1*m2
    bin.stable_factor=1.0;

    //avoid too strong perturbation
    if(bin.fpert>0.01*std::max(bin.m1,bin.m2)/(apo*apo)) return 2;
    
#ifdef STABLE_CHECK_DEBUG
    std::cerr<<"res=1"<<std::endl;
#endif
        
    return 1;
}

template<class Tptcl>
bool stab3check(PtclTree<Tptcl> &bout, PtclTree<Tptcl> &bin, const PS::F64 rbin, const PS::F64 rin, const PS::F64 rout) {
#ifdef STABLE_CHECK_DEBUG
    std::cerr<<"STAB3 bout ax="<<bout.ax<<" ecc="<<bout.ecc<<" m1="<<bout.m1<<" m2="<<bout.m2<<" period="<<bout.peri
             <<" apo="<<bout.ax*(1.0+bout.ecc)<<" peri="<<bout.ax*(1.0-bout.ecc)
             <<std::endl;
#endif
    if(bout.ax<0) {
        bout.tstep=-1.0;
        bout.stable_factor=0;
#ifdef STABLE_CHECK_DEBUG
        std::cerr<<"False: bout.ax<0: "<<bout.ax<<std::endl;
#endif
        return false;
    }
    PS::F64 apo_bout=bout.ax*(1.0+bout.ecc);
    PS::F64 peri_bout=bout.ax*(1.0-bout.ecc);
    PS::F64 apo_bin=bin.ax*(1.0+bin.ecc);
    bout.tstep = bin.tstep;
    if(peri_bout>apo_bin) {
        PS::F64 dr_p_a = peri_bout*peri_bout-apo_bin*apo_bin;
        PS::F64 r_pert = 4.0*peri_bout*apo_bin/(dr_p_a*dr_p_a);
        PS::F64 pert_ratio = bout.mass/bin.mass*r_pert*apo_bin*apo_bin;
        if(pert_ratio>0.01) {
            bout.stable_factor=0;
#ifdef STABLE_CHECK_DEBUG
            std::cerr<<"True: pert_ratio>0.01: "<<pert_ratio<<std::endl;
#endif
            return true;
        }
        else bout.stable_factor=1;
    }
    else {
        bout.stable_factor=0;
#ifdef STABLE_CHECK_DEBUG
            std::cerr<<"True: peri_bout<apo_bin: "<<peri_bout;
#endif
        return true;
    }
    if(apo_bout>rin&&(peri_bout>0.01*rbin||apo_bout>rout)) {
        bout.stable_factor = 1;
#ifdef STABLE_CHECK_DEBUG
        std::cerr<<"False: apo_bout>rin&&(peri_bout>0.01*rbin||apo_bout>rout): apo_bout="<<apo_bout<<" peri_bout=<<"<<std::endl;
#endif
        return false;
    }
#ifdef STABLE_CHECK_DEBUG
    std::cerr<<"True"<<std::endl;
#endif
    return true;
}

template<class Tptcl>
bool stab4check(PtclTree<Tptcl> &bout, PtclTree<Tptcl> &bin1, PtclTree<Tptcl> &bin2, const PS::F64 rbin, const PS::F64 rin, const PS::F64 rout) {
#ifdef STABLE_CHECK_DEBUG
    std::cerr<<"STAB4 bout ax="<<bout.ax<<" ecc="<<bout.ecc<<" m1="<<bout.m1<<" m2="<<bout.m2<<" period="<<bout.peri
             <<" apo="<<bout.ax*(1.0+bout.ecc)<<" peri="<<bout.ax*(1.0-bout.ecc)
             <<std::endl;
#endif
    if(bout.ax<0) {
        bout.tstep=-1.0;
        bout.stable_factor=0;
        return false;
    }
    //PS::F64 apo_bout=bout.ax*(1.0+bout.ecc);
    //if(apo_bout>rout) return -apo_bout;
    PS::F64 apo_bout=bout.ax*(1.0+bout.ecc);
    PS::F64 peri_bout=bout.ax*(1.0-bout.ecc);
    PS::F64 apo_bin1=bin1.ax*(1.0+bin1.ecc);
    PS::F64 apo_bin2=bin2.ax*(1.0+bin2.ecc);
    // the seperation is large for perturbation estimation
    bout.tstep=std::min(bin1.tstep,bin2.tstep);
    if(peri_bout>apo_bin1+apo_bin2) {
        /* Assume all particles are along one line
           peri_bout -> p3
           apo_bin1/2 -> a1/2
           Then strongest perurbation distance is p3-(a1+a2)
           Weakest p3+(a1+a2)
           The perturbation is estimated as
           M[1/2] * 4 p3*(a1+a2)/(p3^2-(a1+a2)^2)
           The inner force
           M1/a1^2 | M2/a2^2
         */
        PS::F64 apo12 = apo_bin1 + apo_bin2;
        PS::F64 dr_p_a = peri_bout*peri_bout-apo12*apo12;
        PS::F64 r_pert = 4.0*peri_bout*apo12/(dr_p_a*dr_p_a);
        PS::F64 pert_ratio = std::max(bin1.mass/bin2.mass*r_pert*apo_bin2*apo_bin2,
                                      bin2.mass/bin1.mass*r_pert*apo_bin1*apo_bin1);
        if(pert_ratio>0.01) {
            bout.stable_factor=0;
            return true;
        }
        else bout.stable_factor=1;
    }
    else{
        // Unstable
        bout.stable_factor=0;
        return true;
    }
    if(apo_bout>rin&&(peri_bout>0.01*rbin||apo_bout>rout)) {
        bout.stable_factor = 1;
        return false;
    }
    return true;
}


template<class Tptcl>
PS::S32 stabilityCheck(PS::ReallocatableArray<PtclTree<Tptcl>*> &nbin, 
                       PtclTree<Tptcl> &bins, const PS::F64 rbin, const PS::F64 rin, const PS::F64 rout) {
    if(bins.member[0]->status!=0) {
        if(bins.member[1]->status!=0) {
            //PS::F64 fs0 = stab4check(bins, *(PtclTree<Tptcl>*)bins.member[0], *(PtclTree<Tptcl>*)bins.member[1], rbin, rin);
            PS::S32 fs1 = stabilityCheck<Tptcl>(nbin, *(PtclTree<Tptcl>*)bins.member[0], rbin, rin, rout);
            PS::S32 fs2 = stabilityCheck<Tptcl>(nbin, *(PtclTree<Tptcl>*)bins.member[1], rbin, rin, rout);
            if(fs1&&fs2) {
                bool fs0 = stab4check(bins, *(PtclTree<Tptcl>*)bins.member[0], *(PtclTree<Tptcl>*)bins.member[1], rbin, rin, rout);
                if(fs0) return 1;
            }
            if(fs1==1) nbin.push_back((PtclTree<Tptcl>*)bins.member[0]);
            if(fs2==1) nbin.push_back((PtclTree<Tptcl>*)bins.member[1]);
            return 0;
        }
        else {
            //PS::F64 fs0 = stab3check(bins, *(PtclTree<Tptcl>*)bins.member[0], rbin, rin);
            PS::S32 fs1 = stabilityCheck<Tptcl>(nbin, *(PtclTree<Tptcl>*)bins.member[0], rbin, rin, rout);
            if(fs1) {
                bool fs0 = stab3check(bins, *(PtclTree<Tptcl>*)bins.member[0], rbin, rin, rout);
                if(fs0) return 1;
            }
            if(fs1==1) nbin.push_back((PtclTree<Tptcl>*)bins.member[0]);
            return 0;
        }
    }
    else {
        if(bins.member[1]->status!=0) {
            //PS::F64 fs0 = stab3check(bins, *(PtclTree<Tptcl>*)bins.member[1], rbin, rin);
            PS::S32 fs1 = stabilityCheck<Tptcl>(nbin, *(PtclTree<Tptcl>*)bins.member[1], rbin, rin, rout);
            if(fs1) {
                bool fs0 = stab3check(bins, *(PtclTree<Tptcl>*)bins.member[1], rbin, rin, rout);
                if(fs0) return 1;
            }
            if(fs1==1) nbin.push_back((PtclTree<Tptcl>*)bins.member[1]);
            return 0;
        }
        else {
            return stab2check(bins, rbin, rout);
        }
    }
}

