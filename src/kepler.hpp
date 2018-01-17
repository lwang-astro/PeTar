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
    PS::F64 ax, ecc, inc, OMG, omg, tperi, peri, ecca, m1, m2, tstep;
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
        bins[i].id = p[0]->id;
        bins[i].tstep = 0.0; 
        bins[i].status = bins[i].id;
        //bins[i].r_search = std::max(p[0]->r_search,p[1]->r_search);
        bins[i].calcRSearch(dt_tree);
    }

#ifdef HARD_DEBUG
    for(int i=0; i<n_members; i++) assert(bin_host[i]==&bins[n_members-2]);
#endif
}

// return integration step estimation if good
template<class Tptcl>
PS::F64 stab2check(const PtclTree<Tptcl> &bin, const PS::F64 rmax) {
#ifdef STABLE_CHECK_DEBUG
    std::cerr<<"STAB2 ax="<<bin.ax<<" ecc="<<bin.ecc<<" m1="<<bin.m1<<" m2="<<bin.m2<<" peri="<<bin.peri<<" res="
             <<(bin.ax<0 ? 0.0 : (rmax<bin.ax*(1.0+bin.ecc) ? -(bin.ax*(1.0+bin.ecc)) : bin.peri*std::sqrt((bin.m1>bin.m2)?bin.m2/bin.m1:bin.m1/bin.m2)*pow(1.0-bin.ecc,0.41666666)))
             <<std::endl;
#endif
    if(bin.ax<0) return 0.0;
    PS::F64 apo = bin.ax*(1.0+bin.ecc);
    if(apo>rmax) return -apo;
    PS::F64 mrate = (bin.m1>bin.m2)?bin.m2/bin.m1:bin.m1/bin.m2;
    return bin.peri*std::sqrt(mrate)*pow(1.0-bin.ecc,0.41666666);
//    return std::max(std::sqrt(std::abs(1.0-bin.ecc)),0.01)*peri;
}

template<class Tptcl>
PS::F64 stab3check(const PtclTree<Tptcl> &bout, const PtclTree<Tptcl> &bin, const PS::F64 rbin, const PS::F64 rin) {
    PS::F64 s2=stab2check(bin, rbin);
    if (s2<=0.0) return s2;
    PS::F64 s1=stab2check(bout, rin);
    if (s1<=0.0) return s1;
    return std::min(s1,s2);
}

template<class Tptcl>
PS::F64 stab4check(const PtclTree<Tptcl> &bout, const PtclTree<Tptcl> &bin1, const PtclTree<Tptcl> &bin2, const PS::F64 rbin, const PS::F64 rin) {
    PS::F64 s2=stab2check(bin1,rbin);
    if (s2<=0.0) return s2;
    PS::F64 s3=stab2check(bin2,rbin);
    if (s3<=0.0) return s3;
    PS::F64 s1=stab2check(bout,rin);
    if (s1<=0.0) return s1;
    return std::min(s2,s3);
}


template<class Tptcl>
PS::F64 stabilityCheck(PS::ReallocatableArray<PtclTree<Tptcl>*> &nbin, 
                       PtclTree<Tptcl> &bins, const PS::F64 rbin, const PS::F64 rin) {
    PS::F64 fstab=0.0;
    if(bins.member[0]->status!=0) {
        if(bins.member[1]->status!=0) {
            //PS::F64 fs0 = stab4check(bins, *(PtclTree<Tptcl>*)bins.member[0], *(PtclTree<Tptcl>*)bins.member[1], rbin, rin);
            PS::F64 fs1 = stabilityCheck<Tptcl>(nbin, *(PtclTree<Tptcl>*)bins.member[0], rbin, rin);
            PS::F64 fs2 = stabilityCheck<Tptcl>(nbin, *(PtclTree<Tptcl>*)bins.member[1], rbin, rin);
            PS::F64 fs0 = stab2check(bins, rin);
            if(fs0<=0) {
                if(fs1>0) {
                    nbin.push_back((PtclTree<Tptcl>*)bins.member[0]);
                    nbin.back()->tstep = fs1;
                }
                if(fs2>0) {
                    nbin.push_back((PtclTree<Tptcl>*)bins.member[1]);
                    nbin.back()->tstep = fs2;
                }
                fstab = -1.0;
            }
            else {
                fstab = std::min(fs1,fs2);
                if(fstab>0) bins.tstep = fstab;
            }
        }
        else {
            //PS::F64 fs0 = stab3check(bins, *(PtclTree<Tptcl>*)bins.member[0], rbin, rin);
            PS::F64 fs1 = stabilityCheck<Tptcl>(nbin, *(PtclTree<Tptcl>*)bins.member[0], rbin, rin);
            PS::F64 fs0 = stab2check(bins, rin);
            if(fs0<=0) {
                if(fs1>0) {
                    nbin.push_back((PtclTree<Tptcl>*)bins.member[0]);
                    nbin.back()->tstep = fs1;
                }
                fstab = -1.0;
            }
            else {
                fstab = fs1;
                if(fstab>0) bins.tstep = fstab;
            }

        }
    }
    else {
        if(bins.member[1]->status!=0) {
            //PS::F64 fs0 = stab3check(bins, *(PtclTree<Tptcl>*)bins.member[1], rbin, rin);
            PS::F64 fs1 = stabilityCheck<Tptcl>(nbin, *(PtclTree<Tptcl>*)bins.member[1], rbin, rin);
            PS::F64 fs0 = stab2check(bins, rin);
            if(fs0<=0) {
                if(fs1>0) {
                    nbin.push_back((PtclTree<Tptcl>*)bins.member[1]);
                    nbin.back()->tstep = fs1;
                }
                fstab = -1.0;
            }
            else {
                fstab = fs1;
                if (fstab>0) bins.tstep = fstab;
            }
        }
        else {
            fstab = stab2check(bins, rbin);
            if(fstab>0) bins.tstep = fstab;
        }
    }
    
    return fstab;
}

