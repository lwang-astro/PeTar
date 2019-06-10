#pragma once
#include "kepler.hpp"

//!initializaton of system parameters
/*! Obtain Radius parameters, consistent with the input help information
  @param[in]     _tsys:   particle system (soft)
  @param[in,out] _r_in:   changeover function inner boundary
  @param[in,out] _r_out:  changeover function outer boundary
  @param[in,out] _r_bin:  arc group radius criterion
  @param[in,out] _r_search_min: minimum searching radius
  @param[in,out] _r_search_max: maximum searching radius
  @param[out] _v_max:  maximum velocity to calculate r_search
  @param[out] _m_average: averaged mass of particles
  @param[out] _m_max: maximum mass of particles
  @param[in,out] _dt_tree: tree time step
  @param[out] _vel_disp: system velocity dispersion 
  @param[in]     _search_factor: coefficient to calculate r_search
  @param[in]     _ratio_r_cut: _r_out/_r_in
  @param[in]     _n_bin: number of binaries
 */
template<class Tpsys>
void GetInitPar(const Tpsys & _tsys,
                PS::F64 &_r_in,
                PS::F64 &_r_out,
                PS::F64 &_r_bin,
                PS::F64 &_r_search_min,
                PS::F64 &_r_search_max,
                PS::F64 &_v_max,
                PS::F64 &_m_average,
                PS::F64 &_m_max,
                PS::F64 &_dt_tree,
                PS::F64 &_vel_disp,
                const PS::F64 _search_factor,
                const PS::F64 _ratio_r_cut,
                const PS::S64 _n_bin,
                const PS::F64 _theta) {

    // local particle number
    const PS::S64 n_loc = _tsys.getNumberOfParticleLocal();

    // local c.m velocity
    PS::F64vec vel_cm_loc = 0.0;
    // local c.m. mass
    PS::F64 mass_cm_loc = 0.0;
    // local maximum mass
    PS::F64 mass_max_loc = 0.0;

    for(PS::S64 i=0; i<n_loc; i++){
        PS::F64 mi = _tsys[i].mass;
        PS::F64vec vi = _tsys[i].vel;

#ifdef MAIN_DEBUG
        assert(mi>0);
#endif
        mass_cm_loc += mi;
        vel_cm_loc += mi * vi;
        mass_max_loc = std::max(mi, mass_max_loc);
    }

    // global c.m. parameters
    PS::F64    mass_cm_glb = PS::Comm::getSum(mass_cm_loc);
    _m_max = PS::Comm::getMaxValue(mass_max_loc);
    PS::F64vec vel_cm_glb  = PS::Comm::getSum(vel_cm_loc);
    vel_cm_glb /= mass_cm_glb;

    // local velocity square
    PS::F64 vel_sq_loc = 0.0;
    PS::S64 n_vel_loc_count = 0;

    // single particle starting index
    PS::S64 single_start_index = 0;
    const PS::S64 bin_last_id = 2*_n_bin;
    if (_tsys[0].id<bin_last_id) {
        single_start_index = std::min(bin_last_id - _tsys[0].id + 1,n_loc);
        if(single_start_index%2!=0) single_start_index--;
    }
    // binary particle starting index
    const PS::S64 binary_start_index = (_tsys[0].id%2==0)?1:0;

    // calculate velocity dispersion
    for (PS::S64 i=binary_start_index; i<single_start_index; i+=2) {
        PS::F64 m1 = _tsys[i].mass;
        PS::F64 m2 = _tsys[i+1].mass;
        PS::F64vec dv = (m1*_tsys[i].vel + m2*_tsys[i+1].vel)/(m1+m2) - vel_cm_glb;
        vel_sq_loc += dv * dv;
        n_vel_loc_count++;
    }
    
    for (PS::S64 i=single_start_index; i<n_loc; i++){
        PS::F64vec dv = _tsys[i].vel - vel_cm_glb;
        vel_sq_loc += dv * dv;
        n_vel_loc_count++;
    }

    const PS::S64    n_vel_glb_count= PS::Comm::getSum(n_vel_loc_count);
    const PS::S64    n_glb          = PS::Comm::getSum(n_loc);
    const PS::F64    vel_sq_glb     = PS::Comm::getSum(vel_sq_loc);
    _vel_disp   = sqrt(vel_sq_glb / 3.0 / (PS::F64)n_vel_glb_count);

    PS::F64 average_mass_glb = mass_cm_glb/(PS::F64)n_glb;
    _m_average = average_mass_glb;

    // flag to check whether r_ous is already defined
    bool r_out_flag = (_r_out>0);
    
    // if r_out is already defined, calculate r_in based on _ratio_r_cut
    if (r_out_flag) _r_in = _r_out * _ratio_r_cut;
    // calculate r_in based on velocity dispersion and averaged mass, calculate r_out by _ratio_r_cut
    else {
        _r_in = 0.5*average_mass_glb / (_vel_disp*_vel_disp);
        _r_out = _r_in / _ratio_r_cut;
    }

    // if tree time step is not defined, calculate tree time step by r_out and velocity dispersion
    if (_dt_tree==0.0) {
        PS::F64 dt_origin = 0.1*_r_out / _vel_disp;
        _dt_tree = 1.0;
        if (dt_origin<1) while (_dt_tree>dt_origin) _dt_tree *= 0.5;
        else {
            while (_dt_tree<=dt_origin) _dt_tree *= 2.0;
            _dt_tree *= 0.5;
        }
        //// if r_out is not defined, calculate r_out based on tree step and velocity dispersion
        //if (!r_out_flag) {
        //    _r_out = 10.0*_dt_tree*_vel_disp;
        //    _r_in = _r_out*_ratio_r_cut;
        //}
    }

    // if r_bin is not defined, set to theta * r_in
    if (_r_bin==0.0) _r_bin = _theta*_r_in;

    // if r_search_min is not defined, calculate by search_factor*velocity_dispersion*tree_time_step + r_out
    if (_r_search_min==0.0) _r_search_min = _search_factor*_vel_disp*_dt_tree + _r_out;
    // if r_search_max is not defined, calcualte by 5*r_out
    if (_r_search_max==0.0) _r_search_max = 5*_r_out;
    // calculate v_max based on r_search_max, tree time step and search_factor
    _v_max = (_r_search_max - _r_out) / _dt_tree / _search_factor;
}

// f(x) = x/sigma^2*exp(-x^2/(2sigma^2))
// F(x) = 1.0 - exp(-x^2/(2sigma^2))
// x = sqrt(-2*sigma^2*ln(1-F))
// 2sigma^2 = <e^2> or <i^2>
// <> means R.M.S.
double RayleighDistribution(const double sigma){
    static PS::MTTS mt;
    static bool first = true;
    if(first){
        mt.init_genrand( PS::Comm::getRank() );
        first = false;
    }
    double F = mt.genrand_res53();
    /*
    double ret = 0.0;
    do{
	ret = sqrt( -2.0*sigma*sigma*log(1.0-F));
    }while(ret >= 1.0);
    return ret;    
    */
    return sqrt( -2.0*sigma*sigma*log(1.0-F));
}

// p(a)da = C*a^-0.5*da
// P(a) = 2*C*(a^0.5 - a_in^0.5)
// a = (P/(2C) + a_in^0.5)^2
double HayashiDistribution(const double a_in, const double a_out){
    static PS::MTTS mt;
    static bool first = true;
    if(first){
	mt.init_genrand( PS::Comm::getRank() );
	first = false;
    }
    const double C = 0.5 / (sqrt(a_out) - sqrt(a_in));
    double P = mt.genrand_res53();
    double ret = P/(2*C) + sqrt(a_in);
    ret *= ret;
    return ret;
}


// p(a)da = C*a^(p+1)*da
// P(a) = 2*C*(a^0.5 - a_in^0.5)
double HayashiDistributionWithIceLine(const double a_in, const double a_out,
				      const double a_ice, const double f_ice=1.0, const double p_sigma=-1.5){
    //std::cout<<"a_in="<<a_in<<std::endl;
    //const PS::F64 p_sigma = -1.5;
    const PS::F64 p_mass = 2.0 + p_sigma;
    //std::cout<<"p_mass="<<p_mass<<std::endl;
    static PS::MTTS mt;
    //std::cout<<"check 0"<<std::endl;
    static bool first = true;
    //std::cout<<"check 1"<<std::endl;
    //std::cout<<"PS::Comm::getRank()="<<PS::Comm::getRank()<<std::endl;
    if(first){
	mt.init_genrand( PS::Comm::getRank() );
	first = false;
    }
    //std::cout<<"first="<<first<<std::endl;
    double P = mt.genrand_res53();
    double ret = 0.0;
    double C = 0.0;
    double P_ice = 0.0;
    if(a_ice <= a_in || a_ice >= a_out){
	//std::cout<<"check a"<<std::endl;
	C = p_mass / (pow(a_out, p_mass) - pow(a_in, p_mass));
	ret = pow(P*p_mass/C + pow(a_in, p_mass), 1.0/p_mass);
    }
    else{
	//std::cout<<"check b"<<std::endl;
	C = p_mass / ( f_ice*(pow(a_out, p_mass)-pow(a_ice, p_mass)) + pow(a_ice, p_mass) - pow(a_in, p_mass) );
	P_ice = C/p_mass*(pow(a_ice, p_mass) - pow(a_in, p_mass));
	//std::cout<<"C="<<C<<std::endl;
	//std::cout<<"P_ice="<<P_ice<<std::endl;
	if(P < P_ice){
	    ret = pow(P*p_mass/C + pow(a_in, p_mass), 1.0/p_mass);
	}
	else{
	    ret = pow( ((P*p_mass/C + pow(a_in, p_mass) - pow(a_ice, p_mass))/f_ice + pow(a_ice, p_mass)), 1.0/p_mass); // beyond the ice line
	}
    }
    return ret;
}


double AU2CM(const double x){
    return x * 1.49597871e13;
}

double CM2AU(const double x){
    return x / 1.49597871e13;
}

// [M] = 1 [g]
// [L] = 1 [cm]
// [T] = 3.871e3[sec]
// [V] = 2.5833118e-4 [cm/sec] = 2.5833118e-9 [km/sec]
// 1yr = 8.146732e3[T]
// or
// G = 1 [L]^3/([M][T]^2)
// [L] = 1[AU] = 1.495978707e8[km]
// [M] = 1[Msun] = 1.989e30[kg]
// [T] = 5.02198050479e6 [sec] = 0.15924595715 [yr]
// [V] = 29.7886203575 [km/sec]
void MakeKeplerDisk(PS::F64 & mass_planet_glb,
                    PS::F64 *& mass,
                    PS::F64vec *& pos,
                    PS::F64vec *& vel,
                    const long long int n_glb,
                    const long long int n_loc,
                    const double a_in, // [AU]
                    const double a_out, // [AU]
                    const double e_rms, // normalized
                    const double i_rms, // normalized
                    const double dens = 10.0, // [g/cm^2]
                    const double mass_sun = 1.0, //[m_sun]
		    const double a_ice = 0.0,
		    const double f_ice = 1.0,
		    const double power = -1.5,
		    const int seed = 0
    ){
    static const double mass_sun_gram = 1.989e33; //[g]
    PS::MTTS mt;
    //mt.init_genrand( PS::Comm::getRank() );
    mt.init_genrand( PS::Comm::getRank()+seed*PS::Comm::getNumberOfProc() );
    static const double PI = atan(1.0) * 4.0;
    const double AU = AU2CM(1.0);
    //mass_planet_glb = 4.0 * PI * dens * ( sqrt(a_out) - sqrt(a_in) ) * AU * AU / mass_sun_gram; // [Msun]
    mass_planet_glb = 2.0 * PI * dens / (2.0+power) * ( pow(a_out, 2.0+power) - pow(a_in, 2.0+power) ) * AU * AU / mass_sun_gram; // [Msun]
    const double m_planet = mass_planet_glb / n_glb;
    mass = new double[n_loc];
    pos = new PS::F64vec[n_loc];
    vel = new PS::F64vec[n_loc];
    const double h = pow(2.0*m_planet / (3.0*mass_sun), 1.0/3.0);
    PS::F64 e_ave = 0.0;
    const PS::F64 e_sigma = sqrt(0.5*e_rms*e_rms); // this is right procedure
    const PS::F64 i_sigma = sqrt(0.5*i_rms*i_rms);
    //const PS::F64 e_sigma = e_rms;
    //const PS::F64 i_sigma = i_rms;
    for(long long int i=0; i<n_loc; i++){
        mass[i] = m_planet;
	double ax = HayashiDistributionWithIceLine(a_in, a_out, a_ice, f_ice, power);
        double ecc = RayleighDistribution(e_sigma) * h;
        double inc = RayleighDistribution(i_sigma) * h;
        PS::F64vec pos_dummy, vel_dummy;
        double omg = 2.0 * PI * mt.genrand_res53();
        double OMG = 2.0 * PI * mt.genrand_res53();
        double l = 2.0 * PI * mt.genrand_res53(); // mean anomayl
        double u = solve_keplereq(l, ecc); // eccentric anomayl
        OrbParam2PosVel(pos_dummy, pos[i],   vel_dummy, vel[i],
                        mass_sun,  m_planet, ax,        ecc,
                        inc,       OMG,      omg,       u);
        e_ave += ecc*ecc;
    }
    PS::F64 e_ave_glb = sqrt(PS::Comm::getSum(e_ave)/n_glb);
    if(PS::Comm::getRank() == 0){
	std::cerr<<"e_ave_glb="<<e_ave_glb<<std::endl;
	std::cerr<<"e_ave_hill="<<e_ave_glb/h<<std::endl;
	std::cerr<<"e_rms="<<e_rms<<std::endl;
    }
}

template<class Tpsys>
void SetParticleKeplerDisk(Tpsys & psys,
                           const PS::S64 n_glb,
                           PS::S32 & n_loc,
                           PS::F64 & t_sys,
                           const PS::F64 ax_in, // [AU]
                           const PS::F64 ax_out, // [AU]
                           const PS::F64 ecc_rms, // normalized
                           const PS::F64 inc_rms, // normalized
                           const PS::F64 dens = 10.0, // [g/cm^2]
                           const PS::F64 mass_sun = 1.0, //[m_sun]
			   const double a_ice = 0.0,
			   const double f_ice = 1.0,
			   const double power = -1.5,
			   const PS::S32 seed = 0
    ){ 
    PS::S32 my_rank = PS::Comm::getRank();
    PS::S32 n_proc = PS::Comm::getNumberOfProc();
#if 1 // for debug
    if(my_rank == 0){
        n_loc = n_glb;
    }
    else{
        n_loc = 0;
    }
    psys.setNumberOfParticleLocal(n_loc);
#else
    // original
    n_loc = n_glb / n_proc; 
    if( n_glb % n_proc > my_rank) n_loc++;
    psys.setNumberOfParticleLocal(n_loc);
#endif
    PS::F64 * mass;
    PS::F64vec * pos;
    PS::F64vec * vel;
    t_sys = 0.0;

    PS::F64 mass_planet_glb;
    MakeKeplerDisk(mass_planet_glb, mass, pos, vel, n_glb, n_loc,
                   ax_in, ax_out, ecc_rms, inc_rms, dens, mass_sun, a_ice, f_ice, power, seed);
    PS::S64 i_h = n_glb/n_proc*my_rank;
    if( n_glb % n_proc  > my_rank) i_h += my_rank;
    else i_h += n_glb % n_proc;
    for(PS::S32 i=0; i<n_loc; i++){
        psys[i].mass = mass[i];
        psys[i].pos = pos[i];
        psys[i].vel = vel[i];
        psys[i].id = i_h + i;
    }
    delete [] mass;
    delete [] pos;
    delete [] vel;
}

void MakePlummerModel(const double mass_glb,
                      const long long int n_glb,
                      const long long int n_loc,
                      double *& mass,
                      PS::F64vec *& pos,
                      PS::F64vec *& vel,
                      const double eng = -0.25,
                      const int seed = 0){

    assert(eng < 0.0);
    static const double PI = atan(1.0) * 4.0;
    const double r_cutoff = 22.8 / (-3.0 * PI * mass_glb * mass_glb / (64.0 * -0.25)); // 22.8 is cutoff in Nbody units
    //const double r_cutoff = 22.8 * 0.25;
    mass = new double[n_loc];
    pos = new PS::F64vec[n_loc];
    vel = new PS::F64vec[n_loc];

    PS::MTTS mt;
    mt.init_genrand( PS::Comm::getRank() );
    for(int i=0; i<n_loc; i++){
        mass[i] = mass_glb / n_glb;
        double r_tmp = 9999.9;
        while(r_tmp > r_cutoff){ 
            double m_tmp = mt.genrand_res53();
            r_tmp = 1.0 / sqrt( pow(m_tmp, (-2.0/3.0)) - 1.0);
        }
        double phi = 2.0 * PI * mt.genrand_res53();
        double cth = 2.0 * (mt.genrand_real2() - 0.5);
        double sth = sqrt(1.0 - cth*cth);
        pos[i][0] = r_tmp * sth * cos(phi);
        pos[i][1] = r_tmp * sth * sin(phi);
        pos[i][2] = r_tmp * cth;
        while(1){
            const double v_max = 0.1;
            const double v_try = mt.genrand_res53();
            const double v_crit = v_max * mt.genrand_res53();
            if(v_crit < v_try * v_try * pow( (1.0 - v_try * v_try), 3.5) ){
                const double ve = sqrt(2.0) * pow( (r_tmp*r_tmp + 1.0), -0.25 );
                phi = 2.0 * PI * mt.genrand_res53();
                cth = 2.0 * (mt.genrand_res53() - 0.5);
                sth = sqrt(1.0 - cth*cth);
                vel[i][0] = ve * v_try * sth * cos(phi);
                vel[i][1] = ve * v_try * sth * sin(phi);
                vel[i][2] = ve * v_try * cth;
                break;
            }
        }
    }

    PS::F64vec cm_pos = 0.0;
    PS::F64vec cm_vel = 0.0;
    double  cm_mass = 0.0;
    for(PS::S32 i=0; i<n_loc; i++){
        cm_pos += mass[i] * pos[i];
        cm_vel += mass[i] * vel[i];
        cm_mass += mass[i];
    }
    cm_pos /= cm_mass;
    cm_vel /= cm_mass;
    for(PS::S32 i=0; i<n_loc; i++){
        pos[i] -= cm_pos;
        vel[i] -= cm_vel;
    }

    const double r_scale = -3.0 * PI * mass_glb * mass_glb / (64.0 * eng);
    const double coef = 1.0 / sqrt(r_scale);
    for(PS::S32 i=0; i<n_loc; i++){
        pos[i] *= r_scale;
        vel[i] *= coef;
    }

    double r_max_sq = -1.0;
    for(int i=0; i<n_loc; i++){
        if(r_max_sq < pos[i] * pos[i]){
            r_max_sq = pos[i] * pos[i];
        }
    }
//    std::cout<<"r_max= "<<sqrt(r_max_sq)<<std::endl;
}


template<class Tpsys>
void SetParticlePlummer(Tpsys & psys,
                        const PS::S64 n_glb,
                        PS::S32 & n_loc){

    PS::S32 my_rank = PS::Comm::getRank();
    PS::S32 n_proc = PS::Comm::getNumberOfProc();
#if 0 // for debug
    if(my_rank == 0){
        n_loc = n_glb;
    }
    else{
        n_loc = 0;
    }
    psys.setNumberOfParticleLocal(n_loc);
#else
// original
    n_loc = n_glb / n_proc; 
    if( n_glb % n_proc > my_rank) n_loc++;
    psys.setNumberOfParticleLocal(n_loc);
#endif
    PS::F64 * mass;
    PS::F64vec * pos;
    PS::F64vec * vel;

    const PS::F64 m_tot = 1.0;
    const PS::F64 eng = -0.25;
    MakePlummerModel(m_tot, n_glb, n_loc, mass, pos, vel, eng);
    PS::S64 i_h = n_glb/n_proc*my_rank;
    if( n_glb % n_proc  > my_rank) i_h += my_rank;
    else i_h += n_glb % n_proc;
#ifdef BINARY
    for(PS::S32 i=0; i<n_loc; i+=2){
        psys[i].mass = mass[i];
        psys[i].pos = pos[i];
        psys[i].vel = vel[i];
        psys[i].id = i_h + i+1;
        psys[i].status.d = 0;
        psys[i+1].mass = mass[i];
        psys[i+1].pos = pos[i] + PS::F64vec(0.02, 0.02, 0.02);
        psys[i+1].vel = vel[i];
        psys[i+1].id = i_h + i+2;
        psys[i+1].status.d = 0;
        
    }
#else
    for(PS::S32 i=0; i<n_loc; i++){
        psys[i].mass = mass[i];
        psys[i].pos = pos[i];
        psys[i].vel = vel[i];
        psys[i].id = i_h + i + 1;
        psys[i].status.d = 0;
    }
#endif
    delete [] mass;
    delete [] pos;
    delete [] vel;

}

#ifdef MUTIL_ROUT

template<class Tpsys>
void SetBinaryRout(Tpsys & psys, const PS::S32 _n_bin, const PS::F64 g_min, const PS::F64 _r_in, const PS::F64 _r_out, const PS::F64 _m_average) {
    double gamma = std::pow(1.0/g_min,0.33333);
    for (PS::S32 i=0; i<2*_n_bin; i+=2) {
        if (psys[i].r_out<=0||psys[i+i].r_out<=0) {
            double a,ecc;
            PosVel2AxEcc(a,ecc,psys[i].pos, psys[i+1].pos, psys[i].vel, psys[i+1].vel, psys[i].mass, psys[i+1].mass);
            double apo = a*(1+ecc);
            if(apo>0&&apo<_r_in) psys[i+1].r_out = psys[i].r_out = std::max(apo*gamma*std::pow((psys[i].mass+psys[i+1].mass)/_m_average,0.3333),r_out);
            else psys[i+1].r_out = psys[i].r_out = _r_out;
        }
    }
}


template<class Tpsys>
void SetSingleRout(Tpsys & psys, const PS::S32 n, const PS::S32 n_off, const PS::F64 _r_out) {
    for (PS::S32 i=n_off; i<n; i++) if(psys[i].r_out<=0) psys[i].r_out = _r_out;
}

//! update r_out function
template<class Tptcl>
void updateRout(Tptcl** p, const PS::S32 n, const PS::F64 _r_in, const PS::F64 _r_out, const PS::F64 gamma, const PS::F64 _m_average) {
    PS::S32 istart=0,icount=0,ioff[n];
    PS::F64 apomax=0;
    Tptcl* plist[n];
    for(PS::S32 i=0; i<n; i++) {
        PS::F64 ax=0.0,ecc,apo=0.0;
        if(i<n-1) {
            PosVel2AxEcc(ax,ecc,
                         p[i]->pos, p[i+1]->pos,
                         p[i]->vel, p[i+1]->vel,
                         p[i]->mass,p[i+1]->mass);
            //apo=ax*(1.0+ecc);
            apo=ax*(1.0+ecc)*std::pow((p[i]->mass+p[i+1]->mass)/_m_average,0.3333);
        }
        if (ax<0.0||apo>_r_in||i==n-1) {
            if (i==istart) {
                p[i]->r_out=_r_out;
                plist[icount]=p[i];
            }
            else {
                Tptcl** ptemp=new Tptcl*[i-istart+1];
                apomax *=gamma;
                for(PS::S32 j=istart; j<=i; j++) {
                    if(apomax>0.0&&(apomax>1.2*p[j]->r_out||apomax<0.8*p[j]->r_out)) p[j]->r_out = std::max(apomax,_r_out);
                    ptemp[j-istart]=p[j];
                }
                apomax = 0.0;
                plist[icount]=new Tptcl;
                calc_center_of_mass(*(plist[icount]),ptemp,i-istart+1);
                delete[] ptemp;
            }
            istart = i+1;
            ioff[icount] = istart;
            icount++;
        }
        else if (apo>apomax) apomax=apo;
    }
    if (n>icount) {
        updateRout(plist, icount, _r_in, _r_out, gamma, _m_average);
        if(plist[0]->r_out>1.2*p[0]->r_out)
            for (PS::S32 j=0; j<ioff[0]; j++) p[j]->r_out = plist[0]->r_out;
        for (PS::S32 i=1; i<icount; i++) {
            if (plist[i]->r_out>1.2*p[ioff[i-1]]->r_out)
                for (PS::S32 j=ioff[i-1];j<ioff[i];j++) {
                    p[j]->r_out = plist[i]->r_out;
#ifdef HARD_DEBUG
                    assert(p[j]->r_out<_r_out);
#endif
                }
            if (ioff[i]-ioff[i-1]>1) delete plist[i];
        }
    }
}

#endif

