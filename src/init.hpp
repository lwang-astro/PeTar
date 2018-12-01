#pragma once

//!initializaton of system parameters
/*! Obtain Radius parameters, consistent with the input help information
  @param[in]     _tsys:   particle system (soft)
  @param[in,out] _r_in:   changeover function inner boundary
  @param[in,out] _r_out:  changeover function outer boundary
  @param[in,out] _r_bin:  arc group radius criterion
  @param[in,out] _r_search_min: minimum searching radius
  @param[in,out] _r_search_max: maximum searching radius
  @param[in,out] _v_max:  maximum velocity to calculate r_search
  @param[in,out] _m_average: averaged mass of particles
  @param[in,out] _dt_tree: tree time step
  @param[in,out] _vel_disp: system velocity dispersion 
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
                PS::F64 &_dt_tree,
                PS::F64 &_vel_disp,
                const PS::F64 _search_factor,
                const PS::F64 _ratio_r_cut,
                const PS::S64 _n_bin) {

    // local particle number
    const PS::S64 n_loc = _tsys.getNumberOfParticleLocal();

    // local c.m velocity
    PS::F64vec vel_cm_loc = 0.0;
    // local c.m. mass
    PS::F64 mass_cm_loc = 0.0;

    for(PS::S64 i=0; i<n_loc; i++){
        PS::F64 mi = _tsys[i].mass;
        PS::F64vec vi = _tsys[i].vel;

#ifdef MAIN_DEBUG
        assert(mi>0);
#endif
        mass_cm_loc += mi;
        vel_cm_loc += mi * vi;
    }

    // global c.m. parameters
    PS::F64    mass_cm_glb = PS::Comm::getSum(mass_cm_loc);
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
        _r_in = average_mass_glb / (_vel_disp*_vel_disp);
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

    // if r_bin is not defined, set to 0.8 * r_in
    if (_r_bin==0.0) _r_bin = 0.8*_r_in;

    // if r_search_min is not defined, calculate by search_factor*velocity_dispersion*tree_time_step + r_out
    if (_r_search_min==0.0) _r_search_min = _search_factor*_vel_disp*_dt_tree + _r_out;
    // if r_search_max is not defined, calcualte by 5*r_out
    if (_r_search_max==0.0) _r_search_max = 5*_r_out;
    // calculate v_max based on r_search_max, tree time step and search_factor
    _v_max = (_r_search_max - _r_out) / _dt_tree / _search_factor;
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

