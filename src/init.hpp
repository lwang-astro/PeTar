#pragma once

//initializaton-------------------------------------------------------
// Obtain Radius parameters
template<class Tpsys>
void GetR(const Tpsys & system_soft,
          PS::F64 &r_in,
          PS::F64 &r_out,
          PS::F64 &r_bin,
          PS::F64 &r_search_min,
          PS::F64 &m_average,
          PS::F64 &dt,
          PS::F64 &vel_disp,
          const PS::F64 search_factor,
          const PS::F64 ratio_r_cut,
          const PS::S64 n_bin,
          const bool restart_flag){
    //const PS::S32 my_rank=PS::Comm::getRank();
    const PS::S64 n_loc = system_soft.getNumberOfParticleLocal();
    PS::F64vec vel_cm_loc = 0.0;
    PS::F64 mass_cm_loc = 0.0;
    //    PS::F64 mmax_loc = system_soft[0].mass;
    //    PS::F64 mmin_loc = mmax_loc;

    for(PS::S64 i=0; i<n_loc; i++){
        PS::F64 mi = system_soft[i].mass;
        PS::F64vec vi = system_soft[i].vel;
        if(restart_flag) {
            if (system_soft[i].id>0&&system_soft[i].status!=0) continue;
            if (system_soft[i].id<0&&system_soft[i].status<=0) continue;
            if (system_soft[i].id<0) {
                mi = system_soft[i].mass_bk;
                vi = system_soft[i].vel + system_soft[i].acc * dt;
            }
        }
#ifdef MAIN_DEBUG
        assert(mi>0);
#endif
        mass_cm_loc += mi;
        vel_cm_loc += mi * vi;
        //      if (system_soft[i].mass>mmax_loc) mmax_loc = system_soft[i].mass;
        //      if (system_soft[i].mass<mmin_loc) mmin_loc = system_soft[i].mass;
    }

    PS::F64    mass_cm_glb = PS::Comm::getSum(mass_cm_loc);
    PS::F64vec vel_cm_glb  = PS::Comm::getSum(vel_cm_loc);
    //    PS::F64    mmax = PS::Comm::getMaxValue(mmax_loc);
    //    PS::F64    mmin = PS::Comm::getMaxValue(mmin_loc);
      
    vel_cm_glb /= mass_cm_glb;
    PS::F64 vel_sq_loc = 0.0;
    PS::S64 n_vel_loc_count = 0;

    if(restart_flag) {
        for(PS::S64 i=0; i<n_loc; i++) {
            if (system_soft[i].id>0&&system_soft[i].status!=0) continue;
            if (system_soft[i].id<0&&system_soft[i].status<=0) continue;
            PS::F64vec dv;
            if (system_soft[i].id<0) dv = system_soft[i].vel - vel_cm_glb;
            else dv = system_soft[i].vel + system_soft[i].acc*dt - vel_cm_glb;
            vel_sq_loc += dv * dv;
            n_vel_loc_count++;
        }
    }
    else {
        PS::S64 is_start = 0;
        const PS::S64 bin_last_id = 2*n_bin;
        if (system_soft[0].id<bin_last_id) {
            is_start = std::min(bin_last_id - system_soft[0].id + 1,n_loc);
            if(is_start%2!=0) is_start--;
        }
        const PS::S64 ib_start = (system_soft[0].id%2==0)?1:0;
        for (PS::S64 i=ib_start; i<is_start; i+=2) {
            PS::F64 m1 = system_soft[i].mass;
            PS::F64 m2 = system_soft[i+1].mass;
            PS::F64vec dv = (m1*system_soft[i].vel + m2*system_soft[i+1].vel)/(m1+m2) - vel_cm_glb;
            vel_sq_loc += dv * dv;
            n_vel_loc_count++;
        }
        
        for (PS::S64 i=is_start; i<n_loc; i++){
            PS::F64vec dv = system_soft[i].vel - vel_cm_glb;
            vel_sq_loc += dv * dv;
            n_vel_loc_count++;
        }
    }

    const PS::S64    n_vel_glb_count= PS::Comm::getSum(n_vel_loc_count);
    const PS::S64    n_glb          = PS::Comm::getSum(n_loc);
    const PS::F64    vel_sq_glb     = PS::Comm::getSum(vel_sq_loc);
    vel_disp   = sqrt(vel_sq_glb / 3.0 / (PS::F64)n_vel_glb_count);

    PS::F64 average_mass_glb = mass_cm_glb/(PS::F64)n_glb;
    m_average = average_mass_glb;

    bool r_out_flag=r_out>0;
    
    if (r_out_flag) {
        r_in = r_out * ratio_r_cut;
    }
    else {
        r_in = average_mass_glb / (vel_disp*vel_disp);
        // remove 2.33 factor
        // r_out = r_in * pow(mmax/average_mass_glb/gmin,1.0/7.0);
        r_out = r_in / ratio_r_cut;
    }
    
    if (dt==0.0) {
        //PS::F64 dt_origin = 0.125*sqrt(r_out*r_out+r_in*r_in) / vel_disp;
        PS::F64 dt_origin = 0.1*r_out / vel_disp;
        dt = 1.0;
        if (dt_origin<1) while (dt>dt_origin) dt *= 0.5;
        else {
            while (dt<=dt_origin) dt *= 2.0;
            dt *= 0.5;
        }
        if (!r_out_flag) {
            r_out = 10.0*dt*vel_disp;
            r_in = r_out*ratio_r_cut;
        }
    }

    if (r_bin==0.0) r_bin = 0.8*r_in;

    r_search_min = search_factor*vel_disp*dt + r_out;
}

#ifdef MUTIL_ROUT

template<class Tpsys>
void SetBinaryRout(Tpsys & psys, const PS::S32 n_bin, const PS::F64 g_min, const PS::F64 r_in, const PS::F64 r_out, const PS::F64 m_average) {
    double gamma = std::pow(1.0/g_min,0.33333);
    for (PS::S32 i=0; i<2*n_bin; i+=2) {
        if (psys[i].r_out<=0||psys[i+i].r_out<=0) {
            double a,ecc;
            PosVel2AxEcc(a,ecc,psys[i].pos, psys[i+1].pos, psys[i].vel, psys[i+1].vel, psys[i].mass, psys[i+1].mass);
            double apo = a*(1+ecc);
            if(apo>0&&apo<r_in) psys[i+1].r_out = psys[i].r_out = std::max(apo*gamma*std::pow((psys[i].mass+psys[i+1].mass)/m_average,0.3333),r_out);
            else psys[i+1].r_out = psys[i].r_out = r_out;
        }
    }
}


template<class Tpsys>
void SetSingleRout(Tpsys & psys, const PS::S32 n, const PS::S32 n_off, const PS::F64 r_out) {
    for (PS::S32 i=n_off; i<n; i++) if(psys[i].r_out<=0) psys[i].r_out = r_out;
}

//! update r_out function
template<class Tptcl>
void updateRout(Tptcl** p, const PS::S32 n, const PS::F64 r_in, const PS::F64 r_out, const PS::F64 gamma, const PS::F64 m_average) {
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
            apo=ax*(1.0+ecc)*std::pow((p[i]->mass+p[i+1]->mass)/m_average,0.3333);
        }
        if (ax<0.0||apo>r_in||i==n-1) {
            if (i==istart) {
                p[i]->r_out=r_out;
                plist[icount]=p[i];
            }
            else {
                Tptcl** ptemp=new Tptcl*[i-istart+1];
                apomax *=gamma;
                for(PS::S32 j=istart; j<=i; j++) {
                    if(apomax>0.0&&(apomax>1.2*p[j]->r_out||apomax<0.8*p[j]->r_out)) p[j]->r_out = std::max(apomax,r_out);
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
        updateRout(plist, icount, r_in, r_out, gamma, m_average);
        if(plist[0]->r_out>1.2*p[0]->r_out)
            for (PS::S32 j=0; j<ioff[0]; j++) p[j]->r_out = plist[0]->r_out;
        for (PS::S32 i=1; i<icount; i++) {
            if (plist[i]->r_out>1.2*p[ioff[i-1]]->r_out)
                for (PS::S32 j=ioff[i-1];j<ioff[i];j++) {
                    p[j]->r_out = plist[i]->r_out;
#ifdef HARD_DEBUG
                    assert(p[j]->r_out<r_out);
#endif
                }
            if (ioff[i]-ioff[i-1]>1) delete plist[i];
        }
    }
}

#endif

