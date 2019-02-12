#pragma once

#include <map>
#include "kepler.hpp"
#include "hard_force.hpp"
#include "AR.h" /// include AR.h (L.Wang)
#include "ptcl.hpp"

#ifdef HARD_DEBUG
#define DEBUG_ENERGY_LIMIT 1e-6
#define ID_PHASE_SHIFT 4
#endif

#ifdef FIX_STEP_DEBUG
#define STEP_DIVIDER 32.0
#endif

#ifndef ARRAY_ALLOW_LIMIT
#define ARRAY_ALLOW_LIMIT 1000000000
#endif
//const PS::F64 SAFTY_FACTOR_FOR_SEARCH_SQ;

//!leap frog kick for single----------------------------------------------
/* modify the velocity of particle in global system
   reset status to zero
   @param[in,out] _sys: particle system
   @param[in]: _dt: tree step
   @param[in]; _adr: address for single particles
 */
template<class Tsys>
void kickOne(Tsys & _sys, 
             const PS::F64 _dt, 
             const PS::ReallocatableArray<PS::S32>& _adr) {
    const PS::S64 n= _adr.size();
#pragma omp parallel for
    for(PS::S32 i=0; i<n; i++){
        const PS::S32 k=_adr[i];
#ifdef KDKDK_4TH
        _sys[k].vel  += _dt*(_sys[k].acc + 9.0/192.0*_dt*_dt*_sys[k].acorr); 
#else
        _sys[k].vel  += _sys[k].acc * _dt;
#endif
        _sys[k].status = 0;
    }
}

//!leap frog kick for clusters------------------------------------------
/* modify the velocity of particle in local, if particle is from remote note and is not group member, do nothing, need MPI receive to update data
   Recover the mass of members for energy calculation
   @param[in,out] _sys: particle system
   @param[in,out] _ptcl: local particle array in system hard
   @param[in]: _dt: tree step
 */
template<class Tsys, class Tptcl>
void kickCluster(Tsys& _sys,
                 PS::ReallocatableArray<Tptcl>& _ptcl,
                 const PS::F64 _dt) {
    const PS::S64 n= _ptcl.size();
#pragma omp parallel for
    for(PS::S32 i=0; i<n; i++) {
        const PS::S64 cm_adr=-_ptcl[i].status; // notice status is negative 
        const PS::S64 i_adr =_ptcl[i].adr_org;
        // if is group member, recover mass and kick due to c.m. force
        if(cm_adr>0) {
#ifdef HARD_DEBUG
            assert(_ptcl[i].mass_bk>0); 
#endif
            _ptcl[i].mass = _ptcl[i].mass_bk;
#ifdef KDKDK_4TH
            _ptcl[i].vel  += _dt*(_sys[cm_adr].acc + 9.0/192.0*_dt*_dt*_sys[cm_adr].acorr); 
#else
            _ptcl[i].vel += _sys[cm_adr].acc * _dt;
#endif
            // Suppressed because thread unsafe
            //_sys[cm_adr].vel += _sys[cm_adr].acc * _dt/_sys[cm_adr].status; // status has total number of members, to avoid duplicate kick. 
        }
        // non-member particle
        else if(i_adr>=0) {
            // not remote particles
#ifdef KDKDK_4TH
            _ptcl[i].vel  += _dt*(_sys[i_adr].acc + 9.0/192.0*_dt*_dt*_sys[i_adr].acorr); 
#else
            _ptcl[i].vel += _sys[i_adr].acc * _dt;
#endif
        }
    }
}

//!leap frog kick for sending list------------------------------------------
/* Kick single particles in sending list 
   @param[in,out] _sys: particle system
   @param[in,out] _ptcl: local particle array in system hard
   @param[in]: _dt: tree step
 */
template<class Tsys>
void kickSend(Tsys& _sys,
              const PS::ReallocatableArray<PS::S32>& _adr_ptcl_send,
              const PS::F64 _dt) {
    const PS::S64 n= _adr_ptcl_send.size();
#pragma omp parallel for
    for(PS::S32 i=0; i<n; i++) {
        const PS::S64 adr = _adr_ptcl_send[i];
        const PS::S64 cm_adr=-_sys[adr].status; // notice status is negative 
        // if it is group member, should not do kick since c.m. particles are on remote nodes;
        if(cm_adr==0)  {
            _sys[adr].vel += _sys[adr].acc * _dt;
#ifdef KDKDK_4TH
            _sys[adr].vel += _dt*_dt* _sys[adr].acorr /48; 
#endif
        }

#ifdef HARD_DEBUG
        if(cm_adr==0) assert(_sys[adr].mass>0);
        else assert(_sys[adr].mass_bk>0);
#endif
    }
}

//! kick for artifical c.m. particles
/*! Kick c.m. velocity
   @param[in,out] _sys: particle system
   @param[in] _adr_cm_start: c.m. particle starting address
   @param[in] _adr_cm_offset: c.m. address offset
   @param[in] _dt: tree step
 */
template<class Tsys>
void kickCM(Tsys& _sys,
            const PS::S32 _adr_cm_start,
            const PS::S32 _adr_cm_offset,
            const PS::F64 _dt) {
    const PS::S64 n_tot= _sys.getNumberOfParticleLocal();
#pragma omp parallel for
    for(PS::S32 i=_adr_cm_start; i<n_tot; i+= _adr_cm_offset) {
        _sys[i].vel += _sys[i].acc * _dt;
#ifdef KDKDK_4TH
        _sys[i].vel += _dt*_dt* _sys[i].acorr /48; 
#endif
#ifdef HARD_DEBUG
        assert(_sys[i].id<0&&_sys[i].status>0);
#endif
    }
}

template<class Tpsys, class Ttree>
void drift(Tpsys & system,
           const Ttree & tree,
           const PS::F64 dt){
    const PS::S32 n = system.getNumberOfParticleLocal();
#pragma omp parallel for
    for(PS::S32 i=0; i<n; i++){
        //if(tree.getForce(i).n_ngb <= 0){
	if(system[i].n_ngb <= 0){
            system[i].pos  += system[i].vel * dt;
        }
    }
}

//Hermite----------------------------------------------
//! Calculate the maximum time step limit for next block step
/*! Get the maximum time step allown for next block step
  Basic algorithm: the integer of time/dt_min is the binary tree for block step, counting from the minimum digital, the last zero indicate the maximum block step level allown for next step
  @param[in] _time: current time
  @param[in] _dt_max: maximum time step allown
  @param[in] _dt_min: minimum time step allown
 */
PS::F64 calcDtLimit(const PS::F64 _time,
                    const PS::F64 _dt_max,
                    const PS::F64 _dt_min){
    // for first step, the maximum time step is OK
    if(_time==0.0) return _dt_max;
    else {
        // the binary tree for current time position in block step 
        PS::U64 bitmap = _time/_dt_min;
//#ifdef __GNUC__ 
//        PS::S64 dts = __builtin_ctz(bitmap) ;
//        PS::U64 c = (1<<dts);
////        std::cerr<<"time = "<<_time<<"  dt_min = "<<_dt_min<<"  bitmap = "<<bitmap<<"  dts = "<<dts<<std::endl;
//#else

        // block step multiply factor 
        PS::U64 c=1;
        // find the last zero in the binary tree to obtain the current block step level
        while((bitmap&1)==0) {
            bitmap = (bitmap>>1);
            c = (c<<1);
        }
//#endif
        // return the maximum step allown
        return std::min(c*_dt_min,_dt_max);
    }
}

//template <class Tptcl>
//void softKickForCM(Tptcl * ptcl_org,
//                   const PS::S32* cm_list,
//                   const PS::S32  n_cm,
//                   const PS::S32* soft_pert_list,
//                   const PS::F64  dt_soft,
//                   const PS::S32  n_split) {
//    PS::S32 offset = 2*n_split;
//    for (PS::S32 i=0; i<n_cm; i++) {
//        Tptcl* pi = &ptcl_org[cm_list[i]];
//        PS::F64vec fi= PS::F64vec(0.0);
//        const PS::S32* isoft = &soft_pert_list[i*offset];
//        PS::F64 micum = 0.0;
//#ifdef TIDAL_TENSOR
//        for (PS::S32 j=8; j<2*n_split; j++) {
//#else
//        for (PS::S32 j=0; j<2*n_split; j++) {
//#endif
//            Tptcl* pj = &ptcl_org[isoft[j]];
//            fi += pj->mass*pj->vel; // here pj->vel store the soft force of fake members
//            micum += pj->mass;
//#ifdef HARD_DEBUG
//            assert(((pj->status)>>ID_PHASE_SHIFT)==-pi->id);
//#endif
//        }
//#ifdef HARD_DEBUG
//        assert(abs(micum-pi->mass)<1e-10);
//#endif
//        pi->vel += fi/micum * dt_soft;
//    }
//}

//template <class Tptcl>
//void softKickForOneGroup(Tptcl * ptcl_org,
//                         const PS::S32  i_cm,
//                         const PS::S32* group_list,
//                         const PS::S32  group_n,
//                         const PS::S32* soft_pert_list,
//                         const PS::F64  dt_soft,
//                         const PS::S32  n_split) {
//    Tptcl* pi = &ptcl_org[i_cm];
//    PS::F64vec fi= PS::F64vec(0.0);
//    PS::F64 micum = 0.0;
//#ifdef TIDAL_TENSOR
//    for (PS::S32 j=8; j<2*n_split; j++) {
//#else
//    for (PS::S32 j=0; j<2*n_split; j++) {
//#endif
//        Tptcl* pj = &ptcl_org[soft_pert_list[j]];
//        fi += pj->mass*pj->vel; // here pj->vel store the soft force of fake members
//        micum += pj->mass;
//#ifdef HARD_DEBUG
//        assert(((pj->status)>>ID_PHASE_SHIFT)==-pi->id);
//#endif
//    }
//    PS::F64vec vkick = fi/micum * dt_soft;
// 
//#ifdef HARD_DEBUG
//    assert(abs(micum-pi->mass)<1e-10);
//#endif
//    for (PS::S32 i=0; i<group_n; i++) {
//        Tptcl* pk = &ptcl_org[group_list[i]];
//        pk->vel += vkick;
//    }
//    pi->vel += vkick;
//}


//#ifdef ISOLATED
////few-body----------------------------------------------
//template<class Tptcl, class ARC_par>
//void Isolated_Multiple_integrator(Tptcl * ptcl_org,
//                                  const PS::S32 n_ptcl,
//                                  const PS::F64 time_end,
//                                  const PS::F64 dt_limit,
//                                  const PS::F64 rout_single,
//                                  const PS::F64 gamma,
//                                  const PS::F64 m_average,
//#ifdef HARD_CHECK_ENERGY
//                                  PS::F64 &ARC_error_relative,
//                                  PS::F64 &ARC_error,
//                                  PS::S32 N_count[20],
//#endif                         
//                                  const ARC::chainpars &ARC_control,
//                                  ARC_par &Int_pars) {
//    // kepler motion test
//    if (n_ptcl==2) {
//        PS::F64 ax=0,ecc;
//        PS::F64 inc,OMG,omg,tperi;
//        PS::F64 ecc_anomaly_old = PosVel2OrbParam(ax, ecc, inc, OMG, omg, tperi,
//                                                  ptcl_org[0].pos, ptcl_org[1].pos, ptcl_org[0].vel, ptcl_org[1].vel, ptcl_org[0].mass, ptcl_org[1].mass);
//#ifdef HARD_DEBUG
//        std::cerr<<"n_ptcl="<<n_ptcl<<"; ax="<<ax<<"; ecc="<<ecc<<"; peri="<<tperi<<"; pid="<<ptcl_org[0].id<<std::endl;
//#endif
//        if (ax>0.0&&2.0*ax<Int_pars.rin) {
//            // center-of-mass
//            Tptcl pcm;
//            calc_center_of_mass(pcm, ptcl_org, n_ptcl, true);
// 
//            DriveKeplerOrbParam(ptcl_org[0].pos, ptcl_org[1].pos, ptcl_org[0].vel, ptcl_org[1].vel,
//                                ptcl_org[0].mass, ptcl_org[1].mass, time_end, ax, ecc, inc, OMG, omg, ecc_anomaly_old);
//          
//          
//            // integration of center-of-mass
//            pcm.pos += pcm.vel * time_end;
// 
//            center_of_mass_correction(pcm, ptcl_org, n_ptcl);
//          
//            //  PosVel2OrbParam(ax,ecc,inc,OMG,omg,tperi,ptcl_org[0].pos, ptcl_org[1].pos, ptcl_org[0].vel, ptcl_org[1].vel, ptcl_org[0].mass, ptcl_org[1].mass);
//            //  std::cerr<<"A:n_ptcl="<<n_ptcl<<"; ax="<<ax<<"; ecc="<<ecc<<"; peri="<<tperi<<"; pid"<<ptcl_org[0].id<<std::endl;
//            //  if (!kout.is_open()) kout.open("kout");
//            //  for (int i=0;i<n_ptcl;i++) kout<<std::setprecision(17)<<time_origin_<<" "<<ptcl_org[i].mass<<" "<<ptcl_org[i].pos<<" "<<ptcl_org[i].vel<<std::endl;
// 
//            PS::F64 peri = ax*(1+ecc)*gamma*std::pow(pcm.mass/m_average,0.3333);
//            if (peri>1.2*ptcl_org[0].r_out || (peri>0 && peri<0.8*ptcl_org[0].r_out) || ptcl_org[0].r_out!= ptcl_org[1].r_out)
//                ptcl_org[0].r_out = ptcl_org[1].r_out = std::max(peri,rout_single);
//            return;
//        }
//    }
//    else if(n_ptcl==3) {
//        // Not yet implementd
//    }
//      
//    ARC::chain<Tptcl> c((std::size_t)n_ptcl);
//    static thread_local PS::F64 time_sys = 0.0;
// 
//    c.addP(n_ptcl,ptcl_org);
// 
//#ifdef HARD_DEBUG
//    for (PS::S32 i=0; i<n_ptcl; i++) 
//        if (ptcl_org[i].r_out<Int_pars.rin) {
//            std::cerr<<"Error, updated p["<<i<<"].rout ("<<ptcl_org[i].r_out<<") < r_in ("<<Int_pars.rin<<")"<<std::endl;
//            abort();
//        }
//#endif
// 
//    //c.link_int_par(Int_pars);
//    c.init(time_sys,ARC_control,Int_pars);
// 
//#ifdef HARD_CHECK_ENERGY
//    PS::F64 ARC_error_once = c.getPot()+c.getEkin();
//    if(n_ptcl<=20) N_count[n_ptcl-1]++;
//    else std::cerr<<"Large cluster formed, n="<<n_ptcl<<std::endl;
//#endif
//      
//    PS::F64 dscoff=1.0;
//    PS::F64 ds_up_limit = 0.25*dt_limit/c.calc_dt_X(1.0);
//    PS::F64 ds_use = c.calc_next_step_custom();
//      
//    if (ds_use>ds_up_limit) ds_use = ds_up_limit;
// 
//    // convergency check
//    PS::S32 converge_count=0;
//    PS::S32 error_count=0;
//    bool modify_step_flag=false;
//    bool final_flag=false;
//      
//    while(time_end-c.getTime()>ARC_control.dterr) {
// 
//        // if (ptcl_org[0].id==8&&time_origin_==0.0078125) {
//        //	 std::cout<<"ds= "<<ds_use<<" toff= "<<time_end<<std::endl;
//    		
//        //	 FILE* fout=fopen("data","w");
//        //	 fwrite(ptcl_org,sizeof(Tptcl),n_ptcl,fout);
//        //	 fclose(fout);
//        // for (PS::S32 i=0; i<n_ptcl; i++) 
//        //	 std::cout<<ptcl_org[i].getMass()<<" "<<ptcl_org[i].getPos()[0]<<" "<<ptcl_org[i].getPos()[1]<<" "<<ptcl_org[i].getPos()[2]<<" "<<ptcl_org[i].getVel()[0]<<" "<<ptcl_org[i].getVel()[1]<<" "<<ptcl_org[i].getVel()[2]<<std::endl;
//        // }
// 
//        PS::F64 dsf=c.extrapolation_integration(ds_use,ARC_control,time_end,Int_pars);
// 
//        // std::cerr<<"Particle=";
//        // for (PS::S32 i=0; i<n_ptcl; i++) std::cerr<<ptcl_org[i].id<<" ";
//        // std::cerr<<"n="<<n_ptcl<<" Time_end="<<time_end<<" ctime="<<c.getTime()<<" diff="<<time_end-c.getTime()<<" ds="<<ds_use<<" dsf="<<dsf<<std::endl;
// 
//        if (dsf<0) {
//            final_flag=true;
//            converge_count++;
//            if (converge_count>5&&time_end-c.getTime()>ARC_control.dterr*100) {
//                std::cerr<<"Error: Time synchronization fails!\nStep size ds: "<<ds_use<<"\nEnding physical time: "<<time_end<<"\nTime difference: "<<time_end-c.getTime()<<"\nR_in: "<<Int_pars.rin<<"\nR_out: "<<Int_pars.rout<<"\n";
//                //  ds_use = 0.1*c.calc_next_step_custom();
//                //  std::cerr<<"New step size: "<<ds_use<<std::endl;
//                //  modify_step_flag=true;
//                //	converge_count=0;
//                c.dump("ARC_dump.dat");
//                ARC_control.dump("ARC_dump.par");
//                c.print(std::cerr);
//                abort();
//            }
//            else ds_use *= -dsf;
//            // debuging
//            // if (ptcl_org[0].id==267) {
//            //		c.dump("ARC_dump.dat");
//            //		ARC_control.dump("ARC_dump.par");
//            //		c.print(std::cerr);
//            //		abort();
//            // }
//        }
//        else if (dsf==0) {
//            //          char collerr[50]="two particle overlap!";
//            c.info->ErrMessage(std::cerr);
//            error_count++;
//            if(error_count>4) {
//                std::cerr<<"Error: Too much error appear!\nStep size ds: "<<ds_use<<"\nEnding physical time: "<<time_end<<"\nTime difference: "<<time_end-c.getTime()<<"\nR_in: "<<Int_pars.rin<<"\nR_out: "<<Int_pars.rout<<"\n";
//                c.dump("ARC_dump.dat");
//                ARC_control.dump("ARC_dump.par");
//                c.print(std::cerr);
//                abort();
//            }
//            if (c.info->status==5) {
//                dscoff = 0.25;
//                ds_use *= dscoff;
//            }
//            //          else if (c.info->status==6) ds_use *= 0.001;
//            else if (c.info->status==4) ds_use = std::min(dscoff*c.calc_next_step_custom(),ds_up_limit);
//            else ds_use *= 0.1;
//            modify_step_flag=true;
//        }
//        else  {
//            if (final_flag) {
//                if (converge_count>6&&time_end-c.getTime()>ARC_control.dterr*100) {
//                    std::cerr<<"Error: Time synchronization fails!\nStep size ds: "<<ds_use<<"\nEnding physical time: "<<time_end<<"\nTime difference: "<<time_end-c.getTime()<<"\nR_in: "<<Int_pars.rin<<"\nR_out: "<<Int_pars.rout<<"\n";
//                    c.dump("ARC_dump.dat");
//                    ARC_control.dump("ARC_dump.par");
//                    c.print(std::cerr);
//                    abort();
//                }
//                converge_count++;
//            }
//            else if (n_ptcl>2||(modify_step_flag&&error_count==0)) {
//                ds_use = std::min(dscoff*c.calc_next_step_custom(),ds_up_limit);
//                modify_step_flag=false;
//            }
//            // reducing error counter if integration success, this is to avoid the significant change of step may cause some issue
//            if(error_count>0) error_count--;
//        }
//    }
// 
//    // update Rout
//    if(n_ptcl>2) {
//        std::size_t* list=new std::size_t[n_ptcl];
//        Tptcl** plist=new Tptcl*[n_ptcl];
//        c.getList(list);
//#ifdef DEBUG_TEMP
//        std::cerr<<"Curren r_out: n="<<n_ptcl;
//        for (PS::S32 i=0; i<n_ptcl; i++) std::cerr<<std::setw(14)<<list[i]<<" ";
//        std::cerr<<std::endl;
//        for (PS::S32 i=0; i<n_ptcl; i++) std::cerr<<std::setw(14)<<ptcl_org[list[i]].r_out<<" ";
//        std::cerr<<std::endl;
//#endif
//        for (PS::S32 i=0; i<n_ptcl; i++) plist[i] = &(ptcl_org[list[i]]);
//        updateRout(plist,n_ptcl,Int_pars.rin,rout_single,gamma,m_average);
//#ifdef DEBUG_TEMP
//        std::cerr<<"new r_out: n="<<n_ptcl;
//        for (PS::S32 i=0; i<n_ptcl; i++) std::cerr<<std::setw(14)<<list[i]<<" ";
//        std::cerr<<std::endl;
//        for (PS::S32 i=0; i<n_ptcl; i++) std::cerr<<std::setw(14)<<ptcl_org[list[i]].r_out<<" ";
//        std::cerr<<std::endl;
//#endif
//#ifdef HARD_DEBUG
//        for (PS::S32 i=0; i<n_ptcl; i++) 
//            if (ptcl_org[list[i]].r_out<Int_pars.rin) {
//                std::cerr<<"Error, updated p["<<list[i]<<"].rout ("<<ptcl_org[list[i]].r_out<<") < r_in ("<<Int_pars.rin<<")"<<std::endl;
//                abort();
//            }
//#endif
//        delete[] list;
//        delete[] plist;
//    }
// 
//    // error record
//#ifdef HARD_CHECK_ENERGY
//    PS::F64 ARC_error_temp = (c.getPot()+c.getEkin()-ARC_error_once);
//    ARC_error += ARC_error_temp;
//    ARC_error_relative += ARC_error_temp/ARC_error_once;
//#endif
//      
//    // integration of center-of-mass
//    c.cm.pos += c.cm.vel * time_end;
// 
//    c.center_shift_inverse();
// 
//    // if (!arout.is_open()) arout.open("arout");
//    // for (int i=0;i<n_ptcl;i++) arout<<std::setprecision(17)<<time_origin_<<" "<<ptcl_org[i].mass<<" "<<ptcl_org[i].pos<<" "<<ptcl_org[i].vel<<std::endl;
//}
//#endif

