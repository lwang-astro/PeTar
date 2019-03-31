#pragma once

#include <map>

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

