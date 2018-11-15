#pragma once

#ifndef TIDAL_TENSOR
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#endif
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
        _sys[k].vel  += _sys[k].acc * _dt;
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
            _ptcl[i].vel += _sys[cm_adr].acc * _dt;
            // Suppressed because thread unsafe
            //_sys[cm_adr].vel += _sys[cm_adr].acc * _dt/_sys[cm_adr].status; // status has total number of members, to avoid duplicate kick. 
        }
        // non-member particle
        else if(i_adr>=0) {
            // not remote particles
            _ptcl[i].vel += _sys[i_adr].acc * _dt;
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
        if(cm_adr==0)  _sys[adr].vel += _sys[adr].acc * _dt;
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
PS::F64 calcDtLimit(const PS::F64 time_sys,
                    const PS::F64 dt_limit_org,
                    const PS::F64 dtmin){
    if(time_sys==0.0) return dt_limit_org;
    else {
        PS::U64 bitmap = time_sys/dtmin;
//#ifdef HARD_DEBUG
////        assert(bitmap*dtmin==time_sys);
//        if(bitmap==0) {
//            std::cerr<<"Error: dt_min_hard not small enough, time = "<<time_sys<<"; dt_min_hard = "<<dtmin<<std::endl;
//            abort();
//        }
//#endif
//#ifdef __GNUC__ 
//        PS::S64 dts = __builtin_ctz(bitmap) ;
//        PS::U64 c = (1<<dts);
////        std::cerr<<"time = "<<time_sys<<"  dtmin = "<<dtmin<<"  bitmap = "<<bitmap<<"  dts = "<<dts<<std::endl;
//#else
        PS::U64 c=1;
        while((bitmap&1)==0) {
            bitmap = (bitmap>>1);
            c = (c<<1);
        }
//#endif
        return std::min(c*dtmin,dt_limit_org);
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

class PtclH4: public Ptcl{
public:
    PS::F64vec acc0;
    PS::F64vec acc1;
    PS::F64 dt;
    PS::F64 time;
#ifdef HARD_DEBUG_ACC
    PS::F64vec acc2; // for debug
    PS::F64vec acc3; // for debug
#endif

    PtclH4(): dt(0.0), time(0.0) {}
    
    template<class Tptcl>
    PtclH4(const Tptcl &_p): Ptcl(_p), dt(0.0), time(0.0) {}

    void dump(FILE *fp) {
        fwrite(this, sizeof(PtclH4),1,fp);
    }

    void read(FILE *fp) {
        size_t rcount = fread(this, sizeof(PtclH4),1,fp);
        if (rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            abort();
        }
    }
};

class PtclPred{
public:
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64 mass;
    PS::F64 r_search;
};

class PtclForce{
public:
    PS::F64vec acc0; //
    PS::F64vec acc1; //
    void clear(){
        acc0 = acc1 = 0.0;
    }

    void dump(FILE *fp) {
        fwrite(this, sizeof(PtclForce),1,fp);
    }

    void read(FILE *fp) {
        size_t rcount = fread(this, sizeof(PtclForce),1,fp);
        if (rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            abort();
        }
    }
};

//! Neighbor information collector
class NeighborInfo{
public:
    PS::S32 r_min_index; // nearest neighbor index for each ptcl
    PS::F64 r_min2;      // nearest neighbor distance square
    PS::F64 min_mass;  // mimimum mass in neighbors

    NeighborInfo(): r_min_index(-1), r_min2(PS::LARGE_FLOAT), min_mass(PS::LARGE_FLOAT) {}
};

template <class Tphard>
class HermiteIntegrator{
private:
    PS::ReallocatableArray<PtclPred> pred_;
    PS::ReallocatableArray<PtclForce> force_;
    PS::ReallocatableArray<PS::S32> adr_dt_sorted_;
    PS::ReallocatableArray<PS::F64> time_next_;

    PS::ReallocatableArray<PtclH4> ptcl_;   // first c.m.; second single
    PS::ReallocatableArray<Tphard*> ptcl_ptr_; // original ptcl address
    PS::ReallocatableArray<PS::S32> Jlist_; // neighbor list 
    PS::ReallocatableArray<PS::S32> Jlist_disp_; // neighbor list offset 
    PS::ReallocatableArray<PS::S32> Jlist_n_;    // number of neighbors
    PS::ReallocatableArray<NeighborInfo> nb_info_; // neighbor information 
    PS::S32 n_nb_off_;

    ParticleBase pcm_; // c.m. of the cluster

    PS::F64 a0_offset_sq_; /// time step parameter
    PS::S32 n_act_;        /// active particle number

    //  import parameter
    // PS::F64 dt_limit_hard_; // time step limit     
    PS::F64 eta_s_;         // time step parameter 
    PS::F64 r_in_;          // force parameter     
    PS::F64 r_out_;         // force parameter     
    PS::F64 r_oi_inv_;      // 1.0/(r_out-r_in)
    PS::F64 r_A_;           // (r_out-r_in)/(r_out+r_in)
    PS::F64 eps_sq_;        // softening parameter 

    PS::F64 CalcDt2nd(const PS::F64vec & acc0, 
                      const PS::F64vec & acc1, 
                      const PS::F64 eta, 
                      const PS::F64 a0_offset_sq=0.0){
        const PS::F64 s0 = acc0 * acc0 + a0_offset_sq;
        const PS::F64 s1 = acc1 * acc1;
        if(s0 == a0_offset_sq || s1 == 0.0){
            return PS::LARGE_FLOAT;
        }
        else{
            return eta * sqrt( s0 / s1 );
        }
    }

    PS::F64 CalcDt4th(const PS::F64vec & acc0,
                      const PS::F64vec & acc1,
                      const PS::F64vec & acc2,
                      const PS::F64vec & acc3,
                      const PS::F64 eta, 
                      const PS::F64 a0_offset_sq=0.0){
        const PS::F64 s0 = acc0 * acc0 + a0_offset_sq;
        const PS::F64 s1 = acc1 * acc1;
        const PS::F64 s2 = acc2 * acc2;
        const PS::F64 s3 = acc3 * acc3;
        if(s0 == a0_offset_sq || s1 == 0.0){
            return PS::LARGE_FLOAT;
        }
        else{
            return eta * sqrt( (sqrt(s0*s2) + s1) / (sqrt(s1*s3) + s2) );
        }
    }


    //! Center-of-mass frame shift for particles
    /*! Shift positions and velocities of particles in #p from their original frame to their center-of-mass frame\n
      Notice the center-of-mass position and velocity use values from #cm
    */
    void calcCMAndShift(ParticleBase &pcm, PtclH4* ptcl, const PS::S32 n) {
        pcm.mass = 0.0;
        pcm.pos = pcm.vel = 0.0;
        for (PS::S32 i=0; i<n; i++) {
            PS::F64 mi=ptcl[i].mass;
            pcm.mass += mi;
            pcm.pos += mi*ptcl[i].pos;
            pcm.vel += mi*ptcl[i].vel;
        }
        PS::F64 invm = 1.0/pcm.mass;
        pcm.pos *= invm;
        pcm.vel *= invm;
        
        for (PS::S32 i=0;i<n;i++) {
//            std::cerr<<std::setprecision(14)<<"cm sbore i "<<i<<" pos "<<ptcl[i].pos<<" vel "<<ptcl[i].vel<<std::endl;
            ptcl[i].pos -= pcm.pos;
            ptcl[i].vel -= pcm.vel;
//#ifdef HARD_DEBUG
//            std::cerr<<"cm shift i "<<i<<" pos "<<ptcl[i].pos<<" vel "<<ptcl[i].vel<<std::endl;
//#endif
        }
    }

    //! shift ptcl from c.m. frame to original frame
    /*!
     */
    void shiftBackCM(const ParticleBase &pcm, PtclH4* ptcl, const PS::S32 n) {
#ifdef HARD_DEBUG
        PS::F64 m=0;
        for (PS::S32 i=0;i<n;i++) m+= ptcl[i].mass;
        assert(abs(m-pcm.mass)<1e-10);
        //PS::F64vec pos=0,vel=0;
        //for (PS::S32 i=0;i<n;i++) {
        //    pos += ptcl[i].mass*ptcl[i].pos;
        //    vel += ptcl[i].mass*ptcl[i].vel;
        //}
        //pos /=m;
        //vel /=m;
        //std::cerr<<std::setprecision(14)<<"cm pos "<<pos<<" vel "<<vel<<std::endl;
#endif
        for (PS::S32 i=0;i<n;i++) {
            //std::cerr<<"cm bbck i "<<i<<" pos "<<ptcl[i].pos<<" vel "<<ptcl[i].vel<<std::endl;
            ptcl[i].pos += pcm.pos;
            ptcl[i].vel += pcm.vel;
            //std::cerr<<"cm back i "<<i<<" pos "<<ptcl[i].pos<<" vel "<<ptcl[i].vel<<std::endl;
        }
    }

    //
    /* \return fail_flag: if the time step < dt_min, return true (failure)
     */
    template <class Tptcl>
    bool CalcBlockDt2ndAct(Tptcl ptcl[],
                           const PtclForce force[],
                           const PS::S32 adr_array[],
                           const PS::S32 n_act, 
                           const PS::F64 eta,
                           const PS::F64 dt_max,
                           const PS::F64 dt_min,
                           const PS::F64 a0_offset_sq){
        bool fail_flag=false;
        for(PS::S32 i=0; i<n_act; i++){
            const PS::S32 adr = adr_array[i];
            const PS::F64vec a0 = force[adr].acc0;
            const PS::F64vec a1 = force[adr].acc1;

#ifdef FIX_STEP_DEBUG
            ptcl[adr].dt = dt_max/STEP_DIVIDER;
#else        
            const PS::F64 dt_ref = CalcDt2nd(a0, a1, eta, a0_offset_sq);
            PS::F64 dt = dt_max;
            while(dt > dt_ref) dt *= 0.5;
            ptcl[adr].dt = dt;
#endif
            if(dt<dt_min) {
                std::cerr<<"Error: Hermite integrator initial step size ("<<dt<<") < dt_min ("<<dt_min<<")!"<<std::endl;
                std::cerr<<"i="<<i<<" adr="<<adr<<" acc="<<a0<<" acc1="<<a1<<" eta="<<eta<<" a0_offset_sq="<<a0_offset_sq<<std::endl;
                fail_flag=true;
            }
        }
        return fail_flag;
    }

//    //! Calculate acc and jerk for active particles from neighbor particles
//    /*!
//       @param[out] _force: acc and jerk array for output
//       @param[out] _nb_info: neighbor information
//       @param[in] _ptcl: particle list
//       @param[in] _Ilist: active i particle index in ptcl_
//       @param[in] _n_act: active particle number
//       @param[in] _Jlist: New neighbor list for each i particle
//       @param[in] _Jlist_disp: Neighbor list offset for each i particle
//       @param[in] _Jlist_n: New number of neighbors for each i particle
//       @param[in] _rin: inner radius of soft-hard changeover function
//       @param[in] _rout: outer radius of soft-hard changeover function
//       @param[in] _r_oi_inv: 1.0/(_rout-_rin);
//       @param[in] _r_A: (_rout-_rin)/(_rout+_rin);
//       @param[in] _eps2: softing eps square
//       @param[in] _Aint: ARC integrator class
//       @param[in] _n_group: number of groups
//     */
//    template <class Tptcl, class ARCint>
//    inline void CalcAcc0Acc1ActNb(PtclForce _force[],
//                                  const Tptcl _ptcl[],
//                                  const PS::S32 _Ilist[],
//                                  const PS::S32 _n_act,
//                                  const PS::S32 _Jlist[],
//                                  const PS::S32 _Jlist_disp[],
//                                  const PS::S32 _Jlist_n[],
//                                  const PS::F64 _rin,
//                                  const PS::F64 _rout,
//                                  const PS::F64 _r_oi_inv,
//                                  const PS::F64 _r_A,
//                                  const PS::F64 _eps2,
//                                  const ARCint *_Aint,
//                                  const PS::S32 _n_group) {
//        // PS::ReallocatableArray< std::pair<PS::S32, PS::S32> > & merge_pair ){
//        // active particles
//        for(PS::S32 i=0; i<_n_act; i++){
//            const PS::S32 iadr = _Ilist[i];
//            _force[iadr].acc0 = _force[iadr].acc1 = 0.0;
//            // all
//            const PS::S32 joff = _Jlist_disp[iadr];
//            const PS::S32 jn = _Jlist_n[iadr];
//            for(PS::S32 j=0; j<jn; j++){
//                PS::S32 jadr = _Jlist[joff+j];
//#ifdef HARD_DEBUG
//                assert(iadr!=jadr);
//                PS::F64 mcmcheck =0.0;
//#endif
//                if (jadr<_n_group) {
//                    if(_Aint->getMask(jadr)) continue; 
//                    const auto* pj = _Aint->getGroupPtcl(jadr);
//                    PS::F64 sd = _Aint->getSlowDown(jadr);
//                    for(PS::S32 k=0; k<_Aint->getGroupN(jadr); k++) {
//                        PS::F64 r2 = 0.0;
//                        CalcAcc0Acc1R2Cutoff(_ptcl[iadr].pos, _ptcl[iadr].vel,
//                                             _force[iadr].acc0, _force[iadr].acc1, r2,
//                                             pj[k].pos, pj[k].vel/sd, pj[k].mass,
//                                             _eps2, _rout, _rin, _r_oi_inv, _r_A);
//#ifdef HARD_DEBUG
//                        mcmcheck += pj[k].mass;
//                        //std::cerr<<k<<" P "<<pj[k].pos<<" v "<<pj[k].vel<<" sd "<<sd<<std::endl;
//#endif
//                        if(r2<_r_min2) {
//                            _r_min2 = r2;
//                            _j_r_min = j;
//                        }
//                    }
//#ifdef HARD_DEBUG
//                    assert(abs(mcmcheck-_ptcl[jadr].mass)<1e-10);
//                    assert(mcmcheck>0.0);
//#endif                    
//                }
//                else {
//                    PS::F64 r2 = 0.0;
//                    //PS::F64 rout = std::max(ptcl[iadr].r_out, ptcl[jadr].r_out);
//                    CalcAcc0Acc1R2Cutoff(_ptcl[iadr].pos, _ptcl[iadr].vel,
//                                         _force[iadr].acc0, _force[iadr].acc1, r2,
//                                         _ptcl[jadr].pos, _ptcl[jadr].vel, _ptcl[jadr].mass,
//                                         _eps2, _rout, _rin, _r_oi_inv, _r_A);
//                    if(r2<_r_min2) {
//                        _r_min2 = r2;
//                        _j_r_min = j;
//                    }
//                }
//                // if(r2 < ((ptcl[adr].r_merge + ptcl[j].r_merge)*(ptcl[adr].r_merge + ptcl[j].r_merge)) && ptcl[j].mass > 0.0){
//                //     merge_pair.push_back( std::make_pair(adr, j) );
//                // }
//            }
//        }
//    }

    //! calculate one particle f and fdot from all j particles 
    /*! Calculate f and fdot from all j particles, and return neighbor map
       @param[out] _force: acc and jerk 
       @param[out] _nb_flag: neighbor flag map with size of j number, if ture, index j is the neighbor
       @param[out] _nb_info: neighbor information collector
       @param[in] _pi: i particle
       @param[in] _iadr: i particle index in _ptcl
       @param[in] _vcmsdi: c.m. of i particle if i is the member of a group, single case it is zero vector
       @param[in] _sdi: slowdown factor for the group which i belong to, single case it is one
       @param[in] _ptcl: j particle array
       @param[in] _n_tot: total j particle number
       @param[in] _n_group: number of groups
       @param[in] _rin:  inner radius of soft-hard changeover function
       @param[in] _rout: outer radius of soft-hard changeover function
       @param[in] _r_oi_inv: 1.0/(_rout-_rin);
       @param[in] _r_A: (_rout-_rin)/(_rout+_rin);
       @param[in] _eps2: softing eps square
       @param[in] _Aint: ARC integrator class
     */
    template <class Tpi, class Tptcl, class ARCint>
    inline void CalcOneAcc0Acc1(PtclForce &_force, 
                                bool _nb_flag[],
                                NeighborInfo& _nb_info,
                                const Tpi &_pi,
                                const PS::S32 _iadr,
                                const PS::F64vec &_vcmsdi,
                                const PS::F64 _sdi,
                                const Tptcl _ptcl[],
                                const PS::S32 _n_tot,
                                const PS::S32 _n_group,
                                const PS::F64 _rin,
                                const PS::F64 _rout,
                                const PS::F64 _r_oi_inv,
                                const PS::F64 _r_A,
                                const PS::F64 _eps2,
                                const ARCint* _Aint = NULL) {
        for(PS::S32 j=0; j<_n_group; j++) {
            if(_iadr==j) continue;
            if(_Aint->getMask(j)) continue;
#ifdef HARD_DEBUG
            PS::F64 mcmcheck =0.0;
#endif
            const auto* pj = _Aint->getGroupPtcl(j);
            PS::F64 sdj = 1.0/_Aint->getSlowDown(j);
            PS::F64vec vcmsdj = (1.0-sdj)*_ptcl[j].vel;
            for(PS::S32 k=0; k<_Aint->getGroupN(j); k++) {
                PS::F64 r2 = 0.0;
                CalcAcc0Acc1R2Cutoff(_pi.pos, _pi.vel*_sdi+_vcmsdi,
                                     _force.acc0, _force.acc1, r2,
                                     pj[k].pos, pj[k].vel*sdj+vcmsdj, pj[k].mass,
                                     _eps2, _rout, _rin, _r_oi_inv, _r_A);
#ifdef HARD_DEBUG
                mcmcheck += pj[k].mass;
#endif
                PS::F64 rs=std::max(_pi.r_search,pj[k].r_search);
                if(r2<=rs*rs) _nb_flag[j] = true;
                if(r2<_nb_info.r_min2) {
                    _nb_info.r_min2 = r2;
                    _nb_info.r_min_index = j;
                }
                _nb_info.min_mass = std::min(_nb_info.min_mass,_ptcl[j].mass);
            }
#ifdef HARD_DEBUG
            assert(abs(mcmcheck-_ptcl[j].mass)<1e-10);
            assert(mcmcheck>0.0);
#endif                    
        }
        for(PS::S32 j=_n_group; j<_n_tot; j++){
            if(_iadr==j) continue;
            //PS::F64 rout = std::max(ptcl[iadr].r_out, ptcl[j].r_out);
            PS::F64 r2 = 0.0;
            CalcAcc0Acc1R2Cutoff(_pi.pos, _pi.vel*_sdi+_vcmsdi,
                                 _force.acc0, _force.acc1, r2,
                                 _ptcl[j].pos, _ptcl[j].vel, _ptcl[j].mass,
                                 _eps2, _rout, _rin, _r_oi_inv, _r_A);
            // if(r2 < ((ptcl[adr].r_merge + ptcl[j].r_merge)*(ptcl[adr].r_merge + ptcl[j].r_merge)) && ptcl[j].mass > 0.0){
            //     merge_pair.push_back( std::make_pair(adr, j) );
            // }
            PS::F64 rs=std::max(_pi.r_search, _ptcl[j].r_search);
            if(r2<=rs*rs) _nb_flag[j] = true;
            if(r2<_nb_info.r_min2) {
                _nb_info.r_min2 = r2;
                _nb_info.r_min_index = j;
            }
            _nb_info.min_mass = std::min(_nb_info.min_mass,_ptcl[j].mass);
        }
    }

    //! Calculate acc and jerk for active particles from all particles and update neighbor lists
    /*! calculate force and update neighbor lists
       @param[out] _force: acc and jerk array
       @param[out] _nb_info: neighbor information collector
       @param[out] _Jlist: New neighbor list for each i particle
       @param[out] _Jlist_n: New number of neighbors for each i particle
       @param[in] _Jlist_disp: Neighbor list offset for each i particle
       @param[in] _ptcl: particle list
       @param[in] _n_tot: total number of particles
       @param[in] _Ilist: active i particle index in ptcl_
       @param[in] _n_act: active particle number
       @param[in] _rin: inner radius of soft-hard changeover function
       @param[in] _rout: outer radius of soft-hard changeover function
       @param[in] _r_oi_inv: 1.0/(_rout-_rin);
       @param[in] _r_A: (_rout-_rin)/(_rout+_rin);
       @param[in] _eps2: softing eps square
       @param[in] _Aint: ARC integrator class
     */
    template <class Tptcl, class ARCint>
    inline void CalcActAcc0Acc1(PtclForce _force[],
                                NeighborInfo _nb_info[],
                                PS::S32 _Jlist[],
                                PS::S32 _Jlist_n[],
                                const PS::S32 _Jlist_disp[],
                                const Tptcl _ptcl[],
                                const PS::S32 _n_tot,
                                const PS::S32 _Ilist[],
                                const PS::S32 _n_act,
                                const PS::F64 _rin,
                                const PS::F64 _rout,
                                const PS::F64 _r_oi_inv,
                                const PS::F64 _r_A,
                                const PS::F64 _eps2,
                                const ARCint* _Aint) {
        PS::S32 n_group = 0;
        if(_Aint!=NULL) n_group = _Aint->getNGroups();
        const PS::F64vec vzero = PS::F64vec(0.0);
        // PS::ReallocatableArray< std::pair<PS::S32, PS::S32> > & merge_pair ){

        // active iparticles loop
        for(PS::S32 i=0; i<_n_act; i++){
            const PS::S32 iadr = _Ilist[i];
            bool nb_flag[_n_tot];
            for (PS::S32 j=0; j<_n_tot; j++) nb_flag[j]=false;   // Neighbor flag for i particle
            _force[iadr].acc0 = _force[iadr].acc1 = 0.0;
            
            _nb_info[iadr].r_min_index = -1;
            _nb_info[iadr].r_min2 = PS::LARGE_FLOAT;
            _nb_info[iadr].min_mass = PS::LARGE_FLOAT;
            
            // for group particle
            if (iadr<n_group) {
#ifdef HARD_DEBUG
                PS::F64 mcmcheck =0.0;
#endif
                const PS::S32 ni = _Aint->getGroupN(iadr);             // number of members
#ifdef HARD_DEBUG
                assert(ni>0);
#endif
                const auto* pi = _Aint->getGroupPtcl(iadr);
                const PS::F64 sdi = 1.0/_Aint->getSlowDown(iadr);      // slowdown factor
                const PS::F64vec vcmsdi = (1.0-sdi)*_ptcl[iadr].vel;   // slowdown velocity
                PtclForce fp[ni];
//                PS::F64 r2min[_n_tot]={PS::LARGE_FLOAT};

                for (PS::S32 j=0; j<ni; j++) {
                    fp[j].acc0 = fp[j].acc1 = 0.0;
                    CalcOneAcc0Acc1(fp[j], nb_flag, _nb_info[iadr], pi[j], iadr, vcmsdi, sdi, _ptcl, _n_tot, n_group, _rin, _rout, _r_oi_inv, _r_A, _eps2, _Aint);
                    // c.m. force
                    _force[iadr].acc0 += pi[j].mass*fp[j].acc0;
                    _force[iadr].acc1 += pi[j].mass*fp[j].acc1;
                    
#ifdef HARD_DEBUG
                    mcmcheck += pi[j].mass;
#endif
                }

#ifdef HARD_DEBUG
                assert(abs(mcmcheck-_ptcl[iadr].mass)<1e-10);
                assert(mcmcheck>0.0);
                assert(_ptcl[iadr].mass>0);
#endif                    
                // c.m. force
                _force[iadr].acc0 /= _ptcl[iadr].mass;
                _force[iadr].acc1 /= _ptcl[iadr].mass;

//                //update perturber list
//                PS::F64 rsearchi = _ptcl[iadr].r_search;
//                for (PS::S32 j=0; j<_n_tot; j++) {
//                    PS::F64 rsearch= std::max(_ptcl[j].r_search, rsearchi);
//                    if (r2min[j] < rsearch*rsearch) {
//                        _Aint.addPert(iadr,_ptcl[j],_force[j]);
//                    }
//                }
            }
            else {
//                PS::F64 r2[_n_tot];
                CalcOneAcc0Acc1(_force[iadr], nb_flag, _nb_info[iadr], _ptcl[iadr], iadr, vzero, 1.0, _ptcl, _n_tot, n_group, _rin, _rout, _r_oi_inv, _r_A, _eps2, _Aint);
            }

            // Update neighbors
            PS::S32 nb_disp = _Jlist_disp[iadr];
            PS::S32 nb=0;
            for (PS::S32 j=0; j<_n_tot; j++) {
                if(nb_flag[j]) _Jlist[nb_disp+nb++] = j;
            }
            _Jlist_n[iadr] = nb;
#ifdef HARD_DEBUG
            assert(nb_flag[iadr]==false);
            assert(nb<=_n_tot);
#endif
        }
    }
    
    class SortAdr{
    public:
        PS::F64 * time;
        SortAdr(PS::F64 * _time): time(_time){}
        bool operator() (const PS::S32 & left, const PS::S32 & right) const {
            return time[left] < time[right];
        }
    };

    void SortAndSelectIp(PS::S32 adr_dt_sorted[],
                         PS::F64 time_next[],
                         PS::S32 & n_act,
                         const PS::S32 n_tot,
                         PS::S32 group_list[],
                         PS::S32 & group_n,
                         const PS::F64 n_limit){
        // const PS::S32 n_tot = time_next.size();
        //std::cerr<<"before sort"<<std::endl;
        /*
          for(PS::S32 ip=0; ip<ni_old; ip++){
          const PS::S32 adr = adr_dt_sorted[ip];
          time_next[adr] += ptcl[adr].dt; // n_act only
          }
        */
        std::sort(adr_dt_sorted, adr_dt_sorted+n_act, SortAdr(time_next));

        const PS::F64 time_ref = time_next[adr_dt_sorted[0]];
        group_n = 0;
        for(n_act=1; n_act<n_tot; n_act++){
            if(time_ref < time_next[adr_dt_sorted[n_act]]) {
                break;
            }
            if(n_act<n_limit) group_list[group_n++] = n_act;
        }
    }

    void SortAndSelectIp(PS::S32 adr_dt_sorted[],
                         PS::F64 time_next[],
                         PS::S32 & n_act,
                         const PS::S32 n_tot){
        // const PS::S32 n_tot = time_next.size();
        //std::cerr<<"before sort"<<std::endl;
        /*
          for(PS::S32 ip=0; ip<ni_old; ip++){
          const PS::S32 adr = adr_dt_sorted[ip];
          time_next[adr] += ptcl[adr].dt; // n_act only
          }
        */
        std::sort(adr_dt_sorted, adr_dt_sorted+std::min(n_act,n_tot), SortAdr(time_next));

        const PS::F64 time_ref = time_next[adr_dt_sorted[0]];
        for(n_act=1; n_act<n_tot; n_act++){
            if(time_ref < time_next[adr_dt_sorted[n_act]]) {
                break;
            }
        }
    }
    
    void PredictAll(PtclPred pred[],
                    const PtclH4 ptcl[],
                    const PS::S32 n_tot,
                    const PS::F64 time_next){
        static thread_local const PS::F64 inv3 = 1.0 / 3.0;
        for(PS::S32 i=0; i<n_tot; i++){
            const PS::F64 dt = time_next - ptcl[i].time;
            pred[i].pos = ptcl[i].pos + dt*(ptcl[i].vel  + 0.5*dt*(ptcl[i].acc0 + inv3*dt*ptcl[i].acc1));
            pred[i].vel = ptcl[i].vel + dt*(ptcl[i].acc0 + 0.5*dt*ptcl[i].acc1);
            // pred[i].r_out = ptcl[i].r_out;
            pred[i].mass = ptcl[i].mass;
        }
    /*
      if(PS::Comm::getRank() == 0){
      for(PS::S32 i=0; i<n_tot; i++){
      std::cerr<<"pred[i].pos= "<<pred[i].pos
      <<" ptcl [i].pos="<<ptcl [i].pos<<std::endl;
      std::cerr<<"pred[i].vel= "<<pred[i].vel
      <<" ptcl [i].vel="<<ptcl [i].vel<<std::endl;
      }
      }
    */
    }

    //! correct ptcl and calculate step 
    /*! Correct ptcl and calculate next time step
      @param[in,out] _ptcl: particle array store the dt, time, previous acc0, acc1, the acc0, acc1 are updated with new ones
      @param[in] _force: predicted forces
      @param[in] _adr_dt_sorted: active particle index array
      @param[in] _n_act: number of active particles
      @param[in] _dt_max: maximum time step
      @param[in] _dt_min: minimum time step for check
      @param[in] _a0_offset_sq: acc0 offset to determine time step
      @param[in] _eta: time step criterion coefficiency
      \return fail_flag: if the time step < dt_min, return true (failure)
     */
    bool CorrectAndCalcDt4thAct(PtclH4 _ptcl[],
                                const PtclForce _force[],
                                const PS::S32 _adr_dt_sorted[], 
                                const PS::S32 _n_act,
                                const PS::F64 _dt_max,
                                const PS::F64 _dt_min,
                                const PS::F64 _a0_offset_sq,
                                const PS::F64 _eta){
        bool fail_flag=false;
        static thread_local const PS::F64 inv3 = 1.0 / 3.0;
        for(PS::S32 i=0; i<_n_act; i++){
            const PS::S32 adr = _adr_dt_sorted[i];
            PtclH4*     pti = &_ptcl[adr];
            const PtclForce* fpi = &_force[adr];

            const PS::F64 dt = pti->dt;
            const PS::F64 h = 0.5 * dt;
            const PS::F64 hinv = 2.0 / dt;
            const PS::F64vec A0p = (fpi->acc0 + pti->acc0);
            const PS::F64vec A0m = (fpi->acc0 - pti->acc0);
            const PS::F64vec A1p = (fpi->acc1 + pti->acc1)*h;
            const PS::F64vec A1m = (fpi->acc1 - pti->acc1)*h;

            const PS::F64vec vel_new = pti->vel + h*( A0p - inv3*A1m );
            pti->pos += h*( (pti->vel + vel_new) + h*(-inv3*A0m));
            pti->vel = vel_new;

            pti->acc0 = fpi->acc0;
            pti->acc1 = fpi->acc1;
            pti->time += dt;

            const PS::F64vec acc3 = (1.5*hinv*hinv*hinv) * (A1p - A0m);
            const PS::F64vec acc2 = (0.5*hinv*hinv) * A1m + h*acc3;
            const PS::F64 dt_ref = CalcDt4th(pti->acc0, pti->acc1, acc2, acc3, _eta, _a0_offset_sq);

            const PS::F64 dt_old = pti->dt;
#ifdef HARD_DEBUG
            // for debug
            assert(!std::isnan(pti->pos[0]));
            assert(!std::isnan(pti->pos[1]));
            assert(!std::isnan(pti->pos[2]));
            assert(!std::isnan(pti->vel[0]));
            assert(!std::isnan(pti->vel[1]));
            assert(!std::isnan(pti->vel[2]));
            //assert(pti->pos[0]==pti->pos[0]);
            //assert(pti->pos[1]==pti->pos[1]);
            //assert(pti->pos[2]==pti->pos[2]);
            //assert(pti->vel[0]==pti->vel[0]);
            //assert(pti->vel[1]==pti->vel[1]);
            //assert(pti->vel[2]==pti->vel[2]);
            assert(dt_old != 0.0);
#ifdef HARD_DEBUG_ACC
            pti->acc2 = acc2;
            pti->acc3 = acc3;
#endif
#endif
            pti->dt = _dt_max;

#ifdef FIX_STEP_DEBUG
            pti->dt /= STEP_DIVIDER;
#else
            while(pti->dt > dt_ref) pti->dt *= 0.5;
            pti->dt = dt_old*2 < pti->dt ?  dt_old*2 : pti->dt;
//            if(int(pti->time/pti->dt)*pti->dt<pti->time) pti->dt *= 0.5;
#endif

#ifdef HARD_DEBUG
            assert(pti->dt != 0.0);
//            assert(pti->dt >1.0e-12);
#endif

            if(pti->dt <_dt_min) {
                std::cerr<<"Error: Hermite integrator step size ("<<pti->dt<<") < dt_min ("<<_dt_min<<")!"<<std::endl;
                std::cerr<<" pti->time="<<pti->time<<" i="<<i<<" adr="<<adr<<" pos="<<pti->pos<<" vel="<<pti->vel<<" acc="<<pti->acc0<<" acc1="<<pti->acc1
#ifdef HARD_DEBUG_ACC
                         <<" acc2="<<pti->acc2
                         <<" acc3="<<pti->acc3
#endif
                         <<std::endl;
                fail_flag=true;
            }
        }
        return fail_flag;
    }



    template<class Teng>
    void CalcEnergy(const PtclH4 ptcl[], const PS::S32 n_tot, Teng & eng, 
                    const PS::F64 r_in, const PS::F64 r_out, const PS::F64 eps_sq = 0.0){
        eng.kin = eng.pot = eng.tot = 0.0;
#ifndef INTEGRATED_CUTOFF_FUNCTION
        PS::F64 r_oi_inv = 1.0/(r_out-r_in);
        PS::F64 r_A = (r_out-r_in)/(r_out+r_in);
        PS::F64 pot_off = cutoff_pot(1.0, r_oi_inv, r_A, r_in)/r_out;
#endif
        for(PS::S32 i=0; i<n_tot; i++){
            eng.kin += 0.5 * ptcl[i].mass * ptcl[i].vel * ptcl[i].vel;

            for(PS::S32 j=i+1; j<n_tot; j++){
                //PS::F64 r_out = std::max(ptcl[i].r_out,ptcl[j].r_out);
                PS::F64vec rij = ptcl[i].pos - ptcl[j].pos;
                PS::F64 dr = sqrt(rij*rij + eps_sq);
#ifdef INTEGRATED_CUTOFF_FUNCTION 
                PS::F64 k  = 1- CalcW(dr/r_out, r_in/r_out);
                eng.pot -= ptcl[j].mass*ptcl[i].mass/dr*k;
#else
                PS::F64 k  = cutoff_pot(dr, r_oi_inv, r_A, r_in);
                if(dr<r_out) eng.pot -= ptcl[j].mass*ptcl[i].mass*(1.0/dr*k - pot_off);
#endif
            }
        }
        eng.tot = eng.kin + eng.pot;
    }

    //! modify list based on the trace table
    /*! modify list based on the trace table
      @param[out] _list: list to modify
      @param[in] _n_list: list size
      @param[in] _trace: moving table for _list, <0: no change, ==_del_key: remove, others: new position in _list
      @param[in] _del_key: indicator for removing
      \return new list size
     */
    PS::S32 modifyList(PS::S32* _list, 
                    const PS::S32 _n_list,
                    const PS::S32* _trace,
                    const PS::S32 _del_key) {
        if(_n_list==0) return 0;
#ifdef HARD_DEBUG
        assert(_n_list>0);
#endif
        PS::S32 list_size=_n_list;
        // moving offset
        PS::S32 imove=0;
        // index 
        PS::S32 i=0;
        while(i<list_size) {
            PS::S32 iadr=_list[i];
            PS::S32 itrace=_trace[iadr];
            // delete case
            if(itrace==_del_key) {
                imove++;
                list_size--;
            }
            else {
                // replace case
                if(itrace>=0) _list[i] = itrace;
                // move to next
                i++;
            }
            if(i>=list_size) break;
#ifdef HARD_DEBUG
            assert(i+imove<_n_list);
#endif
            _list[i] = _list[i+imove];
        }

        return list_size;
    }


public:

//    HermiteIntegrator(const PS::S32 n) {
//#ifdef HARD_DEBUG
//        assert(pred_.size()==0);
//        assert(force_.size()==0);
//        assert(adr_dt_sorted_.size()==0);
//        assert(time_next_.size()==0);
//#endif
//        resizeArray(n);
//    }
    
    //! get whether two pair is leaving each other or coming
    /*!
      @param[in] _i: first index in ptcl_
      @param[in] _j: second index in ptcl_
      \return true: go away, otherwise income
    */
    bool getDirection(const PS::S32 _i, const PS::S32 _j) {
#ifdef HARD_DEBUG
        assert(_i>=0&&_i<ptcl_.size());
        assert(_j>=0&&_j<ptcl_.size());
#endif                
        PS::F64 xv2=(ptcl_[_i].pos-ptcl_[_j].pos)*(ptcl_[_i].vel-ptcl_[_j].vel);
        return (xv2>0);
    }


    //! check the pair with distance below r_crit for ptcl in adr_dt_sorted_
    /*! First check nearest neighbor distance r_min
        If r_min<r_crit, check the direction, if income, accept as group
      @param[out] _new_group_member_index: new group member index in ptcl_
      @param[out] _new_group_member_adr: new group member address, echo two neighbor is one group
      @param[out] _new_group_offset: new group offset, last one show total number of members.
      @param[in] _r_crit2: group distance criterion square
      @param[in] _Aint: ARC integrator class
      \return number of new groups
     */
    template <class ARCint>
    PS::S32 checkNewGroup(PS::S32 _new_group_member_index[], Tphard* _new_group_member_adr[], PS::S32 _new_group_offset[], const PS::F64 _r_crit2, const ARCint* _Aint) {
#ifdef HARD_DEBUG
        assert(nb_info_.size()==ptcl_.size());
#endif
        PS::S32 n_group = 0;
        if (_Aint!=NULL) n_group = _Aint->getNGroups();
        PS::S32 n_new_group=0, offset=0;

        // here ptcl_.size() is used since any ptcl index can be checked!
        PS::S32 used_mask[ptcl_.size()];
        for (PS::S32 k=0; k<ptcl_.size(); k++) used_mask[k] = -1;

        // check adr_dt_sorted_ list (avoid suppressed ptcl)
        PS::S32 n_check = adr_dt_sorted_.size();
        for (PS::S32 k=0; k<n_check; k++) {
            const PS::S32 i = adr_dt_sorted_[k];
            if(nb_info_[i].r_min2<_r_crit2) {
                const PS::S32 j = nb_info_[i].r_min_index;
#ifdef HARD_DEBUG
                assert(j<ptcl_.size());
                assert(ptcl_[i].mass>0.0);
#endif                
                if(j<0) continue;

                if(!(used_mask[i]>=0 && used_mask[j]>=0)) { // avoid double count
                    bool out_flag=getDirection(i, j);
                    if(!out_flag) {
                        PS::F64 sdi=0.0, sdj=0.0;
                        if(i<n_group) sdi = _Aint->getSlowDown(i);
                        if(j<n_group) sdj = _Aint->getSlowDown(j);
                        if(sdi>1.0&&sdj>1.0) continue;
                        if(sdi>1.0&&sdj==0.0) continue;
                        if(sdi==0.0&&sdj>1.0) continue;
                        // to avoid extremely long time integration, for very large slowdown factor, no merge group
                        // current 100.0 is experimental value
                        if(sdi>100.0||sdj>100.0) continue;

                        PS::S32 insert_group=-1, insert_index=-1;
                        if(used_mask[i]>=0) {
                            insert_group = used_mask[i];
                            insert_index = j;
                        }
                        else if(used_mask[j]>=0) {
                            insert_group = used_mask[j];
                            insert_index = i;
                        }
                        // the case of merging group
                        if(insert_group>=0) {
                            // shift the first index in last group to last 
                            PS::S32 last_group_offset = _new_group_offset[n_new_group-1];
                            _new_group_member_index[offset] = _new_group_member_index[last_group_offset];
                            _new_group_member_adr  [offset] = _new_group_member_adr  [last_group_offset];
                            offset++;

                            // in the case insert_group is not the last group
                            if(insert_group<n_new_group-1) {
                                // shift first index in each group to the end to allow to insert j in i_group
                                for (PS::S32 k=n_new_group-1; k>insert_group; k--) {
                                    PS::S32 k_group_offset = _new_group_offset[k];
                                    PS::S32 k0_group_offset = _new_group_offset[k-1];
                                    _new_group_member_index[k_group_offset] = _new_group_member_index[k0_group_offset];
                                    _new_group_member_adr  [k_group_offset] = _new_group_member_adr  [k0_group_offset];
                                    _new_group_offset[k]++;
                                }
                            }
                            // replace the first position of insert_group with insert ptcl
                            PS::S32 insert_group_offset=_new_group_offset[insert_group];
                            _new_group_member_index[insert_group_offset] = insert_index;
                            _new_group_member_adr[insert_group_offset] = ptcl_ptr_[insert_index];
                            used_mask[insert_index] = insert_group;
                        }
                        else {   // new group case
                            _new_group_offset[n_new_group] = offset;
                            _new_group_member_index[offset] = i;
                            _new_group_member_adr[offset] = ptcl_ptr_[i];
                            used_mask[i] = n_new_group;
                            offset++;

                            _new_group_member_index[offset] = j;
                            _new_group_member_adr[offset] = ptcl_ptr_[j];
                            used_mask[j] = n_new_group;
                            offset++;

                            n_new_group++;
                        }
                    }
                }
            }
        }
        // for total number of members
        _new_group_offset[n_new_group] = offset;
#ifdef HARD_DEBUG
        assert(offset<=ptcl_.size());
#endif
        return n_new_group;
    }
    
    void moveCM(const PS::F64 dt) {
        pcm_.pos += pcm_.vel * dt;
        //std::cerr<<"pcm.pos "<<pcm_.pos<<" pcm.vel "<<pcm_.vel<<" dt "<<dt<<std::endl;
    }

    //!shift ptcl back to original frame
    void shiftBackCM() {
        shiftBackCM(pcm_,ptcl_.getPointer(),ptcl_.size());
    }

    //!shift ptcl to c.m. frame and save c.m.
    void shiftToCM() {
        calcCMAndShift(pcm_, ptcl_.getPointer(), ptcl_.size());
    }

    //! Reserve memory
    /*! To reserve memory for use
      @param[in] _n: maximum particle number to reserve
     */
    void reserveMem(const PS::S32 _n) {
        ptcl_.reserve(_n);
        ptcl_ptr_.reserve(_n);
        Jlist_.reserve(_n*n_nb_off_);
        Jlist_disp_.reserve(_n);
        Jlist_n_.reserve(_n);
        nb_info_.reserve(_n);
        pred_.reserve(_n);
        force_.reserve(_n);
        adr_dt_sorted_.reserve(_n);
        time_next_.reserve(_n);
    }

    ////! Copy particle data to local array
    ///* @param[in] _ptcl: particle data array
    //   @param[in] _n_ptcl: number of particle need to be added
    //   @param[in] _ptcl_list: particle address in _ptcl
    // */
    //template <class Tptcl>
    //void addInitPtcl(Tptcl * _ptcl, 
    //             const PS::S32 _n_ptcl, 
    //             const PS::S32* _ptcl_list) {
    //    for (PS::S32 i=0; i<_n_ptcl; i++) {
    //        ptcl_.increaseSize(1);
    //        ptcl_.back().DataCopy(_ptcl[_ptcl_list[i]]);
    //        //ptcl_.pushBackNoCheck(ptcl_org[ptcl_list[i]]);
    //    }
    //}
    // 
    ////! Copy particle data to local array
    ///* @param[in] _ptcl: particle data array
    //   @param[in] _n_ptcl: number of particles need to be added
    // */
    //template <class Tptcl>
    //void addInitPtcl(Tptcl * _ptcl, 
    //             const PS::S32 _n_ptcl) {
    //    for (PS::S32 i=0; i<_n_ptcl; i++) {
    //        ptcl_.increaseSize(1);
    //        ptcl_.back().DataCopy(_ptcl[i]);
    //    }
    //}

    //! Insert one particle
    /*! insert one particle to local array at position i, notice i is added in the front of adr_dt_sorted_.
      @param[in] _i: insert position
      @param[in] _ptcl: particle to be inserted
      @param[in] _n_list_correct: number of ptcl (count from first) need to correct the neighbor list
      @param[in] _time_sys: new ptcl time
     */
    template <class Tptcl>
    void insertOnePtcl(const PS::S32 _i,
                       Tptcl &_ptcl,
                       const PS::S32 _n_list_correct,
                       const PS::F64 _time_sys) {
#ifdef HARD_DEBUG
        assert(_i<=ptcl_.size());
        assert(ptcl_.size()+1<=n_nb_off_);
#endif 
        // increase all data size by 1
        ptcl_.increaseSize(1);
        ptcl_ptr_.increaseSize(1);
        pred_.increaseSize(1);
        force_.increaseSize(1);
        adr_dt_sorted_.increaseSize(1);
        time_next_.increaseSize(1);
        Jlist_disp_.push_back(Jlist_.size());
        Jlist_.increaseSize(n_nb_off_);
        Jlist_n_.increaseSize(1);
        nb_info_.increaseSize(1);

#ifdef HARD_DEBUG
        assert(ptcl_.size()==pred_.size());
        assert(ptcl_.size()==force_.size());
        assert(ptcl_.size()==time_next_.size());
        assert(ptcl_.size()==Jlist_disp_.size());
        assert(ptcl_.size()==Jlist_n_.size());
        assert(ptcl_.size()==nb_info_.size());
#endif

        // last just add
        if (_i==ptcl_.size()-1) {
            ptcl_.back().DataCopy(_ptcl);
            ptcl_.back().acc0 = ptcl_.back().acc1 = force_.back().acc0 = force_.back().acc1 = 0.0;
            ptcl_.back().time = _time_sys;
            ptcl_.back().dt = 0.0;
            ptcl_ptr_.back()=(Tphard*)&_ptcl;

            // shift adr_dt_sorted
            for (PS::S32 i=ptcl_.size()-1; i>0; i--) {
                adr_dt_sorted_[i] = adr_dt_sorted_[i-1];
            }
            adr_dt_sorted_[0] = _i;
        }
        // shift _i to last and add
        else {
            PS::S32 inew = ptcl_.size()-1;
#ifdef HARD_DEBUG
            assert(&(ptcl_.back())==&(ptcl_[inew]));
#endif
            ptcl_.back() = ptcl_[_i];
            ptcl_ptr_.back() = ptcl_ptr_[_i];
            ptcl_[_i].DataCopy(_ptcl);
            ptcl_[_i].acc0 = ptcl_[_i].acc1 = force_[_i].acc0 = force_[_i].acc1 = 0.0;
            ptcl_[_i].time = _time_sys;
            ptcl_[_i].dt = 0.0;
            ptcl_ptr_[_i]=(Tphard*)&_ptcl;
            
            pred_.back()       = pred_[_i];
            force_.back()      = force_[_i];
            time_next_.back()  = time_next_[_i];
            nb_info_.back()    = nb_info_[_i];

            // find sorted index and change
            for (PS::S32 i=inew; i>0; i--) {
                adr_dt_sorted_[i] = adr_dt_sorted_[i-1];
                if(adr_dt_sorted_[i]==_i) {
                    adr_dt_sorted_[i] = inew;
                }
            }
            adr_dt_sorted_[0] = _i;
            n_act_++;

            // shift Jlist
            Jlist_n_.back()    = Jlist_n_[_i];
            Jlist_n_[_i] = 0;
            PS::S32 disp_last = Jlist_disp_.back();
            PS::S32 disp_i = Jlist_disp_[_i];
            for (PS::S32 i=0; i<Jlist_n_.back(); i++) {
                Jlist_[disp_last+i] = Jlist_[disp_i+i];
            }

            // check nb_info
            for (PS::S32 i=0; i<=inew; i++) {
                if(nb_info_[i].r_min_index==_i) 
                    nb_info_[i].r_min_index=inew;
            }

            // correct neighbor list
            for (PS::S32 i=0; i<_n_list_correct; i++) {
                disp_i = Jlist_disp_[i];
                for (PS::S32 j=0; j<Jlist_n_[i]; j++) {
                    PS::S32 k = disp_i+j;
                    if(Jlist_[k]==_i) Jlist_[k]=inew;
                }
            }
        }
    }

    //! Replace one particle
    /*! Replace one particle of _i and shift the index in adr_dt_softed_ to front
      @param[in] _i: modified index in Hint.ptcl_
      @param[in] _ptcl: new ptcl information
      @param[in] _time_sys: new ptcl time
    */
    template <class Tptcl>
    void modOnePtcl(const PS::S32 _i,
                    const Tptcl &_ptcl,
                    const PS::F64 _time_sys ) {
#ifdef HARD_DEBUG
        assert(_i<=ptcl_.size());
#endif 
        ptcl_[_i].DataCopy(_ptcl);
        ptcl_[_i].acc0 = ptcl_[_i].acc1 = force_[_i].acc0 = force_[_i].acc1 = 0;
        ptcl_[_i].time = _time_sys;
        ptcl_[_i].dt = 0.0;
        ptcl_ptr_[_i] = (Tphard*)&_ptcl;

        // find index and shift to front
        const PS::S32 adr_size=adr_dt_sorted_.size();

        PS::S32 i=0;
        while(adr_dt_sorted_[i]!=_i&&i<adr_size) i++;
        if(i>=adr_size) adr_dt_sorted_.increaseSize(1);

        // shift
        for (PS::S32 k=i; k>0; k--) 
            adr_dt_sorted_[k] = adr_dt_sorted_[k-1];
        adr_dt_sorted_[0] = _i;
        if(i>=n_act_) n_act_++;
    }


    //! Add a list of particles at the end
    /*! Add a list of particles at the end of Hint.ptcl_, notice all new particles have index in front of adr_dt_sorted_.
      @param[in] _ptcl_origin: particle array to be added
      @param[in] _list_origin: adding particle index in _ptcl_origin
      @param[in] _n_list: number of adding particles
      @param[in] _n_nb_correct: number of ptcl (count from first): need to correct the neighbor list
      @param[in] _time_sys: new ptcl time
      @param[in] _adr_dt_front_flag: add new ptcl index in front of adr_dt_sorted_ (true) or end (false)
     */
    template <class Tptcl>
    void addPtclList(Tptcl* _ptcl_origin,
                     const PS::S32* _list_origin,
                     const PS::S32 _n_list,
                     const PS::S32 _n_nb_correct,
                     const PS::F64 _time_sys,
                     const bool _adr_dt_front_flag) {

        // quite if no new
        if(_n_list==0) return;
#ifdef HARD_DEBUG
        assert(_n_list>0);
#endif

        // original size
        const PS::S32 n_org = ptcl_.size();
        const PS::S32 n_adr_dt_org = adr_dt_sorted_.size();

        // increase all data size by _n_list
        ptcl_.increaseSize(_n_list);
        ptcl_ptr_.increaseSize(_n_list);
        pred_.increaseSize(_n_list);
        force_.increaseSize(_n_list);
        adr_dt_sorted_.increaseSize(_n_list);
        time_next_.increaseSize(_n_list);
        Jlist_n_.increaseSize(_n_list);
        Jlist_disp_.increaseSize(_n_list);
        Jlist_.increaseSize(_n_list*n_nb_off_);
        nb_info_.increaseSize(_n_list);

#ifdef HARD_DEBUG
        assert(ptcl_.size()<ARRAY_ALLOW_LIMIT);
        assert(ptcl_.size()==ptcl_ptr_.size());
        assert(ptcl_.size()==pred_.size());
        assert(ptcl_.size()==force_.size());
        assert(ptcl_.size()==time_next_.size());
        assert(ptcl_.size()==Jlist_disp_.size());
        assert(ptcl_.size()==Jlist_n_.size());
        assert(ptcl_.size()==nb_info_.size());
        assert(ptcl_.size()*n_nb_off_==Jlist_.size());
#endif

        // add ptcl in order
        if (_list_origin==NULL) {
            for (PS::S32 i=0; i<_n_list; i++) {
                const PS::S32 inew= n_org+i;
                ptcl_[inew].DataCopy(_ptcl_origin[i]);
                ptcl_[inew].time = _time_sys;
                ptcl_[inew].acc0 = ptcl_[inew].acc1 = force_[inew].acc0 = force_[inew].acc1 = 0.0;
                ptcl_[inew].dt = 0.0;
                ptcl_ptr_[inew] = (Tphard*)&_ptcl_origin[i];
                Jlist_disp_[inew] = (n_org+i)*n_nb_off_;
                //_single_index_origin[_n_single+i] = i;
            }
        }
        // add ptcl from list
        else {
            for (PS::S32 i=0; i<_n_list; i++) {
                const PS::S32 iadr= _list_origin[i];
                const PS::S32 inew= n_org+i;
                ptcl_[inew].DataCopy(_ptcl_origin[iadr]);
                ptcl_[inew].time = _time_sys;
                ptcl_[inew].acc0 = ptcl_[inew].acc1 = force_[inew].acc0 = force_[inew].acc1 = 0.0;
                ptcl_[inew].dt = 0.0;
                ptcl_ptr_[inew] = (Tphard*)&_ptcl_origin[iadr];
                Jlist_disp_[inew] = (n_org+i)*n_nb_off_;
                //_single_index_origin[_n_single+i] = iadr;
            }
        }
        // jlist disp

        //_n_single += _n_list;
        
        // add new ptcl in front of adr_dt_sorted
        if(_adr_dt_front_flag) {
            // shift adr_dt_sorted
            for (PS::S32 i=n_adr_dt_org-1; i>=0; i--) 
                adr_dt_sorted_[i+_n_list] = adr_dt_sorted_[i];
            // add new ptcl to adr_dt_sorted
            for (PS::S32 i=0; i<_n_list; i++) 
                adr_dt_sorted_[i] = n_org+i;
        }
        // add at then end of adr_dt_sorted
        else {
            for (PS::S32 i=n_adr_dt_org; i<n_adr_dt_org+_n_list; i++) 
                adr_dt_sorted_[i] = i;
        }
        n_act_ += _n_list;

        // add new ptcl to neighbor list of c.m.
        for (PS::S32 i=0; i<_n_nb_correct; i++) { 
            // escape the suppressed case
            if(ptcl_[i].mass==0&&Jlist_n_[i]==0) continue;
            const PS::S32 ioff=Jlist_disp_[i]+Jlist_n_[i];
            for (PS::S32 j=ioff; j<ioff+_n_list; j++) 
                Jlist_[j] = n_org+j-ioff;
            Jlist_n_[i] += _n_list;
#ifdef HARD_DEBUG
            assert(Jlist_n_[i]<=n_nb_off_);
#endif
        }
    }
                    
    //! Remove a list of particle
    /*! Remove particles from _list, keep first _n_correct in ptcl_ undeleted (but mass to zero) and update their neighbor lists
      @param[in] _list: removing particle index for ptcl_
      @param[in] _n_list: number of removing particles
      @param[in] _n_group: number of c.m. ptcl (count from first): ptcl will not be delected but only remove from adr_dt_sorted_, also the offset of the index between Hint.ptcl_ _single_index_origin
      @param[in] _n_nb_correct: number of ptcl (count from first) to correct neighbor list
     */
    void removePtclList(const PS::S32* _list,
                        const PS::S32 _n_list,
                        const PS::S32 _n_group,
                        const PS::S32 _n_nb_correct) {
        // quit if no deleted ones
        if(_n_list==0) return;
#ifdef HARD_DEBUG
        assert(_n_list>0);
#endif

        const PS::S32 n_org=ptcl_.size();
        // moving trace (initial -1)
        PS::S32 trace[n_org];
        for (PS::S32 i=0; i<n_org; i++) trace[i] = -1;
        PS::S32 ilast = n_org-1;

        // create moving table (trace)
        // record the new position if the ptcl is moved, if deleted, set to ptcl_.size()
        PS::S32 n_decrease=0;
        for (PS::S32 i=0; i<_n_list; i++) {
            PS::S32 k=_list[i];
            // check no delete case
            if(k<_n_group) {
                trace[k]=n_org;
                // set suppressed group mass to zero
                ptcl_[k].mass = 0.0;
                // set status to -20 to identify the suppressed ptcl
                ptcl_[k].status = -20;
                // set neighbor to zero to avoid issue
                Jlist_n_[k] = 0;
                // set r_min_index to -1
                nb_info_[k].r_min_index = -1;
                continue;
            }

            PS::S32 idel = k;
            // set ilast >= idel for safety.
            ilast = std::max(ilast, idel);
            
            // check whether position is already moved, if so, check the moved new position and set current to delete
            if(trace[k]>=0) {
                idel=trace[k];
                trace[k]=n_org;
            }
#ifdef HARD_DEBUG
            assert(k<n_org);
            assert(idel<n_org);
#endif

            // check the last avaiable particle that can be moved to idel
            // If the ilast is already moved, check whether the new moved position trace[ilast] is before the current idel.
            // If trace[ilast] is after the current idel, update the trace[ilast] to idel and set trace[ilast] as new idel to check, until trace[ilast]=-1 or idel >= ilast
            while (idel>=0) {
                PS::S32 itrlast = -1;
                while(trace[ilast]>=0 && (trace[ilast]<idel || trace[ilast]==n_org) && ilast>idel) ilast--;

                // if ilast > idel, move ilast to idel
                if(idel<ilast) {
                    itrlast=trace[ilast];
                    trace[ilast]=idel;
                }
                // idel is already at last, remove it
                trace[idel]=n_org; 
                idel = itrlast;
            }

            n_decrease++;
        }

        // move data
        for (PS::S32 i=0; i<n_org; i++) {
            // no change case
            if(trace[i]<0) continue;
            // remove case
            else if(trace[i]<n_org) {
                const PS::S32 inew=trace[i];
                ptcl_[inew] = ptcl_[i];
                ptcl_ptr_[inew] = ptcl_ptr_[i];
                pred_[inew] = pred_[i];
                force_[inew] = force_[i];
                time_next_[inew] = time_next_[i];
                nb_info_[inew] = nb_info_[i];
                //_single_index_origin[inew-_n_group] = _single_index_origin[i-_n_group];

                // shift Jlist
                Jlist_n_[inew] = Jlist_n_[i];
                Jlist_n_[i] = 0;
                PS::S32 disp_i = Jlist_disp_[i];
                PS::S32 disp_inew = Jlist_disp_[inew];
                for (PS::S32 i=0; i<Jlist_n_[inew]; i++) 
                    Jlist_[disp_inew+i] = Jlist_[disp_i+i];
            }
        }
        ptcl_.decreaseSize(n_decrease);
        ptcl_ptr_.decreaseSize(n_decrease);
        pred_.decreaseSize(n_decrease);
        force_.decreaseSize(n_decrease);
        time_next_.decreaseSize(n_decrease);
        nb_info_.decreaseSize(n_decrease);
        Jlist_n_.decreaseSize(n_decrease);
        Jlist_disp_.decreaseSize(n_decrease);
        Jlist_.decreaseSize(n_decrease*n_nb_off_);
        //_n_single -= n_decrease;

        // check nb_info
        for (PS::S32 i=0; i<ptcl_.size(); i++) {
            const PS::S32 i_min=nb_info_[i].r_min_index;
            if(i_min>=0) {
                const PS::S32 jtr=trace[i_min];
                if(jtr>=0) {
                    if(jtr<n_org) nb_info_[i].r_min_index=jtr;
                    else nb_info_[i].r_min_index=-1;
                }
            }
        }

        // update adr_dt_sorted
        PS::S32 adr_dt_sorted_new_size=modifyList(adr_dt_sorted_.getPointer(), adr_dt_sorted_.size(), trace, n_org);
        adr_dt_sorted_.decreaseSize(_n_list);

        assert(adr_dt_sorted_new_size==adr_dt_sorted_.size());

        
        // correct neighbor list
        for (PS::S32 i=0; i<_n_nb_correct; i++) {
            // escape suppressed one
            if (trace[i]==n_org) continue;
            
            PS::S32 jlist_i_new_size=modifyList(Jlist_.getPointer(Jlist_disp_[i]), Jlist_n_[i], trace, n_org);
            Jlist_n_[i] = jlist_i_new_size;
#ifdef HARD_DEBUG
            assert(Jlist_n_[i]>=0);
#endif
        }
            
    }                       

    //! Write back particle data to ptcl
    /* @param[out] _ptcl: ptcl array to write
       @param[in] _ptcl_list: particle address in _ptcl
       @param[in] _n_ptcl: number of particles need for copy
       @param[in] _i_start: start index in local particle array for copy
     */
    template <class Tptcl>
    void writeBackPtcl(Tptcl * _ptcl, 
                       const PS::S32* _ptcl_list,
                       const PS::S32 _n_ptcl, 
                       const PS::S32 _i_start) {
#ifdef HARD_DEBUG
        assert(_i_start+_n_ptcl<=ptcl_.size());
#endif
        for (PS::S32 i=0; i<_n_ptcl; i++) {
#ifdef HARD_DEBUG
            assert(_ptcl[_ptcl_list[i]].id==ptcl_[i+_i_start].id);
            assert(!std::isnan(ptcl_[i+_i_start].pos[0]));
            assert(!std::isnan(ptcl_[i+_i_start].pos[1]));
            assert(!std::isnan(ptcl_[i+_i_start].pos[2]));
            assert(!std::isnan(ptcl_[i+_i_start].vel[0]));
            assert(!std::isnan(ptcl_[i+_i_start].vel[1]));
            assert(!std::isnan(ptcl_[i+_i_start].vel[2]));
#endif
            _ptcl[_ptcl_list[i]].DataCopy(ptcl_[i+_i_start]);
        }
    }

    //! Write back particle data to original array
    /*!
      @param[in] _i_start: starting index in ptcl_ to copy
    */
    void writeBackPtcl(const PS::S32 _i_start) {
        for (PS::S32 i=_i_start; i<ptcl_.size(); i++) {
#ifdef HARD_DEBUG
            assert(!std::isnan(ptcl_[i].pos[0]));
            assert(!std::isnan(ptcl_[i].pos[1]));
            assert(!std::isnan(ptcl_[i].pos[2]));
            assert(!std::isnan(ptcl_[i].vel[0]));
            assert(!std::isnan(ptcl_[i].vel[1]));
            assert(!std::isnan(ptcl_[i].vel[2]));
//            assert(ptcl_[i].pos[0]==ptcl_[i].pos[0]);
//            assert(ptcl_[i].pos[1]==ptcl_[i].pos[1]);
//            assert(ptcl_[i].pos[2]==ptcl_[i].pos[2]);
//            assert(ptcl_[i].vel[0]==ptcl_[i].vel[0]);
//            assert(ptcl_[i].vel[1]==ptcl_[i].vel[1]);
//            assert(ptcl_[i].vel[2]==ptcl_[i].vel[2]);
            assert(ptcl_ptr_[i]->id==ptcl_[i].id);
#endif
            ptcl_ptr_[i]->DataCopy(ptcl_[i]);
        }
    }

    //! Write back a list of particle data to original array 
    /*!
      @param[in] _ptcl_list: particle index in ptcl_
      @param[in] _n_ptcl: number of particles need for copy
      @param[in] _n_avoid: if particle index is below _n_avoid, no write (to avoid write group c.m.)
    */
    void writeBackPtcl(const PS::S32* _ptcl_list,
                       const PS::S32 _n_ptcl,
                       const PS::S32 _n_avoid) {
        for (PS::S32 i=0; i<_n_ptcl; i++) {
            const PS::S32 k = _ptcl_list[i];
            if(k<_n_avoid) continue;
#ifdef HARD_DEBUG
            assert(k<ptcl_.size()&&k>=0);
            assert(!std::isnan(ptcl_[k].pos[0]));
            assert(!std::isnan(ptcl_[k].pos[1]));
            assert(!std::isnan(ptcl_[k].pos[2]));
            assert(!std::isnan(ptcl_[k].vel[0]));
            assert(!std::isnan(ptcl_[k].vel[1]));
            assert(!std::isnan(ptcl_[k].vel[2]));
//            assert(ptcl_[k].pos[0]==ptcl_[k].pos[0]);
//            assert(ptcl_[k].pos[1]==ptcl_[k].pos[1]);
//            assert(ptcl_[k].pos[2]==ptcl_[k].pos[2]);
//            assert(ptcl_[k].vel[0]==ptcl_[k].vel[0]);
//            assert(ptcl_[k].vel[1]==ptcl_[k].vel[1]);
//            assert(ptcl_[k].vel[2]==ptcl_[k].vel[2]);
            assert(ptcl_ptr_[k]->id==ptcl_[k].id);
#endif            
            ptcl_ptr_[k]->DataCopy(ptcl_[k]);
        }
    }

    //! Write back one particle data to original array 
    /*!
      @param[in] _i: particle index in ptcl_
    */
    void writeBackOnePtcl(const PS::S32 _i) {
#ifdef HARD_DEBUG
        assert(_i<ptcl_.size()&&_i>=0);
        assert(!std::isnan(ptcl_[_i].pos[0]));
        assert(!std::isnan(ptcl_[_i].pos[1]));
        assert(!std::isnan(ptcl_[_i].pos[2]));
        assert(!std::isnan(ptcl_[_i].vel[0]));
        assert(!std::isnan(ptcl_[_i].vel[1]));
        assert(!std::isnan(ptcl_[_i].vel[2]));
//        assert(ptcl_[_i].pos[0]==ptcl_[_i].pos[0]);
//        assert(ptcl_[_i].pos[1]==ptcl_[_i].pos[1]);
//        assert(ptcl_[_i].pos[2]==ptcl_[_i].pos[2]);
//        assert(ptcl_[_i].vel[0]==ptcl_[_i].vel[0]);
//        assert(ptcl_[_i].vel[1]==ptcl_[_i].vel[1]);
//        assert(ptcl_[_i].vel[2]==ptcl_[_i].vel[2]);
        assert(ptcl_ptr_[_i]->id==ptcl_[_i].id);
#endif            
        ptcl_ptr_[_i]->DataCopy(ptcl_[_i]);
    }

    //! Search perturber and neighbor
    /*! Assume the suppressed ptcl has zero mass
      @param[in] _apo_bin: apo-center distance of binary 
      @param[in] _n_bin: number of binaries (locate at begining of ptcl_)
      @param[in] _n_search: number of particles to search perturber
     */
    void searchPerturberBin(PS::F64 _apo_bin[], const PS::S32 _n_bin, const PS::S32 _n_search) {
        PS::S32 n = ptcl_.size();
        
        // find perturber
        for(PS::S32 i=0; i<_n_search; i++) {
            PS::S32 disp=Jlist_disp_[i];
            PS::S32 n_pert=0;
            nb_info_[i].r_min2 = PS::LARGE_FLOAT;
            nb_info_[i].min_mass = PS::LARGE_FLOAT;
            for(PS::S32 j=0; j<_n_bin; j++) {
                PS::F64vec dr = ptcl_[i].pos-ptcl_[j].pos;
                PS::F64 r2 = dr*dr;
                PS::F64 r_search = std::max(ptcl_[i].r_search,ptcl_[j].r_search) + _apo_bin[j];
                PS::F64 r_search2 = r_search*r_search;
                if (r2<r_search2&&i!=j&&ptcl_[j].mass>0.0) {
                    Jlist_[disp+n_pert]=j;
                    n_pert++;
                    if(r2<nb_info_[i].r_min2) {
                        nb_info_[i].r_min2 = r2;
                        nb_info_[i].r_min_index = j;
                    }
                    nb_info_[i].min_mass = std::min(nb_info_[i].min_mass, ptcl_[j].mass);
                }
            }
            for(PS::S32 j=_n_bin; j<n; j++) {
                PS::F64vec dr = ptcl_[i].pos-ptcl_[j].pos;
                PS::F64 r2 = dr*dr;
                PS::F64 r_search = std::max(ptcl_[i].r_search,ptcl_[j].r_search);
                PS::F64 r_search2 = r_search*r_search;
                if (r2<r_search2&&i!=j&&ptcl_[j].mass>0.0) {
                    Jlist_[disp+n_pert]=j;
                    n_pert++;
                    if(r2<nb_info_[i].r_min2) {
                        nb_info_[i].r_min2 = r2;
                        nb_info_[i].r_min_index = j;
                    }
                    nb_info_[i].min_mass = std::min(nb_info_[i].min_mass, ptcl_[j].mass);
                }
            }
//            Jlist_disp_[i] = disp;
            Jlist_n_[i] = n_pert;
#ifdef HARD_DEBUG
            assert(n_pert<=n_nb_off_);
#endif
        }
    }

    //! Search perturber and neighbor
    /*!
      @param[in] _n_search: number of particles to search perturber
     */
    void searchPerturber(const PS::S32 _n_search) {
        for(PS::S32 i=0; i<_n_search; i++) searchPerturberOne(i);
    }

    //! Search perturber and neighbor for one 
    /*! Search perturber and neighbor for one from adr_dt_sorted list.
        In this case, the suppressed particle can be avoid
      @param[in] _i: index of particle to search perturber
     */
    inline void searchPerturberOne(const PS::S32 _i) {
        PS::S32 n = adr_dt_sorted_.size();
        // find perturber
        PS::S32 disp=Jlist_disp_[_i];
        PS::S32 n_pert=0;
        nb_info_[_i].r_min2 = PS::LARGE_FLOAT;
        nb_info_[_i].min_mass = PS::LARGE_FLOAT;
        for(PS::S32 k=0; k<n; k++) {
            PS::S32 j = adr_dt_sorted_[k];
            PS::F64vec dr = ptcl_[_i].pos-ptcl_[j].pos;
            PS::F64 r2 = dr*dr;
            PS::F64 r_search = std::max(ptcl_[_i].r_search,ptcl_[j].r_search);
            PS::F64 r_search2 = r_search*r_search;
            if (r2<r_search2&&_i!=j) {
#ifdef HARD_DEBUG
                assert(ptcl_[j].mass>0);
#endif
                Jlist_[disp+n_pert]=j;
                n_pert++;
                if(r2<nb_info_[_i].r_min2) {
                    nb_info_[_i].r_min2 = r2;
                    nb_info_[_i].r_min_index = j;
                }
                nb_info_[_i].min_mass = std::min(nb_info_[_i].min_mass, ptcl_[j].mass);
            }
        }
//            Jlist_disp_[i] = disp;
        Jlist_n_[_i] = n_pert;
#ifdef HARD_DEBUG
        assert(n_pert<=n_nb_off_);
#endif
    }

    //! update rsearch of ptcl
    /*!
      @param[in] _i_start: start index to calculate
      @param[in] _dt_tree: tree time step
     */
    void updateRSearch(const PS::S32 _i_start, const PS::F64 _dt_tree) {
        for(PS::S32 i=_i_start; i<ptcl_.size(); i++) {
            ptcl_[i].calcRSearch(_dt_tree);
        }
    }

    PS::S32* getPertList(const PS::S32 i) {
        return &Jlist_[Jlist_disp_[i]];
    }

    PS::S32 getPertN(const PS::S32 i) const {
        return Jlist_n_[i];
    }

    PS::S32 getPertListSize() const {
        return Jlist_.size();
    }

    PtclH4* getPtcl() const {
        return ptcl_.getPointer();
    }
    
    PtclForce* getForce() const {
        return force_.getPointer();
    }

    PS::S32 getPtclN() const {
        return ptcl_.size();
    }

    Tphard** getPtclAdr() const {
        return ptcl_ptr_.getPointer();
    }

    const NeighborInfo& getNbInfo(const PS::S32 i) const {
        return nb_info_[i];
    }

#ifdef HARD_DEBUG_PRINT
    void writePtcl(FILE* _fout, const PS::S32 _i_start) const{
        for (PS::S32 i=_i_start; i<ptcl_.size(); i++) {
            ptcl_[i].ParticleBase::writeAscii(_fout);
        }
    }
#endif

    //! Set integrator parameters
    /*!
      @param[in] _eta_s:  time step coefficient (square)
      @param[in] _r_in:   changeover function inner boundary
      @param[in] _r_out:  changeover function outer boundary
      @param[in] _eps_sq: softening 
      @param[in] _n_nb_off: neighbor list offset (maximum number)
     */
    void setParams(const PS::F64 _eta_s, 
                   const PS::F64 _r_in, const PS::F64 _r_out,
                   const PS::F64 _eps_sq,
                   const PS::S32 _n_nb_off) {
        eta_s_    = _eta_s;         
        r_in_     = _r_in;          
        r_out_    = _r_out;
        eps_sq_   = _eps_sq;
        r_oi_inv_ = 1.0/(_r_out-_r_in);
        r_A_      = (_r_out-_r_in)/(_r_out+_r_in);
        n_nb_off_ = _n_nb_off;
        n_act_    = 0;
    }

    //! calculate a0_offset_sq
    /*! calculate a0_offset_sq for timestep determination
     */
    void calcA0offset() {
        PS::F64 mass_min = PS::LARGE_FLOAT;
        for(PS::S32 i=0; i<ptcl_.size(); i++){
            if(mass_min > ptcl_[i].mass)  mass_min = ptcl_[i].mass;
            //if(rout_min > ptcl_[i].r_out) rout_min = ptcl_[i].r_out;
        }
        a0_offset_sq_ = 0.1 * mass_min / (r_out_ * r_out_);
    }

    //! Initial Hermite ptcl force and step 
    /*! Initial f, fdot and step
      @param[in] _list: ptcl index to initialize
      @param[in] _n_list: number of ptcl to initialize
      @param[in] _time_sys: current set time
      @param[in] _dt_max: maximum time step
      @param[in] _dt_min: minimum time step for check
      @param[in,out] _Aint: ARC integrator class (resolve and shift is used)
     */
    template <class ARCint>
    bool initial(const PS::S32* _list,
                 const PS::S32 _n_list,
                 const PS::F64 _time_sys,
                 const PS::F64 _dt_max,
                 const PS::F64 _dt_min,
                 ARCint* _Aint,
                 const bool _pred_flag = true) {

        // if no particle need initial, quit
        if(_n_list==0) return false;
#ifdef HARD_DEBUG
        assert(_n_list>0);
#endif

        const PS::S32* ptcl_list=_list;
        if(ptcl_list==NULL) ptcl_list = adr_dt_sorted_.getPointer();

        for(PS::S32 i=0; i<_n_list; i++){
            PS::S32 iadr = ptcl_list[i];
            pred_[iadr].r_search = ptcl_[iadr].r_search;
            ptcl_[iadr].time = ptcl_[iadr].dt = _time_sys;
            ptcl_[iadr].acc0 = ptcl_[iadr].acc1 = 0.0;
        }

        if(_pred_flag) PredictAll(pred_.getPointer(), ptcl_.getPointer(), ptcl_.size(), _time_sys);

        if(_Aint!=NULL) {
            _Aint->updateCM(ptcl_.getPointer());
            _Aint->resolve();
        }

        // force::acc0,acc1, neighbor list updated
//        if(_calc_full_flag) 
        CalcActAcc0Acc1(force_.getPointer(), 
                        nb_info_.getPointer(),
                        Jlist_.getPointer(), Jlist_n_.getPointer(), 
                        Jlist_disp_.getPointer(), 
                        ptcl_.getPointer(), ptcl_.size(), 
                        ptcl_list, _n_list,
                        r_in_, r_out_, r_oi_inv_, r_A_, eps_sq_, 
                        _Aint);
//        // only neighbor force calculation
//        else CalcAcc0Acc1ActNb(force_.getPointer(), 
//                               ptcl_.getPointer(), 
//                               ptcl_list, _n_list,
//                               Jlist_.getPointer(), Jlist_disp_.getPointer(), Jlist_n_.getPointer(), 
//                               r_in_, r_out_, r_oi_inv_, r_A_, eps_sq_, _Aint, _n_group);
    
        // store predicted force
        for(PS::S32 i=0; i<_n_list; i++){
            PS::S32 iadr = ptcl_list[i];
            ptcl_[iadr].acc0 = force_[iadr].acc0;
            ptcl_[iadr].acc1 = force_[iadr].acc1;
        }

        if(_Aint!=NULL) _Aint->shift2CM();

        bool fail_flag=CalcBlockDt2ndAct(ptcl_.getPointer(), force_.getPointer(), adr_dt_sorted_.getPointer(), _n_list, 0.01*eta_s_, _dt_max, _dt_min, a0_offset_sq_);

        for(PS::S32 i=0; i<_n_list; i++){
            PS::S32 iadr = ptcl_list[i];
            time_next_[iadr] = ptcl_[iadr].time + ptcl_[iadr].dt;
        }

        return fail_flag;
    }
    
    template<class Energy>
    void CalcEnergy(Energy & eng) {
        CalcEnergy(ptcl_.getPointer(), ptcl_.size(), eng, r_in_, r_out_, eps_sq_);
    }

    PS::F64 getNextTime() {
        return time_next_[adr_dt_sorted_[0]];
    }

    //PS::F64* getNextTimeList() const {
    //    return time_next_.getPointer();
    //}
    
    PS::F64 getOneTime(const std::size_t i) const {
        return ptcl_[i].time;
    }

    PS::F64 getOneDt(const std::size_t i) const {
        return ptcl_[i].dt;
    }

    //! Integration active particles
    /*! Integrated to time_sys
      @param[in] _time_sys: integration target time
      @param[in] _dt_max: maximum time step
      @param[in] _dt_min: minimum time step for check
      @param[in,out] _Aint: ARC integrator class (resolve and shift is used)
      \return fail_flag: If step size < dt_min, return true
     */
    template <class ARCint>
    bool integrateOneStepAct(const PS::F64 _time_sys,
                             const PS::F64 _dt_max,
                             const PS::F64 _dt_min,
                             ARCint* _Aint) {
        PS::S32 n_group = 0;
        // pred::mass,pos,vel updated
        PredictAll(pred_.getPointer(), ptcl_.getPointer(), ptcl_.size(), _time_sys);
        if(_Aint!=NULL) {
            n_group = _Aint->getNGroups();
            _Aint->updateCM(pred_.getPointer());
            _Aint->resolve();
        }
        // force::acc0,acc1, neighbor list updated
//        if(_calc_full_flag) 
        CalcActAcc0Acc1(force_.getPointer(), 
                        nb_info_.getPointer(),
                        Jlist_.getPointer(), Jlist_n_.getPointer(), 
                        Jlist_disp_.getPointer(), 
                        pred_.getPointer(), ptcl_.size(), 
                        adr_dt_sorted_.getPointer(), n_act_, 
                        r_in_, r_out_, r_oi_inv_, r_A_, eps_sq_, 
                        _Aint);
//        // only neighbor force calculation
//        else CalcAcc0Acc1ActNb(force_.getPointer(), 
//                               pred_.getPointer(), 
//                               adr_dt_sorted_.getPointer(), n_act_, 
//                               Jlist_.getPointer(), Jlist_disp_.getPointer(), Jlist_n_.getPointer(), 
//                               r_in_, r_out_, r_oi_inv_, r_A_, eps_sq_, _Aint, _n_group);

        // ptcl_org::pos,vel; pred::time,dt,acc0,acc1,acc2,acc3 updated
        bool fail_flag=CorrectAndCalcDt4thAct(ptcl_.getPointer(), force_.getPointer(), adr_dt_sorted_.getPointer(), n_act_, _dt_max, _dt_min, a0_offset_sq_, eta_s_);

        for(PS::S32 i=0; i<n_act_; i++){
            PS::S32 adr = adr_dt_sorted_[i];
            time_next_[adr] = ptcl_[adr].time + ptcl_[adr].dt;
            // update new perturber list
            if(adr<n_group) _Aint->updatePertOneGroup(adr, ptcl_.getPointer(), force_.getPointer(), getPertList(adr), getPertN(adr));
        }

        if(_Aint!=NULL) {
            // shift member to c.m. frame
            _Aint->shift2CM();
            // go back to original time
            _Aint->updateCM(ptcl_.getPointer());
            //Aint.updateCM(Hint.getPtcl(), group_act_list.getPointer(), group_act_n);
        }
        
        return fail_flag;
    }

    //! Integration one particle to targat time (ingore dt)
    /*! Integrated to time_sys but no prediction and no step calculation
      @param[in] _ptcl_list: particle index list to integrate
      @param[in] _n_list: number of particles
      @param[in] _time_sys: integration target time
      @param[in] _dt_max: maximum time step
      @param[in] _dt_min: minimum time step for check
      @param[in,out] _Aint: ARC integrator class (resolve and shift is used)
      \return fail_flag: If step size < dt_min, return true
     */
    template <class ARCint>
    bool integrateOneListNoPred(const PS::S32* _ptcl_list,
                                const PS::S32 _n_list,
                                const PS::F64 _time_sys,
                                const PS::F64 _dt_max,
                                const PS::F64 _dt_min,
                                ARCint* _Aint) {

        if(_n_list==0) return false;
#ifdef HARD_DEBUG
        assert(_n_list>0);
#endif

        PS::S32 n_group = 0;
        if(_Aint!=NULL) n_group = _Aint->getNGroups();

        // adjust dt
        PS::S32 int_list[_n_list];
        PS::S32 group_int_list[n_group+1];
        PS::S32 n_group_int=0;
        PS::S32 n_int=0;
        for (PS::S32 i=0; i<_n_list; i++) {
            PS::S32 iadr = _ptcl_list[i];
            if(ptcl_[iadr].time==_time_sys) continue;
            ptcl_[iadr].dt = _time_sys - ptcl_[iadr].time;
            int_list[n_int++] = iadr;
            if(iadr<n_group) group_int_list[n_group_int++] = iadr;
        }

        if (n_int==0) return false;

        if(_Aint!=NULL) _Aint->resolve();

        // force::acc0,acc1, neighbor list updated
//        if(_calc_full_flag) 
        CalcActAcc0Acc1(force_.getPointer(), 
                        nb_info_.getPointer(),
                        Jlist_.getPointer(), Jlist_n_.getPointer(), 
                        Jlist_disp_.getPointer(), 
                        pred_.getPointer(), ptcl_.size(), 
                        int_list, n_int,
                        r_in_, r_out_, r_oi_inv_, r_A_, eps_sq_, 
                        _Aint);
//        // only neighbor force calculation
//        else CalcAcc0Acc1ActNb(force_.getPointer(), 
//                               pred_.getPointer(), 
//                               int_list, n_int,
//                               Jlist_.getPointer(), Jlist_disp_.getPointer(), Jlist_n_.getPointer(), 
//                               r_in_, r_out_, r_oi_inv_, r_A_, eps_sq_, _Aint, _n_group);


        bool fail_flag=CorrectAndCalcDt4thAct(ptcl_.getPointer(), 
                                              force_.getPointer(), 
                                              int_list, n_int,
                                              _dt_max, _dt_min, a0_offset_sq_, eta_s_);


        // shift member to c.m. frame
        if(_Aint!=NULL) {
            _Aint->shift2CM();
            // update c.m. of _i
            _Aint->updateCM(ptcl_.getPointer(), group_int_list, n_group_int);
        }

        for(PS::S32 i=0; i<n_int; i++){
            PS::S32 iadr=int_list[i];
            // update time table
            time_next_[iadr] = ptcl_[iadr].time + ptcl_[iadr].dt;
        }

        // update new perturber list
        for(PS::S32 i=0; i<n_group_int; i++) {
            PS::S32 iadr=group_int_list[i];
            _Aint->updatePertOneGroup(iadr, ptcl_.getPointer(), force_.getPointer(), getPertList(iadr), getPertN(iadr));
        }

        return fail_flag;
    }

    void SortAndSelectIp(PS::S32 group_act_list[],
                         PS::S32 &group_act_n,
                         const PS::S32 n_groups) {
        SortAndSelectIp(adr_dt_sorted_.getPointer(), time_next_.getPointer(), n_act_, time_next_.size(), group_act_list, group_act_n, n_groups);
        
    }

    void SortAndSelectIp() {
        SortAndSelectIp(adr_dt_sorted_.getPointer(), time_next_.getPointer(), n_act_, adr_dt_sorted_.size());
    }
    
    PS::S32 getNact() const{
        return n_act_;
    }

    const PS::S32* getActList() const{
        return adr_dt_sorted_.getPointer();
    }

#ifdef HARD_DEBUG
    template <class ARCint>
    void checkAdrList(const ARCint& _Aint) {
        const PS::S32 n_ptcl = ptcl_.size();
        const PS::S32 n_adr = adr_dt_sorted_.size();
        PS::S32 ncheck[n_ptcl];
        for (PS::S32 i=0; i<n_ptcl; i++) ncheck[i]=0;
        for (PS::S32 i=0; i<n_adr; i++) {
            PS::S32 k =adr_dt_sorted_[i];
            assert(time_next_[k] == ptcl_[k].time + ptcl_[k].dt);
            PS::F64 tr=PS::S32(ptcl_[k].time/ptcl_[k].dt);
            assert(tr*ptcl_[k].dt==ptcl_[k].time);
            ncheck[adr_dt_sorted_[i]]++;
        }
        const PS::S32 n_group = _Aint.getNGroups();
        for (PS::S32 i=0; i<n_group; i++) {
            if(_Aint.getMask(i)) assert(ncheck[i]==0);
            else assert(ncheck[i]==1);
        }
        for (PS::S32 i=n_group; i<n_ptcl; i++) {
            assert(ncheck[i]==1);
        }
        for (PS::S32 i=0; i<n_adr-1; i++) {
            assert(ptcl_[adr_dt_sorted_[i]].dt<=ptcl_[adr_dt_sorted_[i+1]].dt);
            assert(time_next_[adr_dt_sorted_[i]]<=time_next_[adr_dt_sorted_[i+1]]);
        }
    }

    void printStepHist(){
        std::map<PS::F64, PS::S32> stephist;
        for(int i=0; i<adr_dt_sorted_.size(); i++) {
            PS::S32 k = adr_dt_sorted_[i];
            std::map<PS::F64, PS::S32>::iterator p = stephist.find(ptcl_[k].dt);
            if (p==stephist.end()) stephist[ptcl_[k].dt]=1;
            else stephist[ptcl_[k].dt]++;
        }
        std::cerr<<"Step hist:\n";
        for(auto i=stephist.begin(); i!=stephist.end(); i++) {
            std::cerr<<std::setw(24)<<i->first;
        }
        std::cerr<<std::endl;
        for(auto i=stephist.begin(); i!=stephist.end(); i++) {
            std::cerr<<std::setw(24)<<i->second;
        }
        std::cerr<<std::endl;
    }
#endif

};

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

#ifdef TIDAL_TENSOR
class TidalTensor{
private:
    PS::F64 T2[6];  // 1st (6)
    PS::F64 T3[10]; // 2nd Tensor (10)
    bool use_flag;
public:

    // only initialize flag
    TidalTensor(): use_flag(false) {}

    void dump(FILE *fp){
        fwrite(this, sizeof(TidalTensor),1,fp);
    }

    void read(FILE *fp){
        size_t rcount = fread(this, sizeof(TidalTensor),1,fp);
        if (rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            abort();
        }
    }

    void reset(){
        use_flag=false;
        for(PS::S32 i=0; i<6; i++) T2[i] = 0;
        for(PS::S32 i=0; i<10; i++) T3[i] = 0;
    }

    //! tidal tensor fitting function,
    /* @param[in] _ptcl_tt: tidal tensor measure particles
       @param[in] _bin: binary information, get the scaling factor of distance
       @param[in] _n_split: artifical particle splitting number
     */
    template<class Tptcl>
    void fit(Tptcl* _ptcl_tt, const Binary& _bin, const PS::S32 _n_split) {
        use_flag=true;
        PS::F64vec fi[8];
#ifdef HARD_DEBUG
        assert(_ptcl_tt[12].mass_bk==0);
        assert(_n_split>4);
#endif
        // get acceleration
        for (PS::S32 i=0; i<8; i++) fi[i] = _ptcl_tt[i].acc;
        // get cofficients
        for (PS::S32 i=0; i<8; i++) {
            T2[0] =  0.250000000000000*fi[0][0] + -0.250000000000000*fi[2][0] +  0.250000000000000*fi[4][0] + -0.250000000000000*fi[6][0];
            T2[1] =  0.125000000000000*fi[0][1] +  0.125000000000000*fi[1][0] + -0.125000000000000*fi[2][1] + -0.125000000000000*fi[3][0] 
                +    0.125000000000000*fi[4][1] +  0.125000000000000*fi[5][0] + -0.125000000000000*fi[6][1] + -0.125000000000000*fi[7][0];
            T2[2] = -0.083333333333333*fi[0][0] +  0.083333333333333*fi[0][2] + -0.083333333333333*fi[1][0] + -0.083333333333333*fi[2][0]
                +   -0.083333333333333*fi[2][2] + -0.083333333333333*fi[3][0] +  0.083333333333333*fi[4][0] +  0.083333333333333*fi[4][2]
                +    0.083333333333333*fi[5][0] +  0.083333333333333*fi[6][0] + -0.083333333333333*fi[6][2] +  0.083333333333333*fi[7][0];
            T2[3] =  0.250000000000000*fi[1][1] + -0.250000000000000*fi[3][1] +  0.250000000000000*fi[5][1] + -0.250000000000000*fi[7][1];
            T2[4] = -0.083333333333333*fi[0][1] + -0.083333333333333*fi[1][1] +  0.083333333333334*fi[1][2] + -0.083333333333333*fi[2][1]
                +   -0.083333333333333*fi[3][1] + -0.083333333333333*fi[3][2] +  0.083333333333333*fi[4][1] +  0.083333333333333*fi[5][1]
                +    0.083333333333333*fi[5][2] +  0.083333333333333*fi[6][1] +  0.083333333333333*fi[7][1] + -0.083333333333334*fi[7][2];
            T2[5] = -0.124999999999999*fi[0][2] + -0.125000000000000*fi[1][2] + -0.125000000000001*fi[2][2] + -0.125000000000000*fi[3][2]
                +    0.125000000000001*fi[4][2] +  0.125000000000000*fi[5][2] +  0.124999999999999*fi[6][2] +  0.125000000000000*fi[7][2];

            T3[0] =  0.250000000000001*fi[0][0] +  0.124999999999999*fi[0][2] +  0.250000000000001*fi[2][0] + -0.124999999999999*fi[2][2]
                +    0.249999999999999*fi[4][0] + -0.125000000000001*fi[4][2] +  0.249999999999999*fi[6][0] +  0.125000000000001*fi[6][2];
            T3[1] =  0.250000000000000*fi[0][1] +  0.125000000000000*fi[1][2] +  0.250000000000000*fi[2][1] + -0.125000000000000*fi[3][2]
                +    0.250000000000000*fi[4][1] + -0.125000000000000*fi[5][2] +  0.250000000000000*fi[6][1] +  0.125000000000000*fi[7][2];
            T3[2] = -0.112500000000000*fi[0][0] +  0.025000000000000*fi[0][2] + -0.012500000000000*fi[1][1] + -0.025000000000000*fi[1][2]
                +    0.112500000000000*fi[2][0] +  0.025000000000000*fi[2][2] +  0.012500000000000*fi[3][1] + -0.025000000000000*fi[3][2]
                +    0.112500000000000*fi[4][0] +  0.025000000000000*fi[4][2] +  0.012500000000000*fi[5][1] + -0.025000000000000*fi[5][2]
                +   -0.112500000000000*fi[6][0] +  0.025000000000000*fi[6][2] + -0.012500000000000*fi[7][1] + -0.025000000000000*fi[7][2];
            T3[3] =  0.125000000000000*fi[0][2] +  0.250000000000000*fi[1][0] + -0.125000000000000*fi[2][2] +  0.250000000000000*fi[3][0]
                +   -0.124999999999999*fi[4][2] +  0.250000000000000*fi[5][0] +  0.125000000000000*fi[6][2] +  0.250000000000000*fi[7][0];
            T3[4] = -0.062500000000000*fi[0][1] + -0.062500000000000*fi[1][0] +  0.062500000000000*fi[2][1] +  0.062500000000000*fi[3][0]
                +    0.062500000000000*fi[4][1] +  0.062500000000000*fi[5][0] + -0.062500000000000*fi[6][1] + -0.062500000000000*fi[7][0];
            T3[5] = -0.125000000000000*fi[0][2] +  0.125000000000000*fi[2][2] +  0.125000000000000*fi[4][2] + -0.125000000000000*fi[6][2];
            T3[6] =  0.250000000000000*fi[1][1] +  0.125000000000000*fi[1][2] +  0.250000000000000*fi[3][1] + -0.125000000000000*fi[3][2]
                +    0.250000000000000*fi[5][1] + -0.125000000000000*fi[5][2] +  0.250000000000000*fi[7][1] +  0.125000000000000*fi[7][2];
            T3[7] = -0.012500000000000*fi[0][0] + -0.025000000000000*fi[0][2] + -0.112500000000000*fi[1][1] +  0.025000000000000*fi[1][2]
                +    0.012500000000000*fi[2][0] + -0.025000000000000*fi[2][2] +  0.112500000000000*fi[3][1] +  0.025000000000000*fi[3][2]
                +    0.012500000000000*fi[4][0] + -0.025000000000000*fi[4][2] +  0.112500000000000*fi[5][1] +  0.025000000000000*fi[5][2]
                +   -0.012500000000000*fi[6][0] + -0.025000000000000*fi[6][2] + -0.112500000000000*fi[7][1] +  0.025000000000000*fi[7][2];
            T3[8] = -0.125000000000000*fi[1][2] +  0.125000000000000*fi[3][2] +  0.125000000000000*fi[5][2] + -0.125000000000000*fi[7][2];
            T3[9] =  0.062500000000000*fi[0][0] +  0.125000000000000*fi[0][2] +  0.062500000000000*fi[1][1] +  0.125000000000000*fi[1][2]
                +   -0.062500000000000*fi[2][0] +  0.125000000000000*fi[2][2] + -0.062500000000000*fi[3][1] +  0.125000000000000*fi[3][2]
                +   -0.062500000000000*fi[4][0] +  0.125000000000000*fi[4][2] + -0.062500000000000*fi[5][1] +  0.125000000000000*fi[5][2]
                +    0.062500000000000*fi[6][0] +  0.125000000000000*fi[6][2] +  0.062500000000000*fi[7][1] +  0.125000000000000*fi[7][2];
        }
        // Rescale
        PS::F64 T2S = 1.0/(_bin.semi*(1+_bin.ecc)*0.35);
        PS::F64 T3S = T2S*T2S;
        for (PS::S32 i=0; i<6;  i++) T2[i] *= T2S;
        for (PS::S32 i=0; i<10; i++) T3[i] *= T3S;
    }

    void eval(PS::F64* acc, const PS::F64vec &pos) const {
        if(!use_flag) return;
        /*
          T2:
          [[0 1 2]
          [1 3 4]
          [2 4 5]]

          T3:
          [[[6 7 8]
          [7 9 10]
          [8 10 11]]

          [[7 9 10]
          [9 12 13]
          [10 13 14]]

          [[8 10 11]
          [10 13 14]
          [11 14 15]]]

         */
        PS::F64 x = pos.x;
        PS::F64 y = pos.y;
        PS::F64 z = pos.z;
        PS::F64 x2 = x*x;
        PS::F64 xy = x*y;
        PS::F64 xz = x*z;
        PS::F64 y2 = y*y;
        PS::F64 yz = y*z;
        PS::F64 z2 = z*z;

        PS::F64 acc0=acc[0];
        PS::F64 acc1=acc[1];
        PS::F64 acc2=acc[2];

        acc0 +=  T2[0]*x + T2[1]*y + T2[2]*z 
            +      T3[0]*x2 + 2*T3[1]*xy + 2*T3[2]*xz + T3[3]*y2 + 2*T3[4]*yz + T3[5]*z2;
        acc1 +=  T2[1]*x + T2[3]*y + T2[4]*z
            +      T3[1]*x2 + 2*T3[3]*xy + 2*T3[4]*xz + T3[6]*y2 + 2*T3[7]*yz + T3[8]*z2;
        acc2 +=  T2[2]*x + T2[4]*y + T2[5]*z
            +      T3[2]*x2 + 2*T3[4]*xy + 2*T3[5]*xz + T3[7]*y2 + 2*T3[8]*yz + T3[9]*z2;

        acc[0] = acc0;
        acc[1] = acc1;
        acc[2] = acc2;
    }
};

////! substract c.m. force from component used for tidal tensor. 
///* @param[in]     _ptcl_cm: c.m. particle
//   @param[in,out] _ptcl_tt: tidal tensor artifical particle (8)
//*/
//template<class Tptcl>
//void subtractFcm(Tptcl& _ptcl_cm, Tptcl* _ptcl_tt) {
//    PS::F64vec& acc_cm = _ptcl_cm.acc;
// 
//    // remove center force
//    for (PS::S32 i=0; i<8; i++)  _ptcl_tt[i].acc -= acc_cm;
//}

#else
class keplerSplineFit{
private:
    const gsl_interp_type *t_;
    PS::F64 peri_;
    PS::F64 tperi_;
    gsl_interp_accel *acc_[3];
    gsl_spline *spline_[3];
public:

    keplerSplineFit(): t_(gsl_interp_cspline_periodic), peri_(0.0) {
        acc_[0]=acc_[1]=acc_[2]=NULL;
        spline_[0]=spline_[1]=spline_[2]=NULL;
    }
    
    //! orbital force fitting function,
    /* @param[in] _ptcl_orb: orbital particles
       @param[in] _bin: binary information, get the scaling factor of distance
       @param[in] _n_split: artifical particle splitting number
     */
    template<class Tptcl> 
    void fit(Tptcl* _ptcl_orb, const Binary& _bin, const PS::S32 _n_split) {
#ifdef HARD_DEBUG
        assert(_n_split>=4);
#endif
        const PS::F64 twopi = PI*2.0;
        peri_  = _bin.peri;
        tperi_ = _bin.tperi - twopi;  // make sure tperi_ is negative
        const PS::U32 np = _n_split+1;
        PS::F64 x[np],y[3][np];
        for(PS::U32 i=0; i<np; i++) {
            x[i] = PS::F64(i)/_n_split*twopi;
            x[i] = _bin.peri/twopi * (x[i] - _bin.ecc*std::sin(x[i]));
        }
        for(PS::S32 i=0; i<_n_split; i++) {            
            for(PS::S32 j=0; j<3; j++)
                y[j][i] = _ptcl_orb[2*i+1].acc[j] - _ptcl_orb[2*i].acc[j];  // acc1 - acc0
        }
        for(PS::S32 i=0; i<3; i++) {
            y[i][_n_split] = y[i][0];

            if(acc_[i]!=NULL) gsl_interp_accel_free(acc_[i]);
            acc_[i] = gsl_interp_accel_alloc();
            if(spline_[i]!=NULL) gsl_spline_free(spline_[i]);
            spline_[i] =  gsl_spline_alloc(t_, np);
            
            gsl_spline_init(spline_[i], x, y[i], np);
        }
    }

    void eval(PS::F64* acc, const PS::F64 time) const {
        PS::F64 dt = time - tperi_;
        dt = dt - (PS::S32)(dt/peri_)*peri_;
        if(spline_[0]!=NULL) {
            for(PS::S32 i=0; i<3; i++)
                acc[i] -= gsl_spline_eval(spline_[i],dt,acc_[i]);
        }
    }

    void reset(){
        for(PS::S32 i=0;i<3;i++) {
            if(acc_[i]!=NULL) gsl_interp_accel_free(acc_[i]);
            if(spline_[i]!=NULL) gsl_spline_free(spline_[i]);
            acc_[i]=NULL;
            spline_[i]=NULL;
        }
    }

    ~keplerSplineFit() {
        reset();
    }
};
#endif

//! ARC integrator extra data
/*!
 */
class ARC_int_pars{
public:
    PS::F64 rout;      ///> r out
    PS::F64 rin;       ///> r in
    PS::F64 r_oi_inv;  ///> 1.0/(rout-rin)
    PS::F64 r_A;       ///> (rout-rin)/(rout+rin)
    PS::F64 pot_off;   ///> (1 + r_A)/rout
    PS::F64 eps2;      ///> eps*eps
    
    ARC_int_pars() {}
    ARC_int_pars(const ARC_int_pars& in_) {
        rout     = in_.rout;
        rin      = in_.rin;
        r_oi_inv = in_.r_oi_inv;
        r_A      = in_.r_A;
        pot_off  = in_.pot_off;
        eps2     = in_.eps2;
    }

    void dump(FILE *fp) {
        fwrite(this, sizeof(ARC_int_pars),1,fp);
    }

    void read(FILE *fp) {
        size_t rcount = fread(this, sizeof(ARC_int_pars),1,fp);
        if (rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            abort();
        }
    }
};

#ifdef TIDAL_TENSOR
class ARC_pert_pars: public ARC_int_pars, public TidalTensor{
public:
    ARC_pert_pars() {}
    ARC_pert_pars(const ARC_int_pars& in_): ARC_int_pars(in_) {}

    void dump(FILE *fp) {
        fwrite(this, sizeof(ARC_pert_pars),1,fp);
    }

    void read(FILE *fp) {
        size_t rcount = fread(this, sizeof(ARC_pert_pars),1,fp);
        if (rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            abort();
        }
    }
};
#else
class ARC_pert_pars: public ARC_int_pars, public keplerSplineFit{
public:
    ARC_pert_pars() {}
    ARC_pert_pars(const ARC_int_pars& in_): ARC_int_pars(in_) {}

    void dump(FILE *fp) {
        fwrite(this, sizeof(ARC_pert_pars),1,fp);
    }

    void read(FILE *fp) {
        size_t rcount = fread(this, sizeof(ARC_pert_pars),1,fp);
        if (rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            abort();
        }
    }

};
#endif

template<class TpARC, class Tpert, class Tpforce>
class ARCIntegrator{
private:
    typedef ARC::chain<TpARC> ARChain;
    typedef ARC::chainpars ARControl;
    PS::ReallocatableArray<ARChain> clist_;
    PS::ReallocatableArray<ARC_pert_pars> par_list_;
    PS::ReallocatableArray<Tpert*> pert_;
    PS::ReallocatableArray<Tpforce*> pforce_;
    PS::ReallocatableArray<PS::S32> pert_n_;
    PS::ReallocatableArray<PS::S32> pert_disp_;
    PS::ReallocatableArray<bool> group_mask_map_; 
    PS::ReallocatableArray<PS::S32> group_mask_list_;
    PS::S32 n_pert_off_;

    ARControl *ARC_control_;
    ARC_int_pars *Int_pars_;

public:
    PS::ReallocatableArray<Binary> bininfo;
    //PS::ReallocatableArray<PS::F64> dt;
#ifdef ARC_SYM
    PS::S32 step_count_limit;
#endif

    ARCIntegrator() {};
    ARCIntegrator(ARControl &contr, ARC_int_pars &par): n_pert_off_(0), ARC_control_(&contr), Int_pars_(&par) {}

    void reserveARMem(const PS::S32 _n) {
#ifdef HARD_DEBUG
        assert(_n<ARRAY_ALLOW_LIMIT);
#endif        
        clist_.reserve(_n);
        //clist_.resizeNoInitialize(n);
        par_list_.reserve(_n);
        //par_list_.resizeNoInitialize(n);
        pert_n_.reserve(_n);
        //pert_n_.resizeNoInitialize(n);
        pert_disp_.reserve(_n);
        //pert_disp_.resizeNoInitialize(n);
        bininfo.reserve(_n);
        //dt.reserve(n);
        group_mask_map_.reserve(_n);
        group_mask_list_.reserve(_n);
    }

    //! Reserve memory for perturber list
    /*! reserve memory for perturber and force address list
      @param[in] _n_group: maximum number of groups (can be enlarged)
      @param[in] _n_pert_off: maximum perturber maxinum number (offset for perturber array, cannot be changed)
     */
    void reservePertMem(const PS::S32 _n_group, const PS::S32 _n_pert_off) {
        pert_.reserve(_n_group*_n_pert_off);
        pforce_.reserve(_n_group*_n_pert_off);
        n_pert_off_ = _n_pert_off;
    }
    //void initialize(PS::S32 group_list[];
    //                ReallocatableArray<TpARC> groups[],
    //                const PS::S32 n_groups,
    //                TpARC* ptcl,
    //                PS::S32 adr_cm[],
    //                ReallocatableArray<PS::S32> pertlist[],
    //                ARControl &control,
    //                ARC_int_pars &Int_pars) {
    //    for(int i=0; i<n_groups; i++) {
    //        PS::S32 icm = adr_cm[i];
    //        PS::S32 ig = group_list[i];
    //        PS::S32 ni = groups[ig].size();
    //        clist_.push_back(ARChain(ni,control));
    //        clist_[i].addP(ni,groups[ig].getPointer());
    //        clist_[i].link_int_par(Int_pars);
    //        for(int j=0; j<pertlist[icm].size(); j++) clist_[i].addPext(ptcl[pertlist[icm][j]]);
    //        clist_[i].init(0);
    //    }
    //    ARC_control_ = &control;
    //    Int_pars_ = &Int_pars;
    //}

    //! Add group of particles to ARC class
    /*! Add one group in ARC, use the suppressed group first, if no, add new group
      @param[in] _ptcl: particle data array
      @param[in] _ptcl_list: particle index array for _ptcl (if NULL, read continuelly from 1 to _n_ptcl
      @param[in] _n_ptcl: number of particles
      @param[in] _ptcl_soft_pert: soft perturbation artifical particles (if NULL, keep pert_par same as before)
      @param[in] _n_split: split number for artifical particles
      @param[in] _ptcl_pert: perturber particle array, notice the first _n_group are c.m. which has consistent order of ARC groups
      @param[in] _force_pert: perturber force array
      @param[in] _ptcl_pert_list: perturber particle index in _ptcl_pert
      @param[in] _n_pert: number of perturbers
      \return adding group index
    */
    template <class Tptcl, class Tpsoft>
    PS::S32 addOneGroup(Tptcl* _ptcl,
                        const PS::S32* _ptcl_list,
                        const PS::S32 _n_ptcl,
                        const Tpsoft* _ptcl_soft_pert,
                        const PS::S32 _n_split,
                        Tpert* _ptcl_pert = NULL,
                        Tpforce* _force_pert = NULL,
                        const PS::S32* _ptcl_pert_list = NULL,
                        const PS::S32 _n_pert = 0) {
        // set current group offset
        const PS::S32 ngroup = clist_.size();
        PS::S32 igroup;

        // check suppressed group
        if(group_mask_list_.size()>0) {
            igroup = group_mask_list_.back();
#ifdef HARD_DEBUG
            assert(group_mask_map_[igroup]);
#endif
            group_mask_map_[igroup] = false;
            group_mask_list_.decreaseSize(1);
        }
        else {
            // add new
            igroup = ngroup;
            pert_n_.push_back(0);
            pert_disp_.push_back(pert_.size());
            group_mask_map_.push_back(false);

            pert_.increaseSize(n_pert_off_);
            pforce_.increaseSize(n_pert_off_);
        
            clist_.increaseSize(1);
            bininfo.increaseSize(1);

#ifdef HARD_DEBUG
            assert(pert_disp_.size()==clist_.size());
            assert(pert_n_.size()==clist_.size());
#endif
            par_list_.push_back(ARC_pert_pars(*Int_pars_));
        }

        // Soft perturbation
        if(_ptcl_soft_pert)
            par_list_[igroup].fit(_ptcl_soft_pert, bininfo[igroup], _n_split);

        // allocate memory
        clist_[igroup].allocate(_n_ptcl);
        
        // set current pert_disp
        const PS::S32 i_pert_off = pert_disp_[igroup];


        // Add members to ARC 
        if(_ptcl_list) { // use index from list 
            for(PS::S32 i=0; i<_n_ptcl; i++) {
                clist_[igroup].addP(_ptcl[_ptcl_list[i]]);
            }
        }
        else { // read one by one
            for(PS::S32 i=0; i<_n_ptcl; i++) {
                clist_[igroup].addP(_ptcl[i]);
            }
        }
        
        // c.m. position is in igroup, put the c.m. particle to perturber list thus it can be used for predicting the c.m. position and velocity
        if(_ptcl_pert!=NULL) {
            pert_  [i_pert_off] = &_ptcl_pert [igroup];   
            pforce_[i_pert_off] = &_force_pert[igroup];
            pert_n_[igroup]++;
        }

        // Add perturber
        for(PS::S32 i=0; i<_n_pert; i++) {
            const PS::S32  k = _ptcl_pert_list[i];
            pert_  [i+i_pert_off+1] = &_ptcl_pert[k];
            pforce_[i+i_pert_off+1] = &_force_pert[k];
            pert_n_[igroup]++;
        }
#ifdef HARD_DEBUG
        assert(_n_pert+1<=n_pert_off_);
        assert(par_list_.size()==clist_.size());
#endif

        // recored c.m. inforamtion
        if(_ptcl_pert!=NULL) {
            clist_[igroup].DataCopy(_ptcl_pert[igroup]);
            //clist_.back().pos  = _ptcl_pert[igroup].pos;
            //clist_.back().vel  = _ptcl_pert[igroup].vel;
            //clist_.back().mass = _ptcl_pert[igroup].mass;
#ifdef HARD_DEBUG
            if(igroup==ngroup) assert(clist_[igroup].mass>0.0);
            //assert(clist_.back().mass==_ptcl_pert[igroup].mass);
#endif
        }

        return igroup;
    }

    //! Clear one group
    /*! Clear one group, remove members, relink to c.m., member number to 0
     */
    void clearOneGroup(const PS::S32 _igroup) {
#ifdef HARD_DEBUG
        assert(_igroup<clist_.size());
#endif
        if(!group_mask_map_[_igroup]) {
            clist_[_igroup].clear();
            clist_[_igroup].slowdown.reset();
            bininfo[_igroup].tstep = -1.0;
            pert_n_[_igroup]=0;
            group_mask_list_.push_back(_igroup);
            group_mask_map_[_igroup]=true;
        }
    }

    //! Copy ARC_par
    /*! Copy ARC_pars (soft perturbation) from _i_source to _i_target
      @param[in] _i_target: target to copy
      @param[in] _i_source: source for copy
     */
    void copyParP2P(const PS::S32 _i_target, const PS::S32 _i_source) {
#ifdef HARD_DEBUG
        assert(_i_target>=0||_i_target<par_list_.size());
        assert(_i_source>=0||_i_source<par_list_.size());
#endif
        par_list_[_i_target] = par_list_[_i_source];
    }

    //! Update perturber list
    /*! Update perturber list for group i
      @param[in] _i_group: group index for update perturber
      @param[in] _ptcl_pert: perturber particle array
      @param[in] _force_pert: perturber force array
      @param[in] _ptcl_pert_list: new perturber particle index
      @param[in] _n_pert: number of perturbers
     */
    void updatePertOneGroup(const PS::S32 _igroup,
                            Tpert* _ptcl_pert,
                            Tpforce* _force_pert,
                            const PS::S32* _ptcl_pert_list,
                            const PS::S32 _n_pert) {
#ifdef HARD_DEBUG
        assert(!group_mask_map_[_igroup]);
#endif
        // Add one for c.m.
        pert_n_[_igroup] = _n_pert + 1;
        const PS::S32 i_pert_off = pert_disp_[_igroup];
        
        for (PS::S32 i=0; i<_n_pert; i++) {
            PS::S32 adr = _ptcl_pert_list[i];
            pert_  [i+i_pert_off+1] = &_ptcl_pert[adr];
            pforce_[i+i_pert_off+1] = &_force_pert[adr];
        }
    }

#ifdef HARD_DEBUG
    //! check perturber list
    void checkPert() {
        for (PS::S32 i=0; i<clist_.size(); i++) {
            if(!group_mask_map_[i]){
                const PS::S32 n_member = clist_[i].getN();
                const PS::S32 i_pert_off = pert_disp_[i];
                for (PS::S32 j=1; j<pert_n_[i]; j++) {
                    for (PS::S32 k=0; k<n_member; k++) {
                        assert(clist_[i].getP(k).id!=pert_[j+i_pert_off]->id);
                    }
                    assert(clist_[i].id!=pert_[j+i_pert_off]->id);
                }
                assert(clist_[i].id==pert_[i_pert_off]->id);
            }
        }
    }
#endif


    //! Set initial slowdown parameter for one unpert group 
    /*! 
      @param[in] _i_group: group index
      @param[in] _tend: ending physical time for integration
      @param[in] _sdfactor: slowdown criterion factor
      @param[in] _tp_factor: if minimum factor of integration time interval / (kappa * period).
    */
    void initialOneSlowDownUnPert(const PS::S32 _i_group, const PS::F64 _tend, const PS::F64 _sdfactor = 1.0e-8, const PS::F64 _tp_factor = 0.01) {
        // isolated case
        if (bininfo[_i_group].semi>0&&bininfo[_i_group].stable_factor>=0) {   
            PS::F64 finner = bininfo[_i_group].semi*(1.0+bininfo[_i_group].ecc);
            finner = (bininfo[_i_group].m1-bininfo[_i_group].m2)/(finner*finner);
            PS::F64 finnersq = finner*finner;
            TpARC p[2];
            OrbParam2PosVel(p[0].pos, p[1].pos, p[0].vel, p[1].vel, bininfo[_i_group].m1, bininfo[_i_group].m2, bininfo[_i_group].semi, bininfo[_i_group].ecc, bininfo[_i_group].inc, bininfo[_i_group].OMG, bininfo[_i_group].omg, PI);
            p[0].mass = bininfo[_i_group].m1;
            p[1].mass = bininfo[_i_group].m2;
#ifdef SOFT_PERT
#ifndef TIDAL_TENSOR
            p[0].status = 0;
            p[1].status = 1;
#endif
#endif
            //center_of_mass_correction(*(TpARC*)&clist_[_i_group], p, 2);
            PS::F64 acc[2][3];
            const PS::S32 ipert = pert_disp_[_i_group];
            //Newtonian_extA(acc, bininfo[i].tperi+bininfo[i].peri, p, 2, &pert_[ipert], &pforce_[ipert], pert_n_[i], &par_list_[i]);
            //if(pert_n_[i]>1) Newtonian_extA_pert(acc, 0.0, p, 2, &pert_[ipert], &pforce_[ipert], pert_n_[i], &par_list_[i]);
            Newtonian_extA_soft(acc, 0.0, p, 2, &pert_[ipert], &pforce_[ipert], pert_n_[_i_group], &par_list_[_i_group]);
            PS::F64 fpertsq = 0.0;
            for(PS::S32 k=0; k<3; k++) {
                PS::F64 dacc = acc[0][k]-acc[1][k];
                fpertsq += dacc*dacc;
            }
            clist_[_i_group].slowdown.setSlowDownPars(bininfo[_i_group].peri, _sdfactor);
            clist_[_i_group].slowdown.updatefratiosq(fpertsq/finnersq);
            clist_[_i_group].slowdown.updatekappa(_tend, 1.0, _tp_factor,-1);
        }
    }

    void initialOneSlowDown(const PS::S32 _i_group, const PS::F64 _tend, const PS::F64 _mpert, const PS::F64 _sdfactor = 1.0e-8, const PS::F64 _tp_factor = 1e-4) {
        if (bininfo[_i_group].semi>0&&bininfo[_i_group].stable_factor>=0) {
            clist_[_i_group].slowdown.setSlowDownPars(bininfo[_i_group].peri, _sdfactor);
            clist_[_i_group].slowdown.updatekappa(_tend, clist_[_i_group].mass/_mpert, _tp_factor,-1);
        }
    }

    //! Update slow down factor for one ARC
    /*! Update slowdown for one ARC
      @param[in] _igroup: index of ARC
      @param[in] _tnow: current time of c.m.
      @param[in] _dt: c.m. step size
      @param[in] _dt_limit: step limit
      @param[in] _mpert: nearest perturber mass
      @param[in] _md_factor: slowdown modification limit factor (negative suppress the limit)
     */
    void updateOneSlowDown(const size_t _igroup, const PS::F64 _tnow, const PS::F64 _dt, const PS::F64 _dt_limit, const PS::F64 _mpert, const PS::F64 _md_factor=1.2) {
#ifdef HARD_DEBUG
        assert(!group_mask_map_[_igroup]);
#endif
        //PS::F64 tp_factor = std::max(1e-4,_dt/_dt_limit);
        //PS::F64 tp_factor = std::max(1e-4,_dt/_dt_limit);
        PS::F64 tp_factor = _dt/_dt_limit;
        //std::cerr<<"i "<<_index<<" dt "<<_dt<<" fac "<<tp_factor<<std::endl;
        clist_[_igroup].slowdown.updatekappa(_tnow+_dt, clist_[_igroup].mass/_mpert, tp_factor, _md_factor);
    }

    void adjustSlowDown(const PS::F64 dt) {
        for (PS::S32 i=0; i<clist_.size(); i++) {
            clist_[i].slowdown.adjustkappa(dt);
        }
    }

    void adjustSlowDownPeriod(const PS::F64 dt, PS::S32* np) {
        for (PS::S32 i=0; i<clist_.size(); i++) {
            np[i] = clist_[i].slowdown.adjustkappaPeriod(dt);
        }
    }
    
    void initialSys() {
        for (PS::S32 i=0; i<clist_.size(); i++) {
            const PS::S32 ipert = pert_disp_[i];
            clist_[i].initSys(0.0, *ARC_control_, &(par_list_.back()), &pert_[ipert], &pforce_[ipert], pert_n_[i]);
        }
    }

    //! Initial one group chain integration parameters
    /*! Initial one group chain integration parameters
      @param[in] _i_group: group to initialize
      @param[in] _time_sys: time to initialize
     */
    void initialOneSys(const PS::S32 _i_group, const PS::F64 _time_sys) {
#ifdef HARD_DEBUG
        assert(!group_mask_map_[_i_group]);
#endif
        const PS::S32 ipert = pert_disp_[_i_group];
        clist_[_i_group].initSys(_time_sys, *ARC_control_, &(par_list_[_i_group]),&pert_[ipert], &pforce_[ipert], pert_n_[_i_group]);
    }

    //! Initial one group chain member and c.m.
    void initialOneChain(const PS::S32 _i_group) {
#ifdef HARD_DEBUG
        assert(!group_mask_map_[_i_group]);
#endif
        clist_[_i_group].initChain();
    }


    void dump(const char* fname, const PS::S32 ic, const PS::F64 time_end, const PS::F64 ds_use) {
        std::FILE* fp = std::fopen(fname,"w");
        if (fp==NULL) {
            std::cerr<<"Error: filename "<<fname<<" cannot be open!\n";
            abort();
        }
        fwrite(&time_end,sizeof(PS::F64),1,fp);
        fwrite(&ds_use,sizeof(PS::F64),1,fp);

        par_list_[ic].dump(fp);
        PS::S32 np = pert_n_[ic];
        fwrite(&np,sizeof(PS::S32),1,fp);
        const PS::S32 ipert = pert_disp_[ic];
        for (PS::S32 i=0;i<np;i++) {
            pert_[ipert+i]->dump(fp);
            pforce_[ipert+i]->dump(fp);
        }

        clist_[ic].dump(fp);
        ARC_control_->dump(fp);
        bininfo[ic].dump(fp);

        std::fclose(fp);
        //clist_[ic].print(std::cerr);
    }

    PS::S64 integrateOneStepSymTwo(const PS::S32 ic, const PS::F64 time_end, const PS::S32 kp) {
#ifdef HARD_DEBUG
        assert(!group_mask_map_[ic]);
#endif
        ARChain* c = &clist_[ic];
        ARC_pert_pars* par = &par_list_[ic];
        PS::F64 ds_use=bininfo[ic].tstep;
        const PS::S32 ipert = pert_disp_[ic];
        PS::F64 timetable[8]; // Notice, assuming sym order is -6
#ifdef ARC_OPT_SYM2
        const PS::F64 m1=c->getP(0).getMass();
        const PS::F64 m2=c->getP(1).getMass();
        const PS::F64 m2_mt = m2/(m1+m2);
        const PS::F64 m1_m2_1 = -m1/m2-1.0;
#endif
        const PS::S32 np=8*kp;
        for (PS::S32 i=0; i<np; i++) {
#ifdef ARC_OPT_SYM2
            c->Symplectic_integration_two(ds_use, *ARC_control_, timetable, m2_mt, m1_m2_1, par, &pert_[ipert], &pforce_[ipert], pert_n_[ic]);
#else 
            c->Symplectic_integration(ds_use, *ARC_control_, timetable, par, &pert_[ipert], &pforce_[ipert], pert_n_[ic]);
#endif
        }
//        std::cout<<std::setprecision(16)<<ds_use<<" "<<kp<<" "
//                 <<c->getTime()<<" "<<time_end<<" "<<getSlowDown(0)
//                 <<std::endl; 
//        for (int j=0; j<c->getN(); j++) {
//            std::cout<<c->getP(j).mass<<" "
//                     <<c->getP(j).pos<<" "
//                     <<c->getP(j).vel<<std::endl;
//        }
#ifdef ARC_WARN       
        if((c->getTime()-time_end)/time_end>1e-6) {
            std::cerr<<"Warning! time not synchronized! t(chain)="<<c->getTime()<<" t="<<time_end<<" diff="<<(c->getTime()-time_end)/time_end<<std::endl;
        }
#endif
        return np;
    }

#ifdef ARC_SYM
    PS::S64 integrateOneStepSym(const PS::S32 ic,
                                const PS::F64 time_end,
                                const PS::F64 dt_limit) {
#ifdef HARD_DEBUG
        assert(!group_mask_map_[ic]);
#endif
        ARChain* c = &clist_[ic];
        ARC_pert_pars* par = &par_list_[ic];
        //PS::F64 ds_up_limit = 0.25*dt_limit/c->calc_dt_X(1.0,*ARC_control_);
        //PS::F64 ds_use = 2.0*bininfo[ic].tstep*std::abs(c->getPt());
        PS::F64 ds_use=bininfo[ic].tstep;
        //PS::F64 ds_use = c->calc_next_step_custom(*ARC_control_,par);
        //if (ds_use>ds_up_limit) ds_use = ds_up_limit;

        const PS::S32 ipert = pert_disp_[ic];
        bool fix_step_flag = false;
        // for high-eccentric binary, it is better to fix step to avoid big step drop, for hyperbolic, fix step is risky
        if(c->getN()==2&&bininfo[ic].ecc>0.99&&bininfo[ic].ecc<1.0) {
            fix_step_flag = true;
            PS::F64 korg=c->slowdown.getkappaorg();
            if(korg<1.0) ds_use *= std::max(0.1,korg);
        }

        PS::S64 stepcount = c->Symplectic_integration_tsyn(ds_use, *ARC_control_, time_end, par, &pert_[ipert], &pforce_[ipert], pert_n_[ic],fix_step_flag, step_count_limit);

#ifdef ARC_WARN
        if(c->info!=NULL) {
            c->info->ErrMessage(std::cerr);
            dump("ARC_dump.dat",ic,time_end,ds_use);
            abort();
        }
#endif
        
#ifdef ARC_DEBUG_DUMP
        if(stepcount<0) {
            dump("ARC_dump.dat",ic,time_end,ds_use);
            std::cerr<<"ic = "<<ic<<" N = "<<c->getN()<<" Np = "<<pert_n_[ic]<<" stepcount = "<<stepcount<<std::endl;
            abort();
        }
#endif
        
        return stepcount;
    }
#else

    PS::S64 integrateOneStepExt(const PS::S32 ic,
                                const PS::F64 time_end,
                                const PS::F64 dt_limit) {
#ifdef HARD_DEBUG
        assert(!group_mask_map_[ic]);
#endif
        ARChain* c = &clist_[ic];
        ARC_pert_pars* par = &par_list_[ic];
        PS::F64 dscoff=1.0;
        PS::F64 ds_up_limit = 0.25*dt_limit/c->calc_dt_X(1.0,*ARC_control_);
        PS::F64 ds_use = c->calc_next_step_custom(*ARC_control_,par);
        //PS::F64 ds_use = 0.5*bininfo[ic].tstep*std::abs(c->GetPt());
        
        if (ds_use>ds_up_limit) ds_use = ds_up_limit;

        PS::S64 nstep=0;

        // convergency check
        PS::S32 converge_count=0;
        PS::S32 error_count=0;
        bool modify_step_flag=false;
        bool final_flag=false;

        while(time_end-c->getTime()>ARC_control_->dterr*c->getTime()) {
            const PS::S32 ipert = pert_disp_[ic];
            PS::F64 dsf=c->extrapolation_integration(ds_use, *ARC_control_, time_end, par, &pert_[ipert], &pforce_[ipert], pert_n_[ic]);
            if (dsf<0) {
                final_flag=true;
                converge_count++;
                if (converge_count>10&&time_end-c->getTime()>ARC_control_->dterr*100) {
                    std::cerr<<"Error: Time synchronization fails!\nStep size ds: "<<ds_use<<"\nEnding physical time: "<<time_end<<"\nTime difference: "<<time_end-c->getTime()<<"\nR_in: "<<Int_pars_->rin<<"\nR_out: "<<Int_pars_->rout<<"\n";
                    dump("ARC_dump.dat",ic,time_end,ds_use);
                    abort();
                }
                else ds_use *= -dsf;
            }
            else if (dsf==0) {
                c->info->ErrMessage(std::cerr);
                error_count++;
                if(error_count>4) {
                    std::cerr<<"Error: Too much error appear!\nStep size ds: "<<ds_use<<"\nEnding physical time: "<<time_end<<"\nTime difference: "<<time_end-c->getTime()<<"\nR_in: "<<Int_pars_->rin<<"\nR_out: "<<Int_pars_->rout<<"\n";
                    dump("ARC_dump.dat",ic,time_end,ds_use);
                    abort();
                }
                if (c->info->status==5) {
                    dscoff = 0.25;
                    ds_use *= dscoff;
                }
                else if (c->info->status==4) ds_use = std::min(dscoff*c->calc_next_step_custom(*ARC_control_, par),ds_up_limit);
                else ds_use *= 0.1;
                modify_step_flag=true;
            }
            else  {
                if (final_flag) {
                    if (converge_count>10&&time_end-c->getTime()>ARC_control_->dterr*100) {
                        std::cerr<<"Error: Time synchronization fails!\nStep size ds: "<<ds_use<<"\nEnding physical time: "<<time_end<<"\nTime difference: "<<time_end-c->getTime()<<"\nR_in: "<<Int_pars_->rin<<"\nR_out: "<<Int_pars_->rout<<"\n";
                        dump("ARC_dump.dat",ic,time_end,ds_use);
                        abort();
                    }
                    converge_count++;
                }
                else if (modify_step_flag&&error_count==0) {
                    ds_use = std::min(dscoff*c->calc_next_step_custom(*ARC_control_, par),ds_up_limit);
                    modify_step_flag=false;
                }
                // reducing error counter if integration success, this is to avoid the significant change of step may cause some issue
                if(error_count>0) error_count--;
            }
        }
#ifdef ARC_PROFILE
        nstep = c->profile.itercount;
#endif

        return nstep;
    }
#endif

    //! Integrate active ARC groups
    /* @param[in] _act_list: active ARC group index list
       @param[in] _n_act: number of active groups
       @param[in] _time_end: end of physical integration time
       @param[in] _dt_limit: physical time step upper limit
     */
    PS::S64 integrateOneStepList(PS::S32 _act_list[],
                                 PS::S32 _n_act,
                                 const PS::F64 _time_end,
                                 const PS::F64 _dt_limit) {
        PS::S64 nstep = 0;
        for(PS::S32 i=0; i<_n_act; i++) {
            if(getMask(i)) continue;
#ifdef ARC_SYM
            nstep += integrateOneStepSym(_act_list[i], _time_end, _dt_limit);
#else
            nstep += integrateOneStepExt(_act_list[i], _time_end, _dt_limit);
#endif
        }

        return nstep;
    }

    //! Integrate active ARC groups
    /* @param[in] _time_end: end of physical integration time
       @param[in] _dt_limit: physical time step upper limit
     */
    PS::S64 integrateOneStepList(const PS::F64 _time_end,
                                 const PS::F64 _dt_limit) {
        PS::S64 nstep = 0;
        for(PS::S32 i=0; i<clist_.size(); i++) {
            if(getMask(i)) continue;
#ifdef ARC_SYM
            nstep += integrateOneStepSym(i, _time_end, _dt_limit);
#else
            nstep += integrateOneStepExt(i, _time_end, _dt_limit);
#endif

        }
        return nstep;
    }

    //! Update the c.m. data from the original particle data
    /*!
      @param[in] _ptcl: original particle array
      @param[in] _ptcl_list: particle index list need to be updated, assume _ptcl and ARC group have consistent index
      @param[in] _n_ptcl: number of particles
     */
    template <class Tptcl>
    void updateCM(Tptcl* _ptcl,
                  PS::S32* _ptcl_list,
                  PS::S32 _n_ptcl) {
        for(PS::S32 i=0; i<_n_ptcl; i++) {
            PS::S32 k = _ptcl_list[i];
#ifdef HARD_DEBUG
            assert(k<clist_.size());
#endif
            clist_[k].pos =  _ptcl[k].pos;
            clist_[k].vel =  _ptcl[k].vel;
            clist_[k].mass = _ptcl[k].mass;
        }
    }

    //! Update the c.m. data from the original particle data
    /*!
      @param[in] _ptcl: original particle array
    */
    template <class Tptcl>
    void updateCM(Tptcl _ptcl[]) {
        for(PS::S32 i=0; i<clist_.size(); i++) {
            clist_[i].pos  = _ptcl[i].pos;
            clist_[i].vel  = _ptcl[i].vel;
            clist_[i].mass = _ptcl[i].mass;
        }
    }

    //! update rsearch of components based on c.m.
    /*!
      @param[in] _dt_tree: tree time step
     */
    void updateRSearch(const PS::F64 _dt_tree) {
        for(PS::S32 i=0; i<clist_.size(); i++) {
            if(getMask(i)) continue;
            clist_[i].calcRSearch(_dt_tree);
            TpARC** ipadr=clist_[i].getPAdr();
            for (PS::S32 k=0; k<clist_[i].getN(); k++)
                ipadr[k]->r_search = clist_[i].r_search;
        }
    }

    //! Shift member ptcls to their c.m. frame
    /*! Shift all group members to their c.m. frame for ARC integration
     */
    void shift2CM() {
        for(PS::S32 i=0; i<clist_.size(); i++) {
            clist_[i].center_shift();
        }
    }

    //! resolve all groups' member particles
    /*! shift compotent coordinates to original frame and save data to original particle address
     */
    void resolve() {
        for(PS::S32 i=0; i<clist_.size(); i++) {
            clist_[i].resolve();
        }
    }

    //! resolve a list of groups
    /*! shift compotent coordinates to original frame and save data to original particle address
      @param[in] _group_list: group list to resovle
      @param[in] _n_group: number of group to resolve
     */
    void resolve(const PS::S32 _group_list[], const PS::S32 _n_group) {
        for(PS::S32 i=0; i<_n_group; i++) {
            clist_[_group_list[i]].resolve();
        }
    }

    //! Check break condition
    /*! Check whether it is necessary to break the chain
      1. get maximum distance (rmax) pair
      2. if rmax>r_crit, check whether the pair is go away or go close
      @param[out] _break_group_list: group index list to break
      @param[out] _break_isplit_list: index in chain list to split for corresponding groups
      @param[in] _r_crit2: distance (square) criterion to check whether need to break
      \return n_group_break: number of groups need to break
     */
    PS::S32 checkBreak(PS::S32* _break_group_list,
                       PS::S32* _break_isplit_list,
                       const PS::F64 _r_crit2) {
        PS::F64 r_max2;
        PS::S32 r_max_index;
        PS::S32 n_group_break=0;
        for (PS::S32 i=0; i<clist_.size(); i++) {
            if(getMask(i)) continue;
            // obtain maximum distance pair
            clist_[i].getRmaxIndex(r_max2, r_max_index);
            if(r_max2>_r_crit2) {
                // check whether outcome or income
                bool out_flag=clist_[i].getDirection(r_max_index);
                if(out_flag) {
                    _break_group_list[n_group_break] = i;
                    _break_isplit_list[n_group_break]= r_max_index;
                    n_group_break++;
                }
            }
        }
        return n_group_break;
    }

    //! Get group original address in chain order
    /*! 
      @param[out] _list particle address array to store the results
      @param[in] _igroup: group index to split
     */
    PS::S32 getPtclAdrChain(TpARC* _list[],
                            const PS::S32 _igroup){
        return clist_[_igroup].getPAdrChain(_list);
    }
    

    //! Get number of members in group i
    /*!
      @param[in] _igroup: group index
      \return number of members in group i
     */
    PS::S32 getGroupN(const PS::S32 _igroup) const {
        return clist_[_igroup].getN();
    }

    //! Get number of groups
    /*!
      \return number of groups
     */
    PS::S32 getNGroups() const {
        return clist_.size();
    }

    //! return integration mask
    /*!
      @param[in] _igroup: index of group
      \return true: suppressed for integration
     */
    bool getMask(const PS::S32 _igroup) const {
        return group_mask_map_[_igroup];
    }

    //! return CM particle (cast pointer to TpARC*)
    /*! 
      @param[in] _igroup: index of group
      \return TpARC* c.m. particle pointer of index i
     */
    TpARC* getCM(const PS::S32 _igroup) {
        return (TpARC*)&clist_[_igroup];
    }

    const TpARC* getGroupPtcl(const PS::S32 i) const {
#ifdef HARD_CHECK_ENERGY
        assert(i<clist_.size());
#endif        
        return &clist_[i].getP(0);
    }

    //! return the address array of member particle 
    /*!
      @param[in] _igroup: index of group
      \return TpARC** particle address array of members
     */
    TpARC** getGroupPtclAdr(const PS::S32 i) const {
#ifdef HARD_CHECK_ENERGY
        assert(i<clist_.size());
#endif        
        return clist_[i].getPAdr();
    }
    
    //! Get slow down factor
    /*! get slow down factor kappa
      @param[in] i: index of groups
      \return slowdown kappa
     */
    PS::F64 getSlowDown(const PS::S32 i) const{
        return clist_[i].slowdown.getkappa();
    }

    //! Get slow down original factor
    /*! get slow down original kappa
     */
    PS::F64 getSlowDownOrg(const PS::S32 i) const{
        return clist_[i].slowdown.getkappaorg();
    }

    //! Print slow down parameters
    /*! Print slow down parameters
      @param[in] _os: ofstream for printing
      @param[in] _i: Chain index
      @param[in] _precision: printed precision for one variable
      @param[in] _width: printing width for one variable
     */
    void printSlowDown(std::ostream& _os, const PS::S32 _i, const PS::S32 _precision=15, const PS::S32 _width=23) {
        clist_[_i].slowdown.print(_os,_precision,_width);
    }


//#ifdef ARC_PROFILE
//    const PS::S64 getNsubstep() const{
//        PS::S64 Nsum = 0;
//        for (int i=0; i<clist_.size(); i++) 
//            Nsum += clist_[i].profile.itercount;
//        return Nsum;
//    }
//#endif

#ifdef ARC_DEBUG_PRINT
    void data_dump(std::ostream& os, const PS::S32 i, const PS::F64 dt_limit) const{
        const ARC_pert_pars* par = &par_list_[i];
        os<<std::setprecision(15)<<dt_limit<<" "
          <<clist_[i].getN()+pert_n_[i]<<" "
          <<par->rin<<" "
          <<par->rout<<" "
          <<par->rout<<" "
          <<par->rin*0.1<<" "
          <<dt_limit<<" "
          <<0.05<<" "
          <<std::sqrt(par->eps2)<<std::endl;
        for (PS::S32 j=0; j<clist_[i].getN(); j++) {
            os<<clist_[i].getP(j).mass<<" "
              <<clist_[i].getP(j).pos<<" "
              <<clist_[i].getP(j).vel<<std::endl;
        }
        for (PS::S32 j=1; j<pert_n_[i]; j++) {
            os<<pert_[pert_disp_[i]+j]->mass<<" "
              <<pert_[pert_disp_[i]+j]->pos<<" "
              <<pert_[pert_disp_[i]+j]->vel<<std::endl;
        }
    }

    //! ARC info print
    /* @param[in] _n_group: current total number of groups already integrated
       @param[in] _n_group_in_cluster: number of groups in current cluster
       @param[in] _n_ptcl: number of real particles
       @param[in] _n_hint: number of Hint particles
       @param[in] _dt_limit: hard time step limit
       @param[in] _kp: kepler period number per step
     */
    bool info_print(std::ostream& os, const PS::S64 _n_group, const PS::S64 _n_group_in_cluster, const PS::S64 _n_ptcl, const PS::S64 _n_hint, const PS::F64 _dt_limit, const PS::S32 _kp, const PS::S32 _n_step_limit=10000) const{
        bool dump_flag=false;
        for (PS::S32 i=0; i<clist_.size(); i++) {
            os<<"ARC_info: "
              <<" i_group_tot="<<_n_group+i
              <<" i_group="<<i
              <<" n_ptcl="<<_n_ptcl
              <<" n_groups="<<_n_group_in_cluster
              <<" n_hint="<<_n_hint
              <<" n_member="<<clist_[i].getN()
              <<" n_pert="<<pert_n_[i]
              <<" semi="<<bininfo[i].semi
              <<" ecc="<<bininfo[i].ecc
              <<" period="<<bininfo[i].peri
              <<" tstep="<<bininfo[i].tstep
              <<" sd="<<clist_[i].slowdown.getkappa();
            PS::S64 nstep = 0;
#ifdef ARC_SYM
            if(_kp>0) nstep = _kp*8;
            else nstep = clist_[i].profile.stepcount[0];
#else
            nstep = clist_[i].profile.itercount;
#endif
            os<<" nstep="<<nstep<<std::endl;

            if(nstep>_n_step_limit) {
                os<<"Data dump for hardtest in the case of large nstep "<<nstep<<std::endl;
                data_dump(os, i, _dt_limit);
                dump_flag = true;
            }
        }
        return dump_flag;
    }
#endif
    
#ifdef HARD_CHECK_ENERGY
    template <class Teng>
    void EnergyRecord(Teng &energy, const PS::S32 sdflag=false) {
        energy.kin = energy.pot = energy.tot = 0.0;
        PS::F64 sd;
        for(PS::S32 i=0; i<clist_.size(); i++) {
            if (sdflag) sd = 1.0/clist_[i].slowdown.getkappa();
            else sd = 1.0;
            energy.kin += sd*clist_[i].getEkin();
            energy.pot += sd*clist_[i].getPot();
            energy.tot += sd*clist_[i].getPt();
        }
    }
#endif                         

#ifdef HARD_DEBUG_PRINT
    template <class Tptcl>
    void writePtcl(FILE* _fout) const{
        for (PS::S32 i=0; i<clist_.size(); i++) {
            PS::F64vec pcm_pos = clist_[i].pos;
            PS::F64vec pcm_vel = clist_[i].vel;
            for (PS::S32 j=0; j<clist_[i].getN(); j++) {
                Tptcl pj = clist_[i].getP(j);
                pj.pos += pcm_pos;
                pj.vel += pcm_vel;
                pj.ParticleBase::writeAscii(_fout);
            }
        }
    }
#endif

};
