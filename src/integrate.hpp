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

//const PS::F64 SAFTY_FACTOR_FOR_SEARCH_SQ;

//!leap frog kick for single----------------------------------------------
/* modify the velocity of particle in global system
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
    for(int i=0; i<n; i++){
        const PS::S32 k=_adr[i];
        _sys[k].vel  += _sys[k].acc * _dt;
    }
}

//!leap frog kick for clusters------------------------------------------
/* modify the velocity of particle in local, if remote non-group particle, do nothing, need MPI receive to update data
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
    for(int i=0; i<n; i++) {
        const PS::S64 cm_adr=_ptcl[i].status;
        const PS::S64 i_adr =_ptcl[i].adr_org;
        // if is group member, recover mass and kick due to c.m. force
        if(cm_adr>0) {
            _ptcl[i].mass = _ptcl[i].mass_bk;
            _ptcl[i].vel += _sys[cm_adr].acc * _dt;
        }
        else if(i_adr>=0) { // non-member particle
            // not remote particles
            _ptcl[i].vel += _sys[i_adr].acc * _dt;
        }
    }
}

template<class Tpsys, class Ttree>
void Drift(Tpsys & system,
           const Ttree & tree,
           const PS::F64 dt){
    const PS::S32 n = system.getNumberOfParticleLocal();
#pragma omp parallel for
    for(int i=0; i<n; i++){
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

template <class Tptcl>
void softKickForOneGroup(Tptcl * ptcl_org,
                         const PS::S32  i_cm,
                         const PS::S32* group_list,
                         const PS::S32  group_n,
                         const PS::S32* soft_pert_list,
                         const PS::F64  dt_soft,
                         const PS::S32  n_split) {
    Tptcl* pi = &ptcl_org[i_cm];
    PS::F64vec fi= PS::F64vec(0.0);
    PS::F64 micum = 0.0;
#ifdef TIDAL_TENSOR
    for (PS::S32 j=8; j<2*n_split; j++) {
#else
    for (PS::S32 j=0; j<2*n_split; j++) {
#endif
        Tptcl* pj = &ptcl_org[soft_pert_list[j]];
        fi += pj->mass*pj->vel; // here pj->vel store the soft force of fake members
        micum += pj->mass;
#ifdef HARD_DEBUG
        assert(((pj->status)>>ID_PHASE_SHIFT)==-pi->id);
#endif
    }
    PS::F64vec vkick = fi/micum * dt_soft;

#ifdef HARD_DEBUG
    assert(abs(micum-pi->mass)<1e-10);
#endif
    for (PS::S32 i=0; i<group_n; i++) {
        Tptcl* pk = &ptcl_org[group_list[i]];
        pk->vel += vkick;
    }
    pi->vel += vkick;
}

class PtclH4: public Ptcl{
public:
    PS::F64vec acc0;
    PS::F64vec acc1;
    PS::F64 dt;
    PS::F64 time;
#ifdef HARD_DEBUG
    PS::F64vec acc2; // for debug
    PS::F64vec acc3; // for debug
#endif

    PtclH4(): dt(0.0), time(0.0) {}
    
    template<class Tp>
    PtclH4(const Tp &p): Ptcl(p), dt(0.0), time(0.0) {}

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

template<class Tptcl>
class HermiteIntegrator{
private:
    PS::ReallocatableArray<PtclPred> pred_;
    PS::ReallocatableArray<PtclForce> force_;
    PS::ReallocatableArray<PS::S32> adr_sorted_;
    PS::ReallocatableArray<PS::F64> time_next_;

    PS::ReallocatableArray<PtclH4> ptcl_;   // first c.m.; second single
    PS::ReallocatableArray<PS::S32> Jlist_;
    PS::ReallocatableArray<PS::S32> Jlist_disp_;
    PS::ReallocatableArray<PS::S32> Jlist_n_;

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
        for (int i=0; i<n; i++) {
            PS::F64 mi=ptcl[i].mass;
            pcm.mass += mi;
            pcm.pos += mi*ptcl[i].pos;
            pcm.vel += mi*ptcl[i].vel;
        }
        PS::F64 invm = 1.0/pcm.mass;
        pcm.pos *= invm;
        pcm.vel *= invm;
        
        for (int i=0;i<n;i++) {
//            std::cerr<<std::setprecision(14)<<"cm sbore i "<<i<<" pos "<<ptcl[i].pos<<" vel "<<ptcl[i].vel<<std::endl;
            ptcl[i].pos -= pcm.pos;
            ptcl[i].vel -= pcm.vel;
//#ifdef HARD_DEBUG
//            std::cerr<<"cm shift i "<<i<<" pos "<<ptcl[i].pos<<" vel "<<ptcl[i].vel<<std::endl;
//#endif
        }
    }

    void shiftBackCM(const ParticleBase &pcm, PtclH4* ptcl, const PS::S32 n) {
#ifdef HARD_DEBUG
        PS::F64 m=0;
        for (int i=0;i<n;i++) m+= ptcl[i].mass;
        assert(m==pcm.mass);
        //PS::F64vec pos=0,vel=0;
        //for (int i=0;i<n;i++) {
        //    pos += ptcl[i].mass*ptcl[i].pos;
        //    vel += ptcl[i].mass*ptcl[i].vel;
        //}
        //pos /=m;
        //vel /=m;
        //std::cerr<<std::setprecision(14)<<"cm pos "<<pos<<" vel "<<vel<<std::endl;
#endif
        for (int i=0;i<n;i++) {
            //std::cerr<<"cm bbck i "<<i<<" pos "<<ptcl[i].pos<<" vel "<<ptcl[i].vel<<std::endl;
            ptcl[i].pos += pcm.pos;
            ptcl[i].vel += pcm.vel;
            //std::cerr<<"cm back i "<<i<<" pos "<<ptcl[i].pos<<" vel "<<ptcl[i].vel<<std::endl;
        }
    }

    //
    /* \return fail_flag: if the time step < dt_min, return true (failure)
     */
    template <class Tp>
    bool CalcBlockDt2ndAct(Tp ptcl[],
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

    template <class Tp, class ARCint>
    inline void CalcAcc0Acc1Act(PtclForce force[],
                         const Tp ptcl[],
                         const PS::S32 Ilist[],
                         const PS::S32 n_act,
                         const PS::S32 Jlist[],
                         const PS::S32 Jlist_disp[],
                         const PS::S32 Jlist_n[],
                         const PS::F64 rin,
                         const PS::F64 rout,
                         const PS::F64 r_oi_inv,
                         const PS::F64 r_A,
                         const PS::F64 eps2,
                         const ARCint *Aint) {
        PS::S32 nbin = 0;
        if (Aint!=NULL) nbin = Aint->getN();
        // PS::ReallocatableArray< std::pair<PS::S32, PS::S32> > & merge_pair ){
        // active particles
        for(PS::S32 i=0; i<n_act; i++){
            const PS::S32 iadr = Ilist[i];
            force[iadr].acc0 = force[iadr].acc1 = 0.0;
            // all
            const PS::S32 joff = Jlist_disp[iadr];
            const PS::S32 jn = Jlist_n[iadr];
            for(PS::S32 j=0; j<jn; j++){
                PS::S32 jadr = Jlist[joff+j];
#ifdef HARD_DEBUG
                assert(iadr!=jadr);
                PS::F64 mcmcheck =0.0;
#endif
                if (jadr<nbin) {
                    const Tptcl* pj = Aint->getGroupPtcl(jadr);
                    PS::F64 sd = Aint->getSlowDown(jadr);
                    for(PS::S32 k=0; k<Aint->getGroupN(jadr); k++) {
                        PS::F64 r2 = 0.0;
                        CalcAcc0Acc1R2Cutoff(ptcl[iadr].pos, ptcl[iadr].vel,
                                             force[iadr].acc0, force[iadr].acc1, r2,
                                             pj[k].pos, pj[k].vel/sd, pj[k].mass,
                                             eps2, rout, rin, r_oi_inv, r_A);
#ifdef HARD_DEBUG
                        mcmcheck += pj[k].mass;
                        //std::cerr<<k<<" P "<<pj[k].pos<<" v "<<pj[k].vel<<" sd "<<sd<<std::endl;
#endif
                    }
#ifdef HARD_DEBUG
                    assert(mcmcheck==ptcl[jadr].mass);
                    assert(mcmcheck>0.0);
#endif                    
                }
                else {
                    PS::F64 r2 = 0.0;
                    //PS::F64 rout = std::max(ptcl[iadr].r_out, ptcl[jadr].r_out);
                    CalcAcc0Acc1R2Cutoff(ptcl[iadr].pos, ptcl[iadr].vel,
                                         force[iadr].acc0, force[iadr].acc1, r2,
                                         ptcl[jadr].pos, ptcl[jadr].vel, ptcl[jadr].mass,
                                         eps2, rout, rin, r_oi_inv, r_A);
                }
                // if(r2 < ((ptcl[adr].r_merge + ptcl[j].r_merge)*(ptcl[adr].r_merge + ptcl[j].r_merge)) && ptcl[j].mass > 0.0){
                //     merge_pair.push_back( std::make_pair(adr, j) );
                // }
            }
        }
    }

    template <class Tpi, class Tp, class ARCint>
    inline void CalcAcc0Acc1AllJ(PtclForce &force, 
                          const Tpi &pi,
                          const PS::S32 iadr,
                          const PS::F64vec &vcmsdi,
                          const PS::F64 sdi,
                          const Tp ptcl[],
                          const PS::S32 n_tot,
                          const PS::S32 nbin,
                          const PS::F64 rin,
                          const PS::F64 rout,
                          const PS::F64 r_oi_inv,
                          const PS::F64 r_A,
                          const PS::F64 eps2,
                          const ARCint* Aint = NULL) {

        for(PS::S32 j=0; j<nbin; j++) {
            if(iadr==j) continue;
#ifdef HARD_DEBUG
            PS::F64 mcmcheck =0.0;
#endif
            const Tptcl* pj = Aint->getGroupPtcl(j);
            PS::F64 sdj = 1.0/Aint->getSlowDown(j);
            PS::F64vec vcmsdj = (1.0-sdj)*ptcl[j].vel;
            for(PS::S32 k=0; k<Aint->getGroupN(j); k++) {
                PS::F64 r2 = 0.0;
                CalcAcc0Acc1R2Cutoff(pi.pos, pi.vel*sdi+vcmsdi,
                                     force.acc0, force.acc1, r2,
                                     pj[k].pos, pj[k].vel*sdj+vcmsdj, pj[k].mass,
                                     eps2, rout, rin, r_oi_inv, r_A);
#ifdef HARD_DEBUG
                mcmcheck += pj[k].mass;
#endif
            }
#ifdef HARD_DEBUG
            assert(mcmcheck==ptcl[j].mass);
            assert(mcmcheck>0.0);
#endif                    
        }
        for(PS::S32 j=nbin; j<n_tot; j++){
            if(iadr==j) continue;
            PS::F64 r2 = 0.0;
            //PS::F64 rout = std::max(ptcl[iadr].r_out, ptcl[j].r_out);
            CalcAcc0Acc1R2Cutoff(pi.pos, pi.vel*sdi+vcmsdi,
                                 force.acc0, force.acc1, r2,
                                 ptcl[j].pos, ptcl[j].vel, ptcl[j].mass,
                                 eps2, rout, rin, r_oi_inv, r_A);
                // if(r2 < ((ptcl[adr].r_merge + ptcl[j].r_merge)*(ptcl[adr].r_merge + ptcl[j].r_merge)) && ptcl[j].mass > 0.0){
                //     merge_pair.push_back( std::make_pair(adr, j) );
                // }
        }
        
    }


    template <class Tp, class ARCint>
    void CalcAcc0Acc1Act(PtclForce force[],
                         const Tp ptcl[],
                         const PS::S32 n_tot,
                         const PS::S32 Ilist[],
                         const PS::S32 n_act,
                         const PS::F64 rin,
                         const PS::F64 rout,
                         const PS::F64 r_oi_inv,
                         const PS::F64 r_A,
                         const PS::F64 eps2,
                         const ARCint* Aint=NULL) {
        PS::S32 nbin = 0;
        if (Aint!=NULL) nbin = Aint->getN();
        static thread_local const PS::F64vec vzero = PS::F64vec(0.0);
        // PS::ReallocatableArray< std::pair<PS::S32, PS::S32> > & merge_pair ){
        // active particles
        for(PS::S32 i=0; i<n_act; i++){
            const PS::S32 iadr = Ilist[i];
            force[iadr].acc0 = force[iadr].acc1 = 0.0;
            
            if (iadr<nbin) {
#ifdef HARD_DEBUG
                PS::F64 mcmcheck =0.0;
#endif
                const PS::S32 ni = Aint->getGroupN(iadr);
                const Tptcl* pi = Aint->getGroupPtcl(iadr);
                const PS::F64 sdi = 1.0/Aint->getSlowDown(iadr);
                const PS::F64vec vcmsdi = (1.0-sdi)*ptcl[iadr].vel;
                PtclForce fp[ni];
                for (PS::S32 j=0; j<ni; j++) {
                    fp[j].acc0 = fp[j].acc1 = 0.0;
                    CalcAcc0Acc1AllJ(fp[j], pi[j], iadr, vcmsdi, sdi, ptcl, n_tot, nbin, rin, rout, r_oi_inv, r_A, eps2, Aint);
                    force[iadr].acc0 += pi[j].mass*fp[j].acc0;
                    force[iadr].acc1 += pi[j].mass*fp[j].acc1;
#ifdef HARD_DEBUG
                    mcmcheck += pi[j].mass;
#endif
                }
#ifdef HARD_DEBUG
                assert(mcmcheck==ptcl[iadr].mass);
                assert(mcmcheck>0.0);
#endif                    
                force[iadr].acc0 /= ptcl[iadr].mass;
                force[iadr].acc1 /= ptcl[iadr].mass;
            }
            else CalcAcc0Acc1AllJ(force[iadr], ptcl[iadr], iadr, vzero, 1.0, ptcl, n_tot, nbin, rin, rout, r_oi_inv, r_A, eps2, Aint);
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

    void SortAndSelectIp(PS::S32 adr_sorted[],
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
          const PS::S32 adr = adr_sorted[ip];
          time_next[adr] += ptcl[adr].dt; // n_act only
          }
        */
        std::sort(adr_sorted, adr_sorted+n_act, SortAdr(time_next));

        const PS::F64 time_ref = time_next[adr_sorted[0]];
        group_n = 0;
        for(n_act=1; n_act<n_tot; n_act++){
            if(time_ref < time_next[adr_sorted[n_act]]) {
                break;
            }
            if(n_act<n_limit) group_list[group_n++] = n_act;
        }
    }

    void SortAndSelectIp(PS::S32 adr_sorted[],
                         PS::F64 time_next[],
                         PS::S32 & n_act,
                         const PS::S32 n_tot){
        // const PS::S32 n_tot = time_next.size();
        //std::cerr<<"before sort"<<std::endl;
        /*
          for(PS::S32 ip=0; ip<ni_old; ip++){
          const PS::S32 adr = adr_sorted[ip];
          time_next[adr] += ptcl[adr].dt; // n_act only
          }
        */
        std::sort(adr_sorted, adr_sorted+n_act, SortAdr(time_next));

        const PS::F64 time_ref = time_next[adr_sorted[0]];
        for(n_act=1; n_act<n_tot; n_act++){
            if(time_ref < time_next[adr_sorted[n_act]]) {
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

    /* \return fail_flag: if the time step < dt_min, return true (failure)
     */
    bool CorrectAndCalcDt4thAct(PtclH4 ptcl[],
                                const PtclForce force[],
                                const PS::S32 adr_sorted[], 
                                const PS::S32 n_act,
                                const PS::F64 dt_max,
                                const PS::F64 dt_min,
                                const PS::F64 a0_offset_sq,
                                const PS::F64 eta){
        bool fail_flag=false;
        static thread_local const PS::F64 inv3 = 1.0 / 3.0;
        for(PS::S32 i=0; i<n_act; i++){
            const PS::S32 adr = adr_sorted[i];
            PtclH4*     pti = &ptcl[adr];
            const PtclForce* fpi = &force[adr];

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
            const PS::F64 dt_ref = CalcDt4th(pti->acc0, pti->acc1, acc2, acc3, eta, a0_offset_sq);

            const PS::F64 dt_old = pti->dt;
#ifdef HARD_DEBUG
            // for debug
            assert(dt_old != 0.0);
            pti->acc2 = acc2;
            pti->acc3 = acc3;
#endif
            pti->dt = dt_max;

#ifdef FIX_STEP_DEBUG
            pti->dt /= STEP_DIVIDER;
#else
            while(pti->dt > dt_ref) pti->dt *= 0.5;
            pti->dt = dt_old*2 < pti->dt ?  dt_old*2 : pti->dt;
#endif

#ifdef HARD_DEBUG
            assert(pti->dt != 0.0);
//            assert(pti->dt >1.0e-12);
#endif

            if(pti->dt <dt_min) {
                std::cerr<<"Error: Hermite integrator step size ("<<pti->dt<<") < dt_min ("<<dt_min<<")!"<<std::endl;
                std::cerr<<" pti->time="<<pti->time<<" i="<<i<<" adr="<<adr<<" pos="<<pti->pos<<" vel="<<pti->vel<<" acc="<<pti->acc0<<" acc1="<<pti->acc1
#ifdef HARD_DEBUG
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

public:

//    HermiteIntegrator(const PS::S32 n) {
//#ifdef HARD_DEBUG
//        assert(pred_.size()==0);
//        assert(force_.size()==0);
//        assert(adr_sorted_.size()==0);
//        assert(time_next_.size()==0);
//#endif
//        resizeArray(n);
//    }
    
    void moveCM(const PS::F64 dt) {
        pcm_.pos += pcm_.vel * dt;
        //std::cerr<<"pcm.pos "<<pcm_.pos<<" pcm.vel "<<pcm_.vel<<" dt "<<dt<<std::endl;
    }
    
    void shiftBackCM() {
        shiftBackCM(pcm_,ptcl_.getPointer(),ptcl_.size());
    }

    void resizeArray(const PS::S32 n) {
        ptcl_.reserve(n);
        Jlist_disp_.reserve(n);
        Jlist_n_.reserve(n);
        pred_.reserve(n);
        force_.reserve(n);
        adr_sorted_.reserve(n);
        time_next_.reserve(n);
    }

    void setPtcl(Tptcl * ptcl_org, 
                 const PS::S32 n_org, 
                 const PS::S32 ptcl_list[],
                 const PS::S32 ptcl_n) {
#ifdef HARD_DEBUG
        assert(ptcl_n<=n_org);
#endif
        resizeArray(ptcl_n);
        for (int i=0; i<ptcl_n; i++) {
            ptcl_.increaseSize(1);
            ptcl_.back().DataCopy(ptcl_org[ptcl_list[i]]);
            //ptcl_.pushBackNoCheck(ptcl_org[ptcl_list[i]]);
        }
    }

    void writeBackPtcl(Tptcl * ptcl_org, 
                 const PS::S32 n_org, 
                 const PS::S32 ptcl_list[],
                 const PS::S32 ptcl_n) {
#ifdef HARD_DEBUG
        assert(ptcl_n<=n_org);
#endif
        for (int i=0; i<ptcl_n; i++) {
            ptcl_org[ptcl_list[i]].DataCopy(ptcl_[i]);
        }
    }

    void searchPerturber(PS::F64 dr_search[], const PS::S32 nbin) {
        // find perturber
        PS::S32 n = ptcl_.size();
        PS::S32 nj_tot = 0;
        for(int i=0; i<n; i++) {
#ifdef HARD_DEBUG
            assert(Jlist_disp_.size()==i);
            assert(Jlist_n_.size()==i);
#endif
            Jlist_disp_.pushBackNoCheck(nj_tot);
            Jlist_n_.pushBackNoCheck(0);
            for(int j=0; j<nbin; j++) {
                PS::F64vec dr = ptcl_[i].pos-ptcl_[j].pos;
                PS::F64 r2 = dr*dr;
                PS::F64 r_search = std::max(ptcl_[i].r_search,ptcl_[j].r_search)*SAFTY_FACTOR_FOR_SEARCH + SAFTY_OFFSET_FOR_SEARCH + dr_search[j];
                PS::F64 r_search2 = r_search*r_search;
                if (r2<r_search2&&i!=j) {
                    Jlist_.push_back(j);
                    Jlist_n_[i]++;
                    nj_tot++;
                }
            }
            for(int j=nbin; j<n; j++) {
                PS::F64vec dr = ptcl_[i].pos-ptcl_[j].pos;
                PS::F64 r2 = dr*dr;
                PS::F64 r_search = std::max(ptcl_[i].r_search,ptcl_[j].r_search)*SAFTY_FACTOR_FOR_SEARCH + SAFTY_OFFSET_FOR_SEARCH;
                PS::F64 r_search2 = r_search*r_search;
                if (r2<r_search2&&i!=j) {
                    Jlist_.push_back(j);
                    Jlist_n_[i]++;
                    nj_tot++;
                }
            }
        }
#ifdef HARD_DEBUG
        for(int i=0; i<n-1; i++) {
            assert(Jlist_n_[i]==Jlist_disp_[i+1]-Jlist_disp_[i]);
            assert(Jlist_n_[i]<=n);
        }
        assert(Jlist_disp_[n-1]+Jlist_n_[n-1]==Jlist_.size());
#endif
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

    void setParams(//const PS::F64 dt_limit_hard,  // time step limit
                   const PS::F64 eta_s,          // time step parameter
                   const PS::F64 r_in,           // force parameter
                   const PS::F64 r_out,          // force parameter
                   const PS::F64 eps_sq){        // softening parameter
        // dt_limit_hard_= dt_limit_hard; 
        eta_s_        = eta_s;         
        r_in_         = r_in;          
        r_out_        = r_out;
        eps_sq_       = eps_sq;
        r_oi_inv_     = 1.0/(r_out-r_in);
        r_A_          = (r_out-r_in)/(r_out+r_in);
    }

    template <class ARCint>
    bool initialize(PS::F64 dt_max,
                    PS::F64 dt_min,
                    PS::S32 group_act_list[],
                    PS::S32 &group_act_n,
                    const PS::S32 n_groups,
                    ARCint* Aint = NULL,
                    const bool calc_full_flag = true) {
        PS::S32 n_ptcl = ptcl_.size();
        pred_.resizeNoInitialize(n_ptcl);
        force_.resizeNoInitialize(n_ptcl);
        adr_sorted_.resizeNoInitialize(n_ptcl);
        time_next_.resizeNoInitialize(n_ptcl);
        
        PS::F64 mass_min = PS::LARGE_FLOAT;
        //PS::F64 rout_min = PS::LARGE_FLOAT;
        for(PS::S32 i=0; i<n_ptcl; i++){
            // pred[i].mass = ptcl[i].mass = ptcl_org[i].mass;
            // pred[i].pos  = ptcl[i].pos  = ptcl_org[i].pos;
            // pred[i].vel  = ptcl[i].vel  = ptcl_org[i].vel;
            // ptcl[i].setRMerge();
            // pred[i].r_merge = ptcl[i].r_merge;
            // ptcl[i].id = ptcl_org[i].id;
            pred_[i].r_search = ptcl_[i].r_search;
            adr_sorted_[i] = i;
            ptcl_[i].time = ptcl_[i].dt = 0.0;
            time_next_[i] = 0.0;
            if(mass_min > ptcl_[i].mass)  mass_min = ptcl_[i].mass;
            //if(rout_min > ptcl_[i].r_out) rout_min = ptcl_[i].r_out;
        }

        a0_offset_sq_ = 0.1 * mass_min / (r_out_ * r_out_);
        n_act_ = n_ptcl;

        //shift c.m.
        calcCMAndShift(pcm_, ptcl_.getPointer(), n_ptcl);

        if(Aint!=NULL) {
            Aint->updateCM(ptcl_.getPointer());
            Aint->resolve();
        }

        if(calc_full_flag) CalcAcc0Acc1Act(force_.getPointer(), ptcl_.getPointer(), ptcl_.size(), adr_sorted_.getPointer(), n_act_, r_in_, r_out_, r_oi_inv_, r_A_, eps_sq_, Aint);
        else CalcAcc0Acc1Act(force_.getPointer(), ptcl_.getPointer(), adr_sorted_.getPointer(), n_act_, Jlist_.getPointer(), Jlist_disp_.getPointer(), Jlist_n_.getPointer(), r_in_, r_out_, r_oi_inv_, r_A_, eps_sq_, Aint);
    
        // store predicted force
        for(PS::S32 i=0; i<n_ptcl; i++){
            ptcl_[i].acc0 = force_[i].acc0;
            ptcl_[i].acc1 = force_[i].acc1;
        }

        if(Aint!=NULL) Aint->shift();

        bool fail_flag=CalcBlockDt2ndAct(ptcl_.getPointer(), force_.getPointer(), adr_sorted_.getPointer(), n_act_, 0.01*eta_s_, dt_max, dt_min, a0_offset_sq_);

        for(PS::S32 i=0; i<n_ptcl; i++){
            time_next_[i] = ptcl_[i].time + ptcl_[i].dt;
        }
        //SortAndSelectIp(adr_sorted_.getPointer(), time_next_.getPointer(), n_act_, time_next_.size(), group_act_list, group_act_n, n_groups);
        SortAndSelectIp(adr_sorted_.getPointer(), time_next_.getPointer(), n_act_, time_next_.size());

        return fail_flag;
    }
    
    template<class Energy>
    void CalcEnergy(Energy & eng) {
        CalcEnergy(ptcl_.getPointer(), ptcl_.size(), eng, r_in_, r_out_, eps_sq_);
    }

    PS::F64 getNextTime() {
        return time_next_[adr_sorted_[0]];
    }

    /*
      \return If step size < dt_min, return true
     */
    template <class ARCint>
    bool integrateOneStep(const PS::F64 time_sys,
                          const PS::F64 dt_max,
                          const PS::F64 dt_min,
                          const bool calc_full_flag = true,
                          ARCint* Aint = NULL) {
        // pred::mass,pos,vel updated
        PredictAll(pred_.getPointer(), ptcl_.getPointer(), ptcl_.size(), time_sys);
        if(Aint!=NULL) {
            Aint->updateCM(pred_.getPointer());
            Aint->resolve();
        }
        // force::acc0,acc1 updated
        if(calc_full_flag) CalcAcc0Acc1Act(force_.getPointer(), pred_.getPointer(), ptcl_.size(), adr_sorted_.getPointer(), n_act_, r_in_, r_out_, r_oi_inv_, r_A_, eps_sq_, Aint);
        else CalcAcc0Acc1Act(force_.getPointer(), pred_.getPointer(), adr_sorted_.getPointer(), n_act_, Jlist_.getPointer(), Jlist_disp_.getPointer(), Jlist_n_.getPointer(), r_in_, r_out_, r_oi_inv_, r_A_, eps_sq_, Aint);
        // ptcl_org::pos,vel; pred::time,dt,acc0,acc1,acc2,acc3 updated
        bool fail_flag=CorrectAndCalcDt4thAct(ptcl_.getPointer(), force_.getPointer(), adr_sorted_.getPointer(), n_act_, dt_max, dt_min, a0_offset_sq_, eta_s_);

        for(PS::S32 i=0; i<n_act_; i++){
            PS::S32 adr = adr_sorted_[i];
            time_next_[adr] = ptcl_[adr].time + ptcl_[adr].dt;
        }

        if(Aint!=NULL) {
            Aint->shift();
            //Aint->updateCM(ptcl_.getPointer());
            //Aint.updateCM(Hint.getPtcl(), group_act_list.getPointer(), group_act_n);
        }
        
        return fail_flag;
    }

    void SortAndSelectIp(PS::S32 group_act_list[],
                         PS::S32 &group_act_n,
                         const PS::S32 n_groups) {
        SortAndSelectIp(adr_sorted_.getPointer(), time_next_.getPointer(), n_act_, time_next_.size(), group_act_list, group_act_n, n_groups);
        
    }

    void SortAndSelectIp() {
        SortAndSelectIp(adr_sorted_.getPointer(), time_next_.getPointer(), n_act_, time_next_.size());
    }
    
    PS::S32 getNact() const{
        return n_act_;
    }

    const PS::S32* getActList() const{
        return adr_sorted_.getPointer();
    }

#ifdef HARD_DEBUG
    void printStepHist(){
        std::map<PS::F64, PS::S32> stephist;
        for(int i=0; i<pred_.size(); i++) {
            std::map<PS::F64, PS::S32>::iterator p = stephist.find(ptcl_[i].dt);
            if (p==stephist.end()) stephist[ptcl_[i].dt]=1;
            else stephist[ptcl_[i].dt]++;
        }
        std::cerr<<"Step hist:\n";
        for(auto i=stephist.begin(); i!=stephist.end(); i++) {
            std::cerr<<std::setw(14)<<i->first;
        }
        std::cerr<<std::endl;
        for(auto i=stephist.begin(); i!=stephist.end(); i++) {
            std::cerr<<std::setw(14)<<i->second;
        }
        std::cerr<<std::endl;
    }
#endif

};

#ifdef ISOLATED
//few-body----------------------------------------------
template<class Tptcl, class ARC_par>
void Isolated_Multiple_integrator(Tptcl * ptcl_org,
                                  const PS::S32 n_ptcl,
                                  const PS::F64 time_end,
                                  const PS::F64 dt_limit,
                                  const PS::F64 rout_single,
                                  const PS::F64 gamma,
                                  const PS::F64 m_average,
#ifdef ARC_ERROR
                                  PS::F64 &ARC_error_relative,
                                  PS::F64 &ARC_error,
                                  PS::S32 N_count[20],
#endif                         
                                  const ARC::chainpars &ARC_control,
                                  ARC_par &Int_pars) {
    // kepler motion test
    if (n_ptcl==2) {
        PS::F64 ax=0,ecc;
        PS::F64 inc,OMG,omg,tperi;
        PS::F64 ecc_anomaly_old = PosVel2OrbParam(ax, ecc, inc, OMG, omg, tperi,
                                                  ptcl_org[0].pos, ptcl_org[1].pos, ptcl_org[0].vel, ptcl_org[1].vel, ptcl_org[0].mass, ptcl_org[1].mass);
#ifdef HARD_DEBUG
        std::cerr<<"n_ptcl="<<n_ptcl<<"; ax="<<ax<<"; ecc="<<ecc<<"; peri="<<tperi<<"; pid="<<ptcl_org[0].id<<std::endl;
#endif
        if (ax>0.0&&2.0*ax<Int_pars.rin) {
            // center-of-mass
            Tptcl pcm;
            calc_center_of_mass(pcm, ptcl_org, n_ptcl, true);

            DriveKeplerOrbParam(ptcl_org[0].pos, ptcl_org[1].pos, ptcl_org[0].vel, ptcl_org[1].vel,
                                ptcl_org[0].mass, ptcl_org[1].mass, time_end, ax, ecc, inc, OMG, omg, ecc_anomaly_old);
          
          
            // integration of center-of-mass
            pcm.pos += pcm.vel * time_end;

            center_of_mass_correction(pcm, ptcl_org, n_ptcl);
          
            //  PosVel2OrbParam(ax,ecc,inc,OMG,omg,tperi,ptcl_org[0].pos, ptcl_org[1].pos, ptcl_org[0].vel, ptcl_org[1].vel, ptcl_org[0].mass, ptcl_org[1].mass);
            //  std::cerr<<"A:n_ptcl="<<n_ptcl<<"; ax="<<ax<<"; ecc="<<ecc<<"; peri="<<tperi<<"; pid"<<ptcl_org[0].id<<std::endl;
            //  if (!kout.is_open()) kout.open("kout");
            //  for (int i=0;i<n_ptcl;i++) kout<<std::setprecision(17)<<time_origin_<<" "<<ptcl_org[i].mass<<" "<<ptcl_org[i].pos<<" "<<ptcl_org[i].vel<<std::endl;

            PS::F64 peri = ax*(1+ecc)*gamma*std::pow(pcm.mass/m_average,0.3333);
            if (peri>1.2*ptcl_org[0].r_out || (peri>0 && peri<0.8*ptcl_org[0].r_out) || ptcl_org[0].r_out!= ptcl_org[1].r_out)
                ptcl_org[0].r_out = ptcl_org[1].r_out = std::max(peri,rout_single);
            return;
        }
    }
    else if(n_ptcl==3) {
        // Not yet implementd
    }
      
    ARC::chain<Tptcl> c((std::size_t)n_ptcl);
    static thread_local PS::F64 time_sys = 0.0;

    c.addP(n_ptcl,ptcl_org);

#ifdef HARD_DEBUG
    for (PS::S32 i=0; i<n_ptcl; i++) 
        if (ptcl_org[i].r_out<Int_pars.rin) {
            std::cerr<<"Error, updated p["<<i<<"].rout ("<<ptcl_org[i].r_out<<") < r_in ("<<Int_pars.rin<<")"<<std::endl;
            abort();
        }
#endif

    //c.link_int_par(Int_pars);
    c.init(time_sys,ARC_control,Int_pars);

#ifdef ARC_ERROR
    PS::F64 ARC_error_once = c.getPot()+c.getEkin();
    if(n_ptcl<=20) N_count[n_ptcl-1]++;
    else std::cerr<<"Large cluster formed, n="<<n_ptcl<<std::endl;
#endif
      
    PS::F64 dscoff=1.0;
    PS::F64 ds_up_limit = 0.25*dt_limit/c.calc_dt_X(1.0);
    PS::F64 ds_use = c.calc_next_step_custom();
      
    if (ds_use>ds_up_limit) ds_use = ds_up_limit;

    // convergency check
    PS::S32 converge_count=0;
    PS::S32 error_count=0;
    bool modify_step_flag=false;
    bool final_flag=false;
      
    while(time_end-c.getTime()>ARC_control.dterr) {

        // if (ptcl_org[0].id==8&&time_origin_==0.0078125) {
        //	 std::cout<<"ds= "<<ds_use<<" toff= "<<time_end<<std::endl;
			
        //	 FILE* fout=fopen("data","w");
        //	 fwrite(ptcl_org,sizeof(Tptcl),n_ptcl,fout);
        //	 fclose(fout);
        // for (PS::S32 i=0; i<n_ptcl; i++) 
        //	 std::cout<<ptcl_org[i].getMass()<<" "<<ptcl_org[i].getPos()[0]<<" "<<ptcl_org[i].getPos()[1]<<" "<<ptcl_org[i].getPos()[2]<<" "<<ptcl_org[i].getVel()[0]<<" "<<ptcl_org[i].getVel()[1]<<" "<<ptcl_org[i].getVel()[2]<<std::endl;
        // }

        PS::F64 dsf=c.extrapolation_integration(ds_use,ARC_control,time_end,Int_pars);

        // std::cerr<<"Particle=";
        // for (PS::S32 i=0; i<n_ptcl; i++) std::cerr<<ptcl_org[i].id<<" ";
        // std::cerr<<"n="<<n_ptcl<<" Time_end="<<time_end<<" ctime="<<c.getTime()<<" diff="<<time_end-c.getTime()<<" ds="<<ds_use<<" dsf="<<dsf<<std::endl;

        if (dsf<0) {
            final_flag=true;
            converge_count++;
            if (converge_count>5&&time_end-c.getTime()>ARC_control.dterr*100) {
                std::cerr<<"Error: Time synchronization fails!\nStep size ds: "<<ds_use<<"\nEnding physical time: "<<time_end<<"\nTime difference: "<<time_end-c.getTime()<<"\nR_in: "<<Int_pars.rin<<"\nR_out: "<<Int_pars.rout<<"\n";
                //  ds_use = 0.1*c.calc_next_step_custom();
                //  std::cerr<<"New step size: "<<ds_use<<std::endl;
                //  modify_step_flag=true;
                //	converge_count=0;
                c.dump("ARC_dump.dat");
                ARC_control.dump("ARC_dump.par");
                c.print(std::cerr);
                abort();
            }
            else ds_use *= -dsf;
            // debuging
            // if (ptcl_org[0].id==267) {
            //		c.dump("ARC_dump.dat");
            //		ARC_control.dump("ARC_dump.par");
            //		c.print(std::cerr);
            //		abort();
            // }
        }
        else if (dsf==0) {
            //          char collerr[50]="two particle overlap!";
            c.info->ErrMessage(std::cerr);
            error_count++;
            if(error_count>4) {
                std::cerr<<"Error: Too much error appear!\nStep size ds: "<<ds_use<<"\nEnding physical time: "<<time_end<<"\nTime difference: "<<time_end-c.getTime()<<"\nR_in: "<<Int_pars.rin<<"\nR_out: "<<Int_pars.rout<<"\n";
                c.dump("ARC_dump.dat");
                ARC_control.dump("ARC_dump.par");
                c.print(std::cerr);
                abort();
            }
            if (c.info->status==5) {
                dscoff = 0.25;
                ds_use *= dscoff;
            }
            //          else if (c.info->status==6) ds_use *= 0.001;
            else if (c.info->status==4) ds_use = std::min(dscoff*c.calc_next_step_custom(),ds_up_limit);
            else ds_use *= 0.1;
            modify_step_flag=true;
        }
        else  {
            if (final_flag) {
                if (converge_count>6&&time_end-c.getTime()>ARC_control.dterr*100) {
                    std::cerr<<"Error: Time synchronization fails!\nStep size ds: "<<ds_use<<"\nEnding physical time: "<<time_end<<"\nTime difference: "<<time_end-c.getTime()<<"\nR_in: "<<Int_pars.rin<<"\nR_out: "<<Int_pars.rout<<"\n";
                    c.dump("ARC_dump.dat");
                    ARC_control.dump("ARC_dump.par");
                    c.print(std::cerr);
                    abort();
                }
                converge_count++;
            }
            else if (n_ptcl>2||(modify_step_flag&&error_count==0)) {
                ds_use = std::min(dscoff*c.calc_next_step_custom(),ds_up_limit);
                modify_step_flag=false;
            }
            // reducing error counter if integration success, this is to avoid the significant change of step may cause some issue
            if(error_count>0) error_count--;
        }
    }

    // update Rout
    if(n_ptcl>2) {
        std::size_t* list=new std::size_t[n_ptcl];
        Tptcl** plist=new Tptcl*[n_ptcl];
        c.getList(list);
#ifdef DEBUG_TEMP
        std::cerr<<"Curren r_out: n="<<n_ptcl;
        for (PS::S32 i=0; i<n_ptcl; i++) std::cerr<<std::setw(14)<<list[i]<<" ";
        std::cerr<<std::endl;
        for (PS::S32 i=0; i<n_ptcl; i++) std::cerr<<std::setw(14)<<ptcl_org[list[i]].r_out<<" ";
        std::cerr<<std::endl;
#endif
        for (PS::S32 i=0; i<n_ptcl; i++) plist[i] = &(ptcl_org[list[i]]);
        updateRout(plist,n_ptcl,Int_pars.rin,rout_single,gamma,m_average);
#ifdef DEBUG_TEMP
        std::cerr<<"new r_out: n="<<n_ptcl;
        for (PS::S32 i=0; i<n_ptcl; i++) std::cerr<<std::setw(14)<<list[i]<<" ";
        std::cerr<<std::endl;
        for (PS::S32 i=0; i<n_ptcl; i++) std::cerr<<std::setw(14)<<ptcl_org[list[i]].r_out<<" ";
        std::cerr<<std::endl;
#endif
#ifdef HARD_DEBUG
        for (PS::S32 i=0; i<n_ptcl; i++) 
            if (ptcl_org[list[i]].r_out<Int_pars.rin) {
                std::cerr<<"Error, updated p["<<list[i]<<"].rout ("<<ptcl_org[list[i]].r_out<<") < r_in ("<<Int_pars.rin<<")"<<std::endl;
                abort();
            }
#endif
        delete[] list;
        delete[] plist;
    }

    // error record
#ifdef ARC_ERROR
    PS::F64 ARC_error_temp = (c.getPot()+c.getEkin()-ARC_error_once);
    ARC_error += ARC_error_temp;
    ARC_error_relative += ARC_error_temp/ARC_error_once;
#endif
      
    // integration of center-of-mass
    c.cm.pos += c.cm.vel * time_end;

    c.center_shift_inverse();

    // if (!arout.is_open()) arout.open("arout");
    // for (int i=0;i<n_ptcl;i++) arout<<std::setprecision(17)<<time_origin_<<" "<<ptcl_org[i].mass<<" "<<ptcl_org[i].pos<<" "<<ptcl_org[i].vel<<std::endl;
}
#endif

#ifdef TIDAL_TENSOR
class TidalTensor{
private:
    PS::F64 T2[6];  // 1st (6)
    PS::F64 T3[10]; // 2nd Tensor (10)
public:

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

    template<class Tptcl>
    void fit(Tptcl* data, const PS::S32* list, const Binary& bin, const PS::S32 n_split) {
        PS::F64vec fi[8];
        if(data[list[10]].mass_bk>0) {
            for (PS::S32 i=0; i<8; i+=2) {
                fi[i]   = data[list[i+1]].vel;
                fi[i+1] = data[list[i  ]].vel;
            }
        }
        else {
            for (PS::S32 i=0; i<8; i++) 
                fi[i] = data[list[i]].vel;
        }
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
        PS::F64 T2S = 1.0/(bin.ax*(1+bin.ecc)*0.35);
        PS::F64 T3S = T2S*T2S;
        for (PS::S32 i=0; i<6;  i++) T2[i] *= T2S;
        for (PS::S32 i=0; i<10; i++) T3[i] *= T3S;
    }

    void eval(double* acc, const PS::F64vec &pos) const {
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
    
    template<class Tptcl> 
    void fit(Tptcl* data, const PS::S32* list, const Binary& bin, const PS::S32 n_split) {
#ifdef HARD_DEBUG
        assert(n_split>=4);
#endif
        const PS::F64 twopi = PI*2.0;
        peri_  = bin.peri;
        tperi_ = bin.tperi - twopi;  // make sure tperi_ is negative
        const PS::U32 np = n_split+1;
        PS::F64 x[np],y[3][np];
        for(PS::U32 i=0;i<np;i++) {
            x[i] = PS::F64(i)/n_split*twopi;
            x[i] = bin.peri/twopi * (x[i] - bin.ecc*std::sin(x[i]));
        }
        for(PS::S32 i=0; i<n_split; i++) {            
            for(PS::S32 j=0; j<3; j++)
                y[j][i] = data[list[2*i+1]].vel[j] - data[list[2*i]].vel[j];  // acc1 - acc0
        }
        for(PS::S32 i=0; i<3; i++) {
            y[i][n_split] = y[i][0];

            if(acc_[i]!=NULL) gsl_interp_accel_free(acc_[i]);
            acc_[i] = gsl_interp_accel_alloc();
            if(spline_[i]!=NULL) gsl_spline_free(spline_[i]);
            spline_[i] =  gsl_spline_alloc(t_, np);
            
            gsl_spline_init(spline_[i], x, y[i], np);
        }
    }

    void eval(double* acc, const PS::F64 time) const {
        PS::F64 dt = time - tperi_;
        dt = dt - (int)(dt/peri_)*peri_;
        for(int i=0; i<3; i++)
            acc[i] -= gsl_spline_eval(spline_[i],dt,acc_[i]);
    }

    ~keplerSplineFit() {
        for(int i=0;i<3;i++) {
            if(acc_[i]!=NULL) gsl_interp_accel_free(acc_[i]);
            if(spline_[i]!=NULL) gsl_spline_free(spline_[i]);
        }
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

template<class Tptcl, class Tpert, class Tpforce>
class ARCIntegrator{
private:
    typedef ARC::chain<Tptcl> ARChain;
    typedef ARC::chainpars ARControl;
    PS::ReallocatableArray<ARChain> clist_;
    PS::ReallocatableArray<ARC_pert_pars> par_list_;
    PS::ReallocatableArray<Tpert*> pert_;
    PS::ReallocatableArray<Tpforce*> pforce_;
    PS::ReallocatableArray<PS::S32> pert_n_;
    PS::ReallocatableArray<PS::S32> pert_disp_;

    ARControl *ARC_control_;
    ARC_int_pars *Int_pars_;

public:
    PS::ReallocatableArray<Binary> bininfo;
    //PS::ReallocatableArray<PS::F64> dt;

    ARCIntegrator() {};
    ARCIntegrator(ARControl &contr, ARC_int_pars &par): ARC_control_(&contr), Int_pars_(&par) {}

    void reserveARMem(const PS::S32 n) {
        clist_.reserve(n);
        //clist_.resizeNoInitialize(n);
        par_list_.reserve(n);
        //par_list_.resizeNoInitialize(n);
        pert_n_.reserve(n);
        //pert_n_.resizeNoInitialize(n);
        pert_disp_.reserve(n);
        //pert_disp_.resizeNoInitialize(n);
        bininfo.reserve(n);
        bininfo.resizeNoInitialize(n);
        //dt.reserve(n);
    }

    void reservePertMem(const PS::S32 n) {
        pert_.reserve(n);
        pforce_.reserve(n);
    }
    //void initialize(PS::S32 group_list[];
    //                ReallocatableArray<Tptcl> groups[],
    //                const PS::S32 n_groups,
    //                Tptcl* ptcl,
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

    void addOneGroup(Tptcl* ptcl_org,
                     const PS::S32* group_list,
                     const PS::S32 n_group,
                     const PS::S32* soft_pert_list,
                     const PS::S32 n_split,
                     Tpert* ptcl_pert = NULL,
                     Tpforce* pert_force = NULL,
                     const PS::S32* pert_list = NULL,
                     const PS::S32 n_pert = 0) {
        const PS::S32 igroup = clist_.size();
        const PS::S32 ioff = pert_.size();
        pert_disp_.push_back(ioff);
        pert_n_.push_back(0);
        clist_.increaseSize(1);
        clist_.back().allocate(n_group);
        for(int i=0; i<n_group; i++) {
            clist_.back().addP(ptcl_org[group_list[i]]);
        }

        // c.m.
        if(ptcl_pert!=NULL) {
            pert_.push_back(&ptcl_pert[igroup]);   // c.m. position is in igroup
            pforce_.push_back(&pert_force[igroup]);
            pert_n_[igroup]++;
        }

        // perturber
        for(int i=0; i<n_pert; i++) {
            const PS::S32  k = pert_list[i];
            //const PS::S32 ik = ioff+pert_n_[igroup]++;
            pert_.push_back(&ptcl_pert[k]);
            pforce_.push_back(&pert_force[k]);
            pert_n_[igroup]++;
        }
        par_list_.push_back(ARC_pert_pars(*Int_pars_));
        par_list_.back().fit(ptcl_org,soft_pert_list,bininfo[igroup],n_split);

        if(ptcl_pert!=NULL) {
            clist_.back().pos = ptcl_pert[igroup].pos;
            clist_.back().vel = ptcl_pert[igroup].vel;
            clist_.back().mass = ptcl_pert[igroup].mass;
        }

        // dt.push_back(0.0);

        //clist_.back().init(0.0, *ARC_control_, &(par_list_.back()));
        return;
    }

    void initialSlowDown(const PS::F64 tend, const PS::F64 sdfactor = 1.0e-8) {
        for (int i=0; i<clist_.size(); i++) {
            if (bininfo[i].ax>0&&bininfo[i].stable_factor>0) {
                PS::F64 finner = bininfo[i].ax*(1.0+bininfo[i].ecc);
                finner = clist_[i].mass/(finner*finner);
                finner = finner*finner;
                clist_[i].slowdown.setSlowDownPars(finner, bininfo[i].peri, sdfactor);
                Tptcl p[2];
                OrbParam2PosVel(p[0].pos, p[1].pos, p[0].vel, p[1].vel, bininfo[i].m1, bininfo[i].m2, bininfo[i].ax, bininfo[i].ecc, bininfo[i].inc, bininfo[i].OMG, bininfo[i].omg, PI);
                p[0].mass = bininfo[i].m1;
                p[1].mass = bininfo[i].m2;
#ifdef SOFT_PERT
#ifndef TIDAL_TENSOR
                p[0].status = 0;
                p[1].status = 1;
#endif
#endif
                //center_of_mass_correction(*(Tptcl*)&clist_[i], p, 2);
                PS::F64 acc[2][3];
                const PS::S32 ipert = pert_disp_[i];
                //Newtonian_extA(acc, bininfo[i].tperi+bininfo[i].peri, p, 2, &pert_[ipert], &pforce_[ipert], pert_n_[i], &par_list_[i]);
                if(pert_n_[i]>1) Newtonian_extA_pert(acc, 0.0, p, 2, &pert_[ipert], &pforce_[ipert], pert_n_[i], &par_list_[i]);
                else Newtonian_extA_soft(acc, 0.0, p, 2, &pert_[ipert], &pforce_[ipert], pert_n_[i], &par_list_[i]);
                PS::F64 fpertsq = 0.0;
                for(int k=0; k<3; k++) {
                    PS::F64 dacc = acc[0][k]-acc[1][k];
                    fpertsq += dacc*dacc;
                }
                clist_[i].slowdown.updatefpertsq(fpertsq);
                clist_[i].slowdown.updatekappa(0.0, tend);
            }
        }
    }

    void updateSlowDown(const PS::F64 tend) {
        for (int i=0; i<clist_.size(); i++) {
            clist_[i].slowdown.updatekappa(clist_[i].getTime(),tend);
        }
    }

    void adjustSlowDown(const PS::F64 dt) {
        for (int i=0; i<clist_.size(); i++) {
            clist_[i].slowdown.adjustkappa(dt);
        }
    }

    void adjustSlowDownPeriod(const PS::F64 dt, PS::S32* np) {
        for (int i=0; i<clist_.size(); i++) {
            np[i] = clist_[i].slowdown.adjustkappaPeriod(dt);
        }
    }
    
    void initial() {
        for (int i=0; i<clist_.size(); i++) {
            clist_[i].init(0.0, *ARC_control_, &(par_list_.back()));
        }
    }

    void dump(const char* fname, const PS::S32 ic, const PS::F64 time_end, const PS::F64 ds_use) {
        std::FILE* fp = std::fopen(fname,"w");
        if (fp==NULL) {
            std::cerr<<"Error: filename ARC_dump.dat cannot be open!\n";
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
        for (int i=0; i<np; i++) {
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


    PS::S64 integrateOneStepSym(const PS::S32 ic,
                                const PS::F64 time_end,
                                const PS::F64 dt_limit) {
        ARChain* c = &clist_[ic];
        ARC_pert_pars* par = &par_list_[ic];
        PS::F64 ds_up_limit = 0.25*dt_limit/c->calc_dt_X(1.0,*ARC_control_);
        //PS::F64 ds_use = 2.0*bininfo[ic].tstep*std::abs(c->getPt());
        PS::F64 ds_use=bininfo[ic].tstep;
        //PS::F64 ds_use = c->calc_next_step_custom(*ARC_control_,par);
        if (ds_use>ds_up_limit) ds_use = ds_up_limit;

        const PS::S32 ipert = pert_disp_[ic];
        bool fix_step_flag = false;
        if(c->getN()==2&&bininfo[ic].ecc>0.99) {
            fix_step_flag = true;
            PS::F64 korg=c->slowdown.getkappaorg();
            if(korg<1.0) ds_use *= std::max(0.1,korg);
        }

        PS::S64 stepcount = c->Symplectic_integration_tsyn(ds_use, *ARC_control_, time_end, par, &pert_[ipert], &pforce_[ipert], pert_n_[ic],fix_step_flag);

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
            std::cerr<<"ic = "<<ic<<" N = "<<c->getN()<<" Np = "<<pert_n_[ic]<<std::endl;
            abort();
        }
#endif
        
        return stepcount;
    }

    PS::S64 integrateOneStepExt(const PS::S32 ic,
                             const PS::F64 time_end,
                             const PS::F64 dt_limit) {
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
                    dump(ic,time_end,ds_use);
                    abort();
                }
                else ds_use *= -dsf;
            }
            else if (dsf==0) {
                c->info->ErrMessage(std::cerr);
                error_count++;
                if(error_count>4) {
                    std::cerr<<"Error: Too much error appear!\nStep size ds: "<<ds_use<<"\nEnding physical time: "<<time_end<<"\nTime difference: "<<time_end-c->getTime()<<"\nR_in: "<<Int_pars_->rin<<"\nR_out: "<<Int_pars_->rout<<"\n";
                    dump(ic,time_end,ds_use);
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
                        dump(ic,time_end,ds_use);
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

    PS::S64 integrateOneStepList(PS::S32 act_list[],
                              PS::S32 n_act,
                              const PS::F64 time_end,
                              const PS::F64 dt_limit) {
        PS::S64 nstep = 0;
        for(int i=0; i<n_act; i++) {
#ifdef ARC_SYM
            nstep += integrateOneStepSym(act_list[i], time_end, dt_limit);
#else
            nstep += integrateOneStepExt(act_list[i], time_end, dt_limit);
#endif
        }

        return nstep;
    }

    PS::S64 integrateOneStepList(const PS::F64 time_end,
                              const PS::F64 dt_limit) {
        PS::S64 nstep = 0;
        for(int i=0; i<clist_.size(); i++) {
#ifdef ARC_SYM
            nstep += integrateOneStepSym(i, time_end, dt_limit);
#else
            nstep += integrateOneStepExt(i, time_end, dt_limit);
#endif

        }
        return nstep;
    }

    template <class Tp>
    void updateCM(Tp ptcl[],
                  PS::S32 act_list[],
                  PS::S32 n_act) {
        for(int i=0; i<n_act; i++) {
            PS::S32 iact = act_list[i];
            clist_[iact].pos = ptcl[iact].pos;
            clist_[iact].vel = ptcl[iact].vel;
            clist_[iact].mass = ptcl[iact].mass;
        }
    }

    template <class Tp>
    void updateCM(Tp ptcl[]) {
        for(int i=0; i<clist_.size(); i++) {
            clist_[i].pos = ptcl[i].pos;
            clist_[i].vel = ptcl[i].vel;
#ifdef HARD_DEBUG
            assert(clist_[i].mass==ptcl[i].mass);
#endif
        }
    }

    void shift() {
        for(int i=0; i<clist_.size(); i++) {
            clist_[i].center_shift();
        }
    }

    void resolve() {
        for(int i=0; i<clist_.size(); i++) {
            clist_[i].resolve();
        }
    }

    PS::S32 getGroupN(const PS::S32 i) const {
        return clist_[i].getN();
    }

    PS::S32 getN() const {
        return clist_.size();
    }

    const Tptcl* getGroupPtcl(const PS::S32 i) const {
#ifdef ARC_ERROR
        assert(i<clist_.size());
#endif        
        return &clist_[i].getP(0);
    }

    PS::F64 getSlowDown(const PS::S32 i) const{
        return clist_[i].slowdown.getkappa();
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
        for (int j=0; j<clist_[i].getN(); j++) {
            os<<clist_[i].getP(j).mass<<" "
              <<clist_[i].getP(j).pos<<" "
              <<clist_[i].getP(j).vel<<std::endl;
        }
        for (int j=1; j<pert_n_[i]; j++) {
            os<<pert_[pert_disp_[i]+j]->mass<<" "
              <<pert_[pert_disp_[i]+j]->pos<<" "
              <<pert_[pert_disp_[i]+j]->vel<<std::endl;
        }
    }
    
    void info_print(std::ostream& os, const PS::S64 ngroups, const PS::S64 ngroupi, const PS::S64 nptcl, const PS::S64 ntot, const PS::F64 dt_limit, const PS::S32 kp) const{
        for (int i=0; i<clist_.size(); i++) {
            os<<"ARC_info(i,n_tot,n_ptcl,n_group,n,n_pert,ax,ecc,peri,kappa,n_step): "
              <<ngroups+i<<" "
              <<ntot<<" "
              <<nptcl<<" "
              <<ngroupi<<" "
              <<clist_[i].getN()<<" "
              <<pert_n_[i]<<" "
              <<bininfo[i].ax<<" "
              <<bininfo[i].ecc<<" "
              <<bininfo[i].peri<<" "
              <<clist_[i].slowdown.getkappa()<<" ";
            PS::S64 nstep = 0;
#ifdef ARC_SYM
            if(kp>0) nstep = kp*8;
            else nstep = clist_[i].profile.stepcount[0];
#else
            PS::S64 nstep = clist_[i].profile.itercount;
#endif
            os<<nstep<<std::endl;

            if(nstep>1e4) {
                os<<"Data dump for hardtest in the case of large nstep "<<nstep<<std::endl;
                data_dump(os, i, dt_limit);
            }
        }
    }
#endif
    
#ifdef ARC_ERROR
    template <class Teng>
    void EnergyRecord(Teng &energy, const PS::S32 sdflag=false) {
        energy.kin = energy.pot = energy.tot = 0.0;
        PS::F64 sd;
        for(int i=0; i<clist_.size(); i++) {
            if (sdflag) sd = 1.0/clist_[i].slowdown.getkappa();
            else sd = 1.0;
            energy.kin += sd*clist_[i].getEkin();
            energy.pot += sd*clist_[i].getPot();
            energy.tot += sd*clist_[i].getPt();
        }
    }
#endif                         

};
