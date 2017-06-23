#pragma once
#ifdef USE_INTRINSIC_FOR_X86
#include<immintrin.h>
#endif

#ifndef NAN_CHECK
#define NAN_CHECK(val) assert((val) == (val));
#endif

#include"integrate.hpp"
#include"cstdlib"
#include"ptcl.hpp"
//#include"stdio.h" /// for debug (L.Wang)

template<class T>
void Print(const T str, std::ostream & fout);

//std::ofstream kout;
//std::ofstream arout;

class PtclHard: public Ptcl{
public:
    PS::S32 id_cluster;
    PS::S32 adr_org;
    static PS::F64 search_factor;
    static PS::F64 r_search_min;

    PtclHard() {}

    template<class Tp>
    PtclHard(const Tp &p): Ptcl(p) {}

    template<class Tp>
    PtclHard(const Tp &p, const PS::S32 _id_cluster, const PS::S32 _adr_org): 
        Ptcl(p), id_cluster(_id_cluster), adr_org(_adr_org) {}

};

class ARC_int_pars{
public:
    static PS::F64 rout, rin, eps2;
    
    ARC_int_pars() {}
    ARC_int_pars(const ARC_int_pars& in_) {
        rout = in_.rout;
        rin  = in_.rin;
        eps2 = in_.eps2;
    }
};

class ARC_pert_pars: public ARC_int_pars, public keplerSplineFit{
    ARC_pert_pars(const ARC_int_pars& in_): ARC_int_pars(in_) {}
};

class SystemHard{
public:
    PS::ReallocatableArray<PtclHard> ptcl_hard_;
    ARC::chainpars ARC_control; ///chain controller (L.Wang)
#ifdef ARC_ERROR
    PS::F64 ARC_error_relative;
    PS::F64 ARC_error;  
    PS::S32 N_count[20];  // counting number of particles in one cluster
#endif
private:
    ARC_int_pars Int_pars; /// ARC integration parameters, rout_, rin_ (L.Wang)
    PS::F64 dt_limit_hard_;
    PS::F64 eta_s_;
    //PS::ReallocatableArray<PtclHard> ptcl_hard_;
    PS::ReallocatableArray<PS::S32> n_ptcl_in_cluster_;
    PS::ReallocatableArray<PS::S32> n_ptcl_in_cluster_disp_;
    PS::F64 time_origin_;
    PS::F64 gamma_;
    PS::F64 r_search_single_;
    PS::F64 m_average_;
    PS::S32 n_split_;

    ///////////
    /// functor
    //    struct OPSortClusterID{
    //        template<class T> bool operator() (const T & left, const T & right) const {
    //            return left.id_cluster < right.id_cluster;
    //        }
    //    };
    //    struct OPSortFirst{
    //        template<class T> bool operator() (const T & left, const T & right) const {
    //            return left.first < right.first;
    //        }
    //    };
    //    struct OPSortSecond{
    //        template<class T> bool operator() (const T & left, const T & right) const {
    //            return left.second < right.second;
    //        }
    //    };

    struct OPLessIDCluster{
        template<class T> bool operator() (const T & left, const T & right) const {
            return left.id_cluster < right.id_cluster;
        }
    };


    PS::F64 calcDtLimit(const PS::F64 time_sys,
                        const PS::F64 dt_limit_org,
                        const PS::F64 time_offset = 0.0){
        PS::F64 dt_limit_ret = dt_limit_org;
        PS::F64 s = (time_sys-time_offset) / dt_limit_ret;
        PS::F64 time_head = ((PS::S64)(s)) * dt_limit_org;
        PS::F64 time_shifted = time_sys - time_head;
        while( fmod(time_shifted, dt_limit_ret) != 0.0) dt_limit_ret *= 0.5;
        return dt_limit_ret;
    }

    template<class Tptcl>
    void driveForMultiClusterImpl(Tptcl * ptcl_org,
                                  const PS::S32 n_ptcl,
                                  const PS::F64 time_end,
                                  PS::ReallocatableArray<Tptcl> & ptcl_new) {
#ifdef HERMITE
        if(n_ptcl>5) {
            SearchGroup<Tptcl> group;
            
            group.findGroups(ptcl_org, n_ptcl);
            
            HermiteIntegrator<Tptcl> Hint;
            Hint.setParams(dt_limit_hard_, eta_s_, Int_pars.rin, Int_pars.rout, Int_pars.eps2);
            Hint.setPtcl(ptcl_org,n_ptcl,group.getPtclList(),group.getNPtcl());
            Hint.searchPerturber();

#ifdef HARD_DEBUG
            Energy E0, E1;
            Hint.CalcEnergyHard(E0);
#endif
            
            PS::S32 group_act_n = 0;
            PS::ReallocatableArray<PS::S32> group_act_list; //active group_list act adr
            // ReallocatableArray<PS::S32> group_list;     //group.adr list
            // ReallocatableArray<PS::S32> status;      //ptcl -> group.adr [non cm is -1] (value of Ptcl.status)
            // ReallocatableArray<PS::S32> status_map;  //ptcl -> group_list index [non cm is -1]
            // ReallocatableArray<PS::S32> adr_cm;         //group_list index -> ptcl.cm
            // group.findGroups(group_list, status, status_map,  adr_cm, group_act_n, ptcl_org, n_ptcl);
            group_act_list.resizeNoInitialize(group.getNPtcl());
            
            ARCIntegrator Aint(ARC_control, Int_pars);

            // first particles in Hint.Ptcl are c.m.
            PS::S32 n_groups = group.getNGroups()
            for (int i=0; i<n_groups; i++) 
                Aint.addOneGroup(ptcl_org, group.getGroup(i), group.getGroupN(i), group.getGroupPertList(i,n_split_), Hint.getPtcl(), Hint.getForce(), Hint.getPertList(i), Hint.getPertN(i)); 
            
            PS::S32 time_sys=0.0;
            PS::F64 dt_limit = calcDtLimit(time_sys, dt_limit_hard_);
            Hint.initialize(dt_limit, group_act_list.getPointer(), group_act_n, n_groups);
            while(time_sys<time_end) {
                time_sys = Hint.getNextTime();
                dt_limit = calcDtLimit(time_sys, dt_limit_hard_);
                Aint.integrateOneStepList(group_act_list.getPointer(), group_act_n, time_sys, dt_limit);
                Hint.integrateOneStep(time_sys,dt_limit);
                Aint.updateCM(Hint.getPtcl(), group_act_list.getPointer(), group_act_n);
                Hint.SortAndSelectIp(group_act_list.getPointer(), group_act_n, n_groups);
            }
            Aint.resolve();
            Hint.writeBackPtcl(ptcl_org,n_ptcl,group.getPtclList(),group.getNPtcl());
            
            //group.resolveGroups(ptcl_org, n_ptcl, group_ptcl_glb.getPointer(), group_list.size(), group_list.getPointer(), adr_cm.getPointer());
            group.resolveGroups();
            group.searchAndMerge(ptcl_org, n_ptcl, Int_pars.rin);
            // Kickcorrect(ptcl_org, group.getRoutChangeList());
            group.generatelist(ptcl_org, n_ptcl, ptcl_new);

            // group.reverseCopy(ptcl_org, n_ptcl);
        }
        else 
#endif
            PS::F64 dt_limit = calcDtLimit(0.0, dt_limit_hard_);
            Multiple_integrator(ptcl_org, n_ptcl, time_end, dt_limit,
                                r_search_single_, gamma_, m_average_,
#ifdef ARC_ERROR
                                ARC_error_relative,
                                ARC_error,
                                N_count,
#endif
                                ARC_control, Int_pars);
    }

public:

    SystemHard(){
#ifdef ARC_ERROR
        ARC_error = 0.0;
        ARC_error_relative = 0.0;
        for(PS::S32 i=0;i<20;i++) N_count[i]=0;
#endif
        //        PS::S32 n_threads = PS::Comm::getNumberOfThread();
    }

    /// start set Chainpars (L.Wang)
    ///
    void setARCParam(const PS::F64 energy_error=1e-10, const PS::F64 dterr=1e-9, const PS::F64 dtmin=1e-24, const PS::S32 exp_method=1, const PS::S32 exp_itermax=20, const PS::S32 den_intpmax=20, const PS::S32 exp_fix_iter=0) {
        ARC_control.setA(Newtonian_cut_AW<PtclHard,ARC_pert_pars>,Newtonian_extA<PtclHard,PtclHard*,PtclForce*,ARC_pert_pars>,Newtonian_timescale<ARC_pert_pars>);
        ARC_control.setabg(0,1,0);
        ARC_control.setErr(energy_error,dtmin,dterr);
        ARC_control.setIterSeq(exp_itermax,3,den_intpmax);
        ARC_control.setIntp(exp_method);
        ARC_control.setIterConst((bool)exp_fix_iter);
        ARC_control.setAutoStep(3);
    }
    /// end set Chainpars (L.Wang)

    void initializeForOneCluster(const PS::S32 n){
        ptcl_hard_.resizeNoInitialize(n);
    }

    ////////////////////////
    // for NON-ISOLATED CLUSTER
    template<class Tsys, class Tptcl, class Tmediator>
    void setPtclForConectedCluster(const Tsys & sys,
                                   const PS::ReallocatableArray<Tmediator> & med,
                                   const PS::ReallocatableArray<Tptcl> & ptcl_recv){
        ptcl_hard_.clearSize();
        n_ptcl_in_cluster_.clearSize(); // clear befor break this function
        for(PS::S32 i=0; i<med.size(); i++){
            if(med[i].adr_sys_ < 0) continue;
            if(med[i].rank_send_ != PS::Comm::getRank()) continue;
            const auto & p = sys[med[i].adr_sys_];
            ptcl_hard_.push_back(PtclHard(p, med[i].id_cluster_, med[i].adr_sys_));
        }

        for(PS::S32 i=0; i<ptcl_recv.size(); i++){
            const Tptcl & p = ptcl_recv[i];
            ptcl_hard_.push_back(PtclHard(p, p.id_cluster_, -(i+1)));
        }

        if(ptcl_hard_.size() == 0) return;
        std::sort(ptcl_hard_.getPointer(), ptcl_hard_.getPointer(ptcl_hard_.size()), 
                  OPLessIDCluster());
        PS::S32 n_tot = ptcl_hard_.size();
        PS::S32 id_cluster_ref = -999;
        for(PS::S32 i=0; i<n_tot; i++){
            if(id_cluster_ref != ptcl_hard_[i].id_cluster){
                id_cluster_ref = ptcl_hard_[i].id_cluster;
                n_ptcl_in_cluster_.push_back(0);
            }
            n_ptcl_in_cluster_.back()++;
        }
        PS::S32 n_cluster = n_ptcl_in_cluster_.size();
        n_ptcl_in_cluster_disp_.resizeNoInitialize(n_cluster+1);
        n_ptcl_in_cluster_disp_[0] = 0;
        for(PS::S32 i=0; i<n_cluster; i++){
            n_ptcl_in_cluster_disp_[i+1] = n_ptcl_in_cluster_disp_[i] + n_ptcl_in_cluster_[i];
        }
    }


    // for NON-ISOLATED CLUSTER
    ////////////////////////


    const PS::ReallocatableArray<PtclHard> & getPtcl() const {
        return ptcl_hard_;
    }

    void setTimeOrigin(const PS::F64 _time_origin){
        time_origin_ = _time_origin;
    }

    void setParam(const PS::F64 _rsearch, 
                  const PS::F64 _rin,
                  const PS::F64 _eps,
                  const PS::F64 _dt_limit_hard,
                  const PS::F64 _eta,
                  const PS::F64 _time_origin,
                  const PS::F64 _gmin,
                  const PS::F64 _m_avarage,
                  const PS::S32 _n_split = 8){
        /// Set chain pars (L.Wang)
		Int_pars.rin  = _rin;
        Int_pars.eps2  = _eps*_eps;
        /// Set chain pars (L.Wang)        
        dt_limit_hard_ = _dt_limit_hard;
        eta_s_ = _eta*_eta;
        time_origin_ = _time_origin;
        gamma_ = std::pow(1.0/_gmin,0.33333);
        r_search_single_ = _rsearch; 
        m_average_ = _m_avarage;
        n_split_ = _n_split;
    }


//////////////////
// for one cluster
    template<class Tsys>
    void setPtclForOneCluster(const Tsys & sys, 
                              const PS::ReallocatableArray<PS::S32> & adr_array){
        // for one cluster
        const PS::S32 n = adr_array.size();
        //ptcl_hard_.resizeNoInitialize(n);
        //n_ptcl_in_cluster_.resizeNoInitialize(n);
        for(PS::S32 i=0; i<n; i++){
            PS::S32 adr = adr_array[i];
            ptcl_hard_[i].DataCopy(sys[adr]);
            //n_ptcl_in_cluster_[i] = 1;
        }
    }

    template<class Tsys>
    void setPtclForOneClusterOMP(const Tsys & sys, 
                                 const PS::ReallocatableArray<PS::S32> & adr_array){
        // for one cluster
        const PS::S32 n = adr_array.size();
        //ptcl_hard_.resizeNoInitialize(n);
        //n_ptcl_in_cluster_.resizeNoInitialize(n);
#pragma omp for schedule(dynamic)
        for(PS::S32 i=0; i<n; i++){
            PS::S32 adr = adr_array[i];
            ptcl_hard_[i].DataCopy(sys[adr]);
            //n_ptcl_in_cluster_[i] = 1;
        }
    }

    void driveForOneCluster(const PS::F64 dt){
        const PS::S32 n = ptcl_hard_.size();
        for(PS::S32 i=0; i<n; i++){
            PS::F64vec dr = ptcl_hard_[i].vel * dt;
            ptcl_hard_[i].pos += dr;
            ptcl_hard_[i].r_search = std::max(std::sqrt(dr*dr)*PtclHard::search_factor, PtclHard::r_search_min);
            // ptcl_hard_[i].r_search= r_search_single_;
            /*
              DriveKeplerRestricted(mass_sun_, 
              pos_sun_, ptcl_hard_[i].pos, 
              vel_sun_, ptcl_hard_[i].vel, dt); 
            */
        }

    }
    void driveForOneClusterOMP(const PS::F64 dt){
        const PS::S32 n = ptcl_hard_.size();
#pragma omp for schedule(dynamic)
        for(PS::S32 i=0; i<n; i++){
            PS::F64vec dr = ptcl_hard_[i].vel * dt;
            ptcl_hard_[i].pos += dr;
            ptcl_hard_[i].r_search = std::max(std::sqrt(dr*dr)*PtclHard::search_factor, PtclHard::r_search_min);
            /*
              DriveKeplerRestricted(mass_sun_, 
              pos_sun_, ptcl_hard_[i].pos, 
              vel_sun_, ptcl_hard_[i].vel, dt); 
            */
        }
    }

    template<class Tsys>
    void writeBackPtclForOneCluster(Tsys & sys, 
                                    const PS::ReallocatableArray<PS::S32> & adr_array){
        const PS::S32 n = ptcl_hard_.size();
        for(PS::S32 i=0; i<n; i++){
            PS::S32 adr = adr_array[i];
            assert(sys[adr].id == ptcl_hard_[i].id);
            sys[adr].DataCopy(ptcl_hard_[i]);
        }
    }

    template<class Tsys>
    void writeBackPtclForOneClusterOMP(Tsys & sys, 
                                       const PS::ReallocatableArray<PS::S32> & adr_array){
        const PS::S32 n = ptcl_hard_.size();
#pragma omp for schedule(dynamic)
        for(PS::S32 i=0; i<n; i++){
            PS::S32 adr = adr_array[i];
            assert(sys[adr].id == ptcl_hard_[i].id);
            sys[adr].DataCopy(ptcl_hard_[i]);
        }
    }
// for one cluster
//////////////////


//////////////////
// for isolated multi cluster only
    template<class Tsys>
    void setPtclForIsolatedMultiCluster(const Tsys & sys,
                                        const PS::ReallocatableArray<PS::S32> & _adr_array,
                                        const PS::ReallocatableArray<PS::S32> & _n_ptcl_in_cluster){
        const PS::S32 n_ptcl = _adr_array.size();
        ptcl_hard_.resizeNoInitialize(n_ptcl);
        for(PS::S32 i=0; i<n_ptcl; i++){
            PS::S32 adr = _adr_array[i];
            ptcl_hard_[i].DataCopy(sys[adr]);
            //  ptcl_hard_[i].n_ngb= sys[adr].n_ngb;
        }
        const PS::S32 n_cluster = _n_ptcl_in_cluster.size();
        n_ptcl_in_cluster_.resizeNoInitialize(n_cluster);
        n_ptcl_in_cluster_disp_.resizeNoInitialize(n_cluster+1);
        n_ptcl_in_cluster_disp_[0] = 0;
        for(PS::S32 i=0; i<n_cluster; i++){
            n_ptcl_in_cluster_[i] = _n_ptcl_in_cluster[i];
            n_ptcl_in_cluster_disp_[i+1] = n_ptcl_in_cluster_disp_[i] + n_ptcl_in_cluster_[i];
        }
    }

    void initailizeForIsolatedMultiCluster(const PS::S32 _n_ptcl,
                                           const PS::ReallocatableArray<PS::S32> & _n_ptcl_in_cluster){
        ptcl_hard_.resizeNoInitialize(_n_ptcl);
        const PS::S32 n_cluster = _n_ptcl_in_cluster.size();
        n_ptcl_in_cluster_.resizeNoInitialize(n_cluster);
        n_ptcl_in_cluster_disp_.resizeNoInitialize(n_cluster+1);
        n_ptcl_in_cluster_disp_[0] = 0;
        for(PS::S32 i=0; i<n_cluster; i++){
            n_ptcl_in_cluster_[i] = _n_ptcl_in_cluster[i];
            n_ptcl_in_cluster_disp_[i+1] = n_ptcl_in_cluster_disp_[i] + n_ptcl_in_cluster_[i];
        }
    }

    template<class Tsys>
    void setPtclForIsolatedMultiClusterOMP(const Tsys & sys,
                                           const PS::ReallocatableArray<PS::S32> & _adr_array,
                                           const PS::ReallocatableArray<PS::S32> & _n_ptcl_in_cluster){
        const PS::S32 n_ptcl = _adr_array.size();
#pragma omp for schedule(dynamic)
        for(PS::S32 i=0; i<n_ptcl; i++){
            PS::S32 adr = _adr_array[i];
            ptcl_hard_[i].DataCopy(sys[adr]);
            //  ptcl_hard_[i].n_ngb = sys[adr].n_ngb;
        }
    }

    template<class Tsys>
    void writeBackPtclForMultiCluster(Tsys & sys, 
                                      const PS::ReallocatableArray<PS::S32> & adr_array){
        writeBackPtclForOneCluster(sys, adr_array);
    }
    template<class Tsys>
    void writeBackPtclForMultiClusterOMP(Tsys & sys, 
                                         const PS::ReallocatableArray<PS::S32> & adr_array){
        writeBackPtclForOneClusterOMP(sys, adr_array);
    }
// for isolated multi cluster only
//////////////////

//////////////////
// for multi cluster
    template<class Tptcl>
    void driveForMultiCluster(const PS::F64 dt){
        const PS::S32 n_cluster = n_ptcl_in_cluster_.size();
        /*
          for(PS::S32 ith=0; ith<PS::Comm::getNumberOfThread(); ith++){
          eng_disp_merge_omp_[ith] = 0.0;
          merge_log_omp_[ith].clearSize();
          }
        */
        for(PS::S32 i=0; i<n_cluster; i++){
            const PS::S32 adr_head = n_ptcl_in_cluster_disp_[i];
            const PS::S32 n_ptcl = n_ptcl_in_cluster_[i];
            PS::ReallocatableArray<Tptcl> extra_ptcl;
            driveForMultiClusterImpl(ptcl_hard_.getPointer(adr_head), n_ptcl, dt, extra_ptcl);
            {
                
            }
        }
    }

    template<class Tptcl>
    void driveForMultiClusterOMP(const PS::F64 dt){

        const PS::S32 n_cluster = n_ptcl_in_cluster_.size();
        //	const PS::S32 ith = PS::Comm::getThreadNum();
#pragma omp for schedule(dynamic)
        for(PS::S32 i=0; i<n_cluster; i++){
            const PS::S32 adr_head = n_ptcl_in_cluster_disp_[i];
            const PS::S32 n_ptcl = n_ptcl_in_cluster_[i];
            PS::ReallocatableArray<Tptcl> extra_ptcl;
            driveForMultiClusterImpl(ptcl_hard_.getPointer(adr_head), n_ptcl, dt, extra_ptcl);
#pragma omp critical
            {
                
            }
            
        }
    }


};

