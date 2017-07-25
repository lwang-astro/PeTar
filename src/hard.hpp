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
#include"cluster_list.hpp"
//#include"stdio.h" /// for debug (L.Wang)

//template<class T>
//void Print(const T str, std::ostream & fout);

//std::ofstream kout;
//std::ofstream arout;

class PtclHard: public Ptcl{
public:
    PS::S32 id_cluster;
    PS::S32 adr_org;

    PtclHard() {}

    template<class Tp>
    PtclHard(const Tp &p): Ptcl(p) {}

    template<class Tp>
    PtclHard(const Tp &p, const PS::S32 _id_cluster, const PS::S32 _adr_org): 
        Ptcl(p), id_cluster(_id_cluster), adr_org(_adr_org) {}

    template<class Tp>
    void DataCopy(const Tp& p) {
        Ptcl::DataCopy(p);
    }

    PtclHard& operator = (const PtclHard& p) {
        Ptcl::DataCopy(p);
        id_cluster = p.id_cluster;
        adr_org = p.adr_org;
        return *this;
    }

};

class ARC_int_pars{
public:
    PS::F64 rout, rin, eps2;
    
    ARC_int_pars() {}
    ARC_int_pars(const ARC_int_pars& in_) {
        rout = in_.rout;
        rin  = in_.rin;
        eps2 = in_.eps2;
    }
};

class ARC_pert_pars: public ARC_int_pars, public keplerSplineFit{
public:
    ARC_pert_pars() {}
    ARC_pert_pars(const ARC_int_pars& in_): ARC_int_pars(in_) {}
};

#ifdef HARD_DEBUG
class HardEnergy {
public:
    PS::F64 kin, pot, tot;
};
#endif

class SystemHard{
public:
    PS::ReallocatableArray<PtclHard> ptcl_hard_;
    ARC::chainpars ARC_control_; ///chain controller (L.Wang)
#ifdef ARC_ERROR
    PS::F64 ARC_error_relative;
    PS::F64 ARC_error;  
    PS::S32 N_count[20];  // counting number of particles in one cluster
#endif
#ifdef HARD_DEBUG
    HardEnergy E0, E1;
    HardEnergy AE0, AE1;
    HardEnergy HE0, HE1;
    HardEnergy ESD0, ESD1;
#endif

private:
    ARC_int_pars Int_pars_; /// ARC integration parameters, rout_, rin_ (L.Wang)
    PS::F64 dt_limit_hard_;
    PS::F64 eta_s_;
    //PS::ReallocatableArray<PtclHard> ptcl_hard_;
    PS::ReallocatableArray<PS::S32> n_ptcl_in_cluster_;
    PS::ReallocatableArray<PS::S32> n_ptcl_in_cluster_disp_;
    PS::F64 time_origin_;
    // PS::F64 gamma_;
    // PS::F64 r_search_single_;
    PS::F64 r_bin_;
    // PS::F64 m_average_;
    PS::S64 id_offset_;
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

    void driveForMultiClusterImpl(PtclHard * ptcl_org,
                                  const PS::S32 n_ptcl,
                                  const PS::F64 time_end,
                                  PS::ReallocatableArray<PtclHard> & ptcl_new) {
//#ifdef HERMITE
//        if(n_ptcl>5) {
        SearchGroup<PtclHard> group;
            
        group.findGroups(ptcl_org, n_ptcl, n_split_);

        if(group.getPtclN()==1) {
            PtclHard* pcm = &ptcl_org[group.getPtclList()[0]];
            PS::S32 iact = 0;
            
            ARCIntegrator<PtclHard, PtclH4, PtclForce, ARC_int_pars, ARC_pert_pars> Aint(ARC_control_, Int_pars_);
            Aint.reserveARMem(1);
            // Aint.reservePertMem(1);
            group.getBinPars(Aint.bininfo[0],ptcl_org,0,n_split_);
            Aint.addOneGroup(ptcl_org, group.getGroup(0), group.getGroupN(0),group.getGroupPertList(0,n_split_), n_split_);
            Aint.updateCM(pcm, &iact, 1, true);

            Aint.initialSlowDown(time_end);
            Aint.initial();

            Aint.adjustSlowDown(time_end);

            Aint.integrateOneStep(0, time_end, dt_limit_hard_);
            
            pcm->pos += pcm->vel * time_end;

            Aint.updateCM(pcm, &iact, 1);
            Aint.resolve();
#ifdef HARD_DEBUG
            fprintf(stderr,"Slowdown factor = %e\n", Aint.getSlowDown(0));
#endif

        }
        else {
            
            HermiteIntegrator<PtclHard> Hint;
            Hint.setParams(dt_limit_hard_, eta_s_, Int_pars_.rin, Int_pars_.rout, Int_pars_.eps2);
            Hint.setPtcl(ptcl_org,n_ptcl,group.getPtclList(),group.getPtclN());
            Hint.searchPerturber();

            PS::F64 time_sys=0.0, time_now;
#ifdef FIX_STEP_DEBUG
            PS::F64 dt_limit = dt_limit_hard_;
#else
            PS::F64 dt_limit = calcDtLimit(time_sys, dt_limit_hard_);
#endif
            
            PS::S32 group_act_n = 0;
            PS::ReallocatableArray<PS::S32> group_act_list; //active group_list act adr
            // ReallocatableArray<PS::S32> group_list;     //group.adr list
            // ReallocatableArray<PS::S32> status;      //ptcl -> group.adr [non cm is -1] (value of Ptcl.status)
            // ReallocatableArray<PS::S32> status_map;  //ptcl -> group_list index [non cm is -1]
            // ReallocatableArray<PS::S32> adr_cm;         //group_list index -> ptcl.cm
            // group.findGroups(group_list, status, status_map,  adr_cm, group_act_n, ptcl_org, n_ptcl);

            group_act_list.resizeNoInitialize(group.getPtclN());
            
            ARCIntegrator<PtclHard, PtclH4, PtclForce, ARC_int_pars, ARC_pert_pars> Aint(ARC_control_, Int_pars_);

            // first particles in Hint.Ptcl are c.m.
            PS::S32 n_groups = group.getNumOfGroups();
            Aint.reserveARMem(n_groups);
            Aint.reservePertMem(Hint.getPertListSize());
            for (int i=0; i<n_groups; i++) {
                group.getBinPars(Aint.bininfo[i],ptcl_org,i,n_split_);
                Aint.addOneGroup(ptcl_org, group.getGroup(i), group.getGroupN(i), group.getGroupPertList(i,n_split_), n_split_, Hint.getPtcl(), Hint.getForce(), Hint.getPertList(i), Hint.getPertN(i)); 
            }
            Aint.initialSlowDown(dt_limit);
            Aint.initial();


            Hint.initialize(dt_limit, group_act_list.getPointer(), group_act_n, n_groups, &Aint);


#ifdef HARD_DEBUG
            CalcEnergyHardFull(E0, AE0, HE0, ESD0, Hint, Aint, group);
            PS::ReallocatableArray<PS::F64> slowdownrecord;
            slowdownrecord.resizeNoInitialize(n_groups);
#endif

            while(time_sys<time_end) {
                time_now = time_sys;
                time_sys = Hint.getNextTime();
#ifdef FIX_STEP_DEBUG
                dt_limit = dt_limit_hard_;
#else
                dt_limit = calcDtLimit(time_sys, dt_limit_hard_);
#endif

#ifdef HARD_DEBUG
                assert(time_sys>time_now);
#endif
                PS::F64 dt_h = time_sys-time_now;
                //Aint.updateSlowDown(time_sys);
#ifdef HARD_DEBUG
                for(int k=0; k<n_groups; k++) {
                    slowdownrecord[k] = std::max(slowdownrecord[k], Aint.getSlowDown(k));
                    assert(Aint.getSlowDown(k)>=1.0);
                }
#endif
                //Aint.integrateOneStepList(group_act_list.getPointer(), group_act_n, time_sys, dt_limit);
                Aint.integrateOneStepList(time_sys, std::min(dt_limit,dt_h));
                Hint.integrateOneStep(time_sys,dt_limit,true,&Aint);
                //Hint.SortAndSelectIp(group_act_list.getPointer(), group_act_n, n_groups);
                Hint.SortAndSelectIp();
            }
        
            Aint.updateCM(Hint.getPtcl());
            Aint.resolve();
            Hint.writeBackPtcl(ptcl_org,n_ptcl,group.getPtclList(),group.getPtclN());

#ifdef HARD_DEBUG
            CalcEnergyHardFull(E1, AE1, HE1, ESD1, Hint, Aint, group);
        
#ifdef HARD_DEBUG_PRINT
            fprintf(stderr,"Slowdown factor = ");
            for(int k=0; k<n_groups; k++) 
                fprintf(stderr,"%e; ",slowdownrecord[k]);
            fprintf(stderr,"\n");
            fprintf(stderr,"H4  Energy: init =%e, end =%e, diff =%e, kin =%e pot =%e\nARC Energy: init =%e, end =%e, diff =%e, error = %e\nTot Energy: init =%e, end =%e, diff =%e, kin =%e pot =%e, Tot-H4-ARC =%e\nTSD Energy: init =%e, end =%e, diff =%e, kin =%e pot =%e\n", 
                    HE0.tot, HE1.tot, HE1.tot-HE0.tot, HE1.kin, HE1.pot, 
                    AE0.kin+AE0.pot, AE1.kin+AE1.pot, AE1.kin+AE1.pot-AE0.kin-AE0.pot, (AE1.kin+AE1.pot+AE1.tot-AE0.kin-AE0.pot-AE0.tot)/AE0.tot,
                    E0.tot, E1.tot, E1.tot-E0.tot, E1.kin, E1.pot, E1.tot-HE1.tot-AE1.kin-AE1.pot,
                    ESD0.tot, ESD1.tot, ESD1.tot-ESD0.tot, ESD1.kin, ESD1.pot);
            Hint.printStepHist();
#endif
#endif
        }
            
        //group.resolveGroups(ptcl_org, n_ptcl, group_ptcl_glb.getPointer(), group_list.size(), group_list.getPointer(), adr_cm.getPointer());
        group.resolveGroups();
        updateRSearch(ptcl_org, group.getPtclList(), group.getPtclN(), time_end);

        group.searchAndMerge(ptcl_org, r_bin_);
        // Kickcorrect(ptcl_org, group.getRoutChangeList());
        group.generateList(ptcl_org, ptcl_new, r_bin_, time_end, id_offset_, n_split_);

            // group.reverseCopy(ptcl_org, n_ptcl);
//        }
//        else {
//#endif
//            PS::F64 dt_limit = calcDtLimit(0.0, dt_limit_hard_);
//            Multiple_integrator(ptcl_org, n_ptcl, time_end, dt_limit,
//                                r_search_single_, gamma_, m_average_,
//#ifdef ARC_ERROR
//                                ARC_error_relative,
//                                ARC_error,
//                                N_count,
//#endif
//                                ARC_control_, Int_pars_);
//        }
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
    void setARCParam(const PS::F64 energy_error=1e-10, const PS::F64 dterr=1e-6, const PS::F64 dtmin=1e-24, const PS::S32 exp_method=1, const PS::S32 exp_itermax=20, const PS::S32 den_intpmax=20, const PS::S32 exp_fix_iter=0) {
#ifdef HARD_DEBUG_DEEP_CHECK
        ARC_control_.setA(Newtonian_cut_AW<PtclHard,ARC_pert_pars>,Newtonian_extA_test<PtclHard,PtclH4*,PtclForce*,ARC_pert_pars>,Newtonian_timescale<ARC_pert_pars>);
#else
        ARC_control_.setA(Newtonian_cut_AW<PtclHard,ARC_pert_pars>,Newtonian_extA<PtclHard,PtclH4*,PtclForce*,ARC_pert_pars>,Newtonian_timescale<ARC_pert_pars>);
#endif
        ARC_control_.setabg(0,1,0);
        ARC_control_.setErr(energy_error,dtmin,dterr);
        ARC_control_.setIterSeq(exp_itermax,3,den_intpmax);
        ARC_control_.setIntp(exp_method);
        ARC_control_.setIterConst((bool)exp_fix_iter);
        ARC_control_.setAutoStep(3);
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
            ptcl_hard_.push_back(PtclHard(p, p.id_cluster, -(i+1)));
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

    void setParam(const PS::F64 _rbin,
                  const PS::F64 _rout,
                  const PS::F64 _rin,
                  const PS::F64 _eps,
                  const PS::F64 _dt_limit_hard,
                  const PS::F64 _eta,
                  const PS::F64 _time_origin,
                  // const PS::F64 _gmin,
                  // const PS::F64 _m_avarage,
                  const PS::S64 _id_offset,
                  const PS::S32 _n_split = 8){
        /// Set chain pars (L.Wang)
		Int_pars_.rin  = _rin;
        Int_pars_.rout = _rout;
        Int_pars_.eps2  = _eps*_eps;
        /// Set chain pars (L.Wang)        
        dt_limit_hard_ = _dt_limit_hard;
        eta_s_ = _eta*_eta;
        time_origin_ = _time_origin;
//        gamma_ = std::pow(1.0/_gmin,0.33333);
        // r_search_single_ = _rsearch; 
        r_bin_           = _rbin;
        // m_average_ = _m_avarage;
        n_split_ = _n_split;
        id_offset_ = _id_offset;
    }

    void updateRSearch(PtclHard* ptcl_org,
                       const PS::S32* ptcl_list,
                       const PS::S32 n_ptcl,
                       const PS::F64 dt_tree) {
        for (PS::S32 i=0; i<n_ptcl; i++) {
            ptcl_org[ptcl_list[i]].calcRSearch(dt_tree);
        }
    }

#ifdef HARD_DEBUG
    template<class Teng>
    void CalcEnergyHard(Teng & eng,  const PS::S32* ptcl_list, const PS::S32 ptcl_n, const PS::S32* group_list, const PS::S32 group_n, const PS::S32 nbin){
        eng.kin = eng.pot = eng.tot = 0.0;
        
        for(PS::S32 i=nbin; i<ptcl_n; i++){
            PtclHard* pi = &ptcl_hard_[ptcl_list[i]];
            eng.kin += 0.5 * pi->mass * pi->vel * pi->vel;

            for(PS::S32 j=i+1; j<ptcl_n; j++){
                PtclHard* pj = &ptcl_hard_[ptcl_list[j]];
                PS::F64vec rij = pi->pos - pj->pos;
                PS::F64 dr = sqrt(rij*rij + Int_pars_.eps2);
                eng.pot -= pj->mass*pi->mass/dr*(1.0 - CalcW(dr/Int_pars_.rout, Int_pars_.rin/Int_pars_.rout));
            }

            for(PS::S32 j=0; j<group_n; j++){
                PtclHard* pj = &ptcl_hard_[group_list[j]];
                PS::F64vec rij = pi->pos - pj->pos;
                PS::F64 dr = sqrt(rij*rij + Int_pars_.eps2);
                eng.pot -= pj->mass*pi->mass/dr*(1.0 - CalcW(dr/Int_pars_.rout, Int_pars_.rin/Int_pars_.rout));
            }
        }

        for(PS::S32 i=0; i<group_n; i++){
            PtclHard* pi = &ptcl_hard_[group_list[i]];
            eng.kin += 0.5 * pi->mass * pi->vel * pi->vel;

            for(PS::S32 j=i+1; j<group_n; j++){
                PtclHard* pj = &ptcl_hard_[group_list[j]];
                PS::F64vec rij = pi->pos - pj->pos;
                PS::F64 dr = sqrt(rij*rij + Int_pars_.eps2);
                eng.pot -= pj->mass*pi->mass/dr*(1.0 - CalcW(dr/Int_pars_.rout, Int_pars_.rin/Int_pars_.rout));
            }
        }
        eng.tot = eng.kin + eng.pot;
    }

    template<class Teng, class TH4, class TARC, class Tgroup>
    void CalcEnergyHardFull(Teng& E, Teng& AE, Teng& HE, Teng& ESD, TH4 &Hint, TARC& Aint, const Tgroup& group){
        Hint.CalcEnergy(HE);
        Teng TMP;
        Aint.EnergyRecord(TMP,true);
        Aint.EnergyRecord(AE);
        CalcEnergyHard(E, group.getPtclList(), group.getPtclN(), group.getGroup(0), group.getGroupListSize(), group.getNumOfGroups());
        ESD.tot = (E.tot - AE.kin-AE.pot) + (TMP.kin+TMP.pot);
        ESD.kin = (E.kin - AE.kin) + TMP.kin;
        ESD.pot = (E.pot - AE.pot) + TMP.pot;
    }
#endif

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
            ptcl_hard_[i].calcRSearch(dt);
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
            ptcl_hard_[i].calcRSearch(dt);
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
            // assert(sys[adr].id == ptcl_hard_[i].id);
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
            // assert(sys[adr].id == ptcl_hard_[i].id);
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
    template<class Tsys, class Tsptcl>
    void driveForMultiCluster(const PS::F64 dt, Tsys & sys){
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
            PS::ReallocatableArray<PtclHard> extra_ptcl;
            driveForMultiClusterImpl(ptcl_hard_.getPointer(adr_head), n_ptcl, dt, extra_ptcl);
#ifdef HARD_DEBUG
            if(extra_ptcl.size()>0) fprintf(stderr,"New particle number = %d\n",extra_ptcl.size());
#endif
            for (PS::S32 j=0; j<extra_ptcl.size(); j++) {
                PS::S32 adr = sys.getNumberOfParticleLocal();
                PS::S32 rank = PS::Comm::getRank();
                sys.addOneParticle(Tsptcl(extra_ptcl[j],rank,adr));
            }
        }
    }

    template<class Tsys, class Tsptcl>
    void driveForMultiClusterOMP(const PS::F64 dt, Tsys & sys){
        const PS::S32 n_cluster = n_ptcl_in_cluster_.size();
        //	const PS::S32 ith = PS::Comm::getThreadNum();
#pragma omp for schedule(dynamic)
        for(PS::S32 i=0; i<n_cluster; i++){
            const PS::S32 adr_head = n_ptcl_in_cluster_disp_[i];
            const PS::S32 n_ptcl = n_ptcl_in_cluster_[i];
            PS::ReallocatableArray<PtclHard> extra_ptcl;
            driveForMultiClusterImpl(ptcl_hard_.getPointer(adr_head), n_ptcl, dt, extra_ptcl);
#pragma omp critical
            {
                for (PS::S32 j=0; j<extra_ptcl.size(); j++) {
                    PS::S32 adr = sys.getNumberOfParticleLocal();
                    PS::S32 rank = PS::Comm::getRank();
                    sys.addOneParticle(Tsptcl(extra_ptcl[j],rank,adr));
                }
            }
        }
    }

    template<class Tsys, class Tsptcl>
    void initialMultiClusterOMP(Tsys & sys, const PS::F64 dt_tree){
        const PS::S32 n_cluster = n_ptcl_in_cluster_.size();
        //	const PS::S32 ith = PS::Comm::getThreadNum();
#pragma omp for schedule(dynamic)
        for(PS::S32 i=0; i<n_cluster; i++){
            const PS::S32 adr_head = n_ptcl_in_cluster_disp_[i];
            const PS::S32 n_ptcl = n_ptcl_in_cluster_[i];
            SearchGroup<PtclHard> group;
            group.findGroups(ptcl_hard_.getPointer(adr_head), n_ptcl, n_split_);
            group.searchAndMerge(ptcl_hard_.getPointer(adr_head), r_bin_);
            PS::ReallocatableArray<PtclHard> ptcl_new;
            group.generateList(ptcl_hard_.getPointer(adr_head), ptcl_new, r_bin_, dt_tree, id_offset_, n_split_);
#pragma omp critical
            {
                for (PS::S32 j=0; j<ptcl_new.size(); j++) {
                    PS::S32 adr = sys.getNumberOfParticleLocal();
                    PS::S32 rank = PS::Comm::getRank();
                    sys.addOneParticle(Tsptcl(ptcl_new[j],rank,adr));
                }
            }
            
        }        
    }


};
