#pragma once
#ifdef USE_INTRINSIC_FOR_X86
#include<immintrin.h>
#endif

#ifndef NAN_CHECK
#define NAN_CHECK(val) assert((val) == (val));
#endif

const PS::F64 SAFTY_OFFSET_FOR_SEARCH = 1e-7;

#include"integrate.hpp"
#include"soft.hpp"
#include"cstdlib"
//#include"stdio.h" /// for debug (L.Wang)

template<class T>
void Print(const T str, std::ostream & fout);

//std::ofstream kout;
//std::ofstream arout;

class PtclHard{
public:
    PS::S64 id;
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64 r_out;
    // PS::S32 n_ngb;
    PS::S32 id_cluster;
    PS::S32 adr_org;
    PS::S32 adr_group;
    static PS::F64 r_factor;
    static PS::F64 dens;
    PtclHard():id(-1), mass(-1.0), adr_group(-1) {}
    //    PtclHard(const PtclHard &p) {
    //        id         = p.id;
    //        mass       = p.mass;
    //        pos        = p.pos;
    //        vel        = p.vel;
    //        r_out      = p.r_out;
    //        id_cluster = p.id_cluster;
    //        adr_org    = p.adr_org;
    //    }
    PtclHard(const PS::S64 _id, 
             const PS::F64 _mass, 
             const PS::F64vec & _pos, 
             const PS::F64vec & _vel,
             const PS::F64 _r_out,
             //  const PS::S32 _n_ngb,
             const PS::S32 _id_cluster,
             const PS::S32 _adr_org,
             const PS::S32 _adr_group=-1): id(_id), mass(_mass), pos(_pos), vel(_vel), r_out(_r_out),
                                        id_cluster(_id_cluster), adr_org(_adr_org), adr_group(_adr_group) {}


    /// start Particle member function (L.Wang)
    //! Get mass (required for \ref ARC::chain)
    /*! \return mass
     */
    const PS::F64 getMass() const{
        return mass;
    }
  
    //! Get position (required for \ref ARC::chain)
    /*! \return position vector (PS::F64[3])
     */
    const PS::F64* getPos() const{
        return &pos[0];
    }

    //! Get velocity (required for \ref ARC::chain)
    /*! \return velocity vector (PS::F64[3])
     */
    const PS::F64* getVel() const{
        return &vel[0];
    }

    //!Set position (required for \ref ARC::chain)
    /*! NAN check will be done
      @param [in] x: particle position in x axis
      @param [in] y: particle position in y axis
      @param [in] z: particle position in z axis
    */
    void setPos(const PS::F64 x, const PS::F64 y, const PS::F64 z) {
        NAN_CHECK(x);
        NAN_CHECK(y);
        NAN_CHECK(z);
    
        pos[0] = x;
        pos[1] = y;
        pos[2] = z;
    }
  
    //!Set velocity (required for \ref ARC::chain)
    /*! NAN check will be done
      @param [in] vx: particle velocity in x axis
      @param [in] vy: particle velocity in y axis 
      @param [in] vz: particle velocity in z axis 
    */
    void setVel(const PS::F64 vx, const PS::F64 vy, const PS::F64 vz) {
        NAN_CHECK(vx);
        NAN_CHECK(vy);
        NAN_CHECK(vz);
    
        vel[0] = vx;
        vel[1] = vy;
        vel[2] = vz;
    }

    //!Set mass (required for \ref ARC::chain)
    /*! NAN check will be done
      @param [in] m: particle mass
    */
    void setMass(const PS::F64 m) {
        NAN_CHECK(m);

        mass = m;
    }
    /// end Particle member function (L.Wang)

};

class SystemHard{
public:
    PS::ReallocatableArray<PtclHard> ptcl_hard_;
    ARC::chainpars<PtclHard,ARC_int_pars> ARC_control; ///chain controller (L.Wang)
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
    PS::F64 r_out_single_;
    PS::F64 m_average_;

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


    template<class Tptcl>
    void driveForMultiClusterImpl(Tptcl * ptcl_org,
                                  const PS::S32 n_ptcl,
                                  const PS::F64 time_end,
                                  ReallocatableArray<std::vector<Tptcl>> & group_ptcl_glb,
                                  ReallocatableArray<PS::S32> & group_ptcl_glb_empty_list) {

#ifdef HERMITE
        if(n_ptcl>5) {
            SearchGroup group;
            group.searchPerturber(ptcl_org, n_ptcl);
            
            HermiteIntegrator Hint(n_ptcl);
            Hint.setPtclList(ptcl_org,n_ptcl,group.getPerts());
            Hint.setParams(dt_limit_hard_, eta_s_, Int_pars.rin, Int_pars.eps2);
            PS::S32 group_act_n = 0;
            ReallocatableArray<PS::S32> group_act_list; //active group_list act adr
            ReallocatableArray<PS::S32> group_list;     //group.adr list
            ReallocatableArray<PS::S32> adr_group;      //ptcl -> group.adr [non cm is -1] (value of Ptcl.adr_group)
            ReallocatableArray<PS::S32> adr_group_map;  //ptcl -> group_list index [non cm is -1]
            ReallocatableArray<PS::S32> adr_cm;         //group_list index -> ptcl.cm
            group.findGroups(group_list, adr_group, adr_group_map,  adr_cm, group_act_n, ptcl_org, n_ptcl);
            group_act_list.resizeNoInitialize(group_act_n);
            
            ARCIntegrator Aint();
            Aint.initialize(group_list.getPointer(), group_ptcl_glb.getPointer(), group_act_n, ptcl_org, adr_cm.getPointer(), group.getPerts(), ARC_control, Int_pars); // here act_list used as adr_list
            
            PS::S32 time_sys=0.0;
            PS::F64 dt_limit;
            Hint.initialize(dt_limit, group_act_list.getPointer(), group_act_n, adr_group.getPointer(), adr_group_map.getPointer());
            while(time_sys<time_end) {
                time_sys = Hint.getNextTime();
                Aint.integrateOneStepList(group_act_list.getPointer(), group_act_n, time_sys, dt_limit);
                Hint.integrateOneStep(dt_limit, group_act_list.getPointer(), group_act_n, adr_group.getPointer(), adr_group_map.getPointer());
                Aint.updateCM(ptcl_org, group_act_list.getPointer(), adr_cm.getPointer(), group_act_n);
            }
            Aint.CMshiftreverse();
            
            group.resolveGroups(ptcl_org, n_ptcl, group_ptcl_glb.getPointer(), group_list.size(), group_list.getPointer(), adr_cm.getPointer());
            group.searchaAndMerge(ptcl_org, n_ptcl, Int_pars.rin);
            // Kickcorrect(ptcl_org, group.getRoutChangeList());
            group.generatelist(ptcl_org, n_ptcl, group_list,  group_ptcl_glb, group_ptcl_glb_empty_list);

            // group.reverseCopy(ptcl_org, n_ptcl);
        }
        else 
#endif
            Multiple_integrator(ptcl_org, n_ptcl, time_end, dt_limit_hard_,
                                r_out_single_, gamma_, m_average_,
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
        ARC_control.setA(Newtonian_cut_AW<PtclHard>,Newtonian_cut_Ap<PtclHard>,Newtonian_timescale);
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
            const FPSoft & p = sys[med[i].adr_sys_];
            // ptcl_hard_.push_back(PtclHard(p.id, p.mass, p.pos, p.vel, p.r_out, p.n_ngb,
            ptcl_hard_.push_back(PtclHard(p.id, p.mass, p.pos, p.vel, p.r_out,
                                          med[i].id_cluster_, med[i].adr_sys_));
        }

        for(PS::S32 i=0; i<ptcl_recv.size(); i++){
            const Tptcl & p = ptcl_recv[i];
            //  ptcl_hard_.push_back(PtclHard(p.id_, p.mass_, p.pos_, p.vel_, p.r_out_, p.n_ngb_,
            ptcl_hard_.push_back(PtclHard(p.id_, p.mass_, p.pos_, p.vel_, p.r_out_, 
                                          p.id_cluster_, -(i+1)));
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

    void setParam(const PS::F64 _rout, 
                  const PS::F64 _rin,
                  const PS::F64 _eps,
                  const PS::F64 _dt_limit_hard,
                  const PS::F64 _eta,
                  const PS::F64 _time_origin,
                  const PS::F64 _gmin,
                  const PS::F64 _m_avarage){
        /// Set chain pars (L.Wang)
		Int_pars.rin  = _rin;
        Int_pars.eps2  = _eps*_eps;
        /// Set chain pars (L.Wang)        
        dt_limit_hard_ = _dt_limit_hard;
        eta_s_ = _eta*_eta;
        time_origin_ = _time_origin;
        gamma_ = std::pow(1.0/_gmin,0.33333);
        r_out_single_ = _rout; 
        m_average_ = _m_avarage;
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
            ptcl_hard_[i].id   = sys[adr].id;
            ptcl_hard_[i].mass = sys[adr].mass;
            ptcl_hard_[i].pos  = sys[adr].pos;
            ptcl_hard_[i].vel  = sys[adr].vel;
            ptcl_hard_[i].r_out= sys[adr].r_out;
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
            ptcl_hard_[i].id   = sys[adr].id;
            ptcl_hard_[i].mass = sys[adr].mass;
            ptcl_hard_[i].pos  = sys[adr].pos;
            ptcl_hard_[i].vel  = sys[adr].vel;
            ptcl_hard_[i].r_out= sys[adr].r_out;
            //n_ptcl_in_cluster_[i] = 1;
        }
    }

    void driveForOneCluster(const PS::F64 dt){
        const PS::S32 n = ptcl_hard_.size();
        for(PS::S32 i=0; i<n; i++){
            ptcl_hard_[i].pos += ptcl_hard_[i].vel * dt;
            ptcl_hard_[i].r_out= r_out_single_;
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
            ptcl_hard_[i].pos += ptcl_hard_[i].vel * dt;
            ptcl_hard_[i].r_out= r_out_single_;
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
            sys[adr].mass = ptcl_hard_[i].mass;
            sys[adr].pos  = ptcl_hard_[i].pos;
            sys[adr].vel  = ptcl_hard_[i].vel;
            sys[adr].r_out = ptcl_hard_[i].r_out;
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
            sys[adr].mass = ptcl_hard_[i].mass;
            sys[adr].pos  = ptcl_hard_[i].pos;
            sys[adr].vel  = ptcl_hard_[i].vel;
            sys[adr].r_out = ptcl_hard_[i].r_out;
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
            ptcl_hard_[i].id   = sys[adr].id;
            ptcl_hard_[i].mass = sys[adr].mass;
            ptcl_hard_[i].pos  = sys[adr].pos;
            ptcl_hard_[i].vel  = sys[adr].vel;
            ptcl_hard_[i].r_out= sys[adr].r_out;
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
            ptcl_hard_[i].id    = sys[adr].id;
            ptcl_hard_[i].mass  = sys[adr].mass;
            ptcl_hard_[i].pos   = sys[adr].pos;
            ptcl_hard_[i].vel   = sys[adr].vel;
            ptcl_hard_[i].r_out = sys[adr].r_out;
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
    void driveForMultiCluster(const PS::F64 dt,
                              ReallocatableArray<std::vector<PtclHard> & group_ptcl_glb){
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
            driveForMultiClusterImpl(ptcl_hard_.getPointer(adr_head), n_ptcl, dt,group_ptcl_glb, group_ptcl_glb_empty_list);
        }
    }

    void driveForMultiClusterOMP(const PS::F64 dt,
                                 ReallocatableArray<std::vector<PtclHard>> & group_ptcl_glb,
                                 ReallocatableArray<PS::S32> & group_ptcl_glb_empty_list){

        const PS::S32 n_cluster = n_ptcl_in_cluster_.size();
        //	const PS::S32 ith = PS::Comm::getThreadNum();
#pragma omp for schedule(dynamic)
        for(PS::S32 i=0; i<n_cluster; i++){
            const PS::S32 adr_head = n_ptcl_in_cluster_disp_[i];
            const PS::S32 n_ptcl = n_ptcl_in_cluster_[i];
            driveForMultiClusterImpl(ptcl_hard_.getPointer(adr_head), n_ptcl, dt, group_ptcl_glb, group_ptcl_glb_empty_list);
            
        }
    }


};

template <class Tptcl>
class systemGroup{
public:
    ReallocatableArray<ReallocatableArray<Tptcl>> groups;     //data
    ReallocatableArray<PS::S32> cm_adr;                       //c.m. index in original particle array
    
};
