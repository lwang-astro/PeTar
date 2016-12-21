#pragma once
#ifdef USE_INTRINSIC_FOR_X86
#include<immintrin.h>
#endif

#ifndef NAN_CHECK
#define NAN_CHECK(val) assert((val) == (val));
#endif

#include"kepler.hpp"
#include"force.hpp"
#include"AR.h" /// include AR.h (L.Wang)
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
    PS::S32 id_cluster;
    PS::S32 adr_org;
    static PS::F64 r_factor;
    static PS::F64 dens;
    PtclHard():id(-1), mass(-1.0){}
    PtclHard(const PS::S64 _id, 
             const PS::F64 _mass, 
             const PS::F64vec & _pos, 
             const PS::F64vec & _vel,
             const PS::S32 _id_cluster,
             const PS::S32 _adr_org): id(_id), mass(_mass), pos(_pos), vel(_vel), 
                                     id_cluster(_id_cluster), adr_org(_adr_org){}

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
    ARC::chainpars ARC_control; ///chain controller (L.Wang)
#ifdef ARC_ERROR
    PS::F64 ARC_error_relative;
    PS::F64 ARC_error;  
    PS::S32 N_count[20];  // counting number of particles in one cluster
#endif
private:
    PS::F64 ARC_int_pars[2]; /// ARC integration parameters, rout_, rin_ (L.Wang)
    PS::F64 dt_limit_hard_;
    //PS::ReallocatableArray<PtclHard> ptcl_hard_;
    PS::ReallocatableArray<PS::S32> n_ptcl_in_cluster_;
    PS::ReallocatableArray<PS::S32> n_ptcl_in_cluster_disp_;
    PS::F64 time_origin_;

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

    struct OPLessIDCluster{
        template<class T> bool operator() (const T & left, const T & right) const {
            return left.id_cluster < right.id_cluster;
        }
    };


    template<class Tptcl>
    void driveForMultiClusterImpl(Tptcl * ptcl_org,
                                  const PS::S32 n_ptcl,
                                  const PS::F64 time_end){

      // kepler motion test
      if (n_ptcl==2) {
        PS::F64 ax,ecc,inc,OMG,omg,tperi;
        PosVel2OrbParam(ax,ecc,inc,OMG,omg,tperi,ptcl_org[0].pos, ptcl_org[1].pos, ptcl_org[0].vel, ptcl_org[1].vel, ptcl_org[0].mass, ptcl_org[1].mass);
        //        std::cerr<<"n_ptcl="<<n_ptcl<<"; ax="<<ax<<"; ecc="<<ecc<<"; peri="<<tperi<<"; pid="<<ptcl_org[0].id<<std::endl;
        if (ax>0.0&&2.0*ax<ARC_int_pars[1]) {
          // center-of-mass
          Tptcl pcm;
          calc_center_of_mass(pcm, ptcl_org, n_ptcl, true);
          
          DriveKepler(ptcl_org[0].mass, ptcl_org[1].mass,ptcl_org[0].pos, ptcl_org[1].pos, ptcl_org[0].vel, ptcl_org[1].vel,time_end);
          
          // integration of center-of-mass
          pcm.setPos(pcm.getPos()[0] + pcm.getVel()[0]*time_end,
                     pcm.getPos()[1] + pcm.getVel()[1]*time_end,
                     pcm.getPos()[2] + pcm.getVel()[2]*time_end);

          center_of_mass_correction(pcm, ptcl_org, n_ptcl);
          
          //          PosVel2OrbParam(ax,ecc,inc,OMG,omg,tperi,ptcl_org[0].pos, ptcl_org[1].pos, ptcl_org[0].vel, ptcl_org[1].vel, ptcl_org[0].mass, ptcl_org[1].mass);
          //          std::cerr<<"A:n_ptcl="<<n_ptcl<<"; ax="<<ax<<"; ecc="<<ecc<<"; peri="<<tperi<<"; pid"<<ptcl_org[0].id<<std::endl;
//          if (!kout.is_open()) kout.open("kout");
//          for (int i=0;i<n_ptcl;i++) kout<<std::setprecision(17)<<time_origin_<<" "<<ptcl_org[i].mass<<" "<<ptcl_org[i].pos<<" "<<ptcl_org[i].vel<<std::endl;
          return;
        }
      }
      
      /// start ARC (L.Wang)
      ARC::chain<Tptcl> c((std::size_t)n_ptcl,ARC_control);
      static thread_local PS::F64 time_sys = 0.0;

      c.addP(n_ptcl,ptcl_org);
      c.Int_pars=ARC_int_pars;
      c.init(time_sys);
#ifdef ARC_ERROR
      PS::F64 ARC_error_once = c.getPot()+c.getEkin();
      N_count[n_ptcl-1]++;
#endif
      
      PS::F64 ds_up_limit = calcDtLimit(time_sys, dt_limit_hard_)/c.calc_dt_X(1.0);
      //      PS::F64 dt_low_limit = ds_up_limit/256.0;
      PS::F64 ds_use = c.calc_next_step_peri();
      //      std::cerr<<"ds_use="<<ds_use<<std::endl;
      
      if (ds_use>ds_up_limit) ds_use = ds_up_limit;
      //      if (ds_use<dt_low_limit) ds_use = dt_low_limit;
      
      while(time_end-c.getTime()>ARC_control.dterr) {
//        if (ptcl_org[0].id==8&&time_origin_==0.0078125) {
//          std::cout<<"ds= "<<ds_use<<" toff= "<<time_end<<std::endl;
        
//          FILE* fout=fopen("data","w");
//          fwrite(ptcl_org,sizeof(Tptcl),n_ptcl,fout);
//          fclose(fout);
//          for (PS::S32 i=0; i<n_ptcl; i++) 
//            std::cout<<ptcl_org[i].getMass()<<" "<<ptcl_org[i].getPos()[0]<<" "<<ptcl_org[i].getPos()[1]<<" "<<ptcl_org[i].getPos()[2]<<" "<<ptcl_org[i].getVel()[0]<<" "<<ptcl_org[i].getVel()[1]<<" "<<ptcl_org[i].getVel()[2]<<std::endl;
//        }
        PS::F64 dsf=c.extrapolation_integration(ds_use,time_end);
        //        std::cerr<<"Particle="<<ptcl_org[0].id<<" n="<<n_ptcl<<" Time_end="<<time_end<<" ctime="<<c.getTime()<<" diff="<<time_end-c.getTime()<<" ds="<<ds_use<<" dsf="<<dsf<<std::endl;
        if (dsf<0) ds_use *= -dsf;
        else if (n_ptcl==2) ds_use *= dsf;
        else ds_use = dsf*c.calc_next_step_peri();
      }

      // error record
#ifdef ARC_ERROR
      PS::F64 ARC_error_temp = (c.getPot()+c.getEkin()-ARC_error_once);
      ARC_error += ARC_error_temp;
      ARC_error_relative += ARC_error_temp/ARC_error_once;
#endif
      
      // integration of center-of-mass
      c.cm.setPos(c.cm.getPos()[0] + c.cm.getVel()[0]*time_end,
				  c.cm.getPos()[1] + c.cm.getVel()[1]*time_end,
				  c.cm.getPos()[2] + c.cm.getVel()[2]*time_end);
      

      c.center_shift_inverse();
//      if (!arout.is_open()) arout.open("arout");
//      for (int i=0;i<n_ptcl;i++) arout<<std::setprecision(17)<<time_origin_<<" "<<ptcl_org[i].mass<<" "<<ptcl_org[i].pos<<" "<<ptcl_org[i].vel<<std::endl;
      /// end ARC (L.Wang)
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
    void setARCParam(const PS::F64 energy_error=1e-10, const PS::F64 dterr=1e-9, const PS::F64 dtmin=1e-24, const PS::S32 exp_method=1, const PS::S32 exp_itermax=20, const PS::S32 exp_fix_iter=0) {
      ARC_control.setA(Newtonian_cut_AW,Newtonian_cut_Ap);
      ARC_control.setabg(0,1,0);
      ARC_control.setEXP(energy_error,dtmin,dterr,exp_itermax,exp_method,3,(bool)exp_fix_iter);
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
            ptcl_hard_.push_back(PtclHard(p.id, p.mass, p.pos, p.vel, 
                                          med[i].id_cluster_, med[i].adr_sys_));
        }

        for(PS::S32 i=0; i<ptcl_recv.size(); i++){
            const Tptcl & p = ptcl_recv[i];
            ptcl_hard_.push_back(PtclHard(p.id_, p.mass_, p.pos_, p.vel_, 
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
                  const PS::F64 _dt_limit_hard,
                  const PS::F64 _time_origin){
        /// Set chain pars (L.Wang)
        ARC_int_pars[0] = _rout; 
        ARC_int_pars[1] = _rin;
        /// Set chain pars (L.Wang)        
        dt_limit_hard_ = _dt_limit_hard;
        time_origin_ = _time_origin;
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
            //n_ptcl_in_cluster_[i] = 1;
        }
    }

    void driveForOneCluster(const PS::F64 dt){
        const PS::S32 n = ptcl_hard_.size();
        for(PS::S32 i=0; i<n; i++){
            ptcl_hard_[i].pos += ptcl_hard_[i].vel * dt;
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
            ptcl_hard_[i].id   = sys[adr].id;
            ptcl_hard_[i].mass = sys[adr].mass;
            ptcl_hard_[i].pos  = sys[adr].pos;
            ptcl_hard_[i].vel  = sys[adr].vel;
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
            driveForMultiClusterImpl(ptcl_hard_.getPointer(adr_head), n_ptcl, dt);
        }
    }

    void driveForMultiClusterOMP(const PS::F64 dt){

        const PS::S32 n_cluster = n_ptcl_in_cluster_.size();
        //	const PS::S32 ith = PS::Comm::getThreadNum();
#pragma omp for schedule(dynamic)
	for(PS::S32 i=0; i<n_cluster; i++){
	    const PS::S32 adr_head = n_ptcl_in_cluster_disp_[i];
	    const PS::S32 n_ptcl = n_ptcl_in_cluster_[i];
            driveForMultiClusterImpl(ptcl_hard_.getPointer(adr_head), n_ptcl, dt);
	}
    }
};
