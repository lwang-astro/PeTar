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

std::ofstream fout_debug;

template<class T>
void Print(const T str, std::ostream & fout);



class SortAdr{
public:
    PS::F64 * time;
    SortAdr(PS::F64 * _time): time(_time){}
    bool operator() (const PS::S32 & left, const PS::S32 & right) const {
        return time[left] < time[right];
    }
};

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

class MediatorHard{
public:
    PS::S32 adr_sys_org_;
    PS::S32 rank_org_;
    PS::S32 adr_phard_;
};

class PtclPred{
public:
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64 r_merge;
};

class PtclForce{
public:
    PS::F64vec acc0; // pla + sun
    PS::F64vec acc1; // pla + sun
    PS::F64vec acc0_pla; // pla only
    PS::F64vec acc1_pla; // pla only
    //PS::F64 r2_ngb;
    void clear(){
        acc0 = acc1 = acc0_pla = acc1_pla = 0.0;
        //r2_ngb = PS::LARGE_FLOAT;
    }
};

class PtclH4{
public:
    PS::F64 mass;
    PS::F64 time;
    PS::F64 dt;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc0;
    PS::F64vec acc1;
    PS::F64vec acc0_pla;
    PS::F64vec acc1_pla;
    PS::F64 r_merge;
    PS::S32 n_ngb;
    PS::S32 id;

    PS::F64vec acc2; // for debug
    PS::F64vec acc3; // for debug
    PS::F64vec acc2_pla; // for debug
    PS::F64vec acc3_pla; // for debug
    static PS::F64 r_factor;
    static PS::F64 dens;

    void setRMerge(){
        static const PS::F64 rho_ave = dens * ( (1.49597871e13*1.49597871e13*1.49597871e13) / 1.989e33); // [Msun/AU^3]
        static const PS::F64 PI = 4.0*atan(1.0);
        static const PS::F64 C = 3.0/(4.0*PI*rho_ave);
        r_merge = cbrt(C*mass) * r_factor; // correct
    }

    void merge(PtclH4 & ptcl_del){
        pos = mass*pos + ptcl_del.mass*ptcl_del.pos;
        vel = mass*vel + ptcl_del.mass*ptcl_del.vel;
        mass = mass + ptcl_del.mass;
        pos /= mass;
        vel /= mass;
        dt = PS::LARGE_FLOAT;
        setRMerge();
    }

    void correctWithSun(const PtclForce & force,
                        const PS::F64 eta,
                        const PS::F64 dt_limit,
                        const PS::F64 a0_offset_sq=0.0){
        static const PS::F64 inv3 = 1.0 / 3.0;
        const PS::F64 h = 0.5 * dt;
        const PS::F64 hinv = 2.0 / dt;
        const PS::F64vec A0p = (force.acc0 + this->acc0);
        const PS::F64vec A0m = (force.acc0 - this->acc0);
        const PS::F64vec A1p = (force.acc1 + this->acc1)*h;
        const PS::F64vec A1m = (force.acc1 - this->acc1)*h;
        const PS::F64vec vel_new = this->vel + h*( A0p - inv3*A1m );
        this->pos += h*( (this->vel + vel_new) + h*(-inv3*A0m));
        this->vel = vel_new;
        this->acc0 = force.acc0;
        this->acc1 = force.acc1;
        this->time += dt;
        const PS::F64vec acc3 = (1.5*hinv*hinv*hinv) * (A1p - A0m);
        const PS::F64vec acc2 = (0.5*hinv*hinv) * A1m + h*acc3;
        const PS::F64vec A0m_pla = (force.acc0_pla - this->acc0_pla);
        const PS::F64vec A1p_pla = (force.acc1_pla + this->acc1_pla)*h;
        const PS::F64vec A1m_pla = (force.acc1_pla - this->acc1_pla)*h;
        const PS::F64vec acc3_pla = (1.5*hinv*hinv*hinv) * (A1p_pla - A0m_pla);
        const PS::F64vec acc2_pla = (0.5*hinv*hinv) * A1m_pla + h*acc3_pla;
        const PS::F64 dt_ref = std::min( CalcDt4th(this->acc0_pla, this->acc1_pla, acc2_pla, acc3_pla, eta, a0_offset_sq), CalcDt4th(this->acc0, this->acc1, acc2, acc3, eta) ); 
        // for debug
        this->acc2 = acc2;
        this->acc3 = acc3;
        this->acc2_pla = acc2_pla;
        this->acc3_pla = acc3_pla;
        // for debug
        this->acc0_pla = force.acc0_pla;
        this->acc1_pla = force.acc1_pla;
        const PS::F64 dt_old = this->dt;
        assert(dt_old != 0.0);
        this->dt = dt_limit;
        while(this->dt > dt_ref) this->dt *= 0.5;
        this->dt = dt_old*2 < this->dt ?  dt_old*2 : this->dt;
    }

    void correct(const PtclForce & force,
                 const PS::F64 eta,
                 const PS::F64 dt_limit,
                 const PS::F64 a0_offset_sq=0.0){
        static const PS::F64 inv3 = 1.0 / 3.0;
        const PS::F64 h = 0.5 * dt;
        const PS::F64 hinv = 2.0 / dt;
        const PS::F64vec A0p = (force.acc0 + this->acc0);
        const PS::F64vec A0m = (force.acc0 - this->acc0);
        const PS::F64vec A1p = (force.acc1 + this->acc1)*h;
        const PS::F64vec A1m = (force.acc1 - this->acc1)*h;
        const PS::F64vec vel_new = this->vel + h*( A0p - inv3*A1m );
        this->pos += h*( (this->vel + vel_new) + h*(-inv3*A0m));
        this->vel = vel_new;
        this->acc0 = force.acc0;
        this->acc1 = force.acc1;
        this->time += dt;
        const PS::F64vec acc3 = (1.5*hinv*hinv*hinv) * (A1p - A0m);
        const PS::F64vec acc2 = (0.5*hinv*hinv) * A1m + h*acc3;
        const PS::F64 dt_ref = CalcDt4th(this->acc0, this->acc1, acc2, acc3, eta, a0_offset_sq);
        // for debug
        this->acc2 = acc2;
        this->acc3 = acc3;
        this->acc2_pla = acc2_pla;
        this->acc3_pla = acc3_pla;
        // for debug
        this->acc0_pla = force.acc0_pla;
        this->acc1_pla = force.acc1_pla;
        const PS::F64 dt_old = this->dt;
        assert(dt_old != 0.0);
        this->dt = dt_limit;
        while(this->dt > dt_ref) this->dt *= 0.5;
        this->dt = dt_old*2 < this->dt ?  dt_old*2 : this->dt;
    }
};

PS::F64 PtclH4::r_factor = 1.0;
PS::F64 PtclH4::dens = 2.0; // [g/cc]

class PtclHard{
public:
    PS::S64 id;
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64vec vel;
#ifdef HERMITE  
    PS::F64 time;
    PS::F64 dt;
    PS::F64vec acc0;
    PS::F64vec acc1;
    PS::F64vec acc0_pla;
    PS::F64vec acc1_pla;
    PS::S32 n_ngb;
    //PS::F64 r_merge;
    PS::S32 rank_org;
    PS::S32 adr_fp;
#endif
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

class MergeLog{
public:
    PS::F64 time;
    PtclH4 ptcl_merged;
    PtclH4 ptcl_dead;
    MergeLog(){
        time = -1.0;
    }
    MergeLog(const PS::F64 t, const PtclH4 & p_m, const PtclH4 & p_d){
        time = t;
        ptcl_merged = p_m;
        ptcl_dead = p_d;
    }
    void dump(std::ostream & fout=std::cout){
        fout<<"time= "<<time<<std::endl;
        fout<<"ptcl_merged.id= "      <<ptcl_merged.id<<std::endl;
        fout<<"ptcl_merged.mass= "    <<ptcl_merged.mass<<std::endl;
        fout<<"ptcl_merged.pos= "     <<ptcl_merged.pos<<std::endl;
        fout<<"ptcl_merged.vel= "     <<ptcl_merged.vel<<std::endl;
        fout<<"ptcl_dead.id= "      <<ptcl_dead.id<<std::endl;
        fout<<"ptcl_dead.mass= "    <<ptcl_dead.mass<<std::endl;
        fout<<"ptcl_dead.pos= "     <<ptcl_dead.pos<<std::endl;
        fout<<"ptcl_dead.vel= "     <<ptcl_dead.vel<<std::endl;
    }
};


class MergeCandInfo{
public:
    PS::S32 id_cluster;
    PS::S32 n_cand;
    PS::S32 adr_ptcl;
    PS::S32 adr_merge_pair;
    MergeCandInfo(): id_cluster(-1), n_cand(0), adr_ptcl(-1), adr_merge_pair(-1){}
    MergeCandInfo(const PS::S32 _id, const PS::S32 _n, const PS::S32 _adr_p, const PS::S32 _adr_m):
        id_cluster(_id), n_cand(_n), adr_ptcl(_adr_p), adr_merge_pair(_adr_m){}
    MergeCandInfo(const MergeCandInfo & _mci):
        id_cluster(_mci.id_cluster), n_cand(_mci.n_cand), adr_ptcl(_mci.adr_ptcl), adr_merge_pair(_mci.adr_merge_pair){}
    void increaseNNgb(){ n_cand++; }
    void dump(std::ostream & fout=std::cout){
        fout<<"id_cluster= "<<id_cluster<<std::endl;
        fout<<"n_cand= "<<n_cand<<std::endl;
        fout<<"adr_ptcl= "<<adr_ptcl<<std::endl;
        fout<<"adr_merge_pair= "<<adr_merge_pair<<std::endl;
    }
};

class SystemHard{
public:
    PS::ReallocatableArray<PtclHard> ptcl_hard_;
    ARC::chainpars chain_control; ///chain controller (L.Wang)
    PS::F64 ARC_int_pars[2]; /// ARC integration parameters, rout_, rin_ (L.Wang)
private:
    PS::F64 eps2_;
    PS::F64 rin_;
    PS::F64 rout_;
    PS::F64 mass_sun_;
    PS::F64vec pos_sun_;
    PS::F64vec vel_sun_;
    PS::F64 eta_s_;
    PS::F64 eta_;
    PS::F64 dt_limit_hard_;
    PS::ReallocatableArray<MediatorHard> mediator_;
    //PS::ReallocatableArray<PtclHard> ptcl_hard_;
    PS::ReallocatableArray<PtclForce> ptcl_force_;
    PS::ReallocatableArray<PS::S32> n_ptcl_in_cluster_;
    PS::ReallocatableArray<PS::S32> n_ptcl_in_cluster_disp_;
    PS::ReallocatableArray<MergeLog> * merge_log_omp_;
    PS::F64 * eng_disp_merge_omp_;
    PS::F64 time_origin_;

    ///////////
    /// functor
    struct OPSortClusterID{
        template<class T> bool operator() (const T & left, const T & right) const {
            return left.id_cluster < right.id_cluster;
        }
    };
    struct OPSortFirst{
        template<class T> bool operator() (const T & left, const T & right) const {
            return left.first < right.first;
        }
    };
    struct OPSortSecond{
        template<class T> bool operator() (const T & left, const T & right) const {
            return left.second < right.second;
        }
    };

    ////////////
    // merge particles
    template<class Tptcl>
    void merge2body(PS::ReallocatableArray<Tptcl> & ptcl,
                    const PS::S32 adr0,
                    const PS::S32 adr1,
                    PS::ReallocatableArray<MergeLog> & merge_log,
                    PS::F64 & eng_disp){
        Tptcl & ptcl0 = ptcl[adr0];
        Tptcl & ptcl1 = ptcl[adr1];
        const PS::F64 m0 = ptcl0.mass;
        const PS::F64 m1 = ptcl1.mass;
        const PS::S64 id0 = ptcl0.id;
        const PS::S64 id1 = ptcl1.id;
        Tptcl & ptcl_merge = (m0>m1) ? ptcl0 : ( (m1>m0) ? ptcl1 : ((id0<id1) ? ptcl0 : ptcl1) );
        Tptcl & ptcl_dead = (&ptcl_merge == &ptcl0) ? ptcl1 : ptcl0;
        const PS::F64 merge_time = ptcl0.time + time_origin_;
        merge_log.push_back( MergeLog(merge_time, ptcl_merge, ptcl_dead) );
        PS::F64vec dr01 = pos_sun_ - ptcl_merge.pos;
        PS::F64vec dr02 = pos_sun_ - ptcl_dead.pos;
        PS::F64vec dr12 = ptcl_merge.pos - ptcl_dead.pos;
        PS::F64 eng_ini = -mass_sun_*ptcl_merge.mass / sqrt(dr01*dr01) - mass_sun_*ptcl_dead.mass / sqrt(dr02*dr02) - ptcl_merge.mass*ptcl_dead.mass / sqrt(dr12*dr12);
        ptcl_merge.merge(ptcl_dead);
        ptcl_dead.mass = 0.0;
        dr01 = pos_sun_ - ptcl_merge.pos;
        PS::F64 eng_fin = -mass_sun_*ptcl_merge.mass / sqrt(dr01*dr01);
        eng_disp += eng_fin - eng_ini;
    }
    // merge particles
    ////////////
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


    void setAcc0Acc1WithSunOneParticle(PS::ReallocatableArray<PtclForce> & force,
                                       const PS::ReallocatableArray<PtclPred> & pred,
                                       const PS::S32 adr_i,
                                       const PS::S32 n_tot){ 
        force[adr_i].acc0 = force[adr_i].acc1 = force[adr_i].acc0_pla = force[adr_i].acc1_pla = 0.0;
        //PS::F64 r2_min = PS::LARGE_FLOAT;
        for(PS::S32 j=0; j<n_tot; j++){
            if(adr_i == j) continue;
            PS::F64 r2 = 0.0;
            CalcAcc0Acc1R2Cutoff(pred[adr_i].pos,       pred[adr_i].vel,
                                 force[adr_i].acc0_pla, force[adr_i].acc1_pla, r2,
                                 pred[j].pos, pred[j].vel, pred[j].mass,
                                 eps2_, rout_, rin_);
        }
        //force[adr_i].r2_ngb = r2_min;
        CalcAcc0Acc1(pred[adr_i].pos,   pred[adr_i].vel,
                     force[adr_i].acc0, force[adr_i].acc1,
                     pos_sun_, vel_sun_, mass_sun_);
        force[adr_i].acc0 += force[adr_i].acc0_pla;
        force[adr_i].acc1 += force[adr_i].acc1_pla;
    }

    void setAcc0Acc1(PS::ReallocatableArray<PtclForce> & force,
                     const PS::ReallocatableArray<PtclPred> & pred,
                     const PS::S32 n_act, 
                     const PS::S32 n_tot, 
                     const PS::ReallocatableArray<PS::S32> & adr_array,
                     PS::ReallocatableArray< std::pair<PS::S32, PS::S32> > & merge_pair){
        for(PS::S32 i=0; i<n_act; i++){
            PS::S32 adr = adr_array[i];
            force[adr].acc0 = force[adr].acc1 = force[adr].acc0_pla = force[adr].acc1_pla = 0.0;
            for(PS::S32 j=0; j<n_tot; j++){
                if(adr == j) continue;
                PS::F64 r2 = 0.0;
                CalcAcc0Acc1R2Cutoff(pred[adr].pos,       pred[adr].vel,
                                     force[adr].acc0_pla, force[adr].acc1_pla, r2,
                                     pred[j].pos, pred[j].vel, pred[j].mass,
                                     eps2_, rout_, rin_);
                if(r2 < ((pred[adr].r_merge + pred[j].r_merge)*(pred[adr].r_merge + pred[j].r_merge)) && pred[j].mass > 0.0){
                    merge_pair.push_back( std::make_pair(adr, j) );
                }
            }
            force[adr].acc0 += force[adr].acc0_pla;
            force[adr].acc1 += force[adr].acc1_pla;
        }
    }

    void setAcc0Acc1WithSun(PS::ReallocatableArray<PtclForce> & force,
                            const PS::ReallocatableArray<PtclPred> & pred,
                            const PS::S32 n_act, 
                            const PS::S32 n_tot, 
                            const PS::ReallocatableArray<PS::S32> & adr_array,
                            PS::ReallocatableArray< std::pair<PS::S32, PS::S32> > & merge_pair){
        for(PS::S32 i=0; i<n_act; i++){
            PS::S32 adr = adr_array[i];
            force[adr].acc0 = force[adr].acc1 = force[adr].acc0_pla = force[adr].acc1_pla = 0.0;
            for(PS::S32 j=0; j<n_tot; j++){
                if(adr == j) continue;
                PS::F64 r2 = 0.0;
                CalcAcc0Acc1R2Cutoff(pred[adr].pos,       pred[adr].vel,
                                     force[adr].acc0_pla, force[adr].acc1_pla, r2,
                                     pred[j].pos, pred[j].vel, pred[j].mass,
                                     eps2_, rout_, rin_);
                if(r2 < ((pred[adr].r_merge + pred[j].r_merge)*(pred[adr].r_merge + pred[j].r_merge)) && pred[j].mass > 0.0){
                    merge_pair.push_back( std::make_pair(adr, j) );
                }
            }
            CalcAcc0Acc1(pred[adr].pos, pred[adr].vel,
                         force[adr].acc0, force[adr].acc1,
                         pos_sun_, vel_sun_, mass_sun_);
            force[adr].acc0 += force[adr].acc0_pla;
            force[adr].acc1 += force[adr].acc1_pla;
        }
    }

    template<class Tptcl>
    void setDt2ndOneParticle(PS::ReallocatableArray<Tptcl> & ptcl,
                             const PS::ReallocatableArray<PtclForce> & force,
                             const PS::S32 adr,
                             const PS::F64 eta,
                             const PS::F64 dt_limit,
                             const PS::F64 a0_offset_sq){
        const PS::F64vec a0_pla = force[adr].acc0_pla;
        const PS::F64vec a1_pla = force[adr].acc1_pla;
        const PS::F64vec a0_tot = force[adr].acc0;
        const PS::F64vec a1_tot = force[adr].acc1;
        const PS::F64 dt_ref = std::min( CalcDt2nd(a0_pla, a1_pla, eta, a0_offset_sq), CalcDt2nd(a0_tot, a1_tot, eta) );
        PS::F64 dt = dt_limit;
        while(dt > dt_ref) dt *= 0.5;
        ptcl[adr].dt = dt;
    }

    template<class Tptcl>
    void setDt2nd(PS::ReallocatableArray<Tptcl> & ptcl,
                  const PS::ReallocatableArray<PtclForce> & force,
                  const PS::ReallocatableArray<PS::S32> & adr_array,
                  const PS::S32 n_act, 
                  const PS::F64 eta,
                  const PS::F64 dt_limit,
                  const PS::F64 a0_offset_sq){
        for(PS::S32 i=0; i<n_act; i++){
            const PS::S32 adr = adr_array[i];
            const PS::F64vec a0_pla = force[adr].acc0_pla;
            const PS::F64vec a1_pla = force[adr].acc1_pla;
            const PS::F64vec a0_tot = force[adr].acc0;
            const PS::F64vec a1_tot = force[adr].acc1;
            const PS::F64 dt_ref = std::min( CalcDt2nd(a0_pla, a1_pla, eta, a0_offset_sq), CalcDt2nd(a0_tot, a1_tot, eta) );
            PS::F64 dt = dt_limit;
            while(dt > dt_ref) dt *= 0.5;
            ptcl[adr].dt = dt;
        }
    }

    void sortAndSelectIp(PS::ReallocatableArray<PS::S32> & adr_sorted,
                         PS::ReallocatableArray<PS::F64> & time_next,
                         const PS::ReallocatableArray<PtclH4> & ptcl,
                         PS::S32 & n_act,
                         const bool state_merge = false){
        const PS::S32 ni_old = n_act;
        const PS::S32 n_tot = time_next.size();
        //std::cerr<<"before sort"<<std::endl;
	/*
        for(PS::S32 ip=0; ip<ni_old; ip++){
            const PS::S32 adr = adr_sorted[ip];
            time_next[adr] += ptcl[adr].dt; // n_act only
        }
	*/
        if(state_merge){
            std::sort(adr_sorted.getPointer(), adr_sorted.getPointer(n_tot),  SortAdr(time_next.getPointer()));
        }
        else{
            std::sort(adr_sorted.getPointer(), adr_sorted.getPointer(ni_old), SortAdr(time_next.getPointer()));
        }
        const PS::F64 time_ref = time_next[adr_sorted[0]];
        for(n_act=1; n_act<n_tot; n_act++){
            if(time_ref < time_next[adr_sorted[n_act]]) {
                break;
            }
        }
    }

    void predictAll(PS::ReallocatableArray<PtclPred> & pred,
                    const PS::ReallocatableArray<PtclH4> & ptcl,
                    const PS::F64 time_next){
        static const PS::F64 inv3 = 1.0 / 3.0;
        PS::S32 n_tot = ptcl.size();
        for(PS::S32 i=0; i<n_tot; i++){
            const PS::F64 dt = time_next - ptcl[i].time;
            pred[i].pos = ptcl[i].pos + dt*(ptcl[i].vel  + 0.5*dt*(ptcl[i].acc0 + inv3*dt*ptcl[i].acc1));
            pred[i].vel = ptcl[i].vel + dt*(ptcl[i].acc0 + 0.5*dt*ptcl[i].acc1);
            pred[i].r_merge = ptcl[i].r_merge;
            pred[i].mass = ptcl[i].mass;
            //pred[i].setRMerge();
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

    void correctIp(PS::ReallocatableArray<PtclH4> & ptcl,
                   const PS::ReallocatableArray<PtclForce> & force,
                   const PS::ReallocatableArray<PS::S32> & adr_sorted, 
                   const PS::S32 n_act,
                   const PS::F64 time_sys,
                   const PS::F64 a0_offset_sq){
        const PS::F64 dt_limit = calcDtLimit(time_sys, dt_limit_hard_);
        for(PS::S32 i=0; i<n_act; i++){
            const PS::S32 adr = adr_sorted[i];
            ptcl[adr].correct(force[adr], eta_, dt_limit, a0_offset_sq);
        }
    }

    template<class Tptcl>
    void driveForMultiClusterImpl(Tptcl * ptcl_org,
                                  const PS::S32 n_ptcl,
                                  const PS::F64 time_end){
      
      /// start ARC (L.Wang)
      ARC::chain<Tptcl> c((std::size_t)n_ptcl,chain_control);
      static thread_local PS::F64 time_sys = 0.0;
      static thread_local PS::F64 dt_limit = calcDtLimit(time_sys, dt_limit_hard_);

      c.addP(n_ptcl,ptcl_org);
      c.Int_pars=ARC_int_pars;
      c.init(time_sys);
      PS::F64 dt_use = 0.2*c.calc_next_step_XVA();
      if (dt_use>dt_limit) dt_use = dt_limit;
      
      while(time_end-c.getTime()>chain_control.dterr) {
        //        std::cerr<<"Before ARC: N"<<n_ptcl<<" tend"<<time_end<<std::endl;
        PS::F64 dsf=c.extrapolation_integration(dt_use,time_end);
        //        std::cerr<<"After ARC: N"<<n_ptcl<<" t"<<c.getTime()<<" tend"<<time_end<<std::endl;        
        if (dsf<0) dt_use *= -dsf;
        //        else dt_use *= dsf;
      }

      c.center_shift_inverse();
      /// end ARC (L.Wang)
      
#ifdef HERMITE      
        static thread_local PS::ReallocatableArray<PtclH4> ptcl;
        ptcl.resizeNoInitialize(n_ptcl);
        static thread_local PS::ReallocatableArray<PtclPred> pred;
        pred.resizeNoInitialize(n_ptcl);
        static thread_local PS::ReallocatableArray<PtclForce> force;
        force.resizeNoInitialize(n_ptcl);
        static thread_local PS::ReallocatableArray<PS::S32> adr_sorted;
        adr_sorted.resizeNoInitialize(n_ptcl);
        static thread_local PS::ReallocatableArray<PS::F64> time_next;
        time_next.resizeNoInitialize(n_ptcl);
        static thread_local PS::ReallocatableArray< std::pair<PS::S32, PS::S32> > merge_pair;
        merge_pair.clearSize();
#ifdef MERGE
        static thread_local PS::ReallocatableArray< std::pair<PS::S32, PS::S32> > merge_pair_tmp;
        merge_pair_tmp.clearSize();
        static thread_local PS::ReallocatableArray<PS::S32> adr_merge_cand;
        adr_merge_cand.clearSize();
        static thread_local PS::ReallocatableArray<MergeCandInfo> merge_cand;
        merge_cand.resizeNoInitialize(n_ptcl);
        static thread_local PS::ReallocatableArray<MergeCandInfo> merge_cand_packed;
        merge_cand_packed.resizeNoInitialize(0);
#endif

        PS::F64 mass_min = PS::LARGE_FLOAT;
        for(PS::S32 i=0; i<n_ptcl; i++){
            pred[i].mass = ptcl[i].mass = ptcl_org[i].mass;
            pred[i].pos  = ptcl[i].pos = ptcl_org[i].pos;
            pred[i].vel  = ptcl[i].vel = ptcl_org[i].vel;
            ptcl[i].setRMerge();
            pred[i].r_merge = ptcl[i].r_merge;
            ptcl[i].id = ptcl_org[i].id;
            adr_sorted[i] = i;
            ptcl[i].time = ptcl[i].dt = 0.0;
            time_next[i] = 0.0;
            if(mass_min > pred[i].mass) mass_min = pred[i].mass;
        }

        Energy eng_init, eng_now;
        CalcEnergyHard(ptcl.getPointer(), n_ptcl, eng_init, rout_, rin_, mass_sun_, pos_sun_, vel_sun_, eps2_);
        PS::F64 time_sys = 0.0;
        PS::F64 dt_limit = calcDtLimit(time_sys, dt_limit_hard_);
        const PS::F64 a0_offset_sq = 0.1 * mass_min / (rout_ * rout_);
        PS::S32 n_act = n_ptcl;
        setAcc0Acc1(force, pred, n_act, n_ptcl, adr_sorted, merge_pair);

        for(PS::S32 i=0; i<n_ptcl; i++){
            ptcl[i].acc0 = force[i].acc0;
            ptcl[i].acc1 = force[i].acc1;
            ptcl[i].acc0_pla = force[i].acc0_pla;
            ptcl[i].acc1_pla = force[i].acc1_pla;
        }
        setDt2nd(ptcl, force, adr_sorted, n_act, eta_s_, dt_limit, a0_offset_sq);
        for(PS::S32 i=0; i<n_ptcl; i++){
            time_next[i] = ptcl[i].time + ptcl[i].dt;
        }
        sortAndSelectIp(adr_sorted, time_next, ptcl, n_act);
        PS::S32 n_loop = 0;
        while(time_sys != time_end){
            bool flag_merge = false;
            merge_pair.clearSize();
            n_loop++;
            time_sys = time_next[adr_sorted[0]];
            dt_limit = calcDtLimit(time_sys, dt_limit_hard_);
            predictAll(pred, ptcl, time_sys);
            setAcc0Acc1(force, pred, n_act, n_ptcl, adr_sorted, merge_pair);
            correctIp(ptcl, force, adr_sorted, n_act, time_sys, a0_offset_sq);
            for(PS::S32 i=0; i<n_act; i++){
                PS::S32 adr = adr_sorted[i];
                time_next[adr] = ptcl[adr].time + ptcl[adr].dt;
            }
#ifdef MERGE
            if(merge_pair.size() > 0){
                flag_merge = true;
                merge_pair_tmp.resizeNoInitialize(merge_pair.size());
                for(PS::S32 i=0; i<merge_pair.size(); i++){
                    merge_pair_tmp[i].first = merge_pair[i].second;
                    merge_pair_tmp[i].second = merge_pair[i].first;
                }
                std::sort(merge_pair.getPointer(),     merge_pair.getPointer(merge_pair.size()),         OPSortFirst());
                std::sort(merge_pair_tmp.getPointer(), merge_pair_tmp.getPointer(merge_pair_tmp.size()), OPSortFirst());
                for(PS::S32 i=0; i<merge_pair_tmp.size(); i++){
                    if(ptcl[merge_pair_tmp[i].first].time != time_sys ){
                        merge_pair.push_back(merge_pair_tmp[i]);
                    }
                }
                PS::S32 adr_ptcl_prev = -PS::LARGE_INT;
                adr_merge_cand.clearSize();
                for(PS::S32 i=0; i<merge_pair.size(); i++){
                    PS::S32 adr_ptcl_now = merge_pair[i].first;
                    if( adr_ptcl_now != adr_ptcl_prev){
                        adr_ptcl_prev = adr_ptcl_now;
                        merge_cand[adr_ptcl_now] = MergeCandInfo(ptcl[adr_ptcl_now].id, 0, adr_ptcl_now, i);
                        adr_merge_cand.push_back(adr_ptcl_now);
                    }
                    merge_cand[adr_ptcl_now].n_cand++;
                }
                // correct
                for(PS::S32 i=0; i<adr_merge_cand.size(); i++){
                    const PS::S32 adr = adr_merge_cand[i];
                    if(ptcl[adr].time == time_sys) continue;
                    ptcl[adr].dt = time_sys - ptcl[adr].time;
                    setAcc0Acc1OneParticle(force, pred, adr, n_ptcl);
                    ptcl[adr].correct(force[adr], eta_, dt_limit, a0_offset_sq);
                    time_next[adr] = ptcl[adr].time + ptcl[adr].dt;
                }
                //#ifdef BEFOR_CORRECT
                //fout_debug<<"C: label propagation method"<<std::endl;
                // grouping cluster using label propagating method
                while(1){
                    bool flag_itr = false;
                    for(PS::S32 i=0; i<adr_merge_cand.size(); i++){
                        const PS::S32 adr_i = adr_merge_cand[i];
                        PS::S32 & id_cluster_i = merge_cand[adr_i].id_cluster;
                        for(PS::S32 j=0; j<merge_cand[adr_i].n_cand; j++){
                            //for(PS::S32 j=adr_j; j<adr_j+merge_cand[adr_i].n_cand; j++){
                            const PS::S32 adr_j = merge_pair[merge_cand[adr_i].adr_merge_pair+j].second;
                            //fout_debug<<"adr_i= "<<adr_i<<" adr_j= "<<adr_j<<std::endl;
                            assert(merge_pair[merge_cand[adr_i].adr_merge_pair].first == adr_i);
                            const PS::S32 id_cluster_j = merge_cand[adr_j].id_cluster;
                            if(id_cluster_i > id_cluster_j){
                                id_cluster_i = id_cluster_j;
                                flag_itr = true;
                            }
                        }
                    }
                    if(!flag_itr) break;
                }
                // pack merge_cand_packed
                merge_cand_packed.clearSize();
                for(PS::S32 i=0; i<adr_merge_cand.size(); i++){
                    merge_cand_packed.push_back(merge_cand[adr_merge_cand[i]]);
                }
                std::sort(merge_cand_packed.getPointer(), merge_cand_packed.getPointer(merge_cand_packed.size()), OPSortClusterID());
                fout_debug<<std::endl;
                //fout_debug<<"F: BEFORE MERGE"<<std::endl;
                //#ifdef BEFOR_CORRECT
                if(merge_cand_packed.size() > 0){
                    PS::S32 id_cluster_prev = merge_cand_packed[0].id_cluster;
                    PS::S32 id_max = -PS::LARGE_INT;
                    PS::S32 adr_merged = -PS::LARGE_INT;
                    PS::F64 mass_max = -PS::LARGE_FLOAT;
                    //PS::S32 adr_prev = 0;
                    PS::S32 n_offset = 0;
                    PS::S32 n_cnt = 0;
                    PS::F64 mass_cm = 0.0;
                    PS::F64vec pos_cm = 0.0;
                    PS::F64vec vel_cm = 0.0;
                    PS::S32 ith = PS::Comm::getThreadNum();
                    //#ifdef MERGE_DEBUG
#if 1
                    for(PS::S32 i=0; i<merge_cand_packed.size()+1; i++){
                        if(i == merge_cand_packed.size() || merge_cand_packed[i].id_cluster != id_cluster_prev){
                            for(PS::S32 j0=n_offset; j0<n_offset+n_cnt; j0++){
                                PS::S32 adr0 = merge_cand_packed[j0].adr_ptcl;
                                pred[adr0].mass = ptcl[adr0].mass;
                                pred[adr0].pos  = ptcl[adr0].pos;
                                pred[adr0].vel  = ptcl[adr0].vel;
                            }
                            PS::F64 eng_bef = 0.0;
#ifdef FROM_ONLY_NEIGHOBOR
                            for(PS::S32 j0=n_offset; j0<n_offset+n_cnt; j0++){
                                PS::S32 adr0 = merge_cand_packed[j0].adr_ptcl;
                                //fout_debug<<"adr0= "<<adr0<<" n_offset= "<<n_offset<<" adr_merged= "<<adr_merged<<" n_cnt="<<n_cnt<<std::endl;
                                const PS::F64vec drsun = pos_sun_ - ptcl[adr0].pos;
                                eng_bef += ptcl[adr0].mass*(0.5*ptcl[adr0].vel*ptcl[adr0].vel - mass_sun_/sqrt(drsun*drsun));
                                //fout_debug<<"eng_bef(4)= "<<eng_bef<<std::endl;
                                if(adr0 == adr_merged) continue;
                                merge_log_omp_[ith].push_back( MergeLog(time_sys+time_origin_, ptcl[adr_merged], ptcl[adr0]) );
                                const PS::F64vec drtarget  = ptcl[adr_merged].pos - ptcl[adr0].pos;
                                //fout_debug<<"sqrt(drtarget*drtarget)="<<sqrt(drtarget*drtarget)<<std::endl;
                                eng_bef -= ptcl[adr_merged].mass*ptcl[adr0].mass/(sqrt(drtarget*drtarget));
                                //fout_debug<<"eng_bef(3)= "<<eng_bef<<" ptcl[adr_merged].mass*ptcl[adr0].mass/(sqrt(drtarget*drtarget))= "
                                <<ptcl[adr_merged].mass*ptcl[adr0].mass/(sqrt(drtarget*drtarget))<<std::endl;
                                for(PS::S32 j1=j0+1; j1<n_offset+n_cnt; j1++){
                                    PS::S32 adr1 = merge_cand_packed[j1].adr_ptcl;
                                    if(adr1 == adr_merged) continue;
                                    const PS::F64vec dr01  = ptcl[adr1].pos - ptcl[adr0].pos;
                                    eng_bef -= ptcl[adr1].mass*ptcl[adr0].mass/(sqrt(dr01*dr01));
                                }
                                //fout_debug<<"eng_bef(2)= "<<eng_bef<<std::endl;
                            }
#else //FROM_ONLY_NEIGHOBOR
                            for(PS::S32 j0=n_offset; j0<n_offset+n_cnt; j0++){
                                Energy eng_one_particle;
                                PS::S32 adr0 = merge_cand_packed[j0].adr_ptcl;
                                CalcEnergyHardOneParticle(pred.getPointer(), n_ptcl, eng_one_particle, rout_, rin_, mass_sun_, pos_sun_, vel_sun_, adr0);
                                //fout_debug<<"adr0= "<<adr0<<" eng_one_particle.tot= "<<eng_one_particle.tot<<std::endl;
                                eng_bef += eng_one_particle.tot;
                            }
                            for(PS::S32 j0=n_offset; j0<n_offset+n_cnt; j0++){
                                PS::S32 adr0 = merge_cand_packed[j0].adr_ptcl;
                                for(PS::S32 j1=j0+1; j1<n_offset+n_cnt; j1++){
                                    PS::S32 adr1 = merge_cand_packed[j1].adr_ptcl;
                                    const PS::F64vec dr01  = ptcl[adr1].pos - ptcl[adr0].pos;
                                    const PS::F64 dr01_s = sqrt(dr01*dr01);
                                    //fout_debug<<"dr01_s= "<<dr01_s<<" rin_= "<<rin_<<std::endl;
                                    //eng_bef += ptcl[adr1].mass*ptcl[adr0].mass/(sqrt(dr01*dr01));
                                    eng_bef += ptcl[adr1].mass*ptcl[adr0].mass/dr01_s*(1.0 - CalcW(dr01_s/rout_, rin_/rout_));
                                }
                            }
#endif //FROM_ONLY_NEIGHOBOR
                            pos_cm /= mass_cm;
                            vel_cm /= mass_cm;
                            ///* for debug
                            ptcl[adr_merged].mass = mass_cm;
                            ptcl[adr_merged].pos = pos_cm;
                            ptcl[adr_merged].vel = vel_cm;
                            ptcl[adr_merged].time = time_sys;
                            ptcl[adr_merged].setRMerge();
                            pred[adr_merged].mass = ptcl[adr_merged].mass; // To copy to pred is required to calc force accurately (setAcc0And... use pred[]) 
                            pred[adr_merged].pos  = ptcl[adr_merged].pos;
                            pred[adr_merged].vel  = ptcl[adr_merged].vel;
                            pred[adr_merged].r_merge = ptcl[adr_merged].r_merge;
                            //for debug */
                            for(PS::S32 j0=n_offset; j0<n_offset+n_cnt; j0++){
                                PS::S32 adr0 = merge_cand_packed[j0].adr_ptcl;
                                //fout_debug<<"adr0= "<<adr0<<" adr_merged="<<adr_merged<<std::endl;
                                //fout_debug<<"ptcl[adr0].mass= "<<ptcl[adr0].mass<<" ptcl[adr0].pos= "<<ptcl[adr0].pos<<std::endl;
                                if(adr0 == adr_merged) continue;
                                merge_log_omp_[ith].push_back( MergeLog(time_next[0]+time_origin_, ptcl[adr_merged], ptcl[adr0]) );
                                ///* for debug
                                ptcl[adr0].mass = 0.0;
                                ptcl[adr0].dt = PS::LARGE_FLOAT;
                                time_next[adr0] = PS::LARGE_FLOAT;
                                ptcl[adr0].pos = PS::LARGE_FLOAT;
                                pred[adr0].mass = ptcl[adr0].mass; // To copy to pred is required to calc force accurately (setAcc0And... use pred[]) 
                                pred[adr0].pos  = ptcl[adr0].pos;
                                pred[adr0].vel  = ptcl[adr0].vel;
                                //for debug */
                            }
			    /*
                            pred[adr_merged].mass = ptcl[adr_merged].mass;
                            pred[adr_merged].pos  = ptcl[adr_merged].pos;
                            pred[adr_merged].vel  = ptcl[adr_merged].vel;
			    */
#ifdef FROM_ONLY_NEIGHOBOR
                            PS::F64vec drtarget = ptcl[adr_merged].pos - pos_sun_;
                            PS::F64 eng_fin = ptcl[adr_merged].mass * (-mass_sun_/(sqrt(drtarget*drtarget)) + 0.5*ptcl[adr_merged].vel*ptcl[adr_merged].vel);
#else //FROM_ONLY_NEIGHOBOR
                            Energy eng_one_particle_fin;
                            CalcEnergyHardOneParticle(pred.getPointer(), n_ptcl, eng_one_particle_fin, rout_, rin_, mass_sun_, pos_sun_, vel_sun_, adr_merged);
                            PS::F64 eng_fin = eng_one_particle_fin.tot;
#endif //FROM_ONLY_NEIGHOBOR
                            eng_disp_merge_omp_[ith] += 0.0;
                            ///* for debug
                            eng_disp_merge_omp_[ith] -= eng_fin - eng_bef; // maybe correct
                            setAcc0Acc1OneParticle(force, pred, adr_merged, n_ptcl); // new
                            ptcl[adr_merged].acc0     = force[adr_merged].acc0;
                            ptcl[adr_merged].acc1     = force[adr_merged].acc1;
                            ptcl[adr_merged].acc0_pla = force[adr_merged].acc0_pla;
                            ptcl[adr_merged].acc1_pla = force[adr_merged].acc1_pla;
                            setDt2ndOneParticle(ptcl, force, adr_merged, eta_s_, dt_limit, a0_offset_sq);
                            time_next[adr_merged] = ptcl[adr_merged].time + ptcl[adr_merged].dt;
                            if(i==merge_cand_packed.size()) break;
                            mass_max = -PS::LARGE_FLOAT;
                            n_cnt = 0;
                            mass_cm = 0.0;
                            pos_cm = 0.0;
                            vel_cm = 0.0;
                            id_max = ptcl[merge_cand_packed[i].adr_ptcl].id;
                            n_offset += n_cnt;
                            adr_merged = merge_cand_packed[i].adr_ptcl;
                            id_cluster_prev = merge_cand_packed[i].id_cluster;
                        } //if (i == merge_cand_packed.size() || ...)
                        PS::S32 adr_now = merge_cand_packed[i].adr_ptcl;
                        PS::S32 id_now = ptcl[adr_now].id;
                        if(ptcl[adr_now].mass > mass_max){
                            mass_max = ptcl[adr_now].mass;
                            id_max = id_now;
                            adr_merged = adr_now;
                        }
                        else if(ptcl[adr_now].mass == mass_max){
                            if(id_max > ptcl[adr_now].id){
                                id_max = id_now;
                                adr_merged = adr_now;
                            }
                        }
                        //fout_debug<<"ptcl[adr_now].pos= "<<ptcl[adr_now].pos<<std::endl;
                        mass_cm += ptcl[adr_now].mass;
                        pos_cm  += ptcl[adr_now].mass*ptcl[adr_now].pos;
                        vel_cm  += ptcl[adr_now].mass*ptcl[adr_now].vel;
                        n_cnt++;

                    } //for(merge_cand_packed.size)
                    //fout_debug<<std::endl;
#endif //#if 1

                } // end of if(merge)
            }
#endif // MERGE
            sortAndSelectIp(adr_sorted, time_next, ptcl, n_act, flag_merge);
        }
        CalcEnergyHard(ptcl.getPointer(), n_ptcl, eng_now, rout_, rin_, mass_sun_, pos_sun_, vel_sun_, eps2_);
        for(PS::S32 i=0; i<n_ptcl; i++){
            ptcl_org[i].mass = ptcl[i].mass;
            ptcl_org[i].pos  = ptcl[i].pos;
            ptcl_org[i].vel  = ptcl[i].vel;
        }
#endif
    }

public:

    SystemHard(){
        PS::S32 n_threads = PS::Comm::getNumberOfThread();
        merge_log_omp_ = new PS::ReallocatableArray<MergeLog>[n_threads];
        eng_disp_merge_omp_ = new PS::F64[n_threads];
    }

    /// start set Chainpars (L.Wang)
    ///
    void setARCParam(const PS::F64 energy_error=1e-12, const PS::F64 dterr=1e-10, const PS::F64 dtmin=1e-24, const PS::S32 exp_method=1, const PS::S32 exp_itermax=20, const PS::S32 exp_fix_iter=0) {
      chain_control.setA(Newtonian_cut_AW,Newtonian_cut_Ap);
      chain_control.setabg(0,1,0);
      chain_control.setEXP(energy_error,dtmin,dterr,exp_itermax,exp_method,3,(bool)exp_fix_iter);
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
	//fout_debug<<"(in setPtcl 1)ptcl_hard_.size()= "<<ptcl_hard_.size()<<std::endl;
        for(PS::S32 i=0; i<ptcl_recv.size(); i++){
            const Tptcl & p = ptcl_recv[i];
            ptcl_hard_.push_back(PtclHard(p.id_, p.mass_, p.pos_, p.vel_, 
                                          p.id_cluster_, -(i+1)));
        }
	//fout_debug<<"(in setPtcl 2)ptcl_hard_.size()= "<<ptcl_hard_.size()<<std::endl;
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

    void setSun(const PS::F64 m, const PS::F64vec & p, const PS::F64vec & v){
        mass_sun_ = m;
        pos_sun_ = p;
        vel_sun_ = v;
    }

    void setTimeOrigin(const PS::F64 _time_origin){
        time_origin_ = _time_origin;
    }

    void setParam(const PS::F64 _rout, 
                  const PS::F64 _rin, 
                  const PS::F64 _dt_limit_hard,
                  const PS::F64 _eta,
                  const PS::F64 _eta_s,
                  const PS::F64 _time_origin,
                  const PS::F64 _eps2=0.0){
        rout_ = _rout;
        rin_  = _rin;
        /// Set chain pars (L.Wang)
        ARC_int_pars[0] = rout_; 
        ARC_int_pars[1] = rin_;
        /// Set chain pars (L.Wang)        
        dt_limit_hard_ = _dt_limit_hard;
        eta_ = _eta;
        eta_s_ = _eta_s;
        time_origin_ = _time_origin;
        eps2_ = _eps2;
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
	const PS::S32 ith = PS::Comm::getThreadNum();
	//eng_disp_merge_omp_[ith] = 0.0;
	//merge_log_omp_[ith].clearSize();
#pragma omp for schedule(dynamic)
	for(PS::S32 i=0; i<n_cluster; i++){
	    const PS::S32 adr_head = n_ptcl_in_cluster_disp_[i];
	    const PS::S32 n_ptcl = n_ptcl_in_cluster_[i];
            driveForMultiClusterImpl(ptcl_hard_.getPointer(adr_head), n_ptcl, dt);
	}
    }
// for multi cluster
//////////////////

    void dumpMergeLog(std::ostream & fout){
        const PS::S32 n_th = PS::Comm::getNumberOfThread();
        for(PS::S32 ith=0; ith<n_th; ith++){
            const PS::S32 log_size = merge_log_omp_[ith].size();
            for(PS::S32 ilog=0; ilog<log_size; ilog++){
                merge_log_omp_[ith][ilog].dump(fout);
            }
        }
    }

    PS::F64 getEngDispMerge(){
        const PS::S32 n_th = PS::Comm::getNumberOfThread();
        PS::F64 eng_disp = 0.0;
        for(PS::S32 ith=0; ith<n_th; ith++){
            //fout_debug<<"eng_disp_merge_omp_[ith]="<<eng_disp_merge_omp_[ith]<<std::endl;
            eng_disp += eng_disp_merge_omp_[ith];
        }
        return eng_disp;
    }

};
