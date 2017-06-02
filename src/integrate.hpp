#pragma once

#include "kepler.hpp"
#include "hard_force.hpp"
#include "rsearch.hpp"
#include "AR.h" /// include AR.h (L.Wang)

#ifdef HARD_DEBUG
#define DEBUG_ENERGY_LIMIT 1e-6
#endif

//leap frog----------------------------------------------
template<class Tpsys, class Ttree>
void Kick(Tpsys & system,
          const Ttree & tree,
          const PS::F64 dt){
    const PS::S32 n = system.getNumberOfParticleLocal();
#pragma omp parallel for
    for(int i=0; i<n; i++){
	system[i].vel  += system[i].acc * dt;
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
class PtclH4{
public:
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc0;
    PS::F64vec acc1;
    PS::F64 mass;
    PS::F64 time;
    PS::F64 dt;
    PS::F64 r_out;
#ifdef HARD_DEBUG
    PS::F64vec acc2; // for debug
    PS::F64vec acc3; // for debug
#endif
};

class PtclForce{
public:
    PS::F64vec acc0; //
    PS::F64vec acc1; //
    void clear(){
        acc0 = acc1 = 0.0;
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
class HermiteIntegrator{
private:
    PS::ReallocatableArray<PtclH4> pred_;
    PS::ReallocatableArray<PtclForce> force_;
    PS::ReallocatableArray<PS::S32> adr_sorted_;
    PS::ReallocatableArray<PS::F64> time_next_;

    PS::F64 a0_offset_sq_; /// time step parameter
    PS::F64 time_sys_;     /// current time
    PS::S32 n_act_;        /// active particle number

    //  import parameter
    PS::F64 dt_limit_hard_; // time step limit     
    PS::F64 eta_s_;         // time step parameter 
    PS::F64 r_in_;          // force parameter     
    PS::F64 eps_sq_;        // softening parameter 
    
    // import data
    Tptcl * ptcl_org_;
    PS::S32 n_ptcl_;
    ReallocatableArray<PS::S32> *Jlist_;

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


    template<class Tptcl>
    void CalcBlockDt2ndAct(Tptcl ptcl[],
                           const PtclForce force[],
                           const PS::S32 adr_array[],
                           const PS::S32 n_act, 
                           const PS::F64 eta,
                           const PS::F64 dt_limit,
                           const PS::F64 a0_offset_sq){
        for(PS::S32 i=0; i<n_act; i++){
            const PS::S32 adr = adr_array[i];
            const PS::F64vec a0 = force[adr].acc0;
            const PS::F64vec a1 = force[adr].acc1;
        
            const PS::F64 dt_ref = CalcDt2nd(a0, a1, eta, a0_offset_sq);
            PS::F64 dt = dt_limit;
            while(dt > dt_ref) dt *= 0.5;
            ptcl[adr].dt = dt;
        }
    }

    template<class Tptcl>
    void CalcAcc0Acc1Act(PtclForce force[],
                         const Tptcl ptcl[],
                         const PS::S32 n_act,
                         const PS::S32 Ilist[],
                         const ReallocatableArray<PS::S32> Jlist[],
                         const PS::F64 rin,
                         const PS::F64 eps2) {
        // PS::ReallocatableArray< std::pair<PS::S32, PS::S32> > & merge_pair ){
        // active particles
        for(PS::S32 i=0; i<n_act; i++){
            PS::S32 iadr = Ilist[i];
            force[iadr].acc0 = force[iadr].acc1 = 0.0;
            // all
            for(PS::S32 j=0; j<Jlist[i].size(); j++){
                PS::S32 jadr = Jlist[i][j];
#ifdef HARD_DEBUG
                assert(iadr!=jadr);
#endif
                PS::F64 r2 = 0.0;
                PS::F64 rout = std::max(ptcl[i].r_out, ptcl[j].r_out);
                CalcAcc0Acc1R2Cutoff(ptcl[adr].pos, ptcl[adr].vel,
                                     force[adr].acc0, force[adr].acc1, r2,
                                     ptcl[j].pos, ptcl[j].vel, ptcl[j].mass,
                                     eps2, rout, rin);
                // if(r2 < ((ptcl[adr].r_merge + ptcl[j].r_merge)*(ptcl[adr].r_merge + ptcl[j].r_merge)) && ptcl[j].mass > 0.0){
                //     merge_pair.push_back( std::make_pair(adr, j) );
                // }
            }
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
                         PS::S32 group_list[],
                         PS::S32 & group_n,
                         const PS::S32 adr_group[],
                         const PS::S32 adr_group_map[],
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
        group_n = 0;
        for(n_act=1; n_act<n_tot; n_act++){
            PS::S32 adr = adr_sorted[n_act]
            if(time_ref < time_next[adr]) {
                break;
            }
            if(adr_group[adr]>=0) {
                group_list[group_n++] = adr_group_map[adr];
            }
        }
    }

    template<class Tptcl>
    void PredictAll(PtclH4 pred[],
                    const Tptcl ptcl[],
                    const PS::S32 n_tot,
                    const PS::F64 time_next){
        static const PS::F64 inv3 = 1.0 / 3.0;
        for(PS::S32 i=0; i<n_tot; i++){
            const PS::F64 dt = time_next - pred[i].time;
            pred[i].pos = ptcl[i].pos + dt*(ptcl[i].vel  + 0.5*dt*(pred[i].acc0 + inv3*dt*pred[i].acc1));
            pred[i].vel = ptcl[i].vel + dt*(pred[i].acc0 + 0.5*dt*pred[i].acc1);
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

    template<class Tptcl>
    void CorrectAndCalcDt4thAct(Tptcl ptcl[],
                                PtclH4 pred[],
                                const PtclForce force[],
                                const PS::S32 adr_sorted[], 
                                const PS::S32 n_act,
                                const PS::F64 dt_limit,
                                const PS::F64 a0_offset_sq,
                                const PS::F64 eta){
        static const PS::F64 inv3 = 1.0 / 3.0;
        for(PS::S32 i=0; i<n_act; i++){
            const PS::S32 adr = adr_sorted[i];
            Tptcl*     pti = &ptcl[adr];
            PtclH4*    pri = &pred[adr];
            const PtclForce* fpi = &force[adr];

            const PS::F64 dt = pri->dt;
            const PS::F64 h = 0.5 * dt;
            const PS::F64 hinv = 2.0 / dt;
            const PS::F64vec A0p = (fpi->acc0 + pri->acc0);
            const PS::F64vec A0m = (fpi->acc0 - pri->acc0);
            const PS::F64vec A1p = (fpi->acc1 + pri->acc1)*h;
            const PS::F64vec A1m = (fpi->acc1 - pri->acc1)*h;

            pti->vel = pri->vel + h*( A0p - inv3*A1m );
            pti->pos = pri->pos + h*( (pri->vel + pti->vel) + h*(-inv3*A0m));

            pri->acc0 = fpi->acc0;
            pri->acc1 = fpi->acc1;
            pri->time += dt;

            const PS::F64vec acc3 = (1.5*hinv*hinv*hinv) * (A1p - A0m);
            const PS::F64vec acc2 = (0.5*hinv*hinv) * A1m + h*acc3;
            const PS::F64 dt_ref = CalcDt4th(pri->acc0, pri->acc1, acc2, acc3, eta, a0_offset_sq);

            const PS::F64 dt_old = pri->dt;
#ifdef HARD_DEBUG
            // for debug
            assert(dt_old != 0.0);
            pri->acc2 = acc2;
            pri->acc3 = acc3;
#endif
            pri->dt = dt_limit;
            while(pri->dt > dt_ref) pri->dt *= 0.5;
            pri->dt = dt_old*2 < pri->dt ?  dt_old*2 : pri->dt;
        }
    }


    template<class Teng>
    void CalcEnergyHard(const Tptcl ptcl[], const PS::S32 n_tot, Teng & eng, 
                        const PS::F64 r_in,  const PS::F64 eps_sq = 0.0){
        eng.kin = eng.pot = eng.tot = 0.0;
        for(PS::S32 i=0; i<n_tot; i++){
            eng.kin += 0.5 * ptcl[i].mass * ptcl[i].vel * ptcl[i].vel;

            for(PS::S32 j=i+1; j<n_tot; j++){
                PS::F64 r_out = std::max(ptcl[i].r_out,ptcl[j].r_out);
                PS::F64vec rij = ptcl[i].pos - ptcl[j].pos;
                PS::F64 dr = sqrt(rij*rij + eps_sq);
                eng.pot -= ptcl[i].mass/dr*(1.0 - CalcW(dr/r_out, r_in/r_out));
            }
        }
        eng.tot = eng.kin + eng.pot;
    }

public:

    HermiteIntegrator(const PS::S32 n) {
#ifdef HARD_DEBUG
        assert(pred_.size()==0);
        assert(force_.size()==0);
        assert(adr_sorted_.size()==0);
        assert(time_next_.size()==0);
#endif
        resizeArray(n);
    }

    void resizeArray(const PS::S32 n) {
        pred_.resizeNoInitialize(n);
        force_.resizeNoInitialize(n);
        adr_sorted_.resizeNoInitialize(n);
        time_next_.resizeNoInitialize(n);
    }

    void setPtclList(Tptcl * ptcl_org, 
                     const PS::S32 n_ptcl, 
                     ReallocatableArray<PS::S32> Jlist[]) {
        ptcl_org_ = ptcl_org;
        n_ptcl_   = n_ptcl;
        Jlist_    = Jlist;
        if(n_ptcl_>pred_.size()) {
            std::cerr<<"Error, array size not enough! "<<pred_.size()<<" < "<<n_ptcl_<<std::endl;
            abort();
        }
    }

    void setParams(const PS::F64 dt_limit_hard,  // time step limit
                   const PS::F64 eta_s,          // time step parameter
                   const PS::F64 r_in,           // force parameter
                   const PS::F64 eps_sq){        // softening parameter
        dt_limit_hard_= dt_limit_hard; 
        eta_s_        = eta_s;         
        r_in_         = r_in;          
        eps_sq_       = eps_sq;
    }

    void initialize(PS::F64 &dt_limit,
                    PS::S32 group_act_list[],
                    PS::S32 &group_act_n,
                    const PS::S32 adr_group[],
                    const PS::S32 adr_group_map[]) {
        PS::F64 mass_min = PS::LARGE_FLOAT;
        PS::F64 rout_min = PS::LARGE_FLOAT;
        for(PS::S32 i=0; i<n_ptcl_; i++){
            // pred[i].mass = ptcl[i].mass = ptcl_org[i].mass;
            // pred[i].pos  = ptcl[i].pos  = ptcl_org[i].pos;
            // pred[i].vel  = ptcl[i].vel  = ptcl_org[i].vel;
            // ptcl[i].setRMerge();
            // pred[i].r_merge = ptcl[i].r_merge;
            // ptcl[i].id = ptcl_org[i].id;
            pred_[i].r_out = ptcl_org_[i].r_out;
            adr_sorted_[i] = i;
            pred_[i].time = pred_[i].dt = 0.0;
            time_next_[i] = 0.0;
            if(mass_min > ptcl_org_[i].mass)  mass_min = ptcl_org_[i].mass;
            if(rout_min > ptcl_org_[i].r_out) rout_min = ptcl_org_[i].r_out;
        }

        time_sys_ = 0.0;
        a0_offset_sq_ = 0.1 * mass_min / (rout_min * rout_min);
        n_act_ = n_ptcl_;

        CalcAcc0Acc1Act(force.getPointer(), ptcl_org_, n_act_, adr_sorted.getPointer(), Jlist_, r_in_, eps_sq_);
    
        // store predicted force
        for(PS::S32 i=0; i<n_ptcl; i++){
            pred_[i].acc0 = force_[i].acc0;
            pred_[i].acc1 = force_[i].acc1;
        }

        dt_limit = calcDtLimit(time_sys_, dt_limit_hard_);
        CalcBlockDt2ndAct(pred_.getPointer(), force_.getPointer(), adr_sorted_.getPointer(), n_act_, eta_s_, dt_limit, a0_offset_sq_);

        for(PS::S32 i=0; i<n_ptcl; i++){
            time_next_[i] = pred_[i].time + pred_[i].dt;
        }
        SortAndSelectIp(adr_sorted_.getPointer(), time_next_.getPointer(), n_act_, group_act_list, group_act_n, adr_group, adr_group_map, n_ptcl_);
        
    }
    
    void CalcEnergyHard(Energy & eng) {
        CalcEnergyHard(ptcl_org_, n_ptcl_, eng, r_in_, eps_sq_);
    }

    PS::F64 getNextTime() {
        return time_next_[adr_sorted[0]];
    }

    void integrateOneStep(PS::F64 &dt_limit,
                          PS::S32 group_act_list[],
                          PS::S32 &group_act_n,
                          const PS::S32 adr_group[],
                          const PS::S32 adr_group_map[]) {
        // begin integration loop
        time_sys_ = time_next_[adr_sorted[0]];
        dt_limit = calcDtLimit(time_sys_, dt_limit_hard_);
        // pred::mass,pos,vel updated
        PredictAll(pred_.getPointer(), ptcl_org_, n_ptcl_, time_sys_);
        // force::acc0,acc1 updated
        CalcAcc0Acc1Act(force_.getPointer(), pred_.getPointer(), n_act_, adr_sorted_.getPointer(), Jlist_, r_in_, eps_sq_);
        // ptcl_org::pos,vel; pred::time,dt,acc0,acc1,acc2,acc3 updated
        CorrectAndCalcDt4thAct(ptcl_org_, pred_.getPointer(), force_.getPointer(), adr_sorted_.getPointer(), n_act_, dt_limit, a0_offset_sq_, eta_s_);

        for(PS::S32 i=0; i<n_act; i++){
            PS::S32 adr = adr_sorted_[i];
            time_next_[adr] = pred_[adr].time + pred_[adr].dt;
        }
        SortAndSelectIp(adr_sorted_.getPointer(), time_next_.getPointer(), n_act_, group_act_list, group_act_n, adr_group, adr_group_map, n_ptcl_);

        // for(PS::S32 i=0; i<n_ptcl; i++){
        //     ptcl_org[i].mass = ptcl[i].mass;
        //     ptcl_org[i].pos  = ptcl[i].pos;
        //     ptcl_org[i].vel  = ptcl[i].vel;
        // }
        return time_sys_;
    }

    const PS::S32 getNact() const{
        return N_act_;
    }

    const PS::S32* getActList() const{
        return adr_sorted_.getPointer();
    }

};

//few-body----------------------------------------------
template<class Tptcl, class ARC_par>
void Isolated_Multiple_integrator(Tptcl * ptcl_org,
                                  const PS::S32 n_ptcl,
                                  const PS::F64 time_end,
                                  const PS::F64 dt_limit_hard,
                                  const PS::F64 rout_single,
                                  const PS::F64 gamma,
                                  const PS::F64 m_average,
#ifdef ARC_ERROR
                                  PS::F64 &ARC_error_relative,
                                  PS::F64 &ARC_error,
                                  PS::S32 N_count[20],
#endif                         
                                  const ARC::chainpars<Tptcl, ARC_par> &ARC_control,
                                  ARC_int_pars &Int_pars) {
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
      
    ARC::chain<Tptcl,ARC_int_pars> c((std::size_t)n_ptcl,ARC_control);
    static thread_local PS::F64 time_sys = 0.0;

    c.addP(n_ptcl,ptcl_org);

#ifdef HARD_DEBUG
    for (PS::S32 i=0; i<n_ptcl; i++) 
        if (ptcl_org[i].r_out<Int_pars.rin) {
            std::cerr<<"Error, updated p["<<i<<"].rout ("<<ptcl_org[i].r_out<<") < r_in ("<<Int_pars.rin<<")"<<std::endl;
            abort();
        }
#endif

    c.link_int_par(Int_pars);
    c.init(time_sys);

#ifdef ARC_ERROR
    PS::F64 ARC_error_once = c.getPot()+c.getEkin();
    if(n_ptcl<=20) N_count[n_ptcl-1]++;
    else std::cerr<<"Large cluster formed, n="<<n_ptcl<<std::endl;
#endif
      
    PS::F64 dscoff=1.0;
    PS::F64 ds_up_limit = 0.25*calcDtLimit(time_sys, dt_limit_hard)/c.calc_dt_X(1.0);
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

        PS::F64 dsf=c.extrapolation_integration(ds_use,time_end);

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
#ifdef SAFETY_CHECK
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

template<class Tptcl>
class ARCIntegrator{
private:
    typedef ARC::chain<Tptcl,ARC_int_pars> ARChain;
    ReallocatableArray<ARChain> clist_;

    ARC::chainpars<Tptcl, ARC_par> *ARC_control_;
    ARC_int_pars *Int_pars_;

public:
    void initialize(PS::S32 group_list[];
                    ReallocatableArray<Tptcl> groups[],
                    const PS::S32 n_groups,
                    Tptcl* ptcl,
                    PS::S32 adr_cm[],
                    ReallocatableArray<PS::S32> pertlist[],
                    ARC::chainpars<Tptcl, ARC_par> &control,
                    ARC_int_pars &Int_pars) {
        for(int i=0; i<n_groups; i++) {
            PS::S32 icm = adr_cm[i];
            PS::S32 ig = group_list[i];
            PS::S32 ni = groups[ig].size();
            clist_.push_back(ARChain(ni,control));
            clist_[i].addP(ni,groups[ig].getPointer());
            clist_[i].link_int_par(Int_pars);
            for(int j=0; j<pertlist[icm].size(); j++) clist_[i].addPext(ptcl[pertlist[icm][j]]);
            clist_[i].init(0);
        }
        ARC_control_ = &control;
        Int_pars_ = &Int_pars;
    }

    void IntegrateOneStep(const PS::S32 ic,
                          const PS::F64 time_end,
                          const PS::F64 dt_limit) {
        ARChain* c = &clist_[ic];
        PS::F64 dscoff=1.0;
        PS::F64 ds_up_limit = 0.25*dt_limit/c->calc_dt_X(1.0);
        PS::F64 ds_use = c->calc_next_step_custom();
        
        if (ds_use>ds_up_limit) ds_use = ds_up_limit;

        // convergency check
        PS::S32 converge_count=0;
        PS::S32 error_count=0;
        bool modify_step_flag=false;
        bool final_flag=false;

        while(time_end-c->getTime()>ARC_control_->dterr) {
            PS::F64 dsf=c->extrapolation_integration(ds_use,time_end);
            if (dsf<0) {
                final_flag=true;
                converge_count++;
                if (converge_count>5&&time_end-c->getTime()>ARC_control_->dterr*100) {
                    std::cerr<<"Error: Time synchronization fails!\nStep size ds: "<<ds_use<<"\nEnding physical time: "<<time_end<<"\nTime difference: "<<time_end-c->getTime()<<"\nR_in: "<<Int_pars_->rin<<"\nR_out: "<<Int_pars_->rout<<"\n";
                    c->dump("ARC_dump.dat");
                    ARC_control_->dump("ARC_dump.par");
                    c->print(std::cerr);
                    abort();
                }
                else ds_use *= -dsf;
            }
            else if (dsf==0) {
                c->info->ErrMessage(std::cerr);
                error_count++;
                if(error_count>4) {
                    std::cerr<<"Error: Too much error appear!\nStep size ds: "<<ds_use<<"\nEnding physical time: "<<time_end<<"\nTime difference: "<<time_end-c->getTime()<<"\nR_in: "<<Int_pars_->rin<<"\nR_out: "<<Int_pars_->rout<<"\n";
                    c->dump("ARC_dump.dat");
                    ARC_control_->dump("ARC_dump.par");
                    c->print(std::cerr);
                    abort();
                }
                if (c->info->status==5) {
                    dscoff = 0.25;
                    ds_use *= dscoff;
                }
                else if (c->info->status==4) ds_use = std::min(dscoff*c->calc_next_step_custom(),ds_up_limit);
                else ds_use *= 0.1;
                modify_step_flag=true;
            }
            else  {
                if (final_flag) {
                    if (converge_count>6&&time_end-c->getTime()>ARC_control_->dterr*100) {
                        std::cerr<<"Error: Time synchronization fails!\nStep size ds: "<<ds_use<<"\nEnding physical time: "<<time_end<<"\nTime difference: "<<time_end-c->getTime()<<"\nR_in: "<<Int_pars_->rin<<"\nR_out: "<<Int_pars_->rout<<"\n";
                        c->dump("ARC_dump.dat");
                        ARC_control_->dump("ARC_dump.par");
                        c->print(std::cerr);
                        abort();
                    }
                    converge_count++;
                }
                else if (n_ptcl>2||(modify_step_flag&&error_count==0)) {
                    ds_use = std::min(dscoff*c->calc_next_step_custom(),ds_up_limit);
                    modify_step_flag=false;
                }
                // reducing error counter if integration success, this is to avoid the significant change of step may cause some issue
                if(error_count>0) error_count--;
            }
        }
    }

    void integrateOneStepList(PS::S32 act_list[];
                              PS::S32 n_act;
                              const PS::F64 time_end,
                              const PS::F64 dt_limit) {
        for(int i=0; i<n_act; i++) {
            IntegrateOneStep(act_list[i], time_end, dt_limit);
        }
    }

    void updateCM(Tptcl ptcl[],
                  PS::S32 act_list[];
                  PS::S32 adr_cm[];
                  PS::S32 n_act;) {
        for(int i=0; i<n_act; i++) {
            PS::S32 iact = act_list[i];
            PS::S32 icm  = adr_cm[iact];
            clist_[iact].cm.pos = ptcl[icm].pos;
            clist_[iact].cm.vel = ptcl[icm].vel;
#ifdef HARD_DEBUG
            assert(clist_[iact].cm.mass==ptcl[icm].mass);
#endif
        }
    }

    void CMshiftreverse() {
        for(int i=0; i<clist_.size(); i++) {
            clist_[i].center_shift_inverse();
        }
    }

#ifdef ARC_ERROR
    void ErrorRecord(PS::F64 &ARC_error) {
        ARC_error = 0.0
        for(int i=0; i<clist_.size(); i++) {
            ARC_error += clist_[i].getPot()+clist_[i].getEkin();
        }
    }
#endif                         

};
