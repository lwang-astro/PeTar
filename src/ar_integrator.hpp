#pragma once

#include "hermite_ptcl.hpp"
#include "ar_perturber.hpp"
#include "ar_interaction.hpp"
#include "AR/symplectic_integrator.h"

class ARIntegrator{
private:
    typedef AR::SymplecticIntegrator<Ptcl, PtclH4, ARPerturber, ARInteraction, ARInformation> ARSym;
    SymplecticManager<ARInteraction>* manager_;
    PS::ReallocatableArray<ARSym> groups_;
    PS::ReallocatableArray<bool> group_mask_map_; 
    PS::ReallocatableArray<PS::S32> group_mask_list_;
    
public:
    //! Reserve memory
    void reserveMem(const PS::S32 _n) {
#ifdef HARD_DEBUG
        assert(_n<ARRAY_ALLOW_LIMIT);
#endif        
        groups_.reserve(_n);
        bininfo.reserve(_n);
        group_mask_map_.reserve(_n);
        group_mask_list_.reserve(_n);
    }

    //! Add group of particles to AR class
    /*! Add one group in AR, use the suppressed group first, if no, add new group
      @param[in] _ptcl: particle data array
      @param[in] _ptcl_list: particle index array for _ptcl (if NULL, read continuelly from 1 to _n_ptcl
      @param[in] _n_ptcl: number of particles
      @param[in] _ptcl_soft_pert: soft perturbation artifical particles (if NULL, keep pert_par same as before)
      @param[in] _n_split: split number for artifical particles
      @param[in] _ptcl_pert: perturber particle array, notice the first _n_group are c.m. which has consistent order of AR groups
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
                        PtclH4* _ptcl_pert = NULL,
                        const PS::S32* _ptcl_pert_list = NULL,
                        const PS::S32 _n_pert = 0) {
        // set current group offset
        const PS::S32 ngroup = groups_.size();
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

            group_mask_map_.push_back(false);
            groups_.increaseSize(1);
            bininfo.increaseSize(1);
        }

        ARSym& group = groups_[igroup];

        // Soft perturbation
        if(_ptcl_soft_pert)
            group.perturber.soft_pert.fit(_ptcl_soft_pert, bininfo[igroup], manager_->interaction.r_crit, manager_->interaction.n_split);

        // allocate memory
        group.particles.setMode(AR::LisdMode::copy);
        group.particles.reserveMem(_n_ptcl);
        group.reserveForceMem();
        
        // Add members to AR 
        if(_ptcl_list) { // use index from list 
            for(PS::S32 i=0; i<_n_ptcl; i++) {
                group.addMemberAndAddress(_ptcl[_ptcl_list[i]]);
            }
        }
        else { // read one by one
            for(PS::S32 i=0; i<_n_ptcl; i++) {
                group.addMemberAndAddress(_ptcl[i]);
            }
        }
        
        // calculate the c.m.
        group.particles.calcCenterOfMass();
        // shift to c.m. frame
        group.particles.shiftToCenterOfMassFrame();

        // Add perturber
        group.perturber.neighbor.setMode(AR::ListMode::local);
        for(PS::S32 i=0; i<_n_pert; i++) {
            const PS::S32  k = _ptcl_pert_list[i];
            group.perturber.neighbor.addMember(&_ptcl_pert[k]);
        }

        return igroup;
    }

    //! Clear one group
    /*! Clear one group, remove members, relink to c.m., member number to 0
     */
    void clearOneGroup(const PS::S32 _igroup) {
#ifdef HARD_DEBUG
        assert(_igroup<groups_.size());
#endif
        if(!group_mask_map_[_igroup]) {
            groups_[_igroup].clear();
            bininfo[_igroup].tstep = -1.0;
            group_mask_list_.push_back(_igroup);
            group_mask_map_[_igroup]=true;
        }
    }

    //! Update perturber list
    /*! Update perturber list for group i
      @param[in] _i_group: group index for update perturber
      @param[in] _ptcl_pert: perturber particle array
      @param[in] _ptcl_pert_list: new perturber particle index
      @param[in] _n_pert: number of perturbers
     */
    void updatePertOneGroup(const PS::S32 _igroup,
                            PtclH4* _ptcl_pert,
                            const PS::S32* _ptcl_pert_list,
                            const PS::S32 _n_pert) {
#ifdef HARD_DEBUG
        assert(!group_mask_map_[_igroup]);
#endif
        // Add one for c.m.
        ARSym& group = groups_[_igroup];
        group.perturber.neighbor.resizeNoInitialize(0);
        for (PS::S32 i=0; i<_n_pert; i++) {
            PS::S32 adr = _ptcl_pert_list[i];
            group.perturber.neighbor.addMember(&_ptcl_pert[adr]);
        }
    }

#ifdef HARD_DEBUG
    //! check perturber list
    bool checkPert() {
        for (PS::S32 i=0; i<groups_.size(); i++) {
            if(!group_mask_map_[i]){
                ARSym& group = groups_[i];
                const PS::S32 n_member = group.particles.getSize();
                const PS::S32 n_pert = group.perturber.neighbor.getSize();
                for (PS::S32 j=0; j<n_pert; j++) {
                    for (PS::S32 k=0; k<n_member; k++) {
#ifdef HARD_DEBUG_DUMP
                        if (group.particles[k].id==group.perturber.neighbor[j].id) {
                            std::cerr<<"Error: member id == perturber id\n";
                            return true;
                        }
#else
                        assert(group.particles[k].id!=group.perturber.neighbor[j].id);
#endif
                    }
#ifdef HARD_DEBUG_DUMP
                    if (groups_[i].id==pert_[j+i_pert_off]->id) {
                        std::cerr<<"Error: groups_[i].id!=pert_[j+i_pert_off]->id\n";
                        return true;
                    }
#else
                    assert(groups_[i].id!=pert_[j+i_pert_off]->id);
#endif
                }
#ifdef HARD_DEBUG_DUMP
                if (groups_[i].id!=pert_[i_pert_off]->id) {
                    std::cerr<<"Error: groups_[i].id!=pert_[j+i_pert_off]->id\n";
                    return true;
                }
#else
                assert(groups_[i].id==pert_[i_pert_off]->id);
#endif
            }
        }
        return false;
    }
#endif


    //! Set initial slowdown parameter for one unpert group 
    /*! 
      @param[in] _i_group: group index
      @param[in] _tend: ending physical time for integration
      @param[in] _sdfactor: slowdown criterion factor
    */
    void initialOneSlowDownUnPert(const PS::S32 _i_group, const PS::F64 _tend, const PS::F64 _sdfactor) {
        // isolated case
        if (bininfo[_i_group].semi>0&&bininfo[_i_group].stable_factor>=0) {   
            // estimate inner acceleration diference at apo-center
            PS::F64 finner = bininfo[_i_group].semi*(1.0+bininfo[_i_group].ecc);
            finner = (bininfo[_i_group].m1+bininfo[_i_group].m2)/(finner*finner);
            PS::F64 finnersq = finner*finner;
            // get apo-center position and velocity
            TpAR p[2];
            OrbParam2PosVel(p[0].pos, p[1].pos, p[0].vel, p[1].vel, bininfo[_i_group].m1, bininfo[_i_group].m2, bininfo[_i_group].semi, bininfo[_i_group].ecc, bininfo[_i_group].inc, bininfo[_i_group].OMG, bininfo[_i_group].omg, PI);
            p[0].mass = bininfo[_i_group].m1;
            p[1].mass = bininfo[_i_group].m2;
#ifdef SOFT_PERT
#ifndef TIDAL_TENSOR
            p[0].status = 0;
            p[1].status = 1;
#endif
#endif
            //center_of_mass_correction(*(TpAR*)&groups_[_i_group], p, 2);
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
            groups_[_i_group].slowdown.setSlowDownPars(bininfo[_i_group].peri, _sdfactor/((p[0].mass+p[1].mass)*Ptcl::mean_mass_inv), std::max(1.0, _tend/bininfo[_i_group].peri));
            groups_[_i_group].slowdown.initialFRatioSqRange(fpertsq/finnersq);
            //groups_[_i_group].slowdown.updatekappa(_tend, 1.0, _tp_factor,-1);
            groups_[_i_group].slowdown.updateKappaMin();
            //groups_[_i_group].slowdown.updateKappa();
        }
    }

    //void initialOneSlowDown(const PS::S32 _i_group, const PS::F64 _tend, const PS::F64 _mpert, const PS::F64 _sdfactor, const PS::F64 _tp_factor) {
    void initialOneSlowDown(const PS::S32 _i_group, const PS::F64 _dt_limit_hard, const PS::F64 _sdfactor, bool _set_one_flag=false) {
        if (bininfo[_i_group].semi>0&&bininfo[_i_group].stable_factor>=0) {
            groups_[_i_group].slowdown.setSlowDownPars(bininfo[_i_group].peri, _sdfactor/((bininfo[_i_group].m1+bininfo[_i_group].m2)*Ptcl::mean_mass_inv), std::max(1.0, _dt_limit_hard/bininfo[_i_group].peri));
            //groups_[_i_group].slowdown.updatekappa(_tend, groups_[_i_group].mass/_mpert, _tp_factor,-1);
            if (_set_one_flag) groups_[_i_group].slowdown.setKappa(1.0);
            else groups_[_i_group].slowdown.updateKappaMin();
            //groups_[_i_group].slowdown.updateKappa();
        }
    }

    //! Update slow down factor for one AR
    /*! Update slowdown for one AR
      @param[in] _igroup: index of AR
      @param[in] _md_factor: slowdown modification limit factor (negative suppress the limit)
     */
    void updateOneSlowDown(const size_t _igroup, const PS::F64 _md_factor) {
#ifdef HARD_DEBUG
        assert(!group_mask_map_[_igroup]);
#endif
        //PS::F64 tp_factor = std::max(1e-4,_dt/_dt_limit);
        //PS::F64 tp_factor = std::max(1e-4,_dt/_dt_limit);
        //PS::F64 tp_factor = _dt/_dt_limit;
        //std::cerr<<"i "<<_index<<" dt "<<_dt<<" fac "<<tp_factor<<std::endl;
        groups_[_igroup].slowdown.updateKappaMinPeriod(_md_factor);
        //groups_[_igroup].slowdown.updateKappa();
    }

    void adjustSlowDown(const PS::F64 dt) {
        for (PS::S32 i=0; i<groups_.size(); i++) {
            groups_[i].slowdown.adjustkappa(dt);
        }
    }

    void adjustSlowDownPeriod(const PS::F64 dt, PS::S32* np) {
        for (PS::S32 i=0; i<groups_.size(); i++) {
            np[i] = groups_[i].slowdown.adjustkappaPeriod(dt);
        }
    }
    
    // return fail_flag
    bool initialSys() {
        for (PS::S32 i=0; i<groups_.size(); i++) {
            const PS::S32 ipert = pert_disp_[i];
            groups_[i].initSys(0.0, *AR_control_, &(par_list_.back()), &pert_[ipert], &pforce_[ipert], pert_n_[i]);
#ifdef AR_WARN
            if(groups_[i].info!=NULL) {
                groups_[i].info->ErrMessage(std::cerr);
                return true;
            }
#endif
        }
        return false;
    }

    //! Initial one group chain integration parameters
    /*! Initial one group chain integration parameters
      @param[in] _i_group: group to initialize
      @param[in] _time_sys: time to initialize
      \return fail_flag
     */
    bool initialOneSys(const PS::S32 _i_group, const PS::F64 _time_sys) {
#ifdef HARD_DEBUG
#ifdef HARD_DEBUG_DUMP
        if (group_mask_map_[_i_group]) return true;
#else
        assert(!group_mask_map_[_i_group]);
#endif
#endif
        const PS::S32 ipert = pert_disp_[_i_group];
        groups_[_i_group].initSys(_time_sys, *AR_control_, &(par_list_[_i_group]),&pert_[ipert], &pforce_[ipert], pert_n_[_i_group]);
#ifdef AR_WARN
        if(groups_[_i_group].info!=NULL) {
            groups_[_i_group].info->ErrMessage(std::cerr);
            return true;
        }
#endif
        return false;
    }

    //! Initial one group chain member and c.m.
    void initialOneChain(const PS::S32 _i_group) {
#ifdef HARD_DEBUG
        assert(!group_mask_map_[_i_group]);
#endif
        groups_[_i_group].initChain();
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

        groups_[ic].dump(fp);
        AR_control_->dump(fp);
        bininfo[ic].dump(fp);

        std::fclose(fp);
        //groups_[ic].print(std::cerr);
    }

    PS::S64 integrateOneStepSymTwo(const PS::S32 ic, const PS::F64 time_end, const PS::S32 kp) {
#ifdef HARD_DEBUG
        assert(!group_mask_map_[ic]);
#endif
        ARSym* c = &groups_[ic];
        AR_pert_pars* par = &par_list_[ic];
        PS::F64 ds_use=bininfo[ic].tstep;
        const PS::S32 ipert = pert_disp_[ic];
        PS::F64 timetable[8]; // Notice, assuming sym order is -6
#ifdef AR_OPT_SYM2
        const PS::F64 m1=c->getP(0).getMass();
        const PS::F64 m2=c->getP(1).getMass();
        const PS::F64 m2_mt = m2/(m1+m2);
        const PS::F64 m1_m2_1 = -m1/m2-1.0;
#endif
        const PS::S32 np=8*kp;
        for (PS::S32 i=0; i<np; i++) {
#ifdef AR_OPT_SYM2
            c->Symplectic_integration_two(ds_use, *AR_control_, timetable, m2_mt, m1_m2_1, par, &pert_[ipert], &pforce_[ipert], pert_n_[ic]);
#else 
            c->Symplectic_integration(ds_use, *AR_control_, timetable, par, &pert_[ipert], &pforce_[ipert], pert_n_[ic]);
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
#ifdef AR_WARN       
        if((c->getTime()-time_end)/time_end>1e-6) {
            std::cerr<<"Warning! time not synchronized! t(chain)="<<c->getTime()<<" t="<<time_end<<" diff="<<(c->getTime()-time_end)/time_end<<std::endl;
        }
#endif
        return np;
    }

#ifdef AR_SYM
    //! integration arc with symplectic method
    /*! 
      @param[in] _igroup: the AR group id
      @param[in] _time_end: finishing time
      @param[in] _dt_limit: physical time step limit
      \return stepcount; negative means fail case
     */
    PS::S64 integrateOneStepSym(const PS::S32 _igroup,
                                const PS::F64 _time_end,
                                const PS::F64 _dt_limit) {
#ifdef HARD_DEBUG
        assert(!group_mask_map_[_igroup]);
        assert(_time_end>groups_[_igroup].getTime());
#endif
        ARSym* c = &groups_[_igroup];
        AR_pert_pars* par = &par_list_[_igroup];
        //PS::F64 ds_up_limit = 0.25*_dt_limit/c->calc_dt_X(1.0,*AR_control_);
        //PS::F64 ds_use = 2.0*bininfo[_igroup].tstep*std::abs(c->getPt());
        PS::F64 ds_use=bininfo[_igroup].tstep;
        // in case dt is much less then period, reduce step
        //if (_dt_limit<bininfo[_igroup].peri*c->slowdown.getkappa()) {
        //    ds_use *= _dt_limit/(bininfo[_igroup].peri*c->slowdown.getkappa());
        //}
        if(c->slowdown.isUsed()) {
            PS::F64 korg=c->slowdown.getkappaorg();
            // in strong perturbed case, avoid too large step size
            if(korg<1.0) ds_use *= 1.0/8.0*std::pow(korg,1.0/6.0);
        }
        //PS::F64 ds_use = c->calc_next_step_custom(*AR_control_,par);
        //if (ds_use>ds_up_limit) ds_use = ds_up_limit;

        const PS::S32 ipert = pert_disp_[_igroup];
        // for standard case, not fix step
        PS::S32 fix_step_flag = 0;
        // for two-body or few-body with closed orbit
        if (c->getN()==2||(c->getN()>2&&bininfo[_igroup].semi>0)) fix_step_flag = 1;
        // for high-eccentric binary, it is better to fix step to avoid big step drop, for hyperbolic, fix step is risky
        //if(c->getN()==2&&bininfo[_igroup].ecc>0.99&&bininfo[_igroup].ecc<1.0) {
        //    fix_step_flag = 2;
        //    PS::F64 korg=c->slowdown.getkappaorg();
        //    if(korg<1.0) ds_use *= 1.0/8.0*std::pow(korg,1.0/6.0);
        //}

        PS::S64 stepcount = c->Symplectic_integration_tsyn(ds_use, *AR_control_, _time_end, par, &pert_[ipert], &pforce_[ipert], pert_n_[_igroup],fix_step_flag, step_count_limit);

#ifdef HARD_DEBUG
        assert(_time_end+0.5*_dt_limit>c->getTime());
#endif

#ifdef AR_WARN
        if(c->info!=NULL) {
            c->info->ErrMessage(std::cerr);
        }
#endif
        
#ifdef AR_DEBUG_DUMP
        if(stepcount<0) {
            dump("AR_dump.dat",_igroup,_time_end,ds_use);
            std::cerr<<"Igroup = "<<_igroup<<" N = "<<c->getN()<<" Np = "<<pert_n_[_igroup]<<" stepcount = "<<stepcount<<std::endl;
        }
#endif
        
        return stepcount;
    }
#else

    //! integration arc with extrapolation method
    /*! 
      @param[in] _igroup: the AR group id
      @param[in] _time_end: finishing time
      @param[in] _dt_limit: physical time step limit
      \return stepcount; negative means fail case
     */
    PS::S64 integrateOneStepExt(const PS::S32 _igroup,
                                const PS::F64 _time_end,
                                const PS::F64 _dt_limit) {
#ifdef HARD_DEBUG
        assert(!group_mask_map_[_igroup]);
#endif
        ARSym* c = &groups_[_igroup];
        AR_pert_pars* par = &par_list_[_igroup];
        PS::F64 dscoff=1.0;
        PS::F64 ds_up_limit = 0.25*_dt_limit/c->calc_dt_X(1.0,*AR_control_);
        PS::F64 ds_use = c->calc_next_step_custom(*AR_control_,par);
        //PS::F64 ds_use = 0.5*bininfo[_igroup].tstep*std::abs(c->GetPt());
        
        if (ds_use>ds_up_limit) ds_use = ds_up_limit;

        PS::S64 nstep=0;

        // convergency check
        PS::S32 converge_count=0;
        PS::S32 error_count=0;
        bool modify_step_flag=false;
        bool final_flag=false;

        while(_time_end-c->getTime()>AR_control_->dterr*c->getTime()) {
            const PS::S32 ipert = pert_disp_[_igroup];
            PS::F64 dsf=c->extrapolation_integration(ds_use, *AR_control_, _time_end, par, &pert_[ipert], &pforce_[ipert], pert_n_[_igroup]);
            if (dsf<0) {
                final_flag=true;
                converge_count++;
                if (converge_count>10&&_time_end-c->getTime()>AR_control_->dterr*100) {
                    std::cerr<<"Error: Time synchronization fails!\nStep size ds: "<<ds_use<<"\nEnding physical time: "<<_time_end<<"\nTime difference: "<<_time_end-c->getTime()<<"\nR_in: "<<Int_pars_->rin<<"\nR_out: "<<Int_pars_->rout<<"\n";
                    dump("AR_dump.dat",_igroup,_time_end,ds_use);
                    return -1;
                }
                else ds_use *= -dsf;
            }
            else if (dsf==0) {
                c->info->ErrMessage(std::cerr);
                error_count++;
                if(error_count>4) {
                    std::cerr<<"Error: Too much error appear!\nStep size ds: "<<ds_use<<"\nEnding physical time: "<<_time_end<<"\nTime difference: "<<_time_end-c->getTime()<<"\nR_in: "<<Int_pars_->rin<<"\nR_out: "<<Int_pars_->rout<<"\n";
                    dump("AR_dump.dat",_igroup,_time_end,ds_use);
                    return -1;
                }
                if (c->info->status==5) {
                    dscoff = 0.25;
                    ds_use *= dscoff;
                }
                else if (c->info->status==4) ds_use = std::min(dscoff*c->calc_next_step_custom(*AR_control_, par),ds_up_limit);
                else ds_use *= 0.1;
                modify_step_flag=true;
            }
            else  {
                if (final_flag) {
                    if (converge_count>10&&_time_end-c->getTime()>AR_control_->dterr*100) {
                        std::cerr<<"Error: Time synchronization fails!\nStep size ds: "<<ds_use<<"\nEnding physical time: "<<_time_end<<"\nTime difference: "<<_time_end-c->getTime()<<"\nR_in: "<<Int_pars_->rin<<"\nR_out: "<<Int_pars_->rout<<"\n";
                        dump("AR_dump.dat",_igroup,_time_end,ds_use);
                        return -1;
                    }
                    converge_count++;
                }
                else if (modify_step_flag&&error_count==0) {
                    ds_use = std::min(dscoff*c->calc_next_step_custom(*AR_control_, par),ds_up_limit);
                    modify_step_flag=false;
                }
                // reducing error counter if integration success, this is to avoid the significant change of step may cause some issue
                if(error_count>0) error_count--;
            }
        }
#ifdef AR_PROFILE
        nstep = c->profile.itercount;
#endif

        return nstep;
    }
#endif

    //! Integrate active AR groups
    /* @param[in] _act_list: active AR group index list
       @param[in] _n_act: number of active groups
       @param[in] _time_end: end of physical integration time
       @param[in] _dt_limit: physical time step upper limit
       \return stepcount; if negative, error happen
     */
    PS::S64 integrateOneStepList(PS::S32 _act_list[],
                                 PS::S32 _n_act,
                                 const PS::F64 _time_end,
                                 const PS::F64 _dt_limit) {
        PS::S64 nstep = 0;
        for(PS::S32 i=0; i<_n_act; i++) {
            if(getMask(i)) continue;
            PS::S64 nstep_i;
#ifdef AR_SYM
            nstep_i = integrateOneStepSym(_act_list[i], _time_end, _dt_limit);
#else
            nstep_i = integrateOneStepExt(_act_list[i], _time_end, _dt_limit);
#endif
            if (nstep_i<0) return nstep_i; // error case
            else nstep += nstep_i;
        }

        return nstep;
    }

    //! Integrate active AR groups
    /* @param[in] _time_end: end of physical integration time
       @param[in] _dt_limit: physical time step upper limit
       \return stepcount; if negative, error happen
     */
    PS::S64 integrateOneStepList(const PS::F64 _time_end,
                                 const PS::F64 _dt_limit) {
        PS::S64 nstep = 0;
        for(PS::S32 i=0; i<groups_.size(); i++) {
            if(getMask(i)) continue;
            PS::S64 nstep_i;
#ifdef AR_SYM
            nstep_i = integrateOneStepSym(i, _time_end, _dt_limit);
#else
            nstep_i = integrateOneStepExt(i, _time_end, _dt_limit);
#endif
            if (nstep_i<0) return nstep_i; // error case
            else nstep += nstep_i;
        }
        return nstep;
    }

    //! Update the c.m. data from the original particle data
    /*!
      @param[in] _ptcl: original particle array
      @param[in] _ptcl_list: particle index list need to be updated, assume _ptcl and AR group have consistent index
      @param[in] _n_ptcl: number of particles
     */
    template <class Tptcl>
    void updateCM(Tptcl* _ptcl,
                  PS::S32* _ptcl_list,
                  PS::S32 _n_ptcl) {
        for(PS::S32 i=0; i<_n_ptcl; i++) {
            PS::S32 k = _ptcl_list[i];
#ifdef HARD_DEBUG
            assert(k<groups_.size());
#endif
            groups_[k].pos =  _ptcl[k].pos;
            groups_[k].vel =  _ptcl[k].vel;
            groups_[k].mass = _ptcl[k].mass;
        }
    }

    //! Update the c.m. data from the original particle data
    /*!
      @param[in] _ptcl: original particle array
    */
    template <class Tptcl>
    void updateCM(Tptcl _ptcl[]) {
        for(PS::S32 i=0; i<groups_.size(); i++) {
            groups_[i].pos  = _ptcl[i].pos;
            groups_[i].vel  = _ptcl[i].vel;
            groups_[i].mass = _ptcl[i].mass;
        }
    }

    //! update rsearch of components based on c.m.
    /*!
      @param[in] _dt_tree: tree time step
      @param[in] _v_max: maximum velocity used to calcualte r_search
      \return rsearch maximum
     */
    PS::F64 updateRSearch(const PS::F64 _dt_tree, const PS::F64 _v_max) {
        PS::F64 dt_reduce_factor=1.0;
        for(PS::S32 i=0; i<groups_.size(); i++) {
            if(getMask(i)) continue;
            PS::F64 dt_reduce_fi = groups_[i].calcRSearch(_dt_tree, _v_max);
            dt_reduce_factor = std::max(dt_reduce_fi, dt_reduce_factor);
            TpAR** ipadr=groups_[i].getPAdr();
            for (PS::S32 k=0; k<groups_[i].getN(); k++)
                ipadr[k]->r_search = groups_[i].r_search;
        }
        return dt_reduce_factor;
    }

    //! Shift member ptcls to their c.m. frame
    /*! Shift all group members to their c.m. frame for AR integration
     */
    void shift2CM() {
        for(PS::S32 i=0; i<groups_.size(); i++) {
            groups_[i].center_shift();
        }
    }

    //! resolve all groups' member particles
    /*! shift compotent coordinates to original frame and save data to original particle address
     */
    void resolve() {
        for(PS::S32 i=0; i<groups_.size(); i++) {
            groups_[i].resolve();
        }
    }

    //! resolve a list of groups
    /*! shift compotent coordinates to original frame and save data to original particle address
      @param[in] _group_list: group list to resovle
      @param[in] _n_group: number of group to resolve
     */
    void resolve(const PS::S32 _group_list[], const PS::S32 _n_group) {
        for(PS::S32 i=0; i<_n_group; i++) {
            groups_[_group_list[i]].resolve();
        }
    }

    //! Check break condition
    /*! Check whether it is necessary to break the chain
      Inner distance criterion:
       1. get maximum distance (rmax) pair
       2. if rmax>r_crit, check whether the pair is go away or go close

      Perturbation criterion for closed orbit:
       2. if perturbation / inner force square > 2, break
      @param[out] _break_group_list: group index list to break
      @param[out] _break_isplit_list: index in chain list to split for corresponding groups
      @param[in] _r_crit2: distance (square) criterion to check whether need to break
      @param[in] _sd_factor: slowdown reference factor
      \return n_group_break: number of groups need to break
     */
    PS::S32 checkBreak(PS::S32* _break_group_list,
                       PS::S32* _break_isplit_list,
                       const PS::F64 _r_crit2) {
        PS::F64 r_max2;
        PS::S32 r_max_index;
        PS::S32 n_group_break=0;
        for (PS::S32 i=0; i<groups_.size(); i++) {
            if(getMask(i)) continue;
            // obtain maximum distance pair
            groups_[i].getRmaxIndex(r_max2, r_max_index);
            bool break_flag = false;
            bool out_flag;
            if(r_max2>_r_crit2) {
                // check whether outcome or income
                out_flag=groups_[i].getDirection(r_max_index);
                if(out_flag) break_flag = true;
            }
            // check strong perturbed case
            PS::F64 fp_sq = groups_[i].slowdown.getFratioSq();
            if (fp_sq>2 && bininfo[i].semi>0) break_flag = true;

            // check few-body inner perturbation
            PS::F64 frinsqi=1.0, frinsqj=1.0;
            if(groups_[i].getN()>2) {
                out_flag=groups_[i].getDirection(r_max_index);
                if(out_flag) {
                    groups_[i].getFratioInnerSq(frinsqi, frinsqj, r_max_index, *AR_control_, Int_pars_);
                    PS::F64 sd_factor = groups_[i].slowdown.getSDRef();
                    // if slowdown factor is large, break the group
                    if (std::min(frinsqi,frinsqj)<sd_factor*sd_factor) break_flag = true;
#ifdef ADJUST_GROUP_DEBUG
                    std::cout<<"Check inner fratio, i_group:"<<i<<" left:"<<frinsqi<<" right:"<<frinsqj<<std::endl;
#endif
                }
            }

            if (break_flag) {
#ifdef ADJUST_GROUP_DEBUG
                std::cout<<"Break group, group index: "<<i
                         <<" N_member: "<<groups_[i].getN()
                         <<" break index: "<<r_max_index
                         <<" Out case: "<<out_flag
                         <<" separation square: "<<r_max2
                         <<" r_crit square: "<<_r_crit2
                         <<" pert_ratio_square: "<<fp_sq
                         <<" inner pert ratio square, left: "<<frinsqi
                         <<" right: "<<frinsqj
                         <<std::endl;
#endif
                _break_group_list[n_group_break] = i;
                _break_isplit_list[n_group_break] = r_max_index;
                n_group_break++;
            }
        }
        return n_group_break;
    }

    //! Get group original address in chain order
    /*! 
      @param[out] _list particle address array to store the results
      @param[in] _igroup: group index to split
     */
    PS::S32 getPtclAdrChain(TpAR* _list[],
                            const PS::S32 _igroup){
        return groups_[_igroup].getPAdrChain(_list);
    }
    

    //! Get number of members in group i
    /*!
      @param[in] _igroup: group index
      \return number of members in group i
     */
    PS::S32 getGroupN(const PS::S32 _igroup) const {
        return groups_[_igroup].getN();
    }

    //! Get number of groups
    /*!
      \return number of groups
     */
    PS::S32 getNGroups() const {
        return groups_.size();
    }

    //! return integration mask
    /*!
      @param[in] _igroup: index of group
      \return true: suppressed for integration
     */
    bool getMask(const PS::S32 _igroup) const {
        return group_mask_map_[_igroup];
    }

    //! return CM particle (cast pointer to TpAR*)
    /*! 
      @param[in] _igroup: index of group
      \return TpAR* c.m. particle pointer of index i
     */
    TpAR* getCM(const PS::S32 _igroup) {
        return (TpAR*)&groups_[_igroup];
    }

    const TpAR* getGroupPtcl(const PS::S32 i) const {
#ifdef HARD_CHECK_ENERGY
        assert(i<groups_.size());
#endif        
        return &groups_[i].getP(0);
    }

    //! return the address array of member particle 
    /*!
      @param[in] _igroup: index of group
      \return TpAR** particle address array of members
     */
    TpAR** getGroupPtclAdr(const PS::S32 i) const {
#ifdef HARD_CHECK_ENERGY
        assert(i<groups_.size());
#endif        
        return groups_[i].getPAdr();
    }
    
    //! Get slow down factor
    /*! get slow down factor kappa
      @param[in] i: index of groups
      \return slowdown kappa
     */
    PS::F64 getSlowDown(const PS::S32 i) const{
        return groups_[i].slowdown.getkappa();
    }

    //! Get perturbation factor square
    /*! get current maximum perturbation/inner force square
      @param[in] i: index of groups
      \return fpert square
    */
    PS::F64 getFratioSq(const PS::S32 i) const{
        return groups_[i].slowdown.getFratioSq();
    }

    //! Get slow down original factor
    /*! get slow down original kappa
     */
    PS::F64 getSlowDownOrg(const PS::S32 i) const{
        return groups_[i].slowdown.getkappaorg();
    }

    //! Print slow down parameters
    /*! Print slow down parameters
      @param[in] _os: ofstream for printing
      @param[in] _i: Chain index
      @param[in] _precision: printed precision for one variable
      @param[in] _width: printing width for one variable
     */
    void printSlowDown(std::ostream& _os, const PS::S32 _i, const PS::S32 _precision=15, const PS::S32 _width=23) {
        groups_[_i].slowdown.print(_os,_precision,_width);
    }


//#ifdef AR_PROFILE
//    const PS::S64 getNsubstep() const{
//        PS::S64 Nsum = 0;
//        for (int i=0; i<groups_.size(); i++) 
//            Nsum += groups_[i].profile.itercount;
//        return Nsum;
//    }
//#endif

#ifdef AR_DEBUG_PRINT
    void data_dump(std::ostream& os, const PS::S32 i, const PS::F64 dt_limit) const{
        const AR_pert_pars* par = &par_list_[i];
        os<<std::setprecision(15)<<dt_limit<<" "
          <<groups_[i].getN()+pert_n_[i]<<" "
          <<par->rin<<" "
          <<par->rout<<" "
          <<par->rout<<" "
          <<par->rin*0.1<<" "
          <<dt_limit<<" "
          <<0.05<<" "
          <<std::sqrt(par->eps2)<<std::endl;
        for (PS::S32 j=0; j<groups_[i].getN(); j++) {
            os<<groups_[i].getP(j).mass<<" "
              <<groups_[i].getP(j).pos<<" "
              <<groups_[i].getP(j).vel<<std::endl;
        }
        for (PS::S32 j=1; j<pert_n_[i]; j++) {
            os<<pert_[pert_disp_[i]+j]->mass<<" "
              <<pert_[pert_disp_[i]+j]->pos<<" "
              <<pert_[pert_disp_[i]+j]->vel<<std::endl;
        }
    }

    //! AR info print
    /* @param[in] _n_group: current total number of groups already integrated
       @param[in] _n_group_in_cluster: number of groups in current cluster
       @param[in] _n_ptcl: number of real particles
       @param[in] _n_hint: number of Hint particles
       @param[in] _dt_limit: hard time step limit
       @param[in] _kp: kepler period number per step
     */
    bool info_print(std::ostream& os, const PS::S64 _n_group, const PS::S64 _n_group_in_cluster, const PS::S64 _n_ptcl, const PS::S64 _n_hint, const PS::F64 _dt_limit, const PS::S32 _kp, const PS::S32 _n_step_limit=10000) const{
        bool dump_flag=false;
        for (PS::S32 i=0; i<groups_.size(); i++) {
            if(getMask(i)) continue;
            os<<"AR_info: "
              <<" i_group_tot="<<_n_group+i
              <<" i_group="<<i
              <<" n_ptcl="<<_n_ptcl
              <<" n_groups="<<_n_group_in_cluster
              <<" n_hint="<<_n_hint
              <<" n_member="<<groups_[i].getN()
              <<" n_pert="<<pert_n_[i]
              <<" semi="<<bininfo[i].semi
              <<" ecc="<<bininfo[i].ecc
              <<" period="<<bininfo[i].peri
              <<" tstep="<<bininfo[i].tstep
              <<" sd="<<groups_[i].slowdown.getkappa();
            PS::S64 nstep = 0;
#ifdef AR_SYM
            if(_kp>0) nstep = _kp*8;
            else nstep = groups_[i].profile.stepcount[0];
#else
            nstep = groups_[i].profile.itercount;
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
        for(PS::S32 i=0; i<groups_.size(); i++) {
            if (sdflag) sd = 1.0/groups_[i].slowdown.getkappa();
            else sd = 1.0;
            energy.kin += sd*groups_[i].getEkin();
            energy.pot += sd*groups_[i].getPot();
            energy.tot += sd*groups_[i].getPt();
        }
    }
#endif                         

#ifdef HARD_DEBUG_PRINT
    template <class Tptcl>
    void writePtcl(FILE* _fout) const{
        for (PS::S32 i=0; i<groups_.size(); i++) {
            PS::F64vec pcm_pos = groups_[i].pos;
            PS::F64vec pcm_vel = groups_[i].vel;
            for (PS::S32 j=0; j<groups_[i].getN(); j++) {
                Tptcl pj = groups_[i].getP(j);
                pj.pos += pcm_pos;
                pj.vel += pcm_vel;
                pj.ParticleBase::writeAscii(_fout);
            }
        }
    }
#endif

};
