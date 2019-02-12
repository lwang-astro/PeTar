#pragma once

#include "hermite_ptcl.hpp"
#include "ar_perturber.hpp"
#include "ar_interaction.hpp"
#include "AR/symplectic_integrator.h"

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
    PS::F64 r_bin;     ///> for tidal tensor lscale
    
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
            par_list_[igroup].fit(_ptcl_soft_pert, bininfo[igroup], Int_pars_->r_bin, _n_split);

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
    bool checkPert() {
        for (PS::S32 i=0; i<clist_.size(); i++) {
            if(!group_mask_map_[i]){
                const PS::S32 n_member = clist_[i].getN();
                const PS::S32 i_pert_off = pert_disp_[i];
                for (PS::S32 j=1; j<pert_n_[i]; j++) {
                    for (PS::S32 k=0; k<n_member; k++) {
#ifdef HARD_DEBUG_DUMP
                        if (clist_[i].getP(k).id==pert_[j+i_pert_off]->id) {
                            std::cerr<<"Error: clist_[i].getP(k).id==pert_[j+i_pert_off]->id\n";
                            return true;
                        }
#else
                        assert(clist_[i].getP(k).id!=pert_[j+i_pert_off]->id);
#endif
                    }
#ifdef HARD_DEBUG_DUMP
                    if (clist_[i].id==pert_[j+i_pert_off]->id) {
                        std::cerr<<"Error: clist_[i].id!=pert_[j+i_pert_off]->id\n";
                        return true;
                    }
#else
                    assert(clist_[i].id!=pert_[j+i_pert_off]->id);
#endif
                }
#ifdef HARD_DEBUG_DUMP
                if (clist_[i].id!=pert_[i_pert_off]->id) {
                    std::cerr<<"Error: clist_[i].id!=pert_[j+i_pert_off]->id\n";
                    return true;
                }
#else
                assert(clist_[i].id==pert_[i_pert_off]->id);
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
            clist_[_i_group].slowdown.setSlowDownPars(bininfo[_i_group].peri, _sdfactor/((p[0].mass+p[1].mass)*Ptcl::mean_mass_inv), std::max(1.0, _tend/bininfo[_i_group].peri));
            clist_[_i_group].slowdown.initialFRatioSqRange(fpertsq/finnersq);
            //clist_[_i_group].slowdown.updatekappa(_tend, 1.0, _tp_factor,-1);
            clist_[_i_group].slowdown.updateKappaMin();
            //clist_[_i_group].slowdown.updateKappa();
        }
    }

    //void initialOneSlowDown(const PS::S32 _i_group, const PS::F64 _tend, const PS::F64 _mpert, const PS::F64 _sdfactor, const PS::F64 _tp_factor) {
    void initialOneSlowDown(const PS::S32 _i_group, const PS::F64 _dt_limit_hard, const PS::F64 _sdfactor, bool _set_one_flag=false) {
        if (bininfo[_i_group].semi>0&&bininfo[_i_group].stable_factor>=0) {
            clist_[_i_group].slowdown.setSlowDownPars(bininfo[_i_group].peri, _sdfactor/((bininfo[_i_group].m1+bininfo[_i_group].m2)*Ptcl::mean_mass_inv), std::max(1.0, _dt_limit_hard/bininfo[_i_group].peri));
            //clist_[_i_group].slowdown.updatekappa(_tend, clist_[_i_group].mass/_mpert, _tp_factor,-1);
            if (_set_one_flag) clist_[_i_group].slowdown.setKappa(1.0);
            else clist_[_i_group].slowdown.updateKappaMin();
            //clist_[_i_group].slowdown.updateKappa();
        }
    }

    //! Update slow down factor for one ARC
    /*! Update slowdown for one ARC
      @param[in] _igroup: index of ARC
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
        clist_[_igroup].slowdown.updateKappaMinPeriod(_md_factor);
        //clist_[_igroup].slowdown.updateKappa();
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
    
    // return fail_flag
    bool initialSys() {
        for (PS::S32 i=0; i<clist_.size(); i++) {
            const PS::S32 ipert = pert_disp_[i];
            clist_[i].initSys(0.0, *ARC_control_, &(par_list_.back()), &pert_[ipert], &pforce_[ipert], pert_n_[i]);
#ifdef ARC_WARN
            if(clist_[i].info!=NULL) {
                clist_[i].info->ErrMessage(std::cerr);
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
        clist_[_i_group].initSys(_time_sys, *ARC_control_, &(par_list_[_i_group]),&pert_[ipert], &pforce_[ipert], pert_n_[_i_group]);
#ifdef ARC_WARN
        if(clist_[_i_group].info!=NULL) {
            clist_[_i_group].info->ErrMessage(std::cerr);
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
    //! integration arc with symplectic method
    /*! 
      @param[in] _igroup: the ARC group id
      @param[in] _time_end: finishing time
      @param[in] _dt_limit: physical time step limit
      \return stepcount; negative means fail case
     */
    PS::S64 integrateOneStepSym(const PS::S32 _igroup,
                                const PS::F64 _time_end,
                                const PS::F64 _dt_limit) {
#ifdef HARD_DEBUG
        assert(!group_mask_map_[_igroup]);
        assert(_time_end>clist_[_igroup].getTime());
#endif
        ARChain* c = &clist_[_igroup];
        ARC_pert_pars* par = &par_list_[_igroup];
        //PS::F64 ds_up_limit = 0.25*_dt_limit/c->calc_dt_X(1.0,*ARC_control_);
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
        //PS::F64 ds_use = c->calc_next_step_custom(*ARC_control_,par);
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

        PS::S64 stepcount = c->Symplectic_integration_tsyn(ds_use, *ARC_control_, _time_end, par, &pert_[ipert], &pforce_[ipert], pert_n_[_igroup],fix_step_flag, step_count_limit);

#ifdef HARD_DEBUG
        assert(_time_end+0.5*_dt_limit>c->getTime());
#endif

#ifdef ARC_WARN
        if(c->info!=NULL) {
            c->info->ErrMessage(std::cerr);
        }
#endif
        
#ifdef ARC_DEBUG_DUMP
        if(stepcount<0) {
            dump("ARC_dump.dat",_igroup,_time_end,ds_use);
            std::cerr<<"Igroup = "<<_igroup<<" N = "<<c->getN()<<" Np = "<<pert_n_[_igroup]<<" stepcount = "<<stepcount<<std::endl;
        }
#endif
        
        return stepcount;
    }
#else

    //! integration arc with extrapolation method
    /*! 
      @param[in] _igroup: the ARC group id
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
        ARChain* c = &clist_[_igroup];
        ARC_pert_pars* par = &par_list_[_igroup];
        PS::F64 dscoff=1.0;
        PS::F64 ds_up_limit = 0.25*_dt_limit/c->calc_dt_X(1.0,*ARC_control_);
        PS::F64 ds_use = c->calc_next_step_custom(*ARC_control_,par);
        //PS::F64 ds_use = 0.5*bininfo[_igroup].tstep*std::abs(c->GetPt());
        
        if (ds_use>ds_up_limit) ds_use = ds_up_limit;

        PS::S64 nstep=0;

        // convergency check
        PS::S32 converge_count=0;
        PS::S32 error_count=0;
        bool modify_step_flag=false;
        bool final_flag=false;

        while(_time_end-c->getTime()>ARC_control_->dterr*c->getTime()) {
            const PS::S32 ipert = pert_disp_[_igroup];
            PS::F64 dsf=c->extrapolation_integration(ds_use, *ARC_control_, _time_end, par, &pert_[ipert], &pforce_[ipert], pert_n_[_igroup]);
            if (dsf<0) {
                final_flag=true;
                converge_count++;
                if (converge_count>10&&_time_end-c->getTime()>ARC_control_->dterr*100) {
                    std::cerr<<"Error: Time synchronization fails!\nStep size ds: "<<ds_use<<"\nEnding physical time: "<<_time_end<<"\nTime difference: "<<_time_end-c->getTime()<<"\nR_in: "<<Int_pars_->rin<<"\nR_out: "<<Int_pars_->rout<<"\n";
                    dump("ARC_dump.dat",_igroup,_time_end,ds_use);
                    return -1;
                }
                else ds_use *= -dsf;
            }
            else if (dsf==0) {
                c->info->ErrMessage(std::cerr);
                error_count++;
                if(error_count>4) {
                    std::cerr<<"Error: Too much error appear!\nStep size ds: "<<ds_use<<"\nEnding physical time: "<<_time_end<<"\nTime difference: "<<_time_end-c->getTime()<<"\nR_in: "<<Int_pars_->rin<<"\nR_out: "<<Int_pars_->rout<<"\n";
                    dump("ARC_dump.dat",_igroup,_time_end,ds_use);
                    return -1;
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
                    if (converge_count>10&&_time_end-c->getTime()>ARC_control_->dterr*100) {
                        std::cerr<<"Error: Time synchronization fails!\nStep size ds: "<<ds_use<<"\nEnding physical time: "<<_time_end<<"\nTime difference: "<<_time_end-c->getTime()<<"\nR_in: "<<Int_pars_->rin<<"\nR_out: "<<Int_pars_->rout<<"\n";
                        dump("ARC_dump.dat",_igroup,_time_end,ds_use);
                        return -1;
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
#ifdef ARC_SYM
            nstep_i = integrateOneStepSym(_act_list[i], _time_end, _dt_limit);
#else
            nstep_i = integrateOneStepExt(_act_list[i], _time_end, _dt_limit);
#endif
            if (nstep_i<0) return nstep_i; // error case
            else nstep += nstep_i;
        }

        return nstep;
    }

    //! Integrate active ARC groups
    /* @param[in] _time_end: end of physical integration time
       @param[in] _dt_limit: physical time step upper limit
       \return stepcount; if negative, error happen
     */
    PS::S64 integrateOneStepList(const PS::F64 _time_end,
                                 const PS::F64 _dt_limit) {
        PS::S64 nstep = 0;
        for(PS::S32 i=0; i<clist_.size(); i++) {
            if(getMask(i)) continue;
            PS::S64 nstep_i;
#ifdef ARC_SYM
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
      @param[in] _v_max: maximum velocity used to calcualte r_search
      \return rsearch maximum
     */
    PS::F64 updateRSearch(const PS::F64 _dt_tree, const PS::F64 _v_max) {
        PS::F64 dt_reduce_factor=1.0;
        for(PS::S32 i=0; i<clist_.size(); i++) {
            if(getMask(i)) continue;
            PS::F64 dt_reduce_fi = clist_[i].calcRSearch(_dt_tree, _v_max);
            dt_reduce_factor = std::max(dt_reduce_fi, dt_reduce_factor);
            TpARC** ipadr=clist_[i].getPAdr();
            for (PS::S32 k=0; k<clist_[i].getN(); k++)
                ipadr[k]->r_search = clist_[i].r_search;
        }
        return dt_reduce_factor;
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
        for (PS::S32 i=0; i<clist_.size(); i++) {
            if(getMask(i)) continue;
            // obtain maximum distance pair
            clist_[i].getRmaxIndex(r_max2, r_max_index);
            bool break_flag = false;
            bool out_flag;
            if(r_max2>_r_crit2) {
                // check whether outcome or income
                out_flag=clist_[i].getDirection(r_max_index);
                if(out_flag) break_flag = true;
            }
            // check strong perturbed case
            PS::F64 fp_sq = clist_[i].slowdown.getFratioSq();
            if (fp_sq>2 && bininfo[i].semi>0) break_flag = true;

            // check few-body inner perturbation
            PS::F64 frinsqi=1.0, frinsqj=1.0;
            if(clist_[i].getN()>2) {
                out_flag=clist_[i].getDirection(r_max_index);
                if(out_flag) {
                    clist_[i].getFratioInnerSq(frinsqi, frinsqj, r_max_index, *ARC_control_, Int_pars_);
                    PS::F64 sd_factor = clist_[i].slowdown.getSDRef();
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
                         <<" N_member: "<<clist_[i].getN()
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

    //! Get perturbation factor square
    /*! get current maximum perturbation/inner force square
      @param[in] i: index of groups
      \return fpert square
    */
    PS::F64 getFratioSq(const PS::S32 i) const{
        return clist_[i].slowdown.getFratioSq();
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
            if(getMask(i)) continue;
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
