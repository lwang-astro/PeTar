#pragma once

#include "AR/Float.h"
#include "AR/symplectic_integrator.h"
#include "hermite_ptcl.hpp"
#include "ar_interaction.hpp"
#include "ar_perturber.hpp"
#include "ar_information.hpp"
#include "kepler.hpp"

class ARIntegrator{
private:
    typedef AR::SymplecticIntegrator<Ptcl, PtclH4, ARPerturber, ARInteraction, ARInformation> ARSym;
    SymplecticManager<ARInteraction>* manager_;
    PS::ReallocatableArray<ARSym> groups_;
    PS::ReallocatableArray<bool> group_mask_map_; 
    PS::ReallocatableArray<int> group_mask_list_;
    
public:
    PS::ReallocatableArray<Binary> bininfo;

    //! Reserve memory
    void reserveMem(const int _n) {
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
    int addOneGroup(Tptcl* _ptcl,
                        const int* _ptcl_list,
                        const int _n_ptcl,
                        const Tpsoft* _ptcl_soft_pert,
                        const int _n_split,
                        PtclH4* _ptcl_pert = NULL,
                        const int* _ptcl_pert_list = NULL,
                        const int _n_pert = 0) {
        // set current group offset
        const int ngroup = groups_.size();
        int igroup;

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
            for(int i=0; i<_n_ptcl; i++) {
                group.addMemberAndAddress(_ptcl[_ptcl_list[i]]);
            }
        }
        else { // read one by one
            for(int i=0; i<_n_ptcl; i++) {
                group.addMemberAndAddress(_ptcl[i]);
            }
        }
        
        // calculate the c.m.
        group.particles.calcCenterOfMass();
        // shift to c.m. frame
        group.particles.shiftToCenterOfMassFrame();

        // Add perturber
        group.perturber.neighbor.setMode(AR::ListMode::local);
        for(int i=0; i<_n_pert; i++) {
            const int  k = _ptcl_pert_list[i];
            group.perturber.neighbor.addMember(&_ptcl_pert[k]);
        }

        return igroup;
    }

    //! Clear one group
    /*! Clear one group, remove members, relink to c.m., member number to 0
     */
    void clearOneGroup(const int _igroup) {
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
    void updatePertOneGroup(const int _igroup,
                            PtclH4* _ptcl_pert,
                            const int* _ptcl_pert_list,
                            const int _n_pert) {
#ifdef HARD_DEBUG
        assert(!group_mask_map_[_igroup]);
#endif
        // Add one for c.m.
        ARSym& group = groups_[_igroup];
        group.perturber.neighbor.resizeNoInitialize(0);
        for (int i=0; i<_n_pert; i++) {
            int adr = _ptcl_pert_list[i];
            group.perturber.neighbor.addMember(&_ptcl_pert[adr]);
        }
    }

#ifdef HARD_DEBUG
    //! check perturber list
    bool checkPert() {
        for (int i=0; i<groups_.size(); i++) {
            if(!group_mask_map_[i]){
                ARSym& group = groups_[i];
                const int n_member = group.particles.getSize();
                const int n_pert = group.perturber.neighbor.getSize();
                for (int j=0; j<n_pert; j++) {
                    for (int k=0; k<n_member; k++) {
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
                    if (group.particles.cm.id==group.perturber.neighbor[j].id) {
                        std::cerr<<"Error: cm id== perturber id\n";
                        return true;
                    }
#else
                    assert(group.particles.cm.id!=group.perturber.neighbor[j].id);
#endif
                }
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
    void initialOneSlowDownUnPert(const int _i_group, const Float _tend, const Float _sdfactor) {
        // isolated case
        assert(_i_group<bininfo.size());
        Binary& bin = bininfo[_i_group];
        ARSym& group = groups_[_i_group];

        if (bin.semi>0 && bin.stable_factor>=0) {   
            // apo-center acceleration
            Float apo = bin.semi * (1.0 + bin.ecc);
            Float mtot = (bin.m1 + bin.m2);
            Float acc_apo = mtot/(apo*apo);

            // change bin ecca to apo-center
            Float ecca_bk = bin.ecca;
            bin.ecca = PI;
            // get apo-center position and velocity
            Ptcl p[2];
            p[0].mass = bin.m1;
            p[1].mass = bin.m2;
            OrbParam2PosVel(p[0], p[1], bin);
            bin.ecca = ecca_bk;
#ifdef SOFT_PERT
#ifndef TIDAL_TENSOR
            p[0].status = 0;
            p[1].status = 1;
#endif
#endif
            // soft perturbation force
            AR::Force force[2];
            manager->interaction.calcAccAndSlowDownPert(force, p, 2, group.particles.cm, group.perturber);
            Float acc_pert_sq = 0.0;
            for (int k=0; k<3; k++) {
                Float dacc = force[1].acc_pert[k]-force[0].acc_pert[k];
                fpertsq += dacc*dacc;
            }
            Float acc_pert = sqrt(acc_pert_sq);
            group.initialSlowDownReference(_sdfactor/(mtot*Ptcl::mean_mass_inv), std::max(1.0, _tend/bin.peri));
            group.slowdown.calcSlowDownFactor(acc_apo, acc_pert_sq);
        }
    }

    //! Set initial slowdown parameter for one group
    /*! 
      @param[in] _i_group: group index
      @param[in] _dt_limit_hard: maximum time step for hermite (hard)
      @param[in] _sdfactor: slowdown criterion factor
      @parma[in] _set_one_flag: true: set initial slowdown factor to 1.0
    */
    void initialOneSlowDown(const int _i_group, const Float _dt_limit_hard, const Float _sdfactor, bool _set_one_flag=false) {
        assert(_i_group<bininfo.size());
        Binary& bin = bininfo[_i_group];
        ARSym& group = groups_[_i_group];
        
        if (bin.semi>0 && bin.stable_factor>=0) {
            Float mtot = bin.m1 + bin.m2;
            group.initialSlowDownReference(_sdfactor/(mtot*Ptcl::mean_mass_inv), std::max(1.0, _dt_limit_hard/bin.peri));
            if (_set_one_flag) group.slowdown.setSlowDownFactor(1.0);
        }
    }

    //! Initial one group chain integration parameters
    /*! Initial one group chain integration parameters
      @param[in] _i_group: group to initialize
      @param[in] _time_sys: time to initialize
      \return fail_flag
     */
    bool initialOneGroup(const int _i_group, const Float _time_sys) {
#ifdef HARD_DEBUG
#ifdef HARD_DEBUG_DUMP
        if (group_mask_map_[_i_group]) return true;
#else
        assert(!group_mask_map_[_i_group]);
#endif
#endif
        groups_[_i_group].initial(_time_sys);
        return false;
    }

    //! integration arc with symplectic method
    /*! 
      @param[in] _igroup: the AR group id
      @param[in] _time_end: finishing time
      @param[in] _dt_limit: physical time step limit
      \return fail_flag
     */
    bool integrateOneStep(const int _igroup,
                          const Float _time_end,
                          const Float _dt_limit) {
#ifdef HARD_DEBUG
        assert(!group_mask_map_[_igroup]);
        assert(_time_end>groups_[_igroup].getTime());
#endif
        ARSym& group = groups_[_igroup];
        Float ds_use=bininfo[_igroup].tstep;
        Float kappa_org = group.slowdown.getSlowDownFactorOrigin();
        // in strong perturbed case, avoid too large step size
        if(kappa_org<1.0) {
            ds_use *= 1.0/8.0*std::pow(korg,1.0/6.0);
        }
        // for standard case, not fix step
        int fix_step_option = FixStepOption::none;
        // for two-body or few-body with closed orbit, fix step to avoid energy jump
        const int n_particle = group.particles.getSize();
        if (n_particle==2||(n_particle>2&&bininfo[_igroup].semi>0)) fix_step_option = FixStepOption::later;

        bool fail_flag = group.integrateToTime(ds_use, _time_end, fix_step_option);

#ifdef HARD_DEBUG
        assert(_time_end+0.5*_dt_limit>group.getTime());
#endif

        return fail_flag;
    }

    //! Integrate active AR groups
    /* @param[in] _act_list: active AR group index list
       @param[in] _n_act: number of active groups
       @param[in] _time_end: end of physical integration time
       @param[in] _dt_limit: physical time step upper limit
       \return fail_flag
     */
    bool integrateOneStepList(int _act_list[],
                                 int _n_act,
                                 const Float _time_end,
                                 const Float _dt_limit) {
        long long int nstep = 0;
        for(int i=0; i<_n_act; i++) {
            if(getMask(i)) continue;
            bool fail_flag = integrateOneStep(_act_list[i], _time_end, _dt_limit);
            if (fail_flag) return true; // error case
        }
        return false;
    }

    //! Integrate active AR groups
    /* @param[in] _time_end: end of physical integration time
       @param[in] _dt_limit: physical time step upper limit
       \return fail_flag;
     */
    bool integrateOneStepList(const Float _time_end,
                                 const Float _dt_limit) {
        long long int nstep = 0;
        for(int i=0; i<groups_.size(); i++) {
            if(getMask(i)) continue;
            bool fail_flag = integrateOneStep(i, _time_end, _dt_limit);
            if (fail_flag) return true; // error case
        }
        return false;
    }

    //! Update the c.m. data from the original particle data
    /*!
      @param[in] _ptcl: original particle array
      @param[in] _ptcl_list: particle index list need to be updated, assume _ptcl and AR group have consistent index
      @param[in] _n_ptcl: number of particles
     */
    template <class Tptcl>
    void updateCM(Tptcl* _ptcl,
                  int* _ptcl_list,
                  int _n_ptcl) {
        for(int i=0; i<_n_ptcl; i++) {
            int k = _ptcl_list[i];
#ifdef HARD_DEBUG
            assert(k<groups_.size());
#endif
            groups_[k].particles.cm.pos =  _ptcl[k].pos;
            groups_[k].particles.cm.vel =  _ptcl[k].vel;
            groups_[k].particles.cm.mass = _ptcl[k].mass;
        }
    }

    //! Update the c.m. data from the original particle data
    /*!
      @param[in] _ptcl: original particle array
    */
    template <class Tptcl>
    void updateCM(Tptcl _ptcl[]) {
        for(int i=0; i<groups_.size(); i++) {
            groups_[i].particles.cm.pos  = _ptcl[i].pos;
            groups_[i].particles.cm.vel  = _ptcl[i].vel;
            groups_[i].particles.cm.mass = _ptcl[i].mass;
        }
    }

    //! update rsearch of components based on c.m.
    /*!
      @param[in] _dt_tree: tree time step
      @param[in] _v_max: maximum velocity used to calcualte r_search
      \return rsearch maximum
     */
    Float updateRSearch(const Float _dt_tree, const Float _v_max) {
        Float dt_reduce_factor=1.0;
        for(int i=0; i<groups_.size(); i++) {
            if(getMask(i)) continue;
            Float dt_reduce_fi = groups_[i].calcRSearch(_dt_tree, _v_max);
            dt_reduce_factor = std::max(dt_reduce_fi, dt_reduce_factor);
            TpAR** ipadr=groups_[i].getPAdr();
            for (int k=0; k<groups_[i].getN(); k++)
                ipadr[k]->r_search = groups_[i].r_search;
        }
        return dt_reduce_factor;
    }

    //! Shift member ptcls to their c.m. frame
    /*! Shift all group members to their c.m. frame for AR integration
     */
    void shift2CM() {
        for(int i=0; i<groups_.size(); i++) {
            groups_[i].particles.shiftToCenterOfMassFrame();
        }
    }

    //! resolve all groups' member particles
    /*! shift compotent coordinates to original frame and save data to original particle address
     */
    void resolve() {
        for(int i=0; i<groups_.size(); i++) {
            groups_[i].particles.shiftToOriginFrame();
            groups_[i].writeBackMemberAll();
        }
    }

    //! resolve a list of groups
    /*! shift compotent coordinates to original frame and save data to original particle address
      @param[in] _group_list: group list to resovle
      @param[in] _n_group: number of group to resolve
     */
    void resolve(const int _group_list[], const int _n_group) {
        for(int i=0; i<_n_group; i++) {
            int k = _group_list[i];
            groups_[k].shiftToOriginFrame();
            groups_[k].writeBackMemberAll();
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
    int checkBreak(int* _break_group_list,
                   int* _break_isplit_list,
                   const Float _r_crit2) {
        Float r_max2;
        int r_max_index;
        int n_group_break=0;
        for (int i=0; i<groups_.size(); i++) {
            if(getMask(i)) continue;
            ARSym& group = groups_[i];
            const int n_particle = group.particles.getSize();
            Ptcl* particles = group.particles.getDataAddress();
            // two-body case
            Float dr2, drdv;
            // generate binary tree
            info.generateBinaryTree(particles, n_particle, manager->interaction.dt_tree, manager->interaction.v_max);
            PtclTree<Ptcl>& tree_root = info.getBinaryTreeRoot();
            group.info.getDrDv(dr2, drdv, *tree_root.member[0], *tree_root.member[1]);
                
            // obtain maximum distance pair
            bool break_flag = false;
            // check whether outcome or income
            if(dr2>_r_crit2 && drdv>0.0) break_flag = true;
            // check strong perturbed case
            Float kappa_org = group.slowdown.getSlowDownFactorOrigin();
            if (kappa_org<0.01 && tree_root.semi>0) break_flag = true;

            // check few-body inner perturbation
            Float frinsqi=1.0, frinsqj=1.0;
            if(n_particle>2 && drdv>0.0) {
                AR::SlowDown sd;
                sd.initialSlowDownReference(group.slowdown.getSlowDownFactorReference(),group.slowdown.getSlowDownFactorMax());
                for (int k=0; k<2; k++) {
                    int n_group_sub = tree_root.member[k]->status;
                    if (n_group_sub>1) {
                        PtclTree<Ptcl>* tree_sub = (PtclTree<Ptcl>*) tree_root.member[k];
                        Float semi_db = 2.0*tree_sub->semi;
                        // inner hyperbolic case
                        if(semi_db<0.0) {
                            break_flag = true;
                            break;
                        }
                        Float f_in = tree_sub->m1*tree_sub->m2/(semi_db*semi_db*semi_db);
                        Float f_out = tree_sub->mass*tree_root.member[k-1].mass/(sqrt(dr2)*dr2);
                        Float kappa_in = sd.calcSlowDownFactor(f_in, f_out);
                        // if slowdown factor is large, break the group
                        if (kappa_in>1.0) {
                            break_flag = true;
                            break;
                        }
#ifdef ADJUST_GROUP_DEBUG
                        std::cout<<"Check inner kappa, i_group:"<<i<<" k:"<<k<<" kappa:"<<kappa_in<<std::endl;
#endif
                    }
                }
            }

            if (break_flag) {
#ifdef ADJUST_GROUP_DEBUG
                std::cout<<"Break group, group index: "<<i
                         <<" N_member: "<<n_particle;
                         <<" break index: "<<
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
    int getPtclAdrChain(TpAR* _list[],
                            const int _igroup){
        return groups_[_igroup].getPAdrChain(_list);
    }
    

    //! Get number of members in group i
    /*!
      @param[in] _igroup: group index
      \return number of members in group i
     */
    int getGroupN(const int _igroup) const {
        return groups_[_igroup].getN();
    }

    //! Get number of groups
    /*!
      \return number of groups
     */
    int getNGroups() const {
        return groups_.size();
    }

    //! return integration mask
    /*!
      @param[in] _igroup: index of group
      \return true: suppressed for integration
     */
    bool getMask(const int _igroup) const {
        return group_mask_map_[_igroup];
    }

    //! return CM particle (cast pointer to TpAR*)
    /*! 
      @param[in] _igroup: index of group
      \return TpAR* c.m. particle pointer of index i
     */
    TpAR* getCM(const int _igroup) {
        return (TpAR*)&groups_[_igroup];
    }

    const TpAR* getGroupPtcl(const int i) const {
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
    TpAR** getGroupPtclAdr(const int i) const {
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
    Float getSlowDown(const int i) const{
        return groups_[i].slowdown.getkappa();
    }

    //! Get perturbation factor square
    /*! get current maximum perturbation/inner force square
      @param[in] i: index of groups
      \return fpert square
    */
    Float getFratioSq(const int i) const{
        return groups_[i].slowdown.getFratioSq();
    }

    //! Get slow down original factor
    /*! get slow down original kappa
     */
    Float getSlowDownOrg(const int i) const{
        return groups_[i].slowdown.getkappaorg();
    }

    //! Print slow down parameters
    /*! Print slow down parameters
      @param[in] _os: ofstream for printing
      @param[in] _i: Chain index
      @param[in] _precision: printed precision for one variable
      @param[in] _width: printing width for one variable
     */
    void printSlowDown(std::ostream& _os, const int _i, const int _precision=15, const int _width=23) {
        groups_[_i].slowdown.print(_os,_precision,_width);
    }


//#ifdef AR_PROFILE
//    const long long int getNsubstep() const{
//        long long int Nsum = 0;
//        for (int i=0; i<groups_.size(); i++) 
//            Nsum += groups_[i].profile.itercount;
//        return Nsum;
//    }
//#endif

#ifdef AR_DEBUG_PRINT
    void data_dump(std::ostream& os, const int i, const Float dt_limit) const{
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
        for (int j=0; j<groups_[i].getN(); j++) {
            os<<groups_[i].getP(j).mass<<" "
              <<groups_[i].getP(j).pos<<" "
              <<groups_[i].getP(j).vel<<std::endl;
        }
        for (int j=1; j<pert_n_[i]; j++) {
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
    bool info_print(std::ostream& os, const long long int _n_group, const long long int _n_group_in_cluster, const long long int _n_ptcl, const long long int _n_hint, const Float _dt_limit, const int _kp, const int _n_step_limit=10000) const{
        bool dump_flag=false;
        for (int i=0; i<groups_.size(); i++) {
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
            long long int nstep = 0;
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
    void EnergyRecord(Teng &energy, const int sdflag=false) {
        energy.kin = energy.pot = energy.tot = 0.0;
        Float sd;
        for(int i=0; i<groups_.size(); i++) {
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
        for (int i=0; i<groups_.size(); i++) {
            Floatvec pcm_pos = groups_[i].pos;
            Floatvec pcm_vel = groups_[i].vel;
            for (int j=0; j<groups_[i].getN(); j++) {
                Tptcl pj = groups_[i].getP(j);
                pj.pos += pcm_pos;
                pj.vel += pcm_vel;
                pj.ParticleBase::writeAscii(_fout);
            }
        }
    }
#endif

};
