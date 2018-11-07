#pragma once
#ifdef USE_INTRINSIC_FOR_X86
#include<immintrin.h>
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

    template<class Tptcl>
    PtclHard(const Tptcl& _p, const PS::F64 _r_search, const PS::F64 _mass_bk, const PS::S64 _id, const PS::S64 _status, const PS::S32 _id_cluster, const PS::S32 _adr_org): 
        Ptcl(_p, _r_search, _mass_bk, _id, _status), id_cluster(_id_cluster), adr_org(_adr_org) {}

    template<class Tptcl>
    PtclHard(const Tptcl &_p, const PS::S32 _id_cluster, const PS::S32 _adr_org): 
        Ptcl(_p), id_cluster(_id_cluster), adr_org(_adr_org) {}

    template<class Tptcl>
    PtclHard(const Tptcl &_p) {
        Ptcl::DataCopy(_p);
        id_cluster = _p.id_cluster;
        adr_org = _p.adr_org;
    }

    // notice datacopy ignore id_cluster and adr_org
    template<class Tptcl>
    void DataCopy(const Tptcl& _p) {
        Ptcl::DataCopy(_p);
    }

    template<class Tptcl>
    PtclHard& operator = (const Tptcl& _p) {
        Ptcl::DataCopy(_p);
        id_cluster = _p.id_cluster;
        adr_org = _p.adr_org;
        return *this;
    }

    void dump(FILE *_fout) {
        fwrite(this, sizeof(*this),1,_fout);
    }

    void read(FILE *_fin) {
        size_t rcount = fread(this, sizeof(*this),1,_fin);
        if (rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            abort();
        }
    }
    void print(std::ostream & _fout){
        Ptcl::print(_fout);
        std::cerr<<" id_cluster="<<id_cluster
                 <<" adr_org="<<adr_org;
    }
};


void PtclHardDump(FILE *_fout, PtclHard * _ptcl, const PS::S32 _n) {
    fwrite(&_n, sizeof(PS::S32), 1, _fout);
    for(int i=0; i<_n; i++) {
        _ptcl[i].dump(_fout);
    }
    PS::F64 ptcl_st_dat[3];
    ptcl_st_dat[0] = Ptcl::search_factor;
    ptcl_st_dat[1] = Ptcl::r_search_min;
    ptcl_st_dat[2] = Ptcl::mean_mass_inv;
    fwrite(ptcl_st_dat, sizeof(PS::F64),3,_fout);
}

void PtclHardRead(FILE *_fin, PS::ReallocatableArray<PtclHard> & _ptcl) {
    PS::S32 n;
    size_t rcount = fread(&n, sizeof(PS::S32),1,_fin);
    if (rcount<1) {
        std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
        abort();
    }
    PtclHard ptmp;
    for(int i=0; i<n; i++) {
        ptmp.read(_fin);
        _ptcl.push_back(ptmp);
    }
    PS::F64 ptcl_st_dat[3];
    rcount = fread(ptcl_st_dat, sizeof(PS::F64),3,_fin);
    if (rcount<3) {
        std::cerr<<"Error: Data reading fails! requiring data number is 3, only obtain "<<rcount<<".\n";
        abort();
    }
    Ptcl::search_factor = ptcl_st_dat[0];
    Ptcl::r_search_min  = ptcl_st_dat[1];
    Ptcl::mean_mass_inv = ptcl_st_dat[2];
}


#ifdef HARD_CHECK_ENERGY
class HardEnergy {
public:
    PS::F64 kin, pot, tot;
};
#endif

class SystemHard{
public:
#ifdef PROFILE
    PS::S64 ARC_substep_sum;
    PS::F64 ARC_n_groups;
#endif
#ifdef HARD_CHECK_ENERGY
    PS::F64 hard_dE, hard_dESD;
    PS::F64 hard_dE_limit;
#endif
#ifdef ARC_SYM
    PS::S32 arc_step_count_limit;
#endif

private:
    // Notice: if new variables added, change pardump also
    PS::F64 dt_limit_hard_;
    PS::F64 dt_min_hard_;
    PS::F64 eta_s_;
    PS::F64 time_origin_;
    PS::F64 r_bin_;
    PS::F64 sdfactor_;
    PS::S64 id_offset_;
    PS::S32 n_split_;
    
    PS::ReallocatableArray<PtclHard> ptcl_hard_;
    PS::ReallocatableArray<PS::S32> n_ptcl_in_cluster_;
    PS::ReallocatableArray<PS::S32> n_ptcl_in_cluster_disp_;
    PS::ReallocatableArray<PS::S32> n_group_in_cluster_;
    PS::ReallocatableArray<PS::S32> n_group_in_cluster_offset_;
    PS::ReallocatableArray<PS::S32> adr_first_ptcl_arti_in_cluster_;
    PS::S32 n_group_member_remote_; // number of members in groups but in remote nodes
    ARC::chainpars ARC_control_pert_; ///chain controller for perturbed(L.Wang)
    ARC::chainpars ARC_control_soft_; ///chain controller for no perturber(L.Wang)
    ARC_int_pars Int_pars_; /// ARC integration parameters, rout_, rin_ (L.Wang)

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

    //! Find groups and create aritfical particles to sys
    /* @param[in,out] _sys: global particle system
       @param[in,out] _ptcl_local: local saved particle data (will be reordered due to the groups)
       @param[in]     _n_ptcl_in_cluster: number of particles in one cluster
       @param[in]     _n_ptcl_in_cluster_disp: boundar of particle cluster
       @param[out]    _n_group_in_cluster: number of groups in one cluster
       @param[out]    _n_group_in_cluster_offset: boundary of groups in _adr_first_ptcl_arti_in_cluster
       @param[out]    _adr_first_ptcl_arti_in_cluster: address of the first artifical particle in each groups
       @param[in]     _rbin: binary detection criterion radius
       @param[in]     _rin: inner radius of soft-hard changeover function
       @param[in]     _rout: outer radius of soft-hard changeover function
       @param[in]     _dt_tree: tree time step for calculating r_search
       @param[in]     _id_offset: for artifical particles, the offset of starting id.
       @param[in]     _n_split: split number for artifical particles
     */
    template<class Tsys, class Tptcl>
    void findGroupsAndCreateArtificalParticlesImpl(Tsys & _sys,
                                                   PtclHard* _ptcl_local,
                                                   PS::ReallocatableArray<PS::S32>& _n_ptcl_in_cluster,
                                                   PS::ReallocatableArray<PS::S32>& _n_ptcl_in_cluster_disp,
                                                   PS::ReallocatableArray<PS::S32>& _n_group_in_cluster,
                                                   PS::ReallocatableArray<PS::S32>& _n_group_in_cluster_offset,
                                                   PS::ReallocatableArray<PS::S32>& _adr_first_ptcl_arti_in_cluster,
                                                   const PS::F64 _rbin,
                                                   const PS::F64 _rin,
                                                   const PS::F64 _rout,
                                                   const PS::F64 _dt_tree,
                                                   const PS::S64 _id_offset,
                                                   const PS::S32 _n_split) { 
        const PS::S32 n_cluster = _n_ptcl_in_cluster.size();
#ifdef HARD_DEBUG
        assert(n_cluster<ARRAY_ALLOW_LIMIT);
#endif        
        _n_group_in_cluster.resizeNoInitialize(n_cluster);
        n_group_member_remote_=0;

        const PS::S32 num_thread = PS::Comm::getNumberOfThread();
        PS::ReallocatableArray<PtclHard> ptcl_artifical[num_thread];

#pragma omp parallel for schedule(dynamic)
        for (PS::S32 i=0; i<n_cluster; i++){
            const PS::S32 ith = PS::Comm::getThreadNum();
            PtclHard* ptcl_in_cluster = _ptcl_local + _n_ptcl_in_cluster_disp[i];
            const PS::S32 n_ptcl = _n_ptcl_in_cluster[i];
            // reset status
            for(PS::S32 j=0; j<n_ptcl; j++) ptcl_in_cluster[j].status = 0;
            // search groups
            SearchGroup<PtclHard> group;
            // merge groups
            if (n_ptcl==2) group.searchAndMerge(ptcl_in_cluster, n_ptcl, _rout);
            else group.searchAndMerge(ptcl_in_cluster, n_ptcl, _rin);

            // generate artifical particles,
            group.generateList(i, ptcl_in_cluster, n_ptcl, ptcl_artifical[ith], _n_group_in_cluster[i], _rbin, _rin, _rout, _dt_tree, _id_offset, _n_split);
        }

        // n_group_in_cluster_offset
        _n_group_in_cluster_offset.resizeNoInitialize(n_cluster+1);
        _n_group_in_cluster_offset[0] = 0;
        for (PS::S32 i=0; i<n_cluster; i++) 
            _n_group_in_cluster_offset[i+1] = _n_group_in_cluster_offset[i] + _n_group_in_cluster[i];
#ifdef HARD_DEBUG
        assert(_n_group_in_cluster_offset[n_cluster]<ARRAY_ALLOW_LIMIT);
#endif        
        _adr_first_ptcl_arti_in_cluster.resizeNoInitialize(_n_group_in_cluster_offset[n_cluster]);


        // add artifical particle to particle system
        PS::S32 rank = PS::Comm::getRank();
        // Get the address offset for new artifical ptcl array in each thread in _sys
        PS::S64 sys_ptcl_artifical_thread_offset[num_thread+1];
        sys_ptcl_artifical_thread_offset[0] = _sys.getNumberOfParticleLocal();
        for(PS::S32 i=0; i<num_thread; i++) 
            sys_ptcl_artifical_thread_offset[i+1] = sys_ptcl_artifical_thread_offset[i] + ptcl_artifical[i].size();
        _sys.setNumberOfParticleLocal(sys_ptcl_artifical_thread_offset[num_thread]);
        
#pragma omp parallel for        
        for(PS::S32 i=0; i<num_thread; i++) {
            GroupPars gpar(_n_split);
            const PS::S32 n_artifical_per_group = gpar.n_ptcl_artifical;
            // ptcl_artifical should be integer times of 2*n_split+1
            assert(ptcl_artifical[i].size()%n_artifical_per_group==0);
            // Add particle to ptcl sys
            for (PS::S32 j=0; j<ptcl_artifical[i].size(); j++) {
                PS::S32 adr = j+sys_ptcl_artifical_thread_offset[i];
                ptcl_artifical[i][j].adr_org=adr;
                _sys[adr]=Tptcl(ptcl_artifical[i][j],rank,adr);
            }
            PS::S32 group_offset=0, j_group_recored=-1;
            // Update the status of group members to c.m. address in ptcl sys. Notice c.m. is at the end of an artificial particle group
            for (PS::S32 j=0; j<ptcl_artifical[i].size(); j+=n_artifical_per_group) {
                // obtain group member nember
                gpar.getGroupIndex(&ptcl_artifical[i][j]);
                PS::S32 j_cm = gpar.offset_cm + j;
                PS::S32 n_members = gpar.n_members;
                PS::S32 i_cluster = gpar.i_cluster;
                PS::S32 j_group = gpar.i_group;
                PS::F64 rsearch_member=ptcl_artifical[i][j].r_search;
                // make sure group index increase one by one
                assert(j_group==j_group_recored+1);
                j_group_recored=j_group;
                // update member status
                for (PS::S32 k=0; k<n_members; k++) {
                    PS::S32 kl = _n_ptcl_in_cluster_disp[i_cluster]+group_offset+k;
                    PS::S64 ptcl_k=_ptcl_local[kl].adr_org;
                    if(ptcl_k>=0) {
#ifdef HARD_DEBUG
                        // check whether ID is consistent.
                        if(k==0) assert(_sys[ptcl_k].id==-ptcl_artifical[i][j_cm].id);
#endif
                        // save c.m. address and shift mass to mass_bk, set rsearch
                        _sys[ptcl_k].status = -ptcl_artifical[i][j_cm].adr_org; //save negative address
                        _sys[ptcl_k].mass_bk = _sys[ptcl_k].mass;
                        _sys[ptcl_k].r_search = rsearch_member;
#ifdef SPLIT_MASS
                        _sys[ptcl_k].mass = 0;
#endif
                    }
                    else {
                        // this is remoted member;
                        n_group_member_remote_++;
                    }
#ifdef HARD_DEBUG
                    // check whether ID is consistent.
                    if(k==0) assert(_ptcl_local[kl].id==-ptcl_artifical[i][j_cm].id);
#endif
                    _ptcl_local[kl].status = -ptcl_artifical[i][j_cm].adr_org;
                    _ptcl_local[kl].mass_bk = _ptcl_local[kl].mass;
                    _ptcl_local[kl].r_search = rsearch_member;
#ifdef SPLIT_MASS
                    _ptcl_local[kl].mass = 0;
#endif
                }
                // shift cluster
                if(j_group==_n_group_in_cluster[i_cluster]-1) {
                    group_offset=0;
                    j_group_recored=-1;
                }
                else group_offset += n_members; // group offset in the ptcl list index of one cluster
                // j_group should be consistent with n_group[i_cluster];
                assert(j_group<=_n_group_in_cluster[i_cluster]);

                // save first address of artifical particle
                _adr_first_ptcl_arti_in_cluster[_n_group_in_cluster_offset[i_cluster]+j_group] = ptcl_artifical[i][j].adr_org;
            }
        }
    }

    //! soft force correction use tree neighbor search for one particle
    /*
      @param[in,out] _psoft: particle in global system need to be corrected for acc and pot
      @param[in] _tree: tree for force
      @param[in] _rin: cutoff inner radius;
      @param[in] _rout: cutoff outer radius;
      @param[in] _r_oi_inv: 1.0/(_rout-_rin);
      @param[in] _r_A: (_rout-_rin)/(_rout+_rin);
      @param[in] _eps_sq: softing eps square
     */
    template <class Tpsoft, class Ttree, class Tepj>
    void correctForceWithCutoffTreeNeighborOneParticleImp(Tpsoft& _psoft, 
                                                          Ttree& _tree,
                                                          const PS::F64 _rin,
                                                          const PS::F64 _rout,
                                                          const PS::F64 _r_oi_inv,
                                                          const PS::F64 _r_A,
                                                          const PS::F64 _eps_sq) {
        Tepj * ptcl_nb = NULL;
        PS::S32 n_ngb = _tree.getNeighborListOneParticle(_psoft, ptcl_nb);
#ifdef HARD_DEBUG
        assert(n_ngb >= 1);
#endif
        PS::F64vec& pos_j= _psoft.pos;
        PS::F64vec& acc_j= _psoft.acc;
        PS::F64& pot_j = _psoft.pot_tot;
        PS::S64 stat_j = _psoft.status;

        // self-potential correction 
        // no correction for orbital artifical particles because the potential are not used for any purpose
        // no correction for member particles because their mass is zero during the soft force calculation, the self-potential contribution is also zero.
        if (stat_j==0) pot_j += _psoft.mass/_rout; // single
        //else if (stat_j<0) pot_j += _psoft.mass_bk/_rout; // member

        // loop neighbors
        for(PS::S32 k=0; k<n_ngb; k++){
            if (ptcl_nb[k].id == _psoft.id) continue;
            PS::S32 pot_control_flag;
            PS::S64 stat_k = ptcl_nb[k].status;
            if (stat_k==0) pot_control_flag = 0; //single 
            else if (stat_k<0) pot_control_flag = 1; //member
            else pot_control_flag = 2; //artifical

            CalcAccPotShortWithLinearCutoff(pos_j, acc_j, pot_j, 
                                            ptcl_nb[k].pos, ptcl_nb[k].mass, ptcl_nb[k].mass_bk, 
                                            pot_control_flag, _eps_sq,
                                            _r_oi_inv, _r_A, _rout, _rin);
        }
    }

    //! soft force correction for artifical particles in one cluster
    /* 1. Correct cutoff for artifical particles
       2. The c.m. force is substracted from tidal tensor force
       3. c.m. force is replaced by the averaged force on orbital particles
       @param[in,out] _sys: global particle system, acc is updated
       @param[in] _ptcl_local: particle in systme_hard
       @param[in] _adr_real_start: real particle start address in _ptcl_local
       @param[in] _adr_real_end:   real particle end (+1) address in _ptcl_local
       @param[in] _n_group:  number of groups in cluster
       @param[in] _adr_first_ptcl_arti_in_cluster: address of the first artifical particle in each groups
       @param[in] _rin: cutoff inner radius;
       @param[in] _rout: cutoff outer radius;
       @param[in] _r_oi_inv: 1.0/(_rout-_rin);
       @param[in] _r_A: (_rout-_rin)/(_rout+_rin);
       @param[in] _n_split: artifical particle splitting number
       @param[in] _eps_sq: softing eps square
     */
    template <class Tsys>
    void correctForceWithCutoffArtificalOneClusterImp(Tsys& _sys, 
                                                      const PtclHard* _ptcl_local,
                                                      const PS::S32 _adr_real_start,
                                                      const PS::S32 _adr_real_end,
                                                      const PS::S32 _n_group,
                                                      const PS::S32* _adr_first_ptcl_arti_in_cluster,
                                                      const PS::F64 _rin,
                                                      const PS::F64 _rout,
                                                      const PS::F64 _r_oi_inv,
                                                      const PS::F64 _r_A,
                                                      const PS::S32 _n_split,
                                                      const PS::F64 _eps_sq) {

        GroupPars gpars(_n_split);
        for (int j=0; j<_n_group; j++) {  // j: j_group
            PS::S32 j_start = _adr_first_ptcl_arti_in_cluster[j];
            PS::S32 j_cm = j_start + gpars.offset_cm;

            // loop all artifical particles: tidal tensor, orbital and c.m. particle
            for (int k=j_start; k<=j_cm; k++) {  
                // k: k_ptcl_arti
                PS::F64vec& pos_k= _sys[k].pos;
                PS::F64vec& acc_k= _sys[k].acc;
                PS::F64& pot_k = _sys[k].pot_tot;


                // loop orbital artifical particle
                // group
                for (int kj=0; kj<_n_group; kj++) { // group
                    PS::S32 kj_start = _adr_first_ptcl_arti_in_cluster[kj];
                    PS::S32 kj_cm = kj_start + gpars.offset_cm;

                    // particle arti orbital
                    for (int kk=kj_start+gpars.offset_orb; kk<kj_cm; kk++) {
                        if(kk==k) continue; //avoid same particle
                        auto& ptcl_kk = _sys[kk];
                        CalcAccPotShortWithLinearCutoff(pos_k, acc_k, pot_k, 
                                                        ptcl_kk.pos, ptcl_kk.mass, ptcl_kk.mass_bk, 
                                                        2, _eps_sq,
                                                        _r_oi_inv, _r_A, _rout, _rin);
                    }
                }

                // loop real particle
                for (int kj=_adr_real_start; kj<_adr_real_end; kj++) {
                    const PtclHard* ptcl_kj_ptr = &_ptcl_local[kj];
                    PS::S32 pot_control_flag = ptcl_kj_ptr->status<0? 1: 0;
                    CalcAccPotShortWithLinearCutoff(pos_k, acc_k, pot_k, 
                                                    ptcl_kj_ptr->pos, ptcl_kj_ptr->mass, ptcl_kj_ptr->mass_bk, 
                                                    pot_control_flag, _eps_sq,
                                                    _r_oi_inv, _r_A, _rout, _rin);
                }
            }
            
            // for c.m. particle
            PS::F64vec& acc_cm = _sys[j_cm].acc;

#ifdef TIDAL_TENSOR
            // substract c.m. force (acc) from tidal tensor force (acc)
            for (PS::S32 k=gpars.offset_tt; k<gpars.offset_orb; k++)  _sys[j_start+k].acc -= acc_cm;
#endif
                
            // After c.m. force used, it can be replaced by the averaged force on orbital particles
            acc_cm=PS::F64vec(0.0);
            PS::F64 m_ob_tot = 0.0;

            PS::S32 job_start = j_start+gpars.offset_orb;
            for (PS::S32 k=job_start; k<j_cm; k++) {
                acc_cm += _sys[k].mass*_sys[k].acc; 
                m_ob_tot += _sys[k].mass;
//#ifdef HARD_DEBUG
//                assert(((_sys[k].status)>>ID_PHASE_SHIFT)==-_sys[j_cm].id);
//#endif
            }
            acc_cm /= m_ob_tot;

#ifdef HARD_DEBUG
            assert(abs(m_ob_tot-_sys[j_cm].mass_bk)<1e-10);
#endif

            //PS::F64vec& pos_j= _sys[j_cm].pos;
            //PS::F64vec& acc_j= _sys[j_cm].acc;
            //PS::F64& pot_j = _sys[j_cm].pot_tot;
            // 
            //// loop artifical particle orbital
            //for (int k=0; k<n_group; k++) { // group
            //    PS::S32 k_start = _adr_first_ptcl_arti_in_cluster[_n_group_in_cluster_offset[i]+k];
            //    PS::S32 k_cm = k_start + 2*_n_split;
            //    for (int ki=k_start+_n_split; ki<k_cm; ki++) {
            //        auto& ptcl_k = _sys[ki];
            //        CalcAccPotShortWithLinearCutoff(pos_j, acc_j, pot_j, 
            //                                        ptcl_k.pos, ptcl_k.mass, ptcl_k.mass_bk, 
            //                                        2, _eps_sq,
            //                                        r_oi_inv, r_A, _rout, _rin);
            //    }
            //}
            // 
            //// loop real particle
            //for (int k=_n_ptcl_in_cluster_offset[i]; k<_n_ptcl_in_cluster_offset[i+1]; k++) {
            //    PtclHard* ptcl_k_ptr = &_ptcl_local[k];
            //    PS::S32 pot_control_flag = ptcl_k_ptr->status>0? 1: 0;
            //    CalcAccPotShortWithLinearCutoff(pos_j, acc_j, pot_j, 
            //                                    ptcl_k_ptr->pos, ptcl_k_ptr->mass, ptcl_k_ptr->mass_bk, 
            //                                    pot_control_flag, _eps_sq,
            //                                    r_oi_inv, r_A, _rout, _rin);
            //}
        }
    }

    //! Soft force correction due to different cut-off function
    /* Use cluster member
       1. first correct for artifical particles, then for cluster member. 
       2. The c.m. force is substracted from tidal tensor force
       3. c.m. force is replaced by the averaged force on orbital particles

       @param[in,out] _sys: global particle system, acc is updated
       @param[in] _ptcl_local: particle in systme_hard
       @param[in] _n_ptcl_in_cluster: number of particles in clusters
       @param[in] _n_ptcl_in_cluster_offset: boundary of clusters in _adr_sys_in_cluster
       @parma[in] _n_group_in_cluster: number of groups in clusters
       @param[in] _n_group_in_cluster_offset: boundary of groups in _adr_first_ptcl_arti_in_cluster
       @param[in] _adr_first_ptcl_arti_in_cluster: address of the first artifical particle in each groups
       @param[in] _rin: cutoff inner radius;
       @param[in] _rout: cutoff outer radius;
       @param[in] _n_split: artifical particle splitting number
       @param[in] _eps_sq: softing eps square
    */
    template <class Tsys>
    void correctForceWithCutoffClusterImp(Tsys& _sys, 
                                          const PtclHard* _ptcl_local,
                                          const PS::ReallocatableArray<PS::S32>& _n_ptcl_in_cluster,
                                          const PS::ReallocatableArray<PS::S32>& _n_ptcl_in_cluster_offset,
                                          const PS::ReallocatableArray<PS::S32>& _n_group_in_cluster,
                                          const PS::ReallocatableArray<PS::S32>& _n_group_in_cluster_offset,
                                          const PS::ReallocatableArray<PS::S32>& _adr_first_ptcl_arti_in_cluster,
                                          const PS::F64 _rin,
                                          const PS::F64 _rout,
                                          const PS::S32 _n_split,
                                          const PS::F64 _eps_sq) {
        // cutoff function parameter
        const PS::F64 r_oi_inv = 1.0/(_rout-_rin);
        const PS::F64 r_A = (_rout-_rin)/(_rout+_rin);
        //const PS::S32 ptcl_art_offset = 2*_n_split+1;

        const PS::S32 n_cluster = _n_ptcl_in_cluster.size();
#pragma omp parallel for schedule(dynamic)
        for (int i=0; i<n_cluster; i++) {  // i: i_cluster
            PS::S32 adr_real_start= _n_ptcl_in_cluster_offset[i];
            PS::S32 adr_real_end= _n_ptcl_in_cluster_offset[i+1];
            // artifical particle group number
            PS::S32 n_group = _n_group_in_cluster[i];
            //PS::S32 n_group_offset = _n_group_in_cluster_offset[i];
            const PS::S32* adr_first_ptcl_arti = &_adr_first_ptcl_arti_in_cluster[_n_group_in_cluster_offset[i]];

            // correction for artifical particles
            correctForceWithCutoffArtificalOneClusterImp(_sys, _ptcl_local, adr_real_start, adr_real_end, n_group, adr_first_ptcl_arti, _rin, _rout, r_oi_inv, r_A, _n_split, _eps_sq);

            // obtain correction for real particles in clusters
            GroupPars gpars(_n_split);
            for (int j=adr_real_start; j<adr_real_end; j++) {
                PS::S64 adr = _ptcl_local[j].adr_org;
#ifdef HARD_DEBUG
                assert(_sys[adr].id==_ptcl_local[j].id);
#endif
                PS::F64vec& pos_j= _sys[adr].pos;
                PS::F64vec& acc_j= _sys[adr].acc;
                PS::F64& pot_j = _sys[adr].pot_tot;
            
                PS::S64 stat_j = _sys[adr].status;
                //self-potential correction for non-group member, group member has mass zero, so no need correction
                if(stat_j==0) pot_j += _sys[adr].mass/_rout;
                //else if(stat_j<0) pot_j += _sys[adr].mass_bk/_rout; 

                // cluster member
                for (int k=adr_real_start; k<adr_real_end; k++) {
                    if(k==j) continue;
                    const PtclHard* ptcl_k_ptr = &_ptcl_local[k];
                    PS::S32 pot_control_flag = ptcl_k_ptr->status<0? 1: 0;
                    CalcAccPotShortWithLinearCutoff(pos_j, acc_j, pot_j, 
                                                    ptcl_k_ptr->pos, ptcl_k_ptr->mass, ptcl_k_ptr->mass_bk, 
                                                    pot_control_flag, _eps_sq,
                                                    r_oi_inv, r_A, _rout, _rin);
                }

                // orbital artifical particle
                for (int k=0; k<n_group; k++) {
                    // loop artifical particle orbital
                    PS::S32 k_start = adr_first_ptcl_arti[k];
                    PS::S32 k_cm = k_start + gpars.offset_cm;
                    for (int ki=k_start+gpars.offset_orb; ki<k_cm; ki++) {
                        auto& ptcl_k = _sys[ki];
                        CalcAccPotShortWithLinearCutoff(pos_j, acc_j, pot_j, 
                                                        ptcl_k.pos, ptcl_k.mass, ptcl_k.mass_bk, 
                                                        2, _eps_sq,
                                                        r_oi_inv, r_A, _rout, _rin);
                    }
                }
            
//#ifdef HARD_DEBUG
//                if(stat_j<0) assert(-_sys[-stat_j].id==_ptcl_local[j].id);
//#endif
                //// group member, use c.m. acc
                //if(stat_j>0) acc_j = _sys[stat_j].acc;
            }
        }
    }

//! Soft force correction due to different cut-off function
/* Use tree neighbor search
   1. first correct for artifical particles use cluster information, 
   2. The c.m. force is substracted from tidal tensor force
   3. c.m. force is replaced by the averaged force on orbital particles
   4. then use tree neighbor search for local cluster real member. 
   @param[in] _sys: global particle system, acc is updated
   @param[in] _tree: tree for force
   @param[in] _ptcl_local: particle in systme_hard
   @param[in] _n_ptcl_in_cluster: number of particles in clusters
   @param[in] _n_ptcl_in_cluster_offset: boundary of clusters in _adr_sys_in_cluster
   @parma[in] _n_group_in_cluster: number of groups in clusters
   @param[in] _n_group_in_cluster_offset: boundary of groups in _adr_first_ptcl_arti_in_cluster
   @param[in] _adr_first_ptcl_arti_in_cluster: address of the first artifical particle in each groups
   @param[in] _adr_send: particle in sending list of connected clusters
   @param[in] _rin: cutoff inner radius;
   @param[in] _rout: cutoff outer radius;
   @param[in] _n_split: artifical particle splitting number
   @param[in] _eps_sq: softing eps square
*/
    template <class Tsys, class Tpsoft, class Ttree, class Tepj>
    void correctForceWithCutoffTreeNeighborAndClusterImp(Tsys& _sys,
                                                         Ttree& _tree,
                                                         const PtclHard* _ptcl_local,
                                                         const PS::ReallocatableArray<PS::S32>& _n_ptcl_in_cluster,
                                                         const PS::ReallocatableArray<PS::S32>& _n_ptcl_in_cluster_offset,
                                                         const PS::ReallocatableArray<PS::S32>& _n_group_in_cluster,
                                                         const PS::ReallocatableArray<PS::S32>& _n_group_in_cluster_offset,
                                                         const PS::ReallocatableArray<PS::S32>& _adr_first_ptcl_arti_in_cluster,
                                                         const PS::ReallocatableArray<PS::S32>& _adr_send,
                                                         const PS::F64 _rin,
                                                         const PS::F64 _rout,
                                                         const PS::S32 _n_split,
                                                         const PS::F64 _eps_sq) {

        // cutoff function parameter
        const PS::F64 r_oi_inv = 1.0/(_rout-_rin);
        const PS::F64 r_A = (_rout-_rin)/(_rout+_rin);
        //const PS::S32 ptcl_art_offset = 2*_n_split+1;

        const PS::S32 n_cluster = _n_ptcl_in_cluster.size();

#pragma omp parallel for schedule(dynamic)
        for (int i=0; i<n_cluster; i++) {  // i: i_cluster
            PS::S32 adr_real_start= _n_ptcl_in_cluster_offset[i];
            PS::S32 adr_real_end= _n_ptcl_in_cluster_offset[i+1];
            // artifical particle group number
            PS::S32 n_group = _n_group_in_cluster[i];
            //PS::S32 n_group_offset = _n_group_in_cluster_offset[i];
            const PS::S32* adr_first_ptcl_arti = &_adr_first_ptcl_arti_in_cluster[_n_group_in_cluster_offset[i]];

            // correction for artifical particles
            correctForceWithCutoffArtificalOneClusterImp(_sys, _ptcl_local, adr_real_start, adr_real_end, n_group, adr_first_ptcl_arti, _rin, _rout, r_oi_inv, r_A, _n_split, _eps_sq);

            // obtain correction for real particles in clusters use tree neighbor search
            for (int j=adr_real_start; j<adr_real_end; j++) {
                PS::S64 adr = _ptcl_local[j].adr_org;
                // only do for local particles
#ifdef HARD_DEBUG
                if(adr>=0) assert(_sys[adr].id==_ptcl_local[j].id);
#endif
                if(adr>=0) correctForceWithCutoffTreeNeighborOneParticleImp<Tpsoft, Ttree, Tepj>(_sys[adr], _tree, _rin, _rout, r_oi_inv, r_A, _eps_sq);
            }
        }

        const PS::S32 n_send = _adr_send.size();
#pragma omp parallel for 
        // sending list to other nodes need also be corrected.
        for (int i=0; i<n_send; i++) {
            PS::S64 adr = _adr_send[i];
            correctForceWithCutoffTreeNeighborOneParticleImp<Tpsoft, Ttree, Tepj>(_sys[adr], _tree, _rin, _rout, r_oi_inv, r_A, _eps_sq); 
        }
    }
    
//! soft force correction completely use tree neighbor search
/* @param[in,out] _sys: global particle system, acc is updated
   @param[in] _tree: tree for force
   @param[in] _ptcl_local: particle in systme_hard, only used to get adr_org
   @param[in] _n_ptcl: total number of particles in all clusters
   @param[in] _adr_ptcl_artifical_start: start address of artifical particle in _sys
   @param[in] _rin: cutoff inner radius;
   @param[in] _rout: cutoff outer radius;
   @param[in] _n_split: artifical particle splitting number
   @param[in] _eps_sq: softing eps square
*/
    template <class Tsys, class Tpsoft, class Ttree, class Tepj>
    void correctForceWithCutoffTreeNeighborImp(Tsys& _sys, 
                                               Ttree& _tree, 
                                               const PtclHard* _ptcl_local,
                                               const PS::S32 _n_ptcl,
                                               const PS::S32 _adr_ptcl_artifical_start,
                                               const PS::F64 _rin,
                                               const PS::F64 _rout,
                                               const PS::S32 _n_split,
                                               const PS::F64 _eps_sq) { 
        // cutoff function parameter
        const PS::F64 r_oi_inv = 1.0/(_rout-_rin);
        const PS::F64 r_A = (_rout-_rin)/(_rout+_rin);

        // for real particle
#pragma omp parallel for schedule(dynamic)
        for (int i=0; i<_n_ptcl; i++) {
            PS::S64 adr = _ptcl_local[i].adr_org;
            if(adr>=0) correctForceWithCutoffTreeNeighborOneParticleImp<Tpsoft, Ttree, Tepj>(_sys[adr], _tree, _rin, _rout, r_oi_inv, r_A, _eps_sq);
        }

        // for artifical particle
        const PS::S32 n_tot = _sys.getNumberOfParticleLocal();
#pragma omp parallel for schedule(dynamic)
        for (int i=_adr_ptcl_artifical_start; i<n_tot; i++) 
            correctForceWithCutoffTreeNeighborOneParticleImp<Tpsoft, Ttree, Tepj>(_sys[i], _tree, _rin, _rout, r_oi_inv, r_A, _eps_sq);

        
        GroupPars gpars(_n_split);
#ifdef HARD_DEBUG
        assert((n_tot-_adr_ptcl_artifical_start)%gpars.n_ptcl_artifical==0);
#endif
#pragma omp parallel for schedule(dynamic)
        for (int i=_adr_ptcl_artifical_start; i<n_tot; i+=gpars.n_ptcl_artifical){
            PS::S32 i_cm = i + gpars.offset_cm;
            PS::F64vec& acc_cm = _sys[i_cm].acc;
#ifdef TIDAL_TENSOR
            // substract c.m. force (acc) from tidal tensor force (acc)
            for (PS::S32 k=gpars.offset_tt; k<gpars.offset_orb; k++)  _sys[i+k].acc -= acc_cm;
#endif

            // After c.m. force used, it can be replaced by the averaged force on orbital particles
            acc_cm=PS::F64vec(0.0);
            PS::F64 m_ob_tot = 0.0;

            PS::S32 ob_start = i+gpars.offset_orb;
            for (PS::S32 k=ob_start; k<i_cm; k++) {
                acc_cm += _sys[k].mass*_sys[k].acc; 
                m_ob_tot += _sys[k].mass;
//#ifdef HARD_DEBUG
//                assert(((_sys[k].status)>>ID_PHASE_SHIFT)==-_sys[j_cm].id);
//#endif
            }
            acc_cm /= m_ob_tot;

#ifdef HARD_DEBUG
            assert(abs(m_ob_tot-_sys[i_cm].mass_bk)<1e-10);
#endif
        }
    }

//! soft force correction completely use tree neighbor search for all particles
/* @param[in,out] _sys: global particle system, acc is updated
   @param[in] _tree: tree for force
   @param[in] _adr_ptcl_artifical_start: start address of artifical particle in _sys
   @param[in] _rin: cutoff inner radius;
   @param[in] _rout: cutoff outer radius;
   @param[in] _n_split: artifical particle splitting number
   @param[in] _eps_sq: softing eps square
*/
    template <class Tsys, class Tpsoft, class Ttree, class Tepj>
    void correctForceWithCutoffTreeNeighborImp(Tsys& _sys, 
                                               Ttree& _tree, 
                                               const PS::S32 _adr_ptcl_artifical_start,
                                               const PS::F64 _rin,
                                               const PS::F64 _rout,
                                               const PS::S32 _n_split,
                                               const PS::F64 _eps_sq) { 
        // cutoff function parameter
        const PS::F64 r_oi_inv = 1.0/(_rout-_rin);
        const PS::F64 r_A = (_rout-_rin)/(_rout+_rin);
        // for artifical particle
        const PS::S32 n_tot = _sys.getNumberOfParticleLocal();

#pragma omp parallel for schedule(dynamic)
        for (int i=0; i<n_tot; i++) {
            correctForceWithCutoffTreeNeighborOneParticleImp<Tpsoft, Ttree, Tepj>(_sys[i], _tree, _rin, _rout, r_oi_inv, r_A, _eps_sq);
        }
        GroupPars gpars(_n_split);
#ifdef HARD_DEBUG
        assert((n_tot-_adr_ptcl_artifical_start)%gpars.n_ptcl_artifical==0);
#endif
#pragma omp parallel for schedule(dynamic)
        for (int i=_adr_ptcl_artifical_start; i<n_tot; i+=gpars.n_ptcl_artifical){
            PS::S32 i_cm = i + gpars.offset_cm;
            PS::F64vec& acc_cm = _sys[i_cm].acc;
#ifdef TIDAL_TENSOR
            // substract c.m. force (acc) from tidal tensor force (acc)
            for (PS::S32 k=gpars.offset_tt; k<gpars.offset_orb; k++)  _sys[i+k].acc -= acc_cm;
#endif

            // After c.m. force used, it can be replaced by the averaged force on orbital particles
            acc_cm=PS::F64vec(0.0);
            PS::F64 m_ob_tot = 0.0;

            PS::S32 ob_start = i+gpars.offset_orb;
            for (PS::S32 k=ob_start; k<i_cm; k++) {
                acc_cm += _sys[k].mass*_sys[k].acc; 
                m_ob_tot += _sys[k].mass;
//#ifdef HARD_DEBUG
//                assert(((_sys[k].status)>>ID_PHASE_SHIFT)==-_sys[j_cm].id);
//#endif
            }
            acc_cm /= m_ob_tot;

#ifdef HARD_DEBUG
            assert(abs(m_ob_tot-_sys[i_cm].mass_bk)<1e-10);
#endif
        }
    }

#ifdef ADJUST_GROUP_DEBUG
public:
#endif

    //! Calculate binary info
    /*!
      @param[out] _bininfo: binary information (outer most, tstep choose minimum)
      @param[in] _ptcl_list: ptcl index in _ptcl_origin
      @param[in] _n_ptcl: number of members
      @param[in] _ptcl_origin: original ptcl array
      @param[in] _dt_tree: tree time step for rsearch
    */
    void calcBinPars(Binary& _bininfo,
                     PS::S32* _ptcl_list,
                     const PS::S32 _n_ptcl,
                     PtclHard* _ptcl_origin,
                     const PS::F64 _dt_tree) {

        PtclTree<PtclHard> bins[_n_ptcl-1];
        keplerTreeGenerator(bins, _ptcl_list, _n_ptcl, _ptcl_origin, _dt_tree);

        _bininfo = *(Binary*)&(bins[_n_ptcl-2]);

        for (PS::S32 j=0; j<_n_ptcl-1; j++) 
            _bininfo.tstep = std::min(bins[j].tstep, _bininfo.tstep);
    }


    //! check group and adjust 
    /*! 
      1. check group change:
      First, check Aint situation
      if rmax(i)>rbin, split
      else rmax(i)<=rbin, check rmin(i)
          if perturb(rmax(i),rmin(i)) santisfy slowdown of subgroup >1, split
      
      second, check Hint situation  
      if rmin(i)<rbin,
         if i or j is Aint, check slowdownorg
            if kappaorg(i,j) <1, merge

      2. Integrated selected ptcl to _time_sys

      3. Add/remove Hint

      4. Initialize new ARC group, suppress unused one

      5. Initialize Hint new ptcl

      @param[in,out] _Hint: Hermite integrator class
      @param[in,out] _Aint: ARC integrator class
      @param[in,out] _ptcl_origin: originial ptcl array in current hard cluster
      @param[in,out] _group_member_index_origin: array of group member index in _ptcl_origin 
      @param[in,out] _n_group: number of groups
      @param[in,out] _single_index_origin: array of single index in _ptcl_origin
      @param[in,out] _n_single: number of singles
      @param[in] _time_sys: system time
      @param[in] _dt_max: time step maximum
      @param[in] _dt_min: time step minimum
      @param[in] _dt_tree: tree time step
     */
    template <class Tsoft>
    void adjustGroup(HermiteIntegrator<PtclHard>& _Hint,
                     ARCIntegrator<Ptcl, PtclH4, PtclForce>& _Aint,
                     PtclHard* _ptcl_origin,
                     const PS::F64 _time_sys,
                     const PS::F64 _dt_max,
                     const PS::F64 _dt_tree) {
        PS::S32 n_group=_Aint.getNGroups();
        PS::S32 n_hint = _Hint.getPtclN();
#ifdef HARD_DEBUG
        PS::S32 n_tot = _Hint.getPtclN();
        for (PS::S32 i=0; i<n_group; i++) n_tot += _Aint.getGroupN(i);
        n_tot -= n_group;
#endif
        const PS::F64 r_bin2 = r_bin_*r_bin_;
        // check Aint group to break
        PS::S32 break_group_list[n_group+1];
        PS::S32 break_split_index_list[n_group+1];
        PS::S32 n_group_break=_Aint.checkBreak(break_group_list, break_split_index_list, r_bin2);

        // check Hint to find new group
        PS::S32 hint_size=_Hint.getPtclN();
        PtclHard* new_group_member_adr_origin[hint_size];
        PS::S32 new_group_member_index_hint[hint_size];
        PS::S32 n_group_new = _Hint.checkNewGroup(new_group_member_index_hint, new_group_member_adr_origin, r_bin2, &_Aint);

        if(n_group_break==0&&n_group_new==0) return;

#ifdef ADJUST_GROUP_DEBUG
        for (PS::S32 i=0; i<n_group; i++) {
            PS::S32 n_member= _Aint.getGroupN(i);
            std::cout<<"Group i "<<i<<" N_member:"<<n_member<<std::endl;
            Ptcl** group_ptr= _Aint.getGroupPtclAdr(i);
            for (PS::S32 k=0; k<n_member; k++) {
                std::cout<<" "<<PS::S32((PtclHard*)group_ptr[k]-_ptcl_origin);
            }
            std::cout<<std::endl;
        }
        std::cout<<"N_Hint: "<<n_hint<<std::endl;
        PtclHard** ptcl_ptr = _Hint.getPtclAdr(); 
        for (PS::S32 i=0; i<n_hint; i++) {
            std::cout<<" "<<PS::S32(ptcl_ptr[i]-_ptcl_origin);
        }
        std::cout<<std::endl;
#endif
   
        // Integrate breaking ptcl to current time
        if(n_group_break>0) {
#ifdef ADJUST_GROUP_DEBUG
            for(PS::S32 i=0; i<n_group_break; i++) std::cout<<"Break group: "<<break_group_list[i]<<" group N: "<<_Aint.getGroupN(break_group_list[i])<<" split index: "<<break_split_index_list[i]<<std::endl;
#endif
            _Hint.integrateOneListNoPred(break_group_list, n_group_break, _time_sys, _dt_max, dt_min_hard_, &_Aint);
            //_Hint.writeBackPtcl(break_group_list, n_group_break);
            //_Aint.resolve(break_group_list, n_group_break);
        }

        // Integrate breaking ptcl to current time
        if(n_group_new>0) {
#ifdef ADJUST_GROUP_DEBUG
            std::cout<<"New group index in Hint:\n";
            for(PS::S32 i=0; i<n_group_new; i++) std::cout<<"group "<<i<<": "<<new_group_member_index_hint[2*i]<<" "<<new_group_member_index_hint[2*i+1]<<std::endl;
            std::cout<<"New group address in origin ptcl:\n";
            for(PS::S32 i=0; i<n_group_new; i++) std::cout<<"group "<<i<<": "<<PS::S32(new_group_member_adr_origin[2*i]-_ptcl_origin)<<" "<<PS::S32(new_group_member_adr_origin[2*i+1]-_ptcl_origin)<<std::endl;
#endif
            _Hint.integrateOneListNoPred(new_group_member_index_hint, 2*n_group_new, _time_sys, _dt_max, dt_min_hard_, &_Aint);
            _Hint.writeBackPtcl(new_group_member_index_hint, 2*n_group_new);
        }
        _Aint.resolve();
         
        PS::S32 hint_new_ptcl_index_origin[2*n_group+1]; // new ptcl index from _ptcl_origin
        PS::S32 hint_mod_ptcl_index[n_group+n_hint]; // modified ptcl index from Hint.ptcl_ (notice this is not _ptcl!)
        PS::S32 hint_del_ptcl_index[n_hint+1]; // del ptcl index from Hint.ptcl_. If it is c.m. only suppress the integration (remove index from time step list) 
        PS::S32 n_hint_new=0, n_hint_del=0, n_hint_mod=0;
        bool group_update_mask[n_group+1]={false};
        const PS::S32 n_group_old = n_group;

        if(n_group_break>0) {
            //break groups
            for (PS::S32 i=0; i<n_group_break; i++) {
                PS::S32 i_break = break_group_list[i];
                PS::S32 i_split = break_split_index_list[i];
                group_update_mask[i_break]=true;

                // get index in chain list
                PS::S32 n_member = _Aint.getGroupN(i_break);
//#ifdef HARD_DEBUG
//                assert(n_member==_group_member_index_origin[i_break].size());
//#endif
                // splitted group member index in Aint.chain
                Ptcl* first_adr_origin[n_member];
                Ptcl* second_adr_origin[n_member];
                PS::S32 n_first, n_second;
                _Aint.splitOneGroup(first_adr_origin, n_first, second_adr_origin, n_second, i_break, i_split);

                // current first_list and second_list are index in chain.lst, replace them to index in _ptcl
                PS::S32 first_index_origin[n_member], second_index_origin[n_member];
                if (n_first==1&&n_second>1) {
                    // if first is single, switch to second 
                    PS::S32 i_single = PS::S32((PtclHard*)first_adr_origin[0]-_ptcl_origin);
                    for (PS::S32 j=0; j<n_second; j++) first_index_origin[j] = PS::S32((PtclHard*)second_adr_origin[j] - _ptcl_origin);
                    n_first = n_second;
                    second_index_origin[0] = i_single;
                    n_second = 1;
                }
                else {
                    for (PS::S32 j=0; j<n_first;  j++)  first_index_origin[j] = PS::S32((PtclHard*)first_adr_origin[j] - _ptcl_origin);
                    for (PS::S32 j=0; j<n_second; j++) second_index_origin[j] = PS::S32((PtclHard*)second_adr_origin[j] - _ptcl_origin);
                }

#ifdef ADJUST_GROUP_DEBUG
                std::cout<<"Two group origin index:\n n_first: "<<n_first<<std::endl;
                for (PS::S32 k =0; k<n_first; k++) std::cout<<" "<<first_index_origin[k];
                std::cout<<"\nn_second: "<<n_second<<std::endl;
                for (PS::S32 k =0; k<n_second; k++) std::cout<<" "<<second_index_origin[k];
                std::cout<<std::endl;
#endif
                
                // 2 single case
                if(n_first==1&&n_second==1) {
                    // clear current group
                    _Aint.clearOneGroup(i_break);

                    hint_del_ptcl_index[n_hint_del++] = i_break; // used suppress dt

                    hint_new_ptcl_index_origin[n_hint_new++] = first_index_origin[0];
                    hint_new_ptcl_index_origin[n_hint_new++] = second_index_origin[0];

                    // clear group index list
                    //_group_member_index_origin[i_break].resizeNoInitialize(0);
                }
                // group - single
                else if(n_second==1) {
                    // recreate group
                    _Aint.clearOneGroup(i_break);
                    PS::S32 i_new=_Aint.addOneGroup(_ptcl_origin, first_index_origin, n_first, (Tsoft*)NULL, n_split_, _Hint.getPtcl(), _Hint.getForce());
#ifdef HARD_DEBUG
                    assert(i_new==i_break);
#endif
                    calcBinPars(_Aint.bininfo[i_new], first_index_origin, n_first, _ptcl_origin, _dt_tree);
                    
                    //_Aint.initial(i_break, _time_sys);
                    //_Aint.generateCMfromMembers(i_break);

                    if(i_new!=i_break) hint_del_ptcl_index[n_hint_del++] = i_break;
                    hint_mod_ptcl_index[n_hint_mod++] = i_new; // used to renew c.m.
                    hint_new_ptcl_index_origin[n_hint_new++] = second_index_origin[0];
                    
                    // update group member index 
                    //_group_member_index_origin[i_new].resizeNoInitialize(n_first);
                    //for (int j=0; j<n_first; j++) _group_member_index_origin[i_new][j] = first_index_origin[j];

                    n_group = _Aint.getNGroups();
                }
                // group - group
                else {
                    // recreate first group
                    _Aint.clearOneGroup(i_break);
                    PS::S32 i_new=_Aint.addOneGroup(_ptcl_origin, first_index_origin, n_first, (Tsoft*)NULL, n_split_, _Hint.getPtcl(), _Hint.getForce());
#ifdef HARD_DEBUG
                    assert(i_new==i_break);
#endif
                    calcBinPars(_Aint.bininfo[i_new], first_index_origin, n_first, _ptcl_origin, _dt_tree);

                    //_Aint.initial(i_break, _time_sys);
                    //_Aint.generateCMfromMembers(i_break);

                    // add new group
                    PS::S32 i_new2=_Aint.addOneGroup(_ptcl_origin, second_index_origin, n_second, (Tsoft*)NULL, n_split_, _Hint.getPtcl(), _Hint.getForce()); 

                    calcBinPars(_Aint.bininfo[i_new2], second_index_origin, n_second, _ptcl_origin, _dt_tree);

                    // copy soft perturbation parameters to new group
                    _Aint.copyParP2P(i_new2, i_break);
                    //_Aint.initial(i_group_new, _time_sys);
                    //_Aint.generateCMfromMembers(i_group_new);

                    if(i_new!=i_break) hint_del_ptcl_index[n_hint_del++] = i_break;
                    hint_mod_ptcl_index[n_hint_mod++] = i_new; // renew
                    hint_mod_ptcl_index[n_hint_mod++] = i_new2; // new c.m., put in mod list

                    // update group member index 
                    //_group_member_index_origin[i_new].resizeNoInitialize(n_first);
                    //for (int j=0; j<n_first; j++) _group_member_index_origin[i_new][j] = first_index_origin[j];
                    //_group_member_index_origin[i_new2].resizeNoInitialize(n_second);
                    //for (int j=0; j<n_second; j++) _group_member_index_origin[i_new2][j] = second_index_origin[j];

                    n_group = _Aint.getNGroups();
                }                
            }
        }
        
        if(n_group_new>0) {
            // new groups
            for (PS::S32 i=0; i<n_group_new; i++) {
                PS::S32 i1_hint = new_group_member_index_hint[i*2];
                PS::S32 i2_hint = new_group_member_index_hint[i*2+1];
                PtclHard* i1_adr_origin = new_group_member_adr_origin[i*2];
                PtclHard* i2_adr_origin = new_group_member_adr_origin[i*2+1];

                bool i1_group_flag=i1_hint<n_group_old;
                bool i2_group_flag=i2_hint<n_group_old;

                // if i1/2 are groups and need break, do nothing
                if(i1_group_flag&&group_update_mask[i1_hint]) continue;
                if(i2_group_flag&&group_update_mask[i2_hint]) continue;

                // switch i1/i2 if i1 is single and i2 is group
                if(!i1_group_flag&&i2_group_flag) {
                    i1_group_flag = true;
                    i2_group_flag = false;
                    PS::S32 i_tmp = i1_hint;
                    PtclHard* iadr_tmp = i1_adr_origin;
                    i1_hint = i2_hint;
                    i1_adr_origin = i2_adr_origin;
                    i2_hint = i_tmp;
                    i2_adr_origin = iadr_tmp;
                }

                // single - single
                if(!i1_group_flag&&!i2_group_flag) {
                    // create new group
                    // get ptcl index from _ptcl
                    PS::S32 new_index_origin[2]={PS::S32(i1_adr_origin-_ptcl_origin), PS::S32(i2_adr_origin-_ptcl_origin)};
                    
                    PS::S32 i_new = _Aint.addOneGroup(_ptcl_origin, new_index_origin, 2, (Tsoft*)NULL, n_split_, _Hint.getPtcl(), _Hint.getForce());

                    calcBinPars(_Aint.bininfo[i_new], new_index_origin, 2, _ptcl_origin, _dt_tree);
                    //_Aint.initial(_n_group, _time_sys);
                    //_Aint.generateCMfromMembers(_n_group);

                    hint_del_ptcl_index[n_hint_del++] = i1_hint;
                    hint_del_ptcl_index[n_hint_del++] = i2_hint;
                    hint_mod_ptcl_index[n_hint_mod++] = i_new; // new c.m.
                    
                    //_group_member_index_origin[i_new].resizeNoInitialize(2);
                    //for (int j=0; j<2; j++) _group_member_index_origin[i_new][j] = new_index_origin[j];

                    n_group = _Aint.getNGroups();
                }
                // group -single
                else if(!i2_group_flag) {
                    PS::S32 n_members_i1 = _Aint.getGroupN(i1_hint);
//#ifdef HARD_DEBUG
//                    assert(n_members_i1==_group_member_index_origin[i1_hint].size());
//#endif
                    PS::S32 new_group_member_index_origin[n_members_i1+1];
                    Ptcl** group_member_adr_origin_i1 = _Aint.getGroupPtclAdr(i1_hint);
                    for (PS::S32 j=0; j<n_members_i1; j++) {
                        new_group_member_index_origin[j] = PS::S32((PtclHard*)group_member_adr_origin_i1[j]-_ptcl_origin);
                    }
                    new_group_member_index_origin[n_members_i1] = PS::S32(i2_adr_origin-_ptcl_origin);
                    // increase group member index by 1
                    //_group_member_index_origin[i1_hint].increaseSize(1);
                    //_group_member_index_origin[i1_hint][n_members_i1] = i2_hint;
                    
                    // renew group i1
                    _Aint.clearOneGroup(i1_hint);
                    PS::S32 i_new=_Aint.addOneGroup(_ptcl_origin, new_group_member_index_origin, n_members_i1+1, (Tsoft*)NULL, n_split_, _Hint.getPtcl(), _Hint.getForce());
#ifdef HARD_DEBUG
                    assert(i_new==i1_hint);
#endif
                    calcBinPars(_Aint.bininfo[i_new], new_group_member_index_origin, n_members_i1+1, _ptcl_origin, _dt_tree);
                    //_Aint.initial(i1_hint, _time_sys);
                    //_Aint.generateCMfromMembers(i1_hint);
                    
                    hint_del_ptcl_index[n_hint_del++] = i2_hint;
                    hint_mod_ptcl_index[n_hint_mod++] = i_new; // renew i1
                    if(i_new!=i1_hint) {
                        hint_del_ptcl_index[n_hint_del++] = i1_hint; // suppress i1
                        n_group = _Aint.getNGroups();
                    }
                }
                // group - group
                else {
                    PS::S32 n_members_i1 = _Aint.getGroupN(i1_hint);
                    PS::S32 n_members_i2 = _Aint.getGroupN(i2_hint);
//#ifdef HARD_DEBUG
//                    assert(n_members_i1==_group_member_index_origin[i1_hint].size());
//                    assert(n_members_i2==_group_member_index_origin[i2_hint].size());
//#endif
                    PS::S32 new_group_member_index_origin[n_members_i1+n_members_i2];
                    Ptcl** group_member_adr_origin_i1 = _Aint.getGroupPtclAdr(i1_hint);
                    Ptcl** group_member_adr_origin_i2 = _Aint.getGroupPtclAdr(i2_hint);
                    for (PS::S32 j=0; j<n_members_i1; j++) {
                        new_group_member_index_origin[j] = PS::S32((PtclHard*)group_member_adr_origin_i1[j]-_ptcl_origin);
                    }
                    for (PS::S32 j=0; j<n_members_i2; j++) {
                        new_group_member_index_origin[j+n_members_i1] = PS::S32((PtclHard*)group_member_adr_origin_i2[j]-_ptcl_origin);
                    }

                    // merge second group member index list to first
                    //_group_member_index_origin[i1_hint].increaseSize(n_members_i2);
                    //for (PS::S32 j=0; j<n_members_i2; j++) 
                    //    _group_member_index_origin[i1_hint][n_members_i1+j] = _group_member_index_origin[i2_hint][j];
                    //_group_member_index_origin[i2_hint].resizeNoInitialize(0);
                    
                    // renew group i1
                    _Aint.clearOneGroup(i1_hint);
                    PS::S32 i_new=_Aint.addOneGroup(_ptcl_origin, new_group_member_index_origin, n_members_i1+n_members_i2, (Tsoft*)NULL, n_split_, _Hint.getPtcl(), _Hint.getForce());
#ifdef HARD_DEBUG
                    assert(i_new==i1_hint);
#endif

                    calcBinPars(_Aint.bininfo[i_new], new_group_member_index_origin, n_members_i1+n_members_i2, _ptcl_origin, _dt_tree);
                    //_Aint.initial(i1_hint, _time_sys);
                    //_Aint.generateCMfromMembers(i1_hint);

                    // clear group i2
                    _Aint.clearOneGroup(i2_hint);

                    hint_del_ptcl_index[n_hint_del++] = i2_hint; // suppress i2
                    hint_mod_ptcl_index[n_hint_mod++] = i_new; // renew i1
                    if(i_new!=i1_hint) {
                        hint_del_ptcl_index[n_hint_del++] = i1_hint; // suppress i1
                        n_group = _Aint.getNGroups();
                    }
                }
            }
        }

        // delete ptcl
        _Hint.removePtclList(hint_del_ptcl_index, n_hint_del, n_group_old, n_group_old);
        // renew c.m. ptcl
        for (PS::S32 i=0; i<n_hint_mod; i++) {
            PS::S32 i_mod_hint=hint_mod_ptcl_index[i];
            // initial new group
            _Aint.initialOneChain(i_mod_hint);
            Ptcl* ptcl_cm = _Aint.getCM(i_mod_hint);
            // update research
            ptcl_cm->calcRSearch(_dt_tree);
            ptcl_cm->id = -_Aint.getGroupPtcl(i_mod_hint)->id;
            // update c.m.
            if(i_mod_hint<n_group_old) _Hint.modOnePtcl(i_mod_hint, *ptcl_cm, _time_sys);
            // insert new c.m.
            else _Hint.insertOnePtcl(i_mod_hint, *ptcl_cm, n_group_old, _time_sys);
        }
        // add new ptcl
        _Hint.addPtclList(_ptcl_origin, hint_new_ptcl_index_origin, n_hint_new, n_group_old, _time_sys, true);

        // initial ARC system
        for (PS::S32 i=0; i<n_hint_mod; i++) {
             PS::S32 i_mod_hint=hint_mod_ptcl_index[i];
            _Hint.searchPerturberOne(i_mod_hint);
            _Aint.initialOneSys(i_mod_hint, _time_sys);
            _Aint.initialOneSlowDown(i_mod_hint, _time_sys+dt_limit_hard_, _Hint.getNbInfo(i_mod_hint).min_mass, sdfactor_, 1.0);
#ifdef HARD_DEBUG_PRINT
            fprintf(stderr,"New Group initial Slowdown parameters:");
            std::cerr<<"Aint k="<<i_mod_hint<<std::endl;
            _Aint.printSlowDown(std::cerr,i_mod_hint);
#endif            
        }

        // initial the new ptcl 
        _Hint.initial(NULL, n_hint_mod+n_hint_new, _time_sys, _dt_max, dt_min_hard_, &_Aint);

        // Update Aint perturber list
        for (PS::S32 i=0; i<n_group; i++) {
            if(!_Aint.getMask(i))
                _Aint.updatePertOneGroup(i, _Hint.getPtcl(), _Hint.getForce(), _Hint.getPertList(i), _Hint.getPertN(i));
        }

#ifdef ADJUST_GROUP_DEBUG
        std::cout<<"New group: n group="<<n_group<<"\n";
        for (PS::S32 i=0; i<n_group; i++) {
            PS::S32 n_member= _Aint.getGroupN(i);
            std::cout<<"Group i "<<i<<" N_member:"<<n_member<<std::endl;
            Ptcl** group_ptr= _Aint.getGroupPtclAdr(i);
            for (PS::S32 k=0; k<n_member; k++) {
                std::cout<<" "<<PS::S32((PtclHard*)group_ptr[k]-_ptcl_origin);
            }
            std::cout<<std::endl;
        }
        n_hint = _Hint.getPtclN();
        std::cout<<"New N_Hint: "<<n_hint<<std::endl;
        ptcl_ptr = _Hint.getPtclAdr(); 
        for (PS::S32 i=0; i<n_hint; i++) {
            std::cout<<" "<<PS::S32(ptcl_ptr[i]-_ptcl_origin);
        }
        std::cout<<std::endl;
#endif
        
#ifdef HARD_DEBUG
        // index consistent check
        PS::S32 n_count_check[n_tot] = {0};
        PtclH4* hint_ptcl = _Hint.getPtcl();
        for (PS::S32 k=0; k<n_group; k++) {
            if(!_Aint.getMask(k)) {
                const Ptcl* arc_cm = _Aint.getCM(k);
                assert(arc_cm->mass==hint_ptcl[k].mass);
                assert(arc_cm->pos.x==hint_ptcl[k].pos.x);
                const Ptcl* arc_ptcl = _Aint.getGroupPtcl(k);
                Ptcl** group_ptr= _Aint.getGroupPtclAdr(k);
                for (PS::S32 ig=0; ig<_Aint.getGroupN(k); ig++) {
                    PS::S32 adr= PS::S32((PtclHard*)group_ptr[ig] - _ptcl_origin);
                    assert(arc_ptcl[ig].id==_ptcl_origin[adr].id);
                    assert(adr<n_tot);
                    n_count_check[adr]++;
                }
            }
        }
        PtclHard** hint_ptcl_ptr = _Hint.getPtclAdr(); 
        for (PS::S32 k=n_group; k<_Hint.getPtclN(); k++) {
            PS::S32 adr=PS::S32(hint_ptcl_ptr[k]-_ptcl_origin);
            assert(hint_ptcl[k].id==_ptcl_origin[adr].id);
            assert(adr<n_tot);
            n_count_check[adr]++;
        }
        for (PS::S32 k=0; k<n_tot; k++) 
            assert(n_count_check[k]==1);

        // check perturber
        _Aint.checkPert();
#endif

    }

#ifdef ADJUST_GROUP_DEBUG
private:
#endif    
    //! Hard integration for clusters
    /* The local particle array are integrated. 
       No update of artifical particle pos and vel, eccept the artifical c.m. particle are kicked with acc. 
       The status of local particle in groups are set to 0 for first components if no tidal tensor memthod is used.
       @param[in,out] _ptcl_local: local particle in system_hard for integration
       @param[in] _n_ptcl: particle number in cluster
       @param[in,out] _ptcl_artifical: artifical particle array, c.m. are kicked 
       @param[in] _n_group: group number in cluster
       @param[in] _time_end: integration ending time (initial time is fixed to 0)
     */
    template <class Tsoft>
    void driveForMultiClusterImpl(PtclHard * _ptcl_local,
                                  const PS::S32 _n_ptcl,
                                  Tsoft* _ptcl_artifical,
                                  const PS::S32 _n_group,
                                  const PS::F64 _time_end) {
#ifdef HARD_CHECK_ENERGY
        std::map<PS::S32, PS::S32> N_count;  // counting number of particles in one cluster
        HardEnergy E0, E1;
        HardEnergy AE0, AE1;
        HardEnergy HE0, HE1;
        HardEnergy ESD0, ESD1;
#endif
#ifdef HARD_DEBUG_PROFILE
        N_count[_n_ptcl]++;
#endif
#ifdef HARD_DEBUG_DUMP
        PS::ReallocatableArray<PtclHard> ptcl_bk;
        ptcl_bk.reserve(_n_ptcl);
        for(int i=0; i<_n_ptcl; i++) ptcl_bk.pushBackNoCheck(_ptcl_local[i]);
#ifdef HARD_DEBUG_PRE_DUMP
        dump("hard_pre_dump",_time_end, ptcl_bk.getPointer(), _n_ptcl, _ptcl_artifical, _n_group*(2*n_split_+1), _n_group);
#endif
#endif
        PS::S32 nstepcount = 0;

#ifdef HARD_DEBUG
        if (_n_ptcl>400) {
            std::cerr<<"Large cluster, n_ptcl="<<_n_ptcl<<" n_group="<<_n_group<<std::endl;
            for (PS::S32 i=0; i<_n_ptcl; i++) {
                if(_ptcl_local[i].r_search>10*_ptcl_local[i].r_search_min) {
                    std::cerr<<"i = "<<i<<" ";
                    _ptcl_local[i].print(std::cerr);
                }
            }
        }
#endif

        //** suppressed, use address offset instead
        ///* The index:
        //   group_member_index: store all member index in _ptcl_local registered in ARCint, co-modified when ARC groups changes
        //   group_mask_int:  masker for all groups, if true, this group is not integrated
        //   single_index: store all member index in _ptcl_local registered in Hint, co-modified when Hint ptcl changes
        //   
        //   Notice the Hint first n_group are c.m. particles, all these c.m. particles should have same order as in ARCint all the time.
        // */
        //PS::ReallocatableArray<PS::S32> group_member_index[_n_ptcl];  // Group member index in _ptcl_local array
        //PS::S32 n_group = _n_group; // number of groups, including blend groups

        // prepare initial groups with artifical particles
        PS::S32 adr_first_ptcl[_n_group+1];
        PS::S32 adr_cm_ptcl[_n_group+1];
        PS::S32 n_group_offset[_n_group+1]; // ptcl member offset in _ptcl_local
        n_group_offset[0] = 0;

        //GroupPars gpars[_n_group](n_split_);
        GroupPars gpars[_n_group+1];
        for(int i=0; i<_n_group; i++) {
            gpars[i].init(n_split_);
            adr_first_ptcl[i] = i*gpars[i].n_ptcl_artifical;
            adr_cm_ptcl[i] = adr_first_ptcl[i]+gpars[i].offset_cm;
            gpars[i].getGroupIndex(&_ptcl_artifical[adr_first_ptcl[i]]);
            n_group_offset[i+1] = n_group_offset[i] + gpars[i].n_members;
#ifdef HARD_DEBUG
            assert(gpars[i].id == _ptcl_local[n_group_offset[i]].id);
#endif
            // initialize group_member_index
            //group_member_index[i].reserve(gpars[i].n_members+4);
            //group_member_index[i].resizeNoInitialize(gpars[i].n_members);
            //for (int j=0; j<gpars[i].n_members; j++) {
            //    group_member_index[i][j] = n_group_offset[i] + j;
            //}
        }
#ifdef HARD_DEBUG
        if(_n_group>0) {
            if(n_group_offset[_n_group]<_n_ptcl)
                assert(_ptcl_local[n_group_offset[_n_group]].status==0);
            assert(_ptcl_local[n_group_offset[_n_group]-1].status<0);
        }
#endif

        // single particle start index in _ptcl_local
        PS::S32 i_single_start = n_group_offset[_n_group];
        // number of single particles
        PS::S32 n_single_init = _n_ptcl - i_single_start;
#ifdef HARD_DEBUG
        assert(n_single_init>=0);
#endif

        // recover group member masses
        for(int i=0; i<i_single_start; i++) {
            //_ptcl_local[i].mass = _ptcl_local[i].mass_bk;
#ifdef HARD_DEBUG
            assert(_ptcl_local[i].status<0);
            assert(_ptcl_local[i].mass>0);
#endif
            _ptcl_local[i].mass_bk = 0.0;
        }

#ifndef TIDAL_TENSOR
        // In orbital fitting soft perturbation, the status is used to identify which component the member belong to
        for(int i=0; i<_n_group; i++) {
            // only first component is enough.
            for(int j=0; j<gpars[i].n_members_1st; j++)
                _ptcl_local[n_group_offset[i]+j].status = 0; 
        }
#endif

        // pre-process for c.m. particle,
        for(int i=0; i<_n_group; i++){
            PS::S32 icm = adr_cm_ptcl[i];
            // kick c.m. (not done in previous kick function to avoid multi-kick)
            //_ptcl_artifical[icm].vel += _ptcl_artifical[icm].acc * _time_end; (not do here to avoid half time step issue)
            // recover mass
            _ptcl_artifical[icm].mass = _ptcl_artifical[icm].mass_bk;
#ifdef HARD_DEBUG
            // check id 
            PS::S32 id_mem[2];
            id_mem[0] = _ptcl_local[n_group_offset[i]].id;
            id_mem[1] = _ptcl_local[n_group_offset[i]+gpars[i].n_members_1st].id;
            // id_offset unknown, try to substract id information via calculation between neighbor particles
            for (int j=0; j<gpars[i].n_ptcl_artifical-1; j+=2) {
                // first member
                PS::S32 id_offset_j1 = _ptcl_artifical[adr_first_ptcl[i]+j].id - j/2- id_mem[0]*n_split_;
                // second member
                PS::S32 id_offset_j2 = _ptcl_artifical[adr_first_ptcl[i]+j+1].id - j/2 - id_mem[1]*n_split_;
                assert(id_offset_j1==id_offset_j2);
            }

            // check whether c.m. pos. and vel. are consistent
            PS::F64 mass_cm_check=0.0;
            PS::F64vec vel_cm_check=PS::F64vec(0.0);
            PS::F64vec pos_cm_check=PS::F64vec(0.0);
            
            for(int j=0; j<gpars[i].n_members; j++) {
                PS::S32 k = n_group_offset[i]+j;
                mass_cm_check += _ptcl_local[k].mass;
                vel_cm_check +=  _ptcl_local[k].vel*_ptcl_local[k].mass;
                pos_cm_check +=  _ptcl_local[k].pos*_ptcl_local[k].mass;
            }
            vel_cm_check /= mass_cm_check;
            pos_cm_check /= mass_cm_check;

            assert(abs(mass_cm_check-_ptcl_artifical[icm].mass)<1e-10);
            PS::F64vec dvec = vel_cm_check-_ptcl_artifical[icm].vel;
            PS::F64vec dpos = pos_cm_check-_ptcl_artifical[icm].pos;
            assert(abs(dvec*dvec)<1e-20);
            assert(abs(dpos*dpos)<1e-20);
#endif

        }

        // Only one group with all particles in group
        if(_n_group==1&&n_single_init==0) {
            PS::S32 icm = adr_cm_ptcl[0];
 
 
            // create c.m. particles
            Ptcl pcm(_ptcl_artifical[icm]);
            PS::S32 iact = 0;
 
            ARCIntegrator<Ptcl, PtclH4, PtclForce> Aint(ARC_control_soft_, Int_pars_);
#ifdef ARC_SYM
            Aint.step_count_limit = arc_step_count_limit;
#endif
            Aint.reserveARMem(1);
            Aint.reservePertMem(1,1);
            gpars[0].getBinPars(Aint.bininfo[0], _ptcl_artifical);
#ifdef HARD_DEBUG_PRINT
            std::cerr<<"Hard: one group, n="<<_n_ptcl<<std::endl;
            Aint.bininfo[0].print(std::cerr,20,true);
#endif
#ifdef TIDAL_TENSOR
            PS::S32 i_soft_pert_offset = gpars[0].offset_tt;
#else
            PS::S32 i_soft_pert_offset = gpars[0].offset_orb;
#endif
            Aint.addOneGroup(_ptcl_local, NULL, gpars[0].n_members, &_ptcl_artifical[i_soft_pert_offset], n_split_);
            Aint.updateCM(&pcm, &iact, 1);
 
            Aint.initialOneChain(0);
            Aint.initialOneSys(0, 0.0);
            Aint.initialOneSlowDownUnPert(0, _time_end, sdfactor_, 1.0);
 
#ifdef ARC_SYM_SD_PERIOD
            PS::S32 kp=0;
            Aint.adjustSlowDownPeriod(_time_end, &kp);
#else
            Aint.adjustSlowDown(_time_end);
#endif
 
#ifdef HARD_CHECK_ENERGY
            Aint.EnergyRecord(AE0);
#endif 
 
#ifdef ARC_SYM
#ifdef ARC_SYM_SD_PERIOD
            if(group_n==2&&kp>0) nstepcount +=Aint.integrateOneStepSymTwo(0, _time_end, kp);
            else nstepcount +=Aint.integrateOneStepSym(0, _time_end, dt_limit_hard_);
#else
            nstepcount +=Aint.integrateOneStepSym(0, _time_end, dt_limit_hard_);
#endif
#else 
            nstepcount +=Aint.integrateOneStepExt(0, _time_end, dt_limit_hard_);
#endif
            
            pcm.pos += pcm.vel * _time_end;
 
            Aint.updateCM(&pcm, &iact, 1);
            Aint.resolve();
#ifdef HARD_CHECK_ENERGY
            Aint.EnergyRecord(AE1);
            PS::F64 dEtot = AE1.kin+AE1.pot+AE1.tot-AE0.kin-AE0.pot-AE0.tot;
            hard_dE += dEtot;
#ifdef HARD_DEBUG_PRINT
            fprintf(stderr,"Slowdown factor = %e\n", Aint.getSlowDown(0));
            fprintf(stderr,"ARC Energy: init =%e, end =%e, diff =%e, error = %e\n", 
                    AE0.kin+AE0.pot, AE1.kin+AE1.pot, AE1.kin+AE1.pot-AE0.kin-AE0.pot, (AE1.kin+AE1.pot+AE1.tot-AE0.kin-AE0.pot-AE0.tot)/AE0.tot);
#endif
#ifdef HARD_DEBUG_DUMP
            if(fabs(dEtot)>hard_dE_limit) {
                std::cerr<<"Hard energy significant: "<<dEtot<<std::endl;
                std::cerr<<"Dump data:"<<std::endl;
                dump("hard_dump",_time_end,  ptcl_bk.getPointer(), _n_ptcl, _ptcl_artifical, _n_group*(2*n_split_+1), _n_group);
                abort();
            }
#endif
#endif
#ifdef ARC_DEBUG_PRINT
#ifdef ARC_SYM_SD_PERIOD
            Aint.info_print(std::cerr, ARC_n_groups, 1, _n_ptcl, 0, dt_limit_hard_,kp);
#else
            Aint.info_print(std::cerr, ARC_n_groups, 1, _n_ptcl, 0, dt_limit_hard_, 0, 100000);
//            if (dump_flag) {
//                dump("hard_dump",_time_end,  ptcl_bk.getPointer(), _n_ptcl, _ptcl_artifical, _n_group*(2*n_split_+1), _n_group);
//                abort();
//            }
#endif
#endif
#ifdef PROFILE
            //ARC_substep_sum += Aint.getNsubstep();
            ARC_substep_sum += nstepcount;
            ARC_n_groups += 1;
#endif
        }
        else {
            // Single index in _ptcl_local array
            //PS::S32 single_index[_n_ptcl]; 
            //PS::S32 n_single = 0;

            // integration -----------------------------
            HermiteIntegrator<PtclHard> Hint;
            Hint.setParams(eta_s_, Int_pars_.rin, Int_pars_.rout, Int_pars_.eps2, _n_ptcl);

            Hint.reserveMem(_n_ptcl);

            // add c.m.
            Hint.addPtclList(_ptcl_artifical, adr_cm_ptcl, _n_group, 0, 0.0, false);

            // reset n_single
            //n_single = 0;

            // add single
            Hint.addPtclList(&_ptcl_local[i_single_start], NULL, n_single_init, 0, 0.0, false);
            //for (int i=0; i<n_single; i++) single_index[i] += i_single_start;

            PS::F64 time_sys=0.0, time_now;
#ifdef FIX_STEP_DEBUG
            PS::F64 dt_limit = dt_limit_hard_;
#else
            PS::F64 dt_limit = calcDtLimit(time_sys, dt_limit_hard_, dt_min_hard_);
#endif
            
            //PS::S32 group_act_n = 0;
            //PS::ReallocatableArray<PS::S32> group_act_list; //active group_list act adr
            // ReallocatableArray<PS::S32> group_list;     //group.adr list
            // ReallocatableArray<PS::S32> status;      //ptcl -> group.adr [non cm is -1] (value of Ptcl.status)
            // ReallocatableArray<PS::S32> status_map;  //ptcl -> group_list index [non cm is -1]
            // ReallocatableArray<PS::S32> adr_cm;         //group_list index -> ptcl.cm
            // group.findGroups(group_list, status, status_map,  adr_cm, group_act_n, _ptcl_local, _n_ptcl);

#ifdef HARD_DEBUG
#ifdef HARD_DEBUG_PRINT
            std::cerr<<"Hard: mix, n="<<_n_ptcl<<" n_group= "<<_n_group<<std::endl;
            for(int i=0; i<_n_group; i++) std::cerr<<"i_group= "<<i<<" n_members= "<<gpars[i].n_members<<std::endl;
            std::FILE* fp = std::fopen("hard_debug_ptcl","w");
            if (fp==NULL) {
                std::cerr<<"Error: file hard_debug_ptcl cannot be open!\n";
                abort();
            }
#endif
//            assert(n_hint<ARRAY_ALLOW_LIMIT);
#endif        
            //group_act_list.resizeNoInitialize(n_hint);

            
            // Initial Aint
            ARCIntegrator<Ptcl, PtclH4, PtclForce> Aint(ARC_control_pert_, Int_pars_);
#ifdef ARC_SYM
            Aint.step_count_limit = arc_step_count_limit;
#endif
            // Reserve memory space
            Aint.reserveARMem(_n_group+1);
            // first particles in Hint.Ptcl are c.m. thus add 1
            Aint.reservePertMem(_n_group+1,_n_ptcl+1);

            // Obtain binary parameters
            //PS::F64 apo_bin[_n_group+1];
            for (int i=0; i<_n_group; i++) {
                gpars[i].getBinPars(Aint.bininfo[i], &_ptcl_artifical[adr_first_ptcl[i]]);
#ifdef HARD_DEBUG_PRINT
                std::cerr<<"Group i:"<<i<<" ";
                Aint.bininfo[i].print(std::cerr,20,true);
#endif
                //auto &bini= Aint.bininfo[i];
                // In the new tree search way, this is not necessary
                /* Notice in the neighbor search, the resolved members are used.
                   The two members in different binaries can find each other as neighbors, but the c.m. particle may not find another c.m. 
                   To avoid no perturbers issues, the rsearch should add the maximum distance of components in other binaries (apo-center distance).
                */
                //apo_bin[i] = bini.semi*(bini.ecc+1.0); 
            }            
            Hint.searchPerturber(_n_group);

            for (int i=0; i<_n_group; i++) {
//#ifdef HARD_DEBUG
//                    assert(Hint.getPertN(i)>0);
//#endif
#ifdef TIDAL_TENSOR
                PS::S32 i_soft_pert_offset = adr_first_ptcl[i]+gpars[i].offset_tt;
#else
                PS::S32 i_soft_pert_offset = adr_first_ptcl[i]+gpars[i].offset_orb;
#endif
                
                Aint.addOneGroup(&_ptcl_local[n_group_offset[i]], NULL, gpars[i].n_members, &_ptcl_artifical[i_soft_pert_offset], n_split_, Hint.getPtcl(), Hint.getForce(), Hint.getPertList(i), Hint.getPertN(i));
                Aint.initialOneChain(i);
                Aint.initialOneSys(i,time_sys);
                Aint.initialOneSlowDown(i, dt_limit, Hint.getNbInfo(i).min_mass, sdfactor_, 1.0);
            }

#ifdef HARD_CHECK_ENERGY
            CalcEnergyHardFull(_ptcl_local, _n_ptcl, E0, AE0, HE0, ESD0, Hint, Aint);
#endif
#ifdef HARD_DEBUG_PRINT
            fprintf(stderr,"Initial Slowdown parameters:");
            for(int k=0; k<Aint.getNGroups(); k++) {
                std::cerr<<"Aint k="<<k<<std::endl;
                Aint.printSlowDown(std::cerr,k);
            }
#endif

            Hint.shiftToCM(); // shift ptcl to c.m. frame
            Hint.calcA0offset();
            bool fail_flag=Hint.initial(NULL, n_single_init + _n_group, 0.0, dt_limit, dt_min_hard_, &Aint, false);
            Hint.SortAndSelectIp();

//            for (int i=0; i<_n_group; i++) {
//                Aint.updateOneSlowDown(i, Hint.getOneTime(i), Hint.getOneDt(i), dt_limit, Hint.getNbInfo(i).min_mass, -1);
//            }
// 
//#ifdef HARD_DEBUG_PRINT
//            fprintf(stderr,"First Slowdown parameters:");
//            for(int k=0; k<Aint.getNGroups(); k++) {
//                std::cerr<<"Aint k="<<k<<std::endl;
//                Aint.printSlowDown(std::cerr,k);
//            }
//#endif
            if(fail_flag) {
#ifdef HARD_DEBUG_DUMP
                std::cerr<<"Dump hard data. tend="<<_time_end<<" _n_ptcl="<<_n_ptcl<<"\n";
                dump("hard_dump",_time_end, ptcl_bk.getPointer(), _n_ptcl, _ptcl_artifical, _n_group*gpars[0].n_ptcl_artifical, _n_group);
                abort();
#endif
            }

//#ifdef HARD_CHECK_ENERGY
//                PS::ReallocatableArray<PS::F64> slowdownrecord;
#ifdef HARD_DEBUG
            assert(_n_group<ARRAY_ALLOW_LIMIT);
#endif        
//                slowdownrecord.resizeNoInitialize(n_group);
//#endif

            while(time_sys<_time_end) {
                time_now = time_sys;
                time_sys = Hint.getNextTime();
#ifdef FIX_STEP_DEBUG
                dt_limit = dt_limit_hard_;
#else
                dt_limit = calcDtLimit(time_sys, dt_limit_hard_, dt_min_hard_);
#endif

#ifdef HARD_DEBUG
                assert(time_sys>time_now);
#endif
                PS::F64 dt_h = time_sys-time_now;
                for (int k=0; k<_n_group; k++) {
                    if(!Aint.getMask(k))
                        Aint.updateOneSlowDown(k, Hint.getOneTime(k), Hint.getOneDt(k), dt_limit_hard_, Hint.getNbInfo(k).min_mass);
                }
//#ifdef HARD_CHECK_ENERGY
//                    for(int k=0; k<n_group; k++) {
//                        slowdownrecord[k] = std::max(slowdownrecord[k], Aint.getSlowDown(k));
//                        assert(Aint.getSlowDown(k)>=1.0);
//                    }
//#endif
                nstepcount +=Aint.integrateOneStepList(time_sys, std::min(dt_limit,dt_h));
                fail_flag = Hint.integrateOneStepAct(time_sys,dt_limit,dt_min_hard_, &Aint);
                adjustGroup<Tsoft>(Hint, Aint, _ptcl_local, time_sys, dt_limit, _time_end);

                
                if(fail_flag) {
#ifdef HARD_DEBUG_DUMP
                    std::cerr<<"Dump hard data. tend="<<_time_end<<" _n_ptcl="<<_n_ptcl<<"\n";
                    dump("hard_dump",_time_end, ptcl_bk.getPointer(), _n_ptcl, _ptcl_artifical, _n_group*gpars[0].n_ptcl_artifical, _n_group);
                    abort();
#endif
                }
                Hint.SortAndSelectIp();
#ifdef HARD_DEBUG
                // check time step list
                Hint.checkAdrList(Aint);
#endif


#ifdef HARD_DEBUG_PRINT
                fprintf(stderr,"Time = %g, dt = %g, nstep_ARC = %d \n",time_sys, dt_h, nstepcount);
                fprintf(stderr,"Slowdown parameters:");
                for(int k=0; k<Aint.getNGroups(); k++) {
                    std::cerr<<"Aint k="<<k<<std::endl;
                    Aint.printSlowDown(std::cerr,k);
                }
                // fprintf(stderr,"Slowdownfactor: ");
                // for(int k=0; k<_n_group; k++) fprintf(stderr,"%d: %f; ",k,Aint.getSlowDown(k));
                // fprintf(stderr,"\n");
                Hint.printStepHist();
                fprintf(fp,"%25.14e ",time_sys);
                Aint.writePtcl<Ptcl>(fp);
                Hint.writePtcl(fp,Aint.getNGroups());
                fprintf(fp,"\n");
#endif
            }
        
            Hint.moveCM(_time_end);
            Hint.shiftBackCM();
            Aint.updateCM(Hint.getPtcl());
            Aint.resolve();
            Hint.writeBackPtcl(Aint.getNGroups());

#ifdef ARC_DEBUG_PRINT
            Aint.info_print(std::cerr, ARC_n_groups, _n_group, _n_ptcl, Hint.getPtclN(), dt_limit_hard_,0);
#endif
#ifdef HARD_CHECK_ENERGY
            CalcEnergyHardFull(_ptcl_local, _n_ptcl, E1, AE1, HE1, ESD1, Hint, Aint);
            PS::F64 dEtot = E1.tot - E0.tot;
            hard_dE += dEtot;
            hard_dESD += ESD1.tot - ESD0.tot;
#ifdef HARD_DEBUG_PRINT
            fclose(fp);
            fprintf(stderr,"Slowdown parameters:");
            for(int k=0; k<_n_group; k++) {
                std::cerr<<"Aint k="<<k<<std::endl;
                Aint.printSlowDown(std::cerr,k);
            }
            fprintf(stderr,"\n");
            fprintf(stderr,"H4  Energy: init =%e, end =%e, diff =%e, kini =%e kinf =%e poti =%e potf =%e\nARC Energy: init =%e, end =%e, diff =%e, error = %e\nTot Energy: init =%e, end =%e, diff =%e, kin =%e pot =%e, Tot-H4-ARC =%e\nTSD Energy: init =%e, end =%e, diff =%e, kin =%e pot =%e\n", 
                    HE0.tot, HE1.tot, HE1.tot-HE0.tot, HE0.kin, HE1.kin, HE0.pot, HE1.pot, 
                    AE0.kin+AE0.pot, AE1.kin+AE1.pot, AE1.kin+AE1.pot-AE0.kin-AE0.pot, (AE1.kin+AE1.pot+AE1.tot-AE0.kin-AE0.pot-AE0.tot)/AE0.tot,
                    E0.tot, E1.tot, E1.tot-E0.tot, E1.kin, E1.pot, E1.tot-HE1.tot-AE1.kin-AE1.pot,
                    ESD0.tot, ESD1.tot, ESD1.tot-ESD0.tot, ESD1.kin, ESD1.pot);
            Hint.printStepHist();
#endif
#ifdef HARD_DEBUG_DUMP
            if(fabs(dEtot)>hard_dE_limit) {
                std::cerr<<"Hard energy significant: "<<dEtot<<std::endl;
                std::cerr<<"Dump data:"<<std::endl;
                dump("hard_dump",_time_end, ptcl_bk.getPointer(), _n_ptcl, _ptcl_artifical, _n_group*gpars[0].n_ptcl_artifical, _n_group);
                abort();
            }
#endif
#endif
#ifdef PROFILE
            //ARC_substep_sum += Aint.getNsubstep();
            ARC_substep_sum += nstepcount;
            ARC_n_groups += _n_group;
#endif
        }
//        }            
//        else { // no group
//            // initial single_index
//            n_single = _n_ptcl;
//            for (int i=0; i<_n_ptcl; i++) {
//                single_index[i] = i;
//            }
//            n_group = 0;
// 
//            HermiteIntegrator Hint;
//            Hint.setParams(eta_s_, Int_pars_.rin, Int_pars_.rout, Int_pars_.eps2);
//            Hint.resizeArray(_n_ptcl);
// 
//            // Null pointer for arguments of Hint
//            ARCIntegrator<Ptcl, PtclH4, PtclForce> *Aint_null=NULL;
// 
//            // add single
//            Hint.setPtcl(_ptcl_local, _n_ptcl);
// 
//#ifdef HARD_DEBUG
//            assert(_n_ptcl>1);
//#ifdef HARD_DEBUG_PRINT
//            std::cerr<<"Hard: hermite, n="<<_n_ptcl<<std::endl;
//#endif
//            for(int i=0; i<_n_ptcl; i++) {
//                assert(_ptcl_local[i].status==0);
//                assert(_ptcl_local[i].mass>0);
//                assert(_ptcl_local[i].id>0);
//            }
//#endif
// 
//            PS::F64* dr_search=NULL;
//            Hint.searchPerturber(dr_search,0,_n_ptcl);
//            
//            PS::F64 time_sys=0.0;
//#ifdef HARD_DEBUG
//            PS::F64 time_now;
//#endif
// 
//#ifdef FIX_STEP_DEBUG
//            PS::F64 dt_limit = dt_limit_hard_;
//#else
//            PS::F64 dt_limit = calcDtLimit(time_sys, dt_limit_hard_, dt_min_hard_);
//#endif
// 
//            //PS::S32 group_act_n = 0;
//            //PS::ReallocatableArray<PS::S32> group_act_list; //active group_list act adr
// 
//#ifdef HARD_DEBUG
//            assert(_n_ptcl<ARRAY_ALLOW_LIMIT);
//#endif        
//            //group_act_list.resizeNoInitialize(_n_ptcl);
// 
//#ifdef HARD_CHECK_ENERGY
//            // calculate initial energy
//            Hint.CalcEnergy(HE0);
//#endif
// 
//            bool fail_flag=Hint.initialize(dt_limit, dt_min_hard_, 0, Aint_null);
// 
//            if(fail_flag) {
//#ifdef HARD_DEBUG_DUMP
//                std::cerr<<"Dump hard data. tend="<<_time_end<<" _n_ptcl="<<_n_ptcl<<"\n";
//                dump("hard_dump",_time_end, ptcl_bk.getPointer(), _n_ptcl, _ptcl_artifical, 0, 0);
//                abort();
//#endif
//            }
// 
//            while(time_sys<_time_end) {
//#ifdef HARD_DEBUG
//                time_now = time_sys;
//#endif
//                time_sys = Hint.getNextTime();
//#ifdef FIX_STEP_DEBUG
//                dt_limit = dt_limit_hard_;
//#else
//                dt_limit = calcDtLimit(time_sys, dt_limit_hard_, dt_min_hard_);
//#endif
// 
//#ifdef HARD_DEBUG
//                assert(time_sys>time_now);
//#endif
//                fail_flag = Hint.integrateOneStep(time_sys,dt_limit,dt_min_hard_,true, Aint_null);
//                if(fail_flag) {
//#ifdef HARD_DEBUG_DUMP
//                    std::cerr<<"Dump hard data. tend="<<_time_end<<" _n_ptcl="<<_n_ptcl<<"\n";
//                    dump("hard_dump",_time_end, ptcl_bk.getPointer(), _n_ptcl, _ptcl_artifical, 0, 0);
//                    abort();
//#endif
//                }
//                
//                Hint.SortAndSelectIp();
//            }
//            
//            Hint.moveCM(_time_end);
//            Hint.shiftBackCM();
// 
//            Hint.writeBackPtcl(_ptcl_local, _n_ptcl, single_index, 0);
// 
//#ifdef HARD_CHECK_ENERGY
//            Hint.CalcEnergy(HE1);
//            PS::F64 dEtot = HE1.tot-HE0.tot;
//            hard_dE += dEtot;
//#ifdef HARD_DEBUG_PRINT
//            fprintf(stderr,"H4  Energy: init =%e, end =%e, diff =%e, kini =%e kinf =%e poti =%e potf =%e\n", 
//                    HE0.tot, HE1.tot, HE1.tot-HE0.tot, HE0.kin, HE1.kin, HE0.pot, HE1.pot);
//            Hint.printStepHist();
//#endif
//#ifdef HARD_DEBUG_DUMP
//            if(fabs(dEtot)>hard_dE_limit) {
//                std::cerr<<"Hard energy significant: "<<dEtot<<std::endl;
//                std::cerr<<"Dump data:"<<std::endl;
//                dump("hard_dump",_time_end,  ptcl_bk.getPointer(), _n_ptcl, _ptcl_artifical, 0, 0);
//                abort();
//            }
//#endif
//#endif
//        }

    }

public:

    SystemHard(){
#ifdef HARD_DEBUG_PROFILE
        for(PS::S32 i=0;i<20;i++) N_count[i]=0;
#endif
#ifdef PROFILE
        ARC_substep_sum = 0;
#endif
#ifdef HARD_CHECK_ENERGY
        hard_dE = hard_dESD = 0;
#endif
        //        PS::S32 n_threads = PS::Comm::getNumberOfThread();
    }

    /// start set Chainpars (L.Wang)
    ///
    void setARCParam(const PS::F64 energy_error=1e-10, const PS::F64 dterr_pert=1e-6, const PS::F64 dterr_soft=1e-3, const PS::F64 dtmin=1e-24
#ifndef ARC_SYM
                     ,const PS::S32 exp_method=1, const PS::S32 exp_itermax=20, const PS::S32 den_intpmax=20, const PS::S32 exp_fix_iter=0
#endif
        ) {
#ifdef HARD_DEBUG_DEEP_CHECK
        ARC_control_pert_.setA(Newtonian_AW<Ptcl,ARC_pert_pars>,Newtonian_extA_test<Ptcl,PtclH4*,PtclForce*,ARC_pert_pars>,Newtonian_timescale<ARC_pert_pars>);
        ARC_control_soft_.setA(Newtonian_AW<Ptcl,ARC_pert_pars>,Newtonian_extA_test<Ptcl,PtclH4*,PtclForce*,ARC_pert_pars>,Newtonian_timescale<ARC_pert_pars>);
#else
        ARC_control_pert_.setA(Newtonian_AW<Ptcl,ARC_pert_pars>,Newtonian_extA_pert<Ptcl,PtclH4*,PtclForce*,ARC_pert_pars>,Newtonian_timescale<ARC_pert_pars>);
        ARC_control_soft_.setA(Newtonian_AW<Ptcl,ARC_pert_pars>,Newtonian_extA_soft<Ptcl,PtclH4*,PtclForce*,ARC_pert_pars>,Newtonian_timescale<ARC_pert_pars>);
#endif
        ARC_control_pert_.setabg(0,1,0);
        ARC_control_soft_.setabg(0,1,0);
        
        ARC_control_pert_.setErr(energy_error,dtmin,dterr_pert);
        ARC_control_soft_.setErr(energy_error,dtmin,dterr_soft);
#ifdef ARC_SYM
        ARC_control_pert_.setSymOrder(-6);
        ARC_control_soft_.setSymOrder(-6);
#else
        ARC_control_pert_.setIterSeq(exp_itermax,3,den_intpmax);
        ARC_control_soft_.setIterSeq(exp_itermax,3,den_intpmax);

        ARC_control_pert_.setIntp(exp_method);
        ARC_control_soft_.setIntp(exp_method);

        ARC_control_pert_.setIterConst((bool)exp_fix_iter);
        ARC_control_soft_.setIterConst((bool)exp_fix_iter);

        ARC_control_pert_.setAutoStep(3);
        ARC_control_soft_.setAutoStep(3);
#endif
    }
    /// end set Chainpars (L.Wang)

    void initializeForOneCluster(const PS::S32 n){
#ifdef HARD_DEBUG
        assert(n<ARRAY_ALLOW_LIMIT);
#endif        
        ptcl_hard_.resizeNoInitialize(n);
    }

    ////////////////////////
    // for NON-ISOLATED CLUSTER
    template<class Tsys, class Tptcl, class Tmediator>
    void setPtclForConnectedCluster(const Tsys & sys,
                                   const PS::ReallocatableArray<Tmediator> & med,
                                   const PS::ReallocatableArray<Tptcl> & ptcl_recv){
        ptcl_hard_.clearSize();
        n_ptcl_in_cluster_.clearSize(); // clear befor break this function
        for(PS::S32 i=0; i<med.size(); i++){
            if(med[i].adr_sys_ < 0) continue;
            if(med[i].rank_send_ != PS::Comm::getRank()) continue;
            const auto & p = sys[med[i].adr_sys_];
            ptcl_hard_.push_back(PtclHard(p, med[i].id_cluster_, med[i].adr_sys_));
#ifdef HARD_DEBUG
            assert(med[i].adr_sys_<sys.getNumberOfParticleLocal());
            if(p.id<0&&p.status<0) {
                std::cerr<<"Error: ghost particle is selected! i="<<i<<"; med[i].adr_sys="<<med[i].adr_sys_<<std::endl;
                abort();
            }
#endif
        }

        for(PS::S32 i=0; i<ptcl_recv.size(); i++){
            const Tptcl & p = ptcl_recv[i];
            ptcl_hard_.push_back(PtclHard(p, p.id_cluster, -(i+1)));
#ifdef HARD_DEBUG
            if(p.id<0&&p.status<0) {
                std::cerr<<"Error: receive ghost particle! i="<<i<<std::endl;
                abort();
            }
#endif
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
#ifdef HARD_DEBUG
        assert(n_cluster<ARRAY_ALLOW_LIMIT);
#endif        
        n_ptcl_in_cluster_disp_.resizeNoInitialize(n_cluster+1);
        n_ptcl_in_cluster_disp_[0] = 0;
        for(PS::S32 i=0; i<n_cluster; i++){
#ifdef HARD_DEBUG
            assert(n_ptcl_in_cluster_[i]>1);
#endif
            n_ptcl_in_cluster_disp_[i+1] = n_ptcl_in_cluster_disp_[i] + n_ptcl_in_cluster_[i];
        }
    }


    // for NON-ISOLATED CLUSTER
    ////////////////////////

    PS::S32 getGroupPtclRemoteN() const{
        return n_group_member_remote_;
    }

    PS::ReallocatableArray<PtclHard> & getPtcl() {
        return ptcl_hard_;
    }

    PS::S32 getNCluster() const{
        return n_ptcl_in_cluster_.size();
    }

    PS::S32* getClusterNList(const std::size_t i=0) const{
        return n_ptcl_in_cluster_.getPointer(i);
    }

    PS::S32* getClusterNOffset(const std::size_t i=0) const{
        return n_ptcl_in_cluster_disp_.getPointer(i);
    }

    PS::S32* getGroupNList(const std::size_t i=0) const{
        return n_group_in_cluster_.getPointer(i);
    }

    PS::S32* getGroupNOffset(const std::size_t i=0) const{
        return n_group_in_cluster_offset_.getPointer(i);
    }

    PS::S32* getAdrPtclArtFirstList(const std::size_t i=0) const{
        return adr_first_ptcl_arti_in_cluster_.getPointer(i);
    }

    void setTimeOrigin(const PS::F64 _time_origin){
        time_origin_ = _time_origin;
    }

    void setParam(const PS::F64 _rbin,
                  const PS::F64 _rout,
                  const PS::F64 _rin,
                  const PS::F64 _eps,
                  const PS::F64 _dt_limit_hard,
                  const PS::F64 _dt_min_hard,
                  const PS::F64 _eta,
                  const PS::F64 _time_origin,
                  const PS::F64 _sd_factor,
                  // const PS::F64 _gmin,
                  // const PS::F64 _m_avarage,
                  const PS::S64 _id_offset,
                  const PS::S32 _n_split){
        /// Set chain pars (L.Wang)
		Int_pars_.rin  = _rin;
        Int_pars_.rout = _rout;
        Int_pars_.r_oi_inv = 1.0/(_rout-_rin);
        Int_pars_.r_A      = (_rout-_rin)/(_rout+_rin);
        Int_pars_.pot_off  = (1.0+Int_pars_.r_A)/_rout;
        Int_pars_.eps2  = _eps*_eps;
        /// Set chain pars (L.Wang)        
        dt_limit_hard_ = _dt_limit_hard;
        dt_min_hard_   = _dt_min_hard;
        eta_s_ = _eta*_eta;
        sdfactor_ = _sd_factor;
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

#ifdef HARD_CHECK_ENERGY
    //! check energy
    /* @param[in] _ptcl: particle data array
       @param[out] _energy: energy class
       @param[in] _n_ptcl: particle number
    */
    template<class Teng>
    void CalcEnergyHard(PtclHard* _ptcl, Teng & _energy,  const PS::S32 _n_ptcl) {
        _energy.kin = _energy.pot = _energy.tot = 0.0;
        for(PS::S32 i=0; i<_n_ptcl; i++){
            PtclHard* pi = &_ptcl[i];
            _energy.kin += 0.5 * pi->mass * pi->vel * pi->vel;

            for(PS::S32 j=i+1; j<_n_ptcl; j++){
                PtclHard* pj = &_ptcl[j];
                PS::F64vec rij = pi->pos - pj->pos;
                PS::F64 dr = sqrt(rij*rij + Int_pars_.eps2);
#ifdef INTEGRATED_CUTOFF_FUNCTION
                _energy.pot -= pj->mass*pi->mass/dr*(1.0 - CalcW(dr/Int_pars_.rout, Int_pars_.rin/Int_pars_.rout));  
#else
                if(dr<Int_pars_.rout) _energy.pot -= pj->mass*pi->mass*(1.0/dr*cutoff_pot(dr, Int_pars_.r_oi_inv, Int_pars_.r_A, Int_pars_.rin) - Int_pars_.pot_off);
#endif
            }
        }
        _energy.tot = _energy.kin + _energy.pot;
    }


    //! check energy based on list
    /* @param[in] _ptcl: particle data array
       @param[out] _energy: energy class
       @param[in] _ptcl_single_list: single particle list
       @param[in] _n_ptcl_single: single particle number
       @param[in] _ptcl_group_list: group particle list
       @param[in] _n_ptcl_group: group particle number
    */
    template<class Teng>
    void CalcEnergyHard(PtclHard* _ptcl, Teng & _energy,  const PS::S32* _ptcl_single_list, const PS::S32 _n_ptcl_single, const PS::S32* _ptcl_group_list, const PS::S32 _n_ptcl_group, const PS::S32 _i_single_start){
        _energy.kin = _energy.pot = _energy.tot = 0.0;
        for(PS::S32 i=_i_single_start; i<_n_ptcl_single+_i_single_start; i++){
            PtclHard* pi = &_ptcl[_ptcl_single_list[i]];
            _energy.kin += 0.5 * pi->mass * pi->vel * pi->vel;

            for(PS::S32 j=i+1; j<_n_ptcl_single; j++){
                PtclHard* pj = &_ptcl[_ptcl_single_list[j]];
                PS::F64vec rij = pi->pos - pj->pos;
                PS::F64 dr = sqrt(rij*rij + Int_pars_.eps2);
#ifdef INTEGRATED_CUTOFF_FUNCTION
                _energy.pot -= pj->mass*pi->mass/dr*(1.0 - CalcW(dr/Int_pars_.rout, Int_pars_.rin/Int_pars_.rout));  
#else
                if(dr<Int_pars_.rout) _energy.pot -= pj->mass*pi->mass*(1.0/dr*cutoff_pot(dr, Int_pars_.r_oi_inv, Int_pars_.r_A, Int_pars_.rin) - Int_pars_.pot_off);
#endif
            }

            for(PS::S32 j=0; j<_n_ptcl_group; j++){
                PtclHard* pj = &_ptcl[_ptcl_group_list[j]];
                PS::F64vec rij = pi->pos - pj->pos;
                PS::F64 dr = sqrt(rij*rij + Int_pars_.eps2);
#ifdef INTEGRATED_CUTOFF_FUNCTION
                _energy.pot -= pj->mass*pi->mass/dr*(1.0 - CalcW(dr/Int_pars_.rout, Int_pars_.rin/Int_pars_.rout));  
#else
                if(dr<Int_pars_.rout) _energy.pot -= pj->mass*pi->mass*(1.0/dr*cutoff_pot(dr, Int_pars_.r_oi_inv, Int_pars_.r_A, Int_pars_.rin) - Int_pars_.pot_off);
#endif
            }
        }

        for(PS::S32 i=0; i<_n_ptcl_group; i++){
            PtclHard* pi = &_ptcl[_ptcl_group_list[i]];
            _energy.kin += 0.5 * pi->mass * pi->vel * pi->vel;

            for(PS::S32 j=i+1; j<_n_ptcl_group; j++){
                PtclHard* pj = &_ptcl[_ptcl_group_list[j]];
                PS::F64vec rij = pi->pos - pj->pos;
                PS::F64 dr = sqrt(rij*rij + Int_pars_.eps2);
#ifdef INTEGRATED_CUTOFF_FUNCTION
                _energy.pot -= pj->mass*pi->mass/dr*(1.0 - CalcW(dr/Int_pars_.rout, Int_pars_.rin/Int_pars_.rout));  
#else
                if(dr<Int_pars_.rout) _energy.pot -= pj->mass*pi->mass*(1.0/dr*cutoff_pot(dr, Int_pars_.r_oi_inv, Int_pars_.r_A, Int_pars_.rin) - Int_pars_.pot_off);
#endif
            }
        }
        _energy.tot = _energy.kin + _energy.pot;
    }

    template<class Teng, class THint, class TARC>
    void CalcEnergyHardFull(PtclHard* _ptcl, const PS::S32 _n_ptcl, Teng& _E, Teng& _AE, Teng& _HE, Teng& _ESD, THint &_Hint, TARC& _Aint){
        _Hint.CalcEnergy(_HE);
        Teng TMP;
        _Aint.EnergyRecord(TMP,true);
        _Aint.EnergyRecord(_AE);
        CalcEnergyHard(_ptcl, _E, _n_ptcl);
        _ESD.tot = (_E.tot - _AE.kin - _AE.pot) + (TMP.kin+TMP.pot);
        _ESD.kin = (_E.kin - _AE.kin) + TMP.kin;
        _ESD.pot = (_E.pot - _AE.pot) + TMP.pot;
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
            ptcl_hard_[i].adr_org = adr;
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
#pragma omp parallel for schedule(dynamic)
        for(PS::S32 i=0; i<n; i++){
            PS::S32 adr = adr_array[i];
            ptcl_hard_[i].DataCopy(sys[adr]);
            ptcl_hard_[i].adr_org = adr;
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
#pragma omp parallel for schedule(dynamic)
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
//                                    const PS::ReallocatableArray<PS::S32> & adr_array,
                                    PS::ReallocatableArray<PS::S32> & removelist){
        const PS::S32 n = ptcl_hard_.size();
        //PS::ReallocatableArray<PS::S32> removelist(n);
        for(PS::S32 i=0; i<n; i++){
            //PS::S32 adr = adr_array[i];
            PS::S32 adr = ptcl_hard_[i].adr_org;
#ifdef HARD_DEBUG
            assert(sys[adr].id == ptcl_hard_[i].id);
#endif
            sys[adr].DataCopy(ptcl_hard_[i]);
            if(sys[adr].id<0&&sys[adr].status<0) removelist.push_back(adr);
        }
        //sys.removeParticle(removelist.getPointer(), removelist.size());
    }

    template<class Tsys>
    void writeBackPtclForOneClusterOMP(Tsys & sys){
        const PS::S32 n = ptcl_hard_.size();
#pragma omp parallel for schedule(dynamic)
        for(PS::S32 i=0; i<n; i++){
            PS::S32 adr = ptcl_hard_[i].adr_org;
            //PS::S32 adr = adr_array[i];
#ifdef HARD_DEBUG
            assert(sys[adr].id == ptcl_hard_[i].id);
#endif
            sys[adr].DataCopy(ptcl_hard_[i]);
        }
    }

    template<class Tsys>
    void writeBackPtclLocalOnlyOMP(Tsys & sys) {
        const PS::S32 n = ptcl_hard_.size();
#pragma omp parallel for schedule(dynamic)
        for(PS::S32 i=0; i<n; i++){
            PS::S32 adr = ptcl_hard_[i].adr_org;
            //PS::S32 adr = adr_array[i];
#ifdef HARD_DEBUG
            if(adr>=0) assert(sys[adr].id == ptcl_hard_[i].id);
#endif
            if(adr>=0) sys[adr].DataCopy(ptcl_hard_[i]);
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
        const PS::S32 n_cluster = _n_ptcl_in_cluster.size();
#ifdef HARD_DEBUG
        assert(n_cluster<ARRAY_ALLOW_LIMIT);
#endif        
        n_ptcl_in_cluster_.resizeNoInitialize(n_cluster);
        n_ptcl_in_cluster_disp_.resizeNoInitialize(n_cluster+1);
        n_ptcl_in_cluster_disp_[0] = 0;
        for(PS::S32 i=0; i<n_cluster; i++){
            n_ptcl_in_cluster_[i] = _n_ptcl_in_cluster[i];
#ifdef HARD_DEBUG
            assert(n_ptcl_in_cluster_[i]>1);
#endif
            n_ptcl_in_cluster_disp_[i+1] = n_ptcl_in_cluster_disp_[i] + n_ptcl_in_cluster_[i];
        }
        const PS::S32 n_ptcl = _adr_array.size();
#ifdef HARD_DEBUG
        assert(n_ptcl<ARRAY_ALLOW_LIMIT);
#endif        
        ptcl_hard_.resizeNoInitialize(n_ptcl);
        for(PS::S32 i=0; i<n_ptcl; i++){
            PS::S32 adr = _adr_array[i];
            ptcl_hard_[i].DataCopy(sys[adr]);
            ptcl_hard_[i].adr_org = adr;
            //  ptcl_hard_[i].n_ngb= sys[adr].n_ngb;
        }
    }

    void initailizeForIsolatedMultiCluster(const PS::S32 _n_ptcl,
                                           const PS::ReallocatableArray<PS::S32> & _n_ptcl_in_cluster){
#ifdef HARD_DEBUG
        assert(_n_ptcl<ARRAY_ALLOW_LIMIT);
#endif        
        ptcl_hard_.resizeNoInitialize(_n_ptcl);
        const PS::S32 n_cluster = _n_ptcl_in_cluster.size();
#ifdef HARD_DEBUG
        assert(n_cluster<ARRAY_ALLOW_LIMIT);
#endif        
        n_ptcl_in_cluster_.resizeNoInitialize(n_cluster);
        n_ptcl_in_cluster_disp_.resizeNoInitialize(n_cluster+1);
        n_ptcl_in_cluster_disp_[0] = 0;
        for(PS::S32 i=0; i<n_cluster; i++){
            n_ptcl_in_cluster_[i] = _n_ptcl_in_cluster[i];
#ifdef HARD_DEBUG
            assert(n_ptcl_in_cluster_[i]>1);
#endif
            n_ptcl_in_cluster_disp_[i+1] = n_ptcl_in_cluster_disp_[i] + n_ptcl_in_cluster_[i];
        }
    }

    template<class Tsys>
    void setPtclForIsolatedMultiClusterOMP(const Tsys & sys,
                                           const PS::ReallocatableArray<PS::S32> & _adr_array,
                                           const PS::ReallocatableArray<PS::S32> & _n_ptcl_in_cluster){
        const PS::S32 n_ptcl = _adr_array.size();
#pragma omp parallel for schedule(dynamic)
        for(PS::S32 i=0; i<n_ptcl; i++){
            PS::S32 adr = _adr_array[i];
            ptcl_hard_[i].DataCopy(sys[adr]);
            ptcl_hard_[i].adr_org = adr;
            //  ptcl_hard_[i].n_ngb = sys[adr].n_ngb;
        }
    }

    template<class Tsys>
    void writeBackPtclForMultiCluster(Tsys & sys, 
//                                      const PS::ReallocatableArray<PS::S32> & adr_array,
                                      PS::ReallocatableArray<PS::S32> & removelist){
        writeBackPtclForOneCluster(sys, removelist);
    }

    template<class Tsys>
    void writeBackPtclForMultiClusterOMP(Tsys & sys) { 
        writeBackPtclForOneClusterOMP(sys);
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
            const PS::S32 n_group = n_group_in_cluster_[i];
            const PS::S32 adr_ptcl_artifical = adr_first_ptcl_arti_in_cluster_[n_group_in_cluster_offset_[i]];

            driveForMultiClusterImpl(ptcl_hard_.getPointer(adr_head), n_ptcl, &sys[adr_ptcl_artifical], n_group, dt);
//#ifdef HARD_DEBUG
//            if(extra_ptcl.size()>0) fprintf(stderr,"New particle number = %d\n",extra_ptcl.size());
//#endif
//            for (PS::S32 j=0; j<extra_ptcl.size(); j++) {
//                PS::S32 adr = sys.getNumberOfParticleLocal();
//                PS::S32 rank = PS::Comm::getRank();
//                sys.addOneParticle(Tsptcl(extra_ptcl[j],rank,adr));
//            }
        }
    }

    template<class Tsys, class Tsptcl>
    void driveForMultiClusterOMP(const PS::F64 dt, Tsys & sys){
        const PS::S32 n_cluster = n_ptcl_in_cluster_.size();
        //const PS::S32 num_thread = PS::Comm::getNumberOfThread();
        //PS::ReallocatableArray<PtclHard> extra_ptcl[num_thread];
        //// For test
        //PS::ReallocatableArray<std::pair<PS::S32,PS::S32>> n_sort_list;
        //n_sort_list.resizeNoInitialize(n_cluster);
        //for(PS::S32 i=0; i<n_cluster; i++) {
        //    n_sort_list[i].first = n_ptcl_in_cluster_[i];
        //    n_sort_list[i].second= i;
        //}
        //std::sort(n_sort_list.getPointer(),n_sort_list.getPointer()+n_cluster,[](const std::pair<PS::S32,PS::S32> &a, const std::pair<PS::S32,PS::S32> &b){return a.first<b.first;});
#ifdef OMP_PROFILE        
        PS::ReallocatableArray<PS::F64> time_thread(num_thread);
        PS::ReallocatableArray<PS::S64> num_cluster(num_thread);
        for (PS::S32 i=0; i<num_thread; i++) {
          time_thread[i] = 0;
          num_cluster[i] = 0;
        }
#endif
#pragma omp parallel for schedule(dynamic)
        for(PS::S32 i=0; i<n_cluster; i++){
#ifdef OMP_PROFILE
            const PS::S32 ith = PS::Comm::getThreadNum();
            time_thread[ith] -= PS::GetWtime();
#endif
            //const PS::S32 i   = n_sort_list[k].second;
            const PS::S32 adr_head = n_ptcl_in_cluster_disp_[i];
            const PS::S32 n_ptcl = n_ptcl_in_cluster_[i];
            const PS::S32 n_group = n_group_in_cluster_[i];
            PS::S32 adr_ptcl_artifical;
            if(n_group>0) adr_ptcl_artifical=adr_first_ptcl_arti_in_cluster_[n_group_in_cluster_offset_[i]];
            else adr_ptcl_artifical=0;
#ifdef OMP_PROFILE
            num_cluster[ith] += n_ptcl;
#endif
#ifdef HARD_DEBUG_PROFILE
            PS::F64 tstart = PS::GetWtime();
#endif
            driveForMultiClusterImpl(ptcl_hard_.getPointer(adr_head), n_ptcl, &sys[adr_ptcl_artifical], n_group, dt);
#ifdef OMP_PROFILE
            time_thread[ith] += PS::GetWtime();
#endif
#ifdef HARD_DEBUG_PROFILE
            PS::F64 tend = PS::GetWtime();
            std::cerr<<"HT: "<<i<<" "<<ith<<" "<<n_cluster<<" "<<n_ptcl<<" "<<tend-tstart<<std::endl;
#endif
        }
//        if (n_cluster>0) {
//            PS::S32 rank = PS::Comm::getRank();
//            for(PS::S32 i=0; i<num_thread; i++) {
//#ifdef OMP_PROFILE        
//                std::cerr<<"thread: "<<i<<"  Hard Time="<<time_thread[i]<<"  n_ptcl="<<num_cluster[i]<<std::endl;
//#endif
//                for (PS::S32 j=0; j<extra_ptcl[i].size(); j++) {
//                    PS::S32 adr = sys.getNumberOfParticleLocal();
//                    sys.addOneParticle(Tsptcl(extra_ptcl[i][j],rank,adr));
//#ifdef HARD_DEBUG
//                    if(sys[adr].id==10477) {
//                        std::cerr<<"Add particle adr="<<adr;
//                        sys[adr].print(std::cerr);
//                        std::cerr<<std::endl;
//                        std::cerr<<" original: ";
//                        extra_ptcl[i][j].print(std::cerr);
//                        std::cerr<<std::endl;
//                    }
//                    if(extra_ptcl[i][j].id<0&&extra_ptcl[i][j].status<0) {
//                        std::cerr<<"Error: extra particle list contain ghost particle! i_thread="<<i<<" index="<<j<<" rank="<<rank<<" adr="<<adr<<std::endl;
//                        abort();
//                    }
//#endif
//                }
//            }
//        }
    }

    //! Find groups and create aritfical particles to sys
    /* @param[in,out] _sys: global particle system
       @param[in]     _dt_tree: tree time step for calculating r_search
     */
    template<class Tsys, class Tptcl>
    void findGroupsAndCreateArtificalParticlesOMP(Tsys & _sys, 
                                                  const PS::F64 _dt_tree) {
        // isolated clusters
        findGroupsAndCreateArtificalParticlesImpl<Tsys, Tptcl>(_sys, 
                                                               ptcl_hard_.getPointer(),
                                                               n_ptcl_in_cluster_,
                                                               n_ptcl_in_cluster_disp_,
                                                               n_group_in_cluster_,
                                                               n_group_in_cluster_offset_,
                                                               adr_first_ptcl_arti_in_cluster_,
                                                               r_bin_,
                                                               Int_pars_.rin,     
                                                               Int_pars_.rout,    
                                                               _dt_tree, 
                                                               id_offset_,
                                                               n_split_);

    }

    //! potential correction for single cluster
    /* The force kernel have self-interaction on the potential contribution, need to be excluded. _sys acc is updated
       @param[in,out] _sys: global particle system, acc is updated
       @param[in] _ptcl_list: list of single particle in _ptcl
       @param[in] _n_ptcl: number of single particles
     */
    template <class Tsys> 
    void correctPotWithCutoffOMP(Tsys& _sys, 
                                 const PS::ReallocatableArray<PS::S32>& _ptcl_list) {
        const PS::S32 n_ptcl = _ptcl_list.size();
#pragma omp parallel for 
        for (int i=0; i<n_ptcl; i++) {
            const PS::S32 k =_ptcl_list[i];
            _sys[k].pot_tot += _sys[k].mass / Int_pars_.rout;
#ifdef HARD_DEBUG
            // status may not be zero after binary disrupted
            // assert(_sys[k].status==0);
#endif
        }
    }

    //! Soft force correction due to different cut-off function
    /* Use tree neighbor search for local real particles including sending particles.
       Use cluster information correct artifical particles. 
       c.m. force is replaced by the averaged force on orbital particles
       Tidal tensor particle subtract the c.m. acc
       @param[in] _sys: global particle system, acc is updated
       @param[in] _tree: tree for force
       @param[in] _adr_send: particle in sending list of connected clusters
    */
    template <class Tsys, class Tpsoft, class Ttree, class Tepj>
    void correctForceWithCutoffTreeNeighborAndClusterOMP(Tsys& _sys,
                                                         Ttree& _tree,
                                                         const PS::ReallocatableArray<PS::S32>& _adr_send) {
        correctForceWithCutoffTreeNeighborAndClusterImp<Tsys, Tpsoft, Ttree, Tepj>(_sys, _tree, ptcl_hard_.getPointer(), n_ptcl_in_cluster_, n_ptcl_in_cluster_disp_, n_group_in_cluster_, n_group_in_cluster_offset_, adr_first_ptcl_arti_in_cluster_, _adr_send, Int_pars_.rin, Int_pars_.rout, n_split_, Int_pars_.eps2);
    }

    //! Soft force correction due to different cut-off function
    /* Use cluster member, first correct for artifical particles, then for cluster member
       c.m. force is replaced by the averaged force on orbital particles
       Tidal tensor particle subtract the c.m. acc
       @param[in] _sys: global particle system, acc is updated
    */
    template <class Tsys>
    void correctForceWithCutoffClusterOMP(Tsys& _sys) { 
        correctForceWithCutoffClusterImp(_sys, ptcl_hard_.getPointer(), n_ptcl_in_cluster_, n_ptcl_in_cluster_disp_, n_group_in_cluster_, n_group_in_cluster_offset_, adr_first_ptcl_arti_in_cluster_, Int_pars_.rin, Int_pars_.rout, n_split_, Int_pars_.eps2);
    }

    //! Soft force correction due to different cut-off function
    /* Use tree neighbor search for all particles.
       c.m. force is replaced by the averaged force on orbital particles
       Tidal tensor particle subtract the c.m. acc
       @param[in] _sys: global particle system, acc is updated
       @param[in] _tree: tree for force
       @param[in] _adr_ptcl_artifical_start: start address of artifical particle in _sys
    */
    template <class Tsys, class Tpsoft, class Ttree, class Tepj>
    void correctForceWithCutoffTreeNeighborOMP(Tsys& _sys,
                                               Ttree& _tree,
                                               const PS::S32 _adr_ptcl_artifical_start) {
        
        correctForceWithCutoffTreeNeighborImp<Tsys, Tpsoft, Ttree, Tepj>(_sys, _tree, _adr_ptcl_artifical_start, Int_pars_.rin, Int_pars_.rout, n_split_, Int_pars_.eps2);
    }

    ////! soft force correction for sending particles
    ///* Use tree neighbor search for sending particles
    //   
    // */
    //template <class Tsys, class Tpsoft, class Ttree, class Tepj>
    //void correctForceWithCutoffTreeNeighborSendOMP(Tsys& _sys,
    //                                               Ttree& _tree,
    //                                               const PS::ReallocatableArray<PS::S32> & _adr_send) {
    //    
    //    const PS::S32 n_send = _adr_send.size();
    //    // cutoff function parameter
    //    const PS::F64 r_oi_inv = 1.0/(_rout-_rin);
    //    const PS::F64 r_A = (_rout-_rin)/(_rout+_rin);
//#pragma omp parallel for schedule(dynamic)
    //    for (int i=0; i<n_send; i++) 
    //        const PS::S64 adr  = _adr_send(i);
    //        correctForceWithCutoffTreeNeighborOneParticleImp<Tpsoft, Ttree, Tepj>(_sys[adr], _tree, Int_pars_.rin, Int_pars_.rout, r_oi_inv, r_A, Int_pars_.eps2);
    //}

    //template<class Tsys, class Tsptcl>
    //void initialMultiClusterOMP(Tsys & sys, const PS::F64 dt_tree){
    //    const PS::S32 n_cluster = n_ptcl_in_cluster_.size();
    //    //	const PS::S32 ith = PS::Comm::getThreadNum();
//#pragma omp parallel for schedule(dynamic)
    //    for(PS::S32 i=0; i<n_cluster; i++){
    //        const PS::S32 adr_head = n_ptcl_in_cluster_disp_[i];
    //        const PS::S32 n_ptcl = n_ptcl_in_cluster_[i];
    //        SearchGroup<PtclHard> group;
    //        group.findGroups(ptcl_hard_.getPointer(adr_head), n_ptcl, n_split_);
    //        if (group.getPtclN()==2) group.searchAndMerge(ptcl_hard_.getPointer(adr_head), Int_pars_.rout);
    //        else group.searchAndMerge(ptcl_hard_.getPointer(adr_head), Int_pars_.rin);
    //        //group.searchAndMerge(ptcl_hard_.getPointer(adr_head), Int_pars_.rin);
    //        PS::ReallocatableArray<PtclHard> ptcl_new;
    //        group.generateList(ptcl_hard_.getPointer(adr_head), ptcl_new, r_bin_, Int_pars_.rin, Int_pars_.rout, dt_tree, id_offset_, n_split_);
//#pragma omp critical
    //        {
    //            for (PS::S32 j=0; j<ptcl_new.size(); j++) {
    //                PS::S32 adr = sys.getNumberOfParticleLocal();
    //                PS::S32 rank = PS::Comm::getRank();
    //                sys.addOneParticle(Tsptcl(ptcl_new[j],rank,adr));
    //            }
    //        }
    //        
    //    }        
    //}

#ifdef HARD_DEBUG_DUMP
    //! parameter dump function
    /* list: (S32)
       dt_limit_hard_: 1-2
       dt_min_hard_:   3-4
       eta_s_:         5-6
       time_origin_:   7-8
       r_bin_:         9-10
       sdfactor_:      11-12
       id_offset_:     13-14
       n_split_:       15
       
     */
    void pardump(FILE *_p_file){
        fwrite(&dt_limit_hard_, sizeof(PS::F32), 15, _p_file);
        Int_pars_.dump(_p_file);
        ARC_control_pert_.dump(_p_file);
        ARC_control_soft_.dump(_p_file);
    }
    
    //! Dumping data for debuging
    /* 
       @param[in] fname: file name to write
       @param[in] time_end: time ending
       @param[in] ptcl_bk: hard particle backup
       @param[in] n_ptcl: cluster member number
       @param[in] ptcl_arti_bk: artifical particle backup
       @param[in] n_arti: artifical particle number
       @param[in] n_group: number of groups
     */
    template<class Tpsoft>
    void dump(const char* fname, 
              const PS::F64 time_end, 
              PtclHard* ptcl_bk, 
              const PS::S32 n_ptcl, 
              const Tpsoft* ptcl_arti_bk,
              const PS::S32 n_arti,
              const PS::S32 n_group) {
        
        std::FILE* fp = std::fopen(fname,"w");
        if (fp==NULL) {
            std::cerr<<"Error: filename "<<fname<<" cannot be open!\n";
            abort();
        }
        fwrite(&time_end, sizeof(PS::F64),1,fp);
        PtclHardDump(fp, ptcl_bk, n_ptcl);
        fwrite(&n_arti, sizeof(PS::S32),1,fp);
        fwrite(&n_group, sizeof(PS::S32), 1, fp);
        for (int i=0; i<n_arti; i++) ptcl_arti_bk[i].writeBinary(fp);
        pardump(fp);
        fclose(fp);
    }
#endif
#ifdef HARD_DEBUG
    void parread(FILE *fp){
        size_t rcount = fread(&dt_limit_hard_, sizeof(PS::F32),15,fp);
        if (rcount<15) {
            std::cerr<<"Error: Data reading fails! requiring data number is 15, only obtain "<<rcount<<".\n";
            abort();
        }
        Int_pars_.read(fp);
        ARC_control_pert_.read(fp);
        ARC_control_soft_.read(fp);
#ifdef HARD_DEBUG_DEEP_CHECK
        ARC_control_pert_.setA(Newtonian_AW<Ptcl,ARC_pert_pars>,Newtonian_extA_test<Ptcl,PtclH4*,PtclForce*,ARC_pert_pars>,Newtonian_timescale<ARC_pert_pars>);
        ARC_control_soft_.setA(Newtonian_AW<Ptcl,ARC_pert_pars>,Newtonian_extA_test<Ptcl,PtclH4*,PtclForce*,ARC_pert_pars>,Newtonian_timescale<ARC_pert_pars>);
#else
        ARC_control_pert_.setA(Newtonian_AW<Ptcl,ARC_pert_pars>,Newtonian_extA_pert<Ptcl,PtclH4*,PtclForce*,ARC_pert_pars>,Newtonian_timescale<ARC_pert_pars>);
        ARC_control_soft_.setA(Newtonian_AW<Ptcl,ARC_pert_pars>,Newtonian_extA_soft<Ptcl,PtclH4*,PtclForce*,ARC_pert_pars>,Newtonian_timescale<ARC_pert_pars>);
#endif
#ifdef HARD_DEBUG_PRINT
        std::cerr<<"Parameters:"
                 <<"\nrout: "<<Int_pars_.rout
                 <<"\nrin: "<<Int_pars_.rin
                 <<"\neps2: "<<Int_pars_.eps2
                 <<"\npot_off: "<<Int_pars_.pot_off
                 <<"\ndt_limit_hard: "<<dt_limit_hard_
                 <<"\ndt_min_hard: "<<dt_min_hard_
                 <<"\neta_s: "<<eta_s_
                 <<"\nr_bin: "<<r_bin_
                 <<"\nsdfactor: "<<sdfactor_
                 <<"\nid_offset: "<<id_offset_
                 <<"\nn_split: "<<n_split_
                 <<std::endl;
#endif
    }
    
    template <class Tpsoft>
    void driveForMultiClusterOneDebug(PtclHard* _ptcl, const PS::S32 _n_ptcl, Tpsoft* _ptcl_artifical, const PS::S32 _n_group,  const PS::F64 _time_end) {
        driveForMultiClusterImpl(_ptcl, _n_ptcl, _ptcl_artifical, _n_group, _time_end);
    }

    void set_slowdown_factor(const PS::F64 _slowdown_factor) {
        sdfactor_ = _slowdown_factor;
    }

#endif

};
