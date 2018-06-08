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

//! Soft force correction due to different cut-off function
/* Use cluster member, first correct for artifical particles, then for cluster member, cluster member acc is replaced by c.m. acc
   @param[in] _sys: global particle system, acc is updated
   @param[in] _ptcl_local: particle in systme_hard
   @param[in] _n_arti_start: start index of artifical particle in _sys.
   @param[in] _n_arti_end: artifical particle total number (end index+1)
   @param[in] _n_split: artifical particle splitting number
   @param[in] _n_ptcl_in_cluster: number of particles in clusters
   @param[in] _n_ptcl_in_cluster_offset: boundary of clusters in _adr_sys_in_cluster
   @param[in] _rin: cutoff inner radius;
   @param[in] _rout: cutoff outer radius;
   @param[in] _eps_sq: softing eps square
 */
template <class Tsys>
void CorrectForceWithCutoffClusterOMP(Tsys _sys, 
                                      PtclHard* _ptcl_local,
                                      const PS::S32 _n_arti_start,
                                      const PS::S32 _n_arti_end,
                                      const PS::S32 _n_split,
                                      const PS::ReallocatableArray<PS::S32>& _n_ptcl_in_cluster,
                                      const PS::ReallocatableArray<PS::S32>& _n_ptcl_in_cluster_offset,
                                      const PS::F64 _rin;
                                      const PS::F64 _rout;
                                      const PS::F64 _eps_sq=0.0) {
    // cutoff function parameter
    const PS::F64 r_oi_inv = 1.0/(_rout-_rin);
    const PS::F64 r_A = (_rout-_rin)/(_rout+_rin);
    // first obtain the correction for aritifical particles
    const PS::S32 ptcl_art_offset = 2*_n_split+1;
#pragma omp for schedule(dynamic)
    for (int i=_n_arti_start, i<_n_arti_end, i+=ptcl_art_offset){
        PS::S32 i_cluster=_sys[i+2].status;
        PS::S32 j_cm= i+2*_n_split;
        // loop all artifical particles (tidal tensor particle)
        for (int j=i; j<i+_n_split; j++) {
            PS::F64vec& pos_j= _sys[j].pos;
            PS::F64vec& acc_j= _sys[j].acc;
            PS::F64& pot_j = _sys[j].pot_tot;
            // loop orbital artifical particle
            for (int k=i+_n_split; k<j_cm; k++) {
                auto& ptcl_k = _sys[k];
                CalcAccPotShortWithLinearCutoff(pos_j, acc_j, pot_j, 
                                                ptcl_k.pos, ptcl_k.mass, ptcl_k.mass_bk, 
                                                2, _eps_sq,
                                                r_oi_inv, r_A, _rout, _rin);
            }
            // loop real particle
            for (int k=_n_ptcl_in_cluster_offset[i_cluster]; k<_n_ptcl_in_cluster_offset[i_cluster+1]; k++) {
                PtclHard* ptcl_k_ptr = &_ptcl_local[k];
                PS::S32 pot_control_flag = ptcl_k_ptr->status>0? 1: 0;
                CalcAccPotShortWithLinearCutoff(pos_j, acc_j, pot_j, 
                                                ptcl_k_ptr->pos, ptcl_k_ptr->mass, ptcl_k_ptr->mass_bk, 
                                                pot_control_flag, _eps_sq,
                                                r_oi_inv, r_A, _rout, _rin);
            }
        }
        // for c.m. particle
        PS::F64vec& pos_j= _sys[j_cm].pos;
        PS::F64vec& acc_j= _sys[j_cm].acc;
        PS::F64& pot_j = _sys[j_cm].pot_tot;
        // loop artifical particle orbital
        for (int k=i+_n_split; k<j_cm; k++) {
            auto& ptcl_k = _sys[k];
            CalcAccPotShortWithLinearCutoff(pos_j, acc_j, pot_j, 
                                            ptcl_k.pos, ptcl_k.mass, ptcl_k.mass_bk, 
                                            2, _eps_sq,
                                            r_oi_inv, r_A, _rout, _rin);
        }
        // loop real particle
        for (int k=_n_ptcl_in_cluster_offset[i_cluster]; k<_n_ptcl_in_cluster_offset[i_cluster+1]; k++) {
            PtclHard* ptcl_k_ptr = &_ptcl_local[k];
            PS::S32 pot_control_flag = ptcl_k_ptr->status>0? 1: 0;
            CalcAccPotShortWithLinearCutoff(pos_j, acc_j, pot_j, 
                                            ptcl_k_ptr->pos, ptcl_k_ptr->mass, ptcl_k_ptr->mass_bk, 
                                            pot_control_flag, _eps_sq,
                                            r_oi_inv, r_A, _rout, _rin);
        }
    }

    // then obtain correction for real particles in clusters
    const PS::S32 n_cluster = _n_ptcl_in_cluster.size();
#pragma omp for schedule(dynamic)
    for (int i=0; i<n_cluster; i++) {
        PS::S32 adr_start= _n_ptcl_in_cluster_offset[i];
        PS::S32 adr_end= _n_ptcl_in_cluster_offset[i+1];
        for (int j=adr_start; j<adr_end; j++) {
            PS::S64 adr = _ptcl_local[j].adr_org;
#ifdef HARD_DEBUG
            assert(_sys[adr].id==_ptcl_local[j].id);
#endif
            PS::F64vec& pos_j= _sys[adr].pos;
            PS::F64vec& acc_j= _sys[adr].acc;
            PS::F64& pot_j = _sys[adr].pot_tot;
            
            PS::S64 stat_j = _sys[adr].status;
            //self-potential correction for non-group member 
            if(stat_j==0) pot_j += sys[adr].mass/_rout;

            // cluster member
            for (int k=adr_start; k=adr_end; k++) {
                if(k==j) continue;
                PtclHard* ptcl_k_ptr = &_ptcl_local[k];
                PS::S32 pot_control_flag = ptcl_k_ptr->status>0? 1: 0;
                CalcAccPotShortWithLinearCutoff(pos_j, acc_j, pot_j, 
                                                ptcl_k_ptr->pos, ptcl_k_ptr->mass, ptcl_k_ptr->mass_bk, 
                                                pot_control_flag, _eps_sq,
                                                r_oi_inv, r_A, _rout, _rin);
            }

            // orbital artifical particle
            
#ifdef HARD_DEBUG
            assert(-_sys[stat_j].id==_ptcl_local[j].id);
#endif
            // group member, use c.m. acc
            if(stat_j>0) acc_j = _sys[stat_j].acc;
        }
    }
}

//! Soft force correction due to different cut-off function
/* Use tree neighbor search
   @param[in] _sys: global particle system, acc is updated
   @param[in] _tree: tree for force
   @param[in] _ptcl_local: particle in systme_hard
   @param[in] _n_arti_start: start index of artifical particle in _sys.
   @param[in] _n_arti_end: artifical particle total number (end index+1)
   @param[in] _n_split: artifical particle splitting number
   @param[in] _n_ptcl_in_cluster: number of particles in clusters
   @param[in] _n_ptcl_in_cluster_offset: boundary of clusters in _adr_sys_in_cluster
   @param[in] _rin: cutoff inner radius;
   @param[in] _rout: cutoff outer radius;
   @param[in] _eps_sq: softing eps square
 */

template <class Tsys, class Ttree, class Tepj>
void CorrectForceWithCutoffTreeNeighborOMP(Tsys _sys,
                                           Ttree _tree,
                                           PtclHard* _ptcl_local,
                                           const PS::S32 _n_arti_start,
                                           const PS::S32 _n_arti_end,
                                           const PS::S32 _n_split,
                                           const PS::ReallocatableArray<PS::S32>& _n_ptcl_in_cluster,
                                           const PS::ReallocatableArray<PS::S32>& _n_ptcl_in_cluster_offset,
                                           const PS::F64 _rin;
                                           const PS::F64 _rout;
                                           const PS::F64 _eps_sq=0.0) {

    // cutoff function parameter
    const PS::F64 r_oi_inv = 1.0/(_rout-_rin);
    const PS::F64 r_A = (_rout-_rin)/(_rout+_rin);
    // first obtain the correction for aritifical particles
    const PS::S32 ptcl_art_offset = 2*_n_split+1;
#pragma omp for schedule(dynamic)
    for (int i=_n_arti_start, i<_n_arti_end, i+=ptcl_art_offset){
        PS::S32 i_cluster=_sys[i+2].status;
        // loop all artifical particles (tidal tensor particle)
        for (int j=i; j<i+_n_split; j++) {
            PS::F64vec& pos_j= _sys[j].pos;
            PS::F64vec& acc_j= _sys[j].acc;
            PS::F64& pot_j = _sys[j].pot_tot;
            // loop artifical particle (orbital and c.m.)
            for (int k=i+_n_split; k<i+ptcl_art_offset; k++) {
                auto& ptcl_k = _sys[k];
                CalcAccPotShortWithLinearCutoff(pos_j, acc_j, pot_j, 
                                                ptcl_k.pos, ptcl_k.mass, ptcl_k.mass_bk, 
                                                2, _eps_sq,
                                                r_oi_inv, r_A, _rout, _rin);
            }
            // loop real particle
            for (int k=_n_ptcl_in_cluster_offset[i_cluster]; k<_n_ptcl_in_cluster_offset[i_cluster+1]; k++) {
                PtclHard* ptcl_k_ptr = &_ptcl_local[k];
                PS::S32 pot_control_flag = ptcl_k_ptr->status>0? 1: 0;
                CalcAccPotShortWithLinearCutoff(pos_j, acc_j, pot_j, 
                                                ptcl_k_ptr->pos, ptcl_k_ptr->mass, ptcl_k_ptr->mass_bk, 
                                                pot_control_flag, _eps_sq,
                                                r_oi_inv, r_A, _rout, _rin);
            }
        }
        // for c.m. particle
        PS::S32 j_cm= i+2*_n_split;
        PS::F64vec& pos_j= _sys[j_cm].pos;
        PS::F64vec& acc_j= _sys[j_cm].acc;
        PS::F64& pot_j = _sys[j_cm].pot_tot;
        // loop artifical particle (orbital and c.m.)
        for (int k=i+_n_split; k<i+ptcl_art_offset; k++) {
            auto& ptcl_k = _sys[k];
            CalcAccPotShortWithLinearCutoff(pos_j, acc_j, pot_j, 
                                            ptcl_k.pos, ptcl_k.mass, ptcl_k.mass_bk, 
                                            2, _eps_sq,
                                            r_oi_inv, r_A, _rout, _rin);
        }
        // loop real particle
        for (int k=_n_ptcl_in_cluster_offset[i_cluster]; k<_n_ptcl_in_cluster_offset[i_cluster+1]; k++) {
            PtclHard* ptcl_k_ptr = &_ptcl_local[k];
            PS::S32 pot_control_flag = ptcl_k_ptr->status>0? 1: 0;
            CalcAccPotShortWithLinearCutoff(pos_j, acc_j, pot_j, 
                                            ptcl_k_ptr->pos, ptcl_k_ptr->mass, ptcl_k_ptr->mass_bk, 
                                            pot_control_flag, _eps_sq,
                                            r_oi_inv, r_A, _rout, _rin);
        }
    }

    // then obtain correction for real particles in clusters
    const PS::S32 n_cluster = _n_ptcl_in_cluster.size();
#pragma omp for schedule(dynamic)
    for (int i=0; i<n_cluster; i++) {
        PS::S32 adr_start= _n_ptcl_in_cluster_offset[i];
        PS::S32 adr_end= _n_ptcl_in_cluster_offset[i+1];
        for (int j=adr_start; j<adr_end; j++) {
            PS::S64 adr = _ptcl_local[j].adr_org;
            // only do for local particles
            if(adr>=0) {
                Tepj * nbl = NULL;
                PS::S32 n_ngb = _tree.getNeighborListOneParticle(sys[i], nbl) - 1;
#ifdef HARD_DEBUG
                assert(_sys[i].n_ngb >= 0);
#endif
                PS::S64 stat_j = _sys[adr].status;
                // self-potential correction 
                if (stat_j==0) sys[i].pot_tot += sys[i].mass/r_out;
                else if (stat_j<0) sys[i].pot_tot += sys[i].mass_bk/r_out;

                // loop neighbors
                for(PS::S32 k=0; k<n_ngb+1; k++){
                    
                }
            }
        }
    }
}

class PtclHard: public Ptcl{
public:
    PS::S32 id_cluster;
    PS::S32 adr_org;

    PtclHard() {}

    PtclHard(const Ptcl &p): Ptcl(p) {}
    PtclHard(const ParticleBase &p): Ptcl(p) {}

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
    void dump(FILE *fp) {
        fwrite(this, sizeof(*this),1,fp);
    }

    void read(FILE *fp) {
        size_t rcount = fread(this, sizeof(*this),1,fp);
        if (rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            abort();
        }
    }
    void print(std::ostream & fout){
        Ptcl::print(fout);
        std::cerr<<" id_cluster="<<id_cluster
                 <<" adr_org="<<adr_org;
    }
};


void PtclHardDump(FILE *fp, PtclHard * ptcl, const PS::S32 n) {
    fwrite(&n, sizeof(PS::S32), 1, fp);
    for(int i=0; i<n; i++) {
        ptcl[i].dump(fp);
    }
    PS::F64 ptcl_st_dat[3];
    ptcl_st_dat[0] = Ptcl::search_factor;
    ptcl_st_dat[1] = Ptcl::r_search_min;
    ptcl_st_dat[2] = Ptcl::mean_mass_inv;
    fwrite(ptcl_st_dat, sizeof(PS::F64),3,fp);
}

void PtclHardRead(FILE *fp, PS::ReallocatableArray<PtclHard> & ptcl) {
    PS::S32 n;
    size_t rcount = fread(&n, sizeof(PS::S32),1,fp);
    if (rcount<1) {
        std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
        abort();
    }
    PtclHard ptmp;
    for(int i=0; i<n; i++) {
        ptmp.read(fp);
        ptcl.push_back(ptmp);
    }
    PS::F64 ptcl_st_dat[3];
    rcount = fread(ptcl_st_dat, sizeof(PS::F64),3,fp);
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
    PS::ReallocatableArray<PtclHard> ptcl_hard_;
#ifdef PROFILE
    PS::S64 ARC_substep_sum;
    PS::F64 ARC_n_groups;
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
    
    PS::ReallocatableArray<PS::S32> n_ptcl_in_cluster_;
    PS::ReallocatableArray<PS::S32> n_ptcl_in_cluster_disp_;
    PS::ReallocatableArray<PS::S32> n_group_in_cluster_;
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
                                                   PS::ReallocatableArray<PS::S32> &_n_ptcl_in_cluster,
                                                   PS::ReallocatableArray<PS::S32> &_n_ptcl_in_cluster_disp,
                                                   PS::ReallocatableArray<PS::S32> &_n_group_in_cluster,
                                                   const PS::F64 _rbin,
                                                   const PS::F64 _rin,
                                                   const PS::F64 _rout,
                                                   const PS::F64 _dt_tree,
                                                   const PS::S64 _id_offset,
                                                   const PS::S32 _n_split) { 
        const PS::S32 n_cluster = _n_ptcl_in_cluster.size();
        _n_group_in_cluster.resizeNoInitialize(n_cluster);
        n_group_member_remote_=0;

        const PS::S32 num_thread = PS::Comm::getNumberOfThread();
        PS::ReallocatableArray<PtclHard> ptcl_artifical[num_thread];

#pragma omp for schedule(dynamic)
        for (PS::S32 i=0; i<n_cluster; i++){
            const PS::S32 ith = PS::Comm::getThreadNum();
            PtclHard* ptcl_in_cluster = _ptcl_local + _n_ptcl_in_cluster_disp[i];
            const PS::S32 n_ptcl = _n_ptcl_in_cluster[i];
            // search groups
            SearchGroup<PtclHard> group;
            // merge groups
            if (n_ptcl==2) group.searchAndMerge(ptcl_in_cluster, n_ptcl, _rout);
            else group.searchAndMerge(ptcl_in_cluster, n_ptcl, _rin);

            // generate artifical particles,
            group.generateList(i, ptcl_in_cluster, n_ptcl, ptcl_artifical[ith], _n_group_in_cluster[i], _rbin, _rin, _rout, _dt_tree, _id_offset, _n_split);
        }

        // add artifical particle to particle system
        PS::S32 rank = PS::Comm::getRank();
        const PS::S32 n_artifical_per_group = 2*_n_split+1;
        for(PS::S32 i=0; i<num_thread; i++) {
            // ptcl_artifical should be integer times of 2*n_split+1
            assert(ptcl_artifical[i].size()%n_artifical_per_group==0);
            // Add particle to ptcl sys
            for (PS::S32 j=0; j<ptcl_artifical[i].size(); j++) {
                PS::S32 adr = _sys.getNumberOfParticleLocal();
                ptcl_artifical[i][j].adr_org=adr;
                _sys.addOneParticle(Tptcl(ptcl_artifical[i][j],rank,adr));
            }
            PS::S32 group_offset=0, j_group_recored=-1;
            // Update the status of group members to c.m. address in ptcl sys. Notice c.m. is at the end of an artificial particle group
            for (PS::S32 j=0; j<ptcl_artifical[i].size(); j+=n_artifical_per_group) {
                // obtain group member nember
                PS::S32 j_cm = j+2*_n_split;
                PS::S32 n_members = ptcl_artifical[i][j_cm].status;
                PS::S32 i_cluster = ptcl_artifical[i][j+2].status;
                PS::S32 j_group = ptcl_artifical[i][j+3].status;
                PS::F64 rsearch_member=ptcl_artifical[i][j+2].r_search;
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
                        _sys[ptcl_k].status = ptcl_artifical[i][j_cm].adr_org;
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
                    _ptcl_local[kl].status = ptcl_artifical[i][j_cm].adr_org;
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
            }
        }
    }

    
    void driveForMultiClusterImpl(PtclHard * ptcl_org,
                                  const PS::S32 n_ptcl,
                                  const PS::F64 time_end,
                                  PS::ReallocatableArray<PtclHard> & ptcl_new,
                                  const bool first_int_flag=false) {
#ifdef HARD_CHECK_ENERGY
        std::map<PS::S32, PS::S32> N_count;  // counting number of particles in one cluster
        HardEnergy E0, E1;
        HardEnergy AE0, AE1;
        HardEnergy HE0, HE1;
        HardEnergy ESD0, ESD1;
#endif
#ifdef HARD_DEBUG_PROFILE
        N_count[n_ptcl]++;
#endif
#ifdef HARD_DEBUG_DUMP
        PS::ReallocatableArray<PtclHard> ptcl_bk;
        ptcl_bk.reserve(n_ptcl);
        for(int i=0; i<n_ptcl; i++) ptcl_bk.pushBackNoCheck(ptcl_org[i]);
#endif
        PS::S32 nstepcount = 0;

//#ifdef HERMITE
//        if(n_ptcl>5) {
        SearchGroup<PtclHard> group;
        group.findGroups(ptcl_org, n_ptcl, n_split_);

#ifdef HARD_DEBUG
        // If only ghost particle exist, this assertion happen
        assert(group.getPtclN()>0);
#endif

#ifdef TIDAL_TENSOR
        for (PS::S32 i=0; i<group.getNumOfGroups(); i++) 
            subtractFcmAndRecoverCMVec(ptcl_org, group.getPtclIndex(i), group.getGroup(i), group.getGroupN(i), group.getGroupPertList(i,n_split_));
#endif

#ifdef HARD_CM_KICK
        PS::F64 dt_soft = time_end;
        if(first_int_flag) dt_soft *=0.5;
        softKickForCM(ptcl_org, group.getPtclList(), group.getNumOfGroups(),group.getGroupPertList(0,n_split_), dt_soft, n_split_);
#endif

        if(group.getPtclN()==1) {
#ifdef HARD_DEBUG
            assert(group.getNumOfGroups()==1);
#endif
            PtclHard* pcm = &ptcl_org[group.getPtclIndex(0)];
            PS::S32 iact = 0;
            
            ARCIntegrator<PtclHard, PtclH4, PtclForce> Aint(ARC_control_soft_, Int_pars_);
            Aint.reserveARMem(1);
            // Aint.reservePertMem(1);
            group.getBinPars(Aint.bininfo[0],ptcl_org,0,n_split_);
#ifdef HARD_DEBUG_PRINT            
            if(Aint.bininfo[0].ax>Int_pars_.rout)  Aint.bininfo[0].print(std::cerr,13);
#endif
            const PS::S32 *group_list = group.getGroup(0);
            const PS::S32 group_n = group.getGroupN(0);
            Aint.addOneGroup(ptcl_org, group_list, group_n, group.getGroupPertList(0,n_split_), n_split_);
            Aint.updateCM(pcm, &iact, 1);

            Aint.initialSlowDown(time_end, sdfactor_);
            Aint.initial();

#ifdef ARC_SYM_SD_PERIOD
            PS::S32 kp=0;
            Aint.adjustSlowDownPeriod(time_end, &kp);
#else
            Aint.adjustSlowDown(time_end);
#endif

#ifdef HARD_CHECK_ENERGY
            Aint.EnergyRecord(AE0);
#endif 

#ifdef ARC_SYM
#ifdef ARC_SYM_SD_PERIOD
            if(group_n==2&&kp>0) nstepcount +=Aint.integrateOneStepSymTwo(0, time_end, kp);
            else nstepcount +=Aint.integrateOneStepSym(0, time_end, dt_limit_hard_);
#else
            nstepcount +=Aint.integrateOneStepSym(0, time_end, dt_limit_hard_);
#endif
#else 
            nstepcount +=Aint.integrateOneStepExt(0, time_end, dt_limit_hard_);
#endif
            
            pcm->pos += pcm->vel * time_end;

            Aint.updateCM(pcm, &iact, 1);
            Aint.resolve();
#ifdef HARD_CHECK_ENERGY
            Aint.EnergyRecord(AE1);
#ifdef HARD_DEBUG_PRINT
            fprintf(stderr,"Slowdown factor = %e\n", Aint.getSlowDown(0));
            fprintf(stderr,"ARC Energy: init =%e, end =%e, diff =%e, error = %e\n", 
                    AE0.kin+AE0.pot, AE1.kin+AE1.pot, AE1.kin+AE1.pot-AE0.kin-AE0.pot, (AE1.kin+AE1.pot+AE1.tot-AE0.kin-AE0.pot-AE0.tot)/AE0.tot);
#endif
#endif
#ifdef ARC_DEBUG_PRINT
#ifdef ARC_SYM_SD_PERIOD
            Aint.info_print(std::cerr, ARC_n_groups, 1, 1, group.getGroupN(0),dt_limit_hard_,kp);
#else
            Aint.info_print(std::cerr, ARC_n_groups, 1, 1, group.getGroupN(0),dt_limit_hard_,0);
#endif
#endif
#ifdef PROFILE
            //ARC_substep_sum += Aint.getNsubstep();
            ARC_substep_sum += nstepcount;
            ARC_n_groups += 1;
#endif
        }
        else {
            
            HermiteIntegrator<PtclHard> Hint;
            Hint.setParams(eta_s_, Int_pars_.rin, Int_pars_.rout, Int_pars_.eps2);
            Hint.setPtcl(ptcl_org,n_ptcl,group.getPtclList(),group.getPtclN());

            PS::F64 time_sys=0.0, time_now;
#ifdef FIX_STEP_DEBUG
            PS::F64 dt_limit = dt_limit_hard_;
#else
            PS::F64 dt_limit = calcDtLimit(time_sys, dt_limit_hard_, dt_min_hard_);
#endif
            
            PS::S32 group_act_n = 0;
            PS::ReallocatableArray<PS::S32> group_act_list; //active group_list act adr
            // ReallocatableArray<PS::S32> group_list;     //group.adr list
            // ReallocatableArray<PS::S32> status;      //ptcl -> group.adr [non cm is -1] (value of Ptcl.status)
            // ReallocatableArray<PS::S32> status_map;  //ptcl -> group_list index [non cm is -1]
            // ReallocatableArray<PS::S32> adr_cm;         //group_list index -> ptcl.cm
            // group.findGroups(group_list, status, status_map,  adr_cm, group_act_n, ptcl_org, n_ptcl);

            group_act_list.resizeNoInitialize(group.getPtclN());
            
            // Initial Aint
            PS::S32 n_groups = group.getNumOfGroups();
            ARCIntegrator<PtclHard, PtclH4, PtclForce> Aint(ARC_control_pert_, Int_pars_);
            Aint.reserveARMem(n_groups);
            PS::F64 dr_search[n_groups];
            for (int i=0; i<n_groups; i++) {
                group.getBinPars(Aint.bininfo[i],ptcl_org,i,n_split_);
                auto &bini= Aint.bininfo[i];
                /* Notice when artificial particles exist, the c.m. particle may find neighbor among them from another binary.
                   This can result in merging of two binarties to one group.
                   To avoid no perturbers, the rsearch should add the maximum distance of artifical particles to the c.m., which is apo-center distance
                 */
                dr_search[i] = bini.ax*(bini.ecc+1.0); 
            }            
            Hint.searchPerturber(dr_search,n_groups);

            // first particles in Hint.Ptcl are c.m.
            Aint.reservePertMem(Hint.getPertListSize());
            for (int i=0; i<n_groups; i++) {
#ifdef HARD_DEBUG
                assert(Hint.getPertN(i)>0);
#endif
                Aint.addOneGroup(ptcl_org, group.getGroup(i), group.getGroupN(i), group.getGroupPertList(i,n_split_), n_split_, Hint.getPtcl(), Hint.getForce(), Hint.getPertList(i), Hint.getPertN(i)); 
                
            }
            Aint.initialSlowDown(dt_limit, sdfactor_);
            Aint.initial();

#ifdef HARD_CHECK_ENERGY
            CalcEnergyHardFull(ptcl_org, E0, AE0, HE0, ESD0, Hint, Aint, group);
#endif

            bool fail_flag=Hint.initialize(dt_limit, dt_min_hard_, group_act_list.getPointer(), group_act_n, n_groups, &Aint);

            if(fail_flag) {
#ifdef HARD_DEBUG_DUMP
                std::cerr<<"Dump hard data. tend="<<time_end<<" n_ptcl="<<n_ptcl<<"\n";
                dump("hard_dump",time_end, (PS::S32)first_int_flag, ptcl_bk.getPointer(), n_ptcl);
                abort();
#endif
            }

#ifdef HARD_CHECK_ENERGY
            PS::ReallocatableArray<PS::F64> slowdownrecord;
            slowdownrecord.resizeNoInitialize(n_groups);
#endif

            while(time_sys<time_end) {
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
                //Aint.updateSlowDown(time_sys);
#ifdef HARD_CHECK_ENERGY
                for(int k=0; k<n_groups; k++) {
                    slowdownrecord[k] = std::max(slowdownrecord[k], Aint.getSlowDown(k));
                    assert(Aint.getSlowDown(k)>=1.0);
                }
#endif
                //Aint.integrateOneStepList(group_act_list.getPointer(), group_act_n, time_sys, dt_limit);
                nstepcount +=Aint.integrateOneStepList(time_sys, std::min(dt_limit,dt_h));
                fail_flag = Hint.integrateOneStep(time_sys,dt_limit,dt_min_hard_,true,&Aint);
                
                if(fail_flag) {
#ifdef HARD_DEBUG_DUMP
                    std::cerr<<"Dump hard data. tend="<<time_end<<" n_ptcl="<<n_ptcl<<"\n";
                    dump("hard_dump",time_end, (PS::S32)first_int_flag, ptcl_bk.getPointer(), n_ptcl);
                    abort();
#endif
                }
                //Hint.SortAndSelectIp(group_act_list.getPointer(), group_act_n, n_groups);
                Hint.SortAndSelectIp();
            }
        
            Hint.moveCM(time_end);
            Hint.shiftBackCM();
            Aint.updateCM(Hint.getPtcl());
            Aint.resolve();
            Hint.writeBackPtcl(ptcl_org,n_ptcl,group.getPtclList(),group.getPtclN());

#ifdef ARC_DEBUG_PRINT
            Aint.info_print(std::cerr, ARC_n_groups, n_groups, group.getPtclN(), n_ptcl, dt_limit_hard_,0);
#endif
#ifdef HARD_CHECK_ENERGY
            CalcEnergyHardFull(ptcl_org, E1, AE1, HE1, ESD1, Hint, Aint, group);
#ifdef HARD_DEBUG_PRINT
            fprintf(stderr,"Slowdown factor = ");
            for(int k=0; k<n_groups; k++) 
                fprintf(stderr,"%e; ",slowdownrecord[k]);
            fprintf(stderr,"\n");
            fprintf(stderr,"H4  Energy: init =%e, end =%e, diff =%e, kini =%e kinf =%e poti =%e potf =%e\nARC Energy: init =%e, end =%e, diff =%e, error = %e\nTot Energy: init =%e, end =%e, diff =%e, kin =%e pot =%e, Tot-H4-ARC =%e\nTSD Energy: init =%e, end =%e, diff =%e, kin =%e pot =%e\n", 
                    HE0.tot, HE1.tot, HE1.tot-HE0.tot, HE0.kin, HE1.kin, HE0.pot, HE1.pot, 
                    AE0.kin+AE0.pot, AE1.kin+AE1.pot, AE1.kin+AE1.pot-AE0.kin-AE0.pot, (AE1.kin+AE1.pot+AE1.tot-AE0.kin-AE0.pot-AE0.tot)/AE0.tot,
                    E0.tot, E1.tot, E1.tot-E0.tot, E1.kin, E1.pot, E1.tot-HE1.tot-AE1.kin-AE1.pot,
                    ESD0.tot, ESD1.tot, ESD1.tot-ESD0.tot, ESD1.kin, ESD1.pot);
            Hint.printStepHist();
#endif
#ifdef HARD_DEBUG_DUMP
            PS::F64 dEtot = E1.tot-E0.tot;
            if(fabs(dEtot)>1e-4) {
                std::cerr<<"Hard energy significant: "<<dEtot<<std::endl;
                std::cerr<<"Dump data:"<<std::endl;
                dump("hard_dump",time_end, (PS::S32)first_int_flag, ptcl_bk.getPointer(), n_ptcl);
                abort();
            }
#endif
#endif
#ifdef PROFILE
            //ARC_substep_sum += Aint.getNsubstep();
            ARC_substep_sum += nstepcount;
            ARC_n_groups += n_groups;
#endif
        }
            
        //group.resolveGroups(ptcl_org, n_ptcl, group_ptcl_glb.getPointer(), group_list.size(), group_list.getPointer(), adr_cm.getPointer());
        //group.resolveGroups();
        //updateRSearch(ptcl_org, group.getPtclList(), group.getPtclN(), time_end);

        //if (group.getPtclN()==2) group.searchAndMerge(ptcl_org, Int_pars_.rout);
        //else group.searchAndMerge(ptcl_org, Int_pars_.rin);
        //group.searchAndMerge(ptcl_org, Int_pars_.rout);
        // Kickcorrect(ptcl_org, group.getRoutChangeList());
        //group.generateList(ptcl_org, ptcl_new, r_bin_,Int_pars_.rin, Int_pars_.rout, time_end, id_offset_, n_split_);

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
#ifdef HARD_DEBUG_PROFILE
        for(PS::S32 i=0;i<20;i++) N_count[i]=0;
#endif
#ifdef PROFILE
        ARC_substep_sum = 0;
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
        ARC_control_pert_.setA(Newtonian_AW<PtclHard,ARC_pert_pars>,Newtonian_extA_test<PtclHard,PtclH4*,PtclForce*,ARC_pert_pars>,Newtonian_timescale<ARC_pert_pars>);
        ARC_control_soft_.setA(Newtonian_AW<PtclHard,ARC_pert_pars>,Newtonian_extA_test<PtclHard,PtclH4*,PtclForce*,ARC_pert_pars>,Newtonian_timescale<ARC_pert_pars>);
#else
        ARC_control_pert_.setA(Newtonian_AW<PtclHard,ARC_pert_pars>,Newtonian_extA_pert<PtclHard,PtclH4*,PtclForce*,ARC_pert_pars>,Newtonian_timescale<ARC_pert_pars>);
        ARC_control_soft_.setA(Newtonian_AW<PtclHard,ARC_pert_pars>,Newtonian_extA_soft<PtclHard,PtclH4*,PtclForce*,ARC_pert_pars>,Newtonian_timescale<ARC_pert_pars>);
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
        n_ptcl_in_cluster_disp_.resizeNoInitialize(n_cluster+1);
        n_ptcl_in_cluster_disp_[0] = 0;
        for(PS::S32 i=0; i<n_cluster; i++){
            n_ptcl_in_cluster_disp_[i+1] = n_ptcl_in_cluster_disp_[i] + n_ptcl_in_cluster_[i];
        }
    }


    // for NON-ISOLATED CLUSTER
    ////////////////////////

    PS::S32 getGroupPtclRemoteN() const{
        return n_group_member_remote_;
    }

    const PS::ReallocatableArray<PtclHard> & getPtcl() const {
        return ptcl_hard_;
    }

    PS::S32 getNCluster() const{
        return n_ptcl_in_cluster_.size();
    }

    PS::S32* getClusterNList() const{
        return n_ptcl_in_cluster_.getPointer();
    }

    PS::S32* getClusterNOffset() const{
        return n_ptcl_in_cluster_disp_.getPointer();
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
    template<class Teng>
    void CalcEnergyHard(PtclHard* ptcl, Teng & eng,  const PS::S32* ptcl_list, const PS::S32 ptcl_n, const PS::S32* group_list, const PS::S32 group_n, const PS::S32 nbin){
        eng.kin = eng.pot = eng.tot = 0.0;
        for(PS::S32 i=nbin; i<ptcl_n; i++){
            PtclHard* pi = &ptcl[ptcl_list[i]];
            eng.kin += 0.5 * pi->mass * pi->vel * pi->vel;

            for(PS::S32 j=i+1; j<ptcl_n; j++){
                PtclHard* pj = &ptcl[ptcl_list[j]];
                PS::F64vec rij = pi->pos - pj->pos;
                PS::F64 dr = sqrt(rij*rij + Int_pars_.eps2);
#ifdef INTEGRATED_CUTOFF_FUNCTION
                eng.pot -= pj->mass*pi->mass/dr*(1.0 - CalcW(dr/Int_pars_.rout, Int_pars_.rin/Int_pars_.rout));  
#else
                if(dr<Int_pars_.rout) eng.pot -= pj->mass*pi->mass*(1.0/dr*cutoff_pot(dr, Int_pars_.r_oi_inv, Int_pars_.r_A, Int_pars_.rin) - Int_pars_.pot_off);
#endif
            }

            for(PS::S32 j=0; j<group_n; j++){
                PtclHard* pj = &ptcl[group_list[j]];
                PS::F64vec rij = pi->pos - pj->pos;
                PS::F64 dr = sqrt(rij*rij + Int_pars_.eps2);
#ifdef INTEGRATED_CUTOFF_FUNCTION
                eng.pot -= pj->mass*pi->mass/dr*(1.0 - CalcW(dr/Int_pars_.rout, Int_pars_.rin/Int_pars_.rout));  
#else
                if(dr<Int_pars_.rout) eng.pot -= pj->mass*pi->mass*(1.0/dr*cutoff_pot(dr, Int_pars_.r_oi_inv, Int_pars_.r_A, Int_pars_.rin) - Int_pars_.pot_off);
#endif
            }
        }

        for(PS::S32 i=0; i<group_n; i++){
            PtclHard* pi = &ptcl[group_list[i]];
            eng.kin += 0.5 * pi->mass * pi->vel * pi->vel;

            for(PS::S32 j=i+1; j<group_n; j++){
                PtclHard* pj = &ptcl[group_list[j]];
                PS::F64vec rij = pi->pos - pj->pos;
                PS::F64 dr = sqrt(rij*rij + Int_pars_.eps2);
#ifdef INTEGRATED_CUTOFF_FUNCTION
                eng.pot -= pj->mass*pi->mass/dr*(1.0 - CalcW(dr/Int_pars_.rout, Int_pars_.rin/Int_pars_.rout));  
#else
                if(dr<Int_pars_.rout) eng.pot -= pj->mass*pi->mass*(1.0/dr*cutoff_pot(dr, Int_pars_.r_oi_inv, Int_pars_.r_A, Int_pars_.rin) - Int_pars_.pot_off);
#endif
            }
        }
        eng.tot = eng.kin + eng.pot;
    }

    template<class Teng, class TH4, class TARC, class Tgroup>
    void CalcEnergyHardFull(PtclHard* ptcl, Teng& E, Teng& AE, Teng& HE, Teng& ESD, TH4 &Hint, TARC& Aint, const Tgroup& group){
        Hint.CalcEnergy(HE);
        Teng TMP;
        Aint.EnergyRecord(TMP,true);
        Aint.EnergyRecord(AE);
        CalcEnergyHard(ptcl, E, group.getPtclList(), group.getPtclN(), group.getGroup(0), group.getGroupListSize(), group.getNumOfGroups());
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
#pragma omp for schedule(dynamic)
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
#pragma omp for schedule(dynamic)
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
#pragma omp for schedule(dynamic)
        for(PS::S32 i=0; i<n; i++){
            PS::S32 adr = ptcl_hard_[i].adr_org;
            //PS::S32 adr = adr_array[i];
#ifdef HARD_DEBUG
            assert(sys[adr].id == ptcl_hard_[i].id);
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
        n_ptcl_in_cluster_.resizeNoInitialize(n_cluster);
        n_ptcl_in_cluster_disp_.resizeNoInitialize(n_cluster+1);
        n_ptcl_in_cluster_disp_[0] = 0;
        for(PS::S32 i=0; i<n_cluster; i++){
            n_ptcl_in_cluster_[i] = _n_ptcl_in_cluster[i];
            n_ptcl_in_cluster_disp_[i+1] = n_ptcl_in_cluster_disp_[i] + n_ptcl_in_cluster_[i];
        }
        const PS::S32 n_ptcl = _adr_array.size();
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
    void driveForMultiCluster(const PS::F64 dt, Tsys & sys, const bool first_step_flag=false){
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
            driveForMultiClusterImpl(ptcl_hard_.getPointer(adr_head), n_ptcl, dt, extra_ptcl, first_step_flag);
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
    void driveForMultiClusterOMP(const PS::F64 dt, Tsys & sys, const bool first_step_flag=false){
        const PS::S32 n_cluster = n_ptcl_in_cluster_.size();
        const PS::S32 num_thread = PS::Comm::getNumberOfThread();
        PS::ReallocatableArray<PtclHard> extra_ptcl[num_thread];
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
            const PS::S32 ith = PS::Comm::getThreadNum();
#ifdef OMP_PROFILE
            time_thread[ith] -= PS::GetWtime();
#endif
            //const PS::S32 i   = n_sort_list[k].second;
            const PS::S32 adr_head = n_ptcl_in_cluster_disp_[i];
            const PS::S32 n_ptcl = n_ptcl_in_cluster_[i];
#ifdef OMP_PROFILE
            num_cluster[ith] += n_ptcl;
#endif
#ifdef HARD_DEBUG_PROFILE
            PS::F64 tstart = PS::GetWtime();
#endif
            driveForMultiClusterImpl(ptcl_hard_.getPointer(adr_head), n_ptcl, dt, extra_ptcl[ith], first_step_flag);
#ifdef OMP_PROFILE
            time_thread[ith] += PS::GetWtime();
#endif
#ifdef HARD_DEBUG_PROFILE
            PS::F64 tend = PS::GetWtime();
            std::cerr<<"HT: "<<i<<" "<<ith<<" "<<n_cluster<<" "<<n_ptcl<<" "<<tend-tstart<<std::endl;
#endif
        }
        if (n_cluster>0) {
            PS::S32 rank = PS::Comm::getRank();
            for(PS::S32 i=0; i<num_thread; i++) {
#ifdef OMP_PROFILE        
                std::cerr<<"thread: "<<i<<"  Hard Time="<<time_thread[i]<<"  n_ptcl="<<num_cluster[i]<<std::endl;
#endif
                for (PS::S32 j=0; j<extra_ptcl[i].size(); j++) {
                    PS::S32 adr = sys.getNumberOfParticleLocal();
                    sys.addOneParticle(Tsptcl(extra_ptcl[i][j],rank,adr));
#ifdef HARD_DEBUG
//                    if(sys[adr].id==10477) {
//                        std::cerr<<"Add particle adr="<<adr;
//                        sys[adr].print(std::cerr);
//                        std::cerr<<std::endl;
//                        std::cerr<<" original: ";
//                        extra_ptcl[i][j].print(std::cerr);
//                        std::cerr<<std::endl;
//                    }
                    if(extra_ptcl[i][j].id<0&&extra_ptcl[i][j].status<0) {
                        std::cerr<<"Error: extra particle list contain ghost particle! i_thread="<<i<<" index="<<j<<" rank="<<rank<<" adr="<<adr<<std::endl;
                        abort();
                    }
#endif
                }
            }
        }
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
                                                  r_bin_,
                                                  Int_pars_.rin,     
                                                  Int_pars_.rout,    
                                                  _dt_tree, 
                                                  id_offset_,
                                                  n_split_);

        // connected clusters
    }

    //template<class Tsys, class Tsptcl>
    //void initialMultiClusterOMP(Tsys & sys, const PS::F64 dt_tree){
    //    const PS::S32 n_cluster = n_ptcl_in_cluster_.size();
    //    //	const PS::S32 ith = PS::Comm::getThreadNum();
//#pragma omp for schedule(dynamic)
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
    
    void dump(const char* fname, const PS::F64 time_end, const PS::S32 first_int_flag, PtclHard* ptcl_bk, const PS::S32 n_ptcl) {
        std::FILE* fp = std::fopen(fname,"w");
        if (fp==NULL) {
            std::cerr<<"Error: filename ARC_dump.dat cannot be open!\n";
            abort();
        }
        fwrite(&time_end, sizeof(PS::F64),1,fp);
        fwrite(&first_int_flag, sizeof(PS::S32),1,fp);
        PtclHardDump(fp, ptcl_bk, n_ptcl);
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
        ARC_control_pert_.setA(Newtonian_AW<PtclHard,ARC_pert_pars>,Newtonian_extA_test<PtclHard,PtclH4*,PtclForce*,ARC_pert_pars>,Newtonian_timescale<ARC_pert_pars>);
        ARC_control_soft_.setA(Newtonian_AW<PtclHard,ARC_pert_pars>,Newtonian_extA_test<PtclHard,PtclH4*,PtclForce*,ARC_pert_pars>,Newtonian_timescale<ARC_pert_pars>);
#else
        ARC_control_pert_.setA(Newtonian_AW<PtclHard,ARC_pert_pars>,Newtonian_extA_pert<PtclHard,PtclH4*,PtclForce*,ARC_pert_pars>,Newtonian_timescale<ARC_pert_pars>);
        ARC_control_soft_.setA(Newtonian_AW<PtclHard,ARC_pert_pars>,Newtonian_extA_soft<PtclHard,PtclH4*,PtclForce*,ARC_pert_pars>,Newtonian_timescale<ARC_pert_pars>);
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

    void driveForMultiClusterOneDebug(PtclHard* _ptcl, const PS::S32 _n_ptcl,const PS::F64 _time_end, const PS::S32 _first_step_flag) {
        PS::ReallocatableArray<PtclHard> ptcl_new;
        driveForMultiClusterImpl(_ptcl, _n_ptcl, _time_end, ptcl_new, _first_step_flag);
    }

    void set_slowdown_factor(const PS::F64 _slowdown_factor) {
        sdfactor_ = _slowdown_factor;
    }

#endif

};
