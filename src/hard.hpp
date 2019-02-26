#pragma once
#ifdef USE_INTRINSIC_FOR_X86
#include<immintrin.h>
#endif

#include"cstdlib"

#include"AR/symplectic_integrator.h"
#include"Hermite/hermite_integrator.h"
#include"Hermite/hermite_particle.h"
#include"hard_ptcl.hpp"
#include"hard_force.hpp"
#include"hermite_interaction.hpp"
#include"hermite_information.hpp"
#include"hermite_perturber.hpp"
#include"ar_interaction.hpp"
#include"ar_perturber.hpp"
#include"cluster_list.hpp"


//! Hard integrator parameter manager
class HardManager{
public:
    PS::F64 energy_error_max;
    PS::F64 r_tidal_tensor;
    PS::F64 eps_sq;
    PS::S64 id_offset;
    PS::S32 n_split;
    ChangeOver changeover;
    H4::HermiteManager<HermiteInteraction> h4_manager;
    AR::SymplecticManager<ARInteraction> ar_manager;

    //! constructor
    HardManager(): energy_error_max(-1.0), r_tidal_tensor(-1.0), eps_sq(-1.0), id_offset(-1), n_split(-1), changeover(), h4_manager(), ar_manager() {
        h4_manager.interaction.changeover = &changeover;
        ar_manager.interaction.changeover = &changeover;
    }
    
    //! set softening
    void setEpsSq(const PS::F64 _eps_sq) {
        eps_sq = _eps_sq;
        h4_manager.interaction.eps_sq = _eps_sq;
        ar_manager.interaction.eps_sq = _eps_sq;
    }

    //! set gravitational constant
    void setG(const PS::F64 _g) {
        h4_manager.interaction.G = _g;
        ar_manager.interaction.G = _g;
    }

    //! set time step range
    void setDtRange(const PS::F64 _dt_max, const PS::S32 _dt_min_index) {
        h4_manager.step.setDtRange(_dt_max, _dt_min_index);
        ar_manager.time_step_real_min = h4_manager.step.getDtMin();
        ar_manager.time_error_max_real = 0.25*ar_manager.time_step_real_min;
    }

    //! check paramters
    bool checkParams() {
        ASSERT(energy_error_max>0.0);
        ASSERT(r_tidal_tensor>=0.0);
        ASSERT(eps_sq>=0.0);
        ASSERT(id_offset>0);
        ASSERT(n_split>0);
        ASSERT(changeover.checkParams());
        ASSERT(h4_manager.checkParams());
        ASSERT(ar_manager.checkParams());
        return true;
    }

    //! write class data to file with binary format
    /*! @param[in] _fp: FILE type file for output
     */
    void writeBinary(FILE *_fp) {
        size_t size = sizeof(*this) - sizeof(changeover) - sizeof(h4_manager) - sizeof(ar_manager);
        fwrite(this, size, 1, _fp);
        changeover.writeBinary(_fp);
        h4_manager.writeBinary(_fp);
        ar_manager.writeBinary(_fp);
    }

    //! read class data to file with binary format
    /*! @param[in] _fp: FILE type file for reading
     */
    void readBinary(FILE *_fin) {
        size_t size = sizeof(*this) - sizeof(changeover) - sizeof(h4_manager) - sizeof(ar_manager);
        size_t rcount = fread(this, size, 1, _fin);
        if (rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            abort();
        }
        changeover.readBinary(_fin);
        h4_manager.readBinary(_fin);
        ar_manager.readBinary(_fin);
        h4_manager.interaction.changeover = &changeover;
        ar_manager.interaction.changeover = &changeover;
    }

    //! print parameters
    void print(std::ostream & _fout) const{
        _fout<<"energy_error_max : "<<energy_error_max<<std::endl
             <<"r_tidal_tensor   : "<<r_tidal_tensor<<std::endl
             <<"eps_sq           : "<<eps_sq<<std::endl
             <<"id_offset        : "<<id_offset<<std::endl
             <<"n_split          : "<<n_split<<std::endl;
        changeover.print(_fout);
        h4_manager.print(_fout);
        ar_manager.print(_fout);
    }
};


//! Hard system
class SystemHard{
private:
    typedef H4::ParticleH4<PtclHard> PtclH4;
    // Notice: if new variables added, change pardump also
    PS::F64 time_origin_;
    
    PS::ReallocatableArray<PtclH4> ptcl_hard_;
    PS::ReallocatableArray<PS::S32> n_ptcl_in_cluster_;
    PS::ReallocatableArray<PS::S32> n_ptcl_in_cluster_disp_;
    PS::ReallocatableArray<PS::S32> n_group_in_cluster_;
    PS::ReallocatableArray<PS::S32> n_group_in_cluster_offset_;
    PS::ReallocatableArray<PS::S32> adr_first_ptcl_arti_in_cluster_;
    PS::S32 n_group_member_remote_; // number of members in groups but in remote nodes

    struct OPLessIDCluster{
        template<class T> bool operator() (const T & left, const T & right) const {
            return left.id_cluster < right.id_cluster;
        }
    };

public:
    HardManager* manager;

#ifdef PROFILE
    PS::S64 ARC_substep_sum;
    PS::S64 ARC_tsyn_step_sum;
    PS::F64 ARC_n_groups;
    PS::S64 H4_step_sum;
#endif
#ifdef HARD_CHECK_ENERGY
    PS::F64 hard_dE;
#endif

    //! check paramters
    bool checkParams() {
        ASSERT(manager!=NULL);
        ASSERT(manager->checkParams());
        return true;
    }

private:
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
                                                   PtclH4* _ptcl_local,
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
        PS::ReallocatableArray<PtclH4> ptcl_artifical[num_thread];

#pragma omp parallel for schedule(dynamic)
        for (PS::S32 i=0; i<n_cluster; i++){
            const PS::S32 ith = PS::Comm::getThreadNum();
            PtclH4* ptcl_in_cluster = _ptcl_local + _n_ptcl_in_cluster_disp[i];
            const PS::S32 n_ptcl = _n_ptcl_in_cluster[i];
            // reset status
            for(PS::S32 j=0; j<n_ptcl; j++) ptcl_in_cluster[j].status = 0;
            // search groups
            SearchGroup<PtclH4> group;
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
      @param[in] _acorr_flag: flag to do acorr for KDKDK_4TH case
     */
    template <class Tpsoft, class Ttree, class Tepj>
    void correctForceWithCutoffTreeNeighborOneParticleImp(Tpsoft& _psoft, 
                                                          Ttree& _tree,
                                                          const PS::F64 _rin,
                                                          const PS::F64 _rout,
                                                          const PS::F64 _r_oi_inv,
                                                          const PS::F64 _r_A,
                                                          const PS::F64 _eps_sq,
                                                          const bool _acorr_flag=false) {
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

#ifdef KDKDK_4TH
            if(_acorr_flag) 
                CalcAcorrShortWithLinearCutoff(pos_j, acc_j, _psoft.acorr,
                                             ptcl_nb[k].pos, ptcl_nb[k].acc, ptcl_nb[k].mass,
                                             _eps_sq, _r_oi_inv, _r_A, _rout, _rin);
            else
#endif
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
       @param[in] _acorr_flag: flag to do acorr for KDKDK_4TH case
     */
    template <class Tsys>
    void correctForceWithCutoffArtificalOneClusterImp(Tsys& _sys, 
                                                      const PtclH4* _ptcl_local,
                                                      const PS::S32 _adr_real_start,
                                                      const PS::S32 _adr_real_end,
                                                      const PS::S32 _n_group,
                                                      const PS::S32* _adr_first_ptcl_arti_in_cluster,
                                                      const PS::F64 _rin,
                                                      const PS::F64 _rout,
                                                      const PS::F64 _r_oi_inv,
                                                      const PS::F64 _r_A,
                                                      const PS::S32 _n_split,
                                                      const PS::F64 _eps_sq,
                                                      const bool _acorr_flag) {

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
#ifdef KDKDK_4TH
                    if(_acorr_flag) 
                        CalcAcorrShortWithLinearCutoff(pos_k, acc_k, _sys[k].acorr,
                                                     ptcl_kk.pos, ptcl_kk.acc, ptcl_kk.mass,
                                                     _eps_sq, _r_oi_inv, _r_A, _rout, _rin);
                    else
#endif
                        CalcAccPotShortWithLinearCutoff(pos_k, acc_k, pot_k, 
                                                        ptcl_kk.pos, ptcl_kk.mass, ptcl_kk.mass_bk, 
                                                        2, _eps_sq,
                                                        _r_oi_inv, _r_A, _rout, _rin);
                    }
                }

                // loop real particle
                for (int kj=_adr_real_start; kj<_adr_real_end; kj++) {
                    const PtclH4* ptcl_kj_ptr = &_ptcl_local[kj];
                    PS::S32 pot_control_flag = ptcl_kj_ptr->status<0? 1: 0;
#ifdef KDKDK_4TH
                    if(_acorr_flag) {
                        PS::S64 adr_kj = _ptcl_local[kj].adr_org;
                        PS::F64vec& pos_kj = _sys[adr_kj].pos;
                        PS::F64vec& acc_kj = _sys[adr_kj].acc;
                        // be careful for member particle with zero mass
                        PS::F64& mass_kj = pot_control_flag==1? _sys[adr_kj].mass_bk: _sys[adr_kj].mass;

                        CalcAcorrShortWithLinearCutoff(pos_k, acc_k, _sys[k].acorr,
                                                       pos_kj, acc_kj, mass_kj,
                                                       _eps_sq, _r_oi_inv, _r_A, _rout, _rin);
                    }
                    else
#endif
                        CalcAccPotShortWithLinearCutoff(pos_k, acc_k, pot_k, 
                                                        ptcl_kj_ptr->pos, ptcl_kj_ptr->mass, ptcl_kj_ptr->mass_bk, 
                                                        pot_control_flag, _eps_sq,
                                                        _r_oi_inv, _r_A, _rout, _rin);
                }
            }
            
            // for c.m. particle
            PS::F64vec& acc_cm = _sys[j_cm].acc;

            // substract c.m. force (acc) from tidal tensor force (acc)
            for (PS::S32 k=gpars.offset_tt; k<gpars.offset_orb; k++)  _sys[j_start+k].acc -= acc_cm;
                
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
            //    PtclH4* ptcl_k_ptr = &_ptcl_local[k];
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
       @param[in] _acorr_flag: flag to do acorr for KDKDK_4TH case
    */
    template <class Tsys>
    void correctForceWithCutoffClusterImp(Tsys& _sys, 
                                          const PtclH4* _ptcl_local,
                                          const PS::ReallocatableArray<PS::S32>& _n_ptcl_in_cluster,
                                          const PS::ReallocatableArray<PS::S32>& _n_ptcl_in_cluster_offset,
                                          const PS::ReallocatableArray<PS::S32>& _n_group_in_cluster,
                                          const PS::ReallocatableArray<PS::S32>& _n_group_in_cluster_offset,
                                          const PS::ReallocatableArray<PS::S32>& _adr_first_ptcl_arti_in_cluster,
                                          const PS::F64 _rin,
                                          const PS::F64 _rout,
                                          const PS::S32 _n_split,
                                          const PS::F64 _eps_sq,
                                          const bool _acorr_flag) {
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
            const PS::S32* adr_first_ptcl_arti = n_group>0? &_adr_first_ptcl_arti_in_cluster[_n_group_in_cluster_offset[i]] : NULL;

            // correction for artifical particles
            correctForceWithCutoffArtificalOneClusterImp(_sys, _ptcl_local, adr_real_start, adr_real_end, n_group, adr_first_ptcl_arti, _rin, _rout, r_oi_inv, r_A, _n_split, _eps_sq, _acorr_flag);

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
                    const PtclH4* ptcl_k_ptr = &_ptcl_local[k];
                    PS::S32 pot_control_flag = ptcl_k_ptr->status<0? 1: 0;
#ifdef KDKDK_4TH
                    if(_acorr_flag) {
                        PS::S64 adr_k = _ptcl_local[k].adr_org;
                        PS::F64vec& pos_k = _sys[adr_k].pos;
                        PS::F64vec& acc_k = _sys[adr_k].acc;
                        // be careful for member particle with zero mass
                        PS::F64& mass_k = pot_control_flag==1? _sys[adr_k].mass_bk: _sys[adr_k].mass;
                        CalcAcorrShortWithLinearCutoff(pos_j, acc_j, _sys[adr].acorr,
                                                       pos_k, acc_k, mass_k,
                                                       _eps_sq, r_oi_inv, r_A, _rout, _rin);
                    }
                    else
#endif
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
#ifdef KDKDK_4TH
                    if(_acorr_flag) 
                        CalcAcorrShortWithLinearCutoff(pos_j, acc_j, _sys[adr].acorr,
                                                     ptcl_k.pos, ptcl_k.acc, ptcl_k.mass, 
                                                     _eps_sq, r_oi_inv, r_A, _rout, _rin);
                    else
#endif
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
   @param[in] _acorr_flag: flag to do acorr for KDKDK_4TH case
*/
    template <class Tsys, class Tpsoft, class Ttree, class Tepj>
    void correctForceWithCutoffTreeNeighborAndClusterImp(Tsys& _sys,
                                                         Ttree& _tree,
                                                         const PtclH4* _ptcl_local,
                                                         const PS::ReallocatableArray<PS::S32>& _n_ptcl_in_cluster,
                                                         const PS::ReallocatableArray<PS::S32>& _n_ptcl_in_cluster_offset,
                                                         const PS::ReallocatableArray<PS::S32>& _n_group_in_cluster,
                                                         const PS::ReallocatableArray<PS::S32>& _n_group_in_cluster_offset,
                                                         const PS::ReallocatableArray<PS::S32>& _adr_first_ptcl_arti_in_cluster,
                                                         const PS::ReallocatableArray<PS::S32>& _adr_send,
                                                         const PS::F64 _rin,
                                                         const PS::F64 _rout,
                                                         const PS::S32 _n_split,
                                                         const PS::F64 _eps_sq,
                                                         const bool _acorr_flag=false) {

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
            correctForceWithCutoffArtificalOneClusterImp(_sys, _ptcl_local, adr_real_start, adr_real_end, n_group, adr_first_ptcl_arti, _rin, _rout, r_oi_inv, r_A, _n_split, _eps_sq, _acorr_flag);

            // obtain correction for real particles in clusters use tree neighbor search
            for (int j=adr_real_start; j<adr_real_end; j++) {
                PS::S64 adr = _ptcl_local[j].adr_org;
                // only do for local particles
#ifdef HARD_DEBUG
                if(adr>=0) assert(_sys[adr].id==_ptcl_local[j].id);
#endif
                if(adr>=0) correctForceWithCutoffTreeNeighborOneParticleImp<Tpsoft, Ttree, Tepj>(_sys[adr], _tree, _rin, _rout, r_oi_inv, r_A, _eps_sq, _acorr_flag);
            }
        }

        const PS::S32 n_send = _adr_send.size();
#pragma omp parallel for 
        // sending list to other nodes need also be corrected.
        for (int i=0; i<n_send; i++) {
            PS::S64 adr = _adr_send[i];
            correctForceWithCutoffTreeNeighborOneParticleImp<Tpsoft, Ttree, Tepj>(_sys[adr], _tree, _rin, _rout, r_oi_inv, r_A, _eps_sq, _acorr_flag); 
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
   @param[in] _acorr_flag: flag to do acorr for KDKDK_4TH case
*/
    template <class Tsys, class Tpsoft, class Ttree, class Tepj>
    void correctForceWithCutoffTreeNeighborImp(Tsys& _sys, 
                                               Ttree& _tree, 
                                               const PtclH4* _ptcl_local,
                                               const PS::S32 _n_ptcl,
                                               const PS::S32 _adr_ptcl_artifical_start,
                                               const PS::F64 _rin,
                                               const PS::F64 _rout,
                                               const PS::S32 _n_split,
                                               const PS::F64 _eps_sq,
                                               const bool _acorr_flag=false) { 
        // cutoff function parameter
        const PS::F64 r_oi_inv = 1.0/(_rout-_rin);
        const PS::F64 r_A = (_rout-_rin)/(_rout+_rin);

        // for real particle
#pragma omp parallel for schedule(dynamic)
        for (int i=0; i<_n_ptcl; i++) {
            PS::S64 adr = _ptcl_local[i].adr_org;
            if(adr>=0) correctForceWithCutoffTreeNeighborOneParticleImp<Tpsoft, Ttree, Tepj>(_sys[adr], _tree, _rin, _rout, r_oi_inv, r_A, _eps_sq, _acorr_flag);
        }

        // for artifical particle
        const PS::S32 n_tot = _sys.getNumberOfParticleLocal();
#pragma omp parallel for schedule(dynamic)
        for (int i=_adr_ptcl_artifical_start; i<n_tot; i++) 
            correctForceWithCutoffTreeNeighborOneParticleImp<Tpsoft, Ttree, Tepj>(_sys[i], _tree, _rin, _rout, r_oi_inv, r_A, _eps_sq, _acorr_flag);

        
        GroupPars gpars(_n_split);
#ifdef HARD_DEBUG
        assert((n_tot-_adr_ptcl_artifical_start)%gpars.n_ptcl_artifical==0);
#endif
#pragma omp parallel for schedule(dynamic)
        for (int i=_adr_ptcl_artifical_start; i<n_tot; i+=gpars.n_ptcl_artifical){
            PS::S32 i_cm = i + gpars.offset_cm;
            PS::F64vec& acc_cm = _sys[i_cm].acc;
            // substract c.m. force (acc) from tidal tensor force (acc)
            for (PS::S32 k=gpars.offset_tt; k<gpars.offset_orb; k++)  _sys[i+k].acc -= acc_cm;

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
            // substract c.m. force (acc) from tidal tensor force (acc)
            for (PS::S32 k=gpars.offset_tt; k<gpars.offset_orb; k++)  _sys[i+k].acc -= acc_cm;

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

#ifdef HARD_DEBUG
public:
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
       @param[in] _ithread: omp thread id, default 0
     */
    template <class Tsoft>
    void driveForMultiClusterImpl(PtclH4 * _ptcl_local,
                                  const PS::S32 _n_ptcl,
                                  Tsoft* _ptcl_artifical,
                                  const PS::S32 _n_group,
                                  const PS::F64 _time_end,
                                  const PS::S32 _ithread=0) {
        ASSERT(checkParams());
#ifdef HARD_CHECK_ENERGY
        std::map<PS::S32, PS::S32> N_count;  // counting number of particles in one cluster
        PS::F64 etoti, etotf;
#endif
#ifdef HARD_DEBUG_PROFILE
        N_count[_n_ptcl]++;
#endif
//#ifdef HARD_DEBUG_PRINT
//        PS::ReallocatableArray<PtclH4> ptcl_bk_pt;
//        ptcl_bk_pt.reserve(_n_ptcl);
//        ptcl_bk_pt.resizeNoInitialize(_n_ptcl);
//#endif

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
            gpars[i].init(manager->n_split);
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

        //// In orbital fitting soft perturbation, the status is used to identify which component the member belong to
        //for(int i=0; i<_n_group; i++) {
        //    // only first component is enough.
        //    for(int j=0; j<gpars[i].n_members_1st; j++)
        //        _ptcl_local[n_group_offset[i]+j].status = 0; 
        //}

        // pre-process for c.m. particle,
        for(int i=0; i<_n_group; i++){
            PS::S32 icm = adr_cm_ptcl[i];
            // kick c.m. (not done in previous kick function to avoid multi-kick)
            // Cannot do any kick in drift, because the K/D time step is not necessary same
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
                PS::S32 id_offset_j1 = _ptcl_artifical[adr_first_ptcl[i]+j].id - j/2- id_mem[0]*manager->n_split;
                // second member
                PS::S32 id_offset_j2 = _ptcl_artifical[adr_first_ptcl[i]+j+1].id - j/2 - id_mem[1]*manager->n_split;
                assert(id_offset_j1==id_offset_j2);
            }

            // check whether c.m. pos. and vel. are consistent
            PS::F64 mass_cm_check=0.0;
            // Cannot do velocity check because cm is not kicked
            //PS::F64vec vel_cm_check=PS::F64vec(0.0);
            PS::F64vec pos_cm_check=PS::F64vec(0.0);
            
            for(int j=0; j<gpars[i].n_members; j++) {
                PS::S32 k = n_group_offset[i]+j;
                mass_cm_check += _ptcl_local[k].mass;
                //vel_cm_check +=  _ptcl_local[k].vel*_ptcl_local[k].mass;
                pos_cm_check +=  _ptcl_local[k].pos*_ptcl_local[k].mass;
            }
            //vel_cm_check /= mass_cm_check;
            pos_cm_check /= mass_cm_check;

            assert(abs(mass_cm_check-_ptcl_artifical[icm].mass)<1e-10);
            //PS::F64vec dvec = vel_cm_check-_ptcl_artifical[icm].vel;
            PS::F64vec dpos = pos_cm_check-_ptcl_artifical[icm].pos;
            //assert(abs(dvec*dvec)<1e-20);
            assert(abs(dpos*dpos)<1e-20);
#endif

        }

#ifdef HARD_DEBUG_PRINT
        std::cerr<<"Hard: n_ptcl: "<<_n_ptcl<<" n_group: "<<_n_group<<std::endl;
#endif

        // manager
        H4::HermiteManager<HermiteInteraction>* h4_manager = &(manager->h4_manager);
        AR::SymplecticManager<ARInteraction>* ar_manager = &(manager->ar_manager);

        // Only one group with all particles in group
        if(_n_group==1&&n_single_init==0) {

            AR::SymplecticIntegrator<H4::ParticleAR<PtclHard>, PtclH4, ARPerturber, ARInteraction, H4::ARInformation<PtclHard>> sym_int;
            sym_int.manager = ar_manager;

            sym_int.particles.setMode(COMM::ListMode::copy);
            sym_int.particles.reserveMem(gpars[0].n_members);
            sym_int.info.reserveMem(gpars[0].n_members);
            sym_int.perturber.r_crit_sq = h4_manager->r_neighbor_crit*h4_manager->r_neighbor_crit;
            for (PS::S32 i=0; i<gpars[0].n_members; i++) {
                sym_int.particles.addMemberAndAddress(_ptcl_local[i]);
                sym_int.info.particle_index.addMember(i);
            }
            sym_int.reserveIntegratorMem();
            sym_int.info.generateBinaryTree(sym_int.particles);
            PS::S32 icm = adr_cm_ptcl[0];
            PS::S32 i_soft_pert_offset = gpars[0].offset_tt;
            TidalTensor tt;
            tt.fit(&_ptcl_artifical[i_soft_pert_offset], _ptcl_artifical[icm], manager->r_tidal_tensor, manager->n_split);
            sym_int.perturber.soft_pert=&tt;

            // calculate soft_pert_min
            sym_int.perturber.calcSoftPertMin(sym_int.info.getBinaryTreeRoot());
            
            // initialization 
            sym_int.initialIntegration(0.0);
            sym_int.info.calcDsAndStepOption(sym_int.slowdown.getSlowDownFactorOrigin(), ar_manager->step.getOrder()); 

#ifdef HARD_CHECK_ENERGY
            etoti = sym_int.getEtot();
#endif
            // integration
            sym_int.integrateToTime(sym_int.info.ds, _time_end, sym_int.info.fix_step_option);

            sym_int.particles.cm.pos += sym_int.particles.cm.vel * _time_end;

            // update rsearch
            sym_int.particles.cm.calcRSearch(_time_end);
            // copyback
            sym_int.particles.shiftToOriginFrame();
            sym_int.particles.template writeBackMemberAll<PtclH4>();

            for (PS::S32 i=0; i<gpars[0].n_members; i++) {
                _ptcl_local[i].r_search = sym_int.particles.cm.r_search;
            }

#ifdef PROFILE
            ARC_substep_sum += sym_int.profile.step_count;
            ARC_n_groups += 1;
#endif
#ifdef HARD_CHECK_ENERGY
            etotf  = sym_int.getEtot();
#endif
        }
        else {
            // integration -----------------------------
            H4::HermiteIntegrator<PtclHard, PtclH4, HermitePerturber, ARPerturber, HermiteInteraction, ARInteraction, HermiteInformation> h4_int;
            h4_int.manager = h4_manager;
            h4_int.ar_manager = ar_manager;

            h4_int.particles.setMode(COMM::ListMode::link);
            h4_int.particles.linkMemberArray(_ptcl_local, _n_ptcl);
            h4_int.particles.calcCenterOfMass();
            h4_int.particles.shiftToCenterOfMassFrame();
            
            h4_int.groups.setMode(COMM::ListMode::local);
            h4_int.groups.reserveMem(_n_group+3);
            h4_int.reserveIntegratorMem();

            // initial system 
            h4_int.initialSystemSingle(0.0);

            // Tidal tensor 
            TidalTensor tidal_tensor[_n_group+1];
            PS::S32 n_tt = 0;
            
            // add groups
            if (_n_group>0) {
                ASSERT(n_group_offset[_n_group]>0);
                PS::S32 ptcl_index_group[n_group_offset[_n_group]];
                for (PS::S32 i=0; i<n_group_offset[_n_group]; i++) ptcl_index_group[i] = i;
                h4_int.addGroups(ptcl_index_group, n_group_offset, _n_group);

                for (PS::S32 i=0; i<_n_group; i++) {
                    PS::S32 i_soft_pert_offset = adr_first_ptcl[i]+gpars[i].offset_tt;
                    PS::S32 icm = adr_cm_ptcl[i];
                    // correct pos for t.t. cm
                    _ptcl_artifical[icm].pos -= h4_int.particles.cm.pos;
                    tidal_tensor[i].fit(&_ptcl_artifical[i_soft_pert_offset], _ptcl_artifical[icm], manager->r_tidal_tensor, manager->n_split);
                    n_tt ++;
                    h4_int.groups[i].perturber.soft_pert = &tidal_tensor[i];

                    // calculate soft_pert_min
                    h4_int.groups[i].perturber.calcSoftPertMin(h4_int.groups[i].info.getBinaryTreeRoot());
                }
            }

            // initialization 
            h4_int.initialIntegration(); // get neighbors and min particles
            //h4_int.adjustGroups(true);
            //h4_int.initialIntegration();
            h4_int.sortDtAndSelectActParticle();
            h4_int.info.time = h4_int.getTime();

#ifdef HARD_CHECK_ENERGY
            h4_int.info.calcEnergy(h4_int.particles, h4_manager->interaction, true);
            etoti  = h4_int.info.etot0;
#endif
#ifdef HARD_DEBUG_PRINT
            h4_int.info.printColumnTitle(std::cout, WRITE_WIDTH);
            h4_int.particles.printColumnTitle(std::cout, WRITE_WIDTH);
            std::cout<<std::endl;
#endif

            // integration loop
            while (h4_int.info.time<_time_end) {

                
                h4_int.integrateOneStepAct();
                h4_int.adjustGroups(false);
                
                // check tt
                if (n_tt>0) {
                    const PS::S32* group_index = h4_int.getSortDtIndexGroup();
                    const PS::S32 n_init = h4_int.getNInitGroup();
                    for(int i=0; i<n_init; i++) {
                        auto& groupi = h4_int.groups[group_index[i]];
                        groupi.perturber.findCloseSoftPert(tidal_tensor, _n_group, groupi.particles.cm);
                        // calculate soft_pert_min
                        groupi.perturber.calcSoftPertMin(groupi.info.getBinaryTreeRoot());
                    }
                }
                // initial after groups are modified
                h4_int.initialIntegration();
                h4_int.sortDtAndSelectActParticle();
                h4_int.info.time = h4_int.getTime();

#ifdef HARD_DEBUG_PRINT
                //PS::F64 dt_max = 0.0;
                //PS::S32 n_group = h4_int.getNGroup();
                //PS::S32 n_single = h4_int.getNSingle();
                //if (n_group>0) dt_max = h4_int.groups[h4_int.getSortDtIndexGroup()[n_group-1]].particles.cm.dt;
                //if (n_single>0) dt_max = std::max(dt_max, h4_int.particles[h4_int.getSortDtIndexSingle()[n_single-1]].dt);
                //ASSERT(dt_max>0.0);
                if (fmod(h4_int.info.time, h4_manager->step.getDtMax()/1024)==0.0) {
                    h4_int.writeBackGroupMembers();
                    h4_int.info.calcEnergy(h4_int.particles, h4_manager->interaction, false);
            
                    h4_int.info.printColumn(std::cout, WRITE_WIDTH);
                    h4_int.particles.printColumn(std::cout, WRITE_WIDTH);
                    std::cout<<std::endl;
                }
                if (fmod(h4_int.info.time, h4_manager->step.getDtMax())==0.0) {
                    h4_int.printStepHist();
                }
#endif
            }
        
            h4_int.particles.cm.pos += h4_int.particles.cm.vel * _time_end;
            h4_int.writeBackGroupMembers();
            // update research
            const PS::S32* group_index = h4_int.getSortDtIndexGroup();
            for(PS::S32 i=0; i<h4_int.getNGroup(); i++) {
                const PS::S32 k =group_index[i];
                h4_int.groups[k].particles.cm.calcRSearch(_time_end);
                const PS::S32 n_member = h4_int.groups[k].particles.getSize();
                for (PS::S32 j=0; j<n_member; j++) {
                    h4_int.groups[k].particles.getMemberOriginAddress(j)->r_search = h4_int.groups[k].particles.cm.r_search;
                }
            }
            const PS::S32* single_index = h4_int.getSortDtIndexSingle();
            for (PS::S32 i=0; i<h4_int.getNSingle(); i++) {
                h4_int.particles[single_index[i]].calcRSearch(_time_end);
            }

#ifdef HARD_CHECK_ENERGY
            h4_int.info.calcEnergy(h4_int.particles, h4_manager->interaction, false);
            etotf  = h4_int.info.etot;
#endif
            h4_int.particles.shiftToOriginFrame();

#ifdef PROFILE
            //ARC_substep_sum += Aint.getNsubstep();
            H4_step_sum += h4_int.profile.hermite_single_step_count + h4_int.profile.hermite_group_step_count;
            ARC_substep_sum += h4_int.profile.ar_step_count;
            ARC_tsyn_step_sum += h4_int.profile.ar_step_count_tsyn;
            ARC_n_groups += _n_group;
#endif
#ifdef AR_DEBUG_PRINT
            for (PS::S32 i=0; i<h4_int.getNGroup(); i++) {
                const PS::S32 k= group_index[i];
                auto& groupk = h4_int.groups[k];
                std::cerr<<"Group N:"<<ARC_n_groups
                         <<" k:"<<k
                         <<" N_member: "<<groupk.particles.getSize()
                         <<" Slowdown: "<<groupk.slowdown.getSlowDownFactor()
                         <<" Slowdown(org): "<<groupk.slowdown.getSlowDownFactorOrigin()
                         <<" steps: "<<groupk.profile.step_count
                         <<" steps(tsyn): "<<groupk.profile.step_count_tsyn;
                auto& bin = groupk.info.getBinaryTreeRoot();
                std::cerr<<" semi: "<<bin.semi
                         <<" ecc: "<<bin.ecc
                         <<" NB: "<<groupk.perturber.neighbor_address.getSize()
                         <<std::endl;
            }
#endif
        }
#ifdef HARD_CHECK_ENERGY
        hard_dE = etotf - etoti;
#ifdef HARD_DEBUG_PRINT
        std::cerr<<"Hard Energy: init: "<<etoti<<" end: "<<etotf<<" dE: "<<hard_dE<<std::endl;
#endif        
        if (abs(hard_dE) > manager->energy_error_max) {
            std::cerr<<"Hard energy significant ("<<hard_dE<<") !\n";
            DATADUMP();
            abort();
        }
#endif
    }

public:

    SystemHard(){
        manager = NULL;
#ifdef HARD_DEBUG_PROFILE
        for(PS::S32 i=0;i<20;i++) N_count[i]=0;
#endif
#ifdef PROFILE
        ARC_substep_sum = 0;
        ARC_tsyn_step_sum =0;
        ARC_n_groups = 0;
        H4_step_sum = 0;
#endif
#ifdef HARD_CHECK_ENERGY
        hard_dE = 0;
#endif
        //        PS::S32 n_threads = PS::Comm::getNumberOfThread();
    }

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

    PS::ReallocatableArray<PtclH4> & getPtcl() {
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

    //void setParam(const PS::F64 _rbin,
    //              const PS::F64 _rout,
    //              const PS::F64 _rin,
    //              const PS::F64 _eps,
    //              const PS::F64 _dt_limit_hard,
    //              const PS::F64 _dt_min_hard,
    //              const PS::F64 _eta,
    //              const PS::F64 _time_origin,
    //              const PS::F64 _sd_factor,
    //              const PS::F64 _v_max,
    //              // const PS::F64 _gmin,
    //              // const PS::F64 _m_avarage,
    //              const PS::S64 _id_offset,
    //              const PS::S32 _n_split){
    //    /// Set chain pars (L.Wang)
	//    Int_pars_.rin  = _rin;
    //    Int_pars_.rout = _rout;
    //    Int_pars_.r_oi_inv = 1.0/(_rout-_rin);
    //    Int_pars_.r_A      = (_rout-_rin)/(_rout+_rin);
    //    Int_pars_.pot_off  = (1.0+Int_pars_.r_A)/_rout;
    //    Int_pars_.eps2  = _eps*_eps;
    //    Int_pars_.r_bin = _rbin;
    //    /// Set chain pars (L.Wang)        
    //    dt_limit_hard_ = _dt_limit_hard;
    //    dt_min_hard_   = _dt_min_hard;
    //    eta_s_ = _eta*_eta;
    //    sdfactor_ = _sd_factor;
    //    v_max_ = _v_max;
    //    time_origin_ = _time_origin;
//  //      gamma_ = std::pow(1.0/_gmin,0.33333);
    //    // r_search_single_ = _rsearch; 
    //    //r_bin_           = _rbin;
    //    // m_average_ = _m_avarage;
    //    manager->n_split = _n_split;
    //    id_offset_ = _id_offset;
    //}

//    void updateRSearch(PtclH4* ptcl_org,
//                       const PS::S32* ptcl_list,
//                       const PS::S32 n_ptcl,
//                       const PS::F64 dt_tree) {
//        for (PS::S32 i=0; i<n_ptcl; i++) {
//            ptcl_org[ptcl_list[i]].calcRSearch(dt_tree);
//        }
//    }

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


    //! integrate one isolated particle
    /*! integrate one isolated particle and calculate new r_search
      @param[in] _dt: tree time step
     */
    void driveForOneCluster(const PS::F64 _dt) {
        const PS::S32 n = ptcl_hard_.size();
        for(PS::S32 i=0; i<n; i++){
            PS::F64vec dr = ptcl_hard_[i].vel * _dt;
            ptcl_hard_[i].pos += dr;
            ptcl_hard_[i].calcRSearch(_dt);
            // ptcl_hard_[i].r_search= r_search_single_;
            /*
              DriveKeplerRestricted(mass_sun_, 
              pos_sun_, ptcl_hard_[i].pos, 
              vel_sun_, ptcl_hard_[i].vel, dt); 
            */
        }

    }

    //! integrate one isolated particle
    /*! integrate one isolated particle and calculate new r_search
      @param[in] _dt: tree time step
      @param[in] _v_max: maximum velocity used to calculate r_search
     */
    void driveForOneClusterOMP(const PS::F64 _dt) {
        const PS::S32 n = ptcl_hard_.size();
#pragma omp parallel for schedule(dynamic)
        for(PS::S32 i=0; i<n; i++){
            PS::F64vec dr = ptcl_hard_[i].vel * _dt;
            ptcl_hard_[i].pos += dr;
            ptcl_hard_[i].calcRSearch(_dt);
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
    template<class Tpsoft>
    void driveForMultiCluster(const PS::F64 dt, Tpsoft* _ptcl_soft){
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
#ifndef ONLY_SOFT
            const PS::S32 n_group = n_group_in_cluster_[i];
            Tpsoft* ptcl_artifical_ptr=NULL;
            if(n_group>0) ptcl_artifical_ptr = &(_ptcl_soft[adr_first_ptcl_arti_in_cluster_[n_group_in_cluster_offset_[i]]]);
#ifdef HARD_DUMP
            assert(hard_dump.size>0);
            hard_dump[0].backup(ptcl_hard_.getPointer(adr_head), n_ptcl, ptcl_artifical_ptr, n_group, dt, manager->n_split);
#endif
            driveForMultiClusterImpl(ptcl_hard_.getPointer(adr_head), n_ptcl, ptcl_artifical_ptr, n_group, dt);
#else
            auto* pi = ptcl_hard_.getPointer(adr_head);
            for (PS::S32 j=0; j<n_ptcl; j++) {
                PS::F64vec dr = pi[j].vel * dt;
                pi[j].pos += dr;
                pi[j].status = 0;
                pi[j].mass_bk = 0;
                pi[j].calcRSearch(dt);
            }
#endif
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

    template<class Tpsoft>
    void driveForMultiClusterOMP(const PS::F64 dt, Tpsoft* _ptcl_soft){
        const PS::S32 n_cluster = n_ptcl_in_cluster_.size();
        //PS::ReallocatableArray<PtclH4> extra_ptcl[num_thread];
        //// For test
        //PS::ReallocatableArray<std::pair<PS::S32,PS::S32>> n_sort_list;
        //n_sort_list.resizeNoInitialize(n_cluster);
        //for(PS::S32 i=0; i<n_cluster; i++) {
        //    n_sort_list[i].first = n_ptcl_in_cluster_[i];
        //    n_sort_list[i].second= i;
        //}
        //std::sort(n_sort_list.getPointer(),n_sort_list.getPointer()+n_cluster,[](const std::pair<PS::S32,PS::S32> &a, const std::pair<PS::S32,PS::S32> &b){return a.first<b.first;});
#ifdef OMP_PROFILE        
        const PS::S32 num_thread = PS::Comm::getNumberOfThread();
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
#ifndef ONLY_SOFT
            const PS::S32 n_group = n_group_in_cluster_[i];
            Tpsoft* ptcl_artifical_ptr=NULL;
            if(n_group>0) ptcl_artifical_ptr = &(_ptcl_soft[adr_first_ptcl_arti_in_cluster_[n_group_in_cluster_offset_[i]]]);
#ifdef OMP_PROFILE
            num_cluster[ith] += n_ptcl;
#endif
#ifdef HARD_DUMP
            assert(ith<hard_dump.size);
            hard_dump[ith].backup(ptcl_hard_.getPointer(adr_head), n_ptcl, ptcl_artifical_ptr, n_group, dt, manager->n_split);
#endif

#ifdef HARD_DEBUG_PROFILE
            PS::F64 tstart = PS::GetWtime();
#endif
            driveForMultiClusterImpl(ptcl_hard_.getPointer(adr_head), n_ptcl, ptcl_artifical_ptr, n_group, dt, ith);
#ifdef OMP_PROFILE
            time_thread[ith] += PS::GetWtime();
#endif
#ifdef HARD_DEBUG_PROFILE
            PS::F64 tend = PS::GetWtime();
            std::cerr<<"HT: "<<i<<" "<<ith<<" "<<n_cluster<<" "<<n_ptcl<<" "<<tend-tstart<<std::endl;
#endif
#else
            auto* pi = ptcl_hard_.getPointer(adr_head);
            for (PS::S32 j=0; j<n_ptcl; j++) {
                PS::F64vec dr = pi[j].vel * dt;
                pi[j].pos += dr;
                pi[j].status = 0;
                pi[j].mass_bk = 0;
                pi[j].calcRSearch(dt);
            }
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
                                                               manager->h4_manager.r_break_crit,
                                                               manager->changeover.getRin(),
                                                               manager->changeover.getRout(),
                                                               _dt_tree, 
                                                               manager->id_offset,
                                                               manager->n_split);

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
            _sys[k].pot_tot += _sys[k].mass / manager->changeover.getRout();
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
       @param[in] _acorr_flag: flag to do acorr for KDKDK_4TH case
    */
    template <class Tsys, class Tpsoft, class Ttree, class Tepj>
    void correctForceWithCutoffTreeNeighborAndClusterOMP(Tsys& _sys,
                                                         Ttree& _tree,
                                                         const PS::ReallocatableArray<PS::S32>& _adr_send,
                                                         const bool _acorr_flag=false) {
        correctForceWithCutoffTreeNeighborAndClusterImp<Tsys, Tpsoft, Ttree, Tepj>(_sys, _tree, ptcl_hard_.getPointer(), n_ptcl_in_cluster_, n_ptcl_in_cluster_disp_, n_group_in_cluster_, n_group_in_cluster_offset_, adr_first_ptcl_arti_in_cluster_, _adr_send, manager->changeover.getRin(), manager->changeover.getRout(), manager->n_split, manager->eps_sq, _acorr_flag);
    }

    //! Soft force correction due to different cut-off function
    /* Use cluster member, first correct for artifical particles, then for cluster member
       c.m. force is replaced by the averaged force on orbital particles
       Tidal tensor particle subtract the c.m. acc
       @param[in] _sys: global particle system, acc is updated
       @param[in] _acorr_flag: flag to do acorr for KDKDK_4TH case
*/
    template <class Tsys>
    void correctForceWithCutoffClusterOMP(Tsys& _sys, const bool _acorr_flag=false) { 
        correctForceWithCutoffClusterImp(_sys, ptcl_hard_.getPointer(), n_ptcl_in_cluster_, n_ptcl_in_cluster_disp_, n_group_in_cluster_, n_group_in_cluster_offset_, adr_first_ptcl_arti_in_cluster_, manager->changeover.getRin(), manager->changeover.getRout(), manager->n_split, manager->eps_sq, _acorr_flag);
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
        
        correctForceWithCutoffTreeNeighborImp<Tsys, Tpsoft, Ttree, Tepj>(_sys, _tree, _adr_ptcl_artifical_start, manager->changeover.getRin(), manager->changeover.getRout(), manager->n_split, manager->eps_sq);
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
    //        correctForceWithCutoffTreeNeighborOneParticleImp<Tpsoft, Ttree, Tepj>(_sys[adr], _tree, manager->changeover.getRin(), manager->changeover.getRout(), r_oi_inv, r_A, manager->eps_sq);
    //}

    //template<class Tsys, class Tsptcl>
    //void initialMultiClusterOMP(Tsys & sys, const PS::F64 dt_tree){
    //    const PS::S32 n_cluster = n_ptcl_in_cluster_.size();
    //    //	const PS::S32 ith = PS::Comm::getThreadNum();
//#pragma omp parallel for schedule(dynamic)
    //    for(PS::S32 i=0; i<n_cluster; i++){
    //        const PS::S32 adr_head = n_ptcl_in_cluster_disp_[i];
    //        const PS::S32 n_ptcl = n_ptcl_in_cluster_[i];
    //        SearchGroup<PtclH4> group;
    //        group.findGroups(ptcl_hard_.getPointer(adr_head), n_ptcl, manager->n_split);
    //        if (group.getPtclN()==2) group.searchAndMerge(ptcl_hard_.getPointer(adr_head), manager->changeover.getRout());
    //        else group.searchAndMerge(ptcl_hard_.getPointer(adr_head), manager->changeover.getRin());
    //        //group.searchAndMerge(ptcl_hard_.getPointer(adr_head), manager->changeover.getRin());
    //        PS::ReallocatableArray<PtclH4> ptcl_new;
    //        group.generateList(ptcl_hard_.getPointer(adr_head), ptcl_new, Int_pars_.r_bin, manager->changeover.getRin(), manager->changeover.getRout(), dt_tree, id_offset_, manager->n_split);
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


    //template <class Tpsoft>
    //void driveForMultiClusterOneDebug(PtclH4* _ptcl, const PS::S32 _n_ptcl, Tpsoft* _ptcl_artifical, const PS::S32 _n_group,  const PS::F64 _v_max, const PS::F64 _time_end) {
    //    driveForMultiClusterImpl(_ptcl, _n_ptcl, _ptcl_artifical, _n_group, _v_max, _time_end);
    //}

};

