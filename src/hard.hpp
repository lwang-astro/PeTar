#pragma once
#ifdef USE_INTRINSIC_FOR_X86
#include<immintrin.h>
#endif

#include"cstdlib"
#include <algorithm>

#include"AR/symplectic_integrator.h"
#include"Hermite/hermite_integrator.h"
#include"Hermite/hermite_particle.h"
#include"hard_ptcl.hpp"
#include"soft_ptcl.hpp"
#include"hermite_interaction.hpp"
#include"hermite_information.hpp"
#include"hermite_perturber.hpp"
#include"ar_interaction.hpp"
#include"ar_perturber.hpp"
#include"search_group_candidate.hpp"
#include"artificial_particles.hpp"
#include"stability.hpp"

typedef H4::ParticleH4<PtclHard> PtclH4;
typedef H4::ParticleAR<PtclHard> PtclAR;

//! Hard integrator parameter manager
class HardManager{
public:
    PS::F64 energy_error_max;
    PS::F64 eps_sq;
    PS::F64 r_in_base;
    PS::F64 r_out_base;
    PS::F64 n_step_per_orbit;
    ArtificialParticleManager ap_manager;
    H4::HermiteManager<HermiteInteraction> h4_manager;
    AR::TimeTransformedSymplecticManager<ARInteraction> ar_manager;

    //! constructor
    HardManager(): energy_error_max(-1.0), eps_sq(-1.0), r_in_base(-1.0), r_out_base(-1.0), n_step_per_orbit(-1.0), ap_manager(), h4_manager(), ar_manager() {}
    
    //! set softening
    void setEpsSq(const PS::F64 _eps_sq) {
        eps_sq = _eps_sq;
        h4_manager.interaction.eps_sq = _eps_sq;
        ar_manager.interaction.eps_sq = _eps_sq;
    }

    //! set gravitational constant
    void setGravitationalConstant(const PS::F64 _g) {
        ap_manager.gravitational_constant = _g;
#ifdef ORBIT_SAMPLING
        ap_manager.orbit_manager.gravitational_constant = _g;
#endif
        h4_manager.interaction.gravitational_constant = _g;
        ar_manager.interaction.gravitational_constant = _g;
    }

    //! set time step range
    void setDtRange(const PS::F64 _dt_max, const PS::S32 _dt_min_index) {
        h4_manager.step.setDtRange(_dt_max, _dt_min_index);
        ar_manager.time_step_min = h4_manager.step.getDtMin();
        ar_manager.time_error_max = 0.25*ar_manager.time_step_min;
    }

    //! check paramters
    bool checkParams() {
        ASSERT(energy_error_max>0.0);
        ASSERT(eps_sq>=0.0);
        ASSERT(r_in_base>0.0);
        ASSERT(r_out_base>0.0);
        ASSERT(n_step_per_orbit>0.0);
        ASSERT(ap_manager.checkParams());
        ASSERT(h4_manager.checkParams());
        ASSERT(ar_manager.checkParams());
        return true;
    }

    //! write class data to file with binary format
    /*! @param[in] _fp: FILE type file for output
     */
    void writeBinary(FILE *_fp) {
        size_t size = sizeof(*this) - sizeof(ap_manager) - sizeof(h4_manager) - sizeof(ar_manager);
        fwrite(this, size, 1, _fp);
        ap_manager.writeBinary(_fp);
        h4_manager.writeBinary(_fp);
        ar_manager.writeBinary(_fp);
    }

    //! read class data to file with binary format
    /*! @param[in] _fp: FILE type file for reading
     */
    void readBinary(FILE *_fin) {
        size_t size = sizeof(*this) - sizeof(ap_manager) - sizeof(h4_manager) - sizeof(ar_manager);
        size_t rcount = fread(this, size, 1, _fin);
        if (rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            abort();
        }
        ap_manager.readBinary(_fin);
        h4_manager.readBinary(_fin);
        ar_manager.readBinary(_fin);
    }

    //! print parameters
    void print(std::ostream & _fout) const{
        _fout<<"energy_error_max : "<<energy_error_max<<std::endl
             <<"eps_sq           : "<<eps_sq<<std::endl
             <<"r_in_base        : "<<r_in_base<<std::endl
             <<"r_out_base       : "<<r_out_base<<std::endl;
        ap_manager.print(_fout);
        h4_manager.print(_fout);
        ar_manager.print(_fout);
    }
};

struct HardEnergy{
    PS::F64 de;                 // energy error
    PS::F64 de_sd;              // slowdown energy error
    PS::F64 de_sd_change_cum;   // cumulative Etot_SD change due to the change of slowdown factor
    PS::F64 ekin_sd_correction; // correction from Ekin to Etot_sd
    PS::F64 epot_sd_correction; // correction from Epot to Etot_sd

    HardEnergy() {clear();}

    void clear() {
        de = de_sd = de_sd_change_cum = ekin_sd_correction = epot_sd_correction = 0.0;
    }

    void resetEnergyCorrection() {
        ekin_sd_correction = epot_sd_correction = 0.0;
    }

    HardEnergy& operator +=(const HardEnergy& _energy) {
        de    += _energy.de;
        de_sd += _energy.de_sd;
        de_sd_change_cum   += _energy.de_sd_change_cum;
        ekin_sd_correction += _energy.ekin_sd_correction;
        epot_sd_correction += _energy.epot_sd_correction;
        return *this;
    }

};

//! hard integrator 
class HardIntegrator{
public:
    H4::HermiteIntegrator<PtclHard, PtclH4, HermitePerturber, ARPerturber, HermiteInteraction, ARInteraction, HermiteInformation> h4_int; ///> hermite integrator
    AR::TimeTransformedSymplecticIntegrator<PtclAR, PtclH4, ARPerturber, ARInteraction, H4::ARInformation<PtclHard>> sym_int; ///> AR integrator
    HardManager* manager; ///> hard manager
    PS::ReallocatableArray<TidalTensor> tidal_tensor; ///> tidal tensor array
    PS::F64 time_origin;  ///> origin physical time
    PtclH4* ptcl_origin;  ///> original particle array

    AR::InterruptBinary<PtclAR> interrupt_binary; ///> interrupt binary address

#ifdef HARD_DEBUG_PRINT
    PS::ReallocatableArray<PS::S32> n_group_sub_init; ///> initial sub group number in each groups 
    PS::S32 n_group_sub_tot_init; ///> total sub groups in all groups initially
#endif

#ifdef PROFILE
    PS::S64 ARC_substep_sum;
    PS::S64 ARC_tsyn_step_sum;
    PS::S64 H4_step_sum;
#endif

#ifdef HARD_COUNT_NO_NEIGHBOR
    PS::ReallocatableArray<bool> table_neighbor_exist;
    PS::S32 n_neighbor_zero;
#endif

    bool use_sym_int;  ///> use AR integrator flag
    bool is_initialized; ///> indicator whether initialization is done

#ifdef HARD_CHECK_ENERGY
    HardEnergy energy;
#endif

    //! initializer
    HardIntegrator(): h4_int(), sym_int(), manager(NULL), tidal_tensor(), time_origin(-1.0), ptcl_origin(NULL), 
                      interrupt_binary(), 
#ifdef HARD_DEBUG_PRINT
                      n_group_sub_init(), n_group_sub_tot_init(0),
#endif
#ifdef PROFILE
                      ARC_substep_sum(0), ARC_tsyn_step_sum(0), H4_step_sum(0), 
#endif
#ifdef HARD_COUNT_NO_NEIGHBOR
                      table_neighbor_exist(), n_neighbor_zero(0),
#endif
                      use_sym_int(true), is_initialized(false) {
#ifdef HARD_CHECK_ENERGY
                          energy.clear();
#endif
                      }

    //! check parameters
    bool checkParams() {
        ASSERT(manager!=NULL);
        ASSERT(ptcl_origin!=NULL);
        ASSERT(time_origin>=0.0);
        ASSERT(is_initialized);
        return true;
    }

    //! initial integration
    /*!
       @param[in,out] _ptcl: particle array
       @param[in] _n_ptcl: particle number in cluster
       @param[in,out] _ptcl_artificial: artificial particle array, c.m. are kicked 
       @param[in] _n_group: group number in cluster
       @param[in] _n_member_in_group: number of members in each group
       @param[in] _manager: hard manager
       @param[in] _time_origin: initial physical time 
     */
    template <class Tsoft>
    void initial(PtclH4 * _ptcl,
                 const PS::S32 _n_ptcl,
                 Tsoft* _ptcl_artificial,
                 const PS::S32 _n_group,
                 const PS::S32* _n_member_in_group,
                 HardManager* _manager,
                 const PS::F64 _time_origin) {

        // ensure the integrator is not used
        ASSERT(ptcl_origin==NULL);
        ASSERT(manager==NULL);
        ASSERT(!is_initialized);

        is_initialized = true;

        // set manager
        manager = _manager;
        time_origin = _time_origin;
        ptcl_origin = _ptcl;

#ifdef HARD_DEBUG
        if (_n_ptcl>400) {
            std::cerr<<"Large cluster, n_ptcl="<<_n_ptcl<<" n_group="<<_n_group<<std::endl;
            for (PS::S32 i=0; i<_n_ptcl; i++) {
                if(ptcl_origin[i].r_search>10*ptcl_origin[i].r_search_min) {
                    std::cerr<<"i = "<<i<<" ";
                    ptcl_origin[i].print(std::cerr);
                    std::cerr<<std::endl;
                }
            }
        }
#endif

        // prepare initial groups with artificial particles
        PS::S32 adr_first_ptcl[_n_group+1];
        PS::S32 n_group_offset[_n_group+1]; // ptcl member offset in ptcl_origin
        n_group_offset[0] = 0;
        for(int i=0; i<_n_group; i++) 
            n_group_offset[i+1] = n_group_offset[i] + _n_member_in_group[i];
        
        auto& ap_manager = manager->ap_manager;
        if (_ptcl_artificial!=NULL) {
            for(int i=0; i<_n_group; i++) {
                adr_first_ptcl[i] = i*ap_manager.getArtificialParticleN();
                auto* pi = &(_ptcl_artificial[adr_first_ptcl[i]]);
#ifdef HARD_DEBUG
                assert(ap_manager.getMemberN(pi)==_n_member_in_group[i]);
#endif
                // pre-process for c.m. particle
                auto* pcm = ap_manager.getCMParticles(pi);
                // recover mass
                pcm->mass = pcm->group_data.artificial.getMassBackup();

#ifdef ARTIFICIAL_PARTICLE_DEBUG
                ap_manager.checkConsistence(&ptcl_origin[n_group_offset[i]], &(_ptcl_artificial[adr_first_ptcl[i]]));
#endif
            }

#ifdef HARD_DEBUG
            if(_n_group>0) {
                if(n_group_offset[_n_group]<_n_ptcl)
                    assert(ptcl_origin[n_group_offset[_n_group]].group_data.artificial.isSingle());
                if (_ptcl_artificial!=NULL) assert(ptcl_origin[n_group_offset[_n_group]-1].group_data.artificial.isMember());
                else assert(ptcl_origin[n_group_offset[_n_group]-1].group_data.artificial.isSingle());
            }
#endif
        }
#ifdef HARD_DEBUG
        else if(_n_group>0) assert(_n_ptcl==2); // right now only support isolated binary case without artificial particles
#endif

        // single particle start index in ptcl_origin
        PS::S32 i_single_start = n_group_offset[_n_group];
        // number of single particles
        PS::S32 n_single_init = _n_ptcl - i_single_start;
#ifdef HARD_DEBUG
        assert(n_single_init>=0);
#endif

#ifdef HARD_DEBUG
        // check member consistence
        for(int i=0; i<i_single_start; i++) {
            assert(ptcl_origin[i].group_data.artificial.isMember());
            assert(ptcl_origin[i].mass>0);
        }
#endif

#ifdef HARD_DEBUG_PRINT
        std::cerr<<"Hard: n_ptcl: "<<_n_ptcl<<" n_group: "<<_n_group<<std::endl;
#endif

        // manager
        auto& h4_manager = manager->h4_manager;
        auto& ar_manager = manager->ar_manager;
        
        // Only one group with all particles in group
        if(_n_group==1&&n_single_init==0) {
            use_sym_int = true;

            sym_int.manager = &ar_manager;

            sym_int.particles.setMode(COMM::ListMode::copy);
            const PS::S32 n_members = _n_member_in_group[0];
#ifdef HARD_DEBUG
            ASSERT(n_members == _n_ptcl);
#endif
            sym_int.particles.reserveMem(n_members);
            sym_int.info.reserveMem(n_members);
            //sym_int.perturber.r_crit_sq = h4_manager->r_neighbor_crit*h4_manager->r_neighbor_crit;
            for (PS::S32 i=0; i<n_members; i++) {
                sym_int.particles.addMemberAndAddress(ptcl_origin[i]);
                sym_int.info.particle_index.addMember(i);
                sym_int.info.r_break_crit = std::max(sym_int.info.r_break_crit,ptcl_origin[i].getRGroup());
                Float r_neighbor_crit = ptcl_origin[i].getRNeighbor();
                sym_int.perturber.r_neighbor_crit_sq = std::max(sym_int.perturber.r_neighbor_crit_sq, r_neighbor_crit*r_neighbor_crit);                
            }
            sym_int.reserveIntegratorMem();
            sym_int.info.generateBinaryTree(sym_int.particles, ar_manager.interaction.gravitational_constant);
#ifdef SOFT_PERT
            if (_ptcl_artificial!=NULL) {
                Tsoft* api=&(_ptcl_artificial[adr_first_ptcl[0]]);
                auto* apcm = ap_manager.getCMParticles(api);
                auto* aptt = ap_manager.getTidalTensorParticles(api);

                tidal_tensor.resizeNoInitialize(1);
                auto& tt = tidal_tensor[0];
                tt.fit(aptt, *apcm, ap_manager.r_tidal_tensor);
                sym_int.perturber.soft_pert=&tt;
                // set tt group id to the total number of particles
                sym_int.perturber.soft_pert->group_id = n_members;
            }
#endif
            // calculate soft_pert_min
            sym_int.perturber.calcSoftPertMin(sym_int.info.getBinaryTreeRoot(), ar_manager.interaction.gravitational_constant);
            
            // initialization 
            sym_int.initialIntegration(0.0);
            sym_int.info.calcDsAndStepOption(ar_manager.step.getOrder(),  ar_manager.interaction.gravitational_constant); 

            // calculate c.m. changeover
            auto& pcm = sym_int.particles.cm;
            PS::F64 m_fac = pcm.mass*Ptcl::mean_mass_inv;
            pcm.changeover.setR(m_fac, manager->r_in_base, manager->r_out_base);

#ifdef HARD_DEBUG
            if(_ptcl_artificial==NULL) {
                PS::F64 period = sym_int.info.getBinaryTreeRoot().period;
                PS::F64 sd_factor = sym_int.info.getBinaryTreeRoot().slowdown.getSlowDownFactor();
                PS::F64 sd_tmax = manager->ar_manager.slowdown_timescale_max;
                if (1.01*sd_factor*period<=sd_tmax) {
                    std::cerr<<"Warning: isolated binary SD ("<<sd_factor<<") * period ("<<period<<") = "<<period*sd_factor<<"< SD_timescale_max = dt_tree*n_step_per_orbit ("<<sd_tmax<<")"<<std::endl;
                }
            }
#endif
            //check paramters
            ASSERT(sym_int.info.checkParams());
            ASSERT(sym_int.perturber.checkParams());
        }
        else { // use hermite
            use_sym_int = false;

            h4_int.manager = &h4_manager;
            h4_int.ar_manager = &ar_manager;

            h4_int.particles.setMode(COMM::ListMode::link);
            h4_int.particles.linkMemberArray(ptcl_origin, _n_ptcl);

            h4_int.particles.calcCenterOfMass();
            h4_int.particles.shiftToCenterOfMassFrame();
            
            PS::S32 n_group_size_max = _n_ptcl+_n_group;
            h4_int.groups.setMode(COMM::ListMode::local);
            h4_int.groups.reserveMem(n_group_size_max);
            h4_int.reserveIntegratorMem();

            // initial system 
            h4_int.initialSystemSingle(0.0);

#ifdef SOFT_PERT
            // Tidal tensor 
            tidal_tensor.resizeNoInitialize(_n_group+1);
#endif
#ifdef HARD_COUNT_NO_NEIGHBOR
            table_neighbor_exist.resizeNoInitialize(_n_ptcl);
            for (int k=0; k<_n_ptcl; k++) table_neighbor_exist[k] = false;
#endif
            
            // add groups
            if (_n_group>0) {
                ASSERT(n_group_offset[_n_group]>0);
                PS::S32 ptcl_index_group[n_group_offset[_n_group]];
                for (PS::S32 i=0; i<n_group_offset[_n_group]; i++) ptcl_index_group[i] = i;
                h4_int.addGroups(ptcl_index_group, n_group_offset, _n_group);

                for (PS::S32 i=0; i<_n_group; i++) {
                    auto& groupi = h4_int.groups[i];

#ifdef SOFT_PERT
                    auto* api = &(_ptcl_artificial[adr_first_ptcl[i]]);
                    auto* apcm = ap_manager.getCMParticles(api);
                    auto* aptt = ap_manager.getTidalTensorParticles(api);

                    // correct pos for t.t. cm
                    apcm->pos -= h4_int.particles.cm.pos;

                    // fit tidal tensor
                    tidal_tensor[i].fit(aptt, *apcm, ap_manager.r_tidal_tensor);

                    // set tidal_tensor pointer
                    groupi.perturber.soft_pert = &tidal_tensor[i];

                    // set group id of tidal tensor to the n_members
                    groupi.perturber.soft_pert->group_id = groupi.particles.getSize();

                    // same tidal_tensor id to member particle group_data for identification later
                    for (PS::S32 k=0; k<groupi.particles.getSize(); k++) {
                        groupi.particles[k].setTidalTensorID(i+1);
                        ptcl_origin[n_group_offset[i]+k].setTidalTensorID(i+1);
                    }
#endif
                    // calculate soft_pert_min
                    groupi.perturber.calcSoftPertMin(groupi.info.getBinaryTreeRoot(), ar_manager.interaction.gravitational_constant);

                    // calculate c.m. changeover
                    auto& pcm = groupi.particles.cm;
                    PS::F64 m_fac = pcm.mass*Ptcl::mean_mass_inv;

                    ASSERT(m_fac>0.0);
                    pcm.changeover.setR(m_fac, manager->r_in_base, manager->r_out_base);

#ifdef HARD_DEBUG
                    PS::F64 r_out_cm = pcm.changeover.getRout();
                    for (PS::S32 k=0; k<groupi.particles.getSize(); k++) 
                        ASSERT(abs(groupi.particles[k].changeover.getRout()-r_out_cm)<1e-10);
#endif
                }
            }

            // initialization 
            h4_int.initialIntegration(); // get neighbors and min particles

#ifdef HARD_DEBUG_PRINT
            // AR inner slowdown number
            n_group_sub_init.resizeNoInitialize(_n_group);
            n_group_sub_tot_init=0;
            for (int i=0; i<_n_group; i++) {
#ifdef AR_SLOWDOWN_ARRAY
                n_group_sub_init[i] = h4_int.groups[i].binary_slowdown.getSize();
#else
                n_group_sub_init[i] = h4_int.groups[i].info.binarytree.getSize();
#endif
                n_group_sub_tot_init += n_group_sub_init[i];
            }
#endif

            h4_int.adjustGroups(true);

            const PS::S32 n_init = h4_int.getNInitGroup();
            const PS::S32* group_index = h4_int.getSortDtIndexGroup();
            for(int i=0; i<n_init; i++) {
                auto& groupi = h4_int.groups[group_index[i]];
                // calculate c.m. changeover
                auto& pcm = groupi.particles.cm;
                PS::F64 m_fac = pcm.mass*Ptcl::mean_mass_inv;
                ASSERT(m_fac>0.0);
                pcm.changeover.setR(m_fac, manager->r_in_base, manager->r_out_base);

/*  It is not consistent to use tidal tensor for different group
#ifdef SOFT_PERT                
                // check whether all r_out are same (primoridal or not)
                bool primordial_flag = true;
                PS::F64 r_out_cm = groupi.particles.cm.changeover.getRout();
                for (PS::S32 k=0; k<groupi.particles.getSize(); k++) 
                    if (abs(groupi.particles[k].changeover.getRout()-r_out_cm)>1e-10) {
                        primordial_flag =false;
                        break;
                    }
                if (_n_group>0 && primordial_flag) {
                    // check closed tt and only find consistent changeover 
                    PS::F32 tt_index=groupi.perturber.findCloseSoftPert(tidal_tensor, n_tt, n_group_size_max, groupi.particles.cm, r_out_cm);
                    ASSERT(tt_index<n_tt);
                    // calculate soft_pert_min
                    if (tt_index>=0) 
                        groupi.perturber.calcSoftPertMin(groupi.info.getBinaryTreeRoot(), ar_manager.interaction.gravitational_constant);
#ifdef HARD_DEBUG_PRINT
                    std::cerr<<"Find tidal tensor, group i: "<<group_index[i]<<" pcm.r_out: "<<r_out_cm;
                    std::cerr<<" member.r_out: ";
                    for (PS::S32 k=0; k<groupi.particles.getSize(); k++) 
                        std::cerr<<groupi.particles[k].changeover.getRout()<<" ";
                    std::cerr<<" tidal tensor index: "<<tt_index;
                    std::cerr<<std::endl;
#endif
                    tt_index=0;
                }
#endif
*/
            }

            h4_int.initialIntegration();
            h4_int.sortDtAndSelectActParticle();
            h4_int.info.time_origin = h4_int.getTime() + time_origin;

#ifdef HARD_CHECK_ENERGY
            h4_int.calcEnergySlowDown(true);
#endif

#ifdef HARD_DEBUG_PRINT_TITLE
            h4_int.printColumnTitle(std::cout, WRITE_WIDTH, n_group_sub_init.getPointer(), n_group_sub_init.size(), n_group_sub_tot_init);
            std::cout<<std::endl;
#endif
#ifdef HARD_DEBUG_PRINT
            h4_int.printColumn(std::cout, WRITE_WIDTH, n_group_sub_init.getPointer(), n_group_sub_init.size(), n_group_sub_tot_init);
            std::cout<<std::endl;
#endif
        }
    }

    //! Integrate system to time
    /*!
      @param [in] _time_end: time to integrate
      \return interrupt binary address, if no return NULL
     */
    AR::InterruptBinary<PtclAR>& integrateToTime(const PS::F64 _time_end) {
        ASSERT(checkParams());

        // dump flag
        bool dump_flag=false;

        // integration
        if (use_sym_int) {
            interrupt_binary = sym_int.integrateToTime(_time_end);
            if (manager->ar_manager.interrupt_detection_option==1) interrupt_binary.clear();
#ifdef PROFILE
            if (sym_int.profile.step_count>manager->ar_manager.step_count_max&&!dump_flag) {
                std::cerr<<"Large isolated AR step cluster found (dump): step: "<<sym_int.profile.step_count<<std::endl;
                DATADUMP("dump_large_step");
                dump_flag=true;
            }
#endif                
        }
        else {
#ifdef SOFT_PERT
            PS::S32 n_tt = tidal_tensor.size();
#endif
            // integration loop
            while (h4_int.getTime()<_time_end) {
                // integrate groups
                interrupt_binary = h4_int.integrateGroupsOneStep();
                // when binary is interrupted, break integration loop
                if (interrupt_binary.status!=AR::InterruptStatus::none) break;

#ifdef HARD_COUNT_NO_NEIGHBOR
                checkNeighborExist();
#endif

                // integrate singles
                h4_int.integrateSingleOneStepAct();
                h4_int.adjustGroups(false);

                const PS::S32 n_init_group = h4_int.getNInitGroup();
                const PS::S32* group_index = h4_int.getSortDtIndexGroup();
                for(int i=0; i<n_init_group; i++) {
                    auto& groupi = h4_int.groups[group_index[i]];
                    // calculate c.m. changeover
                    auto& pcm = groupi.particles.cm;
                    PS::F64 m_fac = pcm.mass*Ptcl::mean_mass_inv;
                    ASSERT(m_fac>0.0);
                    pcm.changeover.setR(m_fac, manager->r_in_base, manager->r_out_base);
                    

#ifdef SOFT_PERT                
                    // find corresponding tidal tensor if exist
                    PS::S32 tt_id_member = groupi.particles[0].getTidalTensorID();
#ifdef HARD_DEBUG
                    ASSERT(tt_id_member>=0&&tt_id_member<n_tt);
#endif
                    if (tt_id_member>0&&tt_id_member<n_tt) {
                        PS::S32 n_members = groupi.particles.getSize();
                        TidalTensor* tidal_tensor_i = &tidal_tensor[tt_id_member-1];

                        // check whether n member is consistent
                        if (int(tidal_tensor_i->group_id) == n_members) {
                            // check member group_data to find whether all member has the same tidal tensor id
                            bool tt_consistent = true;
                            for (PS::S32 k=1; k<n_members; k++) {
                                if (tt_id_member != groupi.particles[k].getTidalTensorID()) {
                                    tt_consistent = false;
                                    break;
                                }
                            }

                            // if all match, set group tidal tensor 
                            if (tt_consistent) groupi.perturber.soft_pert = tidal_tensor_i;
                        }
                    }

                    /* not used anymore 
                    // check whether all r_out are same (primoridal or not)
                    bool primordial_flag = true;
                    PS::F64 r_out_cm = groupi.particles.cm.changeover.getRout();
                    for (PS::S32 k=0; k<groupi.particles.getSize(); k++) 
                        if (abs(groupi.particles[k].changeover.getRout()-r_out_cm)>1e-10) {
                            primordial_flag =false;
                            break;
                        }

                    if (n_tt>0 && primordial_flag) {
                        // check closed tt and only find consistent changeover 
                        PS::F32 tt_index=groupi.perturber.findCloseSoftPert(tidal_tensor, n_tt, n_group_size_max, groupi.particles.cm, r_out_cm);
                        ASSERT(tt_index<n_tt);
                        // calculate soft_pert_min
                        if (tt_index>=0) 
                            groupi.perturber.calcSoftPertMin(groupi.info.getBinaryTreeRoot(),ar_manager.interaction.gravitational_constant);
#ifdef HARD_DEBUG_PRINT
                        std::cerr<<"Find tidal tensor, group i: "<<group_index[i]<<" pcm.r_out: "<<groupi.particles.cm.changeover.getRout();
                        std::cerr<<" member.r_out: ";
                        for (PS::S32 k=0; k<groupi.particles.getSize(); k++) 
                            std::cerr<<groupi.particles[k].changeover.getRout()<<" ";
                        std::cerr<<" tidal tensor index: "<<tt_index;
                        std::cerr<<std::endl;
#endif

                    }
                    */
#endif
                }
                ASSERT(n_init_group<=h4_int.getNActGroup());

/* Not consistent, should not shift
#ifdef SOFT_PERT
                // update c.m. for Tidal tensor
                if (n_tt>0) {
                    for(int i=n_init_group; i<n_act_group; i++) {
                        auto& groupi = h4_int.groups[group_index[i]];
                        if (groupi.perturber.soft_pert!=NULL) 
                            groupi.perturber.soft_pert->shiftCM(groupi.particles.cm.pos);
                    }
                }
#endif
*/
                // initial after groups are modified
                h4_int.initialIntegration();
                h4_int.sortDtAndSelectActParticle();
                h4_int.info.time_origin = h4_int.getTime() + time_origin;

#ifdef PROFILE
                if (h4_int.profile.ar_step_count>manager->ar_manager.step_count_max&&!dump_flag) {
                    std::cerr<<"Large H4-AR step cluster found (dump): step: "<<h4_int.profile.ar_step_count<<std::endl;
                    DATADUMP("dump_large_step");
                    dump_flag=true;
                } 
#endif                

#ifdef HARD_DEBUG_PRINT
                //PS::F64 dt_max = 0.0;
                //PS::S32 n_group = h4_int.getNGroup();
                //PS::S32 n_single = h4_int.getNSingle();
                //if (n_group>0) dt_max = h4_int.groups[h4_int.getSortDtIndexGroup()[n_group-1]].particles.cm.dt;
                //if (n_single>0) dt_max = std::max(dt_max, h4_int.particles[h4_int.getSortDtIndexSingle()[n_single-1]].dt);
                //ASSERT(dt_max>0.0);
                auto& h4_manager = manager->h4_manager;
                if (fmod(h4_int.getTime(), h4_manager.step.getDtMax()/HARD_DEBUG_PRINT_FEQ)==0.0) {
                    h4_int.calcEnergySlowDown(false);

                    h4_int.printColumn(std::cout, WRITE_WIDTH, n_group_sub_init.getPointer(), n_group_sub_init.size(), n_group_sub_tot_init);
                    std::cout<<std::endl;

                    PS::F64 de_sd   = h4_int.getEnergyErrorSlowDown();
                    if (abs(de_sd) > manager->energy_error_max) {
                        PS::F64 ekin    = h4_int.getEkin();
                        PS::F64 ekin_sd = h4_int.getEkinSlowDown();
                        PS::F64 epot    = h4_int.getEpot();
                        PS::F64 epot_sd = h4_int.getEpotSlowDown();
                        PS::F64 etot_sd = h4_int.getEtotSlowDownRef();
                        PS::F64 de      = h4_int.getEnergyError();
                        PS::F64 de_sd_change_cum = h4_int.getDESlowDownChangeCum();
                        std::cerr<<"Hard energy significant ("<<de_sd<<") !"
                                 <<"  Ekin: "<<ekin
                                 <<"  Epot: "<<epot
                                 <<"  Ekin_SD: "<<ekin_sd
                                 <<"  Epot_SD: "<<epot_sd
                                 <<"  Etot_SD_ref: "<<etot_sd
                                 <<"  dE: "<<de
                                 <<"  dE_SD: "<<de_sd
                                 <<"  dE_SD_CHANGE: "<<de_sd_change_cum
                                 <<std::endl;
                        DATADUMP("hard_dump");
                        ///abort();
                    }
                }
                if (fmod(h4_int.getTime(), h4_manager.step.getDtMax())==0.0) {
                    h4_int.printStepHist();
                }
#endif
            }

        }
        return interrupt_binary;
    }

#ifdef HARD_COUNT_NO_NEIGHBOR
    //! check whether the closest neighbor is < r_out
    void checkNeighborExist() {
        PS::S32 index_offset_group = h4_int.getIndexOffsetGroup();
        PS::S32 n_act_single = h4_int.getNActSingle();
        PS::S32* act_single_index = h4_int.getSortDtIndexSingle();
        for (int k=0; k<n_act_single; k++) {
            PS::S32 i = act_single_index[k];
            ASSERT(i>=0&&i<h4_int.particles.getSize());
            PS::S32 j = h4_int.neighbors[i].r_min_index;
            ASSERT(j<h4_int.particles.getSizeMax()+h4_int.groups.getSize());
            if (j<0) continue;
            PS::F64 rij2 = h4_int.neighbors[i].r_min_sq;
            PS::F64 r_out_i = h4_int.particles[i].changeover.getRout();
            PS::F64 r_out_j = (j<index_offset_group)? h4_int.particles[j].changeover.getRout() : h4_int.groups[j-index_offset_group].particles.cm.changeover.getRout();
            if (rij2<std::max(r_out_i, r_out_j)) table_neighbor_exist[i] = true;
        }
        PS::S32 n_act_group = h4_int.getNActGroup();
        PS::S32* act_group_index = h4_int.getSortDtIndexGroup();
        for (int k=0; k<n_act_group; k++) {
            PS::S32 i = act_group_index[k];
            ASSERT(i>=0&&i<h4_int.particles.getSizeMax()+h4_int.groups.getSize());
            PS::S32 j = h4_int.groups[i].perturber.r_min_index;
            ASSERT(j<h4_int.particles.getSizeMax()+h4_int.groups.getSize());
            if (j<0) continue;
            PS::F64 rij2 = h4_int.groups[i].perturber.r_min_sq;
            PS::F64 r_out_i = h4_int.groups[i].particles.cm.changeover.getRout();
            PS::F64 r_out_j = (j<index_offset_group)? h4_int.particles[j].changeover.getRout() : h4_int.groups[j-index_offset_group].particles.cm.changeover.getRout();
            if (rij2<std::max(r_out_i, r_out_j)) {
                for (int ki=0; ki<h4_int.groups[i].particles.getSize(); ki++) {
                    PS::S32 ki_index = h4_int.groups[i].info.particle_index[ki];
                    ASSERT(ki_index>=0&&ki_index<table_neighbor_exist.size());
                    table_neighbor_exist[ki_index] = true;
                }
            }
        }
    }
#endif

    //! drift c.m. particle of the cluster record group c.m. in group_data and write back data to original particle array
    /*!
       @param[in] _time_end: integration ending time (initial time is fixed to 0)
     */
    void driftClusterCMRecordGroupCMDataAndWriteBack(const PS::F64 _time_end) {
        ASSERT(checkParams());
#ifdef HARD_CHECK_ENERGY
        PS::F64 ekin, epot, ekin_sd, epot_sd, de, de_sd, de_sd_change_cum;
#endif
        if (use_sym_int) {
            auto& pcm = sym_int.particles.cm;
            pcm.pos += pcm.vel * _time_end;

            // update rsearch
            pcm.Ptcl::calcRSearch(_time_end);
            // copyback
            sym_int.particles.shiftToOriginFrame();
            sym_int.particles.template writeBackMemberAll<PtclH4>();

            PS::S32 n_members = sym_int.particles.getSize();

            for (PS::S32 i=0; i<n_members; i++) {
                auto& pi = ptcl_origin[i];
#ifdef STELLAR_EVOLUTION
                if (pi.mass==0.0) {
                    ASSERT(pi.group_data.artificial.isUnused());
                    continue;
                }

                // shift time interrupt in order to get consistent time for stellar evolution in the next drift
                pi.time_interrupt -= _time_end;
#endif

                pi.r_search = std::max(pcm.r_search, pi.r_search);
#ifdef CLUSTER_VELOCITY
                pi.group_data.cm.mass    = pcm.mass;
                pi.group_data.cm.vel[0]  = pcm.vel[0];
                pi.group_data.cm.vel[1]  = pcm.vel[1];
                pi.group_data.cm.vel[2]  = pcm.vel[2];
#endif
#ifdef HARD_DEBUG
                ASSERT(ptcl_origin[i].r_search>ptcl_origin[i].changeover.getRout());
#endif
            }


#ifdef PROFILE
            ARC_substep_sum += sym_int.profile.step_count;
#endif
#ifdef HARD_CHECK_ENERGY
            ekin    = sym_int.getEkin();
            epot    = sym_int.getEpot();
            de      = sym_int.getEnergyError();
#if (defined AR_SLOWDOWN_ARRAY) || (defined AR_SLOWDOWN_TREE)
            ekin_sd = sym_int.getEkinSlowDown();
            epot_sd = sym_int.getEpotSlowDown();
            de_sd   = sym_int.getEnergyErrorSlowDown();
#else
            ekin_sd = ekin;
            epot_sd = epot;
            de_sd   = de;
#endif
            de_sd_change_cum = sym_int.getDESlowDownChangeCum();
#endif
        }
        else {
#ifdef HARD_CHECK_ENERGY
            // PS: writeBackGroupMembers is done inside calcEnergySlowDown
            h4_int.calcEnergySlowDown(false);
            ekin    = h4_int.getEkin();
            ekin_sd = h4_int.getEkinSlowDown();
            epot    = h4_int.getEpot();
            epot_sd = h4_int.getEpotSlowDown();
            de      = h4_int.getEnergyError();
            de_sd   = h4_int.getEnergyErrorSlowDown();
            de_sd_change_cum = h4_int.getDESlowDownChangeCum();
#else
            h4_int.writeBackGroupMembers();
#endif
            h4_int.particles.cm.pos += h4_int.particles.cm.vel * _time_end;

            h4_int.particles.shiftToOriginFrame();

            // update research and group_data.cm
            auto& h4_pcm = h4_int.particles.cm;
            const PS::S32* group_index = h4_int.getSortDtIndexGroup();
            for(PS::S32 i=0; i<h4_int.getNGroup(); i++) {
                const PS::S32 k =group_index[i];
#ifdef HARD_DEBUG
                ASSERT(h4_int.groups[k].particles.cm.changeover.getRout()>0);
#endif
                //h4_int.groups[k].particles.cm.calcRSearch(_dt);
                auto& pcm = h4_int.groups[k].particles.cm;
                pcm.vel += h4_pcm.vel;

                //pcm.calcRSearch(h4_manager.interaction.G*(h4_pcm.mass-pcm.mass), abs(pcm.pot), h4_pcm.vel, _dt);
                pcm.Ptcl::calcRSearch(_time_end);
                const PS::S32 n_member = h4_int.groups[k].particles.getSize();
                //const PS::S32 id_first = h4_int.groups[k].particles.getMemberOriginAddress(0)->id;
                for (PS::S32 j=0; j<n_member; j++) {
                    auto* pj = h4_int.groups[k].particles.getMemberOriginAddress(j);
                    pj->r_search = std::max(pj->r_search, pcm.r_search);
#ifdef STELLAR_EVOLUTION
                    if (pj->mass==0.0) {
                        ASSERT(pj->group_data.artificial.isUnused());
                        continue;
                    }

                    // shift time interrupt in order to get consistent time for stellar evolution in the next drift
                    pj->time_interrupt -= _time_end;
#endif
#ifdef CLUSTER_VELOCITY
                    // save c.m. velocity and mass for neighbor search
                    pj->group_data.cm.mass    = pcm.mass;
                    pj->group_data.cm.vel[0]  = pcm.vel[0];
                    pj->group_data.cm.vel[1]  = pcm.vel[1];
                    pj->group_data.cm.vel[2]  = pcm.vel[2];
#endif
#ifdef HARD_DEBUG
                    ASSERT(pj->r_search>pj->changeover.getRout());
#endif
                }
            }

            const PS::S32* single_index = h4_int.getSortDtIndexSingle();
            for (PS::S32 i=0; i<h4_int.getNSingle(); i++) {
                auto& pi = h4_int.particles[single_index[i]];
#ifdef STELLAR_EVOLUTION
                if (pi.mass==0.0) {
                    ASSERT(pi.group_data.artificial.isUnused());
                    continue;
                }

                // shift time interrupt in order to get consistent time for stellar evolution in the next drift
                pi.time_interrupt -= _time_end;
#endif
#ifdef CLUSTER_VELOCITY
                // set group_data.cm to 0.0 for singles
                pi.group_data.cm.mass    = 0.0;
                pi.group_data.cm.vel[0]  = 0.0;
                pi.group_data.cm.vel[1]  = 0.0;
                pi.group_data.cm.vel[2]  = 0.0;
#endif
                pi.Ptcl::calcRSearch(_time_end);
//                pi.calcRSearch(h4_manager.interaction.G*(h4_pcm.mass-pi.mass), abs(pi.pot), h4_pcm.vel, _dt);
            }


#ifdef PROFILE
            //ARC_substep_sum += Aint.getNsubstep();
            H4_step_sum += h4_int.profile.hermite_single_step_count + h4_int.profile.hermite_group_step_count;
            ARC_substep_sum += h4_int.profile.ar_step_count;
            ARC_tsyn_step_sum += h4_int.profile.ar_step_count_tsyn;

            if (h4_int.profile.ar_step_count>manager->ar_manager.step_count_max) {
                std::cerr<<"Large AR step cluster found: total step: "<<h4_int.profile.ar_step_count<<std::endl;
                //DATADUMP("dump_large_step");
            } 
#endif

#ifdef HARD_COUNT_NO_NEIGHBOR
            for (PS::S32 i=0; i<table_neighbor_exist.size(); i++) {
                if(!table_neighbor_exist[i]) n_neighbor_zero++;
            }
#endif

#ifdef AR_DEBUG_PRINT
            for (PS::S32 i=0; i<h4_int.getNGroup(); i++) {
                const PS::S32 k= group_index[i];
                auto& groupk = h4_int.groups[k];
                auto& bink = groupk.info.getBinaryTreeRoot();
                std::cerr<<"Group k:"<<std::setw(2)<<k
                         <<" N_member: "<<std::setw(4)<<groupk.particles.getSize()
                         <<" step: "<<std::setw(12)<<groupk.profile.step_count_sum
                         <<" step(tsyn): "<<std::setw(10)<<groupk.profile.step_count_tsyn_sum
//                         <<" step(sum): "<<std::setw(12)<<h4_int.profile.ar_step_count
//                         <<" step_tsyn(sum): "<<std::setw(12)<<h4_int.profile.ar_step_count_tsyn
                         <<" Soft_Pert: "<<std::setw(20)<<groupk.perturber.soft_pert_min
                         <<" Pert_In: "<<std::setw(20)<<bink.slowdown.getPertIn()
                         <<" Pert_Out: "<<std::setw(20)<<bink.slowdown.getPertOut()
                         <<" SD: "<<std::setw(20)<<bink.slowdown.getSlowDownFactor()
                         <<" SD(org): "<<std::setw(20)<<bink.slowdown.getSlowDownFactorOrigin()
                         <<" semi: "<<std::setw(20)<<bink.semi
                         <<" ecc: "<<std::setw(20)<<bink.ecc
                         <<" period: "<<std::setw(20)<<bink.period
                         <<" NB: "<<std::setw(4)<<groupk.perturber.neighbor_address.getSize()
                         <<std::endl;
                if (groupk.profile.step_count_tsyn_sum>10000) {
                    DATADUMP("hard_dump");
                }
            }
#endif
        }

#ifdef HARD_CHECK_ENERGY
        energy.de_sd += de_sd;
        energy.de    += de;
        PS::F64 ekin_sd_correction = ekin_sd - ekin;
        PS::F64 epot_sd_correction = epot_sd - epot;
        energy.ekin_sd_correction += ekin_sd_correction;
        energy.epot_sd_correction += epot_sd_correction;
        energy.de_sd_change_cum   += de_sd_change_cum - (ekin_sd_correction + epot_sd_correction);
#ifdef HARD_DEBUG_PRINT
        std::cerr<<"Hard Energy: "
                 <<"  Ekin: "<<ekin
                 <<"  Ekin_SD: "<<ekin_sd
                 <<"  Epot: "<<epot_sd
                 <<"  Epot_SD: "<<epot_sd
                 <<"  dE: "<<de
                 <<"  dE_SD: "<<de_sd
                 <<"  dE_SD_CHANGE: "<<de_sd_change_cum
                 <<"  H4_step_sum: "<<H4_step_sum
                 <<"  ARC_substep_sum: "<<ARC_substep_sum
                 <<"  ARC_tsyn_step_sum: "<<ARC_tsyn_step_sum
                 <<std::endl;
#endif        
        if (abs(de_sd) > manager->energy_error_max) {
            std::cerr<<"Hard energy significant ("<<de_sd<<") !\n";
            DATADUMP("hard_large_energy");
            //abort();
        }
#endif

    }

    //! print interrupt binary information
    /*!
      @param[out] _fout: std::ostream output object
    */
    void printInterruptBinaryInfo(std::ostream & _fout) const{
        _fout<<"Interrupt condition triggered! Time: ";
        if (use_sym_int) _fout<<"(AR) ";
        else _fout<<"(Hermite) ";
        _fout<<interrupt_binary.time_now<<std::endl;
        interrupt_binary.adr->printColumnTitle(_fout);
        _fout<<std::endl;
        interrupt_binary.adr->printColumn(_fout);
        _fout<<std::endl;
        PtclAR::printColumnTitle(_fout);
        _fout<<std::endl;
        for (int j=0; j<2; j++) {
            interrupt_binary.adr->getMember(j)->printColumn(_fout);
            _fout<<std::endl;
        }
    }

    //! clear function
    void clear() {
        sym_int.clear();
        h4_int.clear();
        manager = NULL;
        tidal_tensor.resizeNoInitialize(0);
        time_origin = 0;
        ptcl_origin = NULL;
        interrupt_binary.clear();
        is_initialized = false;

#ifdef PROFILE
        ARC_substep_sum = 0;
        ARC_tsyn_step_sum = 0;
        H4_step_sum = 0;
#endif
#ifdef HARD_COUNT_NO_NEIGHBOR
        table_neighbor_exist.resizeNoInitialize(0);
        n_neighbor_zero = 0;
#endif
#ifdef HARD_CHECK_ENERGY
        energy.clear();
#endif

    }       

};

//! Hard system
class SystemHard{
private:
    PS::F64 time_origin_; // physical origin time
    PS::F64 time_write_back_; // time of writing back data
    
    PS::ReallocatableArray<PtclH4> ptcl_hard_;                        // particle data
    PS::ReallocatableArray<PS::S32> n_ptcl_in_cluster_;               // number of particles in one cluster
    PS::ReallocatableArray<PS::S32> n_ptcl_in_cluster_disp_;          // boundary of particle cluster
    PS::ReallocatableArray<PS::S32> n_group_in_cluster_;              // number of groups in one cluster
    PS::ReallocatableArray<PS::S32> n_group_in_cluster_offset_;       // boundary of groups in n_member_in_group_ and adr_first_ptcl_arti_in_cluster_
    PS::ReallocatableArray<PS::S32> n_member_in_group_;               // number of members in each group
    PS::ReallocatableArray<PS::S32> n_member_in_group_offset_;        // number of members in group index offest in each cluster (for n_meber_in_group_, notice in each cluster, the first offset is zero)
    PS::ReallocatableArray<PS::S32> adr_first_ptcl_arti_in_cluster_;  // address of the first artificial particle in each groups
    PS::ReallocatableArray<PS::S32> i_cluster_changeover_update_;     // cluster index that has member need changeover update
    PS::S32 n_group_member_remote_; // number of members in groups but in remote nodes

    HardIntegrator* hard_int_; ///> hard integrator array
    PS::S32 n_hard_int_max_; ///> array size of hard_int
    PS::S32 n_hard_int_use_; ///> number of used hard integrator
    PS::ReallocatableArray<HardIntegrator*> interrupt_list_; ///> interrupt integrator list
    PS::F64 interrupt_dt_; ///> time end record for interrupt clusters;

    struct OPLessIDCluster{
        template<class T> bool operator() (const T & left, const T & right) const {
            return left.id_cluster < right.id_cluster;
        }
    };


public:
    PS::ReallocatableArray<COMM::BinaryTree<PtclH4,COMM::Binary>> binary_table;
    HardManager* manager;

#ifdef PROFILE
    PS::S64 ARC_substep_sum;
    PS::S64 ARC_tsyn_step_sum;
    PS::S64 ARC_n_groups;
    PS::S64 ARC_n_groups_iso;
    PS::S64 H4_step_sum;
#endif
#ifdef HARD_COUNT_NO_NEIGHBOR
    PS::S64 n_neighbor_zero;
#endif
#ifdef HARD_CHECK_ENERGY
    HardEnergy energy;
#endif

    //! check paramters
    bool checkParams() {
        ASSERT(manager!=NULL);
        ASSERT(manager->checkParams());
        ASSERT(hard_int_!=NULL);
        ASSERT(n_hard_int_max_>0);
        return true;
    }

private:

    //! collect member particle index, set member parameters
    /*! Set type to member (backup mass) and mass to zero, update changeover and rsearch
      Collect member particle address (index) in _par.group_list;
      backup mass of member particle and set mass to zero
      @param[in,out] _par: parameter container 
      @param[in,out] _ptcl: member particle
     */
    template <class Tchp, class Tptcl>
    static void collectGroupMemberAdrAndSetMemberParametersIter(Tchp& _par, Tptcl*& _ptcl) {
        _par.group_list[_par.n++] = _ptcl - _par.adr_ref;
#ifdef HARD_DEBUG
        assert(_ptcl->mass>0.0);
#endif
        if (_par.pcm_id<0) { // isolated binary case
            _ptcl->group_data.artificial.setParticleTypeToMember(_ptcl->mass);
        }
        else { // with artificial particles
            _ptcl->group_data.artificial.setParticleTypeToMember(_ptcl->mass, - (_par.pcm_id+1));
            _ptcl->mass = 0.0;
        }
        PS::F64 rin_cm = _par.changeover_cm->getRin();
        PS::F64 rin_p  = _ptcl->changeover.getRin();
        if (rin_p!=rin_cm) {
            // avoid round-off error case
            if (abs(rin_p-rin_cm)<1e-10) {
                _ptcl->changeover = *_par.changeover_cm;
            }
            else {
                _ptcl->changeover.r_scale_next = rin_cm/rin_p;
                _ptcl->r_search = std::max(_ptcl->r_search, _par.rsearch_cm);
#ifdef ARTIFICIAL_PARTICLE_DEBUG
                if (_ptcl->r_search<rin_p) {
                    std::cerr<<"Error, r_search<r_out found! \n"
                             <<"id: "<<_ptcl->id
                             <<"bin_changeover: ";
                    _par.changeover_cm->print(std::cerr);
                    std::cerr<<"member changeover: ";
                    _ptcl->changeover.print(std::cerr);
                    std::cerr<<"member mass: "<<_ptcl->mass
                             <<"member r_search: "<<_ptcl->r_search
                             <<"cm research: "<<_par.rsearch_cm
                             <<std::endl;
                    abort();
                }
#endif
                _par.changeover_update_flag = true;
            }
        }
    }


    //! calculate binary parameters
    /*! get new changeover, rsearch, id for c.m.
     */
    template <class Tchp, class Tptcl>
    static PS::S64 calcBinaryIDChangeOverAndRSearchIter (Tchp& _par, const PS::S64& _id1, const PS::S64& _id2, COMM::BinaryTree<Tptcl,COMM::Binary>& _bin) {
        // set bin id as the left member id
        // _id1==-1 is the initial status, once id is obtained, it is bin.id of left member
        if (_id1<0) _bin.id = _bin.getLeftMember()->id;
        else _bin.id = _id1;
        if (_id2<0) _bin.id = std::min(_bin.id, _bin.getRightMember()->id);
        else _bin.id = std::min(_bin.id, _id2);
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(_bin.id>0);
#endif
        //set changeover to the same (root) one
        _bin.changeover.setR(_bin.mass*_par.mean_mass_inv, _par.rin, _par.rout);

        // calculate rsearch
        _bin.Ptcl::calcRSearch(_par.dt_tree);
        
        return _bin.id;
    }

    //! correct force and potential for soft force with changeover function
    /*!
      @param[in,out] _pi: particle for correction
      @param[in] _pj: j particle to calculate correction
     */
    template <class Tpi>
    inline void calcAccPotShortWithLinearCutoff(Tpi& _pi,
                                                const Ptcl& _pj) {
        PS::F64 G = ForceSoft::grav_const;
        const PS::F64vec dr = _pi.pos - _pj.pos;
        const PS::F64 dr2 = dr * dr;
        const PS::F64 dr2_eps = dr2 + manager->eps_sq;
        const PS::F64 drinv = 1.0/sqrt(dr2_eps);
        const PS::F64 gmor = G*_pj.mass * drinv;
        const PS::F64 drinv2 = drinv * drinv;
        const PS::F64 gmor3 = gmor * drinv2;
        const PS::F64 dr_eps = drinv * dr2_eps;
        const PS::F64 k = 1.0 - ChangeOver::calcAcc0WTwo(_pi.changeover, _pj.changeover, dr_eps);

        // linear cutoff 
        const PS::F64 r_out = manager->r_out_base;
        const PS::F64 r_out2 = r_out * r_out;
        const PS::F64 dr2_max = (dr2_eps > r_out2) ? dr2_eps : r_out2;
        const PS::F64 drinv_max = 1.0/sqrt(dr2_max);
        const PS::F64 gmor_max = G*_pj.mass * drinv_max;
        const PS::F64 drinv2_max = drinv_max*drinv_max;
        const PS::F64 gmor3_max = gmor_max * drinv2_max;

        auto& pj_artificial = _pj.group_data.artificial;
#ifdef ONLY_SOFT
        const PS::F64 kpot  = 1.0 - ChangeOver::calcPotWTwo(_pi.changeover, _pj.changeover, dr_eps);
        // single, remove linear cutoff, obtain changeover soft potential
        if (pj_artificial.isSingle()) _pi.pot_tot -= dr2_eps>r_out2? 0.0: (gmor*kpot  - gmor_max);   
        // member, mass is zero, use backup mass
        else if (pj.artificial.isMember()) _pi.pot_tot -= dr2_eps>r_out2? 0.0: (G*pj_artificial.getMassBackup()*drinv*kpot  - gmor_max);   
        // (orbitial) artificial, should be excluded in potential calculation, since it is inside neighbor, gmor_max cancel it to 0.0
        else _pi.pot_tot += gmor_max; 
#else
        // single/member, remove linear cutoff, obtain total potential
        if (pj_artificial.isSingle()) _pi.pot_tot -= (gmor - gmor_max);   
        // member, mass is zero, use backup mass
        else if (pj_artificial.isMember()) _pi.pot_tot -= (G*pj_artificial.getMassBackup()*drinv  - gmor_max);   
        // (orbitial) artificial, should be excluded in potential calculation, since it is inside neighbor, gmor_max cancel it to 0.0
        else _pi.pot_tot += gmor_max; 
#endif
        // correct to changeover soft acceleration
        _pi.acc -= (gmor3*k - gmor3_max)*dr;
    }

    //! correct force and potential for soft force with changeover function
    /*!
      @param[in,out] _pi: particle for correction
      @param[in] _pj: j particle to calculate correction
     */
    template <class Tpi>
    inline void calcAccPotShortWithLinearCutoff(Tpi& _pi,
                                                const EPJSoft& _pj) {
        PS::F64 G = ForceSoft::grav_const;
        const PS::F64vec dr = _pi.pos - _pj.pos;
        const PS::F64 dr2 = dr * dr;
#ifdef HARD_DEBUG
        assert(dr2>0.0);
#endif
        const PS::F64 dr2_eps = dr2 + manager->eps_sq;
        const PS::F64 r_out = manager->r_out_base;
        const PS::F64 r_out2 = r_out * r_out;
        const PS::F64 drinv = 1.0/sqrt(dr2_eps);
        const PS::F64 gmor = G*_pj.mass * drinv;
        const PS::F64 drinv2 = drinv * drinv;
        const PS::F64 gmor3 = gmor * drinv2;
        const PS::F64 dr_eps = drinv * dr2_eps;
        ChangeOver chj;
        chj.setR(_pj.r_in, _pj.r_out);
        const PS::F64 k = 1.0 - ChangeOver::calcAcc0WTwo(_pi.changeover, chj, dr_eps);

        // linear cutoff 
        const PS::F64 dr2_max = (dr2_eps > r_out2) ? dr2_eps : r_out2;
        const PS::F64 drinv_max = 1.0/sqrt(dr2_max);
        const PS::F64 gmor_max = G*_pj.mass * drinv_max;
        const PS::F64 drinv2_max = drinv_max*drinv_max;
        const PS::F64 gmor3_max = gmor_max * drinv2_max;

        auto& pj_artificial = _pj.group_data.artificial;
#ifdef ONLY_SOFT
        const PS::F64 kpot  = 1.0 - ChangeOver::calcPotWTwo(_pi.changeover, chj, dr_eps);
        // single, remove linear cutoff, obtain changeover soft potential
        if (pj_artificial.isSingle()) _pi.pot_tot -= dr2_eps>r_out2? 0.0: (gmor*kpot  - gmor_max);   
        // member, mass is zero, use backup mass
        else if (pj_artificial.isMember()) _pi.pot_tot -= dr2_eps>r_out2? 0.0: (G*pj_artificial.getMassBackup()*drinv*kpot  - gmor_max);   
        // (orbitial) artificial, should be excluded in potential calculation, since it is inside neighbor, gmor_max cancel it to 0.0
        else _pi.pot_tot += gmor_max; 
#else
        // single/member, remove linear cutoff, obtain total potential
        if (pj_artificial.isSingle()) _pi.pot_tot -= (gmor - gmor_max);   
        // member, mass is zero, use backup mass
        else if (pj_artificial.isMember()) _pi.pot_tot -= (G*pj_artificial.getMassBackup()*drinv  - gmor_max);   
        // (orbitial) artificial, should be excluded in potential calculation, since it is inside neighbor, gmor_max cancel it to 0.0
        else _pi.pot_tot += gmor_max; 
#endif
        // correct to changeover soft acceleration
        _pi.acc -= (gmor3*k - gmor3_max)*dr;
    }

    //! correct force and potential for changeover function change
    /*!
      @param[in,out] _pi: particle for correction
      @param[in] _pj: j particle to calculate correction
     */
    template <class Tpi>
    inline void calcAccChangeOverCorrection(Tpi& _pi,
                                            const Ptcl& _pj) {
        PS::F64 G = ForceSoft::grav_const;
        const PS::F64vec dr = _pi.pos - _pj.pos;
        const PS::F64 dr2 = dr * dr;
#ifdef HARD_DEBUG
        assert(dr2>0.0);
#endif
        const PS::F64 dr2_eps = dr2 + manager->eps_sq;
        const PS::F64 drinv = 1.0/sqrt(dr2_eps);
        const PS::F64 gmor = G*_pj.mass * drinv;
        const PS::F64 drinv2 = drinv * drinv;
        const PS::F64 gmor3 = gmor * drinv2;
        const PS::F64 dr_eps = drinv * dr2_eps;

        // old
        const PS::F64 kold = 1.0 - ChangeOver::calcAcc0WTwo(_pi.changeover, _pj.changeover, dr_eps);

        // new
        ChangeOver chinew, chjnew;
        chinew.setR(_pi.changeover.getRin()*_pi.changeover.r_scale_next, _pi.changeover.getRout()*_pi.changeover.r_scale_next);
        chjnew.setR(_pj.changeover.getRin()*_pj.changeover.r_scale_next, _pj.changeover.getRout()*_pj.changeover.r_scale_next);
        const PS::F64 knew = 1.0 - ChangeOver::calcAcc0WTwo(chinew, chjnew, dr_eps);

        // correct to changeover soft acceleration
        _pi.acc -= gmor3*(knew-kold)*dr;
    }

    //! correct force and potential for changeover function change
    /*!
      @param[in,out] _pi: particle for correction
      @param[in] _pj: j particle to calculate correction
     */
    template <class Tpi>
    inline void calcAccChangeOverCorrection(Tpi& _pi,
                                            const EPJSoft& _pj) {
        PS::F64 G = ForceSoft::grav_const;
        const PS::F64vec dr = _pi.pos - _pj.pos;
        const PS::F64 dr2 = dr * dr;
        const PS::F64 dr2_eps = dr2 + manager->eps_sq;
        const PS::F64 drinv = 1.0/sqrt(dr2_eps);
        const PS::F64 gmor = G*_pj.mass * drinv;
        const PS::F64 drinv2 = drinv * drinv;
        const PS::F64 gmor3 = gmor * drinv2;
        const PS::F64 dr_eps = drinv * dr2_eps;

        ChangeOver chjold;
        chjold.setR(_pj.r_in, _pj.r_out);
        // old
        const PS::F64 kold = 1.0 - ChangeOver::calcAcc0WTwo(_pi.changeover, chjold, dr_eps);

        // new
        ChangeOver chinew, chjnew;
        chinew.setR(_pi.changeover.getRin()*_pi.changeover.r_scale_next, _pi.changeover.getRout()*_pi.changeover.r_scale_next);
        chjnew.setR(_pj.r_in*_pj.r_scale_next, _pj.r_out*_pj.r_scale_next);
        const PS::F64 knew = 1.0 - ChangeOver::calcAcc0WTwo(chinew, chjnew, dr_eps);

        // correct to changeover soft acceleration
        _pi.acc -= gmor3*(knew-kold)*dr;
    }

#ifdef KDKDK_4TH
    template <class Tpi>
    inline void calcAcorrShortWithLinearCutoff(Tpi& _pi,
                                               const Ptcl& _pj) {
        PS::F64 G = ForceSoft::grav_const;
        const PS::F64 r_out = manager->changeover.getRout();
        const PS::F64 r_out2 = r_out * r_out;

        const PS::F64vec dr = _pi.pos - _pj.pos;
        const PS::F64vec da = _pi.acc - _pi.acc;
        const PS::F64 dr2 = dr * dr;
        const PS::F64 dr2_eps = dr2 + manager->eps_sq;
        const PS::F64 drda = dr*da;
        const PS::F64 drinv = 1.0/sqrt(dr2_eps);
        const PS::F64 drdadrinv = drda*drinv;
        const PS::F64 gmor = G*_pj.mass * drinv;
        const PS::F64 drinv2 = drinv * drinv;
        const PS::F64 gmor3 = gmor * drinv2;
        const PS::F64 dr_eps = drinv * dr2_eps;

        const PS::F64 k = 1.0 - ChangeOver::calcAcc0WTwo(_pi.changeover, _pj.changeover, dr_eps);
        const PS::F64 kdot = - ChangeOver::calcAcc1WTwo(_pi.changeover, _pj.changeover, dr_eps, drdadrinv);

        const PS::F64 dr2_max = (dr2_eps > r_out2) ? dr2_eps : r_out2;
        const PS::F64 drinv_max = 1.0/sqrt(dr2_max);
        const PS::F64 gmor_max = G*_pj.mass * drinv_max;
        const PS::F64 drinv2_max = drinv_max*drinv_max;
        const PS::F64 gmor3_max = gmor_max * drinv2_max;

        const PS::F64 alpha = drda*drinv2;
        const PS::F64 alpha_max = drda * drinv2_max;
        const PS::F64vec acorr_k = gmor3 * (k*da - (3.0*k*alpha - kdot) * dr);
        const PS::F64vec acorr_max = gmor3_max * (da - 3.0*alpha_max * dr);

        _pi.acorr -= 2.0 * (acorr_k - acorr_max);
        //acci + dt_kick * dt_kick * acorri /48; 
    }

    template <class Tpi>
    inline void calcAcorrShortWithLinearCutoff(Tpi& _pi,
                                               const EPJSoft& _pj) {
        PS::F64 G = ForceSoft::grav_const;
        const PS::F64 r_out = manager->changeover.getRout();
        const PS::F64 r_out2 = r_out * r_out;

        const PS::F64vec dr = _pi.pos - _pj.pos;
        const PS::F64vec da = _pi.acc - _pi.acc;
        const PS::F64 dr2 = dr * dr;
        const PS::F64 dr2_eps = dr2 + manager->eps_sq;
        const PS::F64 drda = dr*da;
        const PS::F64 drinv = 1.0/sqrt(dr2_eps);
        const PS::F64 drdadrinv = drda*drinv;
        const PS::F64 gmor = G*_pj.mass * drinv;
        const PS::F64 drinv2 = drinv * drinv;
        const PS::F64 gmor3 = gmor * drinv2;
        const PS::F64 dr_eps = drinv * dr2_eps;
        ChangeOver chj;
        chj.setR(_pj.r_in, _pj.r_out);
        const PS::F64 k = 1.0 - ChangeOver::calcAcc0WTwo(_pi.changeover, chj, dr_eps);
        const PS::F64 kdot = - ChangeOver::calcAcc1WTwo(_pi.changeover, chj, dr_eps, drdadrinv);

        const PS::F64 dr2_max = (dr2_eps > r_out2) ? dr2_eps : r_out2;
        const PS::F64 drinv_max = 1.0/sqrt(dr2_max);
        const PS::F64 gmor_max = G*_pj.mass * drinv_max;
        const PS::F64 drinv2_max = drinv_max*drinv_max;
        const PS::F64 gmor3_max = gmor_max * drinv2_max;

        const PS::F64 alpha = drda*drinv2;
        const PS::F64 alpha_max = drda * drinv2_max;
        const PS::F64vec acorr_k = gmor3 * (k*da - (3.0*k*alpha - kdot) * dr);
        const PS::F64vec acorr_max = gmor3_max * (da - 3.0*alpha_max * dr);

        _pi.acorr -= 2.0 * (acorr_k - acorr_max);
        //acci + dt_kick * dt_kick * acorri /48; 
    }
#endif

    //! soft force correction use tree neighbor search for one particle
    /*
      @param[in,out] _psoft: particle in global system need to be corrected for acc and pot
      @param[in] _tree: tree for force
      @param[in] _acorr_flag: flag to do acorr for KDKDK_4TH case
     */
    template <class Tpsoft, class Ttree, class Tepj>
    void correctForceWithCutoffTreeNeighborOneParticleImp(Tpsoft& _psoft, 
                                                          Ttree& _tree,
                                                          const bool _acorr_flag=false) {
        Tepj * ptcl_nb = NULL;
        PS::S32 n_ngb = _tree.getNeighborListOneParticle(_psoft, ptcl_nb);
#ifdef HARD_DEBUG
        assert(n_ngb >= 1);
#endif
        // self-potential correction 
        // no correction for orbital artificial particles because the potential are not used for any purpose
        // no correction for member particles because their mass is zero during the soft force calculation, the self-potential contribution is also zero.
        if (_psoft.group_data.artificial.isSingle()) _psoft.pot_tot += _psoft.mass/manager->r_out_base; 
        // for binary without artificial particles, correction is needed.
        if (_psoft.group_data.artificial.isMember() && _psoft.getParticleCMAddress()<0) _psoft.pot_tot += _psoft.mass/manager->r_out_base; 

        // loop neighbors
        for(PS::S32 k=0; k<n_ngb; k++){
            if (ptcl_nb[k].id == _psoft.id) continue;

#ifdef KDKDK_4TH
            if(_acorr_flag) 
                calcAcorrShortWithLinearCutoff(_psoft, ptcl_nb[k]);
            else
#endif
                calcAccPotShortWithLinearCutoff(_psoft, ptcl_nb[k]);
        }
    }

    //! soft force correction for artificial particles in one cluster
    /* 1. Correct cutoff for artificial particles
       2. The c.m. force is substracted from tidal tensor force
       @param[in,out] _sys: global particle system, acc is updated
       @param[in] _ptcl_local: particle in systme_hard
       @param[in] _adr_real_start: real particle start address in _ptcl_local
       @param[in] _adr_real_end:   real particle end (+1) address in _ptcl_local
       @param[in] _n_group:  number of groups in cluster
       @param[in] _adr_first_ptcl_arti_in_cluster: address of the first artificial particle in each groups
       @param[in] _acorr_flag: flag to do acorr for KDKDK_4TH case
     */
    template <class Tsys>
    void correctForceWithCutoffArtificialOneClusterImp(Tsys& _sys, 
                                                      const PtclH4* _ptcl_local,
                                                      const PS::S32 _adr_real_start,
                                                      const PS::S32 _adr_real_end,
                                                      const PS::S32 _n_group,
                                                      const PS::S32* adr_first_ptcl_arti_in_cluster_,
                                                      const bool _acorr_flag) {

        auto& ap_manager = manager->ap_manager;
        for (int j=0; j<_n_group; j++) {  // j: j_group
            PS::S32 j_start = adr_first_ptcl_arti_in_cluster_[j];
#ifdef ARTIFICIAL_PARTICLE_DEBUG
            if (j_start<0) assert(_n_group==1&&_adr_real_end-_adr_real_start==2);
#endif
            if (j_start<0) continue;
            auto* p_arti_j = &(_sys[j_start]);
            auto* pj = ap_manager.getTidalTensorParticles(p_arti_j);

            auto correctionLoop = [&](const PS::S32 k) {
                // loop orbital artificial particle
                for (int kj=0; kj<_n_group; kj++) { // group
                    PS::S32 kj_start = adr_first_ptcl_arti_in_cluster_[kj];
                    if (kj_start<0) continue;
                    auto* porb_kj = ap_manager.getOrbitalParticles(&_sys[kj_start]);

                    // particle arti orbital
                    for (int kk=0; kk<ap_manager.getOrbitalParticleN(); kk++) {
                        if(&porb_kj[kk]==&(pj[k])) continue; //avoid same particle
#ifdef KDKDK_4TH
                        if(_acorr_flag) 
                            calcAcorrShortWithLinearCutoff(pj[k], porb_kj[kk]);
                        else
#endif
                            calcAccPotShortWithLinearCutoff(pj[k], porb_kj[kk]);
                    }
                }

                // loop real particle
                for (int kj=_adr_real_start; kj<_adr_real_end; kj++) {
#ifdef KDKDK_4TH
                    if(_acorr_flag) {
                        PS::S64 adr_kj = _ptcl_local[kj].adr_org;
                        calcAcorrShortWithLinearCutoff(pj[k], _sys[adr_kj]);
                    }
                    else
#endif
                        calcAccPotShortWithLinearCutoff(pj[k], _ptcl_local[kj]);
                }
            };

            // loop all tidal tensor particles 
            for (int k=0; k<ap_manager.getTidalTensorParticleN(); k++) {  
                // k: k_ptcl_arti
                correctionLoop(k);
            }

            // for c.m. particle
            pj = ap_manager.getCMParticles(p_arti_j);
            correctionLoop(0);
            
            ap_manager.correctArtficialParticleForce(p_arti_j);
        }
    }

    //! soft force correction completely use tree neighbor search
    /* @param[in,out] _sys: global particle system, acc is updated
       @param[in] _tree: tree for force
       @param[in] _ptcl_local: particle in systme_hard, only used to get adr_org
       @param[in] _n_ptcl: total number of particles in all clusters
       @param[in] _adr_ptcl_artificial_start: start address of artificial particle in _sys
       @param[in] _acorr_flag: flag to do acorr for KDKDK_4TH case
    */
    template <class Tsys, class Tpsoft, class Ttree, class Tepj>
    void correctForceWithCutoffTreeNeighborImp(Tsys& _sys, 
                                               Ttree& _tree, 
                                               const PtclH4* _ptcl_local,
                                               const PS::S32 _n_ptcl,
                                               const PS::S32 _adr_ptcl_artificial_start,
                                               const bool _acorr_flag=false) { 
        // for real particle
#pragma omp parallel for schedule(dynamic)
        for (int i=0; i<_n_ptcl; i++) {
            PS::S64 adr = _ptcl_local[i].adr_org;
            if(adr>=0) correctForceWithCutoffTreeNeighborOneParticleImp<Tpsoft, Ttree, Tepj>(_sys[adr], _tree, _acorr_flag);
        }

        // for artificial particle
        const PS::S32 n_tot = _sys.getNumberOfParticleLocal();
#pragma omp parallel for schedule(dynamic)
        for (int i=_adr_ptcl_artificial_start; i<n_tot; i++) 
            correctForceWithCutoffTreeNeighborOneParticleImp<Tpsoft, Ttree, Tepj>(_sys[i], _tree, _acorr_flag);

        auto& ap_manager = manager->ap_manager;
        const PS::S32 n_artificial_per_group = ap_manager.getArtificialParticleN();
#ifdef HARD_DEBUG
        assert((n_tot-_adr_ptcl_artificial_start)%n_artificial_per_group==0);
#endif
#pragma omp parallel for schedule(dynamic)
        for (int i=_adr_ptcl_artificial_start; i<n_tot; i+=n_artificial_per_group)
            ap_manager.correctArtficialParticleForce(&(_sys[i]));

    }

public:

    SystemHard(){
        manager = NULL;
        hard_int_ = NULL;
        n_hard_int_max_ = 0;
        n_hard_int_use_ = 0;

#ifdef PROFILE
        ARC_substep_sum = 0;
        ARC_tsyn_step_sum =0;
        ARC_n_groups = 0;
        ARC_n_groups_iso = 0;
        H4_step_sum = 0;
#endif
#ifdef HARD_COUNT_NO_NEIGHBOR
        n_neighbor_zero = 0;
#endif
#ifdef HARD_CHECK_ENERGY
        energy.clear();
#endif
        //        PS::S32 n_threads = PS::Comm::getNumberOfThread();
    }

    //! allocate memorgy for HardIntegrator
    /*! For record interrupt clusters. A array of HardIntegrator is allocated
      @param[in] _n_hard_int: number of HardIntegrator
     */
    void allocateHardIntegrator(const PS::S32 _n_hard_int) {
        assert(_n_hard_int>0);
        assert(n_hard_int_use_==0);
        if (hard_int_!=NULL) delete [] hard_int_;
        const PS::S32 num_thread = PS::Comm::getNumberOfThread();
        n_hard_int_max_ = _n_hard_int + num_thread;
        hard_int_ = new HardIntegrator[n_hard_int_max_];
        n_hard_int_use_ = 0;
    }

    ~SystemHard() {
        if (hard_int_!=NULL) delete [] hard_int_;
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
            if(p.mass==0&&p.group_data.artificial.isUnused()) {
                std::cerr<<"Error: unused particle is selected! i="<<i<<"; med[i].adr_sys="<<med[i].adr_sys_<<std::endl;
                abort();
            }
#endif
        }

        for(PS::S32 i=0; i<ptcl_recv.size(); i++){
            const Tptcl & p = ptcl_recv[i];
            ptcl_hard_.push_back(PtclHard(p, p.id_cluster, -(i+1)));
#ifdef HARD_DEBUG
            if(p.mass==0&&p.group_data.artificial.isUnused()) {
                std::cerr<<"Error: receive usused particle! i="<<i<<std::endl;
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

    PS::S32 getNumberOfClusters() const{
        return n_ptcl_in_cluster_.size();
    }

    PS::S32 getNumberOfInterruptClusters() const {
        return interrupt_list_.size();
    }

    HardIntegrator* getInterruptHardIntegrator(const std::size_t i) {
        return interrupt_list_[i];
    }

    PS::S32* getClusterNumberOfMemberList(const std::size_t i=0) const{
        return n_ptcl_in_cluster_.getPointer(i);
    }

    PS::S32* getClusterNumberOfMemberListOffset(const std::size_t i=0) const{
        return n_ptcl_in_cluster_disp_.getPointer(i);
    }

    PS::S32* getGroupNumberOfMemberList(const std::size_t i=0) const{
        return n_group_in_cluster_.getPointer(i);
    }

    PS::S32* getGroupNumberOfMemberListOffset(const std::size_t i=0) const{
        return n_group_in_cluster_offset_.getPointer(i);
    }

    PS::S32* getAdrPtclArtFirstList(const std::size_t i=0) const{
        return adr_first_ptcl_arti_in_cluster_.getPointer(i);
    }

    PS::S32 getNClusterChangeOverUpdate() const {
        return i_cluster_changeover_update_.size();
    }

    void setTimeOrigin(const PS::F64 _time_origin){
        time_origin_ = _time_origin;
    }

    void updateTimeWriteBack() {
        time_write_back_ = time_origin_;
    }

    PS::F64 getTimeOrigin() const {
        return time_origin_;
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
    void setPtclForOneClusterOMP(const Tsys & sys, 
                                 const PS::ReallocatableArray<PS::S32> & adr_array){
        // for one cluster
        const PS::S32 n = adr_array.size();
        //ptcl_hard_.resizeNoInitialize(n);
        //n_ptcl_in_cluster_.resizeNoInitialize(n);
#pragma omp parallel for 
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
    void driveForOneClusterOMP(const PS::F64 _dt) {
        const PS::S32 n = ptcl_hard_.size();
#pragma omp parallel for
        for(PS::S32 i=0; i<n; i++){
            PS::F64vec dr = ptcl_hard_[i].vel * _dt;
            ptcl_hard_[i].pos += dr;
#ifdef HARD_DEBUG
            // to avoid issue in cluster search with velocity
            auto& pcm = ptcl_hard_[i].group_data.cm;
            assert(pcm.mass==0.0);
            assert(pcm.vel.x==0.0);
            assert(pcm.vel.y==0.0);
            assert(pcm.vel.z==0.0);
#endif
            ptcl_hard_[i].Ptcl::calcRSearch(_dt);
            /*
              DriveKeplerRestricted(mass_sun_, 
              pos_sun_, ptcl_hard_[i].pos, 
              vel_sun_, ptcl_hard_[i].vel, dt); 
            */
        }
        
        time_origin_ += _dt;
    }

    template<class Tsys>
    void writeBackPtclForOneCluster(Tsys & sys, 
                                    PS::ReallocatableArray<PS::S32> & _remove_list){
        const PS::S32 n = ptcl_hard_.size();
        //PS::ReallocatableArray<PS::S32> removelist(n);
        for(PS::S32 i=0; i<n; i++){
            //PS::S32 adr = adr_array[i];
            PS::S32 adr = ptcl_hard_[i].adr_org;
#ifdef HARD_DEBUG
            assert(sys[adr].id == ptcl_hard_[i].id);
#endif
            sys[adr].DataCopy(ptcl_hard_[i]);
            if(sys[adr].mass==0.0&&sys[adr].group_data.artificial.isUnused()) {
                _remove_list.push_back(adr);
            }
        }
        updateTimeWriteBack();
    }

    template<class Tsys>
    void writeBackPtclForOneClusterOMP(Tsys & sys){
        const PS::S32 n = ptcl_hard_.size();
#pragma omp parallel for 
        for(PS::S32 i=0; i<n; i++){
            PS::S32 adr = ptcl_hard_[i].adr_org;
            //PS::S32 adr = adr_array[i];
#ifdef HARD_DEBUG
            assert(sys[adr].id == ptcl_hard_[i].id);
#endif
            sys[adr].DataCopy(ptcl_hard_[i]);
        }
        updateTimeWriteBack();
    }

//    template<class Tsys>
//    void writeBackPtclLocalOnlyOMP(Tsys & sys) {
//        const PS::S32 n = ptcl_hard_.size();
//#pragma omp parallel for schedule(dynamic)
//        for(PS::S32 i=0; i<n; i++){
//            PS::S32 adr = ptcl_hard_[i].adr_org;
//            //PS::S32 adr = adr_array[i];
//#ifdef HARD_DEBUG
//            if(adr>=0) assert(sys[adr].id == ptcl_hard_[i].id);
//#endif
//            if(adr>=0) sys[adr].DataCopy(ptcl_hard_[i]);
//        }
//    }
// for one cluster
//////////////////


//////////////////
// for isolated multi cluster only
    template<class Tsys>
    void setPtclForIsolatedMultiClusterOMP(const Tsys & sys,
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
#pragma omp parallel for 
        for(PS::S32 i=0; i<n_ptcl; i++){
            PS::S32 adr = _adr_array[i];
            ptcl_hard_[i].DataCopy(sys[adr]);
            ptcl_hard_[i].adr_org = adr;
            //  ptcl_hard_[i].n_ngb= sys[adr].n_ngb;
        }

        interrupt_list_.resizeNoInitialize(0);
    }

    template<class Tsys>
    void writeBackPtclForMultiCluster(Tsys & _sys, 
                                      PS::ReallocatableArray<PS::S32> & _remove_list){
        writeBackPtclForOneCluster(_sys, _remove_list);
    }

    template<class Tsys>
    void writeBackPtclForMultiClusterOMP(Tsys & _sys) { 
        writeBackPtclForOneClusterOMP(_sys);
    }

    //! Hard integration for clusters
    /*! Integrate (drift) all clusters with OpenMP
      If interrupt integration exist, record in the interrupt_list_;
       @param[in] _dt: integration ending time (initial time is fixed to 0)
       @param[in] _ptcl_soft: global particle array which contains the artificial particles for constructing tidal tensors.
       \return interrupt cluster number
     */
    template<class Tpsoft>
    int driveForMultiClusterOMP(const PS::F64 dt, Tpsoft* _ptcl_soft){
        assert(n_hard_int_use_==0);

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
        const PS::S32 num_thread = PS::Comm::getNumberOfThread();
        assert(n_hard_int_max_>num_thread);

#ifndef ONLY_SOFT
        HardIntegrator* hard_int_thread[num_thread];
        // set new hard_int front pointer 
        HardIntegrator* hard_int_front_ptr = &hard_int_[num_thread];
        for (PS::S32 i=0; i<num_thread; i++) {
            hard_int_thread[i] = &hard_int_[i];
        }
#endif

#ifdef OMP_PROFILE        
        PS::F64 time_thread[num_thread];
        PS::S64 num_cluster[num_thread];
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
            // Hermite + AR integration
            const PS::S32 n_group = n_group_in_cluster_[i];
            Tpsoft* ptcl_artificial_ptr=NULL;
            PS::S32* n_member_in_group_ptr=NULL;
            if(n_group>0) {
                PS::S32 ptcl_arti_first_index = adr_first_ptcl_arti_in_cluster_[n_group_in_cluster_offset_[i]];
                if (ptcl_arti_first_index>=0) ptcl_artificial_ptr = &(_ptcl_soft[ptcl_arti_first_index]);
#ifdef PROFILE
                else ARC_n_groups_iso += 1;
#endif
                n_member_in_group_ptr = &(n_member_in_group_[n_group_in_cluster_offset_[i]]);
            }
#ifdef OMP_PROFILE
            num_cluster[ith] += n_ptcl;
#endif
#ifdef PROFILE
            ARC_n_groups  += n_group;
#endif

#ifdef HARD_DUMP
            assert(ith<hard_dump.size);
            hard_dump[ith].backup(ptcl_hard_.getPointer(adr_head), n_ptcl, ptcl_artificial_ptr, n_group, n_member_in_group_ptr, dt, manager->ap_manager.getArtificialParticleN());
#endif

#ifdef HARD_DEBUG_PROFILE
            PS::F64 tstart = PS::GetWtime();
#endif

            // if interrupt exist, escape initial
            hard_int_thread[ith]->initial(ptcl_hard_.getPointer(adr_head), n_ptcl, ptcl_artificial_ptr, n_group, n_member_in_group_ptr, manager, time_origin_);

            auto& interrupt_binary = hard_int_thread[ith]->integrateToTime(dt);

            if (interrupt_binary.status!=AR::InterruptStatus::none) {
                #pragma omp atomic capture
                hard_int_thread[ith] = hard_int_front_ptr++;

                assert(hard_int_thread[ith]!=&hard_int_[n_hard_int_max_]);
            }
            else {
                hard_int_thread[ith]->driftClusterCMRecordGroupCMDataAndWriteBack(dt);

#ifdef PROFILE
                ARC_substep_sum    += hard_int_thread[ith]->ARC_substep_sum;
                ARC_tsyn_step_sum  += hard_int_thread[ith]->ARC_tsyn_step_sum;
                H4_step_sum        += hard_int_thread[ith]->H4_step_sum;
#endif
#ifdef HARD_COUNT_NO_NEIGHBOR
                n_neighbor_zero    += hard_int_thread[ith]->n_neighbor_zero;
#endif
#ifdef HARD_CHECK_ENERGY
                energy += hard_int_thread[ith]->energy;
#endif
                
                hard_int_thread[ith]->clear();
            }

#ifdef OMP_PROFILE
            time_thread[ith] += PS::GetWtime();
#endif

#ifdef HARD_DEBUG_PROFILE
            PS::F64 tend = PS::GetWtime();
            std::cerr<<"HT: "<<i<<" "<<ith<<" "<<n_cluster<<" "<<n_ptcl<<" "<<tend-tstart<<std::endl;
#endif

#else
            // Only soft drift
            auto* pi = ptcl_hard_.getPointer(adr_head);
            for (PS::S32 j=0; j<n_ptcl; j++) {
                PS::F64vec dr = pi[j].vel * dt;
                pi[j].pos += dr;
#ifdef CLUSTER_VELOCITY
                auto& pij_cm = pi[j].group_data.cm;
                pij_cm.mass = pij_cm.vel.x = pij_cm.vel.y = pij_cm.vel.z = 0.0;
#endif
                pi[j].calcRSearch(dt);
            }
#endif

        }

        // regist interrupted hard integrator
        assert(interrupt_list_.size()==0);
        for (auto iptr = hard_int_; iptr<hard_int_front_ptr; iptr++) 
            if (iptr->is_initialized) {
                assert(iptr->interrupt_binary.status!=AR::InterruptStatus::none);
#ifdef HARD_INTERRUPT_PRINT
                iptr->printInterruptBinaryInfo(std::cerr);
#endif
                interrupt_list_.push_back(iptr);
            }


        // advance time_origin if all clusters finished
        PS::S32 n_interrupt = interrupt_list_.size();
        if (n_interrupt==0) time_origin_ += dt;
        else interrupt_dt_ = dt;

        return n_interrupt;
    }

    //! Finish interrupt integration
    /*! Finish interrupted integrations, if new interruption appear, record in the interrupt_list and this function need to be called again after modification of interrupt clusters
      If no new interrupt cluster appear, update time_origin_ with drift time.
      @param [in] 
     */
    PS::S32 finishIntegrateInterruptClustersOMP() {
        PS::S32 n_interrupt = interrupt_list_.size();
#pragma omp parallel for schedule(dynamic)
        for (PS::S32 i=0; i<n_interrupt; i++) {
            auto hard_int_ptr = interrupt_list_[i];
            auto& interrupt_binary = hard_int_ptr->integrateToTime(interrupt_dt_);

            if (interrupt_binary.status==AR::InterruptStatus::none) {
                hard_int_ptr->driftClusterCMRecordGroupCMDataAndWriteBack(interrupt_dt_);

#ifdef PROFILE
                ARC_substep_sum    += hard_int_ptr->ARC_substep_sum;
                ARC_tsyn_step_sum  += hard_int_ptr->ARC_tsyn_step_sum;
                H4_step_sum        += hard_int_ptr->H4_step_sum;
#endif
#ifdef HARD_CHECK_ENERGY
                energy += hard_int_ptr->energy;
#endif
                
                hard_int_ptr->clear();
            }
        }
        
        // record new interrupt list
        PS::S32 i_front = 0;
        PS::S32 i_end = n_interrupt;
        while (i_front<i_end) {
            auto hard_int_front_ptr = interrupt_list_[i_front];
            if (hard_int_front_ptr->interrupt_binary.status!=AR::InterruptStatus::none) {
                assert(hard_int_front_ptr->is_initialized);
#ifdef HARD_INTERRUPT_PRINT
                hard_int_front_ptr->printInterruptBinaryInfo(std::cerr);
#endif
                i_front++;
            }
            else {
                interrupt_list_[i_front] = interrupt_list_[--i_end];
                interrupt_list_.decreaseSize(1);
            }
        }

        n_interrupt =  interrupt_list_.size();

        // advance time_origin if all clusters finished
        if (n_interrupt==0) time_origin_ += interrupt_dt_;

        return n_interrupt;
    }

    struct GroupIndexInfo{
        PS::S32 cluster_index, group_index, n_members, n_member_offset, isolated_case;
        GroupIndexInfo() {}
        GroupIndexInfo(PS::S32 _cluster, PS::S32 _group, PS::S32 _member, PS::S32 _offset, PS::S32 _isolated_case):
            cluster_index(_cluster), group_index(_group), n_members(_member), n_member_offset(_offset), isolated_case(_isolated_case) {}
    };

    //! generate artificial particles,
    /*  
        @param[in]     _i_cluster: cluster index
        @param[in,out] _ptcl_in_cluster: particle data
        @param[in]     _n_ptcl: total number of particle in _ptcl_in_cluster.
        @param[out]    _ptcl_artificial: artificial particles that will be added
        @parma[out]    _binary_table: binary information table 
        @param[out]    _n_groups: number of groups in current cluster
        @param[out]    _n_members_in_groups: number of members in each group, (cluster_index, group_index, n_members)
        @param[out]    _changeover_update_list: cluster index list for particles with changeover updates
        @param[in,out] _groups: searchGroupCandidate class, which contain 1-D group member index array, will be reordered by the minimum distance chain for each group
        @param[in]     _dt_tree: tree time step for calculating r_search and set stablility checker period limit
     */
    template <class Tptcl>
    void findGroupsAndCreateArtificialParticlesOneCluster(const PS::S32 _i_cluster,
                                                          Tptcl* _ptcl_in_cluster,
                                                          const PS::S32 _n_ptcl,
                                                          PS::ReallocatableArray<Tptcl> & _ptcl_artificial,
                                                          PS::ReallocatableArray<COMM::BinaryTree<PtclH4,COMM::Binary>> & _binary_table,
                                                          PS::S32 &_n_groups,
                                                          PS::ReallocatableArray<GroupIndexInfo>& _n_member_in_group,
                                                          PS::ReallocatableArray<PS::S32>& _changeover_update_list,
                                                          SearchGroupCandidate<Tptcl>& _groups,
                                                          const PS::F64 _dt_tree) {

        PS::S32 group_ptcl_adr_list[_n_ptcl];
        PS::S32 group_ptcl_adr_offset=0;
        bool changeover_update_flag = false;
        _n_groups = 0;
        auto& ap_manager = manager->ap_manager;

        for (int i=0; i<_groups.getNumberOfGroups(); i++) {
            PS::ReallocatableArray<COMM::BinaryTree<Tptcl,COMM::Binary>> binary_tree;   // hierarch binary tree

            const PS::S32 n_members = _groups.getNumberOfGroupMembers(i);
            binary_tree.reserve(n_members);

#ifdef ARTIFICIAL_PARTICLE_DEBUG
            assert(n_members<ARRAY_ALLOW_LIMIT);
#endif        
            PS::S32* member_list = _groups.getMemberList(i);
            binary_tree.resizeNoInitialize(n_members-1);
            // build hierarch binary tree from the minimum distant neighbors

            COMM::BinaryTree<Tptcl,COMM::Binary>::generateBinaryTree(binary_tree.getPointer(), member_list, n_members, _ptcl_in_cluster, ap_manager.gravitational_constant);

            struct {PS::F64 mean_mass_inv, rin, rout, dt_tree; } changeover_rsearch_pars = {Tptcl::mean_mass_inv, manager->r_in_base, manager->r_out_base, _dt_tree};
            // get new changeover and rsearch for c.m.
            binary_tree.back().processTreeIter(changeover_rsearch_pars, (PS::S64)-1, (PS::S64)-1, calcBinaryIDChangeOverAndRSearchIter);
            
            // stability check and break groups
            Stability<Tptcl> stable_checker;
            
            // find closed Kepler hierarchical systems
            stable_checker.stable_binary_tree.reserve(n_members);
            stable_checker.findClosedTree(binary_tree.back());

            // save group information to binary table.
            for (int k=0; k<stable_checker.stable_binary_tree.size(); k++) {
                auto& closed_binary_tree_k = *stable_checker.stable_binary_tree[k];
                const PS::S32 n_members = closed_binary_tree_k.getMemberN();
                PS::S32 start_index_binary_table = _binary_table.size();
                _binary_table.increaseSize(n_members);
                closed_binary_tree_k.getherBinaryTreeIter(_binary_table.getPointer(start_index_binary_table));
            }


            // be careful, here t_crit should be >= hard slowdown_timescale_max to avoid using slowdown for wide binaries
            stable_checker.t_crit = manager->ar_manager.slowdown_timescale_max;
            stable_checker.findStableTree(binary_tree.back());

            // in isolated binary case, if the slowdown period can be larger than tree step * N_step_per_orbit, it is fine without tidal tensor
            if (_n_ptcl==2&&stable_checker.stable_binary_tree.size()==1) {
                // the criterion is: period < 0.5*pi^2 * slowdown_pert_ratio_pef * R_search/G(m1+m2) / (dt_soft*n_step_per_orbit)
                auto& bin = *stable_checker.stable_binary_tree[i];
                const PS::S32 n_members = bin.getMemberN();
                PS::F64 pot_ch_inv = bin.r_search*bin.r_search*bin.r_search/(ap_manager.gravitational_constant*bin.mass);
                PS::F64 dt_nstep = _dt_tree*manager->n_step_per_orbit;
                if (bin.period < 4.93480220054 * manager->ar_manager.slowdown_pert_ratio_ref * pot_ch_inv / dt_nstep) {
                    // Set member particle type, backup mass, collect member particle index to group_ptcl_adr_list
                    //use _ptcl_in_cluster as the first particle address as reference to calculate the particle index.
                    struct { Tptcl* adr_ref; PS::S32* group_list; PS::S32 n; PS::S64 pcm_id; ChangeOver* changeover_cm; PS::F64 rsearch_cm; bool changeover_update_flag;}
                    group_index_pars = { _ptcl_in_cluster,  &group_ptcl_adr_list[group_ptcl_adr_offset], 0, -1, &bin.changeover, bin.r_search, false};
                    bin.processLeafIter(group_index_pars, collectGroupMemberAdrAndSetMemberParametersIter);
                    if (group_index_pars.changeover_update_flag) changeover_update_flag = true;
#ifdef ARTIFICIAL_PARTICLE_DEBUG
                    assert(group_index_pars.n==n_members);
#ifdef ARTIFICIAL_PARTICLE_DEBUG_PRINT
                    std::cerr<<"Isolated binary case: "
                             <<" period: "<<bin.period
                             <<" semi: "<<bin.semi
                             <<" ecc: "<<bin.ecc
                             <<" r_serach: "<<bin.r_search
                             <<" tree_step: "<<_dt_tree
                             <<" n_step_per_orbit: "<<manager->n_step_per_orbit
                             <<std::endl;
#endif
#endif                
                    _n_member_in_group.push_back(GroupIndexInfo(_i_cluster, _n_groups, n_members, group_ptcl_adr_offset, 1));
                    _n_groups++;
                    group_ptcl_adr_offset += n_members;
                    continue;
                }
            }

            // save cluster and group index in artificial particle data as identification for later collection to global sys
            PS::F64 index_group[2];
            index_group[0] = (PS::F64)(_i_cluster+1);

            for (int i=0; i<stable_checker.stable_binary_tree.size(); i++) {
                index_group[1] = (PS::F64)(_n_groups+1);

                auto& binary_stable_i = *stable_checker.stable_binary_tree[i];
                const PS::S32 n_members = binary_stable_i.getMemberN();

                // Make sure the _ptcl_new will not make new array due to none enough capacity during the following loop, 
                const int n_ptcl_artificial = _ptcl_artificial.size();
                _ptcl_artificial.increaseSize(ap_manager.getArtificialParticleN());
                Tptcl* ptcl_artificial_i = &_ptcl_artificial[n_ptcl_artificial];

                // Set member particle type to member, set mass to zero, backup mass and collect member particle index to group_ptcl_adr_list
                // set rsearch and changeover for members
                //use _ptcl_in_cluster as the first particle address as reference to calculate the particle index.
                struct { Tptcl* adr_ref; PS::S32* group_list; PS::S32 n; PS::S64 pcm_id; ChangeOver* changeover_cm; PS::F64 rsearch_cm; bool changeover_update_flag;}
                group_index_pars = { _ptcl_in_cluster,  &group_ptcl_adr_list[group_ptcl_adr_offset], 0, n_ptcl_artificial, &binary_stable_i.changeover, binary_stable_i.r_search, false};
                binary_stable_i.processLeafIter(group_index_pars, collectGroupMemberAdrAndSetMemberParametersIter);
#ifdef ARTIFICIAL_PARTICLE_DEBUG
                assert(group_index_pars.n==n_members);
#endif                
                if (group_index_pars.changeover_update_flag) changeover_update_flag = true;

                // generate artificial particles
                ap_manager.createArtificialParticles(ptcl_artificial_i, binary_stable_i, index_group, 2);

                // set rsearch and changeover for c.m. particle
                auto* pcm = ap_manager.getCMParticles(ptcl_artificial_i);
                pcm->r_search   = binary_stable_i.r_search;
                pcm->r_search  += binary_stable_i.semi*(1+binary_stable_i.ecc);  // depend on the mass ratio, the upper limit distance to c.m. from all members and artificial particles is apo-center distance
                pcm->changeover = binary_stable_i.changeover;
#ifdef ARTIFICIAL_PARTICLE_DEBUG
                assert(pcm->r_search > pcm->changeover.getRout());
#endif                

                // set rsearch and changeover for tidal tensor particles
                Tptcl* ptcl_tt = ap_manager.getTidalTensorParticles(ptcl_artificial_i);
                // use c.m. r_search and changeover
                for (int j=0; j<ap_manager.getTidalTensorParticleN(); j++) {
                    Tptcl& pj = ptcl_tt[j];
                    pj.r_search   = binary_stable_i.r_search;
                    pj.changeover = binary_stable_i.changeover;

#ifdef ARTIFICIAL_PARTICLE_DEBUG
                    //check rsearch consistence:
                    PS::F64 rsearch_bin = pcm->r_search;
                    PS::F64vec dp = pj.pos - binary_stable_i.pos;
                    PS::F64 dr = dp*dp;
                    assert(dr<=rsearch_bin*rsearch_bin);
#endif
                }

                // set rsearch and changeover for orbitial sample particles
                Tptcl* ptcl_orbit=ap_manager.getOrbitalParticles(ptcl_artificial_i);
                PS::S32 n_orbit = ap_manager.getOrbitalParticleN();
                for (int j=0; j<n_orbit; j++) {
                    // use member changeover, if new changeover is different, record the scale ratio 
                    Tptcl& pj = ptcl_orbit[j];
                    pj.changeover =  binary_stable_i.getMember(j%2)->changeover;

                    PS::F64 pj_r_in  = pj.changeover.getRin();
                    PS::F64 bin_r_in = binary_stable_i.changeover.getRin();
                    if (abs( pj_r_in - bin_r_in)>1e-10) {
                        pj.changeover.r_scale_next = bin_r_in/pj_r_in;
                        pj.r_search = std::max(pj.r_search, binary_stable_i.r_search);
                        // not necessary true since the member changeover may inherient from other binaries which can be larger than the new one here.
                        //assert(_bin.changeover.getRin()>=pj->changeover.getRin());

                        // in case the new changeover is small, r_search may be smaller than the previous r_out, which should be avoided
                        if (pj.r_search < pj.changeover.getRout()) {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
                            assert( pj.r_search > pj.changeover.getRout()*pj.changeover.r_scale_next);
#endif 
                            pj.r_search = pj.changeover.getRout(); 
                        }
                    }
                    else pj.r_search = binary_stable_i.r_search;

#ifdef ARTIFICIAL_PARTICLE_DEBUG
                    //check rsearch consistence:
                    PS::F64 rsearch_bin = pcm->r_search;
                    PS::F64vec dp = pj.pos - binary_stable_i.pos;
                    PS::F64 dr = dp*dp;
                    assert(dr<=rsearch_bin*rsearch_bin);
#endif
                }

                _n_member_in_group.push_back(GroupIndexInfo(_i_cluster, _n_groups, n_members, group_ptcl_adr_offset, 0));
                group_ptcl_adr_offset += n_members;
                _n_groups++;
            }
        }

        assert(group_ptcl_adr_offset<=_n_ptcl);

        // Reorder the ptcl that group member come first
        PS::S32 ptcl_list_reorder[_n_ptcl];
        for (int i=0; i<_n_ptcl; i++) ptcl_list_reorder[i] = i;
 
        // shift single after group members
        PS::S32 i_single_front=group_ptcl_adr_offset;
        PS::S32 i_group = 0;
        while (i_group<group_ptcl_adr_offset) {
            // if single find inside group_ptcl_adr_offset, exchange single with group member out of the offset
            if(_ptcl_in_cluster[i_group].group_data.artificial.isSingle()) {
                while(_ptcl_in_cluster[i_single_front].group_data.artificial.isSingle()) {
                    i_single_front++;
                    assert(i_single_front<_n_ptcl);
                }
                // Swap index
                PS::S32 plist_tmp = ptcl_list_reorder[i_group];
                ptcl_list_reorder[i_group] = ptcl_list_reorder[i_single_front];
                ptcl_list_reorder[i_single_front] = plist_tmp;
                i_single_front++; // avoild same particle be replaced
            }
            i_group++;
        }

#ifdef ARTIFICIAL_PARTICLE_DEBUG
        // check whether the list is correct
        PS::S32 plist_new[group_ptcl_adr_offset];
        for (int i=0; i<group_ptcl_adr_offset; i++) plist_new[i] = group_ptcl_adr_list[i];
        std::sort(plist_new, plist_new+group_ptcl_adr_offset, [](const PS::S32 &a, const PS::S32 &b) {return a < b;});
        std::sort(ptcl_list_reorder, ptcl_list_reorder+group_ptcl_adr_offset, [](const PS::S32 &a, const PS::S32 &b) {return a < b;});
        for (int i=0; i<group_ptcl_adr_offset; i++) assert(ptcl_list_reorder[i]==plist_new[i]);
#endif        

        // overwrite the new ptcl list for group members by reorderd list
        for (int i=0; i<group_ptcl_adr_offset; i++) ptcl_list_reorder[i] = group_ptcl_adr_list[i];

        // templately copy ptcl data
        Tptcl ptcl_tmp[_n_ptcl];
        for (int i=0; i<_n_ptcl; i++) ptcl_tmp[i]=_ptcl_in_cluster[i];

        // reorder ptcl
        for (int i=0; i<_n_ptcl; i++) _ptcl_in_cluster[i]=ptcl_tmp[ptcl_list_reorder[i]];

        // changeover update
        if (changeover_update_flag) _changeover_update_list.push_back(_i_cluster);
    }

    //! Find groups and create aritfical particles to sys
    /* @param[in,out] _sys: global particle system
       @param[in]     _dt_tree: tree time step for calculating r_search
     */
    template<class Tsys, class Tptcl>
    void findGroupsAndCreateArtificialParticlesOMP(Tsys & _sys,
                                                   const PS::F64 _dt_tree) {
        const PS::S32 n_cluster = n_ptcl_in_cluster_.size();
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(n_cluster<ARRAY_ALLOW_LIMIT);
#endif        
        n_group_in_cluster_.resizeNoInitialize(n_cluster);
        n_group_member_remote_=0;

        const PS::S32 num_thread = PS::Comm::getNumberOfThread();
        PS::ReallocatableArray<PtclH4> ptcl_artificial_thread[num_thread];
        PS::ReallocatableArray<COMM::BinaryTree<PtclH4,COMM::Binary>> binary_table_thread[num_thread];
        PS::ReallocatableArray<GroupIndexInfo> n_member_in_group_thread[num_thread];
        PS::ReallocatableArray<PS::S32> i_cluster_changeover_update_threads[num_thread];
        for (PS::S32 i=0; i<num_thread; i++) {
            ptcl_artificial_thread[i].resizeNoInitialize(0);
            binary_table_thread[i].resizeNoInitialize(0);
            n_member_in_group_thread[i].resizeNoInitialize(0);
            i_cluster_changeover_update_threads[i].resizeNoInitialize(0);
        }
        auto& ap_manager = manager->ap_manager;

#pragma omp parallel for schedule(dynamic)
        for (PS::S32 i=0; i<n_cluster; i++){
            const PS::S32 ith = PS::Comm::getThreadNum();
            PtclH4* ptcl_in_cluster = ptcl_hard_.getPointer() + n_ptcl_in_cluster_disp_[i];
            const PS::S32 n_ptcl = n_ptcl_in_cluster_[i];
            // reset particle type to single
            for(PS::S32 j=0; j<n_ptcl; j++) {
                // ensure both hard local and global system have reset status, otherwise singles in global system may have wrong status
                ptcl_in_cluster[j].group_data.artificial.setParticleTypeToSingle();
                PS::S64 adr=ptcl_in_cluster[j].adr_org;
                if(adr>=0) _sys[adr].group_data = ptcl_in_cluster[j].group_data;
            }
            // search group_candidates
            SearchGroupCandidate<PtclH4> group_candidate;
            // merge group_candidates
            group_candidate.searchAndMerge(ptcl_in_cluster, n_ptcl);

            // find groups and generate artificial particles for cluster i
            findGroupsAndCreateArtificialParticlesOneCluster(i, ptcl_in_cluster, n_ptcl, ptcl_artificial_thread[ith], binary_table_thread[ith], n_group_in_cluster_[i], n_member_in_group_thread[ith], i_cluster_changeover_update_threads[ith], group_candidate, _dt_tree);
        }

        // gether binary table
        PS::S32 n_binary_table_offset_thread[num_thread+1];
        n_binary_table_offset_thread[0] = 0;
        for (PS::S32 i=0; i<num_thread; i++) 
            n_binary_table_offset_thread[i+1] = n_binary_table_offset_thread[i] + binary_table_thread[i].size();
        binary_table.resizeNoInitialize(n_binary_table_offset_thread[num_thread]);
#pragma omp parallel for 
        for (PS::S32 i=0; i<num_thread; i++) {
            for (PS::S32 k=n_binary_table_offset_thread[i]; k<n_binary_table_offset_thread[i+1]; k++) {
                binary_table[k] = binary_table_thread[i][k-n_binary_table_offset_thread[i]];
            }
        }

        // n_group_in_cluster_offset
        n_group_in_cluster_offset_.resizeNoInitialize(n_cluster+1);
        n_group_in_cluster_offset_[0] = 0;
        for (PS::S32 i=0; i<n_cluster; i++) {
            n_group_in_cluster_offset_[i+1] = n_group_in_cluster_offset_[i] + n_group_in_cluster_[i];
        }
        n_member_in_group_.resizeNoInitialize(n_group_in_cluster_offset_[n_cluster]);
        n_member_in_group_offset_.resizeNoInitialize(n_group_in_cluster_offset_[n_cluster]);

        // artificial particle first address array
        adr_first_ptcl_arti_in_cluster_.resizeNoInitialize(n_group_in_cluster_offset_[n_cluster]);
        for (PS::S32 i=0; i<adr_first_ptcl_arti_in_cluster_.size(); i++) adr_first_ptcl_arti_in_cluster_[i] = -1;

#ifdef ARTIFICIAL_PARTICLE_DEBUG
        for (PS::S32 i=0; i<n_member_in_group_.size(); i++) n_member_in_group_[i] = 0;
#endif

#pragma omp parallel for 
        for (PS::S32 i=0; i<num_thread; i++) {
            for (PS::S32 k=0; k<n_member_in_group_thread[i].size(); k++) {
                PS::S32 i_cluster = n_member_in_group_thread[i][k].cluster_index;
                PS::S32 j_group = n_member_in_group_thread[i][k].group_index;
                PS::S32 n_members = n_member_in_group_thread[i][k].n_members;
                PS::S32 n_member_offset = n_member_in_group_thread[i][k].n_member_offset;
                n_member_in_group_[n_group_in_cluster_offset_[i_cluster]+j_group] =  n_members;
                n_member_in_group_offset_[n_group_in_cluster_offset_[i_cluster]+j_group] =  n_member_offset;
                // update system soft 
                if (n_member_in_group_thread[i][k].isolated_case) {
                    for (PS::S32 j=0; j<n_members; j++) {
                        auto& p_loc = ptcl_hard_[n_ptcl_in_cluster_disp_[i_cluster]+n_member_offset+j];
#ifdef ARTIFICIAL_PARTICLE_DEBUG
                        assert(p_loc.group_data.artificial.isMember());
                        assert(p_loc.getParticleCMAddress()<0);
                        assert(adr_first_ptcl_arti_in_cluster_[n_group_in_cluster_offset_[i_cluster]+j_group]<0);
#endif
                        const PS::S64 p_glb_adr= p_loc.adr_org;
                        if(p_glb_adr>=0) {
                            auto& p_glb = _sys[p_glb_adr];
                            p_glb.group_data= p_loc.group_data;
                            p_glb.changeover= p_loc.changeover;
                            p_glb.r_search  = p_loc.r_search;
                        }
                        else {
                            // this is remoted member;
                            n_group_member_remote_++;
                        }
                    }
                }
            }
        }


#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(n_group_in_cluster_offset_[n_cluster]<ARRAY_ALLOW_LIMIT);
        for (PS::S32 i=0; i<n_member_in_group_.size(); i++) assert(n_member_in_group_[i] > 0);
#endif        

        // add artificial particle to particle system
        PS::S32 rank = PS::Comm::getRank();
        // Get the address offset for new artificial ptcl array in each thread in _sys
        PS::S64 sys_ptcl_artificial_thread_offset[num_thread+1];
        sys_ptcl_artificial_thread_offset[0] = _sys.getNumberOfParticleLocal();
        for(PS::S32 i=0; i<num_thread; i++) {
            sys_ptcl_artificial_thread_offset[i+1] = sys_ptcl_artificial_thread_offset[i] + ptcl_artificial_thread[i].size();
        }
        _sys.setNumberOfParticleLocal(sys_ptcl_artificial_thread_offset[num_thread]);
        
        const PS::S32 n_artificial_per_group = ap_manager.getArtificialParticleN();
#pragma omp parallel for        
        for(PS::S32 i=0; i<num_thread; i++) {
            // ptcl_artificial should be integer times of n_partificial_per_group 
            assert(ptcl_artificial_thread[i].size()%n_artificial_per_group==0);
            // Add particle to ptcl sys
            for (PS::S32 j=0; j<ptcl_artificial_thread[i].size(); j++) {
                PS::S32 adr = j+sys_ptcl_artificial_thread_offset[i];
                ptcl_artificial_thread[i][j].adr_org=adr;
                _sys[adr]=Tptcl(ptcl_artificial_thread[i][j],rank,adr);
            }
            // Update the status of group members to c.m. address in ptcl sys. Notice c.m. is at the end of an artificial particle group
            PS::S32 group_offset=0, j_group_recored=-1;
            for (PS::S32 j=0; j<ptcl_artificial_thread[i].size(); j+=n_artificial_per_group) {
                auto* pj = &ptcl_artificial_thread[i][j];
                auto* pcm = ap_manager.getCMParticles(pj);
                PS::S32 n_members = ap_manager.getMemberN(pj);
                PS::S32 i_cluster = ap_manager.getStoredData(pj,0,true)-1; 
                PS::S32 j_group = ap_manager.getStoredData(pj,1,true)-1;
#ifdef ARTIFICIAL_PARTICLE_DEBUG
                assert(n_member_in_group_[n_group_in_cluster_offset_[i_cluster]+j_group]==n_members);
#endif                
                // make sure group index increase one by one
                assert(j_group==j_group_recored+1);
                j_group_recored=j_group;

#ifdef ARTIFICIAL_PARTICLE_DEBUG
                // check whether ID is consistent.
                PS::S64 id_cm = ap_manager.getCMID(pj);
                bool id_match_flag=false;
#endif
                for (PS::S32 k=0; k<n_members; k++) {
                    PS::S32 kl = n_ptcl_in_cluster_disp_[i_cluster]+group_offset+k;
                    auto& p_loc = ptcl_hard_.getPointer()[kl];
                    // save c.m. address 
#ifdef ARTIFICIAL_PARTICLE_DEBUG
                    assert(group_offset == n_member_in_group_offset_[n_group_in_cluster_offset_[i_cluster]+j_group]);
                    assert(p_loc.getParticleCMAddress()==j+1);
                    assert(p_loc.group_data.artificial.isMember());
#endif                    
                    p_loc.setParticleCMAddress(pcm->adr_org);

#ifdef ARTIFICIAL_PARTICLE_DEBUG
                    if(p_loc.id == id_cm) id_match_flag=true;
#endif
                    // update global particle to be consistent
                    const PS::S64 p_glb_adr= p_loc.adr_org;
                    if(p_glb_adr>=0) {
                        auto& p_glb = _sys[p_glb_adr];
                        p_glb.mass      = p_loc.mass;
                        p_glb.group_data= p_loc.group_data;
                        p_glb.changeover= p_loc.changeover;
                        p_glb.r_search  = p_loc.r_search;
                    }
                    else {
                        // this is remoted member;
                        n_group_member_remote_++;
                    }

                }

#ifdef ARTIFICIAL_PARTICLE_DEBUG
                assert(id_match_flag);
#endif
                // shift cluster
                if(j_group==n_group_in_cluster_[i_cluster]-1) {
                    group_offset=0;
                    j_group_recored=-1;
                }
                else group_offset += n_members; // group offset in the ptcl list index of one cluster
                // j_group should be consistent with n_group[i_cluster];
                assert(j_group<=n_group_in_cluster_[i_cluster]);

                // save first address of artificial particle
                adr_first_ptcl_arti_in_cluster_[n_group_in_cluster_offset_[i_cluster]+j_group] = ptcl_artificial_thread[i][j].adr_org;
            }
        }
        
        // merge i_cluster_changeover
        i_cluster_changeover_update_.resizeNoInitialize(0);
        for(PS::S32 i=0; i<num_thread; i++) {
            for (PS::S32 j=0; j<i_cluster_changeover_update_threads[i].size();j++) 
                i_cluster_changeover_update_.push_back(i_cluster_changeover_update_threads[i][j]);
        }
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        if (i_cluster_changeover_update_.size()>0) 
            std::cerr<<"Changeover change cluster found: T="<<time_origin_<<" n_cluster_change="<<i_cluster_changeover_update_.size()<<std::endl;
#endif
        /*
        // sort data
        PS::S32 i_cluster_size = i_cluster_changeover_update_.size();
        if (i_cluster_size>0) {
            PS::S32* i_cluster_data = i_cluster_changeover_update_.getPointer();
            std::sort(i_cluster_data, i_cluster_data+i_cluster_size, [] (const PS::S32 &a, const PS::S32 &b) { return a<b; });
            // remove dup
            PS::S32* i_end = std::unique(i_cluster_data, i_cluster_data+i_cluster_size);
#ifdef ARTIFICIAL_PARTICLE_DEBUG
            assert(i_end-i_cluster_data>=0&&i_end-i_cluster_data<=i_cluster_size);
            std::cerr<<"Changeover change cluster found: T="<<time_origin_<<" n_cluster_change="<<int(i_end-i_cluster_data)<<std::endl;
#endif
            i_cluster_changeover_update_.resizeNoInitialize(i_end-i_cluster_data);
        }
        */
        
        Ptcl::group_data_mode = GroupDataMode::artificial;
    }

#ifdef CLUSTER_VELOCITY
    //! clear particle group_data.cm to zero 
    /*! update both local and global 
       @param[in,out] _ptcl_soft: global particle
    */
    template <class Tsoft>
    void resetParticleGroupData(Tsoft& _ptcl_soft) {
        Ptcl::group_data_mode = GroupDataMode::cm;
        const PS::S32 n = ptcl_hard_.size();
#pragma omp parallel for 
        for(PS::S32 i=0; i<n; i++){
            // to avoid issue in cluster search with velocity
            auto& pi_cm = ptcl_hard_[i].group_data.cm;
            pi_cm.mass = pi_cm.vel.x = pi_cm.vel.y = pi_cm.vel.z = 0.0;
            PS::S32 adr = ptcl_hard_[i].adr_org;
#ifdef HARD_DEBUG
            assert(adr>=0);
#endif
            _ptcl_soft[adr].group_data.cm  = pi_cm;
        }
    }

    //! set group member particle group_data.cm to c.m. data for search cluster
    /*! update both local and global 
       @param[in,out] _ptcl_soft: global particle
    */
    template <class Tsoft>
    void setParticleGroupDataToCMData(Tsoft& _ptcl_soft) {
        Ptcl::group_data_mode = GroupDataMode::cm;
        auto& ap_manager = manager->ap_manager;
        const PS::S32 n_cluster = n_ptcl_in_cluster_.size();
#pragma omp parallel for schedule(dynamic)
        for(PS::S32 i=0; i<n_cluster; i++){
            const PS::S32 adr_head = n_ptcl_in_cluster_disp_[i];
            const PS::S32 n_ptcl = n_ptcl_in_cluster_[i];
            const PS::S32 n_group = n_group_in_cluster_[i];
            PtclH4* ptcl_local = ptcl_hard_.getPointer(adr_head);

            PS::S32 n_group_offset_local = 0;
            if(n_group>0) {
                for(int k=0; k<n_group; k++) {
                    PS::S32 n_group_in_cluster_offset_k = n_group_in_cluster_offset_[i]+k;
                    PS::S32 ptcl_artificial_adr = adr_first_ptcl_arti_in_cluster_[n_group_in_cluster_offset_k];
                    PS::S32 n_members = n_member_in_group_[n_group_in_cluster_offset_k];
                    // when artificial particles exist
                    if (ptcl_artificial_adr>=0) {
                        auto* pi = &(_ptcl_soft[ptcl_artificial_adr]);
                        auto* pcm = ap_manager.getCMParticles(pi);
                        PS::F64 pcm_mass = pcm->group_data.artificial.getMassBackup();
#ifdef ARTIFICIAL_PARTICLE_DEBUG
                        assert(n_members == ap_manager.getMemberN(pi));
                        ap_manager.checkConsistence(&ptcl_local[n_group_offset_local], pi);
#endif
                        for (int j=n_group_offset_local; j<n_group_offset_local+n_members; j++) {
                            ptcl_local[j].r_search = std::max(pcm->r_search, ptcl_local[j].r_search);
                            auto& pj_cm = ptcl_local[j].group_data.cm;
                            pj_cm.mass  = pcm_mass;
                            pj_cm.vel.x = pcm->vel[0];
                            pj_cm.vel.y = pcm->vel[1];
                            pj_cm.vel.z = pcm->vel[2];
                            PS::S32 adr = ptcl_local[j].adr_org;
                            if(adr>=0) _ptcl_soft[adr].group_data.cm = pj_cm;
                        }
                    }
                    else {
                        // when no artificial particles, calculate c.m. 
#ifdef ARTIFICIAL_PARTICLE_DEBUG
                        // current case only isolated binary 
                        assert(n_members == 2&&n_group==1);
#endif
                        PS::F32 mass_cm=0.0;
                        PS::F32vec vel_cm=PS::F32vec(0.0);
                        for (int j=n_group_offset_local; j<n_group_offset_local+n_members; j++) {
                            auto& pj = ptcl_local[j];
                            mass_cm += pj.mass;
                            vel_cm.x += pj.mass*pj.vel.x;
                            vel_cm.y += pj.mass*pj.vel.y;
                            vel_cm.z += pj.mass*pj.vel.z;
                        }
                        vel_cm /= mass_cm;
                        for (int j=n_group_offset_local; j<n_group_offset_local+n_members; j++) {
                            auto& pj_cm = ptcl_local[j].group_data.cm;
                            pj_cm.mass = mass_cm;
                            pj_cm.vel  = vel_cm;
                            PS::S32 adr = ptcl_local[j].adr_org;
                            if(adr>=0) _ptcl_soft[adr].group_data.cm = pj_cm;
                        }
                    }
                    n_group_offset_local += n_members;
                }
            }
            for (int j=n_group_offset_local; j<n_ptcl; j++) {
                auto& pj_cm = ptcl_local[j].group_data.cm;
                pj_cm.mass  = pj_cm.vel.x = pj_cm.vel.y = pj_cm.vel.z = 0.0;
                PS::S32 adr = ptcl_local[j].adr_org;
                if(adr>=0) _ptcl_soft[adr].group_data.cm = pj_cm;
            }
        }
    }
#endif

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
            _sys[k].pot_tot += _sys[k].mass / manager->r_out_base;
        }
    }

    //! Soft force correction due to different cut-off function
    /*! Use tree neighbor search
      1. first correct for artificial particles use cluster information, 
      2. The c.m. force is substracted from tidal tensor force
      3. c.m. force is replaced by the averaged force on orbital particles
      4. then use tree neighbor search for local cluster real member. 

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
        const PS::S32 n_cluster = n_ptcl_in_cluster_.size();
        const PtclH4* ptcl_local = ptcl_hard_.getPointer();

#pragma omp parallel for schedule(dynamic)
        for (int i=0; i<n_cluster; i++) {  // i: i_cluster
            PS::S32 adr_real_start= n_ptcl_in_cluster_disp_[i];
            PS::S32 adr_real_end= n_ptcl_in_cluster_disp_[i+1];
            // artificial particle group number
            PS::S32 n_group = n_group_in_cluster_[i];
            //PS::S32 n_group_offset = n_group_in_cluster_offset_[i];
            const PS::S32* adr_first_ptcl_arti = &adr_first_ptcl_arti_in_cluster_[n_group_in_cluster_offset_[i]];

            // correction for artificial particles
            correctForceWithCutoffArtificialOneClusterImp(_sys, ptcl_local, adr_real_start, adr_real_end, n_group, adr_first_ptcl_arti, _acorr_flag);

            // obtain correction for real particles in clusters use tree neighbor search
            for (int j=adr_real_start; j<adr_real_end; j++) {
                PS::S64 adr = ptcl_local[j].adr_org;
                // only do for local particles
#ifdef HARD_DEBUG
                if(adr>=0) assert(_sys[adr].id==ptcl_local[j].id);
#endif
                if(adr>=0) correctForceWithCutoffTreeNeighborOneParticleImp<Tpsoft, Ttree, Tepj>(_sys[adr], _tree, _acorr_flag);
            }
        }

        const PS::S32 n_send = _adr_send.size();
#pragma omp parallel for 
        // sending list to other nodes need also be corrected.
        for (int i=0; i<n_send; i++) {
            PS::S64 adr = _adr_send[i];
            correctForceWithCutoffTreeNeighborOneParticleImp<Tpsoft, Ttree, Tepj>(_sys[adr], _tree, _acorr_flag); 
        }
    }

    //! Soft force correction due to different cut-off function
    /* Use cluster member, 
       1. first correct for artificial particles, then for cluster member. 
       2. The c.m. force is substracted from tidal tensor force
       3. c.m. force is replaced by the averaged force on orbital particles

       @param[in] _sys: global particle system, acc is updated
       @param[in] _acorr_flag: flag to do acorr for KDKDK_4TH case
    */
    template <class Tsys>
    void correctForceWithCutoffClusterOMP(Tsys& _sys, const bool _acorr_flag=false) { 
        assert(Ptcl::group_data_mode == GroupDataMode::artificial);

        const PS::S32 n_cluster = n_ptcl_in_cluster_.size();
        auto& ap_manager = manager->ap_manager;
        const PtclH4* ptcl_local = ptcl_hard_.getPointer();
#pragma omp parallel for schedule(dynamic)
        for (int i=0; i<n_cluster; i++) {  // i: i_cluster
            PS::S32 adr_real_start= n_ptcl_in_cluster_disp_[i];
            PS::S32 adr_real_end= n_ptcl_in_cluster_disp_[i+1];
            // artificial particle group number
            PS::S32 n_group = n_group_in_cluster_[i];
            //PS::S32 n_group_offset = n_group_in_cluster_offset_[i];
            const PS::S32* adr_first_ptcl_arti = n_group>0? &adr_first_ptcl_arti_in_cluster_[n_group_in_cluster_offset_[i]] : NULL;

            // correction for artificial particles
            correctForceWithCutoffArtificialOneClusterImp(_sys, ptcl_hard_.getPointer(), adr_real_start, adr_real_end, n_group, adr_first_ptcl_arti, _acorr_flag);

            // obtain correction for real particles in clusters
            for (int j=adr_real_start; j<adr_real_end; j++) {
                PS::S64 adr = ptcl_local[j].adr_org;
#ifdef HARD_DEBUG
                assert(_sys[adr].id==ptcl_local[j].id);
#endif
                //self-potential correction for non-group member, group member has mass zero, so no need correction
                if(_sys[adr].group_data.artificial.isSingle()) _sys[adr].pot_tot += _sys[adr].mass/manager->r_out_base;
                // for binary without artificial particles, correction is needed.
                if (_sys[adr].group_data.artificial.isMember() && _sys[adr].getParticleCMAddress()<0) _sys[adr].pot_tot += _sys[adr].mass/manager->r_out_base; 

                // cluster member
                for (int k=adr_real_start; k<adr_real_end; k++) {
                    if(k==j) continue;
#ifdef KDKDK_4TH
                    if(_acorr_flag) {
                        PS::S64 adr_k = ptcl_local[k].adr_org;
                        calcAcorrShortWithLinearCutoff(_sys[adr], _sys[adr_k]);
                    }
                    else
#endif
                        calcAccPotShortWithLinearCutoff(_sys[adr], ptcl_local[k]);
                }

                // orbital artificial particle
                for (int k=0; k<n_group; k++) {
                    // loop artificial particle orbital
                    PS::S32 k_start = adr_first_ptcl_arti[k];
#ifdef ARTIFICIAL_PARTICLE_DEBUG
                    if (k_start<0) assert(n_group==1&&adr_real_end-adr_real_start==2);
#endif
                    if (k_start<0) continue;
                    auto* porb_k = ap_manager.getOrbitalParticles(&(_sys[k_start]));
                    for (int ki=0; ki<ap_manager.getOrbitalParticleN(); ki++) {
#ifdef KDKDK_4TH
                        if(_acorr_flag) 
                            calcAcorrShortWithLinearCutoff(_sys[adr], porb_k[ki]);
                        else
#endif
                            calcAccPotShortWithLinearCutoff(_sys[adr], porb_k[ki]);
                    }
                }
            
            }
        }
    }

    //! Correct force due to change over factor change 
    /*!
       @param[in,out] _sys: global particle system, acc is updated
       @param[in] _tree: tree for force, to get neighbor list
       @param[in] _adr_send: for connected case, correct sending list
       @param[in] _n_send: sending particle number
     */
    template <class Tsys, class Ttree, class Tepj>
    void correctForceForChangeOverUpdateOMP(Tsys& _sys, Ttree& _tree, 
                                            const PS::S32*  _adr_send=NULL, const PS::S32 _n_send=0) {
        const PS::S32 n_cluster = i_cluster_changeover_update_.size();
        //auto& ap_manager = manager->ap_manager;
#pragma omp parallel for schedule(dynamic)
        for (int i=0; i<n_cluster; i++) {  // i: i_cluster
            PS::S32 i_cluster = i_cluster_changeover_update_[i];
            PS::S32 adr_real_start= n_ptcl_in_cluster_disp_[i_cluster];
            PS::S32 adr_real_end= n_ptcl_in_cluster_disp_[i_cluster+1];
            
            // correction for artificial particles
            /* c.m. force is not sum of orbital force anymore, suppress this part
            // artificial particle group number
            PS::S32 n_group = n_group_in_cluster_[i_cluster];
            const PS::S32* adr_first_ptcl_arti = n_group>0? &adr_first_ptcl_arti_in_cluster_[n_group_in_cluster_offset_[i_cluster]] : NULL;
            for (int j=0; j<n_group; j++) {  // j: j_group
                PS::S32 j_start = adr_first_ptcl_arti[j];
                auto* pj = &(_sys[j_start]);
                
                // loop orbital particles
                auto* porb_j = ap_manager.getOrbitalParticles(pj);
                const PS::S32 n_orb = ap_manager.getOrbitalParticleN();
                for (int k=0; k<n_orb; k++) { 
                    // k: k_ptcl_arti
                    bool changek = porb_j[k].changeover.r_scale_next!=1.0;

                    // loop orbital artificial particle
                    // group
                    for (int kj=0; kj<n_group; kj++) { // group
                        PS::S32 kj_start = adr_first_ptcl_arti[kj];
                        auto* porb_kj = ap_manager.getOrbitalParticles(&_sys[kj_start]);

                        // particle arti orbital
                        if (porb_kj[0].changeover.r_scale_next!=1.0 || changek) {
                            for (int kk=0; kk<n_orb; kk++) {
                                if(&porb_kj[kk]==&porb_j[k]) continue; //avoid same particle
                                // in rare case, equal mass systems with circular orbit can have two members locating at the same positions, this can result in dr2=0 and NaN
                                calcAccChangeOverCorrection(porb_j[k], porb_kj[kk]);
                            }
                        }
                    }

                    //loop real particle
                    for (int kj=adr_real_start; kj<adr_real_end; kj++) {
                        if (ptcl_hard_[kj].changeover.r_scale_next!=1.0 || changek) {
                            calcAccChangeOverCorrection(porb_j[k], ptcl_hard_[kj]);
                        }
                    }
                }

                // correct c.m. average force
                ap_manager.correctOrbitalParticleForce(pj);
            }
            */

            // correction for real particles
            for (int j=adr_real_start; j<adr_real_end; j++) {
                PS::S64 adr = ptcl_hard_[j].adr_org;
                if(adr>=0) {
                    bool change_i = _sys[adr].changeover.r_scale_next!=1.0;
                    Tepj * ptcl_nb = NULL;
                    PS::S32 n_ngb = _tree.getNeighborListOneParticle(_sys[adr], ptcl_nb);
                    for(PS::S32 k=0; k<n_ngb; k++){
                        if (ptcl_nb[k].id == _sys[adr].id) continue;

                        if (ptcl_nb[k].r_scale_next!=1.0 || change_i) 
                            calcAccChangeOverCorrection(_sys[adr], ptcl_nb[k]);
                    }
                }
                // update changeover
                ptcl_hard_[j].changeover.updateWithRScale();
                if(adr>=0) _sys[adr].changeover.updateWithRScale();
            }
            
        }
#pragma omp parallel for 
        // sending list to other nodes need also be corrected.
        for (int i=0; i<_n_send; i++) {
            PS::S64 adr = _adr_send[i];
            bool change_i = _sys[adr].changeover.r_scale_next!=1.0;
            Tepj * ptcl_nb = NULL;
            PS::S32 n_ngb = _tree.getNeighborListOneParticle(_sys[adr], ptcl_nb);
            for(PS::S32 k=0; k<n_ngb; k++){
                if (ptcl_nb[k].id == _sys[adr].id) continue;
                
                if (ptcl_nb[k].r_scale_next!=1.0 || change_i) 
                    calcAccChangeOverCorrection(_sys[adr], ptcl_nb[k]);
            }
            _sys[adr].changeover.updateWithRScale();
        }
        
    }


    //! Soft force correction due to different cut-off function
    /*! Use tree neighbor search for all particles.
       c.m. force is replaced by the averaged force on orbital particles
       Tidal tensor particle subtract the c.m. acc
       @param[in] _sys: global particle system, acc is updated
       @param[in] _tree: tree for force
       @param[in] _adr_ptcl_artificial_start: start address of artificial particle in _sys
    */
    template <class Tsys, class Tpsoft, class Ttree, class Tepj>
    void correctForceWithCutoffTreeNeighborOMP(Tsys& _sys,
                                               Ttree& _tree,
                                               const PS::S32 _adr_ptcl_artificial_start,
                                               const bool _acorr_flag=false) {
        // for artificial particle
        const PS::S32 n_tot = _sys.getNumberOfParticleLocal();

#pragma omp parallel for schedule(dynamic)
        for (int i=0; i<n_tot; i++) {
            correctForceWithCutoffTreeNeighborOneParticleImp<Tpsoft, Ttree, Tepj>(_sys[i], _tree, _acorr_flag);
        }
        auto& ap_manager = manager->ap_manager;
        const PS::S32 n_artificial = ap_manager.getArtificialParticleN();
#ifdef HARD_DEBUG
        assert((n_tot-_adr_ptcl_artificial_start)%n_artificial==0);
#endif
#pragma omp parallel for schedule(dynamic)
        for (int i=_adr_ptcl_artificial_start; i<n_tot; i+=n_artificial)
            ap_manager.correctArtficialParticleForce(&(_sys[i]));
    }

};

