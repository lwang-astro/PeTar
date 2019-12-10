#pragma once

#define ID_PHASE_SHIFT 4
#define ID_PHASE_MASKER 0xF

#include "search_group.hpp"
#include "stability.hpp"
#include "changeover.hpp"

class ArtificialParticleManager{
    PS::S32 n_split_;    // oribital particle splitting number
    PS::S32 n_artificial_;     // number of artificial particles
    PS::S32 index_offset_tt_;  // tidal tensor particle starting index
    PS::S32 index_offset_orb_; // Orbital particle starting index
    PS::S32 index_cm_;        // center of mass index

    template <class Tptcl>
    struct BinPar {
        Tptcl* adr_ref; PS::S32* group_list; PS::S32 n; ChangeOver* changeover;
    };

    // set binary tree parameters 
    template <class Tptcl>
    static PS::S64 setBinChangeOverIDAndGetMemberAdrIter(BinPar<Tptcl>& _par, const PS::S64& _id1, const PS::S64& _id2, COMM::BinaryTree<Tptcl>& _bin) {
        // set bin id as the left member id
        if (_id1<0) _bin.id = _bin.getLeftMember()->id;
        else _bin.id = _id1;
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(_bin.id>0);
#endif

        // leaf case
        if (_id1<0  && _id2<0) {
            // collect address
            for (int k=0; k<2; k++) {
                Tptcl* member = _bin.getMember(k);
#ifdef ARTIFICIAL_PARTICLE_DEBUG
                assert(member->mass>0.0);
#endif                
                _par.group_list[_par.n++] = member - _par.adr_ref;
                member->status.d = 1; // used for number count later
            }
        }

        //set changeover to the same (root) one
        _bin.changeover = *(_par.changeover);
        
        return _bin.id;
    }

    // get new changeover and rsearch for c.m.
    template <class Tchp, class Tptcl>
    static Tchp* calcBinChangeOverAndRSearchIter (Tchp*& _ch, COMM::BinaryTree<Tptcl>& _bin) {
        _bin.changeover.setR(_bin.mass*_ch->mean_mass_inv, _ch->rin, _ch->rout);
        _bin.Ptcl::calcRSearch(_ch->dt_tree);
        return _ch;
    };

    //! generate kepler sampling artificial particles
    /*  Particle status:
        member particle: 1
        Artificial particles:
        0: left member N
        1: right member N
        2: id_cluster
        3: id_group
        4-2*n_split: index+1
        pcm: n_members
        Particle mass_bk:
        pcm: mass(cm)

        @param[in]     _i_cluster: cluster index
        @param[in]     _i_group: group index in the cluster
        @param[in,out] _ptcl_in_cluster: particle data in local cluster
        @param[out]    _ptcl_new: artificial particles that will be added
        @param[out]    _group_ptcl_adr_list: group member particle index list in _ptcl_in_cluster 
        @param[in]     _bin: binary tree root
     */
    template <class Tptcl>
    void keplerOrbitGenerator(const PS::S32 _i_cluster,
                              const PS::S32 _i_group,
                              Tptcl* _ptcl_in_cluster,
                              PS::ReallocatableArray<Tptcl> & _ptcl_new,
                              PS::S32 *_group_ptcl_adr_list,
                              COMM::BinaryTree<Tptcl> &_bin) {
        const PS::F64 dE = 8.0*atan(1.0)/(n_split_-4);
        //! Set member particle status=1, return orderd particle member index list
        /*  _adr_ref: ptcl_org first particle address as reference to calculate the particle index.
        */
        BinPar<Tptcl> bin_par = {_ptcl_in_cluster, _group_ptcl_adr_list, 0, &_bin.changeover};
        _bin.processTreeIter(bin_par, (PS::S64)-1, (PS::S64)-1, setBinChangeOverIDAndGetMemberAdrIter);
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        const PS::S32 n_members = _bin.getMemberN();
        assert(bin_par.n==n_members);
#endif                

        // Make sure the _ptcl_new will not make new array due to none enough capacity during the following loop, otherwise the p[j] pointer will point to wrong position
        const int np = _ptcl_new.size();
        _ptcl_new.increaseSize(n_artificial_);
        Tptcl* p = &_ptcl_new[np];
        
        // set id and status.d 
        for (int i=0; i<n_split_; i++) {
            for(int j=0; j<2; j++) {
                Tptcl* pj = &p[2*i+j];
                Tptcl* member = _bin.getMember(j);
                pj->id = id_offset + abs(member->id)*n_split_ +i;
                pj->status.d = 2*i+j+1;
            }
        }

        // First 8 is used for tidal tensor points
        TidalTensor::createTidalTensorMeasureParticles(p, *((Tptcl*)&_bin), r_tidal_tensor);

        // use c.m. r_search and changeover
        for (int i=0; i<4; i++) {
            for (int j=0; j<2; j++) {
                Tptcl* pj = &p[2*i+j];
                pj->r_search =_bin.getMember(j)->r_search;
                pj->changeover = _bin.changeover;
            }
        }

        // remaining is used for sample points
        PS::F64 mnormal=0.0;
        for (int i=4; i<n_split_; i++) {
            PS::S32 iph = i-4;

            // center_of_mass_shift(*(Tptcl*)&_bin,p,2);
            // generate particles at different orbitial phase
            _bin.orbitToParticle(p[2*i], p[2*i+1], _bin, dE*iph, G);

            // use velocity to weight mass
            PS::F64vec dvvec= p[2*i].vel - p[2*i+1].vel;
            PS::F64 odv = 1.0/std::sqrt(dvvec*dvvec);

            for (int j=0; j<2; j++) {
                Tptcl* pj = &p[2*i+j];
                Tptcl* member = _bin.getMember(j);

                // set mass
                pj->mass = member->mass * odv;

                // center_of_mass_correction 
                pj->pos += _bin.pos;
                pj->vel += _bin.vel;

#ifdef ARTIFICIAL_PARTICLE_DEBUG
                assert(member->mass>0);
#endif

                // use member changeover, if new changeover is different, record the scale ratio 
                pj->changeover =  member->changeover;
                if (abs(pj->changeover.getRin()-_bin.changeover.getRin())>1e-10) {
                    pj->changeover.r_scale_next = _bin.changeover.getRin()/pj->changeover.getRin();
                    pj->r_search = std::max(pj->r_search, _bin.r_search);
#ifdef ARTIFICIAL_PARTICLE_DEBUG
                    assert(pj->r_search > pj->changeover.getRout());
#endif 
                }
                else pj->r_search = _bin.r_search;

#ifdef ARTIFICIAL_PARTICLE_DEBUG
                //check rsearch consistence:
                PS::F64 rsearch_bin = _bin.r_search+_bin.semi*(1+_bin.ecc);
                PS::F64vec dp = pj->pos-_bin.pos;
                PS::F64 dr = dp*dp;
                assert(dr<=rsearch_bin*rsearch_bin);
#endif
            }

            mnormal += odv;
        }

        // normalized the mass of each particles to keep the total mass the same as c.m. mass
        PS::F64 mfactor = 1.0/mnormal;
        for (int i=4; i<n_split_; i++) 
            for (int j=0; j<2; j++) 
                p[2*i+j].mass *= mfactor;


        // store the component member number 
        for (int j=0; j<2; j++) {
            p[j].status.d = _bin.isMemberTree(j) ? ((COMM::BinaryTree<Tptcl>*)(_bin.getMember(j)))->getMemberN() : 1; 
        }
        // store the i_cluster and i_group for identify artificial particles, +1 to avoid 0 value (status.d>0)
        p[2].status.d = _i_cluster+1;
        p[3].status.d = _i_group+1;

        // last member is the c.m. particle
        Tptcl* pcm;
        pcm = &_ptcl_new.back();
        pcm->mass_bk.d = _bin.mass;
        pcm->mass = 0.0;
        pcm->pos = _bin.pos;
        pcm->vel = _bin.vel;
        pcm->id  = - std::abs(_bin.id);
        pcm->r_search = _bin.r_search;
        pcm->changeover = _bin.changeover;

        pcm->r_search += _bin.semi*(1+_bin.ecc);  // depend on the mass ratio, the upper limit distance to c.m. from all members and artificial particles is apo-center distance

        pcm->status.d = _bin.getMemberN();
    }


public:
    PS::S32 n_split;    // test
    PS::F64 r_tidal_tensor;
    PS::F64 r_in_base;
    PS::F64 r_out_base;
    PS::S64 id_offset;
    PS::F64 G; // gravitational constant

    ArtificialParticleManager(): n_split_(-1), n_artificial_(-1), index_offset_tt_(0), index_offset_orb_(8), index_cm_(-1), r_tidal_tensor(-1.0), r_in_base(-1.0), r_out_base(-1.0), id_offset(-1), G(-1.0) {}

    //! check paramters
    bool checkParams() {
        ASSERT(n_split_>=4);
        ASSERT(n_artificial_>=9);
        ASSERT(index_cm_>=8);
        ASSERT(r_tidal_tensor>=0.0);
        ASSERT(r_in_base>0.0);
        ASSERT(r_out_base>0.0);
        ASSERT(id_offset>0);
        ASSERT(G>0);
        return true;
    }

    //! set orbital particle split number
    void setOrbitalParticleSplitN(const PS::S32& _n_split) {
        if (_n_split>(1<<ID_PHASE_SHIFT)) {
            std::cerr<<"Error! ID_PHASE_SHIFT is too small for phase split! shift bit: "<<ID_PHASE_SHIFT<<" n_split_: "<<_n_split<<std::endl;
            abort();
        }
        ASSERT(_n_split>=4);
        n_split_  = _n_split;
        n_artificial_   = 2*_n_split+1;
        index_offset_tt_  = 0;
        index_offset_orb_ = 8;
        index_cm_ = 2*_n_split;
        n_split = _n_split;
    }
      
    //! generate artificial particles,
    /*  
        @param[in]     _i_cluster: cluster index
        @param[in,out] _ptcl_in_cluster: particle data
        @param[in]     _n_ptcl: total number of particle in _ptcl_in_cluster.
        @param[out]    _ptcl_artificial: artificial particles that will be added
        @param[out]    _n_groups: number of groups in current cluster
        @param[in,out] _groups: searchgroup class, which contain 1-D group member index array, will be reordered by the minimum distance chain for each group
        @param[in,out] _empty_list: the list of _ptcl_in_cluster that can be used to store new artificial particles, reduced when used
        @param[in]     _dt_tree: tree time step for calculating r_search
     */
    template <class Tptcl>
    void createArtificialParticles(const PS::S32 _i_cluster,
                             Tptcl* _ptcl_in_cluster,
                             const PS::S32 _n_ptcl,
                             PS::ReallocatableArray<Tptcl> & _ptcl_artificial,
                             PS::S32 &_n_groups,
                             SearchGroup<Tptcl>& _groups,
                             const PS::F64 _dt_tree) {

        PS::S32 group_ptcl_adr_list[_n_ptcl];
        PS::S32 group_ptcl_adr_offset=0;
        _n_groups = 0;
        for (int i=0; i<_groups.getNumberOfGroups(); i++) {
            PS::ReallocatableArray<COMM::BinaryTree<Tptcl>> bins;   // hierarch binary tree
            bins.reserve(_n_groups);

            const PS::S32 n_members = _groups.getNumberOfGroupMembers(i);
#ifdef ARTIFICIAL_PARTICLE_DEBUG
            assert(n_members<ARRAY_ALLOW_LIMIT);
#endif        
            PS::S32* member_list = _groups.getMemberList(i);
            bins.resizeNoInitialize(n_members-1);
            // build hierarch binary tree from the minimum distant neighbors

            COMM::BinaryTree<Tptcl>::generateBinaryTree(bins.getPointer(), member_list, n_members, _ptcl_in_cluster);

            // reset status to 0
            for (int j=0; j<n_members; j++) _ptcl_in_cluster[member_list[j]].status.d=0;

            struct {PS::F64 mean_mass_inv, rin, rout, dt_tree; } ch = {Tptcl::mean_mass_inv, r_in_base, r_out_base, _dt_tree};
            auto* chp = &ch;
            // get new changeover and rsearch for c.m.
            bins.back().processRootIter(chp, calcBinChangeOverAndRSearchIter);
            
            // stability check and break groups
            Stability<Tptcl> stab;
            // be careful, here t_crit should be >= hard slowdown_timescale_max to avoid using slowdown for wide binaries
            stab.t_crit = _dt_tree;
            stab.stable_binary_tree.reserve(_n_groups);
            stab.findStableTree(bins.back());

            for (int i=0; i<stab.stable_binary_tree.size(); i++) {
                keplerOrbitGenerator(_i_cluster, _n_groups, _ptcl_in_cluster, _ptcl_artificial, &group_ptcl_adr_list[group_ptcl_adr_offset], *stab.stable_binary_tree[i]);
                group_ptcl_adr_offset += stab.stable_binary_tree[i]->getMemberN();
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
            if(_ptcl_in_cluster[i_group].status.d==0) {
                while(_ptcl_in_cluster[i_single_front].status.d==0) {
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

        //for (int i=0; i<_empty_list.size(); i++) {
        //    PS::S32 ik = _empty_list[i];
        //    _ptcl_in_cluster[ik].mass = 0.0;
        //    _ptcl_in_cluster[ik].id = -1;
        //    _ptcl_in_cluster[ik].status.d = -1;
        //}

    }

    //! correct orbitial particles force
    /*!
      replace c.m. force by the averaged force on orbital particles
      @param[in,out] _ptcl_artificial: one group of artificial particles 
    */
    template <class Tptcl>
    void correctOrbitalParticleForce(Tptcl* _ptcl_artificial) {
        auto* pcm = getCMParticles(_ptcl_artificial);
        PS::F64vec& acc_cm = pcm->acc;
        
        acc_cm=PS::F64vec(0.0);
        PS::F64 m_ob_tot = 0.0;

        auto* porb = getOrbitalParticles(_ptcl_artificial);
        const PS::S32 n_orb = getOrbitalParticleN();
        for (int k=0; k<n_orb; k++) {
            acc_cm += porb[k].mass*porb[k].acc; 
            m_ob_tot += porb[k].mass;
        }
        acc_cm /= m_ob_tot;

#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(abs(m_ob_tot-pcm->mass_bk.d)<1e-10);
#endif
    }

    //! correct artificial particles force for furture use
    /*! 
      substract c.m. force (acc) from tidal tensor force (acc)\n
      replace c.m. force by the averaged force on orbital particles
      @param[in,out] _ptcl_artificial: one group of artificial particles 
    */
    template <class Tptcl>
    void correctArtficialParticleForce(Tptcl* _ptcl_artificial) {
        auto* pcm = getCMParticles(_ptcl_artificial);
        // substract c.m. force (acc) from tidal tensor force (acc)
        auto* ptt = getTidalTensorParticles(_ptcl_artificial);
        TidalTensor::subtractCMForce(ptt, *pcm);

        // After c.m. force used, it can be replaced by the averaged force on orbital particles
        correctOrbitalParticleForce(_ptcl_artificial);
    }

    //! get oribital particle list address from a artificial particle array
    template <class Tptcl>
    Tptcl* getOrbitalParticles(Tptcl* _ptcl_list)  {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert((PS::S32)(_ptcl_list[index_offset_orb_].status.d)==index_offset_orb_+1);
#endif
        return &_ptcl_list[index_offset_orb_];
    }

    //! get tidal tensor particle list address from a artificial particle array
    template <class Tptcl>
    Tptcl* getTidalTensorParticles(Tptcl* _ptcl_list) {
        return &_ptcl_list[index_offset_tt_];
    }

    //! get c.m. particle list address from a artificial particle array
    template <class Tptcl>
    Tptcl* getCMParticles(Tptcl* _ptcl_list)  {
        return &_ptcl_list[index_cm_];
    }

    //! get artificial particle total number
    PS::S32 getArtificialParticleN() const {
        return n_artificial_;
    }

    //! get orbitial particle number 
    PS::S32 getOrbitalParticleN() const {
        return n_artificial_ - 9;
    }

    //! get left member number
    template <class Tptcl>
    PS::S32 getLeftMemberN(const Tptcl* _ptcl_list) const {
        return PS::S32(_ptcl_list[0].status.d);
    }

    //! get left member number
    template <class Tptcl>
    PS::S32 getMemberN(const Tptcl* _ptcl_list) const {
        return PS::S32(_ptcl_list[index_cm_].status.d);
    }

    //! get left member number
    template <class Tptcl>
    PS::S32 getRightMemberN(const Tptcl* _ptcl_list) const {
        return PS::S32(_ptcl_list[1].status.d);
    }

    //! get cluster id
    template <class Tptcl>
    PS::S32 getClusterID(const Tptcl* _ptcl_list) const {
        return PS::S32(_ptcl_list[2].status.d)-1;
    }

    //! get group id
    template <class Tptcl>
    PS::S32 getGroupID(const Tptcl* _ptcl_list) const {
        return PS::S32(_ptcl_list[3].status.d)-1;
    }

    //! get first member id
    template <class Tptcl>
    PS::S64 getFirstMemberID(const Tptcl* _ptcl_list) const {
        return -PS::S64(_ptcl_list[index_cm_].id);
    }

    //! write class data to file with binary format
    /*! @param[in] _fp: FILE type file for output
     */
    void writeBinary(FILE *_fp) {
        fwrite(this, sizeof(*this), 1, _fp);
    }    

    //! read class data to file with binary format
    /*! @param[in] _fp: FILE type file for reading
     */
    void readBinary(FILE *_fin) {
        size_t rcount = fread(this, sizeof(*this), 1, _fin);
        if (rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            abort();
        }
    }    

    //! print parameters
    void print(std::ostream & _fout) const{
        _fout<<"r_tidal_tensor   : "<<r_tidal_tensor<<std::endl
             <<"r_in_base        : "<<r_in_base<<std::endl
             <<"r_out_base       : "<<r_out_base<<std::endl
             <<"id_offset        : "<<id_offset<<std::endl
             <<"n_split          : "<<n_split_<<std::endl;
    }    

#ifdef ARTIFICIAL_PARTICLE_DEBUG
    template <class Tptcl, class Tpart>
    void checkConsistence(Tptcl* _ptcl_member, Tpart* _ptcl_artificial) {
        // check id 
        assert(getFirstMemberID(_ptcl_artificial) == _ptcl_member[0].id);
        PS::S32 id_mem[2];
        id_mem[0] = _ptcl_member[0].id;
        id_mem[1] = _ptcl_member[getLeftMemberN(_ptcl_artificial)].id;
        // id_offset unknown, try to substract id information via calculation between neighbor particles
        for (int j=0; j<getArtificialParticleN()-1; j+=2) {
            // first member
            PS::S32 id_offset_j1 = _ptcl_artificial[j].id - j/2- id_mem[0]*(getArtificialParticleN()-1)/2;
            // second member
            PS::S32 id_offset_j2 = _ptcl_artificial[j+1].id - j/2 - id_mem[1]*(getArtificialParticleN()-1)/2;
            assert(id_offset_j1==id_offset_j2);
        }

        // check whether c.m. pos. are consistent
        // Cannot do velocity check because cm is not kicked
        PS::F64 mass_cm_check=0.0;
        PS::F64vec pos_cm_check=PS::F64vec(0.0);
        
        for(int j=0; j<getMemberN(_ptcl_artificial); j++) {
            mass_cm_check += _ptcl_member[j].mass;
            pos_cm_check +=  _ptcl_member[j].pos*_ptcl_member[j].mass;
        }
        
        pos_cm_check /= mass_cm_check;

        auto* pcm = getCMParticles(_ptcl_artificial);
        assert(abs(mass_cm_check-pcm->mass)<1e-10);
        PS::F64vec dpos = pos_cm_check-pcm->pos;
        assert(abs(dpos*dpos)<1e-20);
    }
#endif

//    void checkRoutChange(PS::ReallocatableArray<RCList> & r_out_change_list,
//                         PS::ReallocatableArray<PS::ReallocatableArray<PS::S32>> & group_list,
//                         Tptcl* ptcl){
//        for (int i=0; i<group_list.size(); i++) {
//            PS::S64 r_out_max = 0.0;
//            for (int j=0; j<group_list[i].size(); i++) {
//                PS::S32 k = group_list[i][j];
//                r_out_max = std::max(r_out_max,ptcl[k].r_out);
//            }
//            for (int j=0; j<group_list[i].size(); i++) {
//                PS::S32 k = group_list[i][j];
//                if(ptcl[k].r_out != r_out_max) {
//                    r_out_change_list.push_back(RCList(k,ptcl[k].r_out));
//#ifdef ARTIFICIAL_PARTICLE_DEBUG
//                    std::cerr<<"Rout change detected, p["<<k<<"].r_out: "<<ptcl[k].r_out<<" -> "<<r_out_max<<std::endl;
//#endif
//                    ptcl[k].r_out = r_out_max;
//                }
//            }
//        }
//    }


};
