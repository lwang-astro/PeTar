#pragma once

#define ID_PHASE_SHIFT 4
#define ID_PHASE_MASKER 0xF

class ArtificialParticleManager{
    PS::S32 n_split_;    // oribital particle splitting number
    PS::S32 n_artificial_;     // number of artificial particles
    PS::S32 index_offset_tt_;  // tidal tensor particle starting index
    PS::S32 index_offset_orb_; // Orbital particle starting index
    PS::S32 index_cm_;        // center of mass index


public:
    PS::S32 n_split;    // test
    PS::F64 r_tidal_tensor;
    PS::S64 id_offset;
    PS::F64 G; // gravitational constant

    ArtificialParticleManager(): n_split_(-1), n_artificial_(-1), index_offset_tt_(0), index_offset_orb_(8), index_cm_(-1), r_tidal_tensor(-1.0), id_offset(-1), G(-1.0) {}

    //! check paramters
    bool checkParams() {
        ASSERT(n_split_>=4);
        ASSERT(n_artificial_>=9);
        ASSERT(index_cm_>=8);
        ASSERT(r_tidal_tensor>=0.0);
        ASSERT(id_offset>0);
        ASSERT(G>0);
        return true;
    }

    //! create artificial particles 
    /*! First 8 are tidal tensor particles; 9-2*n_split are orbitial sample particles; 2*n_split+1 is c.m.  

      Particle status:
      member particle: 1
      Artificial particles:
      0: left member N
      1: right member N
      2-2*n_split: _data_to_store; otherwises index+1
      pcm: n_members

      Particle mass_bk:
      pcm: mass(cm)

      @param[out]    _ptcl_new: artificial particles that will be added
      @param[in]     _bin: binary tree root
      @param[in]     _data_to_store: array of data to be stored in the status of artificial particles
      @param[in]     _n_data: number of data
    */
    template <class Tptcl>
    void createArtificialParticles(Tptcl* _ptcl_artificial,
                                   COMM::BinaryTree<Tptcl> &_bin, 
                                   const PS::F64* _data_to_store,
                                   const PS::S32 _n_data) {
        // set id and status.d 
        for (int i=0; i<n_split_; i++) {
            for(int j=0; j<2; j++) {
                Tptcl* pj = &_ptcl_artificial[2*i+j];
                Tptcl* member = _bin.getMember(j);
                pj->id = id_offset + abs(member->id)*n_split_ +i;
                pj->status.d = 2*i+j+1;
            }
        }

        // First 8 is used for tidal tensor points
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(r_tidal_tensor<=_bin.changeover.getRin());
#endif
        TidalTensor::createTidalTensorMeasureParticles(_ptcl_artificial, *((Tptcl*)&_bin), r_tidal_tensor);

        // remaining is used for sample points
        const PS::F64 dE = 8.0*atan(1.0)/(n_split_-4);
        PS::F64 mnormal=0.0;
        for (int i=4; i<n_split_; i++) {
            PS::S32 iph = i-4;

            // center_of_mass_shift(*(Tptcl*)&_bin,p,2);
            // generate particles at different orbitial phase
            _bin.orbitToParticle(_ptcl_artificial[2*i], _ptcl_artificial[2*i+1], _bin, dE*iph, G);

            // use velocity to weight mass
            PS::F64vec dvvec= _ptcl_artificial[2*i].vel - _ptcl_artificial[2*i+1].vel;
            PS::F64 odv = 1.0/std::sqrt(dvvec*dvvec);

            for (int j=0; j<2; j++) {
                Tptcl* pj = &_ptcl_artificial[2*i+j];
                Tptcl* member = _bin.getMember(j);

                // set mass
                pj->mass = member->mass * odv;

                // center_of_mass_correction 
                pj->pos += _bin.pos;
                pj->vel += _bin.vel;

#ifdef ARTIFICIAL_PARTICLE_DEBUG
                assert(member->mass>0);
#endif
            }
            mnormal += odv;
        }

        // normalized the mass of each particles to keep the total mass the same as c.m. mass
        PS::F64 mfactor = 1.0/mnormal;
        for (int i=4; i<n_split_; i++) 
            for (int j=0; j<2; j++) 
                _ptcl_artificial[2*i+j].mass *= mfactor;

        // store the component member number 
        for (int j=0; j<2; j++) {
            _ptcl_artificial[j].status.d = _bin.isMemberTree(j) ? ((COMM::BinaryTree<Tptcl>*)(_bin.getMember(j)))->getMemberN() : 1; 
        }

        // store the additional data (should be positive) 
        // ensure the data size is not overflow
        assert(_n_data+2<2*n_split_);
        for (int j=0; j<_n_data; j++) _ptcl_artificial[j+2].status.d = _data_to_store[j];

        // last member is the c.m. particle
        Tptcl* pcm;
        pcm = &_ptcl_artificial[2*n_split_];
        pcm->mass_bk.d = _bin.mass;
        pcm->mass = 0.0;
        pcm->pos = _bin.pos;
        pcm->vel = _bin.vel;
        pcm->id  = - std::abs(_bin.id);
        pcm->status.d = _bin.getMemberN();
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

    //! get artificial particle total number
    PS::S32 getTidalTensorParticleN() const {
        return 8;
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

    //! get center of mass id
    template <class Tptcl>
    PS::S64 getCMID(const Tptcl* _ptcl_list) const {
        return -PS::S64(_ptcl_list[index_cm_].id);
    }

    //! get stored data 
    template <class Tptcl>
    PS::F64 getStoredData(const Tptcl* _ptcl_list, const PS::S32 _index) const {
        return _ptcl_list[_index+2].status.d;
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
             <<"id_offset        : "<<id_offset<<std::endl
             <<"n_split          : "<<n_split_<<std::endl;
    }    

#ifdef ARTIFICIAL_PARTICLE_DEBUG
    template <class Tptcl, class Tpart>
    void checkConsistence(Tptcl* _ptcl_member, Tpart* _ptcl_artificial) {
        // check id 
        PS::S64 id_min = _ptcl_member[0].id;
        for(int j=0; j<getMemberN(_ptcl_artificial); j++)  id_min = std::min(id_min,_ptcl_member[j].id);
        assert(getCMID(_ptcl_artificial) == id_min);
        // not consistent if first member is single and second is binary
        //PS::S32 id_mem[2];
        //id_mem[0] = _ptcl_member[0].id;
        //id_mem[1] = _ptcl_member[getLeftMemberN(_ptcl_artificial)].id;
        //// id_offset unknown, try to substract id information via calculation between neighbor particles
        //for (int j=0; j<getArtificialParticleN()-1; j+=2) {
        //    // first member
        //    PS::S32 id_offset_j1 = _ptcl_artificial[j].id - j/2- id_mem[0]*(getArtificialParticleN()-1)/2;
        //    // second member
        //    PS::S32 id_offset_j2 = _ptcl_artificial[j+1].id - j/2 - id_mem[1]*(getArtificialParticleN()-1)/2;
        //    assert(id_offset_j1==id_offset_j2);
        //}

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

};
