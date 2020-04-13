#pragma once

#include "Common/binary_tree.h"
#include "tidal_tensor.hpp"

#ifdef ORBIT_SAMPLING
#include "orbit_sampling.hpp"
typedef OrbitalSamplingManager OrbitManager;
#else
#include "pseudoparticle_multipole.hpp"
typedef PseudoParticleMultipoleManager OrbitManager;
#endif

//! class to store necessary information for using artificial particles
/*!
                  single    c.m.             members          initial     artificial
      mass_backup 0         mass      (+)    mass     (+)     -LARGE       default is 0.0 / data stored (-/0)
      status      0         n_members (+)    c.m. adr (-)     -LARGE       position in artificial particle array / data stored (+)
 */
class ArtificialParticleInformation{
private:
    PS::F64 mass_backup;
    PS::F64 status;

public:
    //! initializer
    ArtificialParticleInformation(): mass_backup(-PS::LARGE_FLOAT), status(-PS::LARGE_FLOAT) {}

    //! set particle type to member
    /*! @param[in] _mass: mass to backup
     */
    void setParticleTypeToMember(const PS::F64 _mass = PS::LARGE_FLOAT, const PS::F64 _status = -PS::LARGE_FLOAT) {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(_status<0.0);
#endif
        mass_backup = _mass;
        status = _status;
    };

    //! return whether the particle type is member
    bool isMember() const {
        return (status<0.0);
    }

    //! set particle type to artificial
    /*! @param[in] _status: status to save
     */
    void setParticleTypeToCM(const PS::F64 _mass_backup = PS::LARGE_FLOAT, const PS::F64 _status=PS::LARGE_FLOAT) {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(_status>0.0&&_mass_backup>0.0);
#endif
        status = _status;
        mass_backup = _mass_backup;
    };

    //! return whether the particle type is c.m.
    bool isCM() const {
        return (status>0.0 && mass_backup>0.0);
    }

    //! set backup mass
    void setMassBackup(const PS::F64 _mass) {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert((isMember()||isCM())&&_mass>0.0);
#endif
        mass_backup = _mass;
    }

    //! get backup mass
    PS::F64 getMassBackup() const {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(isMember()||isCM());
#endif
        return mass_backup;
    }

    //! set particle type to single
    void setParticleTypeToSingle() {
        mass_backup = 0.0;
        status = 0.0;
    };

    //! return whether the particle type is single
    bool isSingle() const {
        return (status==0.0 && mass_backup==0.0);
    }

    //! set particle type to artificial
    /*! @param[in] _status: status to save
     */
    void setParticleTypeToArtificial(const PS::F64 _status=PS::LARGE_FLOAT, const PS::F64 _mass_backup = 0.0) {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(_status>0.0);
#endif
        status = _status;
        mass_backup = _mass_backup;
    };

    //! return whether the particle type is artificial
    bool isArtificial() const {
        return (status>0.0);
    }

    //! store one data (only in artificial particles)
    /*! positive data stored in status, otherwise in mass_backup;
     */
    void storeData(const PS::F64 _data) {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(isArtificial());
#endif
        if (_data>0.0) status = _data;
        else      mass_backup = _data;
    }

    //! get stored data (only in artificial particles)
    /*! @param[in] is_positive: indicate whether stored data is positive (true) or negative (false)
     */
    PS::F64 getData(const bool is_positive) const {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(isArtificial());
#endif
        if (is_positive) return status;
        else             return mass_backup;
    }

    //! set status
    void setStatus(const PS::F64 _status) {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert((isMember()&&_status<0.0)||(isArtificial()&&_status>0.0));
#endif        
        status = _status;
    }

    //! get status
    PS::F64 getStatus() const {
        return status;
    }

    //! return whether the particle is unused
    bool isUnused() const {
        return (status<0.0 && mass_backup<0.0);
    }

    //! set particle type to unused
    void setParticleTypeToUnused() {
        status = - NUMERIC_FLOAT_MAX;
        mass_backup = - NUMERIC_FLOAT_MAX;
    }

    //! print titles of class members using column style
    /*! print titles of class members in one line for column style
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
     */
    static void printColumnTitle(std::ostream & _fout, const int _width=20) {
        _fout<<std::setw(_width)<<"mass_bk"
             <<std::setw(_width)<<"status";
    }

    //! print data of class members using column style
    /*! print data of class members in one line for column style. Notice no newline is printed at the end
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
     */
    void printColumn(std::ostream & _fout, const int _width=20){
        _fout<<std::setw(_width)<<mass_backup
             <<std::setw(_width)<<status;
    }

    //! write class data with ASCII format
    /*! @param[in] _fout: file IO for write
     */
    void writeAscii(FILE* _fout) const{
        fprintf(_fout, "%26.17e %26.17e ", 
                this->mass_backup, this->status);
    }

    //! write class data with BINARY format
    /*! @param[in] _fout: file IO for write
     */
    void writeBinary(FILE* _fin) const{
        fwrite(this, sizeof(ArtificialParticleInformation), 1, _fin);
    }

    //! read class data with ASCII format
    /*! @param[in] _fin: file IO for read
     */
    void readAscii(FILE* _fin) {
        PS::S64 rcount=fscanf(_fin, "%lf %lf ",
                              &this->mass_backup, &this->status);
        if (rcount<2) {
            std::cerr<<"Error: Data reading fails! requiring data number is 2, only obtain "<<rcount<<".\n";
            abort();
        }
    }

    //! read class data with BINARY format
    /*! @param[in] _fin: file IO for read
     */
    void readBinary(FILE* _fin) {
        size_t rcount = fread(this, sizeof(ArtificialParticleInformation), 1, _fin);
        if (rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            abort();
        }
    }

    //! print parameters
    void print(std::ostream & _fout) const{
        _fout<<" mass_bk="<<mass_backup
             <<" status="<<status;
    }    
    
};

//! class to organize artificial particles
class ArtificialParticleManager{
public:
    // be careful to make consistent read/writeBinary and operator = if modified
    PS::F64 r_tidal_tensor;  ///> tidal tensor maximum distance of particles
    PS::S64 id_offset;       ///> offset to generate ID for artificial particles
    PS::F64 gravitational_constant; ///> gravitational constant
    OrbitManager orbit_manager; ///> orbit particle method

    //! initializer, require setParticleSplitN later to initialize the size
    ArtificialParticleManager(): r_tidal_tensor(-1.0), id_offset(-1), gravitational_constant(-1.0), orbit_manager() {}

    //! check paramters
    bool checkParams() {
        ASSERT(TidalTensor::getParticleN()%2==0);
        ASSERT(r_tidal_tensor>=0.0);
        ASSERT(id_offset>0);
        ASSERT(gravitational_constant>0);
        ASSERT(orbit_manager.checkParams());
        return true;
    }

    //! create artificial particles 
    /*! First TidalTensor::getParticleN() are tidal tensor particles; others-1 are orbitial sample particles; last is c.m.  
      id: 
      tt/orb: id_offset + abs(member->id)*(n_artificial-1)/2 + index/2;
      cm: - abs(_bin.id)

      status:
      Tidal tensors/orbital
      0: left member N
      1: right member N
      others: index+1 (_data_to_store>0.0)

      c.m.: n_members

      mass_backup:
      Tidial tensors/orbital: 0.0 (_data_to_store <=0.0)
      c.m.:  mass(cm)

      @param[out]    _ptcl_new: artificial particles that will be added
      @param[in]     _bin: binary tree root
      @param[in]     _data_to_store: array of data to be stored in the status of artificial particles
      @param[in]     _n_data: number of data
    */
    template <class Tptcl>
    void createArtificialParticles(Tptcl* _ptcl_artificial,
                                   COMM::BinaryTree<Tptcl,COMM::Binary> &_bin, 
                                   const PS::F64* _data_to_store,
                                   const PS::S32 _n_data) {
        ASSERT(checkParams());
        // set id and status.d except cm particle
        PS::S32 n_artificial = getArtificialParticleN();
        for (int i=0; i<n_artificial-1; i++) {
            Tptcl* pi = &_ptcl_artificial[i];
            Tptcl* binary_member_i = _bin.getMember(i%2);
            pi->id = id_offset + abs(binary_member_i->id)*n_artificial +i;
            auto& pi_artificial = pi->group_data.artificial;
            pi_artificial.setParticleTypeToArtificial(PS::F64(i+1));
#ifdef ARTIFICIAL_PARTICLE_DEBUG
            assert(pi_artificial.isArtificial());
#endif
        }

        // First TidalTensor::getParticleN() is used for tidal tensor points
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(r_tidal_tensor<=_bin.changeover.getRin());
#endif
        TidalTensor::createTidalTensorMeasureParticles(_ptcl_artificial, *((Tptcl*)&_bin), r_tidal_tensor);

        // remaining is for orbital sample particles
        orbit_manager.createSampleParticles(&(_ptcl_artificial[getIndexOffsetOrb()]), _bin);
        
        // store the component member number 
        for (int j=0; j<2; j++) {
            PS::S32 n_members = _bin.isMemberTree(j) ? ((COMM::BinaryTree<Tptcl,COMM::Binary>*)(_bin.getMember(j)))->getMemberN() : 1;
            _ptcl_artificial[j].group_data.artificial.storeData(n_members); 
#ifdef ARTIFICIAL_PARTICLE_DEBUG
            assert(n_members>0);
#endif
        }

        // store the additional data (should be positive) 
        // ensure the data size is not overflow
        assert(_n_data+2<n_artificial-1);
        for (int j=0; j<_n_data; j++) {
            _ptcl_artificial[j+2].group_data.artificial.storeData(_data_to_store[j]);
#ifdef ARTIFICIAL_PARTICLE_DEBUG
            assert(_data_to_store[j]>0);
#endif
        }

        // last member is the c.m. particle
        Tptcl* pcm;
        pcm = &_ptcl_artificial[getIndexOffsetCM()];
        pcm->group_data.artificial.setParticleTypeToCM(_bin.mass, _bin.getMemberN());
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(pcm->group_data.artificial.isCM());
#endif        
        if (getOrbitalParticleN()>0) pcm->mass = 0.0;
        else pcm->mass = _bin.mass;
        pcm->pos = _bin.pos;
        pcm->vel = _bin.vel;
        pcm->id  = - std::abs(_bin.id);
    }

    //! correct orbit-samping/pseudo particles force
    /*!
      replace c.m. force by the averaged force on sample/pseudo particles
      @param[in,out] _ptcl_artificial: one group of artificial particles 
    */
    template <class Tptcl>
    void correctOrbitalParticleForce(Tptcl* _ptcl_artificial) {
        auto* porb = getOrbitalParticles(_ptcl_artificial);
        const PS::S32 n_orb = getOrbitalParticleN();

        if (n_orb>0) {
            auto* pcm = getCMParticles(_ptcl_artificial);
            PS::F64vec& acc_cm = pcm->acc;
        
            acc_cm=PS::F64vec(0.0);
            PS::F64 m_ob_tot = 0.0;

            for (int k=0; k<n_orb; k++) {
                acc_cm += porb[k].mass*porb[k].acc; 
                m_ob_tot += porb[k].mass;
            }
            acc_cm /= m_ob_tot;

#ifdef ARTIFICIAL_PARTICLE_DEBUG
            assert(abs(m_ob_tot-pcm->group_data.artificial.getMassBackup())<1e-10);
#endif
        }
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

        // not consistent
        // After c.m. force used, it can be replaced by the averaged force on orbital particles
        // correctOrbitalParticleForce(_ptcl_artificial);
    }

    //! get oribit/pseudo particle list address from a artificial particle array
    template <class Tptcl>
    Tptcl* getOrbitalParticles(Tptcl* _ptcl_list)  {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(_ptcl_list[getIndexOffsetOrb()].group_data.artificial.isArtificial());
#endif
        return &_ptcl_list[getIndexOffsetOrb()];
    }

    //! get tidal tensor particle list address from a artificial particle array
    template <class Tptcl>
    Tptcl* getTidalTensorParticles(Tptcl* _ptcl_list) {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(_ptcl_list[getIndexOffsetTT()].group_data.artificial.isArtificial());
#endif
        return &_ptcl_list[getIndexOffsetTT()];
    }

    //! get c.m. particle list address from a artificial particle array
    template <class Tptcl>
    Tptcl* getCMParticles(Tptcl* _ptcl_list)  {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(_ptcl_list[getIndexOffsetCM()].group_data.artificial.isCM());
#endif
        return &_ptcl_list[getIndexOffsetCM()];
    }

    //! tidal tensor particle index offset
    PS::S32 getIndexOffsetTT() const {
        return 0;
    }

    //! orbital/pseudo particle index offset
    PS::S32 getIndexOffsetOrb() const {
        return TidalTensor::getParticleN();
    }

    //! CM particle index offset
    PS::S32 getIndexOffsetCM() const {
        return TidalTensor::getParticleN() + orbit_manager.getParticleN();
    }

    //! get artificial particle total number
    PS::S32 getArtificialParticleN() const {
        return TidalTensor::getParticleN() + orbit_manager.getParticleN() + 1;
    }

    //! get artificial particle total number
    PS::S32 getTidalTensorParticleN() const {
        return TidalTensor::getParticleN();
    }

    //! get orbitial particle number 
    PS::S32 getOrbitalParticleN() const {
        return orbit_manager.getParticleN();
    }

    //! get left member number
    template <class Tptcl>
    PS::S32 getLeftMemberN(const Tptcl* _ptcl_list) const {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(_ptcl_list[0].group_data.artificial.isArtificial());
#endif
        return PS::S32(_ptcl_list[0].group_data.artificial.getData());
    }

    //! get member number
    template <class Tptcl>
    PS::S32 getMemberN(const Tptcl* _ptcl_list) const {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(_ptcl_list[getIndexOffsetCM()].group_data.artificial.isCM());
#endif
        return PS::S32(_ptcl_list[getIndexOffsetCM()].group_data.artificial.getData(true));
    }

    //! get right member number
    template <class Tptcl>
    PS::S32 getRightMemberN(const Tptcl* _ptcl_list) const {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(_ptcl_list[1].group_data.artificial.isArtificial());
#endif
        return PS::S32(_ptcl_list[1].group_data.artificial.getData(true));
    }

    //! get center of mass id
    template <class Tptcl>
    PS::S64 getCMID(const Tptcl* _ptcl_list) const {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(_ptcl_list[getIndexOffsetCM()].group_data.artificial.isCM());
#endif
        return -PS::S64(_ptcl_list[getIndexOffsetCM()].id);
    }

    //! get stored data 
    template <class Tptcl>
    PS::F64 getStoredData(const Tptcl* _ptcl_list, const PS::S32 _index, const bool _is_positive) const {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(_ptcl_list[_index+2].group_data.artificial.isArtificial());
#endif
        return _ptcl_list[_index+2].group_data.artificial.getData(_is_positive);
    }

    //! write class data to file with binary format
    /*! @param[in] _fp: FILE type file for output
     */
    void writeBinary(FILE *_fp) {
        fwrite(&r_tidal_tensor, sizeof(PS::F64), 1, _fp);
        fwrite(&id_offset,      sizeof(PS::S64), 1, _fp);
        fwrite(&gravitational_constant, sizeof(PS::F64), 1, _fp);
        orbit_manager.writeBinary(_fp);
    }    

    //! read class data to file with binary format
    /*! @param[in] _fp: FILE type file for reading
     */
    void readBinary(FILE *_fin) {
        size_t rcount = 0;
        rcount += fread(&r_tidal_tensor, sizeof(PS::F64), 1, _fin);
        rcount += fread(&id_offset,      sizeof(PS::S64), 1, _fin);
        rcount += fread(&gravitational_constant, sizeof(PS::F64), 1, _fin);
        if (rcount<3) {
            std::cerr<<"Error: Data reading fails! requiring data number is 2, only obtain "<<rcount<<".\n";
            abort();
        }
        orbit_manager.readBinary(_fin);
    }    

    //! print parameters
    void print(std::ostream & _fout) const{
        _fout<<"r_tidal_tensor   : "<<r_tidal_tensor<<std::endl
             <<"id_offset        : "<<id_offset<<std::endl
             <<"G:               : "<<gravitational_constant<<std::endl;
        orbit_manager.print(_fout);
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
        assert(abs(mass_cm_check-pcm->group_data.artificial.getMassBackup())<1e-10);
        PS::F64vec dpos = pos_cm_check-pcm->pos;
        assert(abs(dpos*dpos)<1e-20);
    }
#endif

};
