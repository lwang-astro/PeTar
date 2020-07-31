#pragma once
#include<particle_simulator.hpp>
#include"particle_base.hpp"
#include"changeover.hpp"
#define ASSERT assert
#include"artificial_particles.hpp"

const PS::F64 SAFTY_FACTOR_FOR_SEARCH = 0.99;
//const PS::F64 SAFTY_FACTOR_FOR_SEARCH_SQ = SAFTY_FACTOR_FOR_SEARCH * SAFTY_FACTOR_FOR_SEARCH;
//const PS::F64 SAFTY_OFFSET_FOR_SEARCH = 1e-7;
//const PS::F64 SAFTY_OFFSET_FOR_SEARCH = 0.0;

//! group data delivery, used for two purpose
/*! artificial is used to store information of artificial particles
    cm is used to store the c.m. particle mass and velocity
 */
union GroupDataDeliver{
    ArtificialParticleInformation artificial;
    struct {PS::F32 mass; PS::F32vec vel;} cm;

    GroupDataDeliver(): artificial() {}
    
    GroupDataDeliver(const GroupDataDeliver& _data): artificial(_data.artificial) {}

    GroupDataDeliver& operator= (const GroupDataDeliver& _data) {
        artificial = _data.artificial;
        return *this;
    }
};

//! group data deliver mode
enum class GroupDataMode{artificial=1, cm=2, none=0};

//! Particle class 
class Ptcl: public ParticleBase{
public:
    PS::F64 r_search;
    PS::S64 id; // positive for single, artificial but not cm, 0 for unused particle
    GroupDataDeliver group_data;
    ChangeOver changeover;
    static PS::F64 search_factor;
    static PS::F64 r_search_min;
    static PS::F64 r_group_crit_ratio;
    static PS::F64 mean_mass_inv;
    static PS::F64vec vel_cm;
    static GroupDataMode group_data_mode;

    Ptcl(): ParticleBase(), r_search(-PS::LARGE_FLOAT), id(0), group_data(), changeover() {}

    template<class Tptcl>
    Ptcl(const Tptcl& _p) { Ptcl::DataCopy(_p);  }

    template<class Tptcl>
    Ptcl(const Tptcl& _p, const PS::F64 _r_search, const PS::S64 _id, const GroupDataDeliver& _group_data, const ChangeOver& _co): ParticleBase(_p), r_search(_r_search), id(_id), group_data(_group_data), changeover(_co) {}

    template<class Tptcl>
    void DataCopy(const Tptcl& _p) {
        ParticleBase::DataCopy(_p);
        r_search   = _p.r_search;
        id         = _p.id;
        group_data = _p.group_data;
        changeover = _p.changeover;
    }

    template<class Tptcl>
    Ptcl& operator = (const Tptcl& _p) {
        Ptcl::DataCopy(_p);
        return *this;
    }

    void print(std::ostream & _fout) const{
        ParticleBase::print(_fout);
        _fout<<" r_search="<<r_search
             <<" id="<<id;
        group_data.artificial.print(_fout);
        changeover.print(_fout);
    }

    //! print titles of class members using column style
    /*! print titles of class members in one line for column style
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
     */
    static void printColumnTitle(std::ostream & _fout, const int _width=20) {
        ParticleBase::printColumnTitle(_fout, _width);
        _fout<<std::setw(_width)<<"r_search"
             <<std::setw(_width)<<"id";
        ArtificialParticleInformation::printColumnTitle(_fout, _width);
        ChangeOver::printColumnTitle(_fout, _width);
    }

    //! print column title with meaning (each line for one column)
    /*! @param[out] _fout: std::ostream output object
      @param[in] _counter: offset of the number counter for each line to indicate the column index (defaulted 0)
      @param[in] _offset: the printing whitespace offset for each line (defaulted 0)
      \return: the total counter of columns
     */
    static int printTitleWithMeaning(std::ostream & _fout, const int _counter=0, const int _offset=0) {
        int counter = _counter;
        counter = ParticleBase::printTitleWithMeaning(_fout, counter, _offset);
        counter++;
        _fout<<std::setw(_offset)<<" "<<counter<<". r_search: neighbor searching radius (0.0)\n";
        counter++;
        _fout<<std::setw(_offset)<<" "<<counter<<". id: identification of particle, should be a positive unique value (>0)\n";
        counter = ArtificialParticleInformation::printTitleWithMeaning(_fout, counter, _offset);
        counter = ChangeOver::printTitleWithMeaning(_fout, counter, _offset);
        return counter;
    }

    //! print data of class members using column style
    /*! print data of class members in one line for column style. Notice no newline is printed at the end
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
     */
    void printColumn(std::ostream & _fout, const int _width=20){
        ParticleBase::printColumn(_fout, _width);
        _fout<<std::setw(_width)<<r_search
             <<std::setw(_width)<<id;
        group_data.artificial.printColumn(_fout, _width);
        changeover.printColumn(_fout, _width);
    }

    //! write class data with ASCII format
    /*! @param[in] _fout: file IO for write
     */
    void writeAscii(FILE* _fout) const{
        ParticleBase::writeAscii(_fout);
        fprintf(_fout, "%26.17e %lld ", 
                this->r_search, this->id);
        group_data.artificial.writeAscii(_fout);
        changeover.writeAscii(_fout);
    }

    //! write class data with BINARY format
    /*! @param[in] _fout: file IO for write
     */
    void writeBinary(FILE* _fout) const{
        ParticleBase::writeBinary(_fout);
        fwrite(&(this->r_search), sizeof(PS::F64), 2, _fout);
        group_data.artificial.writeBinary(_fout);
        changeover.writeBinary(_fout);
    }

    //! read class data with ASCII format
    /*! @param[in] _fin: file IO for read
     */
    void readAscii(FILE* _fin) {
        ParticleBase::readAscii(_fin);
        PS::S64 rcount=fscanf(_fin, "%lf %lld ",
                              &this->r_search, &this->id);
        if (rcount<2) {
            std::cerr<<"Error: Data reading fails! requiring data number is 2, only obtain "<<rcount<<".\n";
            abort();
        }
        group_data.artificial.readAscii(_fin);
        changeover.readAscii(_fin);
    }

    //! read class data with BINARY format
    /*! @param[in] _fin: file IO for read
     */
    void readBinary(FILE* _fin) {
        ParticleBase::readBinary(_fin);
        size_t rcount = fread(&(this->r_search), sizeof(PS::F64), 2, _fin);
        if (rcount<2) {
            std::cerr<<"Error: Data reading fails! requiring data number is 4, only obtain "<<rcount<<".\n";
            abort();
        }
        group_data.artificial.readBinary(_fin);
        changeover.readBinary(_fin);
    }

    //! set status to c.m. particle address 
    void setParticleCMAddress(const PS::S64 _adr) {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(group_data.artificial.isMember()&&_adr>0);
#endif
        group_data.artificial.setStatus(-_adr);
    }

    //! get c.m. particle address
    PS::S64 getParticleCMAddress() const {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(group_data.artificial.isMember());
#endif
        PS::F64 status = group_data.artificial.getStatus();
        if (status==-PS::LARGE_FLOAT) return -1;
        else return PS::S64(-status);
    }

    //! calculate new rsearch
    /*! calculate r_search based on velocity and tree step 
     */
    void calcRSearch(const PS::F64 _dt_tree) {
        PS::F64vec dv = vel- vel_cm;
        PS::F64 v = std::sqrt(dv*dv);
        r_search = std::max(v*_dt_tree*search_factor+changeover.getRout(), r_search_min);
        //r_search = std::max(std::sqrt(vel*vel)*dt_tree*search_factor, std::sqrt(mass*mean_mass_inv)*r_search_min);
#ifdef HARD_DEBUG
        assert(r_search>0);
        assert(r_search_min>0);
        assert(search_factor>0);
#endif
    }

    //! Get neighbor distance criterion 
    PS::F64 getRNeighbor() const {
#ifdef HARD_DEBUG
        assert(r_search>changeover.getRout());
#endif 
        return r_search;
    }

    //! Get group candidate distance criterion
    PS::F64 getRGroupCandidate() const {
        return changeover.getRin();
    }

    //! Get group distance criterion
    PS::F64 getRGroup() const {
#ifdef HARD_DEBUG
        assert(r_group_crit_ratio>0);
#endif
        return changeover.getRin()*r_group_crit_ratio;
    }

    
};

#undef ASSERT
