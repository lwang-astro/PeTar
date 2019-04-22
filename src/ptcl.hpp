#pragma once
#include<particle_simulator.hpp>
#include"usr_define.hpp"
#include"changeover.hpp"

const PS::F64 SAFTY_FACTOR_FOR_SEARCH = 0.99;
//const PS::F64 SAFTY_FACTOR_FOR_SEARCH_SQ = SAFTY_FACTOR_FOR_SEARCH * SAFTY_FACTOR_FOR_SEARCH;
//const PS::F64 SAFTY_OFFSET_FOR_SEARCH = 1e-7;
//const PS::F64 SAFTY_OFFSET_FOR_SEARCH = 0.0;

//! Particle class 
class Ptcl: public ParticleBase{
public:
    /*
                single           c.m.                       members               unused        suppressed c.m.
      id         id          id of first member (-)            id                   -1          id of previous c.m. (-)
      status      0          member number                  c.m. adr (-)            -1            -20 
      mass_bk     0             mass                         mass                 unknown       unknown
                 fake members                                                                            
                id_offset+id*n_split+iphase                                             
                1. first component member number 2. second. 3. i_cluster+1, 4. i_group+1, others: (c.m.id<<ID_PHASE_SHIFT)|i
                  binary parameters                                                 

      PS: mass_bk is used to store perturber force in searchpart
          suppressed c.m. is set in HermiteIntegrator.removePtclList
     */
    PS::F64 r_search;
    PS::F64 mass_bk;
    PS::S64 id;
    PS::S64 status;
    ChangeOver changeover;
    static PS::F64 search_factor;
    static PS::F64 r_search_min;
    static PS::F64 r_group_crit_ratio;
    static PS::F64 mean_mass_inv;

    Ptcl(): r_search(-1.0), id(-10), status(-10), changeover() {}

    template<class Tptcl>
    Ptcl(const Tptcl& _p) { Ptcl::DataCopy(_p);  }

    template<class Tptcl>
    Ptcl(const Tptcl& _p, const PS::F64 _r_search, const PS::F64 _mass_bk, const PS::S64 _id, const PS::S64 _status, const ChangeOver& _co): ParticleBase(_p), r_search(_r_search), mass_bk(_mass_bk), id(_id), status(_status), changeover(_co) {}

    template<class Tptcl>
    void DataCopy(const Tptcl& _p) {
        ParticleBase::DataCopy(_p);
        r_search = _p.r_search;
        mass_bk  = _p.mass_bk;
        id       = _p.id;
        status   = _p.status;
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
             <<" mass_bk="<<mass_bk
             <<" id="<<id
             <<" status="<<status;
        changeover.print(_fout);
    }

    //! print titles of class members using column style
    /*! print titles of class members in one line for column style
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
     */
    void printColumnTitle(std::ostream & _fout, const int _width=20) {
        ParticleBase::printColumnTitle(_fout, _width);
        _fout<<std::setw(_width)<<"r_search"
             <<std::setw(_width)<<"mass_bk"
             <<std::setw(_width)<<"id"
             <<std::setw(_width)<<"status";
        ChangeOver::printColumnTitle(_fout, _width);
    }

    //! print data of class members using column style
    /*! print data of class members in one line for column style. Notice no newline is printed at the end
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
     */
    void printColumn(std::ostream & _fout, const int _width=20){
        ParticleBase::printColumn(_fout, _width);
        _fout<<std::setw(_width)<<r_search
             <<std::setw(_width)<<mass_bk
             <<std::setw(_width)<<id
             <<std::setw(_width)<<status;
        changeover.printColumn(_fout, _width);
    }

    void writeAscii(FILE* _fout) const{
        ParticleBase::writeAscii(_fout);
        fprintf(_fout, "%26.17e %26.17e %lld %lld ", 
                this->r_search, this->mass_bk, this->id, this->status);
        changeover.writeAscii(_fout);
    }

    void writeBinary(FILE* _fin) const{
        ParticleBase::writeBinary(_fin);
        fwrite(&(this->r_search), sizeof(PS::F64), 4, _fin);
        changeover.writeBinary(_fin);
    }

    void readAscii(FILE* _fin) {
        ParticleBase::readAscii(_fin);
        PS::S64 rcount=fscanf(_fin, "%lf %lf %lld %lld ",
                              &this->r_search, &this->mass_bk, &this->id, &this->status);
        if (rcount<4) {
            std::cerr<<"Error: Data reading fails! requiring data number is 4, only obtain "<<rcount<<".\n";
            abort();
        }
        changeover.readAscii(_fin);
    }

    void readBinary(FILE* _fin) {
        ParticleBase::readBinary(_fin);
        size_t rcount = fread(&(this->r_search), sizeof(PS::F64), 4, _fin);
        if (rcount<4) {
            std::cerr<<"Error: Data reading fails! requiring data number is 4, only obtain "<<rcount<<".\n";
            abort();
        }
        changeover.readBinary(_fin);
    }

    //! calculate new rsearch
    /*! calculate r_search based on velocity and tree step 
     */
    void calcRSearch(const PS::F64 _dt_tree) {
        PS::F64 v = std::sqrt(vel*vel);
        r_search = std::max(v*_dt_tree*search_factor+changeover.getRout(), r_search_min);
        //r_search = std::max(std::sqrt(vel*vel)*dt_tree*search_factor, std::sqrt(mass*mean_mass_inv)*r_search_min);
#ifdef HARD_DEBUG
        assert(r_search>0);
#endif
    }

    //! Get neighbor distance criterion 
    PS::F64 getRNeighbor() const {
        return r_search;
#ifdef HARD_DEBUG
        assert(r_search>0);
#endif 
   }

    //! Get neighbor distance criterion 
    PS::F64 getRBreak() const {
        return changeover.getRin()*r_group_crit_ratio;
    }

//    PS::F64 calcRSearch(const PS::F64 _dt_tree, const PS::F64 _v_max) {
//        PS::F64 v = std::sqrt(vel*vel);
//        PS::F64 dt_reduce_factor = 1.0;
//        if (v>_v_max) {
//            dt_reduce_factor = v/_v_max;
//            v = _v_max;
//        }
//        r_search = std::max(v*_dt_tree*search_factor, r_search_min);
//        //r_search = std::max(std::sqrt(vel*vel)*dt_tree*search_factor, std::sqrt(mass*mean_mass_inv)*r_search_min);
//#ifdef HARD_DEBUG
//        assert(r_search>0);
//#endif
//        return dt_reduce_factor;
//    }

};

PS::F64 Ptcl::r_search_min = 0.0;
PS::F64 Ptcl::search_factor= 0.0;
PS::F64 Ptcl::mean_mass_inv= 0.0; // mean mass inverse
PS::F64 Ptcl::r_group_crit_ratio =0.0;
