#pragma once
#include<particle_simulator.hpp>
#include"usr_define.hpp"

const PS::F64 SAFTY_FACTOR_FOR_SEARCH = 0.99;
//const PS::F64 SAFTY_FACTOR_FOR_SEARCH_SQ = SAFTY_FACTOR_FOR_SEARCH * SAFTY_FACTOR_FOR_SEARCH;
//const PS::F64 SAFTY_OFFSET_FOR_SEARCH = 1e-7;
//const PS::F64 SAFTY_OFFSET_FOR_SEARCH = 0.0;

class Ptcl: public ParticleBase{
public:
    /*
                single           c.m.                       members               unused
      id         id          id of first member (-)            id                   -1
      status      0          member number                  c.m. adr (-)            -1
      mass_bk     0             mass                         mass                 unknown
                 fake members                                                                            
                id_offset+id*n_split+iphase                                             
                1. first component member number 2. second. 3. i_cluster+1, 4. i_group+1, others: (c.m.id<<ID_PHASE_SHIFT)|i
                  binary parameters                                                 
     */
    PS::F64 r_search;
    PS::F64 mass_bk;
    PS::S64 id;
    PS::S64 status;
    static PS::F64 search_factor;
    static PS::F64 r_search_min;
    static PS::F64 mean_mass_inv;

    Ptcl(): id(-10), status(-10) {}

    template<class Tptcl>
    Ptcl(const Tptcl& _p) { Ptcl::DataCopy(_p);  }

    template<class Tptcl>
    Ptcl(const Tptcl& _p, const PS::F64 _r_search, const PS::F64 _mass_bk, const PS::S64 _id, const PS::S64 _status): ParticleBase(_p), r_search(_r_search), mass_bk(_mass_bk), id(_id), status(_status)  {}

    template<class Tptcl>
    void DataCopy(const Tptcl& _p) {
        ParticleBase::DataCopy(_p);
        r_search = _p.r_search;
        mass_bk  = _p.mass_bk;
        id       = _p.id;
        status   = _p.status;
    }

    template<class Tptcl>
    Ptcl& operator = (const Tptcl& _p) {
        Ptcl::DataCopy(_p);
        return *this;
    }

    void print(std::ostream & fout){
        ParticleBase::print(fout);
        fout<<" r_search="<<r_search
            <<" mass_bk="<<mass_bk
            <<" id="<<id
            <<" status="<<status;
    }

    void writeAscii(FILE* fp) const{
        ParticleBase::writeAscii(fp);
        fprintf(fp, "%26.17e %26.17e %lld %lld ", 
                this->r_search, this->mass_bk, this->id, this->status);
    }

    void writeBinary(FILE* fp) const{
        ParticleBase::writeBinary(fp);
        fwrite(&(this->r_search), sizeof(PS::F64), 4, fp);
    }

    void readAscii(FILE* fp) {
        ParticleBase::readAscii(fp);
        PS::S64 rcount=fscanf(fp, "%lf %lf %lld %lld ",
                              &this->r_search, &this->mass_bk, &this->id, &this->status);
        if (rcount<4) {
            std::cerr<<"Error: Data reading fails! requiring data number is 4, only obtain "<<rcount<<".\n";
            abort();
        }
    }

    void readBinary(FILE* fp) {
        ParticleBase::readBinary(fp);
        size_t rcount = fread(&(this->r_search), sizeof(PS::F64), 4, fp);
        if (rcount<4) {
            std::cerr<<"Error: Data reading fails! requiring data number is 4, only obtain "<<rcount<<".\n";
            abort();
        }
    }

    void calcRSearch(const PS::F64 dt_tree) {
        r_search = std::max(std::sqrt(vel*vel)*dt_tree*search_factor, r_search_min);
        //r_search = std::max(std::sqrt(vel*vel)*dt_tree*search_factor, std::sqrt(mass*mean_mass_inv)*r_search_min);
#ifdef HARD_DEBUG
        assert(r_search>0);
#endif
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
};

PS::F64 Ptcl::r_search_min = 0.0;
PS::F64 Ptcl::search_factor= 0.0;
PS::F64 Ptcl::mean_mass_inv= 0.0; // mean mass inverse
