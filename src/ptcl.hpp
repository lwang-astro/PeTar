#pragma once
#include<particle_simulator.hpp>
#include"usr_define.hpp"

const PS::F64 SAFTY_FACTOR_FOR_SEARCH = 1.05;
const PS::F64 SAFTY_FACTOR_FOR_SEARCH_SQ = SAFTY_FACTOR_FOR_SEARCH * SAFTY_FACTOR_FOR_SEARCH;
const PS::F64 SAFTY_OFFSET_FOR_SEARCH = 1e-7;
//const PS::F64 SAFTY_OFFSET_FOR_SEARCH = 0.0;

class Ptcl: public ParticleBase{
public:
    /*
                single           c.m.                       members            
      id         id          id of first member (-)            id
      states      0           number of members             fake pert id (-)
      mass_bk   unknown         mass                         mass
                 fake members                  unused
                      id                         -1
                 c.m. id+phase order(+)          -1
                    unknown                    unknown
     */
    PS::F64 r_search;
    PS::F64 mass_bk;
    PS::S64 id;
    PS::S64 status;

    Ptcl(): id(-1), status(-1) {}

    template<class Tp>
    Ptcl(const Tp& p_): ParticleBase(p_), id(-1), status(-1) {}

    template<class Tp>
    Ptcl(const Tp& p_, const PS::F64 r_search_, const PS::F64 mass_bk_, const PS::S64 id_, const PS::S64 status_): ParticleBase(p_), r_search(r_search_), mass_bk(mass_bk_), id(id_), status(status_)  {}

    template<class Tp>
    void DataCopy(const Tp& p) {
        ParticleBase::DataCopy(p);
        r_search = p.r_search;
        mass_bk  = p.mass_bk;
        id       = p.id;
        status   = p.status;
    }

    void dump(std::ofstream & fout){
        ParticleBase::dump(fout);
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

    void readAscii(FILE* fp) {
        ParticleBase::readAscii(fp);
        PS::S64 rcount=fscanf(fp, "%lf %lf %lld %lld ",
                              &this->r_search, &this->mass_bk, &this->id, &this->status);
        if (rcount<4) {
            std::cerr<<"Error: Data reading fails! requiring data number is 8, only obtain "<<rcount<<".\n";
            abort();
        }
    }

};
