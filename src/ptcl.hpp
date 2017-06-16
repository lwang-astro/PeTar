#pragma once
#include<particle_simulator.hpp>
#include"usr_define.hpp"

class Ptcl: public ParticleBase{
public:
    /*
                single           c.m.                       members            
      id         id          id of first member (-)            id
      states      0           number of members             fake pert id (-)
                 fake members                  unused
                      id                         -1
                 c.m. id+phase order(+)          -1
     */
    PS::S64 id;
    PS::S64 status;
    PS::F64 r_search;

    Ptcl(): id(-1), status(-1) {}

    template<class Tp>
    Ptcl(const Tp& p, const PS::S64 id_, const PS::S64 status_, const PS::F64 r_search_):
        id(id_), status(status_), r_search(r_search_) {ParticleBase::DataCopy(p);}

    template<class Tp>
    void DataCopy(const Tp& p) {
        ParticleBase::DataCopy(p);
        id       = p.id;
        status   = p.status;
        r_search = p.r_search;
    }

    void dump(std::ofstream & fout){
        ParticleBase::dump(fout);
        fout<<" id="<<id
            <<" status="<<status
            <<" r_search="<<r_search;
    }
};
