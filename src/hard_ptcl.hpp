#pragma once

#include "ptcl.hpp"

class PtclHard: public Ptcl{
public:
    PS::S32 id_cluster;
    PS::S32 adr_org;

    PtclHard(): Ptcl(), id_cluster(-1), adr_org(-1) {}

    template<class Tptcl>
    PtclHard(const Tptcl& _p, const PS::F64 _r_search, const PS::F64 _mass_bk, const PS::S64 _id, const PS::S64 _status, const ChangeOver& _co, const PS::S32 _id_cluster, const PS::S32 _adr_org): 
        Ptcl(_p, _r_search, _mass_bk, _id, _status, _co), id_cluster(_id_cluster), adr_org(_adr_org) {}

    template<class Tptcl>
    PtclHard(const Tptcl &_p, const PS::S32 _id_cluster, const PS::S32 _adr_org): 
        Ptcl(_p), id_cluster(_id_cluster), adr_org(_adr_org) {}

    template<class Tptcl>
    PtclHard(const Tptcl &_p) {
        Ptcl::DataCopy(_p);
        id_cluster = _p.id_cluster;
        adr_org = _p.adr_org;
    }

    // notice datacopy ignore id_cluster and adr_org
    template<class Tptcl>
    void DataCopy(const Tptcl& _p) {
        Ptcl::DataCopy(_p);
    }

    template<class Tptcl>
    PtclHard& operator = (const Tptcl& _p) {
        Ptcl::DataCopy(_p);
        id_cluster = _p.id_cluster;
        adr_org = _p.adr_org;
        return *this;
    }

    void setTidalTensorID(const PS::S32 _id) {
#ifdef HARD_DEBUG
        assert(_id>0);
#endif
        group_data.artificial.setStatus(PS::F64(-_id));
    }

    PS::S32 getTidalTensorID() const {
        return PS::S32(- group_data.artificial.getStatus());
    }

    // calculate reseach based on potential and velocity
    /*!
      Potential criterion: 
      pot = -Gm/r0
      -Gm/r0 + 1/2 v^2 = -Gm/(r0+dr)
      dr = r0/[2Gm/(r0 v^2) - 1]
      If energy is positive, dr < 0

      velocity criterion:
      v*dt_tree

      Use min of two
     */
    void calcRSearch(const PS::F64 _Gm, const PS::F64 _pot, const PS::F64vec& _vel_cm, const PS::F64 _dt_tree) {
        PS::F64vec vrel = Ptcl::vel-_vel_cm;
        PS::F64 v2rel = vrel*vrel;
        PS::F64 r0 = _pot>0? _Gm/_pot : PS::LARGE_FLOAT;
        PS::F64 q = 2.0*_pot/v2rel;
        if (q>1.0) r0 /= (q-1.0);

        PS::F64 v2 = Ptcl::vel*Ptcl::vel;
        PS::F64 v = std::sqrt(v2);
        r_search = std::max(std::min(v*_dt_tree, r0)*search_factor+Ptcl::changeover.getRout(), r_search_min);
    }

    void dump(FILE *_fout) {
        fwrite(this, sizeof(*this),1,_fout);
    }

    void read(FILE *_fin) {
        size_t rcount = fread(this, sizeof(*this),1,_fin);
        if (rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            abort();
        }
    }
    void print(std::ostream & _fout) const{
        Ptcl::print(_fout);
        std::cerr<<" id_cluster="<<id_cluster
                 <<" adr_org="<<adr_org;
    }

};
