#pragma once

#include "ptcl.hpp"

class PtclHard: public Ptcl{
public:
    PS::S32 id_cluster;
    PS::S32 adr_org;
    static PS::F64 r_group_over_in;
    static PS::F64 r_search_group_over_in;

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

    //! Get group candidate distance criterion
    PS::F64 getRGroupCandidate() const {
        return changeover.getRin()*r_search_group_over_in;
    }

    //! Get group distance criterion
    PS::F64 getRGroup() const {
#ifdef HARD_DEBUG
        assert(r_group_over_in>0);
#endif
        return changeover.getRin()*r_group_over_in;
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

    //! Get neighbor distance criterion 
    PS::F64 getRNeighbor() const {
#ifdef HARD_DEBUG
        // If a binary's velocity is zero. its r_search can be the same as r_out, because r_out can be > r_search_min, then calcRSearch return r_out.
        assert(r_search>=changeover.getRout());
#endif 
        return r_search;
    }

};
