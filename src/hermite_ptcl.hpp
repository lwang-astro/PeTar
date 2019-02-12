#pragma once

#include "ptcl.hpp"

//! Particle type for hermite integrator
class PtclH4: public Ptcl{
public:
    PS::F64vec acc0;
    PS::F64vec acc1;
    PS::F64 dt;
    PS::F64 time;
#ifdef HARD_DEBUG_ACC
    PS::F64vec acc2; // for debug
    PS::F64vec acc3; // for debug
#endif

    PtclH4(): dt(0.0), time(0.0) {}
    
    template<class Tptcl>
    PtclH4(const Tptcl &_p): Ptcl(_p), dt(0.0), time(0.0) {}

    void dump(FILE *fp) {
        fwrite(this, sizeof(PtclH4),1,fp);
    }

    void read(FILE *fp) {
        size_t rcount = fread(this, sizeof(PtclH4),1,fp);
        if (rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            abort();
        }
    }
};

//! Particle predictor type for hermite integrator
class PtclPred{
public:
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64 mass;
    PS::F64 r_search;
};

//! Particle force for hermite integrator
class PtclForce{
public:
    PS::F64vec acc0; //
    PS::F64vec acc1; //
    void clear(){
        acc0 = acc1 = 0.0;
    }

    void dump(FILE *fp) {
        fwrite(this, sizeof(PtclForce),1,fp);
    }

    void read(FILE *fp) {
        size_t rcount = fread(this, sizeof(PtclForce),1,fp);
        if (rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            abort();
        }
    }
};
