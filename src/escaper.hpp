#pragma once
#include<particle_simulator.hpp>

//! escaper checker
class Escaper{
public:
    PS::F64 r_escape_sq; // escape distance criterion square
    bool check_energy_flag; // if true, escaper should have E>0

    Escaper(): r_escape_sq(PS::LARGE_FLOAT), check_energy_flag(true) {}

    //! check escaper based on distance and c.m. of the system
    template <class Tptcl, class Tpcm>
    bool isEscaper(Tptcl& _p, Tpcm& _pcm) {
        //assert(_p.pot_tot!=0.0);
        bool is_escape=true;

        if (check_energy_flag) { 
            // check energy
            PS::F64 kin = _p.vel*_p.vel;
            PS::F64 tot = _p.pot_tot + kin;
            if (tot<0) is_escape=false;
        }

        // check distance
#ifdef RECORD_CM_IN_HEADER
        PS::F64vec dr = _p.pos;
#else
        PS::F64vec dr = _p.pos - _pcm.pos;
#endif
        PS::F64 r2 = dr*dr;
        if (r2<r_escape_sq) is_escape=false;

	return is_escape;
    }
};
