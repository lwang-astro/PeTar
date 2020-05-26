#pragma once
#include<particle_simulator.hpp>

//! escaper checker
class Escaper{
public:
    PS::F64 r_escape; // escape distance criterion

    Escaper(): r_escape(PS::LARGE_FLOAT) {}

    //! check escaper based on distance and c.m. of the system
    template <class Tptcl, class Tpcm>
    bool isEscaper(Tptcl& _p, Tpcm& _pcm) {
        //assert(_p.pot_tot!=0.0);
        PS::F64 kin = _p.vel*_p.vel;
        PS::F64 tot = _p.pot_tot + kin;
        PS::F64vec dr = _p.pos - _pcm.pos;
        PS::F64 r2 = dr*dr;
        if (tot>0&&r2>r_escape*r_escape) return true;
        else return false;
    }
};
