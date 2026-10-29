#pragma once
#include "Common/binary_tree.h"
#include <iostream>

//! PseudoParticle Multipole method for counter force from binaries
class PseudoParticleMultipoleManager{
public:
    //! check paramters
    bool checkParams() {
        return true;
    }

    //! create pseudoparticle multipole particles
    /*! Create three particles that have the same quadrupole moment as orbit-averaged binaries.
        @param[in] _ptcl_artificial: particle array to store the sample particles, 2*n_split_ will be used
        @param[in] _bin: binary orbit 
     */
    template <class Tptcl, class Tbin>
    void createSampleParticles(Tptcl* _ptcl_artificial,
                               Tbin &_bin) {
        PS::F64 m12 = _bin.mass;
        PS::F64 mu = _bin.m1*_bin.m2/m12;
        PS::F64 prefactor = _bin.semi*std::sqrt(mu/m12);
        PS::F64 pmass = m12/3.0;
        PS::F64 ecc2 = _bin.ecc*_bin.ecc;

        PS::F64 beta = prefactor*std::sqrt((1-ecc2)/4.0);
        PS::F64 alpha= prefactor*std::sqrt(3*ecc2+0.75);

        // Canonical orientation (2026-10-29): particleToOrbit stores the orbit
        // under the (member1, member2) convention of the binding chain, so when
        // the member ordering swaps at a group reformation, rot_self shifts by
        // pi and this triangle point-reflects: the quadrupole (the design
        // target) is invariant, but the spurious odd multipoles of the 3-point
        // discretization flip sign and every external particle receives a
        // discontinuous force jump (~(a/d)^3 relative, global). Fix: define the
        // in-plane +x axis to point toward the orbit-averaged position of the
        // HEAVIER member (a physical invariant: <r_heavy> = +(m_light/M) *
        // (3/2) a e * e_hat of the member1-convention orbit). This is enforced
        // by point-reflecting the triangle (equivalent to rot_self += pi) when
        // member1 is the lighter one; masses alone decide for m1 != m2, and
        // the member IDs break the exact tie so the rule is a total canonical
        // ordering. The quadrupole is unchanged under the reflection.
        // Residual corner: node IDs are set only on the AR-info/Hermite
        // refresh trees (checkAndSetBinaryPairID), not on the stable-tree
        // nodes at group formation, so an exactly-equal-AGGREGATE-mass
        // hierarchical pair falls through the tie break there (id 0 == 0)
        // and may point-reflect once at the first refresh.
        PS::F64 flip = 1.0;
        if (_bin.m1 < _bin.m2) flip = -1.0;
        else if (_bin.m1 == _bin.m2) {
            auto* mem0 = _bin.getMember(0);
            auto* mem1 = _bin.getMember(1);
            if (mem0 != NULL && mem1 != NULL && mem0->id > mem1->id) flip = -1.0;
        }

        Tptcl* pi = &(_ptcl_artificial[0]);
        pi->mass = pmass;
        pi->pos = PS::F64vec(flip*0, flip*2*beta, 0 );
        _bin.rotateToOriginalFrame(&(pi->pos.x));
        pi->pos += _bin.pos;
        pi->vel = _bin.vel;
        
        pi = &(_ptcl_artificial[1]);
        pi->mass = pmass;
        pi->pos = PS::F64vec(flip*alpha, flip*(-beta), 0 );
        _bin.rotateToOriginalFrame(&(pi->pos.x));
        pi->pos += _bin.pos;
        pi->vel = _bin.vel;

        pi = &(_ptcl_artificial[2]);
        pi->mass = pmass;
        pi->pos = PS::F64vec(flip*(-alpha), flip*(-beta), 0 );
        _bin.rotateToOriginalFrame(&(pi->pos.x));
        pi->pos += _bin.pos;
        pi->vel = _bin.vel;

#ifdef ARTIFICIAL_PARTICLE_DEBUG
        PS::F64vec dv = (_ptcl_artificial[0].pos + _ptcl_artificial[1].pos + _ptcl_artificial[2].pos)/3 - _bin.pos;
        assert(dv*dv<1e-10);
#endif
    }

    //! get particle number 
    static PS::S32 getParticleN() {
        return 3;
    }

    //! write class data to file with binary format
    /*! @param[in] _fp: FILE type file for output
     */
    void writeBinary(FILE *_fp) const {
    }

    void printColumnBinary(std::ostream& _fout) const {
    }

    //! read class data to file with binary format
    /*! @param[in] _fp: FILE type file for reading
     */
    void readBinary(FILE *_fin) {
    }    

    void readBinary(std::istream& _fin) {
    }

    //! print parameters
    void print(std::ostream & _fout) const{
    }
};
