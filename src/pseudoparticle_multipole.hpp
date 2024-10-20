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
    template <class Tptcl>
    void createSampleParticles(Tptcl* _ptcl_artificial,
                               COMM::BinaryTree<Tptcl,COMM::Binary> &_bin) {
        PS::F64 m12 = _bin.mass;
        PS::F64 mu = _bin.m1*_bin.m2/m12;
        PS::F64 prefactor = _bin.semi*std::sqrt(mu/m12);
        PS::F64 pmass = m12/3.0;
        PS::F64 ecc2 = _bin.ecc*_bin.ecc;

        PS::F64 beta = prefactor*std::sqrt((1-ecc2)/4.0);
        PS::F64 alpha= prefactor*std::sqrt(3*ecc2+0.75);

        Tptcl* pi = &(_ptcl_artificial[0]);
        pi->mass = pmass;
        pi->pos = PS::F64vec(0, 2*beta, 0 );
        _bin.rotateToOriginalFrame(&(pi->pos.x));
        pi->pos += _bin.pos;
        pi->vel = _bin.vel;
#ifdef PETAR_USE_MPFRC
        pi->pos_mp = pi->pos;
#endif        
        
        pi = &(_ptcl_artificial[1]);
        pi->mass = pmass;
        pi->pos = PS::F64vec(alpha, -beta, 0 );
        _bin.rotateToOriginalFrame(&(pi->pos.x));
        pi->pos += _bin.pos;
        pi->vel = _bin.vel;
#ifdef PETAR_USE_MPFRC
        pi->pos_mp = pi->pos;
#endif        

        pi = &(_ptcl_artificial[2]);
        pi->mass = pmass;
        pi->pos = PS::F64vec(-alpha, -beta, 0 );
        _bin.rotateToOriginalFrame(&(pi->pos.x));
        pi->pos += _bin.pos;
        pi->vel = _bin.vel;
#ifdef PETAR_USE_MPFRC
        pi->pos_mp = pi->pos;
#endif        

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
    void writeBinary(FILE *_fp) {
    }

    //! read class data to file with binary format
    /*! @param[in] _fp: FILE type file for reading
     */
    void readBinary(FILE *_fin) {
    }    

    //! print parameters
    void print(std::ostream & _fout) const{
    }
};
