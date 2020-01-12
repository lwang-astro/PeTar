#pragma once

#include "Common/list.h"
#include "Hermite/hermite_particle.h"
#include "hard_ptcl.hpp"
#include "Common/binary_tree.h"
#include "tidal_tensor.hpp"


//! Perturber class for AR integration
class ARPerturber: public H4::Neighbor<PtclHard>{
public:
    typedef H4::Neighbor<PtclHard> NB;
    TidalTensor* soft_pert;  ///> soft perturbation 
    Float soft_pert_min; ///> minimum soft perturbation

    ARPerturber(): NB(), soft_pert(NULL), soft_pert_min(Float(0.0)) {}

    //! clear function
    void clear() {
        NB::clear();
        if (soft_pert!=NULL) {
            ASSERT(soft_pert->group_id>0.0);
            soft_pert->group_id = -soft_pert->group_id;
            soft_pert = NULL;
        }
    }

    //! check parameters status
    bool checkParams() {
        ASSERT(NB::checkParams());
        //ASSERT(soft_pert_min>=0.0);
        //ASSERT(soft_pert!=NULL);
        return true;
    }

    //! find close tidal tensor and if (-) tensor group id is the same as input, initial tidal tensor c.m.
    /*! if the tidal tensor is already in used (group_id>=0), copy a new one after _n_tt
      @param[in,out] _tt: tensor array
      @param[in,out] _n_tt: number of current tensor
      @param[in] _n_max: maximum size of tensor array
      @param[in] _cm: c.m. particle
      @param[in] _gid: group id (not necessary integer)
      \return the tidal tensor index, if no match, return -1
     */
    PS::S32 findCloseSoftPert(TidalTensor* _tt, int& _n_tt, const int _n_max, const H4::ParticleH4<PtclHard>& _cm, const PS::F64 _gid) {
        ASSERT(_gid>0.0);
        const PS::F64vec& pos = _cm.pos;
        PS::F64 r_min2=NUMERIC_FLOAT_MAX;
        PS::S32 r_min_index=-1;
        for (int i=0; i<_n_tt; i++) {
            PS::F64vec dr = pos - _tt[i].pos;
            PS::F64 r2 = dr*dr;
            if (r_min2>r2) {
                r_min_index = i;
                r_min2 = r2;
            }
        }
        // if the close tt is already used, copy a new one
        //if (_tt[r_min_index].group_id>=0) {
        //ASSERT(_n_tt<_n_max);
            
            //_tt[_n_tt] = _tt[r_min_index];
            //soft_pert = &_tt[_n_tt];
            //_n_tt++;
        //}
        if (-_tt[r_min_index].group_id==_gid) {
            soft_pert = &_tt[r_min_index];
            soft_pert->group_id = _gid;
            // update c.m.
            soft_pert->shiftCM(pos);
            return r_min_index;
        }
        else 
            return -1;
    }

    //! calculate soft perturbation 
    /*! calculate soft perturbation for two members of one binary
      @param[in] _p1: member 1 
      @param[in] _p2: member 2
      \return soft perturbation for slowdown
     */
    template <class Tptcl>
    Float calcSoftPertSlowDownBinary(const Tptcl& _p1, const Tptcl& _p2) {
        Float pert = 0.0;
#ifdef SOFT_PERT
        if(soft_pert!=NULL) {
            Float acc_p1[3] = {0.0, 0.0, 0.0};
            Float acc_p2[3] = {0.0, 0.0, 0.0};
            soft_pert->eval(acc_p1, _p1.pos);
            soft_pert->eval(acc_p2, _p2.pos);
            Float dacc[3] = {acc_p1[0]-acc_p2[0], 
                             acc_p1[1]-acc_p2[1],
                             acc_p1[2]-acc_p2[2]};
            pert = dacc[0]*dacc[0] + dacc[1]*dacc[1] + dacc[2]*dacc[2];
        }
#endif
        return pert;
    }

    //! calculate soft_pert_min
    template <class Tptcl>
    void calcSoftPertMin(const COMM::BinaryTree<Tptcl>& _bin, const Float _G) {
        // hyperbolic case
        if(_bin.semi<0.0) soft_pert_min = 0.0;
        else { // close orbit
            ParticleBase p[2];
            _bin.calcParticlesEcca(p[0], p[1], COMM::PI, _G);
            Float dacc_soft = calcSoftPertSlowDownBinary(p[0], p[1]);
            //soft_pert_min = _bin.mass*dacc_soft/(2.0*abs(_bin.semi));
            Float apo = _bin.semi*(1.0+_bin.ecc);
            soft_pert_min = _bin.mass*dacc_soft/(_G*apo*apo);
            //soft_pert_min = _bin.mass*_bin.mass*dacc_soft;
        }
    }
};
