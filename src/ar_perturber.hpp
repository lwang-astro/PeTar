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
#ifdef EXTERNAL_HARD
    PtclHard* global_cm; // c.m. of the hard cluster in the global frame
#endif

    ARPerturber(): NB(), soft_pert(NULL), soft_pert_min(Float(0.0)) 
#ifdef EXTERNAL_HARD
                 , global_cm(NULL)
#endif
    {}

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
#ifdef EXTERNAL_HARD
        ASSERT(global_cm!=NULL);
#endif
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
        PS::F64 r_min2=PS::LARGE_FLOAT;
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

    //! calculate soft_pert_min for slowdown pert_out
    /*! \Delta F = G m_cm m_p (apo) / rp^3
        Pert_out = \Delta F /(G apo)
        @param[in] _bin: binary parameters
        @param[in] _G: gravitatioal constant
     */
    template <class Tptcl>
    void calcSoftPertMin(const AR::BinaryTree<Tptcl>& _bin, const Float _G) {
        soft_pert_min = 0.0;
#ifdef SOFT_PERT
        if(soft_pert!=NULL&&_bin.semi>0.0) {
            ParticleBase p[2];
            _bin.calcParticlesEcca(p[0], p[1], COMM::PI, _G);
            Float dacc_soft = 0.0;
            Float acc_p1[3] = {0.0, 0.0, 0.0};
            Float acc_p2[3] = {0.0, 0.0, 0.0};
            soft_pert->eval(acc_p1, p[0].pos);
            soft_pert->eval(acc_p2, p[1].pos);
            Float dacc[3] = {acc_p1[0]-acc_p2[0], 
                             acc_p1[1]-acc_p2[1],
                             acc_p1[2]-acc_p2[2]};
            dacc_soft = std::sqrt(dacc[0]*dacc[0] + dacc[1]*dacc[1] + dacc[2]*dacc[2]);
            Float apo = _bin.semi*(1.0+_bin.ecc);
            soft_pert_min = _bin.mass*dacc_soft/(_G*apo);
        }
#endif
    }
};
