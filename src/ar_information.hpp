#pragma once

#include "AR/list.h"
#include "AR/Float.h"
#include "ptcl.hpp"
#include "kepler.hpp"

//! contain group information
class ARInformation{
public:
    AR::List<int> particle_index;
    AR::List<PtclTree<Ptcl>> binarytree;
    
    //! reserve memory
    void reserveMem(const int _nmax) {
        particle_index.setMode(ListMode::local);
        binarytree.setMode(ListMode::local);
        binarytree.reserveMem(_nmax);
        particle_index.reserveMem(_nmax);
    }
    
    //! generate binary tree 
    void generateBinaryTree(Ptcl* _particles, const int _n_particle, const Float _dt_tree, const Float _v_max) {
        binarytree.resizeNoInitialize(_n_particle-1);
        particle_index.resizeNoInitialize(_n_particle);
        int* index = particle_index.getDataAddress();
        for (int i=0; i<_n_particle; i++) index[i]=i;
        keplerTreeGenerator(binarytree.getDataAddress(), index, _n_particle, _particles, _dt_tree, _v_max);
    }

    PtclTree<Ptcl>& getBinaryTreeRoot() const {
        int n = binarytree.getSize();
        assert(n>0);
        return binarytree[n-1];
    }

    //! get dr * dv for two particles
    /*!
      @param[out] _dr2: dr*dr
      @param[out] _drdr: dr*dv
      @param[in] _p1: particle 1
      @param[in] _p2: particle 2
     */
    void getDrDv(Float& _dr2, Float& _drdv, const Ptcl& _p1, const Ptcl& _p2) {
        Float dx[3],dv[3];
        Float* pos1 = _p1.pos;
        Float* pos2 = _p2.pos;
        Float* vel1 = _p1.vel;
        Float* vel2 = _p2.vel;
        dx[0] = pos1[0] - pos2[0];
        dx[1] = pos1[1] - pos2[1];
        dx[2] = pos1[2] - pos2[2];

        dv[0] = vel1[0] - vel2[0];
        dv[1] = vel1[1] - vel2[1];
        dv[2] = vel1[2] - vel2[2];
        
        _dr2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
        _drdv= dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2];
    }
};
