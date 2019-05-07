#pragma once

#include "Common/list.h"
#include "Hermite/hermite_particle.h"
#include "hard_ptcl.hpp"
#include "Common/binary_tree.h"

//! Tidal tensor perterbation for AR
class TidalTensor{
private:
    PS::F64 T1[3];     // 0 constant 
    PS::F64 T2[9];  // 1st (9)    general tensor
    //PS::F64 T2[6];  // 1st (6)  symmetry tensor
    PS::F64 T3[10]; // 2nd Tensor (10)
public:
    PS::F64vec pos;  // position of c.m.
    PS::S32 group_id; // indicate which group use the tensor

    TidalTensor(): T1{0.0,0.0,0.0}, T2{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, T3{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, pos(0.0), group_id(-1) {}

    void dump(FILE *fp){
        fwrite(this, sizeof(TidalTensor),1,fp);
    }

    void read(FILE *fp){
        size_t rcount = fread(this, sizeof(TidalTensor),1,fp);
        if (rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            abort();
        }
    }

    void clear(){
        T1[0] = T1[1] = T1[2] = 0;
        for(PS::S32 i=0; i<9; i++) T2[i] = 0;
        for(PS::S32 i=0; i<10; i++) T3[i] = 0;
        pos = 0.0;
        group_id = -1;
    }

    //! tidal tensor fitting function,
    /*! 
       Symmetry T2:
       xx xy xz 0  1  2
       yx yy yz 1  3  4
       zx zy zz 2  4  5

       General T2:
       xx xy xz 0  1  2
       yx yy yz 3  4  5
       zx zy zz 6  7  8

       Symmetry T3:
       xxx xxy xxz  0  1  2
       xyx xyy xyz  1  3  4
       xzx xzy xzz  2  4  5

       yxx yxy yxz  1  3  4
       yyx yyy yyz  3  6  7  
       yzx yzy yzz  4  7  8
      
       zxx zxy zxz  2  4  5
       zyx zyy zyz  4  7  8
       zzx zzy zzz  5  8  9
       
       @param[in] _ptcl_tt: tidal tensor measure particles
       @param[in] _ptcl_cm: tidal tensor measure particle c.m.
       @param[in] _r_bin: particle box size
       @param[in] _n_split: artifical particle splitting number
    */
    template<class Tptcl>
    void fit(Tptcl* _ptcl_tt, Tptcl& _ptcl_cm,  const PS::F64 _r_bin, const PS::S32 _n_split) {
        // get c.m. position
        pos = _ptcl_cm.pos;

        PS::F64vec fi[8];

        // Marked the mass_bk == 0 in generate orbits for consistent check
        ASSERT(_ptcl_tt[12].mass_bk==0);
        ASSERT(_n_split>4);

        // get acceleration
        for (PS::S32 i=0; i<8; i++) fi[i] = _ptcl_tt[i].acc;

        // get cofficients
        // T1, assume input force already remove the c.m.
        T1[0] = T1[1] = T1[2] = 0.0;
        
        // T2, general form
        // 0 1 2
        T2[0] =  0.250000000000000*fi[0][0] + -0.250000000000000*fi[2][0] +  0.250000000000000*fi[4][0] + -0.250000000000000*fi[6][0];
        T2[1] =  0.125000000000000*fi[0][1] +  0.125000000000000*fi[1][0] + -0.125000000000000*fi[2][1] + -0.125000000000000*fi[3][0] 
            +    0.125000000000000*fi[4][1] +  0.125000000000000*fi[5][0] + -0.125000000000000*fi[6][1] + -0.125000000000000*fi[7][0];
        T2[2] = -0.083333333333333*fi[0][0] +  0.083333333333333*fi[0][2] + -0.083333333333333*fi[1][0] + -0.083333333333333*fi[2][0]
            +   -0.083333333333333*fi[2][2] + -0.083333333333333*fi[3][0] +  0.083333333333333*fi[4][0] +  0.083333333333333*fi[4][2]
            +    0.083333333333333*fi[5][0] +  0.083333333333333*fi[6][0] + -0.083333333333333*fi[6][2] +  0.083333333333333*fi[7][0];

        // 3 4 5
        T2[3] =  T2[1];
        T2[4] =  0.250000000000000*fi[1][1] + -0.250000000000000*fi[3][1] +  0.250000000000000*fi[5][1] + -0.250000000000000*fi[7][1];
        T2[5] = -0.083333333333333*fi[0][1] + -0.083333333333333*fi[1][1] +  0.083333333333334*fi[1][2] + -0.083333333333333*fi[2][1]
            +   -0.083333333333333*fi[3][1] + -0.083333333333333*fi[3][2] +  0.083333333333333*fi[4][1] +  0.083333333333333*fi[5][1]
            +    0.083333333333333*fi[5][2] +  0.083333333333333*fi[6][1] +  0.083333333333333*fi[7][1] + -0.083333333333334*fi[7][2];

        // 6 7 8
        T2[6] =  T2[2];
        T2[7] =  T2[4];
        T2[8] = -0.124999999999999*fi[0][2] + -0.125000000000000*fi[1][2] + -0.125000000000001*fi[2][2] + -0.125000000000000*fi[3][2]
            +    0.125000000000001*fi[4][2] +  0.125000000000000*fi[5][2] +  0.124999999999999*fi[6][2] +  0.125000000000000*fi[7][2];

        // T3, symmetry form
        T3[0] =  0.250000000000001*fi[0][0] +  0.124999999999999*fi[0][2] +  0.250000000000001*fi[2][0] + -0.124999999999999*fi[2][2]
            +    0.249999999999999*fi[4][0] + -0.125000000000001*fi[4][2] +  0.249999999999999*fi[6][0] +  0.125000000000001*fi[6][2];
        T3[1] =  0.250000000000000*fi[0][1] +  0.125000000000000*fi[1][2] +  0.250000000000000*fi[2][1] + -0.125000000000000*fi[3][2]
            +    0.250000000000000*fi[4][1] + -0.125000000000000*fi[5][2] +  0.250000000000000*fi[6][1] +  0.125000000000000*fi[7][2];
        T3[2] = -0.112500000000000*fi[0][0] +  0.025000000000000*fi[0][2] + -0.012500000000000*fi[1][1] + -0.025000000000000*fi[1][2]
            +    0.112500000000000*fi[2][0] +  0.025000000000000*fi[2][2] +  0.012500000000000*fi[3][1] + -0.025000000000000*fi[3][2]
            +    0.112500000000000*fi[4][0] +  0.025000000000000*fi[4][2] +  0.012500000000000*fi[5][1] + -0.025000000000000*fi[5][2]
            +   -0.112500000000000*fi[6][0] +  0.025000000000000*fi[6][2] + -0.012500000000000*fi[7][1] + -0.025000000000000*fi[7][2];
        T3[3] =  0.125000000000000*fi[0][2] +  0.250000000000000*fi[1][0] + -0.125000000000000*fi[2][2] +  0.250000000000000*fi[3][0]
            +   -0.124999999999999*fi[4][2] +  0.250000000000000*fi[5][0] +  0.125000000000000*fi[6][2] +  0.250000000000000*fi[7][0];
        T3[4] = -0.062500000000000*fi[0][1] + -0.062500000000000*fi[1][0] +  0.062500000000000*fi[2][1] +  0.062500000000000*fi[3][0]
            +    0.062500000000000*fi[4][1] +  0.062500000000000*fi[5][0] + -0.062500000000000*fi[6][1] + -0.062500000000000*fi[7][0];
        T3[5] = -0.125000000000000*fi[0][2] +  0.125000000000000*fi[2][2] +  0.125000000000000*fi[4][2] + -0.125000000000000*fi[6][2];
        T3[6] =  0.250000000000000*fi[1][1] +  0.125000000000000*fi[1][2] +  0.250000000000000*fi[3][1] + -0.125000000000000*fi[3][2]
            +    0.250000000000000*fi[5][1] + -0.125000000000000*fi[5][2] +  0.250000000000000*fi[7][1] +  0.125000000000000*fi[7][2];
        T3[7] = -0.012500000000000*fi[0][0] + -0.025000000000000*fi[0][2] + -0.112500000000000*fi[1][1] +  0.025000000000000*fi[1][2]
            +    0.012500000000000*fi[2][0] + -0.025000000000000*fi[2][2] +  0.112500000000000*fi[3][1] +  0.025000000000000*fi[3][2]
            +    0.012500000000000*fi[4][0] + -0.025000000000000*fi[4][2] +  0.112500000000000*fi[5][1] +  0.025000000000000*fi[5][2]
            +   -0.012500000000000*fi[6][0] + -0.025000000000000*fi[6][2] + -0.112500000000000*fi[7][1] +  0.025000000000000*fi[7][2];
        T3[8] = -0.125000000000000*fi[1][2] +  0.125000000000000*fi[3][2] +  0.125000000000000*fi[5][2] + -0.125000000000000*fi[7][2];
        T3[9] =  0.062500000000000*fi[0][0] +  0.125000000000000*fi[0][2] +  0.062500000000000*fi[1][1] +  0.125000000000000*fi[1][2]
            +   -0.062500000000000*fi[2][0] +  0.125000000000000*fi[2][2] + -0.062500000000000*fi[3][1] +  0.125000000000000*fi[3][2]
            +   -0.062500000000000*fi[4][0] +  0.125000000000000*fi[4][2] + -0.062500000000000*fi[5][1] +  0.125000000000000*fi[5][2]
            +    0.062500000000000*fi[6][0] +  0.125000000000000*fi[6][2] +  0.062500000000000*fi[7][1] +  0.125000000000000*fi[7][2];

        // Rescale
        //PS::F64 T2S = 1.0/(_bin.semi*(1+_bin.ecc)*0.35);
        PS::F64 T2S = 1.0/(_r_bin*0.16);
        PS::F64 T3S = T2S*T2S;
        for (PS::S32 i=0; i<6;  i++) T2[i] *= T2S;
        for (PS::S32 i=0; i<10; i++) T3[i] *= T3S;
    }

    //! Shift c.m. to new reference position
    /*! Only the 1st order tensor need a correction from 2nd order 
      
      ### 1st order:
      T2:
      xx xy xz 0  1  2
      yx yy yz 3  4  5
      zx zy zz 6  7  8

      ### 2nd order:
      xxx xxy xxz  0  1  2
      xyx xyy xyz  1  3  4
      xzx xzy xzz  2  4  5

      yxx yxy yxz  1  3  4
      yyx yyy yyz  3  6  7  
      yzx yzy yzz  4  7  8
      
      zxx zxy zxz  2  4  5
      zyx zyy zyz  4  7  8
      zzx zzy zzz  5  8  9
      
      @param[in] _pos: new c.m. position
     */
    void shiftCM(const PS::F64vec & _pos) {
        PS::F64vec dr = _pos-pos;

        PS::F64 x = dr.x;
        PS::F64 y = dr.y;
        PS::F64 z = dr.z;
        PS::F64 x2 = x*x;
        PS::F64 xy = x*y;
        PS::F64 xz = x*z;
        PS::F64 y2 = y*y;
        PS::F64 yz = y*z;
        PS::F64 z2 = z*z;

        // T1 += T2^dr + dr^T3^dr
        T1[0] += T2[0]*x + T2[1]*y + T2[2]*z 
            +    T3[0]*x2 + 2*T3[1]*xy + 2*T3[2]*xz + T3[3]*y2 + 2*T3[4]*yz + T3[5]*z2;
        T1[1] += T2[3]*x + T2[4]*y + T2[5]*z
            +    T3[1]*x2 + 2*T3[3]*xy + 2*T3[4]*xz + T3[6]*y2 + 2*T3[7]*yz + T3[8]*z2;
        T1[2] += T2[6]*x + T2[7]*y + T2[8]*z
            +    T3[2]*x2 + 2*T3[4]*xy + 2*T3[5]*xz + T3[7]*y2 + 2*T3[8]*yz + T3[9]*z2;
        
        // T2 += 2*dr^T3
        T2[0] += 2.0*(T3[0]*x + T3[1]*y + T3[2]*z); // xx: xxx*x + xyx*y + xzx*z
        T2[1] += 2.0*(T3[1]*x + T3[3]*y + T3[4]*z); // xy: xxy*x + xyy*y + xzy*z
        T2[2] += 2.0*(T3[2]*x + T3[4]*y + T3[5]*z); // xy: xxz*x + xyz*y + xzz*z

        T2[3] += 2.0*(T3[1]*x + T3[3]*y + T3[4]*z); // yx: yxx*x + yyx*y + yzx*z
        T2[4] += 2.0*(T3[3]*x + T3[6]*y + T3[7]*z); // yy: yxy*x + yyy*y + yzy*z
        T2[5] += 2.0*(T3[4]*x + T3[7]*y + T3[8]*z); // yy: yxz*x + yyz*y + yzz*z

        T2[6] += 2.0*(T3[2]*x + T3[4]*y + T3[5]*z); // zx: zxx*x + zyx*y + zzx*z
        T2[7] += 2.0*(T3[4]*x + T3[7]*y + T3[8]*z); // zy: zxy*x + zyy*y + zzy*z
        T2[8] += 2.0*(T3[5]*x + T3[8]*y + T3[9]*z); // zy: zxz*x + zyz*y + zzz*z

        // update c.m.
        pos = _pos;
    }

    void eval(PS::F64* acc, const PS::F64vec &pos) const {
        /*
          T2:
          [[0 1 2]
           [3 4 5]
           [6 7 8]]

          T3:
          [[[6  7  8 ]
            [7  9  10]
            [8  10 11]]

           [[7  9  10]
            [9  12 13]
            [10 13 14]]

           [[8  10 11]
            [10 13 14]
            [11 14 15]]]

        */
        PS::F64 x = pos.x;
        PS::F64 y = pos.y;
        PS::F64 z = pos.z;
        PS::F64 x2 = x*x;
        PS::F64 xy = x*y;
        PS::F64 xz = x*z;
        PS::F64 y2 = y*y;
        PS::F64 yz = y*z;
        PS::F64 z2 = z*z;

        PS::F64 acc0=acc[0];
        PS::F64 acc1=acc[1];
        PS::F64 acc2=acc[2];

        acc0 +=  T1[0] + T2[0]*x + T2[1]*y + T2[2]*z 
            +      T3[0]*x2 + 2*T3[1]*xy + 2*T3[2]*xz + T3[3]*y2 + 2*T3[4]*yz + T3[5]*z2;
        acc1 +=  T1[1] + T2[3]*x + T2[4]*y + T2[5]*z
            +      T3[1]*x2 + 2*T3[3]*xy + 2*T3[4]*xz + T3[6]*y2 + 2*T3[7]*yz + T3[8]*z2;
        acc2 +=  T1[2] + T2[6]*x + T2[7]*y + T2[8]*z
            +      T3[2]*x2 + 2*T3[4]*xy + 2*T3[5]*xz + T3[7]*y2 + 2*T3[8]*yz + T3[9]*z2;

        acc[0] = acc0;
        acc[1] = acc1;
        acc[2] = acc2;
    }

    PS::F64 evalPot(const PS::F64vec &pos) const {
        PS::F64 x = pos.x;
        PS::F64 y = pos.y;
        PS::F64 z = pos.z;
        PS::F64 x2 = x*x;
        PS::F64 xy = x*y;
        PS::F64 xz = x*z;
        PS::F64 y2 = y*y;
        PS::F64 yz = y*z;
        PS::F64 z2 = z*z;

        PS::F64 acc0 =  T1[0] + T2[0]*x + T2[1]*y + T2[2]*z 
            +      T3[0]*x2 + 2*T3[1]*xy + 2*T3[2]*xz + T3[3]*y2 + 2*T3[4]*yz + T3[5]*z2;
        PS::F64 acc1 =  T1[1] + T2[1]*x + T2[3]*y + T2[4]*z
            +      T3[1]*x2 + 2*T3[3]*xy + 2*T3[4]*xz + T3[6]*y2 + 2*T3[7]*yz + T3[8]*z2;
        PS::F64 acc2 =  T1[2] + T2[2]*x + T2[4]*y + T2[5]*z
            +      T3[2]*x2 + 2*T3[4]*xy + 2*T3[5]*xz + T3[7]*y2 + 2*T3[8]*yz + T3[9]*z2;
        
        return - x*acc0 - y*acc1 - z*acc2;
    }
};


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
            soft_pert->group_id = -1;
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

    //! find close tidal tensor and initial tidal tensor c.m.
    /*! if the tidal tensor is already in used (group_id>=0), copy a new one after _n_tt
      @param[in,out] _tt: tensor array
      @param[in,out] _n_tt: number of current tensor
      @param[in] _n_max: maximum size of tensor array
      @param[in] _cm: c.m. particle
      @param[in] _gid: group id
     */
    void findCloseSoftPert(TidalTensor* _tt, int& _n_tt, const int _n_max, const H4::ParticleH4<PtclHard>& _cm, const PS::S32 _gid) {
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
        if (_tt[r_min_index].group_id>=0) {
            ASSERT(_n_tt<_n_max);
            _tt[_n_tt] = _tt[r_min_index];
            soft_pert = &_tt[_n_tt];
            _n_tt++;
        }
        else soft_pert = &_tt[r_min_index];
        soft_pert->group_id = _gid;
        // update c.m.
        soft_pert->shiftCM(pos);
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
            pert = std::sqrt(dacc[0]*dacc[0] + dacc[1]*dacc[1] + dacc[2]*dacc[2]);
        }
#endif
        return pert;
    }

    //! calculate soft_pert_min
    template <class Tptcl>
    void calcSoftPertMin(const COMM::BinaryTree<Tptcl>& _bin) {
        // hyperbolic case
        if(_bin.semi<0.0) soft_pert_min = 0.0;
        else { // close orbit
            ParticleBase p[2];
            _bin.calcParticlesEcca(p[0], p[1], COMM::PI);
            Float dacc_soft = calcSoftPertSlowDownBinary(p[0], p[1]);
            soft_pert_min = _bin.mass*dacc_soft/(2.0*abs(_bin.semi));
        }
    }
};
