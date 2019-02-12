#pragma once

#include "AR/Float.h"
#include "AR/list.h"
#include "ptclh4.hpp"

#ifndef TIDAL_TENSOR
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#endif

#ifdef TIDAL_TENSOR
//! Tidal tensor perterbation for AR
class TidalTensor{
private:
    PS::F64 T2[6];  // 1st (6)
    PS::F64 T3[10]; // 2nd Tensor (10)
    bool use_flag;
public:

    // only initialize flag
    TidalTensor(): use_flag(false) {}

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

    void reset(){
        use_flag=false;
        for(PS::S32 i=0; i<6; i++) T2[i] = 0;
        for(PS::S32 i=0; i<10; i++) T3[i] = 0;
    }

    //! tidal tensor fitting function,
    /* @param[in] _ptcl_tt: tidal tensor measure particles
       @param[in] _bin: binary information, get the scaling factor of distance
       @param[in] _n_split: artifical particle splitting number
     */
    template<class Tptcl>
    void fit(Tptcl* _ptcl_tt, const Binary& _bin, const PS::F64 _r_bin, const PS::S32 _n_split) {
        use_flag=true;
        PS::F64vec fi[8];
#ifdef HARD_DEBUG
        assert(_ptcl_tt[12].mass_bk==0);
        assert(_n_split>4);
#endif
        // get acceleration
        for (PS::S32 i=0; i<8; i++) fi[i] = _ptcl_tt[i].acc;
        // get cofficients
        for (PS::S32 i=0; i<8; i++) {
            T2[0] =  0.250000000000000*fi[0][0] + -0.250000000000000*fi[2][0] +  0.250000000000000*fi[4][0] + -0.250000000000000*fi[6][0];
            T2[1] =  0.125000000000000*fi[0][1] +  0.125000000000000*fi[1][0] + -0.125000000000000*fi[2][1] + -0.125000000000000*fi[3][0] 
                +    0.125000000000000*fi[4][1] +  0.125000000000000*fi[5][0] + -0.125000000000000*fi[6][1] + -0.125000000000000*fi[7][0];
            T2[2] = -0.083333333333333*fi[0][0] +  0.083333333333333*fi[0][2] + -0.083333333333333*fi[1][0] + -0.083333333333333*fi[2][0]
                +   -0.083333333333333*fi[2][2] + -0.083333333333333*fi[3][0] +  0.083333333333333*fi[4][0] +  0.083333333333333*fi[4][2]
                +    0.083333333333333*fi[5][0] +  0.083333333333333*fi[6][0] + -0.083333333333333*fi[6][2] +  0.083333333333333*fi[7][0];
            T2[3] =  0.250000000000000*fi[1][1] + -0.250000000000000*fi[3][1] +  0.250000000000000*fi[5][1] + -0.250000000000000*fi[7][1];
            T2[4] = -0.083333333333333*fi[0][1] + -0.083333333333333*fi[1][1] +  0.083333333333334*fi[1][2] + -0.083333333333333*fi[2][1]
                +   -0.083333333333333*fi[3][1] + -0.083333333333333*fi[3][2] +  0.083333333333333*fi[4][1] +  0.083333333333333*fi[5][1]
                +    0.083333333333333*fi[5][2] +  0.083333333333333*fi[6][1] +  0.083333333333333*fi[7][1] + -0.083333333333334*fi[7][2];
            T2[5] = -0.124999999999999*fi[0][2] + -0.125000000000000*fi[1][2] + -0.125000000000001*fi[2][2] + -0.125000000000000*fi[3][2]
                +    0.125000000000001*fi[4][2] +  0.125000000000000*fi[5][2] +  0.124999999999999*fi[6][2] +  0.125000000000000*fi[7][2];

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
        }
        // Rescale
        //PS::F64 T2S = 1.0/(_bin.semi*(1+_bin.ecc)*0.35);
        PS::F64 T2S = 1.0/(_r_bin*0.16);
        PS::F64 T3S = T2S*T2S;
        for (PS::S32 i=0; i<6;  i++) T2[i] *= T2S;
        for (PS::S32 i=0; i<10; i++) T3[i] *= T3S;
    }

    void eval(PS::F64* acc, const PS::F64vec &pos) const {
        if(!use_flag) return;
        /*
          T2:
          [[0 1 2]
          [1 3 4]
          [2 4 5]]

          T3:
          [[[6 7 8]
          [7 9 10]
          [8 10 11]]

          [[7 9 10]
          [9 12 13]
          [10 13 14]]

          [[8 10 11]
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

        acc0 +=  T2[0]*x + T2[1]*y + T2[2]*z 
            +      T3[0]*x2 + 2*T3[1]*xy + 2*T3[2]*xz + T3[3]*y2 + 2*T3[4]*yz + T3[5]*z2;
        acc1 +=  T2[1]*x + T2[3]*y + T2[4]*z
            +      T3[1]*x2 + 2*T3[3]*xy + 2*T3[4]*xz + T3[6]*y2 + 2*T3[7]*yz + T3[8]*z2;
        acc2 +=  T2[2]*x + T2[4]*y + T2[5]*z
            +      T3[2]*x2 + 2*T3[4]*xy + 2*T3[5]*xz + T3[7]*y2 + 2*T3[8]*yz + T3[9]*z2;

        acc[0] = acc0;
        acc[1] = acc1;
        acc[2] = acc2;
    }
};

////! substract c.m. force from component used for tidal tensor. 
///* @param[in]     _ptcl_cm: c.m. particle
//   @param[in,out] _ptcl_tt: tidal tensor artifical particle (8)
//*/
//template<class Tptcl>
//void subtractFcm(Tptcl& _ptcl_cm, Tptcl* _ptcl_tt) {
//    PS::F64vec& acc_cm = _ptcl_cm.acc;
// 
//    // remove center force
//    for (PS::S32 i=0; i<8; i++)  _ptcl_tt[i].acc -= acc_cm;
//}

#else
class keplerSplineFit{
private:
    const gsl_interp_type *t_;
    PS::F64 peri_;
    PS::F64 tperi_;
    gsl_interp_accel *acc_[3];
    gsl_spline *spline_[3];
public:

    keplerSplineFit(): t_(gsl_interp_cspline_periodic), peri_(0.0) {
        acc_[0]=acc_[1]=acc_[2]=NULL;
        spline_[0]=spline_[1]=spline_[2]=NULL;
    }
    
    //! orbital force fitting function,
    /* @param[in] _ptcl_orb: orbital particles
       @param[in] _bin: binary information, get the scaling factor of distance
       @param[in] _n_split: artifical particle splitting number
     */
    template<class Tptcl> 
    void fit(Tptcl* _ptcl_orb, const Binary& _bin, const PS::S32 _n_split) {
#ifdef HARD_DEBUG
        assert(_n_split>=4);
#endif
        const PS::F64 twopi = PI*2.0;
        peri_  = _bin.peri;
        tperi_ = _bin.tperi - twopi;  // make sure tperi_ is negative
        const PS::U32 np = _n_split+1;
        PS::F64 x[np],y[3][np];
        for(PS::U32 i=0; i<np; i++) {
            x[i] = PS::F64(i)/_n_split*twopi;
            x[i] = _bin.peri/twopi * (x[i] - _bin.ecc*std::sin(x[i]));
        }
        for(PS::S32 i=0; i<_n_split; i++) {            
            for(PS::S32 j=0; j<3; j++)
                y[j][i] = _ptcl_orb[2*i+1].acc[j] - _ptcl_orb[2*i].acc[j];  // acc1 - acc0
        }
        for(PS::S32 i=0; i<3; i++) {
            y[i][_n_split] = y[i][0];

            if(acc_[i]!=NULL) gsl_interp_accel_free(acc_[i]);
            acc_[i] = gsl_interp_accel_alloc();
            if(spline_[i]!=NULL) gsl_spline_free(spline_[i]);
            spline_[i] =  gsl_spline_alloc(t_, np);
            
            gsl_spline_init(spline_[i], x, y[i], np);
        }
    }

    void eval(PS::F64* acc, const PS::F64 time) const {
        PS::F64 dt = time - tperi_;
        dt = dt - (PS::S32)(dt/peri_)*peri_;
        if(spline_[0]!=NULL) {
            for(PS::S32 i=0; i<3; i++)
                acc[i] -= gsl_spline_eval(spline_[i],dt,acc_[i]);
        }
    }

    void reset(){
        for(PS::S32 i=0;i<3;i++) {
            if(acc_[i]!=NULL) gsl_interp_accel_free(acc_[i]);
            if(spline_[i]!=NULL) gsl_spline_free(spline_[i]);
            acc_[i]=NULL;
            spline_[i]=NULL;
        }
    }

    ~keplerSplineFit() {
        reset();
    }
};
#endif

//! Perturber class for AR integration
class ARPerturber{
public:
    List<PtclH4*> neighbor; //> neighbor perturbers
    
    
    ARPerturber(): neighbor() {}

    void clear() {
        pert_single = NULL;
        pert_group  = NULL;
        nb_single   = NULL;
        nb_group    = NULL;
        pert_soft   = NULL;
        cm_self     = NULL;
    }

};
