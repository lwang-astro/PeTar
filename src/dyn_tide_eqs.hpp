#pragma once

#include "Common/Float.h"
#include "Common/binary_tree.h"

//! Dynamic Tide 
class DynamicTide{
public:
    Float gravitational_constant;
    Float speed_of_light;

    DynamicTide(): gravitational_constant(-1), speed_of_light(-1) {}

    //! check whether parameters values are correct
    /*! \return true: all correct
     */
    bool checkParams() {
        ASSERT(gravitational_constant>0.0);
        ASSERT(speed_of_light>0.0);
        return true;
    }        

    //! evolve hyperbolic orbit based on GW effect
    /*!
      @param[in] _bin: binary data, semi and ecc are updated
     */
    template <class TBinary>
    void evolveOrbitGW(TBinary& _bin) {
        Float e = _bin.ecc;
        Float mtot = _bin.m1+_bin.m2;
        Float m12 = _bin.m1*_bin.m2;
        Float peri = _bin.semi * (1 - _bin.ecc);
        Float e2 = e*e;
        Float c = speed_of_light;
        Float c2 = c*c;
        Float c5 = c2*c2*c;
        Float ge = (24.0*acos(-1/e) * (1+73.0/24.0*e2 + 37.0/96.0*e2*e2) 
                    + sqrt(e2-1)*(301.0/6.0+673.0/12.0*e2)) / pow(1+e, 3.5);
        Float G_over_r = gravitational_constant/peri;
        Float G_over_r3 = G_over_r * G_over_r * G_over_r;
        Float G_over_r7 = G_over_r3 * G_over_r3 * G_over_r;
        Float Etid = 8.0 * sqrt(G_over_r7 * mtot) * m12*m12 * ge /(15.0*c5);

        // assuming angular momentum conserved, the semi-latus rectum is also conserved
        Float pold = _bin.semi *(1.0 - _bin.ecc*_bin.ecc);

        // update binary orbit
        // update semi
        Float GM12 = gravitational_constant * m12;
        Float Ebin = - GM12 / (2.0 * _bin.semi);
        Float Ebin_new = Ebin - Etid;
        _bin.semi = - GM12 / (2.0 * Ebin_new);

        // use semi-latus rectum  to update ecc
        _bin.ecc = sqrt(1.0 - pold/_bin.semi);
        ASSERT(_bin.ecc>=0.0);
        
    }

    //! evolve hyperbolic orbit based on dynamical tide implementation from Alessandro Alberto Trani
    /*!
      @param[in,out] _bin: binary data, semi and ecc are updated
      @param[in] rad1: stellar radius of p1 (should be in the same unit of semi)
      @param[in] rad2: stellar radius of p2
      @param[in] poly_type: polynomial type
     */
    template <class TBinary>
    void evolveOrbitPoly(TBinary& _bin, const Float& rad1, const Float& rad2, const Float& poly_type) {
        ASSERT(_bin.getMemberN()==2);
        ASSERT(_bin.semi<0);
        ASSERT(_bin.ecc>1);

        ASSERT(_bin.m1>0);
        ASSERT(_bin.m2>0);
        ASSERT(poly_type==1.5||poly_type==3.0);

		Float peri = _bin.semi * (1.0 - _bin.ecc);
        ASSERT(peri>rad1+rad2);
        
        //Float rad1 = _bse_manager.getStellarRadius(p1->star);
        //Float rad2 = _bse_manager.getStellarRadius(p2->star);
        Float mtot = _bin.m1 + _bin.m2;

		// Press & Taukolsky 1977
		Float p_over_r1 = peri / rad1;
		Float p_over_r2 = peri / rad2;
		Float r1_over_peri = 1.0 / p_over_r1;
		Float r2_over_peri = 1.0 / p_over_r2;
		Float eta1 = sqrt(_bin.m1 / mtot * p_over_r1 * p_over_r1 * p_over_r1);
		Float eta2 = sqrt(_bin.m2 / mtot * p_over_r2 * p_over_r2 * p_over_r2);

		// Mardling & Aarseth 2001
		Float mard = 0.5 * fabs(eta1 - 2.);
		mard = (0.5 + 0.25 * sqrt(mard * mard * mard));
		Float neweta1 = eta1 * pow(2.0 / (1.0 + _bin.ecc), mard);

		mard = 0.5 * fabs(eta2 - 2.);
		mard = (0.5 + 0.25 * sqrt(mard * mard * mard));
		Float neweta2 = eta2 * pow(2. / (1.0 + _bin.ecc), mard);

		Float Etid = 0;

		// TIDE ON 1
		if ((neweta1 > 0) & (neweta1 < 10)) {
			Etid += calcEtidPolynomicalFit(_bin.m2, r1_over_peri, peri, neweta1, poly_type);
		}

		// TIDE ON 2
		if ((neweta2 > 0) & (neweta2 < 10)) {
			Etid += calcEtidPolynomicalFit(_bin.m1, r2_over_peri, peri, neweta2, poly_type);
		}

        // assuming angular momentum conserved, the semi-latus rectum is also conserved
        Float pold = _bin.semi *(1.0 - _bin.ecc*_bin.ecc);

        // update binary orbit
        // update semi
        Float GM12 = gravitational_constant * _bin.m1 * _bin.m2;
        Float Ebin = - GM12 / (2.0 * _bin.semi);
        Float Ebin_new = Ebin - Etid;
        _bin.semi = - GM12 / (2.0 * Ebin_new);

        // use semi-latus rectum  to update ecc
        _bin.ecc = sqrt(1.0 - pold/_bin.semi);
        ASSERT(_bin.ecc>=0.0);
        //_bin.calcParticles(gravitational_constant);
        //p1->pos += _bin.pos;
        //p2->pos += _bin.pos;
        //p1->vel += _bin.vel;
        //p2->vel += _bin.vel;
    }

    //! Tidal energy based on Polynomical fits from Portegies Zwart et al. 1993
    /*!
      @param[in] mpert: perturber mass
      @param[in] r_peri: stellar radius over peri-center distance
      @param[in] peri: peri-center distance
      @param[in] eta: parameter (Press & Teutoslky 1977)
      @param[in] poly_type: Polynomical type (1.5 or 3.0)
     */
    Float calcEtidPolynomicalFit(const Float& mpert, const Float& r_peri, const Float& peri, const Float& eta, const Float& poly_type) {

        // Tidal fits from Portegies Zwart et al. 1993
        const Float fA_n15 = -0.397;
        const Float fB_n15 = 1.678;
        const Float fC_n15 = 1.277;
        const Float fD_n15 = -12.42;
        const Float fE_n15 = 9.446;
        const Float fF_n15 = -5.550;

        const Float fA_n3 = -1.124;
        const Float fB_n3 = 0.877;
        const Float fC_n3 = -13.37;
        const Float fD_n3 = 21.55;
        const Float fE_n3 = -16.8;
        const Float fF_n3 = 4.124;

        Float fA, fB, fC, fD, fE, fF;
        if (poly_type == 3) {
            fA = fA_n3, fB = fB_n3, fC = fC_n3, fD = fD_n3, fE = fE_n3, fF = fF_n3;
        } else if (poly_type == 1.5) {
            fA = fA_n15, fB = fB_n15, fC = fC_n15, fD = fD_n15, fE = fE_n15, fF = fF_n15;
        } else {
            return 0.0;
        }
	
        // Tidal energy
        Float etid = gravitational_constant * r_peri * r_peri * r_peri * r_peri * r_peri * mpert * mpert / peri;
        // Now computing Teta
        if (eta < 1) { // most likely a disruptive encounter
            etid *= pow(fA, 10);
        } else if (eta < 10) {
            Float let = log10(eta);
            etid *= pow(fA + let * (fB + let * (fC + let * (fD + let * (fE + let * fF)))), 10);
        }

        return etid;
    }

};
