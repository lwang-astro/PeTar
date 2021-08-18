#pragma once

#include "Common/Float.h"
#include "Common/binary_tree.h"

//! Dynamic Tide based on implementation from Alessandro Alberto Trani
class DynamicTide{
public:
    Float gravitational_constant;
    Float poly_type;

    DynamicTide(): gravitational_constant(-1) {}

    //! check whether parameters values are correct
    /*! \return true: all correct
     */
    bool checkParams() {
        ASSERT(gravitational_constant>0.0);
        ASSERT(poly_type==1.5 || poly_type==3.0);
        return true;
    }        

    //! evolve hyperbolic orbit based on dynamical tide
    /*!
      @param[in,out] _bin: binary data
      @param[in] rad1: stellar radius of p1 (should be in the same unit of semi)
      @param[in] rad2: stellar radius of p2
     */
    template <class TBinary>
    void evolveOrbit(TBinary& _bin, const Float& rad1, const Float& rad2) {
        ASSERT(_bin.getMemberN()==2);
        ASSERT(_bin.semi<0);
        ASSERT(_bin.ecc>1);

        auto* p1 = _bin.getLeftMember();
        auto* p2 = _bin.getRightMember();
        ASSERT(p1->mass>0);
        ASSERT(p2->mass>0);

        COMM::Vector3<Float> dr(p2->pos[0] - p1->pos[0], p2->pos[1] - p1->pos[1], p2->pos[2] - p1->pos[2]);
        COMM::Vector3<Float> dv(p2->vel[0] - p1->vel[0], p2->vel[1] - p1->vel[1], p2->vel[2] - p1->vel[2]);
        
		Float peri = _bin.semi * (1.0 - _bin.ecc);
        
        //Float rad1 = _bse_manager.getStellarRadius(p1->star);
        //Float rad2 = _bse_manager.getStellarRadius(p2->star);
        Float mtot = p1->mass + p2->mass;

		// Press & Taukolsky 1977
		Float p_over_r1 = peri / rad1;
		Float p_over_r2 = peri / rad2;
		Float r1_over_peri = 1.0 / p_over_r1;
		Float r2_over_peri = 1.0 / p_over_r2;
		Float eta1 = sqrt(p1->mass / mtot * p_over_r1 * p_over_r1 * p_over_r1);
		Float eta2 = sqrt(p2->mass / mtot * p_over_r2 * p_over_r2 * p_over_r2);

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
			Etid += calcEtidPolynomicalFit(p2->mass, r1_over_peri, peri, neweta1, poly_type);
		}

		// TIDE ON 2
		if ((neweta2 > 0) & (neweta2 < 10)) {
			Etid += calcEtidPolynomicalFit(p1->mass, r2_over_peri, peri, neweta2, poly_type);
		}

        // update binary orbit
        // Keep peri-center distance the same, update semi and then ecc
        Float GM12 = gravitational_constant * p1->mass * p2->mass;
        Float Ebin = - GM12 / (2.0 * _bin.semi);
        Float Ebin_new = Ebin - Etid;
        _bin.semi = - GM12 / (2.0 * Ebin_new);
        _bin.ecc  = 1.0 - peri/_bin.semi;
        ASSERT(_bin.ecc>=0.0);
        _bin.calcParticles(gravitational_constant);
        p1->pos += _bin.pos;
        p2->pos += _bin.pos;
        p1->vel += _bin.vel;
        p2->vel += _bin.vel;
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
        Float etid = r_peri * r_peri * r_peri * r_peri * r_peri * mpert * mpert / peri;
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
