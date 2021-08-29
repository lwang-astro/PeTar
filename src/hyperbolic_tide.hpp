#pragma once

#include "Common/Float.h"
#include "Common/binary_tree.h"

//! HyperbolicTide 
/*! calculate orbital change due to dynamical tide or gravitational wave radiation for hyperbolic encounters
 */
class HyperbolicTide{
public:
    Float gravitational_constant;  
    Float speed_of_light;  // should have the same unit set as that of G

    HyperbolicTide(): gravitational_constant(-1), speed_of_light(-1) {}

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
      \return Etid: energy loss
     */
    template <class TBinary>
    void evolveOrbitGW(TBinary& _bin, Float& _Etid, Float& _Ltid) {

        // calculate energy loss (Hansen 1972, PRD, 5, 1021; correction from Turner 1977, ApJ, 216, 610)

        Float e = _bin.ecc;
        Float a = _bin.semi;
        Float m1 = _bin.m1;
        Float m2 = _bin.m2;

        ASSERT(e>1.0);
        ASSERT(m1>0);
        ASSERT(m2>0);
        ASSERT(a<0);

        Float c = speed_of_light;
        Float c2 = c*c;
        Float c5 = c2*c2*c;

        Float e2 = e*e;
        Float theta0 = acos(1/e);
        Float ge = (COMM::PI - theta0) * (96.0 + 292*e2 + 37*e2*e2) + sqrt(e2 - 1)/3.0 * (602.0 + 673*e2);

        Float mtot = m1 + m2;
        Float m12  = m1 * m2;
        Float m12sq = m12*m12;
        Float p = a * (1 - e2);
        Float G_over_r = gravitational_constant/p;
        Float G_over_r3 = G_over_r * G_over_r * G_over_r;
        Float G_over_r7 = G_over_r3 * G_over_r3 * G_over_r;
        Float Etid = 2.0 * sqrt(G_over_r7 * mtot) * m12sq * ge / (15.0 * c5);

        // calculate the angular momentum loss (Hansen 1972)
        ge = (COMM::PI - theta0) * (8.0 + 7*e2) + sqrt(e2 - 1) * (13 + 2*e2);
        Float Ltid = 8.0 * G_over_r3 * p * m12sq * ge / (15.0 * c5);


        // update binary orbit
        // update semi
        Float GM12 = gravitational_constant * m12;
        Float Ebin = - GM12 / (2.0 * a);
        Float Ebin_new = Ebin - Etid;
        _bin.semi = - GM12 / (2.0 * Ebin_new);


        // update ecc
        Float mfac = GM12 * m12 / mtot;
        Float Lbin = sqrt(mfac * p);
        Float Lbin_new = Lbin - Ltid;
        ASSERT(Lbin_new>=0);
        _bin.ecc = sqrt(1.0 - Lbin_new*Lbin_new / (mfac*_bin.semi));
        ASSERT(_bin.ecc>=0.0);
        
        _Etid = Etid;
        _Ltid = Ltid;
    }

    //! evolve hyperbolic orbit based on dynamical tide implementation from Alessandro Alberto Trani
    /*!
      @param[in,out] _bin: binary data, semi and ecc are updated
      @param[in] rad1: stellar radius of p1 (should be in the same unit of semi)
      @param[in] rad2: stellar radius of p2
      @param[in] poly_type: polynomial type

      \return Etid: energy loss
     */
    template <class TBinary>
    Float evolveOrbitPoly(TBinary& _bin, const Float& rad1, const Float& rad2, const Float& poly_type) {
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

        return Etid;
    }

    //! Tidal energy based on Polynomical fits from Portegies Zwart et al. 1993
    /*!
      @param[in] _mpert: perturber mass
      @param[in] _rad_over_peri: stellar radius over peri-center distance
      @param[in] _peri: peri-center distance
      @param[in] _eta: parameter (Press & Teutoslky 1977)
      @param[in] _poly_type: Polynomical type (1.5 or 3.0)
     */
    Float calcEtidPolynomicalFit(const Float& _mpert, const Float& _rad_over_peri, const Float& _peri, const Float& _eta, const Float& _poly_type) {
        // Tidal fits from Portegies Zwart et al. 1993
        Float fA=0, fB=0, fC=0, fD=0, fE=0, fF=0;
        if (_poly_type == 1.5) {
            fA = -0.397;
            fB = 1.678;
            fC = 1.277;
            fD = -12.42;
            fE = 9.446;
            fF = -5.550;
        }
        else if (_poly_type == 3.0) {
            fA = -1.124;
            fB = 0.877;
            fC = -13.37;
            fD = 21.55;
            fE = -16.8;
            fF = 4.124;
        } 
        else {
            std::cerr<<"Error, polynomical types should be 1.5, 3.0, given "<<_poly_type<<std::endl;
            abort();
        }
	
        // Tidal energy
        Float r_over_p2 =  _rad_over_peri * _rad_over_peri;
        Float r_over_p5 = r_over_p2 * r_over_p2 * _rad_over_peri;
        Float etid = gravitational_constant * r_over_p5 * _mpert * _mpert / _peri;

        // Now computing T(eta)
        if (_eta < 1) { // most likely a disruptive encounter
            etid *= pow(10, fA);
        } else if (_eta < 10) {
            Float let = log10(_eta);
            etid *= pow(10, fA + let * (fB + let * (fC + let * (fD + let * (fE + let * fF)))));
        }

        return etid;
    }

};
