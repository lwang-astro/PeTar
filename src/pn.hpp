#pragma once

#include <cmath>
#include "Common/Float.h"
/***************************************************************************/
/*
  Coded by       : Peter Berczik (on the base of Gabor Kupi original PN code)
  Version number : 1.0
  Last redaction : 2011.V.6. 11:00

  Modified by Long Wang (2024), fixed several errors of formula 
  1. All PN jerk formula use A[K] which exclude Newton acceleration. This is not consistent, use AT[k] instead.
  2. PN2 jerk calculation A2D 3rd term is not correct, V2R -> V2R*rdot*rdot*rdot
*/

class PostNewtonian{
public:
    Float speed_of_light;
    Float gravitational_constant;
    Float precession_criterion; // precession angle to determine which order of PN should be on
    
    PostNewtonian(): speed_of_light(1.0), gravitational_constant(1.0), precession_criterion(1.0e-8) {}
    
    //! check whether parameters values are correct
    /*! \return true: all correct
     */
    bool checkParams() {
        ASSERT(speed_of_light>0.0);
        ASSERT(gravitational_constant>0.0);
        return true;
    }

    //! print parameters
    void print(std::ostream & _fout) const{
        _fout<<"speed_of_light: "<<speed_of_light<<std::endl
             <<"G: "<<gravitational_constant<<std::endl
             <<"precession_criterion: "<<precession_criterion<<std::endl;
    }

    //! determine which PN terms should be switched on based on an estimator due to Einstein drift
    /*!
      If satisfy, switch on PN1, 2, 2.5, others are not yet tested.
      @param[out] used_pn_orders: determine which PN orders are calculated: set in a bool array: [PN1, PN2, PN2.5, PN3, PN3.5, SPIN]
      @param[in] r: distance between two particles
      @param[in] m: total mass of two particles
      
      return: if use PN: true
    */
    bool setUsedPNOrders(bool used_pn_orders[6], Float r, Float m) {
        const Float c = speed_of_light;
        Float c2 = c*c;
        Float theta = 4*gravitational_constant*m/(r*c2); // precession angle 
        bool used_pn = false;
        // PN1
        if (theta > precession_criterion) {
            used_pn_orders[0] = true;
            used_pn_orders[1] = true;
            used_pn_orders[2] = true;
            used_pn_orders[3] = false;
            used_pn_orders[4] = false;
            used_pn_orders[5] = false;
            used_pn = true;
        }
        else {
            for (int i=0; i<6; i++) used_pn_orders[i] = false;
        }

        return used_pn;
    }

    //! PN acceleration and adot between two particles
    /*!
      @param[out] a_pn1[6][3]: (3D) Newton and PN acceleration of p1, each PN term are separately saved in array: [Newton, PN1, PN2, PN2.5, PN3, PN3.5]
      @param[out] a_pn2[6][3]: Newton and PN acceleration of p2 for each PN order 
      @param[out] adot_pn1[6][3]: Newton and PN jerk of p1 for each PN order
      @param[out] adot_pn2[6][3]: Newton PN jerk of p2 for each PN order
      @param[out] dspin1[3]: spin change rate of p1
      @param[out] dspin2[3]: spin change rate of p2
      @param[in] p1: particle receiving PN force
      @param[in] p2: particle giving PN force
      @param[in] used_pn_orders: determine which PN orders are calculated: set in a bool array: [PN1, PN2, PN2.5, PN3, PN3.5, SPIN]
      @param[in] calc_adot: if true, calculated jerk, (default: true)

      Return: Maximum PN orders indices [-1: not used; 0: PN1; 1: PN2; 2: PN2.5; 3: PN3; 4: PN3.5; 5: Spin]
    */
    template<class Tp> 
    int calcAccJerkPN(Float a_pn1[][3], Float a_pn2[][3], 
                      Float adot_pn1[][3], Float adot_pn2[][3],
                      Float dspin1[], Float dspin2[],
                      const Tp& p1, const Tp& p2, 
                      bool used_pn_orders[6], const bool calc_adot = true) {

        // used_order maximum
        int used_order_max = -1;

        // interface to particle
        Float m1 = p1.mass;
        Float m2 = p2.mass;
        const auto& x1 = p1.pos;
        const auto& x2 = p2.pos;
        const auto& v1 = p1.vel;
        const auto& v2 = p2.vel;

        const Float& G = gravitational_constant;

        // Final acc and adot for each order of PN
        /*a_pn1   [0 - PN0; 1 - PN1; 2 - PN2; 3 - PN2.5, 4 - PN3, 5 - PN3.5] [3]	for p1
          adot_pn1[0 - PN0; 1 - PN1; 2 - PN2; 3 - PN2.5, 4 - PN3, 5 - PN3.5] [3]	for p1

          a_pn2   [0 - PN0; 1 - PN1; 2 - PN2; 3 - PN2.5, 4 - PN3, 5 - PN3.5] [3]	for p2
          adot_pn2[0 - PN0; 1 - PN1; 2 - PN2; 3 - PN2.5, 4 - PN3, 5 - PN3.5] [3]	for p2
        */
        Float SPIN[3][2];

        bool used_spin = used_pn_orders[5];
        if (used_spin) {
            for(int k=0;k<3;k++) {
                // Currently 3D Spin is not implemented in particle
                //SPIN[k][0] = p1.spin[k];
                //SPIN[k][1] = p2.spin[k];
                SPIN[k][0] = SPIN[k][1] = 0.0; 
            }
        }
        else {
            for(int k=0;k<3;k++) SPIN[k][0] = SPIN[k][1] = 0.0;
        }

        const Float PI = 4.0*atan(1.0);
        const Float PI2 = PI*PI;

        // Speed of light "c" and its powers
        Float c_1 = speed_of_light;

        Float c_2 = c_1*c_1;
        Float c_4 = c_2*c_2;
        Float c_5 = c_4*c_1;
        Float c_6 = c_5*c_1;
        Float c_7 = c_6*c_1;

        // Mass parameters
        Float M   = m1+m2;
        Float eta = m1*m2/(M*M);

        // relative position and velocity
        Float x[3], v[3]; 
        for (int k=0; k<3; k++) {
            x[k] = x1[k] - x2[k];
            v[k] = v1[k] - v2[k];
        }

        Float r2 = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
        Float r  = std::sqrt(r2);
        Float r3 = r2*r;

        Float GMOR = G*M/r;
        Float dv2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
        //Float VWHOLE = std::sqrt(dv2);
        Float rdot     = (x[0]*v[0] + x[1]*v[1] + x[2]*v[2])/r;

        // Newton accelerations ~1/c^0
        Float N[3];
        for (int k=0; k<3; k++) {
            N[k] = x[k]/r;

            a_pn1[0][k] = -G*m2*x[k]/r3;
            a_pn2[0][k] =  G*m1*x[k]/r3;
        }

        Float A[3] = {0.0};
        Float AK2 = 0.0;
        Float BK2 = 0.0;
        // PN1 ~1/c^2
        if(used_pn_orders[0]) {
            Float A1 = 2.0*(2.0+eta)*GMOR-(1.0+3.0*eta)*dv2 +1.5*eta*rdot*rdot;
            Float B1 = 2.0*(2.0-eta)*rdot;

            AK2 = A1/c_2;
            BK2 = B1/c_2;

            for (int k=0; k<3; k++) {
                Float CK2 = AK2*N[k] + BK2*v[k];
                a_pn1[1][k] =  G*CK2/r2*m2;
                a_pn2[1][k] = -G*CK2/r2*m1;
                
                A[k] += GMOR*CK2/r;
            }
            used_order_max = 0;
        }

        // PN2 ~1/c^4
        Float AK4 = 0.0;
        Float BK4 = 0.0;
        if(used_pn_orders[1]) {
            Float A2 = -0.75*(12.0+29.0*eta)*GMOR*GMOR-eta*(3.0-4.0*eta)*dv2*dv2-1.875*eta*(1.0-3.0*eta)*rdot*rdot*rdot*rdot+0.5*eta*(13.0-4.0*eta)*GMOR*dv2+(2.0+25.0*eta+2.0*eta*eta)*GMOR*rdot*rdot+1.5*eta*(3.0-4.0*eta)*dv2*rdot*rdot;
            Float B2 = -0.5*rdot*((4.0+41.0*eta+8.0*eta*eta)*GMOR-eta*(15.0+4.0*eta)*dv2+3.0*eta*(3.0+2.0*eta)*rdot*rdot);

            AK4 = A2/c_4;
            BK4 = B2/c_4;

            for (int k=0; k<3; k++) {
                Float CK4 = AK4*N[k] + BK4*v[k];
                a_pn1[2][k] =  G*CK4/r2*m2;
                a_pn2[2][k] = -G*CK4/r2*m1;

                A[k] += GMOR*CK4/r;
            }
            used_order_max = 1;
        }

        // PN2.5 ~1/c^5
        Float AK5 = 0.0;
        Float BK5 = 0.0;
        if (used_pn_orders[2]) {
            Float A2_5 = 1.6*eta*GMOR*rdot*(17.0*GMOR/3.0+3.0*dv2);
            Float B2_5 = -1.6*eta*GMOR*(3.0*GMOR+dv2);

            AK5 = A2_5/c_5;
            BK5 = B2_5/c_5;

            for (int k=0; k<3; k++) {
                Float CK5 = AK5*N[k] + BK5*v[k];
                a_pn1[3][k] =  G*CK5/r2*m2;
                a_pn2[3][k] = -G*CK5/r2*m1;

                A[k] += GMOR*CK5/r;
            }
            used_order_max = 2;
        }

        // PN3 ~1/c^6
        Float AK6 = 0.0;
        Float BK6 = 0.0;
        if(used_pn_orders[3]) {
            Float A3 = GMOR*GMOR*GMOR*(16.0+(1399.0/12.0-41.0*PI2/16.0)*eta+
                                    71.0*eta*eta/2.0)+eta*(20827.0/840.0+123.0*PI2/64.0-eta*eta)
                *GMOR*GMOR*dv2-(1.0+(22717.0/168.0+615.0*PI2/64.0)*eta+
                                 11.0*eta*eta/8.0-7.0*eta*eta*eta)*GMOR*GMOR*rdot*rdot-
                0.25*eta*(11.0-49.0*eta+52.0*eta*eta)*dv2*dv2*dv2+
                35.0*eta*(1.0-5.0*eta+5.0*eta*eta)*rdot*rdot*rdot*rdot*rdot*rdot/16.0-
                0.25*eta*(75.0+32.0*eta-40.0*eta*eta)*GMOR*dv2*dv2-
                0.5*eta*(158.0-69.0*eta-60.0*eta*eta)*GMOR*rdot*rdot*rdot*rdot+
                eta*(121.0-16.0*eta-20.0*eta*eta)*GMOR*dv2*rdot*rdot+
                3.0*eta*(20.0-79.0*eta+60.0*eta*eta)*dv2*dv2*rdot*rdot/8.0-
                15.0*eta*(4.0-18.0*eta+17.0*eta*eta)*dv2*rdot*rdot*rdot*rdot/8.0;

            Float B3 = rdot*((4.0+(5849.0/840.0+123.0*PI2/32.0)*eta-25.0*eta*eta-
                            8.0*eta*eta*eta)*GMOR*GMOR+eta*(65.0-152.0*eta-48.0*eta*eta)*
                           dv2*dv2/8.0+15.0*eta*(3.0-8.0*eta-2.0*eta*eta)*rdot*rdot*rdot*rdot/8.0+
                           eta*(15.0+27.0*eta+10.0*eta*eta)*GMOR*dv2-eta*(329.0+177.0*eta+
                                                                            108.0*eta*eta)*GMOR*rdot*rdot/6.0-
                           3.0*eta*(16.0-37.0*eta-16.0*eta*eta)*dv2*rdot*rdot/4.0);

            AK6 = A3/c_6;
            BK6 = B3/c_6;
            
            for (int k=0; k<3; k++) {
                Float CK6 = AK6*N[k] + BK6*v[k];
                a_pn1[4][k] =  G*CK6/r2*m2;
                a_pn2[4][k] = -G*CK6/r2*m1;

                A[k] += GMOR*CK6/r;
            }
            used_order_max = 3;
        }

		// PN3.5 ~1/c^7
        Float AK7 = 0.0;
        Float BK7 = 0.0;
        if(used_pn_orders[4]) {
            Float A3_5 = GMOR*eta*(dv2*dv2*(-366.0/35.0-12.0*eta)+dv2*rdot*rdot*(114.0+12.0*eta)-112.0*rdot*rdot*rdot*rdot+GMOR*(dv2*(-692.0/35.0+724.0*eta/15.0)+rdot*rdot*(-294.0/5.0-376.0*eta/5.0)+GMOR*(-3956.0/35.0-184.0*eta/5.0)));
            Float B3_5 = 8.0*eta*GMOR*((1325.0+546.0*eta)*GMOR*GMOR/42.0+(313.0+42.0*eta)*dv2*dv2/28.0+75.0*rdot*rdot*rdot*rdot-(205.0+777.0*eta)*GMOR*dv2/42.0+(205.0+424.0*eta)*GMOR*rdot*rdot/12.0-3.0*(113.0+2.0*eta)*dv2*rdot*rdot/4.0)/5.0;


            AK7 = A3_5/c_7;
            BK7 = B3_5/c_7;
            
            for (int k=0; k<3; k++) {
                Float CK7 = AK7*N[k] + BK7*v[k];
                a_pn1[5][k] =  G*CK7/r2*m2;
                a_pn2[5][k] = -G*CK7/r2*m1;

                A[k] += GMOR*CK7/r;
            }
            used_order_max = 4;
        }

        // Spin accelerations
        // PN accelerations
        Float S1[3], S2[3], KSS[3], KSSIG[3], VDS, VDSIG, NDS, NDSIG, NDV;
        if (used_spin) {
            Float DM = m1 - m2;

            Float XS[3], XA[3];
            for(int k=0;k<3;k++) {
                S1[k] = SPIN[k][0]*G*G*m1*m1/c_1;			// fizikai spin 
                S2[k] = SPIN[k][1]*G*G*m2*m2/c_1;
                KSS[k] = S1[k]+S2[k];
                KSSIG[k] = M*(S2[k]/m2-S1[k]/m1);
                XS[k]  = 0.5*(SPIN[k][0]+SPIN[k][1]);
                XA[k]  = 0.5*(SPIN[k][0]-SPIN[k][1]);
            }

            //NCV crossproduct of N[k] and relative v = N[k]Xv[j]
            Float NCV[3];
            NCV[0] = N[1]*v[2] - N[2]*v[1];  
            NCV[1] = N[2]*v[0] - N[0]*v[2];  
            NCV[2] = N[0]*v[1] - N[1]*v[0];  

            //NCS crossproduct of N[k] and KSS = N[k]XKSS
            Float NCS[3];
            NCS[0] = N[1]*KSS[2] - N[2]*KSS[1];  
            NCS[1] = N[2]*KSS[0] - N[0]*KSS[2];  
            NCS[2] = N[0]*KSS[1] - N[1]*KSS[0];  

            //NCSIG crossproduct of N[k] and KSSIG = N[k]XKSSIG
            Float NCSIG[3];
            NCSIG[0] = N[1]*KSSIG[2] - N[2]*KSSIG[1];  
            NCSIG[1] = N[2]*KSSIG[0] - N[0]*KSSIG[2];  
            NCSIG[2] = N[0]*KSSIG[1] - N[1]*KSSIG[0];  

            //VCS crossproduct of v[k] and KSS = v[k]XKSS
            Float VCS[3];
            VCS[0] = v[1]*KSS[2] - v[2]*KSS[1];  
            VCS[1] = v[2]*KSS[0] - v[0]*KSS[2];  
            VCS[2] = v[0]*KSS[1] - v[1]*KSS[0];  

            //VCSIG crossproduct of v[k] and KSSIG = v[k]XKSSIG
            Float VCSIG[3];
            VCSIG[0] = v[1]*KSSIG[2] - v[2]*KSSIG[1];  
            VCSIG[1] = v[2]*KSSIG[0] - v[0]*KSSIG[2];  
            VCSIG[2] = v[0]*KSSIG[1] - v[1]*KSSIG[0];  

            Float SDNCV = KSS[0]*NCV[0]+KSS[1]*NCV[1]+KSS[2]*NCV[2];

            Float SIGDNCV = KSSIG[0]*NCV[0]+KSSIG[1]*NCV[1]+KSSIG[2]*NCV[2];

            NDV = N[0]*v[0] + N[1]*v[1] + N[2]*v[2];

            Float XS2 = XS[0]*XS[0]+XS[1]*XS[1]+XS[2]*XS[2];
            Float XA2 = XA[0]*XA[0]+XA[1]*XA[1]+XA[2]*XA[2];

            Float NXA = N[0]*XA[0]+N[1]*XA[1]+N[2]*XA[2];
            Float NXS = N[0]*XS[0]+N[1]*XS[1]+N[2]*XS[2];

            VDS = v[0]*KSS[0]+v[1]*KSS[1]+v[2]*KSS[2];
            VDSIG = v[0]*KSSIG[0]+v[1]*KSSIG[1]+v[2]*KSSIG[2];

            NDS = N[0]*KSS[0]+N[1]*KSS[1]+N[2]*KSS[2];
            NDSIG = N[0]*KSSIG[0]+N[1]*KSSIG[1]+N[2]*KSSIG[2];

            Float C1_5[3],C2[3], C2_5[3];
            //Float C1_5D[3],C2D[3], C2_5D[3];

            for(int k=0; k<3; k++) {
                C1_5[k] = (N[k]*(12.0*SDNCV+6.0*DM*SIGDNCV/M)+9.0*NDV*NCS[k]+3.0*DM*NDV*NCSIG[k]/M -7.0*VCS[k]-3.0*DM*VCSIG[k]/M)/r3;
                C2[k] = -GMOR*GMOR*GMOR/r*3.0*eta*(N[k]*(XS2-XA2-5.0*NXS*NXS+5.0*NXA*NXA)+2.0*(XS[k]*NXS-XA[k]*NXA));
                C2_5[k] = (N[k]*(SDNCV*(-30.0*eta*NDV*NDV+24.0*eta*dv2-GMOR*(38.0+25.0*eta))+DM/M*SIGDNCV*(-15.0*eta*NDV*NDV+12.0*eta*dv2-GMOR*(18.0+14.5*eta)))+NDV*v[k]*(SDNCV*(-9.0+9.0*eta)+DM/M* SIGDNCV*(-3.0+6.0*eta))+NCV[k]*(NDV*VDS*(-3.0+3.0*eta)-8.0*GMOR*eta*NDS-DM/M*(4.0*GMOR*eta*NDSIG+3.0*NDV*VDSIG))+NDV*NCS[k]*(-22.5*eta*NDV*NDV+21.0*eta*dv2-GMOR*(25.0+15.0*eta)) +DM/M*NDV*NCSIG[k]*(-15.0*eta*NDV*NDV+12.0*eta*dv2-GMOR*(9.0+8.5*eta))+VCS[k]*(16.5*eta*NDV*NDV+GMOR*(21.0+9.0*eta)-14.0*eta*dv2)+DM/M*VCSIG[k]*(9.0*eta*NDV*NDV-7.0*eta*dv2+GMOR*(9.0+4.5*eta)))/r3;

                A[k] += C1_5[k]/c_2 + C2[k]/c_4 + C2_5[k]/c_4;
            } 

            used_order_max = 5;
        } /* if Spin is used*/   

        if (calc_adot) {
            // PN jerks

            Float AT[3];
            AT[0] = A[0] - GMOR*N[0]/r;		
            AT[1] = A[1] - GMOR*N[1]/r;
            AT[2] = A[2] - GMOR*N[2]/r;

            /*
              AT[0] = A[0];
              AT[1] = A[1];
              AT[2] = A[2];
            */

            Float V2R = dv2/r + AT[0]*N[0] + AT[1]*N[1] + AT[2]*N[2] - rdot*rdot/r;
            Float VA = AT[0]*v[0] + AT[1]*v[1] + AT[2]*v[2];

            // Newton ~1/c^0
            //Float NDOT[3];
            for(int k=0;k<3;k++) {
                //NDOT[k] = (v[k]-N[k]*rdot)/r;
                adot_pn1[0][k] = -G*m2*(v[k]/r3 - 3.0*rdot*x[k]/r2/r2);
                adot_pn2[0][k] =  G*m1*(v[k]/r3 - 3.0*rdot*x[k]/r2/r2);
            }

            /*
            Float NVDOT = NDOT[0]*v[0]+NDOT[1]*v[1]+NDOT[2]*v[2]+N[0]*AT[0]+N[1]*AT[1]+N[2]*AT[2];

            //NDOTCV crossproduct of NDOT[k] and relative v = NDOT[k]Xv[j]
            Float NDOTCV[3];
            NDOTCV[0] = NDOT[1]*v[2] - NDOT[2]*v[1];  
            NDOTCV[1] = NDOT[2]*v[0] - NDOT[0]*v[2];  
            NDOTCV[2] = NDOT[0]*v[1] - NDOT[1]*v[0];  

            //NCA crossproduct of N and AT = N[k]XAT[j]
            Float NCA[3];
            NCA[0] = N[1]*AT[2] - N[2]*AT[1];  
            NCA[1] = N[2]*AT[0] - N[0]*AT[2];  
            NCA[2] = N[0]*AT[1] - N[1]*AT[0];  

            */

            // PN1 ~1/c^2
            if(used_pn_orders[0]) {
                Float A1D = -2.0*(2.0+eta)*GMOR*rdot/r - 2.0*(1.0+3.0*eta)*VA + 3.0*eta*rdot*V2R;
                Float B1D =  2.0*(2.0-eta)*V2R;

                Float ADK2 = A1D/c_2;
                Float BDK2 = B1D/c_2;

                for (int k=0; k<3; k++) {
                    Float CK2 = AK2*N[k] + BK2*v[k];
                    Float CDK2 = ADK2*N[k]+BDK2*v[k];
                    adot_pn1[1][k] =  (-2.0*GMOR*rdot*CK2/r2 + GMOR*CDK2/r + GMOR*(AK2*(v[k]-N[k]*rdot)/r+BK2*AT[k])/r)*m2/M;
                    //adot_pn2[1][k] = -(-2.0*GMOR*rdot*CK2/r2 + GMOR*CDK2/r + GMOR*(AK2*(v[k]-N[k]*rdot)/r+BK2*A[k])/r)*m1/M;
                    adot_pn2[1][k] =  -adot_pn1[1][k]*m1/m2;
                }
            }

            // PN2 ~1/c^4
            if(used_pn_orders[1]) {
                Float A2D =  1.5*(12.0+29.0*eta)*GMOR*GMOR*rdot/r -eta*(3.0-4.0*eta)*4.0*dv2*VA - 7.5*eta*(1.0-3.0*eta)*V2R*rdot*rdot*rdot -0.5*eta*(13.0-4.0*eta)*GMOR*rdot*dv2/r+eta*(13.0-4.0*eta)*GMOR*VA -(2.0+25.0*eta+2.0*eta*eta)*GMOR*rdot*rdot*rdot/r+2.0*(2.0+25.0*eta+2.0*eta*eta)*GMOR*rdot*V2R + 3.0*eta*(3.0-4.0*eta)*VA*rdot*rdot + 3.0*eta*(3.0-4.0*eta)*dv2*rdot*V2R;
                Float B2D = -0.5*V2R*((4.0+41.0*eta+8.0*eta*eta)*GMOR - eta*(15.0+4.0*eta)*dv2+3.0*eta*(3.0+2.0*eta)*rdot*rdot) - 0.5*rdot*(-(4.0+41.0*eta+8.0*eta*eta)*GMOR*rdot/r - 2.0*eta*(15.0+4.0*eta)*VA + 6.0*eta*(3.0+2.0*eta)*rdot*V2R);

                Float ADK4 = A2D/c_4;
                Float BDK4 = B2D/c_4;
                
                for (int k=0; k<3; k++) {
                    Float CK4 = AK4*N[k] + BK4*v[k];
                    Float CDK4 = ADK4*N[k]+BDK4*v[k];
                    adot_pn1[2][k] =  (-2.0*GMOR*rdot*CK4/r2 + GMOR*CDK4/r + GMOR*(AK4*(v[k]-N[k]*rdot)/r+BK4*AT[k])/r)*m2/M;
                    //adot_pn2[2][k] = -G*(-2.0*GMOR*rdot*CK4/r2 + GMOR*CDK4/r + GMOR*(AK4*(v[k]-N[k]*rdot)/r+BK4*AT[k])/r)*m1/M;
                    adot_pn2[2][k] = -adot_pn1[2][k]*m1/m2;
                }
            }

            // PN2.5 ~1/c^5
            if (used_pn_orders[2]) {
                Float A2_5D = -1.6*eta*GMOR*rdot*rdot*(17.0/3.0*GMOR+3.0*dv2)/r +1.6*eta*GMOR*V2R*(17.0/3.0*GMOR+3.0*dv2)+1.6*eta*GMOR*rdot*(-17.0*GMOR*rdot/3.0/r+6.0*VA);
                Float B2_5D = 1.6*eta*GMOR*rdot*(3.0*GMOR+dv2)/r - 1.6*eta*GMOR*(-3.0*GMOR*rdot/r+2.0*VA);

                Float ADK5 = A2_5D/c_5;
                Float BDK5 = B2_5D/c_5;

                for (int k=0; k<3; k++) {
                    Float CK5 = AK5*N[k] + BK5*v[k];
                    Float CDK5 = ADK5*N[k]+BDK5*v[k];
                    adot_pn1[3][k] =  (-2.0*GMOR*rdot*CK5/r2 + GMOR*CDK5/r + GMOR*(AK5*(v[k]-N[k]*rdot)/r+BK5*AT[k])/r)*m2/M;
                    //adot_pn2[3][k] = -(-2.0*GMOR*rdot*CK5/r2 + GMOR*CDK5/r + GMOR*(AK5*(v[k]-N[k]*rdot)/r+BK5*AT[k])/r)*m1/M;
                    adot_pn2[3][k] = -adot_pn1[3][k]*m1/m2;
                }
            }

            // PN3 ~1/c^6
            if (used_pn_orders[3]) {
                Float A3D =  6.0*eta*rdot*rdot*rdot*rdot*rdot*V2R*(35.0-175.0*eta+175.0*eta*eta)/16.0 + eta*(4.0*rdot*rdot*rdot*V2R*dv2 + 2.0*rdot*rdot*rdot*rdot*VA)*(-15.0+135.0*eta/2.0-255.0*eta*eta/4.0)/2.0 + eta*(2.0*rdot*V2R*dv2*dv2+4.0*rdot*rdot*dv2*VA)/2.0*(15.0-237.0*eta/2.0+45.0*eta*eta) + 6.0*dv2*dv2*VA*eta*(-11.0/4.0-49.0*eta/4.0-13.0*eta*eta) + GMOR*(4.0*rdot*rdot*rdot*V2R*eta*(-79.0+69.0/2.0*eta+30.0*eta*eta) + eta*(2.0*rdot*V2R*dv2+2.0*rdot*rdot*VA)*(121.0-16.0*eta-20.0*eta*eta)+4.0*dv2*VA*eta*(-75.0/4.0-8.0*eta+10.0*eta*eta)) - GMOR*rdot*((-79.0+69.0*eta/2.0+30.0*eta*eta)*rdot*rdot*rdot*rdot*eta+eta*rdot*rdot*dv2*(121.0-16.0*eta-20.0*eta*eta)+eta*dv2*dv2*(-75.0/4.0-8.0*eta+10.0*eta*eta))/r - 2.0*GMOR*GMOR*rdot*(rdot*rdot*((-1.0-615.0*PI2*eta/64.0)-22717.0*eta/168.0-11.0*eta*eta/8.0+7.0*eta*eta*eta)+eta*dv2*((20827.0/840.0+123.0*PI2/64.0)-eta*eta))/r + GMOR*GMOR*(2.0*rdot*V2R*((-1.0-615*PI2*eta/64.0)-22717.0*eta/168.0-11.0*eta*eta/8.0+7*eta*eta*eta)+2.0*eta*VA*((20827.0/840.0 +123.0*PI2/64.0)-eta*eta)) - 3.0*GMOR*GMOR*GMOR*rdot*(16.0+(1399.0/12.0-41.0*PI2/16.0)*eta+71.0*eta*eta/2.0)/r;
                Float B3D = 75.0*rdot*rdot*rdot*rdot*V2R*eta*(3.0/8.0-eta-.25*eta*eta)+eta*(3.0*rdot*rdot*V2R*dv2+2.0*rdot*rdot*rdot*VA)*(-12.0+111.0*eta/4.0+12.0*eta*eta)+eta*(V2R*dv2*dv2+4.0*rdot*dv2*VA)*(65.0/8.0-19.0*eta-6.0*eta*eta)-GMOR*rdot*(rdot*rdot*rdot*eta*(-329.0/6.0-59.0*eta/2.0-18.0*eta*eta)+rdot*dv2*eta*(15.0+27.0*eta+10.0*eta*eta))/r+GMOR*(3.0*rdot*rdot*V2R*eta*(-329.0/6.0-59.0*eta/2.0-18.0*eta*eta)+eta*(V2R*dv2+2.0*rdot*VA)*(15.0+27.0*eta+10.0*eta*eta))-2.0*GMOR*GMOR*rdot*(rdot*((4.0+123.0*PI2*eta/32.0)+5849.0*eta/840.0-25.0*eta*eta-8.0*eta*eta*eta))/r+GMOR*GMOR*(V2R*((4.0+123.0*PI2*eta/32.0)+5849.0/840.0*eta-25.0*eta*eta-8.0*eta*eta*eta));

                Float ADK6 = A3D/c_6;
                Float BDK6 = B3D/c_6;
                
                for (int k=0; k<3; k++) {
                    Float CK6 = AK6*N[k] + BK6*v[k];
                    Float CDK6= ADK6*N[k]+BDK6*v[k];
                    adot_pn1[4][k] =  (-2.0*GMOR*rdot*CK6/r2 + GMOR*CDK6/r + GMOR*(AK6*(v[k]-N[k]*rdot)/r+BK6*AT[k])/r)*m2/M;
                    //adot_pn2[4][k] = -(-2.0*GMOR*rdot*CK6/r2 + GMOR*CDK6/r + GMOR*(AK6*(v[k]-N[k]*rdot)/r+BK6*AT[k])/r)*m1/M;
                    adot_pn2[4][k] = -adot_pn1[4][k]*m1/m2;
                }

            }

            // PN3.5 ~1/c^7
            if (used_pn_orders[4]) {
                Float A3_5D = GMOR*eta*(-rdot*(dv2*dv2*(-366.0/35.0-12.0*eta)+dv2*rdot*rdot*(114.0+12.0*eta)+rdot*rdot*rdot*rdot*(-112.0))/r+4.0*dv2*VA*(-366.0/35.0-12.0*eta)+2.0*(VA*rdot*rdot+rdot*V2R*dv2)*(114.0+12.0*eta)+4.0*rdot*rdot*rdot*V2R*(-112.0)+GMOR*(2.0*VA*(-692.0/35.0+724.0*eta/15.0)+2.0*rdot*V2R*(-294.0/5.0-376.0*eta/5.0)-2.0*rdot*(dv2*(-692.0/35.0+724.0*eta/15.0)+rdot*rdot*(-294.0/5.0-376.0*eta/5.0))/r-3.0*GMOR*rdot*(-3956.0/35.0-184.0*eta/5.0)/r));
                Float B3_5D = GMOR*eta*(4.0*dv2*VA*(626.0/35.0+12.0*eta/5.0)+2.0*(VA*rdot*rdot+dv2*rdot*V2R)*(-678.0/5.0-12.0*eta/5.0)+4.0*rdot*rdot*rdot*V2R*120.0-rdot*(dv2*dv2*(626.0/35.0+12.0*eta/5.0)+dv2*rdot*rdot*(-678.0/5.0-12.0*eta/5.0)+120.0*rdot*rdot*rdot*rdot)/r+GMOR*(2.0*VA*(-164.0/21.0-148.0*eta/5.0)+2*rdot*V2R*(82.0/3.0+848.0*eta/15.0)-2.0*rdot*(dv2*(-164.0/21-148.0*eta/5.0)+rdot*rdot*(82.0/3.0+848.0*eta/15.0))/r-3.0*GMOR*rdot*(1060.0/21.0+104.0*eta/5.0)/r));

                Float ADK7 = A3_5D/c_7;
                Float BDK7 = B3_5D/c_7;
                    
                for (int k=0; k<3; k++) {
                    Float CK7 = AK7*N[k] + BK7*v[k];
                    Float CDK7 = ADK7*N[k]+BDK7*v[k];
                    adot_pn1[5][k] =  (-2.0*GMOR*rdot*CK7/r2 + GMOR*CDK7/r + GMOR*(AK7*(v[k]-N[k]*rdot)/r+BK7*AT[k])/r)*m2/M;
                    //adot_pn2[5][k] = -(-2.0*GMOR*rdot*CK7/r2 + GMOR*CDK7/r + GMOR*(AK7*(v[k]-N[k]*rdot)/r+BK7*AT[k])/r)*m1/M;                    
                    adot_pn2[5][k] = -adot_pn1[5][k]*m1/m2;

                }
            }

        }

        // Spin
        if (used_spin) {
            Float L[3];
            //L crossproduct of x[k] and relative v = x[k]Xv[j]
            L[0] = x[1]*v[2] - x[2]*v[1];  
            L[1] = x[2]*v[0] - x[0]*v[2];  
            L[2] = x[0]*v[1] - x[1]*v[0];  

            Float LABS = std::sqrt(L[0]*L[0]+L[1]*L[1]+L[2]*L[2]);

            Float LU[3];
            LU[0] = L[0]/LABS;
            LU[1] = L[1]/LABS;
            LU[2] = L[2]/LABS;

            Float S1DLU = S1[0]*LU[0]+S1[1]*LU[1]+S1[2]*LU[2];
            Float S2DLU = S2[0]*LU[0]+S2[1]*LU[1]+S2[2]*LU[2];

            Float SU1[3], SV1[3], SS1[3], SS2[3], SU2[3], SV2[3]; 

            Float DM = m1 - m2;

            for (int k=0; k<3; k++) {
                SU1[k] = GMOR*eta*(N[k]*(-4.0*VDS-2.0*DM/M*VDSIG)+ v[k]*(3.0*NDS+DM/M*NDSIG)+NDV*(2.0*KSS[k]+DM/M*KSSIG[k])) /r;
                SV1[k] = GMOR*(N[k]*(VDSIG*(-2.0+4.0*eta)-2.0*DM/M*VDS)+ v[k]*(NDSIG*(1.0-eta)+DM/M*NDS)+NDV*(KSSIG[k]*(1.0- 2.0*eta)+ DM/M*KSS[k]))/r;

                SS1[k] = 0.5*(L[k]*(4.0+3.0*(m2/m1))+ (S2[k]-3.0*S2DLU*LU[k]))/r3;
                SS2[k] = 0.5*(L[k]*(4.0+3.0*(m1/m2))+ (S1[k]-3.0*S1DLU*LU[k]))/r3;

                SU2[k] = GMOR*eta/r*(N[k]*(VDS*(-2.0*dv2+3.0*NDV*NDV- 6.0*eta*NDV*NDV+7.0*GMOR-8.0*eta*GMOR)-14.0*GMOR*NDS*NDV+ DM/M*VDSIG*eta*(-3.0*NDV*NDV-4.0*GMOR)+DM/M*GMOR*NDSIG*NDV* (2.0-eta/2.))+v[k]*(NDS*(2.0*dv2-4.0*eta*dv2-3.0*NDV* NDV+7.5*eta*NDV*NDV+4.0*GMOR-6.0*eta*GMOR)+VDS*NDV*(2.0- 6.0*eta)+ DM/M*NDSIG*(-1.5*eta*dv2+3.0*eta*NDV*NDV-GMOR-3.5*eta* GMOR)-3.0*DM/M*VDSIG*NDV*eta)+KSS[k]*NDV*(dv2-2.0*eta* dv2-1.5*NDV*NDV+3.0*eta*NDV*NDV-GMOR+2.0*eta*GMOR)+ DM/M*KSSIG[k]*NDV*(-eta*dv2+1.5*eta*NDV*NDV+ (eta-1.)*GMOR));
                SV2[k] = GMOR/r*(N[k]*(VDSIG*eta*(-2.0*dv2+6.0*eta*NDV* NDV+(3.0+8.0*eta)*GMOR)+GMOR*NDSIG*NDV*(2.0-22.5*eta+2.0* eta*eta)+ DM/M*VDS*eta*(-3.0*NDV*NDV-4.0*GMOR)+DM/M*GMOR*NDS*NDV*(2.0- 0.5*eta))+v[k]*(NDSIG*(0.5*eta*dv2+2.0*eta*eta*dv2- 4.5*eta*eta*NDV*NDV+(4.5*eta-1.0+8.0*eta*eta)*GMOR)+VDSIG*NDV* eta*(6.0*eta-1.)-3.0*DM/M*VDS*NDV*eta+DM/M*NDS*(-1.5* eta*dv2+ 3.0*eta*NDV*NDV-(1.0+3.5*eta)*GMOR))+KSSIG[k]*NDV*(2.0*eta*eta* dv2-3.0*eta*eta*NDV*NDV+(-1.0+4.0*eta-2.0*eta*eta)*GMOR)+ DM/M*KSS[k]*NDV*(-eta*dv2+1.5*eta*NDV*NDV+(-1.0+eta)* GMOR));
            }

            //SS1 crossproduct of SS1 and S1 = SS1[k]XS1[j]
            Float SS1aux[3],SS2aux[3],SU[3],SV[3];
            SS1aux[0] =   SS1[1]*S1[2] - SS1[2]*S1[1];  
            SS1aux[1] =   SS1[2]*S1[0] - SS1[0]*S1[2];  
            SS1aux[2] =   SS1[0]*S1[1] - SS1[1]*S1[0];  

            SS1[0] = SS1aux[0];  
            SS1[1] = SS1aux[1];  
            SS1[2] = SS1aux[2];  

            //SS2 crossproduct of SS2 and S2 = SS2[k]XS2[j]
            SS2aux[0] =   SS2[1]*S2[2] - SS2[2]*S2[1];  
            SS2aux[1] =   SS2[2]*S2[0] - SS2[0]*S2[2];  
            SS2aux[2] =   SS2[0]*S2[1] - SS2[1]*S2[0];  

            SS2[0] = SS2aux[0];  
            SS2[1] = SS2aux[1];  
            SS2[2] = SS2aux[2];  

            // new values of the spins, returned back to the main program... 
            if ((dspin1==NULL) || (dspin2==NULL)) {
                std::cerr<<"Error: spin array is NULL, cannot save spin data!\n";
                abort();
            }

            for(int k=0; k<3; k++) {
                SU[k] = SU1[k]/c_2 + SU2[k]/c_4 + (SS1[k] + SS2[k])/c_2;
                SV[k] = SV1[k]/c_2 + SV2[k]/c_4+M*(SS2[k]/m2-SS1[k]/ m1)/c_2;

                //KSS[k] = KSS[k] + SU[k]*dt_bh;						// integrate for dt_bh timestep 
                //KSSIG[k] = KSSIG[k] + SV[k]*dt_bh;

                //dspin1[k] = m1*(M*KSS[k]-m2*KSSIG[k])/M/M/m1/m1*c_1 ;       
                //dspin2[k] = m2*(M*KSS[k]+m1*KSSIG[k])/M/M/m2/m2*c_1;
                
                Float mmc = c_1/(G*G*M*M);
                dspin1[k] = mmc*(M*SU[k]-m2*SV[k])/m1;
                dspin2[k] = mmc*(M*SU[k]+m1*SV[k])/m2;
                //XAD[k] = 0.5/(M*M*m1*m2)*(-SU[k]*M*DM-SV[k]*(m1*m1+m2*m2));
                //XSD[k] = 0.5/(M*M*m1*m2)*(SU[k]*M*M+SV[k]*(m1*m1-m2*m2));

            }

            /*
            //NDOTCS crossproduct of NDOT and KSS = NDOT[k]XKSS[j]
            Float NDOTCS[3], NCSU[3], NDOTCSIG[3], NCSV[3], ACS[3], VCSU[3], ACSIG[3], VCSV[3];
            NDOTCS[0] =   NDOT[1]*KSS[2] - NDOT[2]*KSS[1];  
            NDOTCS[1] =   NDOT[2]*KSS[0] - NDOT[0]*KSS[2];  
            NDOTCS[2] =   NDOT[0]*KSS[1] - NDOT[1]*KSS[0];  
            //NCSU crossproduct of N and SU = N[k]XSU[j]
            NCSU[0] =   N[1]*SU[2] - N[2]*SU[1];  
            NCSU[1] =   N[2]*SU[0] - N[0]*SU[2];  
            NCSU[2] =   N[0]*SU[1] - N[1]*SU[0];  
            //NDOTCSIG crossproduct of NDOT and KSSIG = NDOT[k]XKSSIG[j]
            NDOTCSIG[0] =   NDOT[1]*KSSIG[2] - NDOT[2]*KSSIG[1];  
            NDOTCSIG[1] =   NDOT[2]*KSSIG[0] - NDOT[0]*KSSIG[2];  
            NDOTCSIG[2] =   NDOT[0]*KSSIG[1] - NDOT[1]*KSSIG[0];  
            //NCSV crossproduct of N and SV = N[k]XSV[j]
            NCSV[0] =   N[1]*SV[2] - N[2]*SV[1];  
            NCSV[1] =   N[2]*SV[0] - N[0]*SV[2];  
            NCSV[2] =   N[0]*SV[1] - N[1]*SV[0];  
            //ACS crossproduct of AT and KSS = AT[k]XKSS[j]
            ACS[0] =   AT[1]*KSS[2] - AT[2]*KSS[1];  
            ACS[1] =   AT[2]*KSS[0] - AT[0]*KSS[2];  
            ACS[2] =   AT[0]*KSS[1] - AT[1]*KSS[0];  
            //VCSU crossproduct of relative v and SU = v[k]XSU[j]
            VCSU[0] =   v[1]*SU[2] - v[2]*SU[1];  
            VCSU[1] =   v[2]*SU[0] - v[0]*SU[2];  
            VCSU[2] =   v[0]*SU[1] - v[1]*SU[0];  
            //ACSIG crossproduct of AT and KSSIG = AT[k]XKSSIG[j]
            ACSIG[0] =   AT[1]*KSSIG[2] - AT[2]*KSSIG[1];  
            ACSIG[1] =   AT[2]*KSSIG[0] - AT[0]*KSSIG[2];  
            ACSIG[2] =   AT[0]*KSSIG[1] - AT[1]*KSSIG[0];  
            //VCSV crossproduct of relative v and SV = v[k]XSV[j]
            VCSV[0] =   v[1]*SV[2] - v[2]*SV[1];  
            VCSV[1] =   v[2]*SV[0] - v[0]*SV[2];  
            VCSV[2] =   v[0]*SV[1] - v[1]*SV[0];  

            Float SNVDOT = SU[0]*NCV[0]+SU[1]*NCV[1]+SU[2]*NCV[2]+ KSS[0]*NDOTCV[0]+KSS[1]*NDOTCV[1]+KSS[2]*NDOTCV[2]+ KSS[0]*NCA[0]+KSS[1]*NCA[1]+KSS[2]*NCA[2];

            Float SIGNVDOT = SV[0]*NCV[0]+SV[1]*NCV[1]+SV[2]*NCV[2]+ KSSIG[0]*NDOTCV[0]+KSSIG[1]*NDOTCV[1]+KSSIG[2]*NDOTCV[2]+ KSSIG[0]*NCA[0]+KSSIG[1]*NCA[1]+KSSIG[2]*NCA[2];
                
            Float NSDOT = NDOT[0]*KSS[0]+NDOT[1]*KSS[1]+NDOT[2]*KSS[2]+ N[0]*SU[0]+N[1]*SU[1]+N[2]*SU[2];
            Float NSIGDOT = NDOT[0]*KSSIG[0]+NDOT[1]*KSSIG[1]+NDOT[2]*KSSIG[2]+ N[0]*SV[0]+N[1]*SV[1]+N[2]*SV[2];
            Float VSDOT = AT[0]*KSS[0]+AT[1]*KSS[1]+AT[2]*KSS[2]+ v[0]*SU[0]+v[1]*SU[1]+v[2]*SU[2];
            Float VSIGDOT = AT[0]*KSSIG[0]+AT[1]*KSSIG[1]+AT[2]*KSSIG[2]+ v[0]*SV[0]+v[1]*SV[1]+v[2]*SV[2];

            Float NXSDOT = NDOT[0]*XS[0]+NDOT[1]*XS[1]+NDOT[2]*XS[2]+ N[0]*XSD[0]+N[1]*XSD[1]+N[2]*XSD[2];
            Float NXADOT = NDOT[0]*XA[0]+NDOT[1]*XA[1]+NDOT[2]*XA[2]+ N[0]*XAD[0]+N[1]*XAD[1]+N[2]*XAD[2];

            Float C1_5D[3], C2D[3], C2_5D[3];
            for(int k=0; k<3; k++) { 
            C1_5D[k] = -3.0*rdot/r*C1_5[k]+(NDOT[k]*(12.0*SDNCV+6.0*DM/M* SIGDNCV)+N[k]*(12.0*SNVDOT+6.0*DM/M*SIGNVDOT)+9.0*NVDOT* NCS[k]+9.0*NDV*(NDOTCS[k]+NCSU[k])+3.0*DM/M*(NVDOT*NCSIG[k]+ NDV*(NDOTCSIG[k]+NCSV[k]))-7.0*(ACS[k]+VCSU[k])-3.0*DM/M* (ACSIG[k]+VCSV[k]))/(r3);
            C2D[k] = -4.0*rdot/r*C2[k]-GMOR*GMOR*GMOR*3.0*eta/r*(NDOT[k]* (XS2-XA2-5.0*NXS*NXS+5.0*NXA*NXA)+N[k]*(2.0*(XS[0]*XSD[0]+ XS[1]*XSD[1]+XS[2]*XSD[2]-XA[0]*XAD[0]-XA[1]*XAD[1]- XA[2]*XAD[2])-10.0*NXS*NXSDOT+10.0*NXA*NXADOT)+2.0*(XSD[k]* NXS+XS[k]*NXSDOT-XAD[k]*NXA-XA[k]*NXADOT));
            C2_5D[k] = -3.0*rdot/r*C2_5[k]+(NDOT[k]*(SDNCV*(-30.0*eta* NDV*NDV+24.0*eta*dv2-GMOR*(38.0+25.0*eta))+DM/M*SIGDNCV* (-15.0*eta*NDV*NDV+12.0*eta*dv2-GMOR*(18.0+14.5*eta)))+ N[k]*(SNVDOT*(-30.0*eta*NDV*NDV+24.0*eta*dv2-GMOR* (38.0+25.0*eta))+SDNCV*(-60.0*eta*NDV*NVDOT+48.0*eta*VA+ GMOR*rdot/r*(38.0+25.0*eta))+DM/M*SIGNVDOT*(-15.0*eta*NDV* NDV+12.0*eta*dv2-GMOR*(18.0+14.5*eta))+DM/M*SIGDNCV* (-30.0*eta*NDV*NVDOT+24.0*eta*VA+GMOR*rdot/r*(18.0+14.5*eta)))+ (NVDOT*v[k]+NDV*AT[k])*(SDNCV*(-9.0+9.0*eta)+DM/M*SIGDNCV* (-3.0+6.0*eta))+NDV*v[k]*(SNVDOT*(-9.0+9.0*eta)+DM/M* SIGNVDOT*(-3.0+6.0*eta))+(NDOTCV[k]+NCA[k])*(NDV*VDS*(-3.0+ 3.0*eta)-8.0*GMOR*eta*NDS-DM/M*(4.0*GMOR*eta*NDSIG+3.0*NDV*VDSIG) )+NCV[k]*((NVDOT*VDS+NDV*VSDOT)*(-3.0+3.0*eta)-8.0*eta*GMOR* (NSDOT-rdot/r*NDS)-DM/M*(4.0*eta*GMOR*(NSIGDOT-rdot/r*NDSIG)+ 3.0*(NVDOT*VDSIG+NDV*VSIGDOT)))+(NVDOT*NCS[k]+NDV* (NDOTCS[k]+NCSU[k]))*(-22.5*eta*NDV*NDV+21.0*eta*dv2- GMOR*(25.0+15.0*eta))+NDV*NCS[k]*(-45.0*eta*NDV*NVDOT+42.0*eta* VA+GMOR*rdot/r*(25.0+15.0*eta))+DM/M*(NVDOT*NCSIG[k]+NDV* (NDOTCSIG[k]+NCSV[k]))*(-15.0*eta*NDV*NDV+12.0*eta*dv2- GMOR*(9.0+8.5*eta))+DM/M*NDV*NCSIG[k]*(-30.0*eta*NDV*NVDOT+ 24.0*eta*VA+GMOR*rdot/r*(9.0+8.5*eta))+(ACS[k]+VCSU[k])* (16.5*eta*NDV*NDV+GMOR*(21.0+9.0*eta)-14.0*eta*dv2)+ VCS[k]*(33.0*eta*NDV*NVDOT-GMOR*rdot/r*(21.0+9.0*eta)- 28.0*eta*VA)+DM/M*(ACSIG[k]+VCSV[k])*(9.0*eta*NDV*NDV- 7.0*eta*dv2+GMOR*(9.0+4.5*eta))+DM/M*VCSIG[k]*(18.0* eta*NDV*NVDOT-14.0*eta*VA-GMOR*rdot/r*(9.0+4.5*eta)))/ (r3);
            }
            */

        } // spin

        /*
          Float ADK = ADK2+ADK4+ADK5+ADK6+ADK7;
          Float BDK = BDK2+BDK4+BDK5+BDK6+BDK7;

          Float KSAK = AK2+AK4+AK5+AK6+AK7;
          Float KSBK = BK2+BK4+BK5+BK6+BK7;

          Float AD[3] = {0.0};
          for (int k=0; k<3; k++) {
          AD[k] = -2.0*GMOR*rdot*(KSAK*N[k]+KSBK*v[k])/r2 + GMOR*(ADK*N[k]+BDK*v[k])/r + GMOR*(KSAK*(v[k]-N[k]*rdot)/r+KSBK*AT[k])/r + C1_5D[k]/c_2 + C2D[k]/c_4 +C2_5D[k]/c_4;
          }
        */


        // Check RS_DIST conditions !!!
        //Float RS_DIST = 4.0*(2.0*m1/c_2 + 2.0*m2/c_2);

        return used_order_max;

    }

    //! Sum different orders of PN acc and jerk to acc and adot
    /*!
      @param[out] acc[3]: acceleration, to be added (not reset to zero)
      @param[out] adot[3]: jerk, to be added (not reset to zero)
      @param[in] a_pn[6][3]: (3D) Newton and PN acceleration of particle, each PN term are separately saved in array: [Newton, PN1, PN2, PN2.5, PN3, PN3.5]
      @param[in] adot_pn[6][3]: Newton and PN jerk of particle for each PN order
     */
    void sumAccJerkPN(Float acc[3], Float adot[3], Float a_pn[][3], Float adot_pn[][3]) {
        for (int i=0; i<6; i++) {
            for (int k=0; k<3; k++) {
                acc[k] += a_pn[i][k];
                adot[k] += adot_pn[i][k];
            }
        }
    }

    //! write class data to file with binary format
    /*! @param[in] _fp: FILE type file for output
     */
    void writeBinary(FILE *_fp) const {
        fwrite(this, sizeof(*this),1,_fp);
    }

    //! read class data to file with binary format
    /*! @param[in] _fp: FILE type file for reading
     */
    void readBinary(FILE *_fin) {
        size_t rcount = fread(this, sizeof(*this), 1, _fin);
        if (rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            abort();
        }
    }    

};
