#pragma once

#include <array>
#include <tuple>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <random>
#include <algorithm>  
#include "../parallel-random/rand.hpp"

typedef double Float;

//! Gravitational wave recoil kick calculater
class GWKick{
public:    
    Float vscale;          ///> velocity scale (km/s to pc/Myr)

    GWKick(): vscale(-1) {}

    //! check whether parameters values are correct
    /*! \return true: all correct
     */
    bool checkParams() {
        assert(vscale>0.0);
        return true;
    }        

    // base calculation

    // Function to calculate the norm of a vector
    /*! @param[in] vec: vector
        \return the norm of the vector
     */
    Float norm(const std::array<Float, 3>& vec) {
        return std::sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
    }
    
    // Function to normalize a vector
    /*! @param[in] vec: vector
        \return the normalized vector
     */
    std::array<Float,3> normalize(const std::array<Float,3>& vec){
        return {vec[0]/norm(vec),vec[1]/norm(vec),vec[2]/norm(vec)};
    }

    // Function to calculate the dot product of two vectors
    /*! @param[in] a: vector a
        @param[in] b: vector b
        \return the dot product of the two vectors
     */
    Float dot(const std::array<Float, 3>& a, const std::array<Float, 3>& b) {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    // Function to calculate the projuction  of a to b
    /*! @param[in] a: vector a
        @param[in] b: vector b
        \return the proection of a to b
     */
    Float dot_p(const std::array<Float, 3>&a, const std::array<Float, 3>&b) {
        return dot(a, b) / norm(b);
    }

    // Function to calculate the cross product of two vectors
    /*! @param[in] a: vector a
        @param[in] b: vector b
        \return the cross product of the two vectors
     */
    std::array<Float, 3> cross(const std::array<Float, 3>& a, const std::array<Float, 3>& b) {
        return std::array<Float,3>{ a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0] };
    }

    // Function to calculate the multiplication product of a vector and a scalar
    /*! @param[in] vec: vector
        @param[in] scalar: scalar
        \return the multiplication product of a vector and a scalar
     */
    std::array<Float, 3> multiply(const std::array<Float, 3>& vec, Float scalar) {
        return { vec[0] * scalar, vec[1] * scalar, vec[2] * scalar };
    }
    
    // Function to calculate the subtraction of two vectors
    /*! @param[in] a: vector a
        @param[in] b: vector b
        \return the subtraction of the two vectors 
     */
    std::array<Float, 3> subtract(const std::array<Float, 3>& a, const std::array<Float, 3>& b) {
        return { a[0] - b[0], a[1] - b[1], a[2] - b[2] };
    }

    // Function to calculate the addition of two vectors
    /*! @param[in] a: vector a
        @param[in] b: vector b
        \return the addition of the two vectors
     */
    std::array<Float, 3> add(const std::array<Float, 3>& a, const std::array<Float, 3>& b) {
        return { a[0] + b[0], a[1] + b[1], a[2] + b[2] };
    }

    // Function to calculate the division of a vector and a scalar
    /*! @param[in] vec: vector
        @param[in] scalar: scalar
        \return the division of a vector and a scalar
     */
    std::array<Float, 3> divide(const std::array<Float, 3>& vec, Float scalar) {
        return { vec[0] / scalar, vec[1] / scalar, vec[2] / scalar };
    }

    // Function to calculate the matmultiplication product of a matrix and a vector
    /*! @param[in] mat: matrix
        @param[in] vec: vector
        \return the matmultiplication product of a matrix and a vector
     */
    std::array<Float, 3> matmul(const std::array<std::array<Float, 3>, 3>& mat, const std::array<Float, 3>& vec) {
        std::array<Float, 3> result = { 0.0, 0.0, 0.0 };
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                result[i] += mat[i][j] * vec[j];
            }
        }
        return result;
    }

    std::array<Float, 3> randomVectorWithMagnitude(const Float& magnitude) {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_real_distribution<Float> dist(0.0, 1.0);

        Float theta = 2 * M_PI * dist(gen); // 随机角度在 0 到 2π 之间
        Float phi = std::acos(2 * dist(gen) - 1); // 随机角度在 0 到 π 之间

        Float x = magnitude * std::sin(phi) * std::cos(theta);
        Float y = magnitude * std::sin(phi) * std::sin(theta);
        Float z = magnitude * std::cos(phi);

        return {x, y, z};
    }
    
    //Function about the GW kick
    
    // Function to change the axis of the GW kick to the axis of petar
    /*! @param[in] L: orbital angular momentum
        @param[in] dr: separation vector
        @param[in] vkick: kick velocity: [vx, vy, vz]
        \return the kick velocity in the axis of petar
    */
    std::array<Float, 3> calcInverterAxisL(
        const std::array<Float, 3>& L,
        const std::array<Float, 3>& dr,
        const std::array<Float, 3>& vkick) {

        Float drm = norm(dr);
        Float Lm = norm(L);

        // Define unit vectors
        std::array<Float, 3> e1 = divide(dr, drm);
        std::array<Float, 3> e3 = divide(L, Lm);
        std::array<Float, 3> e2 = cross(e3, e1);

        // Construct matrix E
        std::array<std::array<Float, 3>, 3> E = { e1, e2, e3 };

        // Multiply matrix E by vector vkick
        std::array<Float, 3> vkick_new = matmul(E, vkick);

        return vkick_new;
    }

    // Function to change the axis of petar to the axis of the GW kick
    /*! @param[in] Chi1: spin of the primary body
        @param[in] Chi2: spin of the secondary body
        @param[in] L: orbital angular momentum
        @param[in] dr: separation vector
        \return the spin of the primary and secondary body in the GWkick axis
    */
    std::tuple<std::array<Float, 3>, std::array<Float, 3>> calcAxisL(
        const std::array<Float, 3>& L,
        const std::array<Float, 3>& dr,
        const std::array<Float, 3>& Chi1,
        const std::array<Float, 3>& Chi2) {

        Float drm = norm(dr);
        Float Lm = norm(L);

        // Define unit vectors
        std::array<Float, 3> e1 = divide(dr, drm);
        std::array<Float, 3> e3 = divide(L, Lm);
        std::array<Float, 3> e2 = cross(e3, e1);

        std::array<Float, 3> Chi1_new = { dot_p(Chi1, e1), dot_p(Chi1, e2), dot_p(Chi1, e3) };
        std::array<Float, 3> Chi2_new = { dot_p(Chi2, e1), dot_p(Chi2, e2), dot_p(Chi2, e3) };

        return std::make_tuple(Chi1_new, Chi2_new);
    }


    // Function to calculate the delta which  is (q*chi2 - chi1)/(1+q) in petar axis
    /*! @param[in] Chi1: spin of the primary body
        @param[in] Chi2: spin of the secondary body
        @param[in] L: orbital angular momentum
        @param[in] dr: separation vector
        @param[in] q: mass ratio
        \return the delta which  is (q*chi2 - chi1)/(1+q) in petar axis
     */
    std::array<Float, 3> calcDelta(const std::array<Float, 3>& Chi1, 
        const std::array<Float, 3>& Chi2, 
        const std::array<Float, 3>& L,
        const std::array<Float, 3>& dr,
        const Float& q) {
        
        std::array<Float, 3> Chi1_new = Chi1;
        std::array<Float, 3> Chi2_new = Chi2;
        std::tie(Chi1_new, Chi2_new) = calcAxisL(L, dr, Chi1, Chi2);
        
        Float chi1 = norm(Chi1_new);
        Float chi2 = norm(Chi2_new);
        std::array<Float, 3> hatS1 = normalize(Chi1_new);
        std::array<Float, 3> hatS2 = normalize(Chi2_new);

        std::array<Float, 3> Delta = subtract(multiply(hatS2, q * chi2), multiply(hatS1, chi1));
        Delta = multiply(Delta, -1.0 / (1.0 + q));

        return Delta;
    }

    // Function to generate the chit which is (q**2*chi2 + chi1)/(1+q)**2 in petar axis
    /*! @param[in] _Chi1: spin of the primary body
        @param[in] _Chi2: spin of the secondary body
        @param[in] _L: orbital angular momentum
        @param[in] _dr: separation vector
        @param[in] q: mass ratio
        \return the chit which is (q**2*chi2 + chi1)/(1+q)**2 in petar axis
     */
    std::array<Float, 3> calcChit(const std::array<Float, 3>& Chi1, 
        const std::array<Float, 3>& Chi2, 
        const std::array<Float, 3>& L,
        const std::array<Float, 3>& dr,
        const Float& q) {
        
        std::array<Float, 3> Chi1_new = Chi1;
        std::array<Float, 3> Chi2_new = Chi2;
        std::tie(Chi1_new, Chi2_new) = calcAxisL(L, dr, Chi1, Chi2);
        
        Float chi1 = norm(Chi1_new);
        Float chi2 = norm(Chi2_new);
        std::array<Float, 3> hatS1 = normalize(Chi1_new);
        std::array<Float, 3> hatS2 = normalize(Chi2_new);

        std::array<Float, 3> chit = add(multiply(hatS2, q * q * chi2), multiply(hatS1, chi1));
        chit = multiply(chit, 1.0 / pow(1.0 + q, 2));

        return chit;
    }

    // Function to generate a uniformly distributed point inside a sphere
    /*! @param[in] _radius: radius of the sphere
        \return a point inside the sphere
     */
    std::array<Float, 3> uniformPointsInsideSphere(const Float& radius) {
        Float x, y, z;
        while (true) {
            // Generate a random point within the bounding cube
            x = (-2*rand_f64() + 1)*radius;
            y = (-2*rand_f64() + 1)*radius;
            z = (-2*rand_f64() + 1)*radius;

            // Check if the point is inside the sphere
            if (std::sqrt(x * x + y * y + z * z) <= 0.8) {
                std::array<Float, 3> result = {x,y,z};
                return result;
            }
        }
    }

    // Function to generate kick velocity of a GW merger with spins
    /*! @param[out] vkick: kick velocity: [vx, vy, vz]
        @param[in] _Chi1: spin of the primary body
        @param[in] _Chi2: spin of the secondary body
        @param[in] _L: orbital angular momentum
        @param[in] _dr: separation vector
        @param[in] q: mass ratio
        @param[in] maxkick: whether to calculate the maximum kick (false)
        @param[in] inverteraxisl: whether to invert the axis of L (true)
    */
    void calcKickVel(Float vkick[], 
        const Float _Chi1[], 
        const Float _Chi2[], 
        const Float _L[],
        const Float _dr[],
        const Float& q, 
        bool maxkick = false,
        bool inverteraxisl = true ) {

        const std::array<Float, 3> L = {_L[0], _L[1], _L[2]};
        const std::array<Float,3> dr = {_dr[0], _dr[1], _dr[2]};
        std::array<Float, 3> Chi1{ _Chi1[0], _Chi1[1],  _Chi1[2] };
        std::array<Float, 3> Chi2{ _Chi2[0], _Chi2[1],  _Chi2[2] };
        std::array<Float, 3> hatL = { 0.0, 0.0, 1.0 };

        Float eta = q * pow(1.0 + q, -2);   // Symmetric mass ratio

        // Constants
        const Float A = 1.2e4; // km/s
        const Float B = -0.93;
        const Float H = 6.9e3; // km/s
        const Float V11 = 3677.76; // km/s
        const Float VA = 2481.21; // km/s
        const Float VB = 1792.45; // km/s
        const Float VC = 1506.52; // km/s
        const Float C2 = 1140.0; // km/s
        const Float C3 = 2481.0; // km/s
        const Float M_pi = 3.14159265358979323846;
        //std::array<Float, 3> hatS1 = normalize(Chi1);
        //std::array<Float, 3> hatS2 = normalize(Chi2);

        std::array<Float, 3> Delta = calcDelta(Chi1, Chi2, L, dr,q);
        Float Delta_par = dot(Delta, hatL);
        Float Delta_perp = norm(cross(Delta, hatL));
        std::array<Float, 3> chit = calcChit(Chi1, Chi2, L, dr,q);
        Float chit_par = dot(chit, hatL);
        Float chit_perp = norm(cross(chit, hatL));
        Float zeta = 145.0 * M_pi / 180.0; // convert degrees to radians

        Float bigTheta;
        if (maxkick) {
            bigTheta = 0.0;
        }
        else {
            bigTheta = rand_f64() * 2.0 * M_pi;
        }

        Float vm = A * eta * eta * (1.0 + B * eta) * (1.0 - q) / (1.0 + q);
        Float vperp = H * eta * eta * Delta_par;
        Float vpar = 16.0 * eta * eta * (Delta_perp * (V11 + 2.0 * VA * chit_par + 4.0 * VB * pow(chit_par, 2) + 8.0 * VC * pow(chit_par, 3)) + chit_perp * Delta_par * (2.0 * C2 + 4.0 * C3 * chit_par)) * std::cos(bigTheta);
        
        std::array<Float, 3> vkick_new;
        if (inverteraxisl)
         {
            std::array<Float, 3> vkick_L = { vm + vperp * std::cos(zeta), vperp * std::sin(zeta), vpar };
            //for (int i = 0; i < 3; i++)
            //   std::cout << vkick_L[i] << std::endl;
            vkick_new =  calcInverterAxisL(L, dr, vkick_L); 
        }
        else {
            vkick_new = { vm + vperp * std::cos(zeta), vperp * std::sin(zeta), vpar };
        }
        for (int i = 0; i < 3; i++) vkick[i] = vkick_new[i]/vscale;
    }

    // Function to generate final mass of a GW merger with spins
    /*! @param[out] Mfin: final mass
        @param[in] _Chi1: spin of the primary body
        @param[in] _Chi2: spin of the secondary body
        @param[in] _L: orbital angular momentum
        @param[in] _dr: separation vector
        @param[in] q: mass ratio
    */
    void calcFinalMass(Float& Mfin,
        const Float _Chi1[], 
        const Float _Chi2[], 
        const Float _L[],
        const Float _dr[],
        const Float& q){
        
        const std::array<Float, 3> L = {_L[0], _L[1], _L[2]};
        const std::array<Float,3> dr = {_dr[0], _dr[1], _dr[2]};
        std::array<Float, 3> Chi1{ _Chi1[0], _Chi1[1],  _Chi1[2] };
        std::array<Float, 3> Chi2{ _Chi2[0], _Chi2[1],  _Chi2[2] };

        std::array<Float, 3> hatL = { 0.0, 0.0, 1.0 };

        Float eta = q * pow(1.0 + q, -2);
        std::array<Float, 3> chit = calcChit(Chi1, Chi2, L, dr,q);
        Float chit_par = dot(chit, hatL);
        const Float p0 = 0.04827;
        const Float p1 = 0.01707;
        Float Z1 = 1 + pow((1-pow(chit_par, 2)),1/3)* (pow((1+chit_par),1/3)+pow((1-chit_par),1/3));
        Float Z2 = pow((3*pow(chit_par,2) + pow(Z1,2)),1/2);
        Float risco = 3 + Z2 - std::copysign(1.0,chit_par)*pow((3-Z1)*(3+Z1+2*Z2),1/2);
        Float Eisco = pow((1-2/(3*risco)),1/2);
        Float Erad = eta*(1-Eisco) + 4*pow(eta,2)*(4*p0+16*p1*chit_par*(chit_par+1)+Eisco -1);
        Mfin = 1 - Erad;
    }

    // Function to calculate the angles theta1, theta2, and theta12
    /*! @param[in] _Chi1: spin of the primary body
        @param[in] _Chi2: spin of the secondary body
        @param[in] _L: orbital angular momentum
        @param[in] _dr: separation vector
        @param[in] q: mass ratio
        \return the angles theta1, theta2, and theta12
    */
    std::tuple<Float, Float, Float> calcAngle( 
        const std::array<Float, 3>& Chi1, 
        const std::array<Float, 3>& Chi2,
        const std::array<Float, 3>& L,
        const std::array<Float, 3>& dr) {
        std::array<Float, 3> hatL = {0, 0, 1};
        std::array<Float, 3> Chi1_new = Chi1;
        std::array<Float, 3> Chi2_new = Chi2;
        std::tie(Chi1_new, Chi2_new) = calcAxisL(L, dr, Chi1, Chi2);

        Float dot_chi1L = dot(Chi1_new, hatL);
        Float dot_chi2L = dot(Chi2_new, hatL);
        Float dot_chi12 = dot(Chi1_new, Chi2_new);

        Float chi1_norm = norm(Chi1_new);
        Float chi2_norm = norm(Chi2_new);
        Float hatL_norm = norm(hatL);

        Float theta1 = std::acos(std::clamp(dot_chi1L / (chi1_norm * hatL_norm), -1.0, 1.0));
        Float theta2 = std::acos(std::clamp(dot_chi2L / (chi2_norm * hatL_norm), -1.0, 1.0));
        Float theta12 = std::acos(std::clamp(dot_chi12 / (chi2_norm * chi1_norm), -1.0, 1.0));

        return {theta1, theta2, theta12};
    }

    // Function to generate final spin of a GW merger with spins 
    /*! @param[out] spin: spin velocity
        @param[in] _Chi1: spin of the primary body
        @param[in] _Chi2: spin of the secondary body
        @param[in] _L: orbital angular momentum
        @param[in] _dr: separation vector
        @param[in] q: mass ratio
        @param[in] which: which model to use
        \return the kick spin of the GW merger
        */
    void calcFinalSpin(Float spin[],
        const Float _Chi1[], 
        const Float _Chi2[], 
        const Float _L[],
        const Float _dr[],
        const Float& q,
        const std::string& which = "HBR16_34corr") {
        
        const std::array<Float, 3> L = {_L[0], _L[1], _L[2]};
        const std::array<Float,3> dr = {_dr[0], _dr[1], _dr[2]};
        std::array<Float, 3> Chi1{ _Chi1[0], _Chi1[1],  _Chi1[2] };
        std::array<Float, 3> Chi2{ _Chi2[0], _Chi2[1],  _Chi2[2] };

        std::array<Float, 3> Chi1_new;
        std::array<Float, 3> Chi2_new;
        std::tie(Chi1_new, Chi2_new) = calcAxisL(L, dr, Chi1, Chi2);
        
        const std::array<Float, 3> hatL = {0, 0, 1};

        Float eta = q * std::pow(1.0 + q, -2);
        const Float M_pi = 3.14159265358979323846;

        Float norms_chi1 = norm(Chi1_new);
        Float norms_chi2 = norm(Chi2_new);
        std::vector<std::vector<Float>> kfit;
        Float xifit;

        if (which == "HBR16_12" || which == "HBR16_12corr") {
            kfit = {{NAN, -1.2019, -1.20764}, {3.79245, 1.18385, 4.90494}};
            xifit = 0.41616;
        } else if (which == "HBR16_33" || which == "HBR16_33corr") {
            kfit = {{NAN, 2.87025, -1.53315, -3.78893}, {32.9127, -62.9901, 10.0068, 56.1926}, {-136.832, 329.32, -13.2034, -252.27}, {210.075, -545.35, -3.97509, 368.405}};
            xifit = 0.463926;
        } else if (which == "HBR16_34" || which == "HBR16_34corr") {
            kfit = {{NAN, 3.39221, 4.48865, -5.77101, -13.0459}, {35.1278, -72.9336, -86.0036, 93.7371, 200.975}, {-146.822, 387.184, 447.009, -467.383, -884.339}, {223.911, -648.502, -697.177, 753.738, 1166.89}};
            xifit = 0.474046;
        } else {
            throw std::invalid_argument("`which` needs to be one of the following: `HBR16_12`, `HBR16_12corr`, `HBR16_33`, `HBR16_33corr`, `HBR16_34`, `HBR16_34corr`.");
        }

        // Calculate K00 from Eq 11
        size_t n = kfit.size();
        size_t m = kfit[0].size();
        Float sum = 0.0;
        for (size_t i = 1; i < n; ++i) {
            sum += kfit[i][0] / std::pow(4.0, 3 + (i - 1));
        }
        kfit[0][0] = std::pow(4.0, 2) * (0.68646 - sum - std::sqrt(3.0) / 2.0);

        auto [theta1, theta2, theta12] = calcAngle(Chi1, Chi1,L,dr);

        // Eq. 18
        if (which.find("corr") != std::string::npos) {
            Float eps1 = 0.024;
            Float eps2 = 0.024;
            Float eps12 = 0;
            theta1 = theta1 + eps1 * std::sin(theta1);
            theta2 = theta2 + eps2 * std::sin(theta2);
            theta12 = theta12 + eps12 * std::sin(theta12);
        }
        // Eq. 14 - 15
        Float atot = (norms_chi1 * std::cos(theta1) + norms_chi2 * std::cos(theta2) * q * q) / std::pow(1 + q, 2);
    
        Float aeff = atot + xifit * eta * (norms_chi1 * std::cos(theta1) + norms_chi2 * std::cos(theta2));
        
        // Eq. 2 - 6 evaluated at aeff, as specified in Eq. 11
        Float Z1 = 1 + std::pow(1 - aeff * aeff, 1.0 / 3.0) * (std::pow(1 + aeff, 1.0 / 3.0) + std::pow(1 - aeff, 1.0 / 3.0));
        Float Z2 = std::sqrt(3 * aeff * aeff + Z1 * Z1);
        Float risco = 3 + Z2 - std::copysign(1.0, aeff) * std::sqrt((3 - Z1) * (3 + Z1 + 2 * Z2));
        Float Eisco = std::sqrt(1 - 2 / (3 * risco));
        Float Lisco = (2 / (3 * std::sqrt(3))) * (1 + 2 * std::sqrt(3 * risco - 2));
        
        // Eq. 13

        // calculate etatoi
        std::vector<Float> etatoi(n);
        for (size_t i = 0; i < n; ++i) {
            etatoi[i] = std::pow(eta, 1 + i);
        }

        // calculate innersum
        std::vector<Float> innersum(m, 0.0);
        for (size_t j = 0; j < m; ++j) {
            for (size_t i = 0; i < n; ++i) {
                innersum[j] += kfit[i][j] * etatoi[i];
            }
        }

        // calculate aefftoj
        std::vector<Float> aefftoj(m);
        for (size_t j = 0; j < m; ++j) {
            aefftoj[j] = std::pow(aeff, j);
        }

        // calculate sumell
        Float sumell = 0.0;
        for (size_t j = 0; j < m; ++j) {
            sumell += innersum[j] * aefftoj[j];
        }

        Float ell = std::abs(Lisco - 2.0 * atot * (Eisco - 1.0) + sumell);

        // Eq. 16
        Float chifin = 1 / std::pow(1 + q, 2) * std::sqrt(
            (std::pow(norms_chi1, 2) + std::pow(norms_chi2, 2) * std::pow(q, 4) + 
            2 * norms_chi1 * norms_chi2 * std::pow(q, 2) * std::cos(theta12) + 
            2 * (norms_chi1 * std::cos(theta1) + norms_chi2 * std::pow(q, 2) * std::cos(theta2)) * ell * q + 
            std::pow(ell * q, 2)));
        chifin = std::clamp(chifin, 0.0, 1.0);

        // Calculate direction of spin
        Float theta_xy = rand_f64() * 2.0 * M_pi;
        std::array<Float, 3> J = {hatL[0] + Chi1_new[0] + Chi2_new[0], hatL[1] + Chi1_new[1] + Chi2_new[1], hatL[2] + Chi1_new[2] + Chi2_new[2]};
        Float theta_z = std::acos(dot(J, hatL) / (norm(J) * norm(hatL)));

        std::array<Float, 3> chikick = {chifin * std::sin(theta_z) * std::cos(theta_xy), chifin * std::sin(theta_z) * std::sin(theta_xy), chifin * std::cos(theta_z)};
        
        std::array<Float,3> chikick_new =  calcInverterAxisL(L, dr, chikick);
        for (int i = 0; i < 3; i++) spin[i] = chikick_new[i];
    }

};