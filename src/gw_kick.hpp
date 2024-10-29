#pragma once

#include <array>
#include <tuple>
#include "../parallel-random/rand.hpp"

#define Float double

//! Gravitational wave recoil kick calculater
class GWKick{
public:    
    Float speed_of_light;  ///> speed of light (km/s)
    Float vscale;          ///> velocity scale (km/s to pc/Myr)

    GWKick(): speed_of_light(-1), vscale(-1) {}

    //! check whether parameters values are correct
    /*! \return true: all correct
     */
    bool checkParams() {
        ASSERT(speed_of_light>0.0);
        ASSERT(vscale>0.0);
        return true;
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
        //const Float& c = speed_of_light;

        auto norm = [&](const std::array<Float, 3>& vec) {
            return std::sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
            };
        auto normalize = [&](const std::array<Float,3>& vec){
            return std::array<Float,3>{vec[0]/norm(vec),vec[1]/norm(vec),vec[2]/norm(vec)};
            }; 
        auto dot = [&](const std::array<Float, 3>& a, const std::array<Float, 3>& b) {
            return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
            };
        auto mod = [&](const std::array<Float, 3>& a) {
            return std::sqrt(dot(a, a));
            };
        auto dot_p = [&](const std::array<Float, 3>&a, const std::array<Float, 3>&b) {
            return dot(a, b) / mod(b);
            };
        auto cross = [&](const std::array<Float, 3>& a, const std::array<Float, 3>& b) {
            return std::array<Float,3>{ a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0] };
            };
        auto multiply = [&](const std::array<Float, 3>& vec, Float scalar) {
            return std::array<Float, 3>{ vec[0] * scalar, vec[1] * scalar, vec[2] * scalar };
            };
        auto subtract = [&](const std::array<Float, 3>& a, const std::array<Float, 3>& b) {
            return std::array<Float, 3>{ a[0] - b[0], a[1] - b[1], a[2] - b[2] };
            };
        auto add = [&](const std::array<Float, 3>& a, const std::array<Float, 3>& b) {
            return std::array<Float, 3>{ a[0] + b[0], a[1] + b[1], a[2] + b[2] };
            };
        auto divide = [&](const std::array<Float, 3>& vec, Float scalar) {
            return std::array<Float, 3>{ vec[0] / scalar, vec[1] / scalar, vec[2] / scalar };
            };
        auto matmul = [&](const std::array<std::array<Float, 3>, 3>& mat, const std::array<Float, 3>& vec) {
            std::array<Float, 3> result = { 0.0, 0.0, 0.0 };
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    result[i] += mat[i][j] * vec[j];
                }
            }
            return result;
        };

        // Coordinate transformation
        auto calcAxisL = [&](
            const std::array<Float, 3>& L,
            const std::array<Float, 3>& dr,
            const std::array<Float, 3>& Chi1,
            const std::array<Float, 3>& Chi2) {

                Float drm = mod(dr);
                Float Lm = mod(L);

                // Define unit vectors
                std::array<Float, 3> e1 = divide(dr, drm);
                std::array<Float, 3> e3 = divide(L, Lm);
                std::array<Float, 3> e2 = cross(e3, e1);

                // Transformations
                std::array<Float, 3> L_new = { dot_p(L, e1), dot_p(L, e2), dot_p(L, e3) };
                L_new = divide(L_new, Lm);

                std::array<Float, 3> Chi1_new = { dot_p(Chi1, e1), dot_p(Chi1, e2), dot_p(Chi1, e3) };
                std::array<Float, 3> Chi2_new = { dot_p(Chi2, e1), dot_p(Chi2, e2), dot_p(Chi2, e3) };

                return std::make_tuple(Chi1_new, Chi2_new);
            };

        auto calcInverterAxisL = [&](
            const std::array<Float, 3>& L,
            const std::array<Float, 3>& dr,
            const std::array<Float, 3>& vkick) {

                Float drm = mod(dr);
                Float Lm = mod(L);

                // Define unit vectors
                std::array<Float, 3> e1 = divide(dr, drm);
                std::array<Float, 3> e3 = divide(L, Lm);
                std::array<Float, 3> e2 = cross(e3, e1);

                // Construct matrix E
                std::array<std::array<Float, 3>, 3> E = { e1, e2, e3 };

                // Multiply matrix E by vector vkick
                std::array<Float, 3> vkick_new = matmul(E, vkick);
    
                return vkick_new;
            };


        std::array<Float, 3> Chi1{ _Chi1[0], _Chi1[1],  _Chi1[2] };
        std::array<Float, 3> Chi2{ _Chi2[0], _Chi2[1],  _Chi2[2] };

        std::tie(Chi1, Chi2) = calcAxisL(L, dr, Chi1, Chi2);

        // Normalization m1, m2, and eta;
        Float m1 = 1.0 / (1.0 + q);   // Primary mass
        Float m2 = q / (1.0 + q);     // Secondary mass
        Float eta = q * pow(1.0 + q, -2);   // Symmetric mass ratio

        // Turn vectors into scalars
        Float chi1 = norm(Chi1);
        Float chi2 = norm(Chi2);

        // Generate s1 and s2 (scalars)
        Float s1 = chi1 * m1 * m1; // Primary spin magnitude
        Float s2 = chi2 * m2 * m2; // Secondary spin magnitude

        // Spins here are defined in a frame with L along z and S1 in xz
        std::array<Float, 3> hatL = { 0.0, 0.0, 1.0 };

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
        std::array<Float, 3> hatS1 = normalize(Chi1);
        std::array<Float, 3> hatS2 = normalize(Chi2);

        std::array<Float, 3> Delta = subtract(multiply(hatS2, q * chi2), multiply(hatS1, chi1));
        Delta = multiply(Delta, -1.0 / (1.0 + q));
        Float Delta_par = dot(Delta, hatL);
        Float Delta_perp = norm(cross(Delta, hatL));
        std::array<Float, 3> chit = add(multiply(hatS2, q * q * chi2), multiply(hatS1, chi1));
        chit = multiply(chit, 1.0 / pow(1.0 + q, 2));
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
};