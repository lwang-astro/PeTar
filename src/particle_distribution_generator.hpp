#pragma once
#include "Common/binary_tree.h"

//! Particle distribution generator
class ParticleDistributionGenerator {
    // f(x) = x/sigma^2*exp(-x^2/(2sigma^2))
    // F(x) = 1.0 - exp(-x^2/(2sigma^2))
    // x = sqrt(-2*sigma^2*ln(1-F))
    // 2sigma^2 = <e^2> or <i^2>
    // <> means R.M.S.
    static double RayleighDistribution(const double sigma){
        static PS::MTTS mt;
        static bool first = true;
        if(first){
            mt.init_genrand( PS::Comm::getRank() );
            first = false;
        }
        double F = mt.genrand_res53();
        /*
          double ret = 0.0;
          do{
          ret = sqrt( -2.0*sigma*sigma*log(1.0-F));
          }while(ret >= 1.0);
          return ret;    
        */
        return sqrt( -2.0*sigma*sigma*log(1.0-F));
    }

    // p(a)da = C*a^-0.5*da
    // P(a) = 2*C*(a^0.5 - a_in^0.5)
    // a = (P/(2C) + a_in^0.5)^2
    static double HayashiDistribution(const double a_in, const double a_out){
        static PS::MTTS mt;
        static bool first = true;
        if(first){
            mt.init_genrand( PS::Comm::getRank() );
            first = false;
        }
        const double C = 0.5 / (sqrt(a_out) - sqrt(a_in));
        double P = mt.genrand_res53();
        double ret = P/(2*C) + sqrt(a_in);
        ret *= ret;
        return ret;
    }


    // p(a)da = C*a^(p+1)*da
    // P(a) = 2*C*(a^0.5 - a_in^0.5)
    static double HayashiDistributionWithIceLine(const double a_in, const double a_out,
                                          const double a_ice, const double f_ice=1.0, const double p_sigma=-1.5){
        //std::cout<<"a_in="<<a_in<<std::endl;
        //const PS::F64 p_sigma = -1.5;
        const PS::F64 p_mass = 2.0 + p_sigma;
        //std::cout<<"p_mass="<<p_mass<<std::endl;
        static PS::MTTS mt;
        //std::cout<<"check 0"<<std::endl;
        static bool first = true;
        //std::cout<<"check 1"<<std::endl;
        //std::cout<<"PS::Comm::getRank()="<<PS::Comm::getRank()<<std::endl;
        if(first){
            mt.init_genrand( PS::Comm::getRank() );
            first = false;
        }
        //std::cout<<"first="<<first<<std::endl;
        double P = mt.genrand_res53();
        double ret = 0.0;
        double C = 0.0;
        double P_ice = 0.0;
        if(a_ice <= a_in || a_ice >= a_out){
            //std::cout<<"check a"<<std::endl;
            C = p_mass / (pow(a_out, p_mass) - pow(a_in, p_mass));
            ret = pow(P*p_mass/C + pow(a_in, p_mass), 1.0/p_mass);
        }
        else{
            //std::cout<<"check b"<<std::endl;
            C = p_mass / ( f_ice*(pow(a_out, p_mass)-pow(a_ice, p_mass)) + pow(a_ice, p_mass) - pow(a_in, p_mass) );
            P_ice = C/p_mass*(pow(a_ice, p_mass) - pow(a_in, p_mass));
            //std::cout<<"C="<<C<<std::endl;
            //std::cout<<"P_ice="<<P_ice<<std::endl;
            if(P < P_ice){
                ret = pow(P*p_mass/C + pow(a_in, p_mass), 1.0/p_mass);
            }
            else{
                ret = pow( ((P*p_mass/C + pow(a_in, p_mass) - pow(a_ice, p_mass))/f_ice + pow(a_ice, p_mass)), 1.0/p_mass); // beyond the ice line
            }
        }
        return ret;
    }


    static double AU2CM(const double x){
        return x * 1.49597871e13;
    }

    static double CM2AU(const double x){
        return x / 1.49597871e13;
    }

public:
    // [M] = 1 [g]
    // [L] = 1 [cm]
    // [T] = 3.871e3[sec]
    // [V] = 2.5833118e-4 [cm/sec] = 2.5833118e-9 [km/sec]
    // 1yr = 8.146732e3[T]
    // or
    // G = 1 [L]^3/([M][T]^2)
    // [L] = 1[AU] = 1.495978707e8[km]
    // [M] = 1[Msun] = 1.989e30[kg]
    // [T] = 5.02198050479e6 [sec] = 0.15924595715 [yr]
    // [V] = 29.7886203575 [km/sec]
    static void makeKeplerDisk(PS::F64 & mass_planet_glb,
                               PS::F64 *& mass,
                               PS::F64vec *& pos,
                               PS::F64vec *& vel,
                               const long long int n_glb,
                               const long long int n_loc,
                               const double a_in, // [AU]
                               const double a_out, // [AU]
                               const double e_rms, // normalized
                               const double i_rms, // normalized
                               const double dens = 10.0, // [g/cm^2]
                               const double mass_sun = 1.0, //[m_sun]
                               const double a_ice = 0.0,
                               const double f_ice = 1.0,
                               const double power = -1.5,
                               const int seed = 0
        ){
        static const double mass_sun_gram = 1.989e33; //[g]
        PS::MTTS mt;
        //mt.init_genrand( PS::Comm::getRank() );
        mt.init_genrand( PS::Comm::getRank()+seed*PS::Comm::getNumberOfProc() );
        static const double PI = atan(1.0) * 4.0;
        const double AU = AU2CM(1.0);
        //mass_planet_glb = 4.0 * PI * dens * ( sqrt(a_out) - sqrt(a_in) ) * AU * AU / mass_sun_gram; // [Msun]
        mass_planet_glb = 2.0 * PI * dens / (2.0+power) * ( pow(a_out, 2.0+power) - pow(a_in, 2.0+power) ) * AU * AU / mass_sun_gram; // [Msun]
        const double m_planet = mass_planet_glb / n_glb;
        mass = new double[n_loc];
        pos = new PS::F64vec[n_loc];
        vel = new PS::F64vec[n_loc];
        const double h = pow(2.0*m_planet / (3.0*mass_sun), 1.0/3.0);
        PS::F64 e_ave = 0.0;
        const PS::F64 e_sigma = sqrt(0.5*e_rms*e_rms); // this is right procedure
        const PS::F64 i_sigma = sqrt(0.5*i_rms*i_rms);
        //const PS::F64 e_sigma = e_rms;
        //const PS::F64 i_sigma = i_rms;
        for(long long int i=0; i<n_loc; i++){
            mass[i] = m_planet;
            COMM::Binary bin;
            bin.semi = HayashiDistributionWithIceLine(a_in, a_out, a_ice, f_ice, power);
            bin.ecc  = RayleighDistribution(e_sigma) * h;
            bin.incline = RayleighDistribution(i_sigma) * h;
            bin.rot_horizon = 2.0 * PI * mt.genrand_res53();
            bin.rot_self    = 2.0 * PI * mt.genrand_res53();
            double l = 2.0 * PI * mt.genrand_res53(); // mean anomayl
            bin.ecca = bin.calcEccAnomaly(l, bin.ecc); // eccentric anomayl
            bin.m1 = mass_sun;
            bin.m2 = m_planet;
            ParticleBase sun,planet;
            bin.calcParticles(sun, planet);
            pos[i]=planet.pos;
            vel[i]=planet.vel;
            e_ave += bin.ecc*bin.ecc;
        }
        PS::F64 e_ave_glb = sqrt(PS::Comm::getSum(e_ave)/n_glb);
        if(PS::Comm::getRank() == 0){
            std::cerr<<"e_ave_glb="<<e_ave_glb<<std::endl;
            std::cerr<<"e_ave_hill="<<e_ave_glb/h<<std::endl;
            std::cerr<<"e_rms="<<e_rms<<std::endl;
        }
    }

    static void makePlummerModel(const double mass_glb,
                                 const long long int n_glb,
                                 const long long int n_loc,
                                 double *& mass,
                                 PS::F64vec *& pos,
                                 PS::F64vec *& vel,
                                 const double eng = -0.25,
                                 const int seed = 0){

        assert(eng < 0.0);
        static const double PI = atan(1.0) * 4.0;
        const double r_cutoff = 22.8 / (-3.0 * PI * mass_glb * mass_glb / (64.0 * -0.25)); // 22.8 is cutoff in Nbody units
        //const double r_cutoff = 22.8 * 0.25;
        mass = new double[n_loc];
        pos = new PS::F64vec[n_loc];
        vel = new PS::F64vec[n_loc];

        PS::MTTS mt;
        mt.init_genrand( PS::Comm::getRank() );
        for(int i=0; i<n_loc; i++){
            mass[i] = mass_glb / n_glb;
            double r_tmp = 9999.9;
            while(r_tmp > r_cutoff){ 
                double m_tmp = mt.genrand_res53();
                r_tmp = 1.0 / sqrt( pow(m_tmp, (-2.0/3.0)) - 1.0);
            }
            double phi = 2.0 * PI * mt.genrand_res53();
            double cth = 2.0 * (mt.genrand_real2() - 0.5);
            double sth = sqrt(1.0 - cth*cth);
            pos[i][0] = r_tmp * sth * cos(phi);
            pos[i][1] = r_tmp * sth * sin(phi);
            pos[i][2] = r_tmp * cth;
            while(1){
                const double v_max = 0.1;
                const double v_try = mt.genrand_res53();
                const double v_crit = v_max * mt.genrand_res53();
                if(v_crit < v_try * v_try * pow( (1.0 - v_try * v_try), 3.5) ){
                    const double ve = sqrt(2.0) * pow( (r_tmp*r_tmp + 1.0), -0.25 );
                    phi = 2.0 * PI * mt.genrand_res53();
                    cth = 2.0 * (mt.genrand_res53() - 0.5);
                    sth = sqrt(1.0 - cth*cth);
                    vel[i][0] = ve * v_try * sth * cos(phi);
                    vel[i][1] = ve * v_try * sth * sin(phi);
                    vel[i][2] = ve * v_try * cth;
                    break;
                }
            }
        }

        PS::F64vec cm_pos = 0.0;
        PS::F64vec cm_vel = 0.0;
        double  cm_mass = 0.0;
        for(PS::S32 i=0; i<n_loc; i++){
            cm_pos += mass[i] * pos[i];
            cm_vel += mass[i] * vel[i];
            cm_mass += mass[i];
        }
        cm_pos /= cm_mass;
        cm_vel /= cm_mass;
        for(PS::S32 i=0; i<n_loc; i++){
            pos[i] -= cm_pos;
            vel[i] -= cm_vel;
        }

        const double r_scale = -3.0 * PI * mass_glb * mass_glb / (64.0 * eng);
        const double coef = 1.0 / sqrt(r_scale);
        for(PS::S32 i=0; i<n_loc; i++){
            pos[i] *= r_scale;
            vel[i] *= coef;
        }

        double r_max_sq = -1.0;
        for(int i=0; i<n_loc; i++){
            if(r_max_sq < pos[i] * pos[i]){
                r_max_sq = pos[i] * pos[i];
            }
        }
    }

};
