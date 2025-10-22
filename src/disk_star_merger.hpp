#pragma once
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <string>
#include <getopt.h>
#include "io.hpp"
#include <cassert>
#include "Common/Float.h"
#ifndef ASSERT
#define ASSERT assert
#endif
#include "Common/binary_tree.h"
#include "../parallel-random/rand.hpp"

//! Class for managing mergers
enum class StarType:int {none = -1, smbh = 0, bh = 1, star = 2, seed = 3, star_remnant = 4, bh_remnant = 5};

//! IO parameters manager for external perturbation in hard integration
/*! For initializing the COMMON block variables from the commander option.
  The description of each parameter is also provided.
 */
class IOParamsDiskStarMerger{
public:
    IOParamsContainer input_par_store;
    IOParams<double> merger_mass_loss_rate; //!< mass loss rate after merger
    IOParams<double> stellar_radius_power_index; //!< stellar radius power index
    IOParams<double> stellar_radius_scale; //!< stellar radius scale
    IOParams<double> merger_time_delay; //!< time delay for merger to increase mass and change radius
    IOParams<double> stellar_seed_mass; //!< initial mass of star seed
    IOParams<double> initial_equilbrium_mass; //!< initial equilbrium mass of star;
    IOParams<double> lambda0; //!< fraction of star's intrinsic luminosity over the Eddington luminosity without merger
    IOParams<double> helium_fraction_disk; //!< helium fraction in the disk, used to calculate the equilbrium mass
    IOParams<double> salpeter_timescale; //!< salpeter timescale, for the star to reach equilbrium if NUMERIC_FLOAT_MAX, no growth
    IOParams<double> epsilon_helium; //!< helium enrichment efficiency
    IOParams<double> epsilon_bh; //!< the kenetic energy to radiation conversion efficiency of Eddington-limited accretion for BH
    IOParams<double> gravitational_constant; //!< gravitational constant
    IOParams<double> speed_of_light; //!< speed of light
    IOParams<long long int> redistribute_star_mode; //!< redistribute star mode, 0: no redistribute; 1: redistribute star position and velocity to opposite side of the center
    IOParams<std::string> fname_par;

    bool print_flag; //!< print flag
    //! Constructor
    IOParamsDiskStarMerger(): input_par_store(),
                              merger_mass_loss_rate(input_par_store, 0.0, "merger-mass-loss-rate", "mass loss rate for merger"),
                              stellar_radius_power_index(input_par_store, 0.6, "stellar-radius-power", "stellar radius power index 'n', rs = s M^n"),
                              stellar_radius_scale(input_par_store, 0.0046, "stellar-radius-scale", "stellar radius scale 's', rs = s M^n"),
                              merger_time_delay(input_par_store, 0.0, "merger-time-delay", "time delay for merger to increase mass"),
                              stellar_seed_mass(input_par_store, 10.0, "stellar-seed-mass", "initial mass of star seed"),  
                              initial_equilbrium_mass(input_par_store, 253.3124306069483, "initial-equilbrium-mass", "initial equilbrium mass of star"), 
                              lambda0(input_par_store, 0.75, "lambda0", "fraction of star's intrinsic luminosity over the Eddington luminosity without merger"),   
                              helium_fraction_disk(input_par_store, 0.28, "helium-fraction-disk", "helium fraction in the disk, used to calculate the equilbrium mass"),
                              salpeter_timescale(input_par_store, 0, "salpeter-timescale", "salpeter timescale, for the star to reach equilbrium, if zero, no stellar evolution"),
                              epsilon_helium(input_par_store, 0.006, "epsilon-helium", "helium enrichment efficiency, used to calculate helium enrichment timescale"),
                              epsilon_bh(input_par_store, 0.06, "epsilon-bh", "the kenetic energy to radiation conversion efficiency of Eddington-limited accretion for BH"),
                              gravitational_constant(input_par_store, 1.0, "G", "gravitational constant"),
                              speed_of_light(input_par_store, 1.0, "speed-of-light", "speed of light"),
                              redistribute_star_mode(input_par_store, 1, "redistribute-star-mode", "redistribute star mode, 0: no redistribute; 1: redistribute star position and velocity in random position along a circular orbit with the semi-major axis being the distance to the center; 2: redistribute star by choosing next type 3 star"),
                              fname_par    (input_par_store, "input.par", "p", "Input parameter file for external force (this option should be used first before any other options)",NULL,false),
                              print_flag(false) {}

    //! reading parameters from GNU option API
    /*!
      @param[in] argc: number of options
      @param[in] argv: string of options
      @param[in] opt_used_pre: already used option number from previous reading, use to correctly count the remaining argument number
      \return -1 if help is used; else the used number of argv
     */
    int read(int argc, char *argv[], const bool print_format_info=true, const int opt_used_pre=0) {
        static int merger_flag=-1;
        const struct option long_options[] = {
            {merger_mass_loss_rate.key, required_argument, &merger_flag, 0},  
            {stellar_radius_power_index.key, required_argument, &merger_flag, 1},  
            {stellar_radius_scale.key, required_argument, &merger_flag, 2},
            {merger_time_delay.key, required_argument, &merger_flag, 3},
            {stellar_seed_mass.key, required_argument, &merger_flag, 4},
            {initial_equilbrium_mass.key, required_argument, &merger_flag, 5},
            {lambda0.key, required_argument, &merger_flag, 6},
            {helium_fraction_disk.key, required_argument, &merger_flag, 7},
            {salpeter_timescale.key, required_argument, &merger_flag, 8},
            {epsilon_helium.key, required_argument, &merger_flag, 9},
            {epsilon_bh.key, required_argument, &merger_flag, 10},
            {speed_of_light.key, required_argument, &merger_flag, 11},
            {redistribute_star_mode.key, required_argument, &merger_flag, 12},
            {"help",      no_argument,       0, 'h'},
            {0,0,0,0}
        };

        int opt_used=opt_used_pre;
        int copt;
        int option_index;
        optind = 0;
        while ((copt = getopt_long(argc, argv, "-G:p:h", long_options, &option_index)) != -1) 
            switch (copt) {
            case 0:
                switch (merger_flag) {
                case 0:
                    merger_mass_loss_rate.value = atof(optarg);
                    if(print_flag) merger_mass_loss_rate.print(std::cout);
                    opt_used+=2;
                    break;            
                case 1:
                    stellar_radius_power_index.value = atof(optarg);
                    if(print_flag) stellar_radius_power_index.print(std::cout);
                    opt_used+=2;
                    break;            
                case 2:
                    stellar_radius_scale.value = atof(optarg);
                    if(print_flag) stellar_radius_scale.print(std::cout);
                    opt_used+=2;
                    break;
                case 3:
                    merger_time_delay.value = atof(optarg);
                    if(print_flag) merger_time_delay.print(std::cout);
                    opt_used+=2;
                    break;
                case 4:
                    stellar_seed_mass.value = atof(optarg);
                    if(print_flag) stellar_seed_mass.print(std::cout);
                    opt_used+=2;
                    break;
                case 5:
                    initial_equilbrium_mass.value = atof(optarg);
                    if(print_flag) initial_equilbrium_mass.print(std::cout);
                    opt_used+=2;
                    break;
                case 6:
                    lambda0.value = atof(optarg);
                    if(print_flag) lambda0.print(std::cout);
                    opt_used+=2;
                    break;
                case 7:
                    helium_fraction_disk.value = atof(optarg);
                    if(print_flag) helium_fraction_disk.print(std::cout);
                    opt_used+=2;
                    break;
                case 8:
                    salpeter_timescale.value = atof(optarg);
                    if(print_flag) salpeter_timescale.print(std::cout);
                    opt_used+=2;
                    break;
                case 9:
                    epsilon_helium.value = atof(optarg);
                    if(print_flag) epsilon_helium.print(std::cout);
                    opt_used+=2;
                    break;
                case 10:
                    epsilon_bh.value = atof(optarg);
                    if(print_flag) epsilon_bh.print(std::cout);
                    opt_used+=2;
                    break;
                case 11:
                    speed_of_light.value = atof(optarg);
                    if(print_flag) speed_of_light.print(std::cout);
                    opt_used+=2;
                    break;
                case 12:
                    redistribute_star_mode.value = atoi(optarg);
                    if(print_flag) redistribute_star_mode.print(std::cout);
                    opt_used+=2;
                    break;
                default:
                    break;
                }
                break;
            case 'G':
                gravitational_constant.value = atof(optarg);
                if(print_flag) gravitational_constant.print(std::cout);
                opt_used+=2;
                break;
            case 'p':
                fname_par.value = optarg;
                if(print_flag) {
                    std::string fgalpy_par = fname_par.value+".disk_star_merger"; 
                    FILE* fpar_in;
                    if( (fpar_in = fopen(fgalpy_par.c_str(),"r")) == NULL) {
                        fprintf(stderr,"Error: Cannot open file %s.\n", fgalpy_par.c_str());
                        abort();
                    }
                    input_par_store.readAscii(fpar_in);
                    fclose(fpar_in);
                }
                opt_used+=2;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
                input_par_store.mpi_broadcast();
                PS::Comm::barrier();
#endif
                break;
            case 'h':
                if(print_flag){
                    std::cout<<"Disk star merger parameters, options:"<<std::endl;
                    input_par_store.printHelp(std::cout, print_format_info);
                }
                return -1;
            case '?':
                opt_used +=2;
                break;
            default:
                break;
            }
        
        if(print_flag) std::cout<<"----- Finish reading input options of disk star merger parameters -----\n";

        return opt_used;
    }
};

class DiskStarMergerManager {
public:
    Float merger_mass_loss_rate; //!< mass loss rate
    Float stellar_radius_power_index; //!< radius amplifier
    Float stellar_radius_scale; //!< radius scale
    Float merger_time_delay; //!< time delay for merger to increase mass
    Float stellar_seed_mass; //!< initial mass of star seed
    Float initial_equlibrium_mass; //!< initial equilbrium mass of star;
    Float lambda0; //!< fraction of star's intrinsic luminosity over the Eddington luminosity
    Float helium_fraction_disk; //!< helium fraction in the disk, used to calculate the equilbrium mass
    Float salpeter_timescale; //!< Salpeter timescale, for the star to reach equilbrium if NUMERIC_FLOAT_MAX, no growth
    Float epsilon_helium; //!< helium enrichment efficiency, used to calculate helium enrichment timescale
    Float epsilon_bh; //!< the kenetic energy to radiation conversion efficiency of Eddington-limited accretion for BH
    Float gravitational_constant; //!< gravitational constant
    Float speed_of_light; //!< speed of light
    int redistribute_star_mode; //!< redistribute star mode, 0: no redistribute; 1: redistribute star position and velocity to opposite side of the center
    
    //! Constructor
    DiskStarMergerManager(): merger_mass_loss_rate(-1.0), 
                             stellar_radius_power_index(-1.0), 
                             stellar_radius_scale(0.0), 
                             merger_time_delay(-1.0), 
                             stellar_seed_mass(0.0), 
                             initial_equlibrium_mass(0.0), 
                             lambda0(0.0),
                             helium_fraction_disk(0.0),
                             salpeter_timescale(NUMERIC_FLOAT_MAX),
                             epsilon_helium(0.0),
                             epsilon_bh(0.0), 
                             gravitational_constant(0.0),
                             speed_of_light(0.0),
                             redistribute_star_mode(-1) {}

    //! (Necessary) check whether publicly initialized parameters are correctly set
    /*! \return true: all parmeters are correct. In this case no parameters, return true;
     */
    bool checkParams() {
        assert(merger_mass_loss_rate>=0.0);
        assert(stellar_radius_power_index>=0.0);
        assert(stellar_radius_scale>0.0);
        assert(merger_time_delay>=0.0);
        assert(stellar_seed_mass>0.0);
        assert(initial_equlibrium_mass>0.0);
        assert(lambda0>0.0);
        assert(helium_fraction_disk>0.0);
        assert(salpeter_timescale>=0.0);
        assert(epsilon_helium>0.0 && epsilon_helium<=1.0);
        assert(epsilon_bh>0.0 && epsilon_bh<=1.0);
        assert(gravitational_constant>0.0);
        assert(speed_of_light>0.0);
        assert(redistribute_star_mode>=0 && redistribute_star_mode<=2);
        return true;
    }

    //! print parameters
    void print(std::ostream & _fout) const{
        _fout<<"merger_mass_loss_rate : "<<merger_mass_loss_rate<<std::endl
             <<"stellar_radius_power_index : "<<stellar_radius_power_index<<std::endl
             <<"stellar_radius_scale : "<<stellar_radius_scale<<std::endl
             <<"merger_time_delay : "<<merger_time_delay<<std::endl
             <<"stellar_seed_mass : "<<stellar_seed_mass<<std::endl
             <<"initial_equlibrium_mass : "<<initial_equlibrium_mass<<std::endl
             <<"lambda0 : "<<lambda0<<std::endl
             <<"helium_fraction_disk : "<<helium_fraction_disk<<std::endl
             <<"salpeter_timescale : "<<salpeter_timescale<<std::endl
             <<"epsilon_helium : "<<epsilon_helium<<std::endl
             <<"epsilon_bh : "<<epsilon_bh<<std::endl
             <<"gravitational_constant : "<<gravitational_constant<<std::endl
             <<"speed_of_light : "<<speed_of_light<<std::endl
             <<"redistribute_star_mode : "<<redistribute_star_mode<<std::endl;
    }

    //! initial parameters for disk star mergers
    /*!
      @param[in] _input: input parameter
      @param[in] _print_flag: printing flag
     */
    void initial(const IOParamsDiskStarMerger& _input, const bool _print_flag=false) {
        merger_mass_loss_rate = _input.merger_mass_loss_rate.value;
        stellar_radius_power_index = _input.stellar_radius_power_index.value;
        stellar_radius_scale = _input.stellar_radius_scale.value;
        merger_time_delay = _input.merger_time_delay.value;
        stellar_seed_mass = _input.stellar_seed_mass.value;
        initial_equlibrium_mass = _input.initial_equilbrium_mass.value;
        lambda0 = _input.lambda0.value;
        helium_fraction_disk = _input.helium_fraction_disk.value;
        salpeter_timescale = _input.salpeter_timescale.value;
        epsilon_helium = _input.epsilon_helium.value;
        epsilon_bh = _input.epsilon_bh.value;
        gravitational_constant = _input.gravitational_constant.value;
        speed_of_light = _input.speed_of_light.value;
        redistribute_star_mode = _input.redistribute_star_mode.value;
    }

    //! calcMergerProperties
    /*! Calculate the new properties of the merger
        @param[in,out] p1: particle 1, will be merger 
        @param[in] p2: particle 2, will be zero mass particle
        @param[in] time: current time
        */
    template <class TParticle>
    void calcMergerProperties(TParticle* p1, TParticle* p2, const Float& time) {

        TParticle *pm, *p0; // set final merger to star or first BH
        if (p1->star.getType()==StarType::bh && p2->star.getType()==StarType::star) {
            pm = p2;
            p0 = p1;
        }
        else if(p1->star.getType()==StarType::seed) {
            pm = p2;
            p0 = p1;
        }
        else if(p1->star.getType()!=StarType::smbh && p2->star.getType()==StarType::smbh) {
            pm = p2;
            p0 = p1;
        }
        else {
            pm = p1;
            p0 = p2;
        }

        Float mcm = p1->mass + p2->mass;
        for (int k=0; k<3; k++) {
            pm->pos[k] = (p1->mass*p1->pos[k] + p2->mass*p2->pos[k])/mcm;
            pm->vel[k] = (p1->mass*p1->vel[k] + p2->mass*p2->vel[k])/mcm;
            p0->pos[k] = -pm->pos[k];
            p0->vel[k] = -pm->vel[k];
        }
        
        // only increase mass and change radius after time delay
        if (time > merger_time_delay) {
            Float new_mass = mcm * (1 - merger_mass_loss_rate);    
            if (pm->star.getType() == StarType::star) {
                pm->radius = stellar_radius_scale * std::pow(new_mass, stellar_radius_power_index);
                if (p0->star.getType() == StarType::star) {
                    pm->star.helium_fraction = (p1->star.helium_fraction*p1->mass + p2->star.helium_fraction*p2->mass)/mcm;
                }
            }
            else if (pm->star.getType() == StarType::bh) {
                pm->radius = gravitational_constant * new_mass / (speed_of_light * speed_of_light);
            }
            pm->mass = new_mass;
            pm->dm += new_mass - pm->mass;
        }

        p0->dm -= p0->mass;
        p0->mass = 0.0;
        p0->radius = 0.0;

        if (p0->star.getType() == StarType::bh) {
            pm->star.n_merger_bh++;
            p0->star.setType(StarType::bh_remnant);
        }
        else if (p0->star.getType() == StarType::star) {
            pm->star.n_merger_star++;
            pm->star.n_merger_bh += p0->star.n_merger_bh;
            p0->star.setType(StarType::star_remnant);
        }

        pm->star.last_merger_time = time;
        pm->star.last_mass_change_time = time;
        p0->star.last_mass_change_time = time;
        p0->star.n_merger_star = 0;
        p0->star.n_merger_bh = 0;
    }

    //! calculate mass change
    /*! Calculate mass change, for mass < target mass, increase mass; for mass > equilbrium mass, decrease mass
        Increase mass formula:  
        dM/dt = c M^2
        c: mass growth factor
        @param[in] p: particle
        @param[in] time: current time

        \return 0: no change; 1: modified mass
        */
    template <class TParticle>
    int calcMassChange(TParticle* p, const Float& time) {
        int return_flag = 0;

        // no mass change if salpeter_timescale is 0        
        if (salpeter_timescale == 0) return 0; 

        // if type is star, evolve mass to equilbrium mass    
        if (p->star.getType()==StarType::star) {
            Float dt = time - p->star.last_mass_change_time;
            if (dt>0) {
                // Helium growth
                p->star.helium_fraction += lambda0 / (epsilon_helium * salpeter_timescale) * dt;
                
                // Helium fraction should be between 0 and 1, if 1, evolve to BH
                if (p->star.helium_fraction > 1.0) {

                    // calculate time to reach BH                       
                    Float time_bh_form = epsilon_helium * salpeter_timescale * (1.0 - p->star.helium_fraction) / lambda0 + time;
                    p->star.helium_fraction = 1.0;
                    
                    // if helium fraction is 1, then star evolve to post-main sequence and eventually become a BH
                    Float new_mass = initial_equlibrium_mass * std::pow(helium_fraction_disk, 2.5);
                    p->dm += new_mass - p->mass;
                    p->mass = new_mass;

                    // set type to BH
                    p->star.setType(StarType::bh);
                    p->star.last_mass_change_time = time_bh_form;
                    // set Swartzchild radius
                    p->radius = gravitational_constant * p->mass / (speed_of_light * speed_of_light);

                    return_flag = 1;
                }
                else {
                    // calculate equilbrium mass
                    Float equilbrium_mass = initial_equlibrium_mass * std::pow(helium_fraction_disk / p->star.helium_fraction, 2.5);

                    // Eddington accretion rate m'_edd = r/r_grav * m / tau_salpeter = r c^2 / (G tau_salpeter)
                    Float mdot_eddington = p->radius * speed_of_light * speed_of_light / (gravitational_constant * salpeter_timescale);

                    // factor of star's intrinsic luminosity over the Eddington luminosity
                    Float m_frac_8 = std::pow(p->mass / equilbrium_mass, 8.0);
                    Float s_fb = ( 1 - m_frac_8/(1 + m_frac_8));
                    s_fb = s_fb * s_fb;

                    // mass accretion rate m'_acc = S_feedback * m'_edd
                    // wind mass loss rate m'_wind = lambda*(1 - S_feedback)/2 * m'_edd
                    // net mass change rate m'_net = m'_acc - m'_wind
                    Float new_mass = p->mass + (s_fb - lambda0 * (1 - s_fb)/2) * mdot_eddington * dt;

                    // update parameters
                    p->dm += new_mass - p->mass;
                    p->mass = new_mass;

                    p->radius = stellar_radius_scale * std::pow(new_mass, stellar_radius_power_index);
                    p->star.last_mass_change_time = time;

                    return_flag = 1;
                }
            }         
        }

        // if type is black hole, considering accretion
        if (p->star.getType() == StarType::bh) {
            Float dt = time - p->star.last_mass_change_time;
            if (dt>0) {
                // Eddington accretion 
                Float mdot = p->mass /(epsilon_bh * salpeter_timescale);
                Float new_mass = p->mass + mdot * dt;

                p->dm += new_mass - p->mass;
                p->mass = new_mass;

                p->radius = gravitational_constant * p->mass / (speed_of_light * speed_of_light);
                p->star.last_mass_change_time = time;

                return_flag = 1;
            }
        }
        return return_flag;
    }

    //! Set the remnant orbit to the center of mass
    /*! 
        @param[in,out] p: particle array
        @param[in] p_cm: center of mass particle
    */
    template <class TParticle, class Tcm>
    void setRemnantOrbitToCM(TParticle& p, const Tcm& p_cm) {
        if (p.star.getType() == StarType::star_remnant) {
            for (int k=0; k<3; k++) {
                p.pos[k] = p_cm.pos[k];
                p.vel[k] = p_cm.vel[k];
            }
        }
    }

    //! redistribute star position and velocity
    /*!
        redistribute star position and velocity to opposite side of the center     
        @param[in,out] p: particle to redistribute
        @param[in] center: center particle
        @param[in] system: particle system for search star seed
        @param[in] n_system: number of particles in system
        @return 0: no redistribute or need to delete star; 1: redistribute star
    */
    template <class TParticle>
    int redistributeStar(TParticle* p, TParticle* center, TParticle* system, const int n_system) {
        // if type is star, redistribute 
        if (p->star.getType() == StarType::star_remnant) {
            if (redistribute_star_mode == 1) {
                COMM::Binary bin;
                bin.calcOrbit(*p,*center, gravitational_constant);
                if (bin.r>0) {
                    // redistribute star position and velocity to random position assuming a circular orbit with semi-major axis r
                    bin.semi = bin.r;
                    bin.ecc = 0.0;
                    bin.ecca = 2 * COMM::PI * rand_f64();
                    TParticle cm;
                    bin.calcParticles(*p, cm, gravitational_constant);
                    
                    for (int k=0; k<3; k++) {
                        p->pos[k] += center->pos[k] - cm.pos[k];
                        p->vel[k] += center->vel[k] - cm.vel[k];
                    }
                    p->star.setType(StarType::star);
                    p->mass = stellar_seed_mass;
                    p->dm = 0.0;
                    p->radius = stellar_radius_scale * std::pow(stellar_seed_mass, stellar_radius_power_index);
                    return 1;
                }
            }
            else if (redistribute_star_mode == 2) {
                // redistribute star by choosing a star seed with a similar distance to center
                Float pos[3] = {p->pos[0] - center->pos[0], 
                                p->pos[1] - center->pos[1], 
                                p->pos[2] - center->pos[2]};
                Float r = std::sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
                // Search next 10 star seeds and pickup one with the closest r
                int index_seed = -1;
                Float dr = NUMERIC_FLOAT_MAX;
                int n_found = 0;
                for (int i=0; i<n_system; i++) {
                    auto pi = &system[i];
                    if (pi->star.getType() == StarType::seed) {
                        Float pos_seed[3] = {pi->pos[0] - center->pos[0], 
                                             pi->pos[1] - center->pos[1], 
                                             pi->pos[2] - center->pos[2]};
                        Float ri = std::sqrt(pos_seed[0]*pos_seed[0] + pos_seed[1]*pos_seed[1] + pos_seed[2]*pos_seed[2]);
                        Float dri = std::abs(r - ri);
                        if (dri<dr) {
                            dr = dri;
                            index_seed = i;
                        }
                        n_found++;
                        if (n_found>=10) break;
                    }
                }
                if (index_seed >= 0) {
                    auto pi = &system[index_seed];
                    pi->star.setType(StarType::star);
                    pi->mass = stellar_seed_mass;
                    pi->dm = 0.0;
                    pi->radius = stellar_radius_scale * std::pow(stellar_seed_mass, stellar_radius_power_index);
                    pi->star.last_mass_change_time = p->star.last_mass_change_time;
                        
                    pi->star.n_merger_star = 0;
                    pi->star.n_merger_bh = 0;
                    return 0;
                }
            }
        }
        return 0;
    } 

};


//! class for disk star merger parameters of individual stars
class StarParameter{
public:
    long long int type; //!< type of object; 0: supermassive black hole; 1: black hole; 2: star; 3: star seed; 4: star zero mass remnant; 5: black hole zero mass remnant
    long long int n_merger_star; //!< times of merger with star
    long long int n_merger_bh; //!< times of merger with black hole
    Float last_mass_change_time; //!< last mass change time
    Float last_merger_time; //!< last merger time
    Float helium_fraction; //!< helium fraction

    //! Constructor
    StarParameter():  type(-1),
                      n_merger_star(0),
                      n_merger_bh(0),
                      last_mass_change_time(0.0),
                      last_merger_time(0.0),
                      helium_fraction(0.0)
                      {}

    //! set type
    /*!
      @param[in] _type: type of object; 0: supermassive black hole; 1: black hole; 2: star; 3: star seed; 4: star zero mass remnant; 5: black hole zero mass remnant
     */
    void setType(const StarType _type) {
        type = static_cast<long long int>(_type);
    }

    //! get type
    /*!
      @return type of object; 0: supermassive black hole; 1: black hole; 2: star; 3: star seed; 4: star zero mass remnant; 5: black hole zero mass remnant
     */
    StarType getType() const {
        return static_cast<StarType>(type);
    }

    //! initial parameters for disk star merger
    /*!
      @param[in] _type: type of object; 0: supermassive black hole; 1: black hole; 2: star; 3: star seed; 4: star zero mass remnant; 5: black hole zero mass remnant
      @param[in] _last_mass_change_time: time delay for mass approach target
     */
    void initial(const StarType _type, const Float& _last_mass_change_time=0.0) {
        type = static_cast<long long int>(_type);
        n_merger_star = 0;
        n_merger_bh = 0;
        last_mass_change_time = _last_mass_change_time;
        last_merger_time = 0;
        helium_fraction = 0.0;
    }

    //! write class data with ASCII format
    void writeAscii(FILE* fp) const{
        fprintf(fp, "%lld %lld %lld %26.17e %26.17e %26.17e\n", type, n_merger_star, n_merger_bh, last_mass_change_time, last_merger_time, helium_fraction);
    }

    //! read class data with ASCII format
    void readAscii(FILE* fp) {
        int rcount=fscanf(fp, "%lld %lld %lld %lf %lf %lf ", 
                          &type, &n_merger_star, &n_merger_bh, &last_mass_change_time, &last_merger_time, &helium_fraction);
        if(rcount<6) {
            std::cerr<<"Error: Data reading fails! requiring data number is 6, only obtain "<<rcount<<".\n";
            abort();
        }
    }

    //! for print in one line
    void print(std::ostream & fout) const{
        fout<<" type= "<<type
            <<" n_merger_star= "<<n_merger_star
            <<" n_merger_bh= "<<n_merger_bh
            <<" last_mass_change_time= "<<last_mass_change_time
            <<" last_merger_time= "<<last_merger_time
            <<" helium_fraction= "<<helium_fraction;
    }

    //! print titles of class members using column style
    /*! print titles of class members in one line for column style
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
     */
    static void printColumnTitle(std::ostream & _fout, const int _width=20) {
        _fout<<std::setw(_width)<<"type"
             <<std::setw(_width)<<"n_merger_star"
             <<std::setw(_width)<<"n_merger_bh"
             <<std::setw(_width)<<"last_mass_change_time"
             <<std::setw(_width)<<"last_merger_time"
             <<std::setw(_width)<<"helium_fraction";
    }

    //! print data of class members using column style
    /*! print data of class members in one line for column style. Notice no newline is printed at the end
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
     */
    void printColumn(std::ostream & _fout, const int _width=20) const{
        _fout<<std::setw(_width)<<type
             <<std::setw(_width)<<n_merger_star
             <<std::setw(_width)<<n_merger_bh
             <<std::setw(_width)<<last_mass_change_time
             <<std::setw(_width)<<last_merger_time
             <<std::setw(_width)<<helium_fraction;
    }

    //! print column title with meaning (each line for one column)
    /*! @param[out] _fout: std::ostream output object
      @param[in] _counter: offset of the number counter for each line to indicate the column index (defaulted 0)
      @param[in] _offset: the printing whitespace offset for each line (defaulted 0)
      \return: the total counter of columns
     */
    static int printTitleWithMeaning(std::ostream & _fout, const int _counter=0, const int _offset=0) {
        int counter = _counter;
        _fout<<std::setw(_offset)<<" "<<++counter<<". type: type of object; 0: SMBH; 1: stellar-mass BH; 2: star; 3: star seed; 4: star remnant; 5: BH remnant\n";
        _fout<<std::setw(_offset)<<" "<<++counter<<". n_merger_star: times of merger with star\n";
        _fout<<std::setw(_offset)<<" "<<++counter<<". n_merger_bh: times of merger with black hole\n";
        _fout<<std::setw(_offset)<<" "<<++counter<<". last_mass_change_time: last mass change time\n";
        _fout<<std::setw(_offset)<<" "<<++counter<<". last_merger_time: last merger time\n";
        _fout<<std::setw(_offset)<<" "<<++counter<<". helium_fraction: helium fraction\n";
        return counter;
    }
    
};
