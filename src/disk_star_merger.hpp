#pragma once
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <string>
#include <getopt.h>
#include "io.hpp"
#include "Common/Float.h"

//! IO parameters manager for external perturbation in hard integration
/*! For initializing the COMMON block variables from the commander option.
  The description of each parameter is also provided.
 */
class IOParamsDiskStarMerger{
public:
    IOParamsContainer input_par_store;
    IOParams<double> merger_mass_loss_rate; //!< mass loss rate after merger
    IOParams<double> merger_radius_amplifier_rate; //!< radius amplifier rate
    IOParams<double> merger_time_delay; //!< time delay for merger to increase mass and change radius
    IOParams<double> equalibrium_mass; //!< equilibrium mass
    IOParams<double> initial_mass; //!< initial mass
    IOParams<double> growth_time; //!< growth time to equalibrium mass    

    bool print_flag; //!< print flag
    //! Constructor
    IOParamsDiskStarMerger(): input_par_store(),
                              merger_mass_loss_rate(input_par_store, 0.0, "merger-mass-loss-rate", "mass loss rate for merger"),
                              merger_radius_amplifier_rate(input_par_store, 1.0, "merger-radius-rate", "particle radius amplifier rate after merger"),
                              merger_time_delay(input_par_store, 0.0, "merger-time-delay", "time delay for merger to increase mass"),
                              equalibrium_mass(input_par_store, 0.0, "equalibrium-mass", "mass for star approaching equilibrium"),
                              initial_mass(input_par_store, 0.0, "initial-mass", "initial mass for star growth"),
                              growth_time(input_par_store, 0.0, "growth-time", "time scale for mass growth"),
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
            {merger_radius_amplifier_rate.key, required_argument, &merger_flag, 1},  
            {merger_time_delay.key, required_argument, &merger_flag, 2},  
            {equalibrium_mass.key, required_argument, &merger_flag, 3},
            {initial_mass.key, required_argument, &merger_flag, 4},
            {growth_time.key, required_argument, &merger_flag, 5},
            {"help",      no_argument,       0, 'h'},
            {0,0,0,0}
        };

        int opt_used=opt_used_pre;
        int copt;
        int option_index;
        optind = 0;
        while ((copt = getopt_long(argc, argv, "-z:p:h", long_options, &option_index)) != -1) 
            switch (copt) {
            case 0:
                switch (merger_flag) {
                case 0:
                    merger_mass_loss_rate.value = atof(optarg);
                    if(print_flag) merger_mass_loss_rate.print(std::cout);
                    opt_used+=2;
                    break;            
                case 1:
                    merger_radius_amplifier_rate.value = atof(optarg);
                    if(print_flag) merger_radius_amplifier_rate.print(std::cout);
                    opt_used+=2;
                    break;            
                case 2:
                    merger_time_delay.value = atof(optarg);
                    if(print_flag) merger_time_delay.print(std::cout);
                    opt_used+=2;
                    break;
                case 3:
                    equalibrium_mass.value = atof(optarg);
                    if(print_flag) equalibrium_mass.print(std::cout);
                    opt_used+=2;
                    break;
                case 4:
                    initial_mass.value = atof(optarg);
                    if(print_flag) initial_mass.print(std::cout);
                    opt_used+=2;
                    break;
                case 5:
                    growth_time.value = atof(optarg);
                    if(print_flag) growth_time.print(std::cout);
                    opt_used+=2;
                    break;
                default:
                    break;
                }
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

//! Class for managing mergers
class DiskStarMergerManager {
public:
    Float merger_mass_loss_rate; //!< mass loss rate
    Float merger_radius_amplifier_rate; //!< radius amplifier
    Float merger_time_delay; //!< time delay for merger to increase mass
    Float equalibrium_mass; //!< equilibrium mass
    Float initial_mass; //!< initial mass
    Float growth_time; //!< growth timescale
    
    //! Constructor
    DiskStarMergerManager(): merger_mass_loss_rate(0.0), merger_radius_amplifier_rate(1.0), merger_time_delay(0.0), equalibrium_mass(0.0), initial_mass(0.0) {}

    //! (Necessary) check whether publicly initialized parameters are correctly set
    /*! \return true: all parmeters are correct. In this case no parameters, return true;
     */
    bool checkParams() {
        ASSERT(merger_mass_loss_rate>=0.0);
        ASSERT(merger_radius_amplifier_rate>=0.0);
        ASSERT(merger_time_delay>=0.0);
        ASSERT(equalibrium_mass>=0.0);
        ASSERT(initial_mass>=0.0);
        ASSERT(growth_time>=0.0);
        return true;
    }

    //! print parameters
    void print(std::ostream & _fout) const{
        _fout<<"merger_mass_loss_rate : "<<merger_mass_loss_rate<<std::endl
             <<"merger_radius_amplifier_rate : "<<merger_radius_amplifier_rate<<std::endl
             <<"merger_time_delay : "<<merger_time_delay<<std::endl
             <<"equalibrium_mass : "<<equalibrium_mass<<std::endl
             <<"initial_mass : "<<initial_mass<<std::endl
             <<"growth_time : "<<growth_time<<std::endl;
    }

    //! initial parameters for disk star mergers
    /*!
      @param[in] _input: input parameter
      @param[in] _print_flag: printing flag
     */
    void initial(const IOParamsDiskStarMerger& _input, const bool _print_flag=false) {
        merger_mass_loss_rate = _input.merger_mass_loss_rate.value;
        merger_radius_amplifier_rate = _input.merger_radius_amplifier_rate.value;
        merger_time_delay = _input.merger_time_delay.value;
        equalibrium_mass = _input.equalibrium_mass.value;
        initial_mass = _input.initial_mass.value;
        growth_time = _input.growth_time.value;
    }

    //! calcMergerProperties
    /*! Calculate the new properties of the merger
        @param[in,out] p1: particle 1, will be merger 
        @param[in] p2: particle 2, will be zero mass particle
        @param[in] time: current time
        */
    template <class TParticle>
    void calcMergerProperties(TParticle* p1, TParticle* p2, const Float& time) {

        Float mcm = p1->mass + p2->mass;
        for (int k=0; k<3; k++) {
            p1->pos[k] = (p1->mass*p1->pos[k] + p2->mass*p2->pos[k])/mcm;
            p1->vel[k] = (p1->mass*p1->vel[k] + p2->mass*p2->vel[k])/mcm;
        }

        // only increase4mass and change radius after time delay
        if (time>merger_time_delay) {
            Float new_mass = mcm*(1-merger_mass_loss_rate);    
            p1->dm += new_mass - p1->mass;
            p1->mass = new_mass;
            p1->radius = merger_radius_amplifier_rate*p1->radius;
        }

        p2->dm -= p2->mass;
        p2->mass = 0.0;
        p2->radius = 0.0;
    }

    //! calculate mass change
    /*! Calculate mass change, for mass < equalibrium mass, increase mass; for mass > equilibrium mass, decrease mass
        Increase mass formula:  
        dM/dt = c M^2,   c = (Me-Mi)/(2 Me Mi tf),   M(t) = Mi/( 1 - (Me-Mi)/(Me tf) t)
        Mi: initial mass; 
        Me: equilibrium mass; 
        tf: growth timescale;
        @param[in] p: particle
        @param[in] time: current time
        @return new mass
        */
    template <class TParticle>
    void calcMassChange(TParticle* p, const Float& time) {
        Float new_mass = p->mass;
        if (p->stat.type==2) {
            Float dt = time - p->time_record;
            if (dt>0) {
                if (p->mass < equalibrium_mass) {
                    new_mass = initial_mass/(1 - (equalibrium_mass-initial_mass)/(equalibrium_mass*growth_time)*dt);
                    if (new_mass>equalibrium_mass) {
                        new_mass = equalibrium_mass;
                    }
                }
            }
            p->dm += new_mass - p->mass;
            p->mass = new_mass;
        }
    }
    
};

//! class for disk star merger parameters of individual stars
class StarParameter{
public:
    int type; //!< type of star; 0: black hole; 1: star no growth; 2: star with growth
    int merger_star_times; //!< times of merger with star
    int merger_bh_times; //!< times of merger with black hole
    Float growth_time_start; //!< time delay for mass approach equalibrium
    Float last_merger_time; //!< last merger time

    //! Constructor
    StarParameter():  type(-1),
                      merger_star_times(0),
                      merger_bh_times(0),
                      growth_time_start(0.0),
                      last_merger_time(0.0)                    
                      {}


    //! initial parameters for disk star merger
    /*!
      @param[in] _type: type of star; 0: black hole; 1: star no growth; 2: star with growth
      @param[in] _growth_time_start: time delay for mass approach equalibrium
     */
    void initial(const int _type, const Float& _growth_time_start=0.0) {
        type = _type;
        merger_star_times = 0;
        merger_bh_times = 0;
        growth_time_start = _growth_time_start;
        last_merger_time = 0;
    }

        

};
