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
    IOParams<double> mass_loss_rate; //!< mass loss rate
    IOParams<double> radius_amplifier_rate; //!< radius amplifier rate
    IOParams<double> time_delay; //!< time delay for merger to increase mass and change radius

    bool print_flag; //!< print flag
    //! Constructor
    IOParamsDiskStarMerger(): input_par_store(),
                              mass_loss_rate(input_par_store, 0.0, "merger-mass-loss-rate", "mass loss rate for merger"),
                              radius_amplifier_rate(input_par_store, 1.0, "merger-radius-rate", "particle radius amplifier rate after merger"),
                              time_delay(input_par_store, 0.0, "merger-time-delay", "time delay for merger to increase mass"),
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
            {mass_loss_rate.key, required_argument, &merger_flag, 0},  
            {radius_amplifier_rate.key, required_argument, &merger_flag, 1},  
            {time_delay.key, required_argument, &merger_flag, 2},  
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
                    mass_loss_rate.value = atof(optarg);
                    if(print_flag) mass_loss_rate.print(std::cout);
                    opt_used+=2;
                    break;            
                case 1:
                    radius_amplifier_rate.value = atof(optarg);
                    if(print_flag) radius_amplifier_rate.print(std::cout);
                    opt_used+=2;
                    break;            
                case 2:
                    time_delay.value = atof(optarg);
                    if(print_flag) time_delay.print(std::cout);
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
    Float mass_loss_rate; //!< mass loss rate
    Float radius_amplifier_rate; //!< radius amplifier
    Float time_delay; //!< time delay for merger to increase mass
    
    //! Constructor
    DiskStarMergerManager(): mass_loss_rate(0.0), radius_amplifier_rate(1.0), time_delay(0.0) {}

    //! (Necessary) check whether publicly initialized parameters are correctly set
    /*! \return true: all parmeters are correct. In this case no parameters, return true;
     */
    bool checkParams() {
        ASSERT(mass_loss_rate>=0.0);
        ASSERT(radius_amplifier_rate>=0.0);
        ASSERT(time_delay>=0.0);
        return true;
    }

    //! print parameters
    void print(std::ostream & _fout) const{
        _fout<<"mass_loss_rate : "<<mass_loss_rate<<std::endl
             <<"radius_amplifier_rate : "<<radius_amplifier_rate<<std::endl
             <<"time_delay : "<<time_delay<<std::endl;
    }

    //! initial parameters for disk star mergers
    /*!
      @param[in] _input: input parameter
      @param[in] _print_flag: printing flag
     */
    void initial(const IOParamsDiskStarMerger& _input, const bool _print_flag=false) {
        mass_loss_rate = _input.mass_loss_rate.value;
        radius_amplifier_rate = _input.radius_amplifier_rate.value;
        time_delay = _input.time_delay.value;
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

        // only increase mass and change radius after time delay
        if (time>time_delay) {
            Float new_mass = mcm*(1-mass_loss_rate);    
            p1->dm += new_mass - p1->mass;
            p1->mass = new_mass;
            p1->radius = radius_amplifier_rate*p1->radius;
        }

        p2->dm -= p2->mass;
        p2->mass = 0.0;
        p2->radius = 0.0;
    }
    
};