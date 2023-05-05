#pragma once

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <string>
#include <getopt.h>
#include "io.hpp"
#include "Common/Float.h"
#include "static_variables.hpp"

//! IO parameters manager for external perturbation in hard integration
/*! For initializing the COMMON block variables from the commander option.
  The description of each parameter is also provided.
 */
class IOParamsExternalHard{
public:
    IOParamsContainer input_par_store;
    IOParams<long long int> mode; // option to switch perturbation
    IOParams<double> gas_density; 
    IOParams<double> decay_time;
    IOParams<double> sound_speed; 
    IOParams<double> coulomb_log; 

    bool print_flag;
    IOParamsExternalHard(): input_par_store(),
                            mode  (input_par_store, 0,          "ext-hard-mode", "external hard mode: 0, not used; 1, gas dynamical friction (Ostriker 1999, Rozner et al. 2022)"),
                            gas_density  (input_par_store, 1.0, "ext-gas-density",  "gas density in units of PeTar input"),
                            decay_time   (input_par_store, 0.0, "ext-decay-time",  "gas density decay time scale in units of PeTar input"),
                            sound_speed  (input_par_store, 1.0, "ext-sound-speed",  "sound speed in units of PeTar input"),
                            coulomb_log  (input_par_store, 3.1, "ext-coulomb-log",  "coulomb logarithm"),
                            print_flag(false) {}

    //! reading parameters from GNU option API
    /*!
      @param[in] argc: number of options
      @param[in] argv: string of options
      @param[in] opt_used_pre: already used option number from previous reading, use to correctly count the remaining argument number
      \return -1 if help is used; else the used number of argv
     */
    int read(int argc, char *argv[], const int opt_used_pre=0) {
        static int ext_flag=-1;
        const struct option long_options[] = {
            {mode.key,    required_argument, &ext_flag, 0},  
            {gas_density.key, required_argument, &ext_flag, 1},  
            {decay_time.key,  required_argument, &ext_flag, 2},  
            {sound_speed.key, required_argument, &ext_flag, 3},  
            {coulomb_log.key, required_argument, &ext_flag, 4},  
            {"help",      no_argument,       0, 'h'},
            {0,0,0,0}
        };

        int opt_used=opt_used_pre;
        int copt;
        int option_index;
        std::string fname_par;
        optind = 0;
        while ((copt = getopt_long(argc, argv, "-z:p:h", long_options, &option_index)) != -1) 
            switch (copt) {
            case 0:
                switch (ext_flag) {
                case 0:
                    mode.value = atoi(optarg);
                    if(print_flag) mode.print(std::cout);
                    opt_used+=2;
                    break;            
                case 1:
                    gas_density.value = atof(optarg);
                    if(print_flag) gas_density.print(std::cout);
                    opt_used+=2;
                    break;            
                case 2:
                    decay_time.value = atof(optarg);
                    if(print_flag) decay_time.print(std::cout);
                    opt_used+=2;
                    break;            
                case 3:
                    sound_speed.value = atof(optarg);
                    if(print_flag) sound_speed.print(std::cout);
                    opt_used+=2;
                    break;            
                case 4:
                    coulomb_log.value = atof(optarg);
                    if(print_flag) coulomb_log.print(std::cout);
                    opt_used+=2;
                    break;            
                default:
                    break;
                }
                break;
            case 'p':
                fname_par = optarg;
                if(print_flag) {
                    std::string fgalpy_par = fname_par+".exthard"; 
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
                    std::cout<<"External perturbation for hard integration, options:"<<std::endl;
                    input_par_store.printHelp(std::cout, 2, 10, 23);
                    std::cout<<std::endl;
                }
                return -1;
            case '?':
                opt_used +=2;
                break;
            default:
                break;
            }
        
        if(print_flag) std::cout<<"----- Finish reading input options of external perturbation for hard integration -----\n";

        return opt_used;
    }    
};

class ExternalHardForce{
public:
    bool is_used; // indicator whether external force is used
    Float gas_density; 
    Float gas_density_init; 
    Float decay_time;
    Float sound_speed; 
    Float coulomb_log;
    Float time;

    ExternalHardForce(): is_used(false), gas_density(1.0), gas_density_init(1.0), decay_time(0.0), sound_speed(1.0), coulomb_log(3.1), time(0.0) {}

    //! initial parameters for perturbation
    void initial(const IOParamsExternalHard& _input, const Float _time, const bool _print_flag=false) {
        is_used = bool(_input.mode.value);
        gas_density_init = _input.gas_density.value;
        decay_time  = _input.decay_time.value;
        sound_speed = _input.sound_speed.value;
        coulomb_log = _input.coulomb_log.value;
        updateTime(_time);
    }

    //! update time and gas density
    void updateTime(const Float _time) {
        time = _time;
        gas_density = gas_density_init * exp(-time/decay_time);
    }

    //! External force for one particle in hard part
    /*!
      Gas dynamical friction
      Due to the acceleration dependence, this function must be used at the end of acceleration calculation
      (Ostriker 1999, https://ui.adsabs.harvard.edu/abs/1999ApJ...513..252O, 
      Rozner 2022, https://arxiv.org/abs/2212.00807)

      @param[out] _force: acceleration and jerk
      @param[in] _particle: particle data
      
      Return: the next integration time step (default: maximum floating point number)
    */
    template<class Tf, class Tp> 
    Float calcAccJerkExternal(Tf& _force, const Tp& _particle){
        
        auto& mass = _particle.mass;
        auto& vel = _particle.vel;
        Float v2 = vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2];
        Float v = std::sqrt(v2);
        Float v3 = v2*v;

        //auto& pos = _particle.pos;
        //Float r2 = pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2] + sound_speed;
        //Float c1 = alpha/std::pow(r2, 0.5*beta);

        //Float c2 = beta*(pos[0]*vel[0]+pos[1]*vel[1]+pos[2]*vel[2])*c1/r2;
        //_force.acc1[0] += c2*vel[0] - c1*_force.acc0[0];
        //_force.acc1[1] += c2*vel[1] - c1*_force.acc0[1];
        //_force.acc1[2] += c2*vel[2] - c1*_force.acc0[2];
        
        const Float PI = 4.0*atan(1.0);
        Float G2 = ForceSoft::grav_const*ForceSoft::grav_const;
        Float mach = v/sound_speed;
        Float Ifunc, dIfunc;
        if (mach<0.9) {
            Float mach2 = mach*mach;
            Ifunc = 0.5*std::log((1.0+mach)/(1.0-mach)) - mach;
            dIfunc = mach2/(1-mach2);
        }
        else if (mach>=0.9 && mach<1.1) {
            // 2nd order derivative Hermite interpolation
            Float mach2 = mach*mach;
            Float mach3 = mach2*mach;
            Float mach4 = mach2*mach2;
            Float mach5 = mach4*mach;
            Ifunc = 8670.66512394438*mach5 - 43353.850000357*mach4 + 86337.1029009788*mach3 - 85594.6848365398*mach2 + 42251.2454762294*mach - 8309.08207013687;
            dIfunc = 43353.3256197219*mach4 - 173415.400001428*mach3 + 259011.308702936*mach2 - 171189.36967308*mach + 42251.2454762294;
        }
        else if (mach>=1.1) {
            Float mach2 = mach*mach;
            Ifunc = 0.5*std::log(1-1/mach2) + coulomb_log;
            dIfunc = 1/(mach2*mach - mach);
        }

        Float c1 = -4*PI*G2*mass*gas_density/v3*Ifunc;

        auto& acc0 = _force.acc0;
        //Float acc0[3] = {c1*vel[0], c1*vel[1], c1*vel[2]};
        acc0[0] += c1*vel[0];
        acc0[1] += c1*vel[1];
        acc0[2] += c1*vel[2];

        Float vdota = vel[0]*acc0[0] + vel[1]*acc0[1] + vel[2]*acc0[2];
        // d(1/v^3)/dt = -3/v^5 v dot a
        Float c2 = -3*c1/v2*vdota;
        // d(v/ds)/dt  = v dot a / (v*ds) 
        Float c3 = -c1*dIfunc*vdota/(v*sound_speed);
        
        //_force.acc0[0] += acc0[0];
        //_force.acc0[1] += acc0[1];
        //_force.acc0[2] += acc0[2];

        _force.acc1[0] += (c2+c3)*vel[0] + c1*acc0[0];
        _force.acc1[1] += (c2+c3)*vel[1] + c1*acc0[1];
        _force.acc1[2] += (c2+c3)*vel[2] + c1*acc0[2];

        return NUMERIC_FLOAT_MAX;
    }

    //! check whether parameters values are correct
    /*! \return true: all correct
     */
    bool checkParams() {
        return true;
    }    

    //! print parameters
    void print(std::ostream & _fout) const{
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
