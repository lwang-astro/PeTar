#pragma once

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <string>
#include <getopt.h>
#include "io.hpp"
#include "Common/Float.h"
#include "static_variables.hpp"
#ifdef GALPY
#include "galpy_interface.h"
#include "status.hpp"
#endif

//! IO parameters manager for external perturbation in hard integration
/*! For initializing the COMMON block variables from the commander option.
  The description of each parameter is also provided.
 */
class IOParamsExternalHard{
public:
    IOParamsContainer input_par_store;
    IOParams<long long int> mode; // option to switch perturbation
    IOParams<double> sound_speed; 
    IOParams<double> coulomb_log; 
    IOParams<double> polytropic_constant; // K
    IOParams<double> polytropic_exponent; // gamma
#ifdef GALPY
    IOParams<long long int> galpy_gaspot_index;
    IOParams<double> scale_density;
#else
    IOParams<double> gas_density; 
    IOParams<double> decay_time;
#endif
    IOParams<double> gravitational_constant;
    IOParams<std::string> fname_par;
    
    bool print_flag;
    IOParamsExternalHard(): input_par_store(),
                            mode  (input_par_store, 0,          "ext-hard-mode", "external force for hard integration; 0, not used; 1, gas dynamical friction (Ostriker 1999, Rozner et al. 2022); 2, gas dynamical friction only in radial direction"),
                            sound_speed  (input_par_store, 0.0, "ext-sound-speed",  "sound speed in units of PeTar input, if given 0, calculate by sqrt(P/rho), based on hydrostatic equilibrium"),
                            coulomb_log  (input_par_store, 3.1, "ext-coulomb-log",  "coulomb logarithm"),
                            polytropic_constant (input_par_store, 1.0, "ext-K", "Polytropic constant in units of PeTar input, used to evaluate Pressure P = K rho^gamma"),
                            polytropic_exponent (input_par_store, 4.0/3.0, "ext-gamma", "Polytropic exponent, used to evaluate Pressure P = K rho^gamma"),
#ifdef GALPY
                            galpy_gaspot_index(input_par_store, -1, "ext-gaspot-index",  "galpy potential set index for gas component, used for obtaining gas density", "None"),
                            scale_density(input_par_store, 1/G_ASTRO, "ext-scale-density", "scale factor for galpy potential density","1/G"),
#else
                            gas_density  (input_par_store, 1.0, "ext-gas-density",  "gas density in units of PeTar input"),
                            decay_time   (input_par_store, 0.0, "ext-decay-time",  "gas density decay time scale in units of PeTar input"),
#endif
                            gravitational_constant (input_par_store, 1.0, "G", "Gravitational constant", NULL, false),
                            fname_par    (input_par_store, "input.par", "p", "Input parameter file for external force (this option should be used first before any other options)",NULL,false),
    
                            print_flag(false) {}

    //! reading parameters from GNU option API
    /*!
      @param[in] argc: number of options
      @param[in] argv: string of options
      @param[in] print_format_info: if true, print the format information
      @param[in] opt_used_pre: already used option number from previous reading, use to correctly count the remaining argument number
      \return -1 if help is used; else the used number of argv
     */
    int read(int argc, char *argv[], const bool print_format_info=true, const int opt_used_pre=0) {
        static int ext_flag=-1;
        const struct option long_options[] = {
            {mode.key,    required_argument, &ext_flag, 0},  
#ifdef GALPY
            {galpy_gaspot_index.key, required_argument, &ext_flag, 1},  
            {scale_density.key, required_argument, &ext_flag, 2},  
#else
            {gas_density.key, required_argument, &ext_flag, 1},  
            {decay_time.key,  required_argument, &ext_flag, 2},  
#endif
            {sound_speed.key, required_argument, &ext_flag, 3},  
            {coulomb_log.key, required_argument, &ext_flag, 4},  
            {polytropic_constant.key, required_argument, &ext_flag, 5},
            {polytropic_exponent.key, required_argument, &ext_flag, 6},
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
                switch (ext_flag) {
                case 0:
                    mode.value = atoi(optarg);
                    if(print_flag) mode.print(std::cout);
                    opt_used+=2;
                    break;            
#ifdef GALPY
                case 1:
                    galpy_gaspot_index.value = atoi(optarg);
                    if(print_flag) galpy_gaspot_index.print(std::cout);
                    opt_used+=2;
                    break;            
                case 2:
                    scale_density.value = atof(optarg);
                    if(print_flag) scale_density.print(std::cout);
                    opt_used+=2;
                    break;            
#else
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
#endif
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
                case 5:
                    polytropic_constant.value = atof(optarg);
                    if(print_flag) polytropic_constant.print(std::cout);
                    opt_used+=2;
                    break;
                case 6:
                    polytropic_exponent.value = atof(optarg);
                    if(print_flag) polytropic_exponent.print(std::cout);
                    opt_used+=2;
                    break;
                default:
                    break;
                }
                break;
            case 'G':
                gravitational_constant.value = atof(optarg);
                if(print_flag) gravitational_constant.print(std::cout);
                opt_used += 2;
                assert(gravitational_constant.value>0.0);
                break;
            case 'p':
                fname_par.value = optarg;
                if(print_flag) {
                    std::string fgalpy_par = fname_par.value+".exthard"; 
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
                    std::cout<<"----- External perturbation for hard integration options: -----"<<std::endl;
                    input_par_store.printHelp(std::cout, print_format_info);
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
    int mode; 
#ifdef GALPY
    int galpy_gaspot_index; 
    GalpyManager* galpy_manager;
    Status* status;
    Float scale_density;
#else
    Float gas_density; 
    Float gas_density_init; 
    Float decay_time;
    Float time;
#endif
    Float sound_speed; 
    Float coulomb_log;
    Float polytropic_constant;
    Float polytropic_exponent;
    Float gravitational_constant;
    bool calc_sound_speed;

    ExternalHardForce(): mode(0),
#ifdef GALPY
                         galpy_gaspot_index(-1), galpy_manager(NULL), status(NULL), scale_density(1.0),
#else
                         gas_density(1.0), gas_density_init(1.0), decay_time(0.0), time(0.0),
#endif
                         sound_speed(0.0), coulomb_log(3.1), polytropic_constant(1.0), polytropic_exponent(4.0/3.0), gravitational_constant(1.0), calc_sound_speed(true) {}

#ifdef GALPY
    //! initial parameters for perturbation
    /*!
      @param[in] _input: input parameter
      @param[in] _galpy_manager: galpy manager pointer
      @param[in] _status: system status for information of time and pcm position and velocity offsets used for converting to galactic frame
      @param[in] _print_flag: printing flag
     */
    void initial(const IOParamsExternalHard& _input, GalpyManager& _galpy_manager, Status& _status, const bool _print_flag=false) {
        mode = _input.mode.value;
        galpy_gaspot_index = _input.galpy_gaspot_index.value;
        galpy_manager = &_galpy_manager;
        status = &_status;
        scale_density = _input.scale_density.value;
        sound_speed = _input.sound_speed.value;
        coulomb_log = _input.coulomb_log.value;
        polytropic_constant = _input.polytropic_constant.value;
        polytropic_exponent = _input.polytropic_exponent.value;
        gravitational_constant = _input.gravitational_constant.value;
        if (sound_speed>0.0) calc_sound_speed = false;
        else calc_sound_speed = true;
    }

#else
    //! initial parameters for perturbation
    void initial(const IOParamsExternalHard& _input, const Float _time, const bool _print_flag=false) {
        mode = _input.mode.value;
        gas_density_init = _input.gas_density.value;
        decay_time  = _input.decay_time.value;
        sound_speed = _input.sound_speed.value;
        coulomb_log = _input.coulomb_log.value;
        polytropic_constant = _input.polytropic_constant.value;
        polytropic_exponent = _input.polytropic_exponent.value;
        gravitational_constant = _input.gravitational_constant.value;
        if (sound_speed>0.0) calc_sound_speed = false;
        else calc_sound_speed = true;
        updateTime(_time);
    }

    //! update time and gas density
    void updateTime(const Float _time) {
        time = _time;
        gas_density = gas_density_init * exp(-time/decay_time);
    }
#endif

    //! External force for one particle in hard part
    /*!
      Gas dynamical friction
      Due to the acceleration dependence, this function must be used at the end of acceleration calculation
      (Ostriker 1999, https://ui.adsabs.harvard.edu/abs/1999ApJ...513..252O, 
      Rozner 2022, https://arxiv.org/abs/2212.00807)

      @param[out] _acc: acceleration
      @param[in] _particle: particle data
      
      Return: the next integration time step (default: maximum floating point number)
    */
    template<class Tp> 
    Float calcAccExternal(Float* _acc, const Tp& _particle){
        if (mode==0) 
            return NUMERIC_FLOAT_MAX;

        auto& mass = _particle.mass;
        auto& vel = _particle.vel;
        Float v2 = vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2];
        Float v = std::sqrt(v2);
        Float v3 = v2*v;

        const Float PI = 4.0*atan(1.0);
        const Float G2 = gravitational_constant*gravitational_constant;

#ifdef GALPY
        auto& pos = _particle.pos;
        Float pos_g[3] = {pos[0] + status->pcm.pos[0], 
                          pos[1] + status->pcm.pos[1], 
                          pos[2] + status->pcm.pos[2]};
        Float gas_density = scale_density*galpy_manager->calcSetDensity(galpy_gaspot_index, status->time, pos_g, &pos[0]);
#else
        auto& pos_g = _particle.pos;
#endif       
        if (calc_sound_speed) 
            sound_speed = std::sqrt(polytropic_constant*std::pow(gas_density, polytropic_exponent-1));

        Float mach = v/sound_speed;
        Float Ifunc;
        if (mach<0.9) {
            Ifunc = 0.5*std::log((1.0+mach)/(1.0-mach)) - mach;
        }
        else if (mach>=0.9 && mach<1.1) {
            // 2nd order derivative Hermite interpolation
            Float mach2 = mach*mach;
            Float mach3 = mach2*mach;
            Float mach4 = mach2*mach2;
            Float mach5 = mach4*mach;
            Ifunc = 8670.66512394438*mach5 - 43353.850000357*mach4 + 86337.1029009788*mach3 - 85594.6848365398*mach2 + 42251.2454762294*mach - 8309.08207013687;
        }
        else{
            Float mach2 = mach*mach;
            Ifunc = 0.5*std::log(1-1/mach2) + coulomb_log;
        }

        Float c1 = -4*PI*G2*mass*gas_density/v3*Ifunc;

        if (mode==1) {
            // GDF force
            _acc[0] += c1*vel[0];
            _acc[1] += c1*vel[1];
            _acc[2] += c1*vel[2];
        }
        else if (mode==2) {
            // GDF force only in radial direction
            Float r_g = std::sqrt(pos_g[0]*pos_g[0] + pos_g[1]*pos_g[1] + pos_g[2]*pos_g[2]);
            _acc[0] += c1*vel[0]*pos_g[0]/r_g;
            _acc[1] += c1*vel[1]*pos_g[1]/r_g;
            _acc[2] += c1*vel[2]*pos_g[2]/r_g;
        }

        return NUMERIC_FLOAT_MAX;
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
        if (mode==0) 
            return NUMERIC_FLOAT_MAX;

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
        const Float G2 = gravitational_constant*gravitational_constant;

#ifdef GALPY
        auto& pos = _particle.pos;
        Float pos_g[3] = {pos[0] + status->pcm.pos[0], 
                          pos[1] + status->pcm.pos[1], 
                          pos[2] + status->pcm.pos[2]};
        Float gas_density = scale_density*galpy_manager->calcSetDensity(galpy_gaspot_index, status->time, pos_g, &pos[0]);
#else
        auto& pos_g = _particle.pos;
#endif

        if (calc_sound_speed) 
            sound_speed = std::sqrt(polytropic_constant*std::pow(gas_density, polytropic_exponent-1));

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
        else{
            Float mach2 = mach*mach;
            Ifunc = 0.5*std::log(1-1/mach2) + coulomb_log;
            dIfunc = 1/(mach2*mach - mach);
        }

        Float r_g;
        Float c1 = -4*PI*G2*mass*gas_density/v3*Ifunc;

        auto& acc0 = _force.acc0;
        //Float acc0[3] = {c1*vel[0], c1*vel[1], c1*vel[2]};

        if (mode==1) {
            // GDF force
            acc0[0] += c1*vel[0];
            acc0[1] += c1*vel[1];
            acc0[2] += c1*vel[2];
        }
        else{
            // GDF force only in radial direction
            r_g = std::sqrt(pos_g[0]*pos_g[0] + pos_g[1]*pos_g[1] + pos_g[2]*pos_g[2]);
            acc0[0] += c1*vel[0]*pos_g[0]/r_g;
            acc0[1] += c1*vel[1]*pos_g[1]/r_g;
            acc0[2] += c1*vel[2]*pos_g[2]/r_g;
        }

        Float vdota = vel[0]*acc0[0] + vel[1]*acc0[1] + vel[2]*acc0[2];
        // d(1/v^3)/dt = -3/v^5 v dot a
        Float c2 = -3*c1/v2*vdota;
        // d(v/ds)/dt  = v dot a / (v*ds) 
        Float c3 = -c1*dIfunc*vdota/(v*sound_speed);
        
        //_force.acc0[0] += acc0[0];
        //_force.acc0[1] += acc0[1];
        //_force.acc0[2] += acc0[2];

        if (mode==1) {
            // GDF force derivative
            _force.acc1[0] += (c2+c3)*vel[0] + c1*acc0[0];
            _force.acc1[1] += (c2+c3)*vel[1] + c1*acc0[1];
            _force.acc1[2] += (c2+c3)*vel[2] + c1*acc0[2];
        }
        else if (mode==2) {
            // GDF force derivative only in radial direction, assuming pos_g is constant. For time-dependent pos_g, additional term of time derivative of pos_g should be added.
            _force.acc1[0] += ((c2+c3)*vel[0] + c1*acc0[0])*pos_g[0]/r_g;
            _force.acc1[1] += ((c2+c3)*vel[1] + c1*acc0[1])*pos_g[1]/r_g;
            _force.acc1[2] += ((c2+c3)*vel[2] + c1*acc0[2])*pos_g[2]/r_g;
        }

        return NUMERIC_FLOAT_MAX;
    }

    //! check whether parameters values are correct
    /*! \return true: all correct
     */
    bool checkParams() {
        ASSERT(mode>=0 && mode<=2);
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
