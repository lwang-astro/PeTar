#pragma once

#include <getopt.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include "io.hpp"
#include "Common/Float.h"
#include "static_variables.hpp"
#ifdef GAS_DRAG
#include "gas_drag.hpp"
#endif

//! IO parameters manager for external perturbation in hard integration
/*! Stub implementation when no specific external hard-force model is enabled. */
class IOParamsExternalHard{
public:
    IOParamsContainer input_par_store;
    IOParams<long long int> switcher;
    IOParams<std::string> fname_par;
#ifdef GAS_DRAG
    IOParamsGasDrag gas_drag;
#endif
    bool print_flag;

    //! constructor
    IOParamsExternalHard() :
        input_par_store(),
        switcher(input_par_store, 1, "ext-hard-switch", "switch of external hard force; 0: off, 1: on"),
        fname_par(input_par_store, "input.par", "p", "Input parameter file for external force (this option should be used first before any other options)",NULL,false),
#ifdef GAS_DRAG
        gas_drag(),
#endif
        print_flag(false)
    {}

    //! append parameter containers for undefined-option checking
    /*! 
      @param[in,out] _par_list: list of parameter containers
     */
    void appendInputParamStores(std::vector<IOParamsContainer*>& _par_list) {
        _par_list.push_back(&input_par_store);
#ifdef GAS_DRAG
        _par_list.push_back(&gas_drag.input_par_store);
#endif
    }

    //! set gravitational constant for model parameters
    /*! 
      @param[in] _grav_const: gravitational constant
     */
    void setGravitationalConstant(const Float _grav_const) {
#ifdef GAS_DRAG
        gas_drag.gravitational_constant.value = _grav_const;
#endif
    }

    //! get active external hard-force feature name
    /*! 
      
eturn feature name string
     */
    const char* getFeatureName() const {
#ifdef GAS_DRAG
        return "gasdrag";
#else
        return "unknown";
#endif
    }

    //! write external hard parameter sets in ASCII format
    /*! 
      @param[in] _filename_prefix: output filename prefix; the actual filename will be _filename_prefix.exthard
    */
    void writeModelParamsAscii(const std::string& _filename_prefix) {
        std::string fexthard_par = _filename_prefix + ".exthard";
        if (print_flag) std::cout << "Save external_hard_parameters to file " << fexthard_par << std::endl;
        FILE* fpar_out;
        if ((fpar_out = fopen(fexthard_par.c_str(), "w")) == NULL) {
            fprintf(stderr, "Error: Cannot open file %s.\n", fexthard_par.c_str());
            abort();
        }
        input_par_store.writeAscii(fpar_out);
        fclose(fpar_out);
#ifdef GAS_DRAG
        gas_drag.writeModelParamsAscii(_filename_prefix);
#endif
    }

    //! reading parameters from GNU option API
    /*! 
        @param[in] argc: number of options
        @param[in] argv: string of options
        @param[in] print_format_info: if true, print the format information
        @param[in] opt_used_pre: already used option number from previous reading
        
        @return -1 if help is used; else the used number of argv
    */
    int read(int argc, char* argv[], const bool print_format_info=true, const int opt_used_pre=0) {
        int opt_used = opt_used_pre;
        static int ext_flag=-1;
        const struct option long_options[] = {
            {switcher.key, required_argument, &ext_flag, 0},
            {fname_par.key, required_argument, 0, 'p'},
            {"help", no_argument, 0, 'h'},
            {0,0,0,0}
        };

        bool helf_option_used = false;
        int copt;
        int option_index;
        optind = 0;
        while ((copt = getopt_long(argc, argv, "-p:h", long_options, &option_index)) != -1) {
            switch (copt) {
            case 0:
                if (ext_flag==0) {
                    switcher.value = atoll(optarg);
                    if(print_flag) switcher.print(std::cout);
                    opt_used += 2;
                }
                break;
            case 'p':
                fname_par.value = optarg;
                if(print_flag) {
                    std::string fpar_name = fname_par.value+".exthard";
                    FILE* fpar_in;
                    if( (fpar_in = fopen(fpar_name.c_str(),"r")) == NULL) {
                        fprintf(stderr,"Error: Cannot open file %s.\n", fpar_name.c_str());
                        abort();
                    }
                    input_par_store.readAscii(fpar_in);
                    fclose(fpar_in);
                }
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                input_par_store.mpi_broadcast();
                PS::Comm::barrier();
#endif
                break;
            case 'h':
                if (print_flag) {
                    std::cout<<"----- External perturbation for hard integration options: -----"<<std::endl;
                    input_par_store.printHelp(std::cout, print_format_info);
                }
                helf_option_used = true;
                break;
            case '?':
                opt_used += 2;
                break;
            default:
                break;
            }
        }

        if (helf_option_used) opt_used = -1;

#ifdef GAS_DRAG
        gas_drag.print_flag = print_flag;
        opt_used = gas_drag.read(argc, argv, print_format_info, opt_used);
#endif
        
        return opt_used;
    }
};

//! External force manager in hard integration
/*! Stub implementation when no specific external hard-force model is enabled. */
class ExternalHardForce{
public:
    bool enabled;
    FPSoft* center;
    PS::S64* center_id;
#ifdef GAS_DRAG
    GasDragForce gas_drag;
#endif

    //! constructor
    ExternalHardForce()
        : enabled(false)
        , center(NULL)
        , center_id(NULL)
#ifdef GAS_DRAG
        , gas_drag()
#endif
    {}

    //! bind center information used by external hard force
    /*! 
      @param[in] _center: center particle reference
      @param[in] _center_id: center particle id reference
     */
    void bindCenter(FPSoft& _center, PS::S64& _center_id) {
        center = &_center;
        center_id = &_center_id;
    }

#ifdef GALPY
    //! initialize external hard force with galpy support
    /*! 
      @param[in] _input: external hard input parameters
      @param[in] _galpy_manager: galpy manager
      @param[in] _status: system status
      @param[in] _center: center particle
      @param[in] _center_id: center particle id
      @param[in] _print_flag: print initialization information
     */
    template<class TGalpyManager, class TStatus>
    void initial(const IOParamsExternalHard& _input, TGalpyManager& _galpy_manager, TStatus& _status, FPSoft& _center, PS::S64& _center_id, bool _print_flag=false) {
        bindCenter(_center, _center_id);
        enabled = (_input.switcher.value>0);
#ifdef GAS_DRAG
        gas_drag.initial(_input.gas_drag, _galpy_manager, _status, _print_flag);
        enabled = enabled && (gas_drag.mode>0);
#else
        enabled = false;
#endif
    }
#else
    //! initialize external hard force without galpy support
    /*! 
      @param[in] _input: external hard input parameters
      @param[in] _time: current time
      @param[in] _center: center particle
      @param[in] _center_id: center particle id
      @param[in] _print_flag: print initialization information
     */
    void initial(const IOParamsExternalHard& _input, const Float _time, FPSoft& _center, PS::S64& _center_id, const bool _print_flag=false) {
        bindCenter(_center, _center_id);
        enabled = (_input.switcher.value>0);
#ifdef GAS_DRAG
        gas_drag.initial(_input.gas_drag, _time, _print_flag);
        enabled = enabled && (gas_drag.mode>0);
#else
        enabled = false;
#endif
        checkParams();
    }

    //! update time-dependent model state
    /*! 
      @param[in] _time: current time
     */
    void updateTime(const Float _time) {
#ifdef GAS_DRAG
        if (enabled) gas_drag.updateTime(_time);
#endif
    }
#endif

    //! check whether external hard force is enabled at runtime
    bool isEnabled() const {
        return enabled;
    }

    //! calculate external acceleration and jerk
    /*! 
      @param[out] _acc0: acceleration
      @param[out] _acc1: jerk
      @param[in] _particle: target particle
      @param[in] _calc_acc1: if true, calculate jerk
      
eturn time step criterion from model; NUMERIC_FLOAT_MAX if disabled
     */
    template<class Tp>
    Float calcAccJerkExternal(Float* _acc0, Float* _acc1, const Tp& _particle, const bool _calc_acc1) {
        if (!enabled) return NUMERIC_FLOAT_MAX;
        assert(center!=NULL);
        assert(center_id!=NULL);
#ifdef GAS_DRAG
        return gas_drag.calcAccJerkExternal(_acc0, _acc1, _particle, _calc_acc1, *center, *center_id);
#else
        return NUMERIC_FLOAT_MAX;
#endif
    }

    //! check external hard force parameter consistency
    bool checkParams() {
        assert(center!=NULL);
        assert(center_id!=NULL);
#ifdef GAS_DRAG
        return gas_drag.checkParams();
#else
        return true;
#endif
    }

    //! print runtime status and model information
    /*! 
      @param[in,out] _fout: output stream
     */
    void print(std::ostream & _fout) const{
        _fout<<"external hard enabled: "<<enabled<<std::endl;
#ifdef GAS_DRAG
        if (enabled) gas_drag.print(_fout);
        if (enabled && center!=NULL && center_id!=NULL) {
            _fout<<"center id: "<<(*center_id)<<std::endl
                 <<"center mass: "<<center->mass<<std::endl;
        }
#endif
    }

    //! write binary data
    /*! 
      @param[in] _fp: output file pointer
     */
    void writeBinary(FILE *_fp) const {
        fwrite(&enabled, sizeof(enabled), 1, _fp);
#ifdef GAS_DRAG
        gas_drag.writeBinary(_fp);
#endif
    }

    //! write binary columns to stream
    /*! 
      @param[in,out] _fout: output stream
     */
    void printColumnBinary(std::ostream& _fout) const {
        _fout.write(reinterpret_cast<const char*>(&enabled), sizeof(enabled));
#ifdef GAS_DRAG
        gas_drag.printColumnBinary(_fout);
#endif
    }

    //! read binary data from C file API
    /*! 
      @param[in] _fin: input file pointer
     */
    void readBinary(FILE *_fin) {
        size_t rcount = fread(&enabled, sizeof(enabled), 1, _fin);
        if (rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            abort();
        }
#ifdef GAS_DRAG
        gas_drag.readBinary(_fin);
#endif
    }

    //! read binary data from C++ stream API
    /*! 
      @param[in,out] _fin: input stream
     */
    void readBinary(std::istream& _fin) {
        _fin.read(reinterpret_cast<char*>(&enabled), sizeof(enabled));
        if (!_fin) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1.\n";
            abort();
        }
#ifdef GAS_DRAG
        gas_drag.readBinary(_fin);
#endif
    }
};
