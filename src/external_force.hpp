#pragma once

#include "Common/Float.h"
#ifdef GAS_DRAG
#include "gas_drag.hpp"
#endif

//! IO parameters manager for external perturbation in hard integration
/*! For initializing the COMMON block variables from the commander option.
  The description of each parameter is also provided.
 */
class IOParamsExternalHard{
public:
#ifdef GAS_DRAG
    IOParamsGasDrag gas_io;
#endif    
    bool print_flag;

    IOParamsExternalHard(): 
#ifdef GAS_DRAG
        gas_io(), 
#endif
        print_flag(false) {}

    //! reading parameters from GNU option API
    /*!
      @param[in] argc: number of options
      @param[in] argv: string of options
      @param[in] opt_used_pre: already used option number from previous reading, use to correctly count the remaining argument number
      \return -1 if help is used; else the used number of argv
     */
    int read(int argc, char *argv[], const int opt_used_pre=0) {
        int opt_used = 0;
#ifdef GAS_DRAG
        gas_io.print_flag = print_flag;
        opt_used = gas_io.read(argc, argv, std::max(opt_used,0));
#endif
        return opt_used;
    }

    //! write all input_par_stores;
    void writeParameters(FILE *_fp) {
#ifdef GAS_DRAG
        gas_io.input_par_store.writeAscii(_fp);
#endif
    }

};

//! External force for one particle used in Hard part
class ExternalHardForce{
public:
    int mode; // Indicate which external force is used. 1: gas drag
#ifdef GAS_DRAG
    GasDrag gas_drag;
#endif

    ExternalHardForce(): mode(0)
#ifdef GAS_DRAG
                       , gas_drag() 
#endif
    {}

    //! initial parameters for perturbation
    /*!
      @param[in] _input: input parameter
      @param[in] _print_flag: printing flag
     */
    void initial(const IOParamsExternalHard& _input, const bool _print_flag=false) {
#ifdef GAS_DRAG
        gas_drag.initial(_input.gas_io, _print_flag);
        // determine mode
        if (gas_drag.mode>0) mode = 1;
#endif
    }

    //! External force for one particle in hard part
    /*!
      call external module to calculate the acceleration

      @param[out] _acc: acceleration
      @param[in] _particle: particle data

      Return: the next integration time step (default: maximum floating point number)
    */
    template<class Tp> 
    Float calcAccExternal(Float* _acc, const Tp& _particle){
        Float dt = NUMERIC_FLOAT_MAX;
#ifdef GAS_DRAG
        dt = gas_drag.calcAccExternal(_acc, _particle);
#endif
        return dt;
    }


    //! External force for one particle in hard part
    /*!
      @param[out] _force: acceleration and jerk
      @param[in] _particle: particle data
      
      Return: the next integration time step (default: maximum floating point number)
    */
    template<class Tf, class Tp> 
    Float calcAccJerkExternal(Tf& _force, const Tp& _particle){
        Float dt = NUMERIC_FLOAT_MAX;
#ifdef GAS_DRAG
        dt = gas_drag.calcAccJerkExternal(_force, _particle);
#endif
        return dt;
    }

    //! check whether parameters values are correct
    /*! \return true: all correct
     */
    bool checkParams() {
#ifdef GAS_DRAG
        ASSERT(gas_drag.checkParams());
#endif
        return true;
    }    

    //! print parameters
    void print(std::ostream & _fout) const{
#ifdef GAS_DRAG
        gas_drag.print(_fout);
#endif
    }    

    //! write class data to file with binary format
    /*! @param[in] _fp: FILE type file for output
     */
    void writeBinary(FILE *_fp) const {
#ifdef GAS_DRAG
        gas_drag.writeBinary(_fp);
#endif
    }

    //! read class data to file with binary format
    /*! @param[in] _fp: FILE type file for reading
     */
    void readBinary(FILE *_fin) {
#ifdef GAS_DRAG
        gas_drag.readBinary(_fin);
#endif
    }        

};
