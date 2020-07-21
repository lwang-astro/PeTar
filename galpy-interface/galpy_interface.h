#include <iostream>
#include <iomanip>
#include <cassert>
#include <string>
#include <vector>
#include <getopt.h>
#include "../src/io.hpp"

#include <integrateFullOrbit.h>

//! IO parameters for Galpy manager
class IOParamsGalpy{
public:
    IOParamsContainer input_par_store;
    IOParams<std::string> type_args; // potential type and argument list
    IOParams<std::string> pre_define_type; // pre defined type name
    IOParams<std::string> pos_offset; // offset of the system position
    
    bool print_flag;

    IOParamsGalpy(): input_par_store(),
                     type_args (input_par_store, "", "Description of potential type and arguments"),
                     pre_define_type (input_par_store, "", "Add additional Pre-defined potential type to potential list, options are: MWPotential2014"),
                     pos_offset(input_par_store, "0:0:0", "3D position offset to obtain the particle position in the potential field, x, y, z are separated by ':'"),
                     print_flag(false) {}

    //! reading parameters from GNU option API
    /*!
      @param[in] argc: number of options
      @param[in] argv: string of options
      @param[in] opt_used_pre: already used option number from previous reading, use to correctly count the remaining argument number
      \return -1 if help is used; else the used number of argv
     */
    int read(int argc, char *argv[], const int opt_used_pre=0) {
        static int galpy_flag=-1;
        const struct option long_options[] = {
            {"type-arg", required_argument, &galpy_flag, 0},  
            {"pos-offset", required_argument, &galpy_flag, 1}, 
            {"pre-define-type", required_argument, &galpy_flag, 2}, 
            {"help", no_argument, 0, 'h'},
            {0,0,0,0}
        };
        
        int opt_used=0;
        int copt;
        int option_index;
        std::string fname_par;
        optind = 0;
        while ((copt = getopt_long(argc, argv, "-p:h", long_options, &option_index)) != -1) 
            switch (copt) {
            case 0:
                switch (galpy_flag) {
                case 0:
                    type_args.value = optarg;
                    if (print_flag) type_args.print(std::cout);
                    opt_used+=2;
                    break;
                case 1:
                    pos_offset.value = optarg;
                    if (print_flag) pos_offset.print(std::cout);
                    opt_used+=2;
                    break;
                case 2:
                    pre_define_type.value = optarg;
                    if (print_flag) pre_define_type.print(std::cout);
                    opt_used+=2;
                    break;
                default:
                    break;
                }
                break;
            case 'p':
                fname_par = optarg;
                if(print_flag) {
                    std::string fgalpy_par = fname_par+".galpy"; 
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
                    std::cout<<"Galpy options:"<<std::endl;
                    std::cout<<"       Option defaulted values are shown after ':'"<<std::endl;
                    std::cout<<"         --type-arg:  [S] "<<type_args<<std::endl
                             <<"         --pos-offset:[S] "<<pos_offset<<std::endl
                             <<"         --pre-define-type:[S] "<<pre_define_type<<std::endl
                             <<"  Potential type index and corresponding arguments can be found by the Python commander:\n"
                             <<"\n"
                             <<"    from galpy.orbit.integrateFullOrbit import _parse_pot\n"
                             <<"    npot, pot_type, pot_args= _parse_pot(pot_name)\n"
                             <<"  where pot_name is the name of potential, see details of Galpy document\n"
                             <<"  The format of a combination of types and arguments is (no space in the middle):\n"
                             <<"    [type1 index]:[arg1],[arg2],**|[type2 index]:[arg1],[arg2],**\n"
                             <<"  For example:\n"
                             <<"    --type-arg 15:0.0299946,1.8,0.2375|5:0.7574802,0.375,0.035|9:4.85223053,2.0\n"
                             <<"  generate the MWPotential2014 (same to --pre-define-type MWPotential2014) from Galpy\n";
                }
                return -1;
            case '?':
                break;
            default:
                break;
            }
        
        if(print_flag) std::cout<<"----- Finish reading input options of Galpy -----\n";

        return opt_used;
    }    
};

//! A class to manager the API to Galpy
class GalpyManager{
public:
    struct potentialArg* potential_args;
    int npot; // number of potential models 
    double pos_offset[3]; // position offset

    GalpyManager(): potential_args(NULL), npot(0), pos_offset{0.0,0.0,0.0} {}

    //! initialization function
    void initial(const IOParamsGalpy& _input, const bool _print_flag=false) {
        std::size_t istart = 0;
        const std::string& pos_offset_str=_input.pos_offset.value;
        std::size_t inext = pos_offset_str.find_first_of(":");
        if (inext==std::string::npos) {
            std::cerr<<"Error: position offset delimiter ':' is not found in input string "<<pos_offset_str<<std::endl;
            abort();
        }
        pos_offset[0] = std::stod(pos_offset_str.substr(0,inext));
        istart = inext + 1;
        inext = pos_offset_str.find_first_of(":",istart);
        if (inext==std::string::npos) {
            std::cerr<<"Error: 2nd position offset delimiter ':' is not found in input string "<<pos_offset_str<<std::endl;
            abort();
        }
        pos_offset[1] = std::stod(pos_offset_str.substr(istart,inext-istart));
        pos_offset[2] = std::stod(pos_offset_str.substr(inext+1));
        if (_print_flag) std::cout<<"Position offset: "<<pos_offset[0]<<" "<<pos_offset[1]<<" "<<pos_offset[2]<<std::endl;
        
        std::string type_args = _input.type_args.value;
        if (_input.pre_define_type.value=="MWPotential2014") {
            if (type_args.size()==0) type_args="15:0.0299946,1.8,0.2375|5:0.7574802,0.375,0.035|9:4.85223053,2.0";
            else type_args += "|15:0.0299946,1.8,0.2375|5:0.7574802,0.375,0.035|9:4.85223053,2.0";
        }
        std::vector<std::string> type_args_pair;

        std::vector<int> pot_type;
        std::vector<double> pot_args;

        istart = 0;
        inext = type_args.find_first_of("|");
        while (inext!=std::string::npos) {
            type_args_pair.push_back(type_args.substr(istart,inext-istart));
            istart = inext+1;
            inext = type_args.find_first_of("|",istart);
        }
        if (istart!=type_args.size()) type_args_pair.push_back(type_args.substr(istart));

        for (std::size_t i=0; i<type_args_pair.size(); i++) {
            istart = 0;
            inext = type_args_pair[i].find_first_of(":");
            if (inext==std::string::npos) {
                std::cerr<<"Error: potential type index delimiter ':' is not found in input string "<<type_args_pair[i]<<std::endl;
                abort();
            }
            int type_i=std::stoi(type_args_pair[i].substr(0,inext));
            pot_type.push_back(type_i);
            if (_print_flag) std::cout<<"Type index: "<<type_i<<"  args:";
            
            istart = inext+1;
            inext = type_args_pair[i].find_first_of(",",istart);
            while (inext!=std::string::npos) {
                double arg_i = std::stod(type_args_pair[i].substr(istart, inext-istart));
                pot_args.push_back(arg_i);
                if (_print_flag) std::cout<<" "<<arg_i;
                istart = inext+1;
                inext = type_args_pair[i].find_first_of(",",istart);
            }
            if(istart!=type_args.size()) {
                double arg_i = std::stod(type_args_pair[i].substr(istart));
                pot_args.push_back(arg_i);
                if (_print_flag) std::cout<<" "<<arg_i;
            }
            if (_print_flag) std::cout<<std::endl;
        }
        npot = type_args_pair.size();
        
        potential_args = new struct potentialArg[npot];
        int* pot_type_ptr = pot_type.data();
        double* pot_args_ptr = pot_args.data();

        parse_leapFuncArgs_Full(npot, potential_args, &pot_type_ptr, &pot_args_ptr);
    }

    //! calculate acceleration and potential at give position
    /*!
      @param[out] acc: [3] acceleration to return
      @param[out] pot: potential to return 
      @param[in] pos: position of particles
     */
    void calcAccPot(double* acc, double& pot, const double t, const double* pos) {
        if (npot>0) {
            double x = pos[0] + pos_offset[0];
            double y = pos[1] + pos_offset[1];
            double z = pos[2] + pos_offset[2];
        
            double rxy = std::sqrt(x*x+y*y);
            double phi = std::acos(x/rxy);
            double sinphi = y/rxy;
            double cosphi = x/rxy;
            double acc_rxy = calcRforce(rxy, z, phi, t, npot, potential_args);
            double acc_z   = calczforce(rxy, z, phi, t, npot, potential_args);
            double acc_phi = calcPhiforce(rxy, z, phi, t, npot, potential_args);
            acc[0] = cosphi*acc_rxy - sinphi*acc_phi/rxy;
            acc[1] = sinphi*acc_rxy + cosphi*acc_phi/rxy;
            acc[2] = acc_z;
        
            pot = acc[0]*x+acc[1]*y+acc[2]*z;
        }
        else {
            acc[0] = acc[1] = acc[2] = pot = 0.0;
        }
    }

    void clear() {
        if (potential_args!=NULL) {
            free_potentialArgs(npot, potential_args);
            free(potential_args);
        }
    }
};

