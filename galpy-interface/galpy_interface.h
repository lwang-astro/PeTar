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
    IOParams<std::string> config_filename; // configure file name
    IOParams<double> rscale; 
    IOParams<double> tscale; 
    IOParams<double> vscale; 
    IOParams<double> fscale; 
    IOParams<double> pscale; 
    
    bool print_flag;

    IOParamsGalpy(): input_par_store(),
                     type_args (input_par_store, "", "Description of potential type and arguments"),
                     pre_define_type (input_par_store, "", "Add additional Pre-defined potential type to the potential list, options are: MWPotential2014"),
                     pos_offset(input_par_store, "0:0:0", "3D position offset in unit of the input particle data (IN) to obtain the particle position in the potential field, x, y, z are separated by ':'"),
                     config_filename(input_par_store, "", "A configure file that store the types and arguments of potential, will be added to the potential list"),
                     rscale(input_par_store, 1.0, "Radius scale factor from unit of the input particle data (IN) to Galpy distance unit (r[Galpy]=r[IN]*rscale)"),
                     tscale(input_par_store, 1.0, "Time scale factor from unit of the input particle data (IN) to Galpy time (time[Galpy]=time[IN]*tscale)"),
                     vscale(input_par_store, 1.0, "Velocity scale factor from unit of the input particle data (IN) to Galpy velocity unit (v[Galpy]=v[IN]*vscale)"),
                     fscale(input_par_store, 1.0, "Acceleration scale factor from unit of the input particle data (IN) to Galpy acceleration unit (v[Galpy]=v[IN]*vscale)"),
                     pscale(input_par_store, 1.0, "Potential factor from unit of the input particle data (IN) to Galpy potential unit (v[Galpy]=v[IN]*vscale)"),
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
            {"galpy-type-arg",   required_argument, &galpy_flag, 0}, 
            {"galpy-pos-offset", required_argument, &galpy_flag, 1}, 
            {"galpy-set",        required_argument, &galpy_flag, 2}, 
            {"galpy-rscale",     required_argument, &galpy_flag, 3}, 
            {"galpy-tscale",     required_argument, &galpy_flag, 4}, 
            {"galpy-vscale",     required_argument, &galpy_flag, 5}, 
            {"galpy-fscale",     required_argument, &galpy_flag, 6}, 
            {"galpy-pscale",     required_argument, &galpy_flag, 7}, 
            {"galpy-conf-file",  required_argument, &galpy_flag, 8}, 
            {"help", no_argument, 0, 'h'},
            {0,0,0,0}
        };
        
        int opt_used=opt_used_pre;
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
                case 3:
                    rscale.value = atof(optarg);
                    if(print_flag) rscale.print(std::cout);
                    opt_used+=2;
                    break;
                case 4:
                    tscale.value = atof(optarg);
                    if(print_flag) tscale.print(std::cout);
                    opt_used+=2;
                    break;
                case 5:
                    vscale.value = atof(optarg);
                    if(print_flag) vscale.print(std::cout);
                    opt_used+=2;
                    break;
                case 6:
                    fscale.value = atof(optarg);
                    if(print_flag) fscale.print(std::cout);
                    opt_used+=2;
                    break;
                case 7:
                    pscale.value = atof(optarg);
                    if(print_flag) pscale.print(std::cout);
                    opt_used+=2;
                    break;
                case 8:
                    config_filename.value = optarg;
                    if(print_flag) config_filename.print(std::cout);
                    opt_used+=2;
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
                    std::cout<<"         --galpy-type-arg:  [S] "<<type_args<<std::endl
                             <<"             The format of a combination of types and arguments have two styles: (no space in the middle):\n"
                             <<"             1) [type1]:[arg1-1],[arg1-2],**|[type2]:[arg2-1],[arg2-2],**\n"
                             <<"             2) [type1],[type2],..:[arg1-1],[arg1-2],[arg2-1],**\n"
                             <<"             where '|' split different potential types;\n"
                             <<"                   ':' separates the type indices and the argument list\n"
                             <<"                   ',' separates the items of types or arguments in their lists\n"
                             <<"             For example:\n"
                             <<"                1) --galpy-type-arg 15:0.0299946,1.8,0.2375|5:0.7574802,0.375,0.035|9:4.85223053,2.0\n"
                             <<"                2) --galpy-type-arg 15,5,9:0.0299946,1.8,0.2375,0.7574802,0.375,0.035,4.85223053,2.0\n"
                             <<"                both can generate the MWPotential2014 from Galpy (same to --galpy-set MWPotential2014)\n"
                             <<"             The defaults values of types and arguments for supported potentials can be found by using the commander, petar.galpy.help.\n"
                             <<"         --galpy-pos-offset:[S] "<<pos_offset<<std::endl
                             <<"         --galpy-set:       [S] "<<pre_define_type<<std::endl
                             <<"         --galpy-conf-file  [S] "<<config_filename<<std::endl
                             <<"         --galpy-rscale:    [F] "<<rscale<<std::endl
                             <<"         --galpy-tscale:    [F] "<<tscale<<std::endl
                             <<"         --galpy-vscale:    [F] "<<vscale<<std::endl
                             <<"         --galpy-fscale:    [F] "<<fscale<<std::endl
                             <<"         --galpy-pscale:    [F] "<<pscale<<std::endl
                             <<"              The default arguments in potentials are using the Galpy (natural) unit (velocity in [220 km/s], distance in [8 kpc]).\n"
                             <<"              If the input particle data have different units, the scaling factor should be properly set."<<std::endl;

                }
                return -1;
            case '?':
                opt_used +=2;
                break;
            default:
                break;
            }
        
        if(print_flag) std::cout<<"----- Finish reading input options of Galpy -----\n";

        return opt_used;
    }    

    //! set nature unit scaling to match Galpy unit and astronomical unit set (Myr, Msun, pc)
    void setStdUnit(const bool print_flag=true) {
        double kms_pcmyr=1.0227121655399913;
        double vbase=220.0;
        double rbase=8.0;
        double vb_pcmyr = vbase*kms_pcmyr;
        rscale.value = 0.001/rbase; // pc to solar position in Milkyway
        vscale.value = 1.0/vb_pcmyr; // pc/Myr to solar velocity in Milkyway
        tscale.value = rscale.value/vscale.value; // myr to time unit in galpy
        fscale.value = vscale.value*vscale.value/rscale.value; // pc/myr^2 to acc unit in galpy
        pscale.value = vscale.value*vscale.value; // pc^2/myr^2 to pot unit in galpy

        if (print_flag) {
            std::cout<<"Unit conversion for Galpy:\n"
                     <<"rscale = "<<rscale.value<<" pc -> solar position (8 kpc) in Milkyway\n"
                     <<"vscale = "<<vscale.value<<" pc/myr -> solar velocity (220 km/s) in Milkyway\n"
                     <<"tscale = "<<tscale.value<<" myr -> solar orbital period (225 Myr)/(2 pi) in Milkyway\n"
                     <<"fscale = "<<fscale.value<<" pc/myr^2 -> galpy acceleration unit\n"
                     <<"pscale = "<<pscale.value<<" pc^2/myr^2 -> galpy potential unit\n";
        }
    }
};

//! A class to manager the API to Galpy
class GalpyManager{
public:
    struct potentialArg* potential_args;
    int npot; // number of potential models 
    double pos_offset[3]; // position offset
    double rscale;
    double tscale;
    double vscale;
    double fscale;
    double pscale;

    GalpyManager(): potential_args(NULL), npot(0), pos_offset{0.0,0.0,0.0}, rscale(1.0), tscale(1.0), vscale(1.0), fscale(1.0), pscale(1.0) {}

    //! initialization function
    void initial(const IOParamsGalpy& _input, const bool _print_flag=false) {
        // unit scale
        rscale = _input.rscale.value;
        vscale = _input.vscale.value;
        tscale = _input.tscale.value;
        fscale = _input.fscale.value;
        pscale = _input.pscale.value;

        // position offset
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
        if (_print_flag) std::cout<<"Central position offset: "<<pos_offset[0]<<" "<<pos_offset[1]<<" "<<pos_offset[2]<<std::endl;

        // add pre-defined type-argu groups
        std::string type_args = _input.type_args.value;
        if (_input.pre_define_type.value=="MWPotential2014") {
            if (type_args.size()==0) type_args="15:0.0299946,1.8,0.2375|5:0.7574802,0.375,0.035|9:4.85223053,2.0";
            else type_args += "|15:0.0299946,1.8,0.2375|5:0.7574802,0.375,0.035|9:4.85223053,2.0";
        }

        std::vector<std::string> type_args_pair;
        std::vector<int> pot_type;
        std::vector<double> pot_args;

        // split type-arg groups
        istart = 0;
        inext = type_args.find_first_of("|");
        while (inext!=std::string::npos) {
            type_args_pair.push_back(type_args.substr(istart,inext-istart));
            istart = inext+1;
            inext = type_args.find_first_of("|",istart);
        }
        if (istart!=type_args.size()) type_args_pair.push_back(type_args.substr(istart));

        if (_print_flag) std::cout<<"Potential combination list:\n";
        // loop group
        for (std::size_t i=0; i<type_args_pair.size(); i++) {
            // get types
            istart = 0;
            inext = type_args_pair[i].find_first_of(":");
            if (inext==std::string::npos) {
                std::cerr<<"Error: potential type index delimiter ':' is not found in input string "<<type_args_pair[i]<<std::endl;
                abort();
            }
            std::string type_str = type_args_pair[i].substr(0,inext);
            std::string arg_str = type_args_pair[i].substr(inext+1);
            
            inext = type_str.find_first_of(",");
            if (_print_flag) std::cout<<"Type index: ";
            while (inext!=std::string::npos) {
                int type_i=std::stoi(type_str.substr(istart, inext-istart));
                pot_type.push_back(type_i);
                if (_print_flag) std::cout<<type_i<<" ";
                istart = inext+1;
                inext = type_str.find_first_of(",",istart);
            }
            if (istart!=type_str.size()) {
                int type_i=std::stoi(type_str.substr(istart));
                pot_type.push_back(type_i);
                if (_print_flag) std::cout<<type_i<<" ";
            }
            if (_print_flag) std::cout<<"args:";
            
            // get arguments
            istart = 0;
            inext = arg_str.find_first_of(",",istart);
            while (inext!=std::string::npos) {
                double arg_i = std::stod(arg_str.substr(istart, inext-istart));
                pot_args.push_back(arg_i);
                if (_print_flag) std::cout<<" "<<arg_i;
                istart = inext+1;
                inext = arg_str.find_first_of(",",istart);
            }
            if(istart!=arg_str.size()) {
                double arg_i = std::stod(arg_str.substr(istart));
                pot_args.push_back(arg_i);
                if (_print_flag) std::cout<<" "<<arg_i;
            }
            if (_print_flag) std::cout<<std::endl;
        }
        npot = pot_type.size();

        // add type arguments from configure file if exist
        if (_input.config_filename.value!="") {
            std::ifstream fin;
            fin.open(_input.config_filename.value.c_str(), std::ifstream::in);
            int nadd=0;
            fin>>nadd;
            assert(nadd>0);
            npot += nadd;
            if (_print_flag) std::cout<<"Type index: ";
            for (int i=0; i<nadd; i++) {
                int type_i;
                fin>>type_i;
                if (fin.eof()) {
                    std::cerr<<"Read type fails! required number is "<<nadd<<" only get "<<i;
                    abort();
                }
                pot_type.push_back(type_i);
                if (_print_flag) std::cout<<type_i<<" ";
            }
            if (_print_flag) std::cout<<"args:";
            while (true) {
                double arg_i;
                fin>>arg_i;
                if (fin.eof()) break;
                pot_args.push_back(arg_i);
                if (_print_flag) std::cout<<" "<<arg_i;
            }
            if (_print_flag) std::cout<<std::endl;
        }
        
        // generate galpy potential argument 
        potential_args = new struct potentialArg[npot];
        int* pot_type_ptr = pot_type.data();
        double* pot_args_ptr = pot_args.data();

        parse_leapFuncArgs_Full(npot, potential_args, &pot_type_ptr, &pot_args_ptr);

        if(_print_flag) std::cout<<"----- Finish initialize Galpy potential -----\n";

    }

    //! calculate acceleration and potential at give position
    /*!
      @param[out] acc: [3] acceleration to return
      @param[out] pot: potential to return 
      @param[in] time: time in input unit
      @param[in] pos: position of particles in input unit
     */
    void calcAccPot(double* acc, double& pot, const double time, const double* pos) {
        if (npot>0) {
            
            double x = (pos[0] + pos_offset[0])*rscale;
            double y = (pos[1] + pos_offset[1])*rscale;
            double z = (pos[2] + pos_offset[2])*rscale;
            double t = time*tscale;
        
            double rxy = std::sqrt(x*x+y*y);
            double phi = std::acos(x/rxy);
            double sinphi = y/rxy;
            double cosphi = x/rxy;
            double acc_rxy = calcRforce(rxy, z, phi, t, npot, potential_args);
            double acc_z   = calczforce(rxy, z, phi, t, npot, potential_args);
            double acc_phi = calcPhiforce(rxy, z, phi, t, npot, potential_args);
            pot = evaluatePotentials(rxy, z, npot, potential_args)/pscale;
            acc[0] = (cosphi*acc_rxy - sinphi*acc_phi/rxy)/fscale;
            acc[1] = (sinphi*acc_rxy + cosphi*acc_phi/rxy)/fscale;
            acc[2] = acc_z/fscale;
        
            //pot = acc[0]*x+acc[1]*y+acc[2]*z;
        }
        else {
            acc[0] = acc[1] = acc[2] = pot = 0.0;
        }
    }

    void clear() {
        if (potential_args!=NULL) {
            free_potentialArgs(npot, potential_args);
            free(potential_args);
            npot = 0;
        }
        pos_offset[0] = pos_offset[1] = pos_offset[2] = 0.0;
    }

    ~GalpyManager() {
        clear();
    }
};

