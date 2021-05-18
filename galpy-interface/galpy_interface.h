#include <iostream>
#include <iomanip>
#include <cassert>
#include <string>
#include <vector>
#include <getopt.h>
#include <cmath>
#include "../src/io.hpp"

#include <integrateFullOrbit.h>

//! IO parameters for Galpy manager
class IOParamsGalpy{
public:
    IOParamsContainer input_par_store;
    IOParams<std::string> type_args; // potential type and argument list
    IOParams<std::string> pre_define_type; // pre defined type name
    IOParams<std::string> config_filename; // configure file name
    IOParams<double> rscale; 
    IOParams<double> tscale; 
    IOParams<double> vscale; 
    IOParams<double> fscale; 
    IOParams<double> pscale; 
    
    bool print_flag;

    IOParamsGalpy(): input_par_store(),
                     type_args (input_par_store, "__NONE__", "galpy-type-arg", "Add potential types and argumentsto the potential list in the center of the galactic reference frame"),
                     pre_define_type (input_par_store, "__NONE__", "galpy-set", "Add a pre-defined potential set to the potential list, options are: MWPotential2014"),
                     config_filename(input_par_store, "__NONE__", "galpy-conf-file", "A configure file of the time- and space-dependent potentials"),
                     rscale(input_par_store, 1.0, "galpy-rscale", "Radius scale factor from unit of the input particle data (IN) to Galpy distance unit (r[8kpc]=r[IN]*rscale)"),
                     tscale(input_par_store, 1.0, "galpy-tscale", "Time scale factor (rscale/vscale) from unit of the input particle data (IN) to Galpy time (time[Solar orbital period/2 pi]=time[IN]*tscale)"),
                     vscale(input_par_store, 1.0, "galpy-vscale", "Velocity scale factor from unit of the input particle data (IN) to Galpy velocity unit (v[220 km/s]=v[IN]*vscale)"),
                     fscale(input_par_store, 1.0, "galpy-fscale", "Acceleration scale factor (vscale^2/rscale) from unit of the input particle data (IN) to Galpy acceleration unit (acc[Galpy]=acc[IN]*fscale)"),
                     pscale(input_par_store, 1.0, "galpy-pscale", "Potential scale factor (vscale^2) from unit of the input particle data (IN) to Galpy potential unit (pot[Galpy]=pot[IN]*pscale)"),
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
            {type_args.key,       required_argument, &galpy_flag, 0}, 
            {config_filename.key, required_argument, &galpy_flag, 1}, 
            {pre_define_type.key, required_argument, &galpy_flag, 2}, 
            {rscale.key,     required_argument, &galpy_flag, 3}, 
            {tscale.key,     required_argument, &galpy_flag, 4}, 
            {vscale.key,     required_argument, &galpy_flag, 5}, 
            {fscale.key,     required_argument, &galpy_flag, 6}, 
            {pscale.key,     required_argument, &galpy_flag, 7}, 
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
                    config_filename.value = optarg;
                    if(print_flag) config_filename.print(std::cout);
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
                    input_par_store.printHelp(std::cout, 2, 10, 23);
                    std::cout<<"*** PS: for --galpy-type-arg: the parameters to set up a set of potentials fixed at the center of the galactic reference frame.\n"
                             <<"        This is for a quick set-up of potentials. For a more complex and flexible set of potentials, use --galpy-conf-file.\n"
                             <<"        The format for a combination of the types and arguments of potentials have two styles: (no space in the middle):\n"
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
                             <<"             The default arguments in potentials are using the Galpy (natural) unit (velocity in [220 km/s], distance in [8 kpc]).\n"
                             <<"             If the input particle data have different units, the scaling factor should be properly set."
                             <<"             The potentials defined in --galpy-type-arg and --galpy-set are both added.\n"
                             <<"             Thus, don't repeat the same potential sets in the two options.\n"
                             <<"       for --galpy-conf-file: the configure file containing the parameters for time- and space-dependent potentials.\n"
                             <<"             Users can add arbitrary number of potentials in two levels depending on time and coordinates, respectively;\n"
                             <<"             Users can use Galpy (_parse_pot function, see its online manual) to find and configure the types and arguments of potentials.\n"
                             <<"             The format of the configure file is like:\n"
                             <<"                 Update_time[Galpy unit] N_sets\n"
                             <<"                 [set 1]\n"
                             <<"                 [set 2]\n"
                             <<"                 ...\n"
                             <<"             This pattern can repeat with an increasing order of Update_time.\n"
                             <<"             PeTar will read the configure file during the running time once the Update_time is reached.\n"
                             <<"             Thus, don't rename, delete or modify the configure file before the simulation ends.\n"
                             <<"             The time interval of update_time should be integer times of the tree time step used in PeTar.\n"
                             <<"             The last group will be used until the end of the simulation.\n"
                             <<"\n"
                             <<"             N_sets indicates the number of potential sets, each set shares the same central position and velocity.\n"
                             <<"             If N_sets = 0, the external potential is switched off.\n"
                             <<"             Each set contains three lines: \n"
                             <<"                 1st: N_types Mode x y z vx vy vz\n"
                             <<"                 2nd: type1 type2 ...\n"
                             <<"                 3rd: arg1-1 arg1-2 ... arg2-1 ...\n"
                             <<"             The definitions of each items:\n"
                             <<"                 N_types: number of potential types\n"
                             <<"                       If there is no argument, the 3rd line is still necessary (empty line).\n"
                             <<"                 Mode: an index to indicate the reference frame of the following position and velocity with two choices: 0 and 1\n"
                             <<"                       0: The galactic frame\n"
                             <<"                       1: The frame of the particle system (following the motion of the particle system)\n"
                             <<"                 x,y,z,vx,vy,vz: the position and velocity of the potential center [Galpy unit] in the reference frame defined in the Mode.\n"
                             <<"             2nd and 3rd lines are types and arguments of each potential, similar to those in --galpy-type-arg.\n"
                             <<"             The number of items in the 2nd line must be the same as Number_of_types.\n"
                             <<"\n"
                             <<"             For example, to add the MWPotential2014 fixed at the galactic center and a Plummer potential following the motion of the particle system,\n"
                             <<"             and then, remove the Plummer potential after 1.0 Galpy time unit, the configure file looks like:\n"
                             <<"                 0.0 2                           #[Update_time N_sets]\n"
                             <<"                 3 0 0.0 0.0 0.0 0.0 0.0 0.0     #[N_types, Mode, x, y, z, vx, vy, vz]\n"
                             <<"                 15 5 9                          #[type1, type2, type3 (MWPotential2014)]\n"
                             <<"                 0.029994597188218 1.8 0.2375 0.75748020193716 0.375 0.035 4.852230533528 2.0 #[arg1-1, arg1-2, ...]\n"
                             <<"                 1 1 0.0 0.0 0.0 0.0 0.0 0.0     #[N_types, Mode, x, y, z, vx, vy, vz]\n"
                             <<"                 17                              #[type1 (Plummer)]\n"
                             <<"                 1.11072675e-8 0.000125          #[arg1-1, arg1-2]\n"
                             <<"                 1.0 1                           #[Update_time N_sets]\n"
                             <<"                 3 G 0.0 0.0 0.0 0.0 0.0 0.0     #[N_types, Mode, x, y, z, vx, vy, vz]\n"
                             <<"                 15 5 9                          #[type1, type2, type3 (MWPotential2014)]\n"
                             <<"                 0.029994597188218 1.8 0.2375 0.75748020193716 0.375 0.035 4.852230533528 2.0 #[arg1-1, arg1-2, ...]\n"
                             <<"             Here the G*M and distance scaling factors are 2.4692087520131e-09 [galpy GM unit] / [pc^3/Myr^2] and 0.000125 [8 kpc] / [pc], respectively;\n"
                             <<"             The plummer sphere has a total mass of 1000 Msun and a scale radius of 1 pc at time zero.\n"
                             <<"             Notice that the comments after the symbol # is for the reference here, they cannot appear in the configure file.\n"
                             <<"       Users can either use --galpy-type-arg and --galpy-set or --galpy-conf-file. But if both are used, the error will appear."
                             <<std::endl;
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
        double kms_pcmyr=1.022712165045695; // pc = 30856775814913673 m, yr = 365.25 days
        double vbase=220.0;
        double rbase=8.0;
        double vb_pcmyr = vbase*kms_pcmyr;
        rscale.value = 0.001/rbase; // pc to solar position in Milkyway
        vscale.value = 1.0/vb_pcmyr; // pc/Myr to solar velocity in Milkyway
        tscale.value = rscale.value/vscale.value; // myr to time unit in galpy
        fscale.value = vscale.value*vscale.value/rscale.value; // pc/Myr^2 to acc unit in galpy
        pscale.value = vscale.value*vscale.value; // pc^2/Myr^2 to pot unit in galpy

        if (print_flag) {
            double GMscale = pscale.value*rscale.value; // pc^3/Myr^2
            std::cout<<"----- Unit conversion for Galpy -----\n"
                     <<"rscale  = "<<rscale.value<<"  [8 kpc] / [pc]\n"
                     <<"vscale  = "<<vscale.value<<"  [220 km/s] / [pc/Myr]\n"
                     <<"tscale  = "<<tscale.value<<"  [Solar orbital period/2 pi] / Myr\n"
                     <<"fscale  = "<<fscale.value<<"  [galpy acceleration unit] / [pc/Myr^2]\n"
                     <<"pscale  = "<<pscale.value<<"  [galpy potential unit] / [pc^2/Myr^2]\n"
                     <<"GMscale = "<<GMscale     <<"  [galpy gravitational constant * mass unit] / [pc^3/Myr^2]\n";
        }
    }
};

//! set of Galpy potentials sharing the same position and velocity of the origin
struct PotentialSet{
    int mode; // mode of origin: 0: galactic frame; 1: local particle system frame
    double pos[3]; // position
    double vel[3]; // velocity
    int npot; // number of potential models 
    struct potentialArg* arguments; //potential arguments array for Gaply

    PotentialSet(): mode(-1), pos{0.0,0.0,0.0}, vel{0.0,0.0,0.0}, npot(0), arguments(NULL) {}

    //! set position and velocity
    void setOrigin(const int _mode, const double* _pos=NULL, const double* _vel=NULL) {
        assert(_mode==0||_mode==1);
        mode = _mode;
        if (_pos!=NULL) {
            pos[0] = _pos[0];
            pos[1] = _pos[1];
            pos[2] = _pos[2];
        }
        if (_vel!=NULL) {
            vel[0] = _vel[0];
            vel[1] = _vel[1];
            vel[2] = _vel[2];
        }
    }

    //! generate Galpy potential arguments 
    /*! 
      @param[in] _npot: number of potential types
      @param[in] _type: the array of type indices
      @param[in] _arg: the array of arguments for each type
     */ 
    void generatePotentialArgs(const int _npot, int* _type, double* _arg) {
        assert(_npot>0);
        npot = _npot;
        arguments = new struct potentialArg[npot]; 
        parse_leapFuncArgs_Full(npot, arguments, &_type, &_arg);
    }

    void clear() {
        if (arguments!=NULL) {
            assert(npot>0);
            free_potentialArgs(npot, arguments);
            free(arguments);
        }
        mode = -1;
        pos[0] = pos[1] = pos[2] = 0.0;
        vel[0] = vel[1] = vel[2] = 0.0;
        npot = 0;
        arguments = NULL;
    }
};

//! A class to manager the API to Galpy
class GalpyManager{
public:
    std::vector<PotentialSet> potential_sets;
    double update_time;
    double rscale;
    double tscale;
    double vscale;
    double fscale;
    double pscale;
    std::ifstream fconf;

    GalpyManager(): potential_sets(), update_time(0.0), rscale(1.0), tscale(1.0), vscale(1.0), fscale(1.0), pscale(1.0), fconf() {}

    //! initialization function
    /*!
      @param[in] _input: input parameters 
      @param[in] _print_flag: if true, printing information to std::cout
     */
    void initial(const IOParamsGalpy& _input, const bool _print_flag=false) {
        // unit scale
        rscale = _input.rscale.value;
        vscale = _input.vscale.value;
        tscale = _input.tscale.value;
        fscale = _input.fscale.value;
        pscale = _input.pscale.value;

        // add pre-defined type-argu groups
        std::string type_args = _input.type_args.value;
        if (_input.pre_define_type.value=="MWPotential2014") {
            if (type_args=="__NONE__") type_args="15:0.0299946,1.8,0.2375|5:0.7574802,0.375,0.035|9:4.85223053,2.0";
            else type_args += "|15:0.0299946,1.8,0.2375|5:0.7574802,0.375,0.035|9:4.85223053,2.0";
        }
        if (type_args!="__NONE__"&&_input.config_filename.value!="__NONE__")  {
            std::cerr<<"Galpy Error: both --galpy-type-arg|--galpy-set and --galpy-conf-file are used, please choose one of them."<<std::endl;
            abort();
        }
        if (type_args!="__NONE__") {

            std::vector<std::string> type_args_pair;
            std::vector<int> pot_type;
            std::vector<double> pot_args;

            // split type-arg groups
            std::size_t istart = 0;
            std::size_t inext = type_args.find_first_of("|");
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
            int npot = pot_type.size();

            // generate galpy potential argument 
            potential_sets.push_back(PotentialSet());
            auto& pset = potential_sets.back();
            pset.setOrigin(0);
            pset.generatePotentialArgs(npot, pot_type.data(), pot_args.data());
        }

        // add type arguments from configure file if exist
        if (_input.config_filename.value!="__NONE__") {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
            int my_rank = PS::Comm::getRank();
            if (my_rank==0) {
#endif
                fconf.open(_input.config_filename.value.c_str(), std::ifstream::in);
                fconf>>update_time;
                if(fconf.eof()) fconf.close();
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
            }
            PS::Comm::broadcast(&update_time, 1, 0);
#endif
        }

        if(_print_flag) std::cout<<"----- Finish initialize Galpy potential -----\n";

    }

    //! Update types and arguments from configure file if update_time <= particle system time
    /*! 
      @param[in] _system_time: the time of globular particley system in PeTar. The time is transferred based on the time scaling factor
      @param[in] _print_flag: if true, print the read types and arguments.
     */
    void updateTypesAndArgsFromFile(const double& _system_time, const bool _print_flag) {
        double time_scaled= _system_time*tscale;
        int nset=-1;
        std::vector<int> mode;
        std::vector<double> origin;
        std::vector<int> pot_type_offset;
        std::vector<int> pot_type;
        std::vector<int> pot_args_offset;
        std::vector<double> pot_args;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
        int my_rank = PS::Comm::getRank();
        if (my_rank==0) {
#endif
            if (fconf.is_open())  {
                if (time_scaled>=update_time) {

                    double update_time_next;

                    while(true) {
                        fconf>>nset;
                        if(fconf.eof()) {
                            std::cerr<<"Galpy config: reading number of sets fails! File reaches EOF.";
                            abort();
                        }
                        // clear
                        mode.resize(0);
                        origin.resize(0);
                        pot_type_offset.resize(0);
                        pot_type.resize(0);
                        pot_args_offset.resize(0);
                        pot_args.resize(0);

                        if (nset>0) {
                            pot_type_offset.push_back(0);
                            pot_args_offset.push_back(0);
                            for (int k=0; k<nset; k++) {
                                int n_pot_k=-1;
                                fconf>>n_pot_k;
                                if(fconf.eof()) {
                                    std::cerr<<"Galpy config: reading number of types for Set "<<k+1<<" fails! File reaches EOF.";
                                    abort();
                                }
                                assert(n_pot_k>0);
                                // mode
                                int mode_k;
                                fconf>>mode_k;
                                if(fconf.eof()) {
                                    std::cerr<<"Galpy config: reading mode for Set "<<k+1<<" fails! File reaches EOF.";
                                    abort();
                                }
                                mode.push_back(mode_k);
                                // origin
                                double posvel[6];
                                fconf>>posvel[0]>>posvel[1]>>posvel[2]>>posvel[3]>>posvel[4]>>posvel[5];
                                if(fconf.eof()) {
                                    std::cerr<<"Galpy config: reading position and velocity of origin for Set "<<k+1<<" fails! File reaches EOF.";
                                    abort();
                                }
                                for (int i=0; i<6; i++) origin.push_back(posvel[i]);
                                // types
                                for (int i=0; i<n_pot_k; i++) {
                                    int type_i;
                                    fconf>>type_i;
                                    if (fconf.eof()) {
                                        std::cerr<<"Galpy config: reading potential types for Set "<<k+1<<" fails! required number is "<<n_pot_k<<" only get "<<i;
                                        abort();
                                    }
                                    pot_type.push_back(type_i);
                                }
                                pot_type_offset.push_back(pot_type.size());

                                // arguments
                                std::string line;
                                std::getline(fconf, line); // skip 2nd line
                                std::getline(fconf, line);
                                std::istringstream fin(line);
                                double arg_i;
                                while (fin>>arg_i) {
                                    pot_args.push_back(arg_i);
                                }
                                pot_args_offset.push_back(pot_args.size());
                            }
                        }

                        fconf>>update_time_next;
                        if(fconf.eof()) {
                            fconf.close();
                            break;
                        }
                        if (update_time_next>time_scaled) {
                            update_time = update_time_next;
                            break;
                        }
                    };

                    if (_print_flag) {
                        std::cout<<"Galpy time: "<<time_scaled;
                        if (fconf.is_open()) std::cout<<" Next update time: "<<update_time;
                        else std::cout<<" Next update time: End_of_simulation";
                        std::cout<<std::endl;
                        for (int k=0; k<nset; k++) {
                            std::cout<<"Set "<<k+1<<" Mode: "<<mode[k]
                                     <<" Pos: "<<origin[6*k]<<" "<<origin[6*k+1]<<" "<<origin[6*k+2]
                                     <<" Vel: "<<origin[6*k+3]<<" "<<origin[6*k+4]<<" "<<origin[6*k+5]
                                     <<" Type index: ";
                            for (int i=pot_type_offset[k]; i<pot_type_offset[k+1]; i++) std::cout<<pot_type[i]<<" ";
                            std::cout<<" Arguments: ";
                            for (int i=pot_args_offset[k]; i<pot_args_offset[k+1]; i++) std::cout<<pot_args[i]<<" ";
                            std::cout<<std::endl;
                        }
                    }
                }
            }
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
        }
        PS::Comm::broadcast(&nset, 1, 0);
#endif
        // update potentials
        if (nset>=0) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            PS::Comm::broadcast(&update_time, 1, 0);
            if (my_rank==0) {
                assert((int)mode.size()==nset);
                assert((int)origin.size()==6*nset);
                assert((int)pot_type_offset.size()==nset+1);
                assert((int)pot_args_offset.size()==nset+1);
            }
            else {
                mode.resize(nset);
                origin.resize(6*nset);
                pot_type_offset.resize(nset+1);
                pot_args_offset.resize(nset+1);
            }
            PS::Comm::broadcast(mode.data(), mode.size(), 0);
            PS::Comm::broadcast(origin.data(), origin.size(), 0);
            PS::Comm::broadcast(pot_type_offset.data(), pot_type_offset.size(), 0);
            PS::Comm::broadcast(pot_args_offset.data(), pot_args_offset.size(), 0);
            
            if (my_rank==0) {
                assert((int)pot_type.size()==pot_type_offset.back());
                assert((int)pot_args.size()==pot_args_offset.back());
            }
            else {
                pot_type.resize(pot_type_offset.back());
                pot_args.resize(pot_args_offset.back());
            }
            PS::Comm::broadcast(pot_type.data(), pot_type.size(), 0);
            if (pot_args.size()>0) PS::Comm::broadcast(pot_args.data(), pot_args.size(), 0);
#endif
            // generate galpy potential argument 
            freePotentialArgs();
            for (int k=0; k<nset; k++) {
                potential_sets.push_back(PotentialSet());
                auto& pset = potential_sets.back();
                pset.setOrigin(mode[k],&origin[6*k],&origin[6*k+3]);
                pset.generatePotentialArgs(pot_type_offset[k+1]-pot_type_offset[k], &(pot_type[pot_type_offset[k]]), &(pot_args[pot_args_offset[k]]));
            }
        }
    }

    //! print reference to cite
    static void printReference(std::ostream & fout, const int offset=4) {
        for (int i=0; i<offset; i++) fout<<" ";
        fout<<"Galpy: Bovy J., 2015, ApJS, 216, 29"<<std::endl;
    }    

    //! calculate acceleration and potential at give position
    /*!
      @param[out] acc: [3] acceleration to return
      @param[out] pot: potential to return 
      @param[in] time: time in input unit
      @param[in] pos_g: position of particles in the galactic frame [input unit]
      @param[in] pos_l: position of particles in the particle system frame [input unit]
     */
    void calcAccPot(double* acc, double& pot, const double time, const double* pos_g, const double* pos_l) {
        if (potential_sets.size()>0) {
            
            double t = time*tscale;

            double x[2] = {pos_g[0]*rscale, pos_l[0]*rscale};
            double y[2] = {pos_g[1]*rscale, pos_l[1]*rscale};
            double z[2] = {pos_g[2]*rscale, pos_l[2]*rscale};

            double rxy[2], phi[2], sinphi[2], cosphi[2];
            for (int i=0; i<2; i++) {
                rxy[i]= std::sqrt(x[i]*x[i]+y[i]*y[i]);
                phi[i]= std::acos(x[i]/rxy[i]);
                sinphi[i] = y[i]/rxy[i];
                cosphi[i] = x[i]/rxy[i];
            }

            pot = 0;
            acc[0] = acc[1] = acc[2] = 0.0;

            for (size_t k=0; k<potential_sets.size(); k++) {
                int i = potential_sets[k].mode;
                assert(i==0||i==1);
                int npot = potential_sets[k].npot;
                auto& potential_args = potential_sets[k].arguments;
                double acc_rxy = calcRforce(rxy[i], z[i], phi[i], t, npot, potential_args);
                double acc_z   = calczforce(rxy[i], z[i], phi[i], t, npot, potential_args);
                double acc_phi = calcPhiforce(rxy[i], z[i], phi[i], t, npot, potential_args);
                pot += evaluatePotentials(rxy[i], z[i], npot, potential_args);
                if (rxy[i]>0.0) {
                    acc[0] += (cosphi[i]*acc_rxy - sinphi[i]*acc_phi/rxy[i]);
                    acc[1] += (sinphi[i]*acc_rxy + cosphi[i]*acc_phi/rxy[i]);
                    acc[2] += acc_z;
                }
            }
            pot /= pscale;
            acc[0] /= fscale;
            acc[1] /= fscale;
            acc[2] /= fscale;

            //pot = acc[0]*x+acc[1]*y+acc[2]*z;
        }
        else {
            acc[0] = acc[1] = acc[2] = pot = 0.0;
        }
    }

    void freePotentialArgs() {
        if (!potential_sets.empty()) {
            for (size_t i=0; i<potential_sets.size(); i++) potential_sets[i].clear();
            potential_sets.resize(0);
        }
    }

    void clear() {
        freePotentialArgs();
        update_time = 0;
        if (fconf.is_open()) fconf.close();
    }

    ~GalpyManager() {
        clear();
    }
};


