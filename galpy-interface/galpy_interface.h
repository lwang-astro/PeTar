#include <iostream>
#include <sstream>
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
                             <<"             If the input particle data have different units, the scaling factor should be properly set.\n"
                             <<"             For example, when the particle mass is in Msun, then the argument 'amp' in a Galpy potential can be calculated by G*M*GMscale\n"
                             <<"             where M is in Msun; G = 0.0044983099795944 [Msun, pc, Myr]; and GMscale = 2.4692087520131e-09  [galpy GM unit] / [pc^3/Myr^2].\n"
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
#ifdef GALPY_VERSION_1_7_1
        parse_leapFuncArgs_Full(npot, arguments, &_type, &_arg);
#else
        parse_leapFuncArgs_Full(npot, arguments, &_type, &_arg, NULL);
#endif
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

//! FRW cosmological model for calculating scale factor (redshift)
class FRWModel{
public:
    double time; // current time (Myr)
    double dt_max; // maximum time step to integrate a (Myr)
    double a; //scale factor
    double H0; // Hubble constant (km s^-1 Mpc^-1)
    double H0_myr; // scaled Hubble constant (Myr^-1)
    double omega_energy; // normalized dark energy density
    double omega_radiation; // normalized radiation density
    double omega_matter; // normalized matter density
    double dt_scale; // if -1, evolve a backwards with -dt, if 1, evolve forwards with dt, but time integration is not multiplied by dt_scale

    inline double getMatterDensity(const double _a) const {
        return omega_matter/_a;
    }

    inline double getRadiationDensity(const double _a) const {
        return omega_radiation/(_a*_a);
    }

    inline double getEnergyDensity(const double _a) const {
        return omega_energy*(_a*_a);
    }

    inline double getCurvatureDensity() const {
        return 1 - omega_energy - omega_radiation - omega_matter;
    }

    inline double getDensity(const double _a) const {
        return getMatterDensity(_a) + getRadiationDensity(_a) + getEnergyDensity(_a) + getCurvatureDensity();
    }

    inline void calcH0Myr() {
        double kms_pcmyr=1.022712165045695; // pc = 30856775814913673 m, yr = 365.25 days
        H0_myr = H0*kms_pcmyr*1e-6;
    }

    //! calculate da/dt
    /*!
      @param[in] _a: scale factor
      \return da/dt
     */
    inline double calcAdot(const double _a) const {
        return H0_myr*sqrt(getDensity(_a));
    }

    //! update scale factor by time step _dt (first order Euler method)
    /*! Iterate 10 times to get high accuracy for small a
      @param[in] _dt: time step
     */
    void updateA(const double _dt) {
        double adot = calcAdot(a);
        double a_tmp = a + adot*dt_scale*_dt;
        
        for (int k=0; k<10; k++) {
            double adot_new = calcAdot(a_tmp);
            a_tmp = a + 0.5*(adot_new+adot)*dt_scale*_dt;
        }
        a = a_tmp;
        time += _dt;
    }

    //! update scale factor with maximum time step dt_max (first order Euler method)
    /*! Iterate 10 times to get high accuracy for small a
      @param[in] _dt: time step
     */
    void updateAwithDtmin(const double _dt) {
        assert(_dt>=0);
        assert(dt_max>0);

        if (_dt==0) return;

        double dt = _dt;
        int nstep = 1;
        if (dt>dt_max) {
            dt = dt_max;
            nstep = int(_dt/dt_max);
        }
        for (int i=0; i<nstep; i++) updateA(dt);
        if (dt*nstep<_dt) {
            dt = _dt - dt*nstep;
            updateA(dt);
        }
    }

    //! get Hubble constant at current a
    double getH() const {
        return calcAdot(a)/a;
    }

};

//! data for MWPotential evolve
/*! The potential set refers to MWPotential2014 in Galpy
 */
class MWPotentialEvolve{
public:
    FRWModel frw;
    struct {
        double m_vir_halo;
        double c_halo; // concentration
        double ac_halo; // formation epoch of halo (Wechsler 2002, Zhao 2003)
        double m_disk;
        double ra_disk;
        double rb_disk;
        double m_bulge;
        double alpha_bulge;
        double rcut_bulge;
    } init;
    struct {
        double m_vir_halo;
        double r_vir_halo;
        double rho_c_halo;
        double gm_halo;
        double rs_halo;
        double gm_disk;
        double ra_disk;
        double rb_disk;
        double gm_bulge;
        double alpha_bulge;
        double rcut_bulge;
    } now;

    //! initial potential parameter from a file
    /*!
      @param [in] _filename: the file contain parameters: 
                 FRWModel: Time[Myr], a, H0[km/s/Mpc], Omega_energy, Omega_radiation, Omega_matter, dt_scale
                 NFWHalo: M_vir[Msun], c, ac
                 Disk: M[Msun], Ra[pc], Rb[pc]
                 Bulge: M[Msun], alpha, rcut[pc]
      @param [in] _time: current physical time [Myr]
     */
    void initialFromFile(const std::string& _filename, const double& _time) {
        const int npar=17;
        double pars[npar];
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
        int my_rank = PS::Comm::getRank();
        if (my_rank==0) {
#endif
            std::ifstream fconf;
            fconf.open(_filename.c_str(), std::ifstream::in);
            if (!fconf.is_open()) {
                std::cerr<<"Error: MWPotentialEvolve configure file "<<_filename.c_str()<<" cannot be open!"<<std::endl;
                abort();
            }
            for (int i=0; i<npar; i++) {
                fconf>>pars[i];
                if (fconf.eof()) {
                    std::cerr<<"Error: MWPotentialEvolve require "<<npar<<" parameters, only "<<i<<" given!"<<std::endl;
                    abort();
                }
            }
            fconf.close();
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
        }
        PS::Comm::broadcast(pars, npar, 0);
#endif

        frw.time = pars[0];
        if (frw.time!=_time) {
            std::cerr<<"Error: the time recorded in the MWpotential parameter file, "<<frw.time<<", is not equal to the current system time, "<<_time<<"!"<<std::endl;
            abort();
        }
        frw.a = pars[1];
        frw.H0 = pars[2];
        frw.calcH0Myr();
        frw.omega_energy = pars[3];
        frw.omega_radiation = pars[4];
        frw.omega_matter = pars[5];
        frw.dt_scale = pars[6];
        frw.dt_max = pars[7];

        init.m_vir_halo = pars[8];
        init.c_halo = pars[9];
        init.ac_halo = pars[10];
        init.m_disk = pars[11];
        init.ra_disk = pars[12];
        init.rb_disk = pars[13];
        init.m_bulge = pars[14];
        init.alpha_bulge = pars[15];
        init.rcut_bulge = pars[16];
    }

    //! write data for restart
    /*! 
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
    */
    void writeData(std::ofstream & _fout) {
        _fout<<frw.time<<" "
             <<frw.a<<" "
             <<frw.H0<<" "
             <<frw.omega_energy<<" "
             <<frw.omega_radiation<<" "
             <<frw.omega_matter<<" "
             <<frw.dt_scale<<" "
             <<frw.dt_max<<" "
             <<init.m_vir_halo<<" "
             <<init.c_halo<<" "
             <<init.ac_halo<<" "
             <<init.m_disk<<" "
             <<init.ra_disk<<" "
             <<init.rb_disk<<" "
             <<init.m_bulge<<" "
             <<init.alpha_bulge<<" "
             <<init.rcut_bulge
             <<std::endl;
    }

    //! calculate MW potential based on the input time
    /*! @param[in] _time: physical time [Myr]
     */
    void calcMWPotentialParameters(const double& _time) {

        double dt = _time-frw.time;
        frw.updateAwithDtmin(dt);

        const double G_astro = 0.0044983099795944;  // Msun, pc, Myr
        const double pi = 3.141592653589793;

        double z = 1/frw.a-1; // redshift
        double Ht = frw.getH(); // Myr^-1
        double H0 = frw.H0_myr; // Myr^-1
        double mass_density_m1 = frw.getMatterDensity(frw.a)/frw.getDensity(frw.a) - 1;
        double mass_density_m1_z0 = frw.getMatterDensity(1.0) - 1;
        double delta = 18*pi*pi + 82*mass_density_m1 - 39*mass_density_m1*mass_density_m1;
        double delta0 = 18*pi*pi + 82*mass_density_m1_z0 - 39*mass_density_m1_z0*mass_density_m1_z0;
        now.rho_c_halo = 3*Ht*Ht/(8*pi*G_astro); // astro unit
        double rho_c_halo0 = 3*H0*H0/(8*pi*G_astro); // astro unit

        // halo, update Mvir, Rvir, c (Wechsler 2002, Zhao 2003) in astro unit
        double mass_z_factor = std::exp(-2*init.ac_halo*z);  
        now.m_vir_halo = init.m_vir_halo*mass_z_factor; // Msun

        now.r_vir_halo = std::pow(3*now.m_vir_halo/(4*pi*delta*now.rho_c_halo),1.0/3.0);
        double rvir_z_factor = std::pow(mass_z_factor*delta0*rho_c_halo0/(delta*now.rho_c_halo), 1.0/3.0);

        // Galpy parameter for halo in Galpy unit
        double c0 = init.c_halo;  // present-day concentration r_vir/r_s
        double c = c0/(z+1); // concentration at reshift z
        double c_factor = 1.0/(std::log(1+c) - c/(1+c));

        now.gm_halo = G_astro*now.m_vir_halo*c_factor;
        now.rs_halo = now.r_vir_halo/c;

        // disk & budge, update mass and radial scale factors (Bullock & Johnston 2005)
        now.gm_disk = G_astro*init.m_disk*mass_z_factor;
        now.ra_disk = init.ra_disk*rvir_z_factor;
        now.rb_disk = init.rb_disk*rvir_z_factor;

        now.gm_bulge = G_astro*init.m_bulge*mass_z_factor;
        now.alpha_bulge = init.alpha_bulge;
        now.rcut_bulge = init.rcut_bulge*rvir_z_factor;

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
    double gmscale;
    std::ifstream fconf;
    std::string set_name;
    std::string set_parfile;
    MWPotentialEvolve mw_evolve;

    GalpyManager(): potential_sets(), update_time(0.0), rscale(1.0), tscale(1.0), vscale(1.0), fscale(1.0), pscale(1.0), gmscale(1.0), fconf(), set_name(), set_parfile(), mw_evolve() {}

    //! initialization function
    /*!
      @param[in] _input: input parameters 
      @param[in] _time: current system time
      @param[in] _print_flag: if true, printing information to std::cout
     */
    void initial(const IOParamsGalpy& _input, const double _time, const bool _print_flag=false) {
        // unit scale
        rscale = _input.rscale.value;
        vscale = _input.vscale.value;
        tscale = _input.tscale.value;
        fscale = _input.fscale.value;
        pscale = _input.pscale.value;
        gmscale = pscale*rscale;

        // add pre-defined type-argu groups
        std::string type_args = _input.type_args.value;
        // Update types and arguments from type-args string
        bool initial_flag = addTypesAndArgsFromString(type_args, true, _print_flag);

        set_name = _input.pre_define_type.value;
        std::size_t ipar = set_name.find_first_of(":");
        if (ipar!=std::string::npos) {
            set_parfile = set_name.substr(ipar+1);
            set_name = set_name.substr(0, ipar);
        }

        if (set_name=="MWPotential2014") {
            int npot = 3;
            int pot_type[3] = {15,5,9};
            double pot_args[8] = {0.0299946, 1.8, 0.2375,  
                                  0.7574802, 0.375, 0.035, 
                                  4.85223053, 2.0};

            if (_print_flag) {
                std::cout<<"Galpy MWPotential2014 combination list:\n"
                         <<"Type index: 15 args: 0.0299946, 1.8, 0.2375\n"
                         <<"Type index: 5 args: 0.7574802, 0.375, 0.035\n"
                         <<"Type index: 9 args: 4.85223053, 2.0\n";
            }

            // generate galpy potential argument 
            potential_sets.push_back(PotentialSet());
            auto& pset = potential_sets.back();
            pset.setOrigin(0);
            pset.generatePotentialArgs(npot, pot_type, pot_args);

            initial_flag = true;
        }

        if (set_name=="MWPotentialEvolve") {
            /*
            int icount = 0;
            std::size_t istart = 0;
            std::size_t inext = set_par_str.find_first_of(",");
            while (inext!=std::string::npos) {
                set_par.push_back(std::stod(set_par_str.substr(istart, inext - istart)));
                istart = inext + 1;
                inext = set_par_str.find_first_of(",",istart);
                icount++;
            }
            if (istart!=set_par_str.size()) {
                icount++;
                set_par.push_back(std::stod(set_par_str.substr(istart)));
            }

            if (icount!=5) {
                std::cerr<<"Galpy Error: MWPotentialEvolve require 5 parameters, given only "<<icount<<"!"<<std::endl;
                abort();
            }
            */
            
            mw_evolve.initialFromFile(set_parfile, _time);

            updateMWPotentialEvolve(0, _print_flag);

            initial_flag = true;
        }        

        if (initial_flag && _input.config_filename.value!="__NONE__")  {
            std::cerr<<"Galpy Error: both --galpy-type-arg|--galpy-set and --galpy-conf-file are used, please choose one of them."<<std::endl;
            abort();
        }

        // add type arguments from configure file if exist
        if (_input.config_filename.value!="__NONE__") {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
            int my_rank = PS::Comm::getRank();
            if (my_rank==0) {
#endif
                fconf.open(_input.config_filename.value.c_str(), std::ifstream::in);
                if (!fconf.is_open()) {
                    std::cerr<<"Error: Galpy configure file "<<_input.config_filename.value.c_str()<<" cannot be open!"<<std::endl;
                    abort();
                }
                fconf>>update_time;
                if(fconf.eof()) fconf.close();
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
            }
            PS::Comm::broadcast(&update_time, 1, 0);
#endif
        }

        if(_print_flag) std::cout<<"----- Finish initialize Galpy potential -----\n";

    }


    //! Update types and arguments from type-args string
    /*! 
      @param[in] _type_args: potential type and argment strings 
      @param[in] _reset_flag: if true, reset potential set, otherwise add to current existing set
      @param[in] _print_flag: if true, print the read types and arguments.

      \return initial flag: if potential is updated, return true; else return false
     */
    bool addTypesAndArgsFromString(const std::string& _type_args, const bool _reset_flag, const bool _print_flag) {
        if (_type_args!="__NONE__") {

            std::vector<std::string> type_args_pair;
            std::vector<int> pot_type;
            std::vector<double> pot_args;

            // split type-arg groups
            std::size_t istart = 0;
            std::size_t inext = _type_args.find_first_of("|");
            while (inext!=std::string::npos) {
                type_args_pair.push_back(_type_args.substr(istart,inext-istart));
                istart = inext+1;
                inext = _type_args.find_first_of("|",istart);
            }
            if (istart!=_type_args.size()) type_args_pair.push_back(_type_args.substr(istart));

            if (_print_flag) std::cout<<"Galpy Potential combination list:\n";
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

            if (_reset_flag) freePotentialArgs();
            // generate galpy potential argument 
            potential_sets.push_back(PotentialSet());
            auto& pset = potential_sets.back();
            pset.setOrigin(0);
            pset.generatePotentialArgs(npot, pot_type.data(), pot_args.data());

            return true;
        }
        return false;
    }

    //! Update Potential types and arguments based on time
    /*!
      @param[in] _system_time: the time of globular particley system in PeTar. The time is transferred based on the time scaling factor
      @param[in] _print_flag: if true, print the read types and arguments.
     */
    void updatePotential(const double& _system_time, const bool _print_flag) {
        
        if (set_name=="MWPotentialEvolve") updateMWPotentialEvolve(_system_time, _print_flag);

        updateTypesAndArgsFromFile(_system_time, _print_flag);
    }

    //! write potential parameters for restart
    /*! Precision is set to 14
      @param[in] _fout std:ofstream to write data
      @param[in] _time current time
     */
    void writePotentialPars(std::ofstream & _fout, const double& _system_time) {
        _fout<<std::setprecision(14);
        if (set_name=="MWPotentialEvolve") {
            assert(mw_evolve.frw.time==_system_time);
            mw_evolve.writeData(_fout);
        }
    }

    //! Update MWpotential evolve (Gomez et al. 2010)
    /*! Assume time unit is Myr!
      @param[in] _system_time: the time of globular particley system in PeTar. The time is transferred based on the time scaling factor
      @param[in] _print_flag: if true, print the read types and arguments.
     */
    void updateMWPotentialEvolve(const double& _system_time, const bool _print_flag) {

        mw_evolve.calcMWPotentialParameters(_system_time);

        int npot = 3;
        int pot_type[3] = {15,5,9};
        auto& mwpot = mw_evolve.now;
        double pot_args[8] = {mwpot.gm_bulge*gmscale, mwpot.alpha_bulge, mwpot.rcut_bulge*rscale,
                              mwpot.gm_disk*gmscale, mwpot.ra_disk*rscale, mwpot.rb_disk*rscale,
                              mwpot.gm_halo*gmscale, mwpot.rs_halo*rscale};
        
        if (_print_flag) {
            std::cout<<std::setprecision(14)
                     <<"Galpy time[Myr]: "<<_system_time
                     <<" FRW a: "<<mw_evolve.frw.a
                     <<" NFW Halo: M_vir[Msun]: "<<mwpot.m_vir_halo<<" R_vir[pc]: "<<mwpot.r_vir_halo<<std::endl
                     <<"Update evolving MWPotential2014:\n"
                     <<"  Type index: 15 args: "<<pot_args[0]<<" "<<pot_args[1]<<" "<<pot_args[2]<<std::endl
                     <<"  Type index: 5  args: "<<pot_args[3]<<" "<<pot_args[4]<<" "<<pot_args[5]<<std::endl
                     <<"  Type index: 9  args: "<<pot_args[6]<<" "<<pot_args[7]<<std::endl;
        }

        freePotentialArgs();
        potential_sets.push_back(PotentialSet());
        auto& pset = potential_sets.back();
        pset.setOrigin(0);
        pset.generatePotentialArgs(npot, pot_type, pot_args);
        
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

                        pot_type_offset.push_back(0);
                        pot_args_offset.push_back(0);
                        if (nset>0) {
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


