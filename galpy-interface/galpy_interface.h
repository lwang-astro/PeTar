#pragma once
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <string>
#include <vector>
#include <getopt.h>
#include <cmath>
#include <algorithm>
#include "../src/io.hpp"

#include <integrateFullOrbit.h>

//! IO parameters for Galpy manager
class IOParamsGalpy{
public:
    IOParamsContainer input_par_store;
    IOParams<std::string> type_args; // potential type and argument list
    IOParams<std::string> pre_define_type; // pre defined type name
    IOParams<std::string> config_filename; // configure file name
    //IOParams<std::string> unit_set; // unscale, bovy 
    IOParams<double> rscale; 
    //IOParams<double> tscale; 
    IOParams<double> vscale; 
    //IOParams<double> fscale; 
    //IOParams<double> pscale; 
    
    bool print_flag;

    IOParamsGalpy(): input_par_store(),
                     type_args (input_par_store, "__NONE__", "galpy-type-arg", "Add potential types and argumentsto the potential list in the center of the galactic reference frame"),
                     pre_define_type (input_par_store, "__NONE__", "galpy-set", "Add a pre-defined potential set to the potential list, options are: MWPotential2014"),
                     config_filename(input_par_store, "__NONE__", "galpy-conf-file", "A configure file of the time- and space-dependent potentials"),
                     //unit_set(input_par_store, "unscale", "galpy-units", "Units conversion set: 'unscale': no conversion; all scale factors are 1.0; 'bovy': radial unit: 8 kpc; velocity unit: 220 km/s"),
                     rscale(input_par_store, 1.0, "galpy-rscale", "Radius scale factor from unit of the input particle data (IN) to Galpy distance unit (1.0)"),
                     //tscale(input_par_store, 1.0, "galpy-tscale", "Time scale factor (rscale/vscale) from unit of the input particle data (IN) to Galpy time (1.0)"),
                     vscale(input_par_store, 1.0, "galpy-vscale", "Velocity scale factor from unit of the input particle data (IN) to Galpy velocity unit (1.0)"),
                     //fscale(input_par_store, 1.0, "galpy-fscale", "Acceleration scale factor (vscale^2/rscale) from unit of the input particle data (IN) to Galpy acceleration unit (1.0)"),
                     //pscale(input_par_store, 1.0, "galpy-pscale", "Potential scale factor (vscale^2) from unit of the input particle data (IN) to Galpy potential unit (1.0)"),
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
            //{tscale.key,     required_argument, &galpy_flag, 4}, 
            {vscale.key,     required_argument, &galpy_flag, 5}, 
            //{fscale.key,     required_argument, &galpy_flag, 6}, 
            //{pscale.key,     required_argument, &galpy_flag, 7}, 
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
                //case 4:
                //    tscale.value = atof(optarg);
                //    if(print_flag) tscale.print(std::cout);
                //    opt_used+=2;
                //    break;
                case 5:
                    vscale.value = atof(optarg);
                    if(print_flag) vscale.print(std::cout);
                    opt_used+=2;
                    break;
                //case 6:
                //    fscale.value = atof(optarg);
                //    if(print_flag) fscale.print(std::cout);
                //    opt_used+=2;
                //    break;
                //case 7:
                //    pscale.value = atof(optarg);
                //    if(print_flag) pscale.print(std::cout);
                //    opt_used+=2;
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
                    std::cout<<"***PS: use petar.galpy.help to check how to setup --galpy-type-arg and --galpy-conf-file."
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

    //! set Bovy unit scale factor (solar distance and velocityr referring to the Galactic center) converge from astronomical unit set (Myr, Msun, pc)
    void setBovyUnit(const bool print_flag=true) {
        double kms_pcmyr=1.022712165045695; // pc = 30856775814913673 m, yr = 365.25 days
        double vbase=220.0;
        double rbase=8.0;
        double vb_pcmyr = vbase*kms_pcmyr;
        rscale.value = 0.001/rbase; // pc to solar position in Milkyway
        vscale.value = 1.0/vb_pcmyr; // pc/Myr to solar velocity in Milkyway

        if (print_flag) {
            double tscale = rscale.value/vscale.value; // myr to time unit in galpy
            double fscale = vscale.value*vscale.value/rscale.value; // pc/Myr^2 to acc unit in galpy
            double pscale = vscale.value*vscale.value; // pc^2/Myr^2 to pot unit in galpy
            double GMscale = pscale*rscale.value; // pc^3/Myr^2
            std::cout<<"----- Unit conversion for Galpy -----\n"
                     <<"rscale  = "<<rscale.value<<"  [8 kpc] / [pc]\n"
                     <<"vscale  = "<<vscale.value<<"  [220 km/s] / [pc/Myr]\n"
                     <<"tscale  = "<<tscale <<"  [Solar orbital period/2 pi] / Myr\n"
                     <<"fscale  = "<<fscale <<"  [galpy acceleration unit] / [pc/Myr^2]\n"
                     <<"pscale  = "<<pscale <<"  [galpy potential unit] / [pc^2/Myr^2]\n"
                     <<"GMscale = "<<GMscale<<"  [galpy gravitational constant * mass unit] / [pc^3/Myr^2]\n";
        }
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
/*! The evolving Milkyway potential referring to MWPotential2014 in Galpy
     Reference: Gomez F. A., Helmi A., Brown A. G. A., Li Y. S., 2010, MNRAS, 408, 935. doi:10.1111/j.1365-2966.2010.17225.x, https://ui.adsabs.harvard.edu/abs/2010MNRAS.408..935G
 */
class MWPotentialEvolve{
private:
    template <class Tstream>
    void labelCheck(Tstream& fconf, const char* match) {
        std::string label;
        fconf>>label;
        if (label!=match) {
            std::cerr<<"Galpy MWPotentialEvolve config: reading label error, should be "<<match<<" given "<<label<<std::endl;
            abort();
        }
    }

    void eofCheck(std::ifstream& fconf, const char* message) {
        if(fconf.eof()) {
            std::cerr<<"Galpy MWPotentialEvolve config: reading "<<message<<" fails! File reaches EOF."<<std::endl;
            abort();
        }
    }

public:
    FRWModel frw;
    struct {
        double m_vir_halo;
        double c_halo; // concentration
        double ac_halo; // formation epoch of halo (Wechsler 2002, Zhao 2003)
        double m_disk;
        double ra_disk;
        double rb_disk;
        double rho_bulge;
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
        double grho_bulge;
        double alpha_bulge;
        double rcut_bulge;
    } now;

    //! initial potential parameter from a file
    /*!
      Parameters in configure file:
                 Time [value(Myr)] a [value] H0 [value(km/s/Mpc)]
                 Omega [energy radiation matter]
                 dt_scale [value] dt_max [value(Myr)]
                 Halo [m_vir(Msun) c ac]                 
                 Disk [M(Msun) ra(pc) rb(pc)]
                 Bulge [rho(Msun/pc^(3-alpha)) alpha rcut(pc)]
      a (scale factor), H0 (Hubble constant) and Omega are used to obtain a at a given time by integrating adot(t)
      Omega are three normazied densities in FRW cosmological model
      dt_scale: scale of time to integrate a, if -1, evolve backwards, if 1, evolve forwards
      dt_max: maximum time step for integrating adot
      NFW Halo (3): virial mass, concentration, formation epoch of halo (Wechsler 2002, Zhao 2003)
      Miyamoto Nagai disk (3): mass, a, b
      Power Spherical bulge with expotential cutoff (3): rho0, power law index, cutoff radius

      @param [in] _filename: the file contain parameters 
      @param [in] _time: current physical time [Myr]
     */
    void initialFromFile(const std::string& _filename, const double& _time) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
        int my_rank = PS::Comm::getRank();
        if (my_rank==0) {
#endif
            std::ifstream fconf;
            fconf.open(_filename.c_str(), std::ifstream::in);
            if (!fconf.is_open()) {
                std::cerr<<"Galpy MWPotentialEvolve config: configure file "<<_filename.c_str()<<" cannot be open!"<<std::endl;
                abort();
            }

            labelCheck(fconf, "Time");
            fconf>>frw.time;
            eofCheck(fconf, "current time (Myr)");
            labelCheck(fconf, "a");
            fconf>>frw.a;
            eofCheck(fconf, "scale factor");
            labelCheck(fconf, "H0");
            fconf>>frw.H0;
            eofCheck(fconf, "Hubble constant");
            labelCheck(fconf, "Omega");
            fconf>>frw.omega_energy>>frw.omega_radiation>>frw.omega_matter;
            eofCheck(fconf, "cosmological densities");
            labelCheck(fconf, "dt_scale");
            fconf>>frw.dt_scale;
            eofCheck(fconf, "time scale");
            labelCheck(fconf, "dt_max");
            fconf>>frw.dt_max;

            eofCheck(fconf, "maximum time step (Myr)");
            labelCheck(fconf, "Halo");
            fconf>>init.m_vir_halo>>init.c_halo>>init.ac_halo;
            eofCheck(fconf, "Halo parameters");
            labelCheck(fconf, "Disk");
            fconf>>init.m_disk>>init.ra_disk>>init.rb_disk;
            eofCheck(fconf, "Disk parameters");
            labelCheck(fconf, "Bulge");
            fconf>>init.rho_bulge>>init.alpha_bulge>>init.rcut_bulge;
            eofCheck(fconf, "Bulge parameters");

            fconf.close();
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
        }
        PS::Comm::broadcast(&frw, 1, 0);
        PS::Comm::broadcast(&init, 1, 0);
#endif

        if (frw.time!=_time) {
            std::cerr<<"Error: the time recorded in the MWpotential parameter file, "<<frw.time<<", is not equal to the current system time, "<<_time<<"!"<<std::endl;
            abort();
        }

        frw.calcH0Myr();
    }

    //! read potential parameter from a file
    /*!
      @param [in] _filename: the file contain parameters (no label): 
                 FRWModel: Time[Myr], a, H0[km/s/Mpc], Omega_energy, Omega_radiation, Omega_matter, dt_scale
                 NFWHalo: M_vir[Msun], c, ac
                 Disk: M[Msun], Ra[pc], Rb[pc]
                 Bulge: M[Msun], alpha, rcut[pc]
      @param [in] _time: current physical time [Myr]
     */
    void readDataFromFile(const std::string& _filename, const double& _time) {
        const int npar=17;
        double pars[npar];
 #ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
         int my_rank = PS::Comm::getRank();
         if (my_rank==0) {
 #endif
             std::ifstream fconf;
             fconf.open(_filename.c_str(), std::ifstream::in);
             if (!fconf.is_open()) {
                 std::cerr<<"Galpy MWPotentialEvolve config: configure file "<<_filename.c_str()<<" cannot be open!"<<std::endl;
                 abort();
             }

             for (int i=0; i<npar; i++) {
                 fconf>>pars[i];
                 if (fconf.eof()) {
                     std::cerr<<"Galpy MWPotentialEvolve config: reading number of parameters is "<<npar<<" parameters, only "<<i<<" given!"<<std::endl;
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
         init.rho_bulge = pars[14];
         init.alpha_bulge = pars[15];
         init.rcut_bulge = pars[16];
    }

    //! write data for restart
    /*! 
      @param[in] _filename: file to save data
    */
    void writeDataToFile(const std::string& _filename) {
        std::ofstream fout;
        fout.open(_filename.c_str(), std::ifstream::out);
        if (!fout.is_open()) {
            std::cerr<<"Error: Galpy potential parameter file to write, "<<_filename<<", cannot be open!"<<std::endl;
            abort();
        }
        fout<<std::setprecision(14)
            <<frw.time<<" "
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
            <<init.rho_bulge<<" "
            <<init.alpha_bulge<<" "
            <<init.rcut_bulge
            <<std::endl;
        fout.close();
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

        now.grho_bulge = G_astro*init.rho_bulge*mass_z_factor/std::pow(rvir_z_factor,3-init.alpha_bulge);
        now.alpha_bulge = init.alpha_bulge;
        now.rcut_bulge = init.rcut_bulge*rvir_z_factor;

    }
};

//! mode, position, velocity and acceleration for set of Galpy potentials
struct PotentialSetPar{
    int mode; // mode of origin: 0: galactic frame; 1: local particle system frame
    double gm; // G*mass of potential (used for computing acc from particle system)
    double pos[3]; // position
    double vel[3]; // velocity
    double acc[3]; // acceleration

    PotentialSetPar(): mode(-1), gm(0), pos{0.0,0.0,0.0}, vel{0.0,0.0,0.0}, acc{0.0,0.0,0.0} {}

    //! print data in one line
    void writeData(std::ostream& fout) {
        fout<<mode<<" "
            <<gm<<" ";
        for (std::size_t i=0; i<3; i++) fout<<pos[i]<<" ";
        for (std::size_t i=0; i<3; i++) fout<<vel[i]<<" ";
        for (std::size_t i=0; i<3; i++) fout<<acc[i]<<" ";
    }

    void readData(std::ifstream& fin) {
        fin>>mode>>gm;
        fin>>pos[0]>>pos[1]>>pos[2];
        fin>>vel[0]>>vel[1]>>vel[2];
        fin>>acc[0]>>acc[1]>>acc[2];
    }

    //! set position,  velocity and acceleration
    void setOrigin(const int _mode, const double _gm=0, const double* _pos=NULL, const double* _vel=NULL, const double* _acc=NULL) {
        assert(_mode>=0&&_mode<=2);
        mode = _mode;
        gm = _gm;
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
        if (_acc!=NULL) {
            acc[0] = _acc[0];
            acc[1] = _acc[1];
            acc[2] = _acc[2];
        }
    }

    void clear() {
        mode = -1;
        gm = 0.0;
        pos[0] = pos[1] = pos[2] = 0.0;
        vel[0] = vel[1] = vel[2] = 0.0;
        acc[0] = acc[1] = acc[2] = 0.0;
    }
    
};

//! set of Galpy potentials sharing the same position and velocity of the origin
struct PotentialSet{
    int npot; // number of potential models 
    struct potentialArg* arguments; //potential arguments array for Gaply

    PotentialSet(): npot(0), arguments(NULL) {}

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
        npot = 0;
        arguments = NULL;
    }

    ~PotentialSet() { 
        clear(); 
    }
};

//! changing argument parameter
struct ChangeArgument{
    int index; // index of argument for change
    int mode; // change mode (1: linear; 2: exponential)
    double rate; // change rate: linear: cofficient a in a*t; exp: cofficient a in exp(a*t)

    //! print data in one line
    void writeData(std::ostream& fout) {
        fout<<index<<" "<<mode<<" "<<rate<<" ";
    }

    //! print data in one line
    void readData(std::ifstream& fin) {
        fin>>index>>mode>>rate;
    }
};

//! A class to manager the API to Galpy
class GalpyManager{
private:
    template <class Tstream>
    void labelCheck(Tstream& fconf, const char* match) {
        std::string label;
        fconf>>label;
        if (label!=match) {
            std::cerr<<"Galpy config: reading label error, should be "<<match<<" given "<<label<<std::endl;
            abort();
        }
    }

    void eofCheck(std::ifstream& fconf, const char* message) {
        if(fconf.eof()) {
            std::cerr<<"Galpy config: reading "<<message<<" fails! File reaches EOF."<<std::endl;
            abort();
        }
    }

    void gtZeroCheck(const int n, const char* message) {
        if (n<=0) {
            std::cerr<<"Galpy config: "<<message<<" <=0!"<<std::endl;
            abort();
        }
    }
    
    void geZeroCheck(const int n, const char* message) {
        if (n<0) {
            std::cerr<<"Galpy config: "<<message<<" <0!"<<std::endl;
            abort();
        }
    }


    //! resize array by insert or erase elements for given offset index
    /*!
      @param[in] n_diff: size difference after update
      @param[in] index: index of array offset for change
      @param[in,out] array: array of data
      @param[in,out] array_offset: offset of differert data groups 
     */
    template <class ttype>
    void resizeArray(const int n_diff, const int index, std::vector<ttype>& array, std::vector<int>& array_offset) {
        int offset = array_offset[index];
        if (n_diff!=0) {
            if (n_diff>0) {
                std::vector<ttype> data(n_diff,ttype());
                array.insert(array.begin()+offset, data.begin(), data.end());
            }
            else 
                array.erase(array.begin()+offset, array.begin()+offset-n_diff);
            for (size_t i=index+1; i<array_offset.size(); i++) 
                array_offset[i] += n_diff;
        }
    }

    //! erase array for one offset index
    /*!
      @param[in] index: index of array offset for remove
      @param[in,out] array: array of data
      @param[in,out] array_offset: offset of differert data groups 
    */
    template <class ttype>
    void eraseArray(const int index, std::vector<ttype>& array, std::vector<int>& array_offset) {
        int n = array_offset[index+1] - array_offset[index];
        int offset = array_offset[index];
        if (n>0) {
            array.erase(array.begin()+offset, array.begin()+offset+n);
        }
        for (size_t i=index+1; i<array_offset.size(); i++) 
            array_offset[i] = array_offset[i+1]-n;
        array_offset.pop_back();
        assert(array_offset.back()==int(array.size()));
    }

        
public:

    // for MPI communication, data IO and initialization
    std::vector<int> pot_type_offset; //  set offset
    std::vector<int> pot_type;    // types for each set
    std::vector<int> pot_args_offset; // arguments of pot for each set
    std::vector<double> pot_args;    // set offset
    double time;   // current time, update in evolveChangingArguments
    std::vector<ChangeArgument> change_args;  // changing argument index to evolve
    std::vector<int> change_args_offset; // set offset
    // galpy potential arguments
    std::vector<PotentialSetPar> pot_set_pars; // potential parameters for each set
    std::vector<PotentialSet> pot_sets; // potential arguments for each set
    double update_time;
    // unit scaling
    double rscale;
    double vscale;
    double tscale;
    double fscale;
    double pscale;
    double gmscale;
    // IO
    std::ifstream fconf;
    std::string set_name;
    std::string set_parfile;
    // evolving Milkyway potential
    MWPotentialEvolve mw_evolve;

    GalpyManager(): pot_type_offset(), pot_type(), 
                    pot_args_offset(), pot_args(), change_args(), change_args_offset(),
                    pot_set_pars(), pot_sets(), update_time(0.0), rscale(1.0), vscale(1.0), tscale(1.0), fscale(1.0), pscale(1.0), gmscale(1.0), fconf(), set_name(), set_parfile(), mw_evolve() {}

    //! print current potential data
    void printData(std::ostream& fout) {
        fout<<"Galpy parameters, time: "<<time;
        fout<<" Next update time: "<<update_time;
        fout<<std::endl;
        int nset = pot_set_pars.size();
        for (int k=0; k<nset; k++) {
            auto& pot_set_par_k = pot_set_pars[k];
            fout<<"Potential set "<<k+1<<" Mode: "<<pot_set_par_k.mode
                <<" GM: "<<pot_set_par_k.gm
                <<" Pos: "<<pot_set_par_k.pos[0]<<" "<<pot_set_par_k.pos[1]<<" "<<pot_set_par_k.pos[2]
                <<" Vel: "<<pot_set_par_k.vel[0]<<" "<<pot_set_par_k.vel[1]<<" "<<pot_set_par_k.vel[2]
                <<" Acc: "<<pot_set_par_k.acc[0]<<" "<<pot_set_par_k.acc[1]<<" "<<pot_set_par_k.acc[2]
                <<"\nPotential type indice: ";
            for (int i=pot_type_offset[k]; i<pot_type_offset[k+1]; i++) 
                fout<<pot_type[i]<<" ";
            fout<<"\nPotential arguments: ";
            for (int i=pot_args_offset[k]; i<pot_args_offset[k+1]; i++) 
                fout<<pot_args[i]<<" ";
            fout<<"\nChange argument [index mode rate]:";
            for (int i=change_args_offset[k]; i<change_args_offset[k+1]; i++) {
                fout<<"["<<change_args[i].index
                    <<" "<<change_args[i].mode
                    <<" "<<change_args[i].rate<<"] ";
            }
            fout<<std::endl;

            if (set_name=="MWPotentialEvolve") 
                fout<<"Scale factor a: "<<mw_evolve.frw.a<<std::endl;
        }
    }        

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
    void broadcastDataMPI() {

        int nset;
        int my_rank = PS::Comm::getRank();
        if (my_rank==0) nset = pot_set_pars.size();
        PS::Comm::broadcast(&nset, 1, 0);

        // update potentials
        if (nset>=0) {

            PS::Comm::broadcast(&update_time, 1, 0);
            if (my_rank==0) {
                assert((int)pot_type_offset.size()==nset+1);
                assert((int)pot_args_offset.size()==nset+1);
                assert((int)change_args_offset.size()==nset+1);
            }
            else {
                pot_set_pars.resize(nset);
                pot_type_offset.resize(nset+1);
                pot_args_offset.resize(nset+1);
                change_args_offset.resize(nset+1);
            }
            PS::Comm::broadcast(pot_set_pars.data(), pot_set_pars.size(), 0);
            PS::Comm::broadcast(pot_type_offset.data(), pot_type_offset.size(), 0);
            PS::Comm::broadcast(pot_args_offset.data(), pot_args_offset.size(), 0);
            PS::Comm::broadcast(change_args_offset.data(), change_args_offset.size(), 0);
            
            if (my_rank==0) {
                assert((int)pot_type.size()==pot_type_offset.back());
                assert((int)pot_args.size()==pot_args_offset.back());
                assert((int)change_args.size()==change_args_offset.back());
            }
            else {
                pot_type.resize(pot_type_offset.back());
                pot_args.resize(pot_args_offset.back());
                change_args.resize(change_args_offset.back());
            }
            PS::Comm::broadcast(pot_type.data(), pot_type.size(), 0);
            if (pot_args.size()>0) PS::Comm::broadcast(pot_args.data(), pot_args.size(), 0);
            if (change_args.size()>0) PS::Comm::broadcast(change_args.data(), change_args.size(), 0);
        }
    }

#endif

    //! initialization function
    /*!
      @param[in] _input: input parameters 
      @param[in] _time: current system time
      @param[in] _conf_name: galpy configure file name for restart
      @param[in] _restart_flag: if true, read the configure file for restart
      @param[in] _print_flag: if true, printing information to std::cout
     */
    void initial(const IOParamsGalpy& _input, const double _time, const std::string _conf_name, const bool _restart_flag, const bool _print_flag=false) {
        // unit scale
        rscale = _input.rscale.value;
        vscale = _input.vscale.value;
        tscale = rscale/vscale;
        fscale = vscale*vscale/rscale;
        pscale = vscale*vscale;
        gmscale = pscale*rscale;

        // initial offset
        pot_type_offset.push_back(0);
        pot_args_offset.push_back(0);
        change_args_offset.push_back(0);

        set_name = _input.pre_define_type.value;
        std::size_t ipar = set_name.find_first_of(":");
        if (ipar!=std::string::npos) {
            set_parfile = set_name.substr(ipar+1);
            set_name = set_name.substr(0, ipar);
        }

        // if restart, read the configure file from _conf_name instead of _input.pre_define_type.
        if (_restart_flag) set_parfile = _conf_name;

        if (set_name=="MWPotential2014") {
                
            pot_set_pars.emplace_back();
            auto& pot_set_par_k = pot_set_pars.back();
            pot_set_par_k.setOrigin(0);

            pot_type.insert(pot_type.end(), {15,5,9});
            pot_type_offset.push_back(pot_type.size());
            
            // Bovy unit
            //double pot_args[8] = {0.0299946, 1.8, 0.2375,  
            //                      0.7574802, 0.375, 0.035, 
            //                      4.85223053, 2.0};

            // Astro unit
            // PowerSphericalPotentialwithCutoff bulge 
            //   rho(r) = (3-a) G*M_g / (4 pi R_g^(3-a)) * (r_1/r)^a * exp ((-r/rc)^2)
            //      args(3): (3-a) G*M_g / (4 pi R_g^(3-a)), a, rc (r_1 = 1)
            //      unit: GMscale/rscale^(3-a), 1, rscale
            // MiyamotoNagai disk
            //   Phi(R,z) = -G*M_g / sqrt(R^2+(a+sqrt(z^2+b^2))^2)
            //      args(3): G*M_g, a, b
            //      unit: GMscale, rscale, rscale
            // NFW halo
            //   rho(r) = G*M_g /(4 pi a^3) * 1/((r/a)*(1+r/a)^2)
            //      args(2): G*M_g, a
            //      unit: GMscale, rscale
            pot_args.insert(pot_args.end(), {251.63858935563147, 1.8, 1899.9999999999998, 
                                             306770418.38588977, 3000.0, 280.0,  
                                             1965095308.192175, 16000.0});
            pot_args_offset.push_back(pot_args.size());

            change_args_offset.push_back(0);

        }

        // MWpotentialEvolve should be the first set
        if (set_name=="MWPotentialEvolve") {
            
            if (_restart_flag) 
                mw_evolve.readDataFromFile(set_parfile, _time);
            else
                mw_evolve.initialFromFile(set_parfile, _time);
            updateMWPotentialEvolve(_time, _print_flag, true);
        }        

        // add pre-defined type-argu groups
        std::string type_args = _input.type_args.value;
        // Update types and arguments from type-args string
        addPotentialFromString(type_args, true, _print_flag);

        // add type arguments from configure file if exist
        if (_input.config_filename.value!="__NONE__") {
            set_name="configure";

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
            int my_rank = PS::Comm::getRank();
            if (my_rank==0) {
#endif
                fconf.open(_input.config_filename.value.c_str(), std::ifstream::in);
                if (!fconf.is_open()) {
                    std::cerr<<"Error: Galpy configure file "<<_input.config_filename.value.c_str()<<" cannot be open!"<<std::endl;
                    abort();
                }
                labelCheck(fconf, "Time");
                if(fconf.eof()) fconf.close();
                fconf>>update_time;

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
            }
            PS::Comm::broadcast(&update_time, 1, 0);
#endif
            // if restart, read configure file at the current time
            if (_restart_flag) {
                // first escape past time to avoid overwritting restart configure file
                updatePotentialFromFile(_time, false);
                readDataFromFile(set_parfile, _print_flag);
            }
        }

        updatePotentialSet();
        
        if(_print_flag) {
            printData(std::cout);
            std::cout<<"----- Finish initialize Galpy potential -----\n";
        }

    }

    //! generate potentail args
    /*!
       clear the old potentialset first
    */
    void updatePotentialSet() {
        freePotentialArgs();
        int nset = pot_set_pars.size(); 
        pot_sets.resize(nset);
        for (int k=0; k<nset; k++) {
            pot_sets[k].generatePotentialArgs(pot_type_offset[k+1]-pot_type_offset[k], &(pot_type[pot_type_offset[k]]), &(pot_args[pot_args_offset[k]]));
        }
    }

    //! Update types and arguments from type-args string
    /*! 
      @param[in] _type_args: potential type and argment strings 
      @param[in] _reset_flag: if true, reset potential set, otherwise add to current existing set
      @param[in] _print_flag: if true, print the read types and arguments.

      \return initial flag: if potential is updated, return true; else return false
     */
    bool addPotentialFromString(const std::string& _type_args, const bool _reset_flag, const bool _print_flag) {
        if (_type_args!="__NONE__") {
            
            pot_set_pars.emplace_back();
            auto& pot_set_par_k = pot_set_pars.back();
            pot_set_par_k.setOrigin(0);

            std::vector<std::string> type_args_pair;

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
            pot_type_offset.push_back(pot_type.size());
            pot_args_offset.push_back(pot_args.size());
            change_args_offset.push_back(0);

            return true;
        }
        return false;
    }

    //! Update MWpotential evolve (Gomez et al. 2010)
    /*! Assume time unit is Myr!
      @param[in] _system_time: the time of globular particley system in PeTar. The time is transferred based on the time scaling factor
      @param[in] _print_flag: if true, print the read types and arguments.
      @param[in] _initial_flag: if true, this is the first time, potential parameter array will be generated.
      \return true if potential is updated (always) else false
     */
    bool updateMWPotentialEvolve(const double& _system_time, const bool _print_flag, const bool _initial_flag) {

        mw_evolve.calcMWPotentialParameters(_system_time);
        auto& mwpot = mw_evolve.now;
        double pot_arg_mw[8] = {mwpot.grho_bulge*gmscale, mwpot.alpha_bulge, mwpot.rcut_bulge*rscale,
                            mwpot.gm_disk*gmscale, mwpot.ra_disk*rscale, mwpot.rb_disk*rscale,
                            mwpot.gm_halo*gmscale, mwpot.rs_halo*rscale};

        if (_initial_flag) {
            pot_set_pars.emplace_back();
            auto& pot_set_par_k = pot_set_pars.back();
            pot_set_par_k.setOrigin(0);

            pot_type.insert(pot_type.end(), {15,5,9});
            pot_type_offset.push_back(pot_type.size());

            pot_args.insert(pot_args.end(), pot_arg_mw, pot_arg_mw+8);
            pot_args_offset.push_back(pot_args.size());
        
            change_args_offset.push_back(0);
        }
        else{
            assert(pot_type[0]==15);
            assert(pot_type[1]==5);
            assert(pot_type[2]==9);
            for (int i=0; i<8; i++) pot_args[i] = pot_arg_mw[i];
        }
        
        if (_print_flag) {
            std::cout<<std::setprecision(14)
                     <<" Evolve MilkyWay potential: FRW a: "<<mw_evolve.frw.a
                     <<" NFW Halo: M_vir[Msun]: "<<mwpot.m_vir_halo<<" R_vir[pc]: "<<mwpot.r_vir_halo<<std::endl;
        }

        return true;
    }

    //! Update types and arguments from configure file if update_time <= particle system time
    /*! 
      @param[in] _system_time: the time of globular particley system in PeTar. The time is transferred based on the time scaling factor
      @param[in] _print_flag: if true, print the read types and arguments.

      \return true if potential is updated; else false
     */
    bool updatePotentialFromFile(const double& _system_time, const bool _print_flag) {
        double time_scaled= _system_time*tscale;
        int nset=-1;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
        int my_rank = PS::Comm::getRank();
        if (my_rank==0) {
#endif
            if (fconf.is_open())  {
                if (time_scaled>=update_time) {

                    double update_time_next;

                    while(true) {
                        labelCheck(fconf, "Task");
                        std::string task;
                        fconf>>task;
                        eofCheck(fconf, "task");                                             int n_set_old = pot_set_pars.size();
                        if (task=="add") {
                            labelCheck(fconf, "Nset");
                            fconf>>nset;
                            eofCheck(fconf, "number of new potential sets");
                            gtZeroCheck(nset, "number of new potential sets");

                            for (int k=0; k<nset; k++) {
                                pot_set_pars.push_back(PotentialSetPar());
                                auto& pot_set_par_k = pot_set_pars.back();
                                // set 
                                labelCheck(fconf, "Set");
                                int set_index;
                                fconf>>set_index;
                                eofCheck(fconf, "set index");
                                if (set_index!=k) {
                                    std::cerr<<"Galpy config: reading set index not match the expected one, should be "<<k<<" given "<<set_index<<std::endl;
                                    abort();
                                }
                                // Ntype
                                labelCheck(fconf, "Ntype");
                                int n_pot_k=-1;
                                fconf>>n_pot_k;
                                eofCheck(fconf, "number of potential types");
                                gtZeroCheck(n_pot_k, "number of potential types");

                                // mode
                                labelCheck(fconf, "Mode");
                                int mode_k;
                                fconf>>mode_k;
                                eofCheck(fconf, "mode");
                                if (mode_k<0 || mode_k>2) {
                                    std::cerr<<"Galpy config: mode should be 0, 1, 2, given "<<mode_k<<std::endl;
                                    abort();
                                }
                                pot_set_par_k.mode = mode_k;

                                // Mass
                                labelCheck(fconf, "GM");
                                fconf>>pot_set_par_k.gm;
                                eofCheck(fconf, "g*m");
                                // Pos
                                labelCheck(fconf, "Pos");
                                for (int j=0; j<3; j++) fconf>>pot_set_par_k.pos[j];
                                eofCheck(fconf, "position");
                                // Vel
                                labelCheck(fconf, "Vel");
                                for (int j=0; j<3; j++) fconf>>pot_set_par_k.vel[j];
                                eofCheck(fconf, "velocity");

                                // type, arg and change_args are not directly save to potential set because of MPI communication
                                // types
                                labelCheck(fconf, "Type");
                                for (int i=0; i<n_pot_k; i++) {
                                    int type_i;
                                    fconf>>type_i;
                                    eofCheck(fconf, "potential type");
                                    pot_type.push_back(type_i);
                                }
                                pot_type_offset.push_back(pot_type.size());

                                // arguments
                                std::string line;
                                std::getline(fconf, line); // skip 2nd line
                                std::getline(fconf, line);
                                std::istringstream fin(line);
                                double arg_i;
                                int n_arg=0;
                                labelCheck(fin, "Arg");
                                while (fin>>arg_i) {
                                    pot_args.push_back(arg_i);
                                    n_arg++;
                                }
                                pot_args_offset.push_back(pot_args.size());

                                // time-dependent arguments
                                labelCheck(fconf, "Nchange");
                                int n_change;
                                fconf>>n_change;
                                eofCheck(fconf, "number of changing arguments");
                                geZeroCheck(n_change, "number of changing arguments");
                                int offset = change_args.size();
                                if (n_change>0) {
                                    change_args.resize(offset+n_change);
                                    labelCheck(fconf, "Index");
                                    for (int i=0; i<n_change; i++) {
                                        int index_i;
                                        fconf>>index_i;
                                        eofCheck(fconf, "index of the changing argument");
                                        if (index_i<-1 || index_i>=n_arg) {
                                            std::cerr<<"Galpy config: index of the changing argument should be >=-1 and <"<<n_arg<<", given "<<index_i<<std::endl;
                                            abort();
                                        }
                                        change_args[offset+i].index = index_i;
                                    }
                                    labelCheck(fconf, "ChangeMode");
                                    for (int i=0; i<n_change; i++) {
                                        int change_mode_i;
                                        fconf>>change_mode_i;
                                        eofCheck(fconf, "changing mode");
                                        if (!(change_mode_i==1||change_mode_i==2)) {
                                            std::cerr<<"Galpy config: change mode must be 1 or 2, given "<<change_mode_i<<std::endl;
                                            abort();
                                        }
                                        change_args[offset+i].mode = change_mode_i;
                                    }
                                    labelCheck(fconf, "ChangeRate");
                                    for (int i=0; i<n_change; i++) {
                                        double rate_i;
                                        fconf>>rate_i;
                                        eofCheck(fconf, "changing rate");
                                        change_args[offset+i].rate = rate_i;
                                    }
                                }
                                change_args_offset.push_back(change_args.size());
                            }
                            nset += n_set_old;
                        }
                        else if (task=="update") {
                            labelCheck(fconf, "Nset");
                            fconf>>nset;
                            eofCheck(fconf, "number of update potential sets");
                            gtZeroCheck(nset, "number of update potential sets");

                            for (int i=0; i<nset; i++) {
                                // set 
                                labelCheck(fconf, "Set");
                                int k;
                                fconf>>k;
                                eofCheck(fconf, "set index");
                                if (k<0 || k>=n_set_old) {
                                    std::cerr<<"Galpy config: changing set index incorrect, should be from 0 to "<<n_set_old-1<<", given "<<k<<std::endl;
                                    abort();
                                }
                                auto& pot_set_par_k = pot_set_pars[k];

                                // Ntype
                                labelCheck(fconf, "Ntype");
                                int n_pot_k=-1;
                                fconf>>n_pot_k;
                                eofCheck(fconf, "number of potential types");
                                gtZeroCheck(n_pot_k, "number of potential types");

                                // mode
                                labelCheck(fconf, "Mode");
                                int mode_k;
                                fconf>>mode_k;
                                eofCheck(fconf, "mode");
                                if (mode_k<-2 || mode_k==-1 || mode_k>2) {
                                    std::cerr<<"Galpy config: mode should be 0, 1, 2 or -2, given "<<mode_k<<std::endl;
                                    abort();
                                }

                                // mass
                                labelCheck(fconf, "GM");
                                fconf>>pot_set_par_k.gm;
                                eofCheck(fconf, "gm");

                                // update position and velocity for mode 2
                                if (mode_k!=-2) {
                                    labelCheck(fconf, "Pos");
                                    for (int j=0; j<3; j++) fconf>>pot_set_par_k.pos[j];
                                    eofCheck(fconf, "position");
                                    labelCheck(fconf, "Vel");
                                    for (int j=0; j<3; j++) fconf>>pot_set_par_k.vel[j];
                                    eofCheck(fconf, "velocity");
                                }
                                else {
                                    double dtmp;
                                    labelCheck(fconf, "Pos");
                                    for (int j=0; j<3; j++) fconf>>dtmp;
                                    eofCheck(fconf, "position");
                                    labelCheck(fconf, "Vel");
                                    for (int j=0; j<3; j++) fconf>>dtmp;
                                    eofCheck(fconf, "velocity");
                                    mode_k = 2;
                                }
                                pot_set_par_k.mode = mode_k;
                                
                                // types
                                // modify array offset if number of pot in the set changes
                                int n_pot_diff = n_pot_k - (pot_type_offset[k+1]-pot_type_offset[k]);
                                resizeArray(n_pot_diff, k, pot_type, pot_type_offset);
                                int offset = pot_type_offset[k];
                                
                                labelCheck(fconf, "Type");
                                for (int j=0; j<n_pot_k; j++) {
                                    int type_j;
                                    fconf>>type_j;
                                    eofCheck(fconf, "potential type");
                                    pot_type[offset+j] = type_j;
                                }

                                // arguments
                                std::string line;
                                std::getline(fconf, line); // skip 2nd line
                                std::getline(fconf, line);
                                std::istringstream fin(line);
                                labelCheck(fin, "Arg");
                                std::vector<double> arg_k;
                                double arg_tmp;
                                while (fin>>arg_tmp) {
                                    arg_k.push_back(arg_tmp);
                                }
                                int n_arg_diff = int(arg_k.size()) - (pot_args_offset[k+1] - pot_args_offset[k]);
                                resizeArray(n_arg_diff, k, pot_args, pot_args_offset);
                                offset = pot_args_offset[k];
                                for (size_t j=0; j<arg_k.size(); j++) {
                                    pot_args[offset+j] = arg_k[j];
                                }

                                // time-dependent arguments
                                labelCheck(fconf, "Nchange");
                                int n_change;
                                fconf>>n_change;
                                eofCheck(fconf, "number of changing arguments");
                                geZeroCheck(n_change, "changing arguments");
                                int n_change_diff = n_change - (change_args_offset[k+1] - change_args_offset[k]);

                                resizeArray(n_change_diff, k, change_args, change_args_offset);
                                if (n_change>0) {
                                    offset = change_args_offset[k];
                                    labelCheck(fconf, "Index");
                                    for (int j=0; j<n_change; j++) {
                                        int index_j;
                                        fconf>>index_j;
                                        eofCheck(fconf, "index of the changing argument");
                                        if (index_j<-1 || index_j>=int(arg_k.size())) {
                                            std::cerr<<"Galpy config: index of the changing argument should be >=-1 and <"<<(arg_k.size())<<", given "<<index_j<<std::endl;
                                            abort();
                                        }
                                        change_args[offset+j].index = index_j;
                                    }
                                    labelCheck(fconf, "ChangeMode");
                                    for (int j=0; j<n_change; j++) {
                                        int change_mode_j;
                                        fconf>>change_mode_j;
                                        eofCheck(fconf, "changing mode");
                                        if (!(change_mode_j==1||change_mode_j==2)) {
                                            std::cerr<<"Galpy config: change mode must be 1 or 2, given "<<change_mode_j<<std::endl;
                                            abort();
                                        }
                                        change_args[offset+j].mode = change_mode_j;
                                    }
                                    labelCheck(fconf, "ChangeRate");
                                    for (int j=0; j<n_change; j++) {
                                        double rate_j;
                                        fconf>>rate_j;
                                        eofCheck(fconf, "changing rate");
                                        change_args[offset+j].rate = rate_j;
                                    }                                    
                                }
                            }
                            nset = n_set_old;
                        }
                        else if (task=="remove") {
                            labelCheck(fconf, "Nset");
                            fconf>>nset;
                            eofCheck(fconf, "number of removing potential sets");
                            gtZeroCheck(nset, "number of removing potential sets");

                            labelCheck(fconf, "Index");
                            std::vector<int> rm_indices;
                            for (int i=0; i<nset; i++) {
                                int k;
                                fconf>>k;
                                eofCheck(fconf, "removing set index");
                                rm_indices.push_back(k);
                            }
                            std::sort(rm_indices.begin(), rm_indices.end());
                            for (int i=0; i<nset; i++) {
                                int idel = rm_indices[i]-i; // each time of remove, the following indices reduce by one, total reducing number is i
                                pot_set_pars.erase(pot_set_pars.begin()+idel);
                                eraseArray(idel, pot_type, pot_type_offset);
                                eraseArray(idel, pot_args, pot_args_offset);
                                eraseArray(idel, change_args, change_args_offset);
                            }
                            nset = n_set_old - nset;
                            geZeroCheck(nset, "number of potential sets after remove");
                        }
                        else {
                            std::cerr<<"Galpy config: Task must be add, update or remove, given "<<task<<std::endl;
                            abort();
                        }

                        std::string label;
                        fconf>>label;
                        if(fconf.eof()) {
                            fconf.close();
                            break;
                        }
                        if(label!="Time") {
                            std::cerr<<"Galpy config: reading label error, should be Time given "<<label<<std::endl;
                            abort();
                        }
                        fconf>>update_time_next;
                        if (update_time_next>time_scaled) {
                            update_time = update_time_next;
                            break;
                        }
                    };
                }
            }
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
        }
        broadcastDataMPI();
#endif        
        return (nset>=0);
    }

    //! Update Potential types and arguments based on time
    /*!
      @param[in] _system_time: the time of globular particley system in PeTar. The time is transferred based on the time scaling factor
      @param[in] _print_flag: if true, print the read types and arguments.
     */
    void updatePotential(const double& _system_time, const bool _print_flag) {

        int nset = pot_set_pars.size();
        assert((int)pot_type_offset.size()==nset+1);
        assert((int)pot_args_offset.size()==nset+1);
        assert((int)change_args_offset.size()==nset+1);
        assert((int)pot_type.size()==pot_type_offset.back());
        assert((int)pot_args.size()==pot_args_offset.back());
        assert((int)change_args.size()==change_args_offset.back());

        bool update_flag = false;
        if (set_name=="MWPotentialEvolve") {
            update_flag = updateMWPotentialEvolve(_system_time, false, false);
        }
        else if (set_name=="configure") {
            update_flag = updatePotentialFromFile(_system_time, _print_flag);
            // update time-dependent argument
            bool change_flag = evolveChangingArguments(_system_time*tscale, false);
            if (update_flag && _print_flag) printData(std::cout);
            update_flag = (update_flag || change_flag);
        }
        
        // generate galpy potential argument 
        if (update_flag) updatePotentialSet();
    }

    //! write potential parameters for restart
    /*! Precision is set to 14
      @param[out] _filename: file to save data
      @param[in] _time current time
     */
    void writePotentialPars(const std::string& _filename, const double& _system_time) {
        if (set_name=="MWPotentialEvolve") {
            assert(mw_evolve.frw.time==_system_time);
            mw_evolve.writeDataToFile(_filename);
        }
        else if (set_name=="configure") {
            assert(time==_system_time*tscale);
            writeDataToFile(_filename);
        }
    }

    //! print reference to cite
    static void printReference(std::ostream & fout, const int offset=4) {
        for (int i=0; i<offset; i++) fout<<" ";
        fout<<"Galpy: Bovy J., 2015, ApJS, 216, 29"<<std::endl;
    }    

    //! reset acceleraltion of potential 
    void resetPotAcc() {
        for (size_t k=0; k<pot_set_pars.size(); k++) {
            double* acc_pot = pot_set_pars[k].acc;
            acc_pot[0] = acc_pot[1] = acc_pot[2] = 0.0;
        }        
    }

    //! kick velocity of moving potential (mode 2)
    /*!
      @param[in] dt: time step for kick [input unit]
    */
    void kickMovePot(const double dt) {
        double dt_scale = dt*tscale;
        for (size_t k=0; k<pot_set_pars.size(); k++) {
            int mode_k = pot_set_pars[k].mode;
            assert(mode_k>=0||mode_k<=2);
            if (mode_k==2) {
                double* vel_k = pot_set_pars[k].vel;
                double* acc_k = pot_set_pars[k].acc;
                vel_k[0] += acc_k[0]*dt_scale;
                vel_k[1] += acc_k[1]*dt_scale;
                vel_k[2] += acc_k[2]*dt_scale;
            }        
        }
    }

    //! drift position of moving potential (mode 2)
    /*!
      @param[in] dt: time step for kick [input unit]
    */
    void driftMovePot(const double dt) {
        double dt_scale = dt*tscale;
        for (size_t k=0; k<pot_set_pars.size(); k++) {
            int mode_k = pot_set_pars[k].mode;
            assert(mode_k>=0||mode_k<=2);
            if (mode_k==2) {
                double* vel_k = pot_set_pars[k].vel;
                double* pos_k = pot_set_pars[k].pos;
                pos_k[0] += vel_k[0]*dt_scale;
                pos_k[1] += vel_k[1]*dt_scale;
                pos_k[2] += vel_k[2]*dt_scale;
            }        
        }
    }


    //! calculate acceleration for moving potentials from other potentials
    /*!
      @param[in] _time: current time [input unit] of particle system
      @param[in] pos_pcm: position of particle system in the galactic frame [input unit] for mode 1 potential
     */
    void calcMovePotAccFromPot(const double _time, const double* pos_pcm) {
        assert(pot_sets.size()==pot_set_pars.size());
        double t = _time*tscale;

        for (size_t k=0; k<pot_set_pars.size(); k++) {
            int mode_k = pot_set_pars[k].mode;
            assert(mode_k>=0||mode_k<=2);
            if (mode_k==2) {
                double* pos_k = pot_set_pars[k].pos; // galactic frame
                double* acc_pot = pot_set_pars[k].acc;
                for (size_t j=0; j<pot_set_pars.size(); j++) {
                    if (k==j) continue;
                    int mode_j = pot_set_pars[j].mode;
                    double* pos_j = pot_set_pars[j].pos;
                    double dx,dy,dz;
                    if (mode_j==1) {
                        // should be consistent in galactic frame
                        dx = pos_k[0]-pos_j[0]-pos_pcm[0];
                        dy = pos_k[1]-pos_j[1]-pos_pcm[1];
                        dz = pos_k[2]-pos_j[2]-pos_pcm[2];
                    }
                    else {
                        dx = pos_k[0]-pos_j[0];
                        dy = pos_k[1]-pos_j[1];
                        dz = pos_k[2]-pos_j[2];
                    }
                    double rxy= std::sqrt(dx*dx+dy*dy);
                    double phi= std::acos(dx/rxy);
                    double sinphi = dy/rxy;
                    double cosphi = dx/rxy;
                    auto& pot_args = pot_sets[j].arguments;
                    int npot = pot_sets[j].npot;
                    double acc_rxy = calcRforce(rxy, dz, phi, t, npot, pot_args);
                    double acc_z   = calczforce(rxy, dz, phi, t, npot, pot_args);
#if (defined GALPY_VERSION_1_7_9) || (defined GALPY_VERSION_1_7_1)
                    double acc_phi = calcPhiforce(rxy, dz, phi, t, npot, pot_args);
#else
                    double acc_phi = calcphitorque(rxy, dz, phi, t, npot, pot_args);
#endif
                    if (rxy>0.0) {
                        acc_pot[0] += (cosphi*acc_rxy - sinphi*acc_phi/rxy);
                        acc_pot[1] += (sinphi*acc_rxy + cosphi*acc_phi/rxy);
                        acc_pot[2] += acc_z;
                    }
                }
            }
        }
    }


    //! calculate acceleration and potential at give position
    /*!
      @param[out] acc: [3] acceleration to return
      @param[out] pot: potential to return 
      @param[in] _time: time in input unit
      @param[in] gm: G*mass of particles [input unit]
      @param[in] pos_g: position of particles in the galactic frame [input unit]
      @param[in] pos_l: position of particles in the particle system frame [input unit]
     */
    void calcAccPot(double* acc, double& pot, const double _time, const double gm, const double* pos_g, const double* pos_l) {
        assert(pot_sets.size()==pot_set_pars.size());
        int nset = pot_sets.size();
        if (nset>0) {
            double t = _time*tscale;

            // galactic frame and rest frame of particle system
            double x[2] = {pos_g[0]*rscale, pos_l[0]*rscale};
            double y[2] = {pos_g[1]*rscale, pos_l[1]*rscale};
            double z[2] = {pos_g[2]*rscale, pos_l[2]*rscale};

            pot = 0;
            acc[0] = acc[1] = acc[2] = 0.0;



            for (int k=0; k<nset; k++) {
                int mode_k = pot_set_pars[k].mode;
                assert(mode_k>=0||mode_k<=2);
                int npot = pot_sets[k].npot;
                double* pos_k = pot_set_pars[k].pos;
                int i = (mode_k & 1); // get first bit to select frame (0: galactic; 1: rest)
                // frame is consistent 
                double dx = x[i]-pos_k[0];
                double dy = y[i]-pos_k[1];
                double dz = z[i]-pos_k[2];
                double rxy= std::sqrt(dx*dx+dy*dy);
                double phi= std::acos(dx/rxy);
                double sinphi = dy/rxy;
                double cosphi = dx/rxy;

                auto& pot_args = pot_sets[k].arguments;
                double acc_rxy = calcRforce(rxy, dz, phi, t, npot, pot_args);
                double acc_z   = calczforce(rxy, dz, phi, t, npot, pot_args);
#if (defined GALPY_VERSION_1_7_9) || (defined GALPY_VERSION_1_7_1)
                double acc_phi = calcPhiforce(rxy, dz, phi, t, npot, pot_args);
#else
                double acc_phi = calcphitorque(rxy, dz, phi, t, npot, pot_args);
#endif
                double pot_i = evaluatePotentials(rxy, dz, npot, pot_args);
                double gm_pot = pot_set_pars[k].gm;
                if (rxy>0.0) {
                    assert(!std::isinf(acc_rxy));
                    assert(!std::isnan(acc_rxy));
                    assert(!std::isinf(acc_phi));
                    assert(!std::isnan(acc_phi));
                    assert(!std::isinf(pot));
                    assert(!std::isnan(pot));
                    pot += pot_i;
                    double acc_x = (cosphi*acc_rxy - sinphi*acc_phi/rxy);
                    double acc_y = (sinphi*acc_rxy + cosphi*acc_phi/rxy);
                    acc[0] += acc_x;
                    acc[1] += acc_y;
                    acc[2] += acc_z;
                    if (mode_k==2) {
                        double* acc_pot = pot_set_pars[k].acc;
                        acc_pot[0] -= gm*acc_x/gm_pot; // anti-acceleration to potential set origin
                        acc_pot[1] -= gm*acc_y/gm_pot;
                        acc_pot[2] -= gm*acc_z/gm_pot;
                    }
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

    //! calculate density for given potential set
    /*!
      @param[in] _iset: potential set index
      @param[in] _time: time in input unit
      @param[in] pos_g: position of particles in the galactic frame [input unit]
      @param[in] pos_l: position of particles in the particle system frame [input unit]
     */
    double calcSetDensity(const int _iset, const double _time, double* pos_g, const double* pos_l) {
        assert(pot_sets.size()==pot_set_pars.size());
        assert(_iset<int(pot_sets.size()));
        double t = _time*tscale;

        // galactic frame and rest frame of particle system
        double x[2] = {pos_g[0]*rscale, pos_l[0]*rscale};
        double y[2] = {pos_g[1]*rscale, pos_l[1]*rscale};
        double z[2] = {pos_g[2]*rscale, pos_l[2]*rscale};
        
        int mode_k = pot_set_pars[_iset].mode;
        assert(mode_k>=0||mode_k<=2);
        int npot = pot_sets[_iset].npot;
        double* pos_k = pot_set_pars[_iset].pos;
        int i = (mode_k & 1); // get first bit to select frame (0: galactic; 1: rest)        
        // frame is consistent 
        double dx = x[i]-pos_k[0];
        double dy = y[i]-pos_k[1];
        double dz = z[i]-pos_k[2];
        double rxy= std::sqrt(dx*dx+dy*dy);
        double phi= std::acos(dx/rxy);

        auto& pot_args = pot_sets[_iset].arguments;
        double density = calcDensity(rxy, dz, phi, t, npot, pot_args);
        return density;
    }

    //! get number of set
    int getNSet() const {
        return pot_set_pars.size();
    }

    //! write data for restart
    /*! 
      Write data sctructure:
      Time N_set Mode[:]  # N_set: set mode
          After first line, the format is (size, array)
      pot_set_par:        # mode, gm, pos, vel per line
      pot_type_offset     # N_set+1: number of pots per set 
      pot_type            # N_pot: pot type list
      pot_args_offset     # N_set+1: number of arguments per set
      pot_args            # N_pot_arg pot argument list
      change_args_offset  # N_set+1: number of changing arguments per set
      change_args         # N_change: changing arg index, mode, rate 

      @param[in] _filename: file to save data
    */
    void writeDataToFile(const std::string& _filename) {
        std::ofstream fout;
        fout.open(_filename.c_str(), std::ifstream::out);
        if (!fout.is_open()) {
            std::cerr<<"Error: Galpy potential parameter file to write, "<<_filename<<", cannot be open!"<<std::endl;
            abort();
        }
        fout<<std::setprecision(14);

        std::size_t n_set = pot_set_pars.size();
        std::size_t n_pot = pot_type.size();
        std::size_t n_pot_offset = pot_type_offset.size();
        std::size_t n_arg = pot_args.size();
        std::size_t n_arg_offset = pot_args_offset.size();
        std::size_t n_change = change_args.size();
        std::size_t n_change_offset = change_args_offset.size();
        fout<<time<<" "
            <<n_set<<" "
            <<std::endl;
        for (std::size_t i=0; i<n_set; i++) {
            pot_set_pars[i].writeData(fout);
            fout<<std::endl;
        }
        fout<<std::endl<<n_pot_offset<<" ";
        for (std::size_t i=0; i<n_pot_offset; i++) {
            fout<<pot_type_offset[i]<<" ";
        }
        fout<<std::endl<<n_pot<<" ";
        for (std::size_t i=0; i<n_pot; i++) {
            fout<<pot_type[i]<<" ";
        }
        fout<<std::endl<<n_arg_offset<<" ";
        for (std::size_t i=0; i<n_arg_offset; i++) {
            fout<<pot_args_offset[i]<<" ";
        }
        fout<<std::endl<<n_arg<<" ";
        for (std::size_t i=0; i<n_arg; i++) {
            fout<<pot_args[i]<<" ";
        }
        fout<<std::endl<<n_change_offset<<" ";
        for (std::size_t i=0; i<n_change_offset; i++) {
            fout<<change_args_offset[i]<<" ";
        }
        fout<<std::endl<<n_change<<" ";
        for (std::size_t i=0; i<n_change; i++) {
            change_args[i].writeData(fout);
        }
        fout<<std::endl;

        fout.close();
    }

    //! read configure data for restart
    /*
      @param[in] _filename: file to save data
      @param[in] _print_flag: if true, print reading parameters
    */
    void readDataFromFile(const std::string& _filename, const bool _print_flag) {

        std::size_t n_set, n_pot, n_pot_offset, n_arg, n_arg_offset, n_change, n_change_offset;

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
        int my_rank = PS::Comm::getRank();

        if (my_rank==0) {
#endif
            std::ifstream fin;
            fin.open(_filename.c_str(), std::ifstream::in);
            if (!fin.is_open()) {
                std::cerr<<"Error: Galpy configure file for restart "<<_filename.c_str()<<" cannot be open!"<<std::endl;
                abort();
            }

            fin>>time>>n_set;
            eofCheck(fin, "time and number of potential sets");
            geZeroCheck(n_set, "number of potential sets");
            pot_set_pars.resize(n_set);
            for (std::size_t i=0; i<n_set; i++) {
                pot_set_pars[i].readData(fin);
                eofCheck(fin, "potential set parameter");
            }
            fin>>n_pot_offset;
            assert(n_pot_offset==n_set+1);
            pot_type_offset.resize(n_pot_offset);
            for (std::size_t i=0; i<n_pot_offset; i++) {
                fin>>pot_type_offset[i];
            }
            eofCheck(fin, "potential type array offset");

            fin>>n_pot;
            geZeroCheck(n_pot, "number of potential types");
            pot_type.resize(n_pot);
            for (std::size_t i=0; i<n_pot; i++) {
                fin>>pot_type[i];
            }
            eofCheck(fin, "potential type");

            fin>>n_arg_offset;
            pot_args_offset.resize(n_arg_offset);
            assert(n_arg_offset==n_set+1);
            for (std::size_t i=0; i<n_arg_offset; i++) {
                fin>>pot_args_offset[i];
            }
            eofCheck(fin, "potential argument array offset");

            fin>>n_arg;
            geZeroCheck(n_arg, "number of potential arguments");
            pot_args.resize(n_arg);
            for (std::size_t i=0; i<n_arg; i++) {
                fin>>pot_args[i];
            }
            eofCheck(fin, "potential arguments");

            fin>>n_change_offset;
            change_args_offset.resize(n_change_offset);
            assert(n_change_offset==n_set+1);
            for (std::size_t i=0; i<n_change_offset; i++) {
                fin>>change_args_offset[i];
            }
            eofCheck(fin, "changing argument offset");

            fin>>n_change;
            geZeroCheck(n_change, "number of changing arguments");
            change_args.resize(n_change);
            for (std::size_t i=0; i<n_change; i++) {
                change_args[i].readData(fin);
            }
            eofCheck(fin, "changing arguments");

            fin.close();

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
        }

        broadcastDataMPI();
#endif

        if (_print_flag) printData(std::cout);
    }

    //! calculate new arguments for changing potentials
    /*! 
      Change index counts from 0 inside the local potential set
      If change index is -1, change the gm of potential set instead.
      To find the correct potential argument, the pot_args_offset is needed.

      @param[in] _time: new time (Galpy Unit)
      @param[in] _print_flag: if true, print updated argument
      \return if changed, return true
    */
    bool evolveChangingArguments(const double& _time, const bool _print_flag) {
        if (change_args.size()>0) {
            double dt = _time - time;
            if (dt>0) {
                if (_print_flag) std::cout<<"Galpy update argument, time [Galpy unit]: "<<_time<<" dt: "<<dt<<std::endl;
                std::size_t n_set = change_args_offset.size()-1;
                for (std::size_t i=0; i<n_set; i++) {
                    for (int k=change_args_offset[i]; k<change_args_offset[i+1]; k++) {
                        int index_local = change_args[k].index;
                        double* change_arg_ptr = NULL;
                        if (index_local==-1)  // change gm instead of arguments
                            change_arg_ptr = &(pot_set_pars[i].gm);
                        else {
                            int index = index_local + pot_args_offset[i]; // change args index counting for individual set, need to obtain the global index of pot_args.
                            change_arg_ptr = &(pot_args[index]);
                        }
                        if (change_args[k].mode==1) {
                            *change_arg_ptr += change_args[k].rate*dt;
                            if (_print_flag) std::cout<<"Set "<<i+1<<" Index: "<<change_args[k].index<<" "<<" new_argument(linear): "<<*change_arg_ptr;
                        }
                        else if(change_args[k].mode==2) {
                            *change_arg_ptr *= std::exp(change_args[k].rate*dt);
                            if (_print_flag) std::cout<<"Set "<<i+1<<" Index: "<<change_args[k].index<<" "<<" new_argument(exp): "<<*change_arg_ptr;
                        }
                        if (_print_flag) std::cout<<std::endl;
                    }
                }
                time = _time;
                return true;
            }
        }
        time = _time;
        return false;
    }

    void freePotentialArgs() {
        if (!pot_sets.empty()) {
            for (size_t i=0; i<pot_sets.size(); i++) pot_sets[i].clear();
            pot_sets.resize(0);
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


