#pragma once

#ifdef P3T_64BIT
#define CALC_EP_64bit
#define CALC_SP_64bit
#define RSQRT_NR_EPJ_X4
#define RSQRT_NR_SPJ_X4

#elif P3T_MIXBIT
#define CALC_EP_64bit
#define RSQRT_NR_EPJ_X4

#else
#define RSQRT_NR_EPJ_X2
//#define RSQRT_NR_SPJ_X2
#endif 

#if defined(INTRINSIC_K) || defined(INTRINSIC_X86)
#define INTRINSIC
#endif

#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<sstream>
//#include<unistd.h>
#include<getopt.h>

#ifdef MPI_DEBUG
#include <mpi.h>
// Send recv debug
int MPI_Isend(void* buffer, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request* req)
{
   int ret;

   int size;
   MPI_Type_size(datatype, &size);
   std::cerr<<"MPI_ISend count "<<count<<" dest "<<dest<<" datatype "<<size<<" tag "<<tag<<std::endl;

   ret = PMPI_Isend(buffer, count, datatype, dest, tag, comm, req);

   return ret;
}

// Debug
int MPI_Irecv(void* buffer, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request* req)
{
   int ret;

   int size;
   MPI_Type_size(datatype, &size);
   std::cerr<<"MPI_IRecv count "<<count<<" dest "<<dest<<" datatype "<<size<<" tag "<<tag<<std::endl;

   ret = PMPI_Irecv(buffer, count, datatype, dest, tag, comm, req);

   return ret;
}
#endif

#include<particle_simulator.hpp>
#include"hard_assert.hpp"
#include"soft_ptcl.hpp"
#include"soft_force.hpp"
#ifdef USE_GPU
#include"force_gpu_cuda.hpp"
#endif
#include"energy.hpp"
#include"hard.hpp"
#include"io.hpp"
#include"status.hpp"
#include"particle_distribution_generator.hpp"
#include"domain.hpp"
#include"cluster_list.hpp"
#include"kickdriftstep.hpp"
#ifdef PROFILE
#include"profile.hpp"
#endif
#include"static_variables.hpp"

//! IO parameters for Petar
class IOParamsPeTar{
public:
    // IO parameters
    IOParamsContainer input_par_store;
    IOParams<PS::F64> ratio_r_cut;
    IOParams<PS::F64> theta;
    IOParams<PS::S32> n_leaf_limit;
#ifdef USE__AVX512
    IOParams<PS::S32> n_group_limit;
#else
    IOParams<PS::S32> n_group_limit;
#endif
    IOParams<PS::S32> n_interrupt_limit;
    IOParams<PS::S32> n_smp_ave;
    IOParams<PS::S32> n_split;
    IOParams<PS::S64> n_bin;
    IOParams<PS::F64> n_step_per_orbit;
    IOParams<PS::F64> time_end;
    IOParams<PS::F64> eta;
    IOParams<PS::F64> gravitational_constant;
    IOParams<PS::S64> n_glb;
    IOParams<PS::S64> id_offset;
    IOParams<PS::F64> dt_soft;
    IOParams<PS::F64> dt_snp;
    IOParams<PS::F64> search_vel_factor;
    IOParams<PS::F64> search_peri_factor;
    IOParams<PS::F64> dt_limit_hard_factor;
    IOParams<PS::S32> dt_min_hermite_index;
    IOParams<PS::S32> dt_min_arc_index;
    IOParams<PS::F64> dt_err_pert;
    IOParams<PS::F64> dt_err_soft;
    IOParams<PS::F64> e_err_arc;
#ifdef HARD_CHECK_ENERGY
    IOParams<PS::F64> e_err_hard;
#endif
    IOParams<PS::S32> step_limit_arc;
    IOParams<PS::F64> eps;
    IOParams<PS::F64> r_out;
    IOParams<PS::F64> r_bin;
    IOParams<PS::F64> r_search_max;
    IOParams<PS::F64> r_search_min;
    IOParams<PS::F64> sd_factor;
    IOParams<PS::S32> data_format;
    IOParams<PS::S32> write_style;
    IOParams<PS::S32> interrupt_detection_option;
    IOParams<std::string> fname_snp;
    IOParams<std::string> fname_par;
    IOParams<std::string> fname_inp;

    // flag
    bool app_flag; // appending data flag
    bool print_flag; 

    IOParamsPeTar(): input_par_store(), 
                     ratio_r_cut  (input_par_store, 0.1,  "r_in / r_out"),
                     theta        (input_par_store, 0.3,  "Openning angle theta"),
                     n_leaf_limit (input_par_store, 20,   "Tree leaf number limit", "optimized value shoudl be slightly >=11+N_bin_sample (20)"),
#ifdef USE__AVX512
                     n_group_limit(input_par_store, 1024, "Tree group number limit", "optimized for x86-AVX512 (1024)"),    
#else
                     n_group_limit(input_par_store, 512,  "Tree group number limit", "optimized for x86-AVX2 (512)"),
#endif
                     n_interrupt_limit(input_par_store, 512, "Interrupt binary number limit"),
                     n_smp_ave    (input_par_store, 100,  "Average target number of sample particles per process"),
                     n_split      (input_par_store, 8,    "Number of binary sample points for tree perturbation force"),
                     n_bin        (input_par_store, 0,    "Number of binaries used for initialization (assume binaries ID=1,2*n_bin)"),
                     n_step_per_orbit(input_par_store, 8, "Number of steps per slow-down binary orbits (binary period/tree timestep) for isolated binaries; also the maximum criterion for switching on tidal tensor method"),
                     time_end     (input_par_store, 10.0, "Finishing time"),
                     eta          (input_par_store, 0.1,  "Hermite time step coefficient eta"),
                     gravitational_constant(input_par_store, 1.0,  "Gravitational constant"),
                     n_glb        (input_par_store, 0,    "Total number of particles, only used to generate particles if needed"),
                     id_offset    (input_par_store, -1,   "Starting id for artificial particles, total number of real particles must be always smaller than this","n_glb+1"),
                     dt_soft      (input_par_store, 0.0,  "Tree timestep","0.1*r_out/sigma_1D"),
                     dt_snp       (input_par_store, 1.0,  "Output time interval of particle dataset"),
                     search_vel_factor (input_par_store, 3.0,  "Neighbor searching coefficient for velocity check (v*dt)"),
                     search_peri_factor(input_par_store, 1.5,  "Neighbor searching coefficient for peri-center check"),
                     dt_limit_hard_factor(input_par_store, 4.0,  "Limit of tree time step/hard time step"),
                     dt_min_hermite_index(input_par_store, 40,   "Power index n for the smallest time step (0.5^n) allowed in Hermite integrator"),
                     dt_min_arc_index    (input_par_store, 64,   "Power index n for the smallest time step (0.5^n) allowed in ARC integrator, suppressed"),
                     dt_err_pert  (input_par_store, 1e-6, "Time synchronization maximum (relative) error for perturbed ARC integrator, suppressed"),
                     dt_err_soft  (input_par_store, 1e-3, "Time synchronization maximum (relative) error for no-perturber (only soft perturbation) ARC integrator, suppressed"),
                     e_err_arc    (input_par_store, 1e-8,"Maximum energy error allown for ARC integrator"),
#ifdef HARD_CHECK_ENERGY
                     e_err_hard   (input_par_store, 1e-4, "Maximum energy error allown for hard integrator"),
#endif
                     step_limit_arc(input_par_store, 1000000, "Maximum step allown for ARC sym integrator"),
                     eps          (input_par_store, 0.0,  "Softerning eps"),
                     r_out        (input_par_store, 0.0,  "Transit function outer boundary radius", "<m>/sigma_1D^2/ratio_r_cut"),
                     r_bin        (input_par_store, 0.0,  "Tidal tensor box size and binary radius criterion", "theta*r_in"),
                     r_search_max (input_par_store, 0.0,  "Maximum search radius criterion", "5*r_out"),
                     r_search_min (input_par_store, 0.0,  "Minimum search radius  value","auto"),
                     sd_factor    (input_par_store, 1e-4, "Slowdown perturbation criterion"),
                     data_format  (input_par_store, 1,    "Data read(r)/write(w) format BINARY(B)/ASCII(A): r-B/w-A (3), r-A/w-B (2), rw-A (1), rw-B (0)"),
                     write_style  (input_par_store, 1,    "File Writing style: 0, no output; 1. write snapshots, status and profile separately; 2. write snapshot and status in one line per step (no MPI support); 3. write only status and profile"),
                     interrupt_detection_option(input_par_store, 0, "detect interruption", "no"),
                     fname_snp(input_par_store, "data","Prefix filename of dataset: [prefix].[File ID]"),
                     fname_par(input_par_store, "input.par", "Input parameter file (this option should be used first before any other options)"),
                     fname_inp(input_par_store, "", "Input data file"),
                     app_flag(false), print_flag(false) {}

    
    //! reading parameters from GNU option API
    /*!
      @param[in] argc: number of options
      @param[in] argv: string of options
      @param[in] print_flag: true: print input
      \return -1 if help is used
     */
    int read(int argc, char *argv[]) {
        static struct option long_options[] = {
            {"number-split", required_argument, 0, 0},        
            {"search-vel-factor", required_argument, 0, 1},  
            {"dt-max-factor", required_argument, 0, 2},  
            {"dt-min-hermite", required_argument, 0, 3}, 
            {"number-group-limit", required_argument, 0, 4},
            {"number-leaf-limit", required_argument, 0, 5},
            {"number-sample-average", required_argument, 0, 6},
            {"energy-err-arc", required_argument, 0, 7}, 
            {"soft-eps", required_argument, 0, 8},       
            {"slowdown-factor", required_argument, 0, 9},
            {"r-ratio", required_argument, 0, 10},       
            {"r-bin",   required_argument, 0, 11},       
            {"search-peri-factor", required_argument, 0, 12}, 
            {"hermite-eta", required_argument, 0, 13}, 
            {"help",no_argument, 0, 'h'},                
#ifdef HARD_CHECK_ENERGY
            {"energy-err-hard", required_argument, 0, 14},  
#endif
            {"step-limit-arc", required_argument, 0, 15},   
            {"disable-print-info", no_argument, 0, 16},
            {"number-interrupt-limt",required_argument, 0, 17},
            {"detect-interrupt", no_argument, 0, 18},
            {"number-step-tt", required_argument, 0, 19},
            {0,0,0,0}
        };

        int copt;
        int option_index;
        int n_opt=0;
        while ((copt = getopt_long(argc, argv, "i:at:s:o:r:b:n:G:T:f:p:w:h", long_options, &option_index)) != -1) 
            switch (copt) {
            case 0:
                n_split.value = atoi(optarg);
                if(print_flag) n_split.print(std::cout);
                assert(n_split.value>=4);
                n_opt+=2;
                break;
            case 1:
                search_vel_factor.value = atof(optarg);
                if(print_flag) search_vel_factor.print(std::cout);
                assert(search_vel_factor.value>0.0);
                n_opt+=2;
                break;
            case 2:
                dt_limit_hard_factor.value = atof(optarg);
                if(print_flag) dt_limit_hard_factor.print(std::cout);
                assert(dt_limit_hard_factor.value > 0.0);
                n_opt+=2;
                break;
            case 3:
                dt_min_hermite_index.value = atoi(optarg);
                if(print_flag) dt_min_hermite_index.print(std::cout);
                assert(dt_min_hermite_index.value > 0);
                n_opt+=2;
                break;
            case 4:
                n_group_limit.value = atoi(optarg);
                if(print_flag) n_group_limit.print(std::cout);
                assert(n_group_limit.value>0);
                n_opt+=2;
                break;
            case 5:
                n_leaf_limit.value = atoi(optarg);
                if(print_flag) n_leaf_limit.print(std::cout);
                assert(n_leaf_limit.value>0);
                n_opt+=2;
                break;
            case 6:
                n_smp_ave.value = atoi(optarg);
                if(print_flag) n_smp_ave.print(std::cout);
                assert(n_smp_ave.value>0.0);
                n_opt+=2;
                break;
            case 7:
                e_err_arc.value = atof(optarg);
                if(print_flag) e_err_arc.print(std::cout);
                assert(e_err_arc.value > 0.0);
                n_opt+=2;
                break;
            case 8:
                eps.value = atof(optarg);
                if(print_flag) eps.print(std::cout);
                assert(eps.value>=0.0);
                n_opt+=2;
                break;
            case 9:
                sd_factor.value = atof(optarg);
                if(print_flag) sd_factor.print(std::cout);
                assert(sd_factor.value>0.0);
                n_opt+=2;
                break;
            case 10:
                ratio_r_cut.value = atof(optarg);
                if(print_flag) ratio_r_cut.print(std::cout);
                assert(ratio_r_cut.value>0.0);
                assert(ratio_r_cut.value<1.0);
                n_opt+=2;
                break;
            case 11:
                r_bin.value = atof(optarg);
                if(print_flag) r_bin.print(std::cout);
                assert(r_bin.value>0.0);
                n_opt+=2;
                break;
            case 12:
                search_peri_factor.value = atof(optarg);
                if(print_flag) search_peri_factor.print(std::cout);
                assert(search_peri_factor.value>=1.0);
                n_opt+=2;
                break;
            case 13:
                eta.value = atof(optarg);
                if(print_flag) eta.print(std::cout);
                assert(eta.value>0.0);
                n_opt+=2;
                break;
#ifdef HARD_CHECK_ENERGY
            case 14:
                e_err_hard.value = atof(optarg);
                if(print_flag) e_err_hard.print(std::cout);
                n_opt+=2;
                break;
#endif
            case 15:
                step_limit_arc.value = atoi(optarg);
                if(print_flag) step_limit_arc.print(std::cout);
                n_opt+=2;
                break;
            case 16:
                print_flag = false;
                n_opt++;
                break;
            case 17:
                n_interrupt_limit.value = atoi(optarg);
                if(print_flag) n_interrupt_limit.print(std::cout);
                assert(n_interrupt_limit.value>0);
                n_opt+=2;
                break;
            case 18:
                interrupt_detection_option.value = 1;
                if(print_flag) interrupt_detection_option.print(std::cout);
                n_opt++;
                break;
            case 19:
                n_step_per_orbit.value = atof(optarg);
                if(print_flag) n_step_per_orbit.print(std::cout);
                assert(n_step_per_orbit.value>=1.0);
                n_opt+=2;
                break;
            case 'i':
                data_format.value = atoi(optarg);
                if(print_flag) data_format.print(std::cout);
                assert(data_format.value>=0&&data_format.value<=3);
                n_opt+=2;
                break;
            case 'a':
                app_flag=true;
                n_opt++;
                break;
            case 't':
                time_end.value = atof(optarg);
                if(print_flag) time_end.print(std::cout);
                assert(time_end.value>=0.0);
                n_opt+=2;
                break;
            case 's':
                dt_soft.value = atof(optarg);
                if(print_flag) dt_soft.print(std::cout);
                assert(dt_soft.value>0.0);
                n_opt+=2;
                break;
            case 'o':
                dt_snp.value = atof(optarg);
                if(print_flag) dt_snp.print(std::cout);
                assert(dt_snp.value>0.0);
                n_opt+=2;
                break;
            case 'r':
                r_out.value = atof(optarg);
                if(print_flag) r_out.print(std::cout);
                assert(r_out.value>0.0);
                n_opt+=2;
                break;
            case 'b':
                n_bin.value = atoi(optarg);
                if(print_flag) n_bin.print(std::cout);
                assert(n_bin.value>=0);
                n_opt+=2;
                break;
            case 'n':
                n_glb.value = atol(optarg);
                if(print_flag) n_glb.print(std::cout);
                assert(n_glb.value>0);
                n_opt+=2;
                break;
            case 'G':
                gravitational_constant.value = atof(optarg);
                if(print_flag) gravitational_constant.print(std::cout);
                assert(gravitational_constant.value>0.0);
                n_opt+=2;
                break;
            case 'T':
                theta.value = atof(optarg);
                if(print_flag) theta.print(std::cout);
                assert(theta.value>=0.0);
                n_opt+=2;
                break;
            case 'f':
                fname_snp.value = optarg;
                if(print_flag) fname_snp.print(std::cout);
                n_opt+=2;
                break;
            case 'p':
                fname_par.value = optarg;
                if(print_flag) {
                    fname_par.print(std::cout);
                    FILE* fpar_in;
                    if( (fpar_in = fopen(fname_par.value.c_str(),"r")) == NULL) {
                        fprintf(stderr,"Error: Cannot open file %s.\n", fname_par.value.c_str());
                        abort();
                    }
                    input_par_store.readAscii(fpar_in);
                    fclose(fpar_in);
                }
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
                input_par_store.mpi_broadcast();
                PS::Comm::barrier();
#endif
                n_opt+=2;
                break;
            case 'w':
                write_style.value = atoi(optarg);
                if(print_flag) write_style.print(std::cout);
                n_opt+=2;
                break;
            case 'h':
                if(print_flag){
                    std::cout<<"Usage: nbody.out [option] [filename]"<<std::endl;
                    std::cout<<"       Option defaulted values are shown after ':'\n"<<std::endl;
                    std::cout<<"  -i: [I] "<<data_format<<std::endl;
                    std::cout<<"          File content:\n"
                             <<"            First line: \n"
                             <<"             1. File_ID: 0 for initialization, else for restarting\n"
                             <<"             2. N_particle \n"
                             <<"             3. Time\n"
                             <<"            Following lines:\n";
                    std::cout<<"             ";
                    FPSoft::printColumnTitle(std::cout,10);
                    std::cout<<std::endl;
/*
                             <<"             1. mass\n"
                             <<"             2. position[3]\n"
                             <<"             5. velocity[3]\n"
                             <<"             8. r_search (0.0)\n"
                             <<"             9. mass_backup (0.0)\n"
                             <<"             10. ID (>0,unique)\n"
                             <<"             11. status (0)\n"
                             <<"             12. r_in (0.0)\n"
                             <<"             13. r_out (0.0)\n"
                             <<"             14. Acceleration[3] (0.0)\n"
                             <<"             17. Potential (0.0)\n"
                             <<"             18. N_neighbor (0)\n"
                             <<"             in () show initialization values which should be used together with FILE_ID = 0"<<std::endl;
*/
                    std::cout<<"  -a:     data output style (except snapshot) becomes appending, defaulted: replace"<<std::endl;
                    std::cout<<"  -t: [F] "<<time_end<<std::endl;
                    std::cout<<"  -s: [F] "<<dt_soft<<std::endl;
                    std::cout<<"  -o: [F] "<<dt_snp<<std::endl;
                    std::cout<<"        --detect-interrupt:      "<<interrupt_detection_option<<std::endl;
                    std::cout<<"        --dt-max-factor:     [F] "<<dt_limit_hard_factor<<std::endl;
                    std::cout<<"        --dt-min-hermite:    [I] "<<dt_min_hermite_index<<std::endl;
                    std::cout<<"        --disable-print-info:  "<<"Do not print information"<<std::endl;
                    std::cout<<"        --disable-write-info:  "<<"Do not write information"<<std::endl;
                    std::cout<<"  -r: [F] "<<r_out<<std::endl;
                    std::cout<<"        --r-ratio:           [F] "<<ratio_r_cut<<std::endl;
                    std::cout<<"        --r-bin:             [F] "<<r_bin<<std::endl;
                    std::cout<<"  -b: [I] "<<n_bin<<std::endl;
                    std::cout<<"  -n: [I] "<<n_glb<<std::endl;
                    std::cout<<"        --number-split:           [I] "<<n_split<<std::endl;
                    std::cout<<"        --number-group-limit:     [I] "<<n_group_limit<<std::endl;
                    std::cout<<"        --number-leaf-limit:      [I] "<<n_leaf_limit<<std::endl;
                    std::cout<<"        --number-interrupt-limit: [I] "<<n_interrupt_limit<<std::endl;
                    std::cout<<"        --number-sample-average:  [I] "<<n_smp_ave<<std::endl;
                    std::cout<<"        --number-step-tt:         [F] "<<n_step_per_orbit<<std::endl;
                    std::cout<<"  -G: [F] "<<gravitational_constant<<std::endl;
                    std::cout<<"  -T: [F] "<<theta<<std::endl;
                    std::cout<<"        --hermite-eta:       [F] "<<eta<<std::endl;
                    std::cout<<"        --search-vel-factor: [F] "<<search_vel_factor<<std::endl;
                    std::cout<<"        --search-peri-factor:[F] "<<search_peri_factor<<std::endl;
                    std::cout<<"        --energy-err-arc:    [F] "<<e_err_arc<<std::endl;
#ifdef HARD_CHECK_ENERGY
                    std::cout<<"        --energy-err-hard:   [F] "<<e_err_hard<<std::endl;
#endif
                    std::cout<<"        --step-limit-arc:    [F] "<<step_limit_arc<<std::endl;
                    std::cout<<"        --slowdown-factor:   [F] "<<sd_factor<<std::endl;
                    std::cout<<"        --soft-eps:          [F] "<<eps<<std::endl;
                    std::cout<<"  -f: [S] "<<fname_snp<<std::endl;
                    std::cout<<"  -p: [S] "<<fname_par<<std::endl;
                    std::cout<<"  -w: [I] "<<write_style<<std::endl;
                    std::cout<<"  -h(--help):               print help"<<std::endl;
                    std::cout<<"*** PS: r_in : transit function inner boundary radius\n"
                             <<"        r_out: transit function outer boundary radius\n"
                             <<"        sigma: half-mass radius velocity dispersion\n"
                             <<"        n_bin: number of primordial binaries\n"
                             <<"        <m>  : averaged mass"<<std::endl;
                }
                return -1;
            }
        
        if (argc-n_opt>1) {
            fname_inp.value =argv[argc-1];
            if(print_flag) std::cout<<"Reading data file name: "<<fname_inp.value<<std::endl;
        }        

        if(print_flag) std::cout<<"----- Finish reading input options -----\n";

        return 0;
    }

    //! check paramters
    bool checkParams() {
        assert(n_split.value>=4);
        assert(search_vel_factor.value>0.0);
        assert(dt_limit_hard_factor.value > 0.0);
        assert(dt_min_hermite_index.value > 0);
        assert(e_err_arc.value > 0.0);
        assert(eps.value>=0.0);
        assert(sd_factor.value>0.0);
        assert(ratio_r_cut.value>0.0);
        assert(ratio_r_cut.value<1.0);
        assert(r_bin.value>0.0);
        assert(search_peri_factor.value>=1.0);
        assert(data_format.value>=0||data_format.value<=3);
        assert(time_end.value>=0.0);
        assert(dt_soft.value>0.0);
        assert(dt_snp.value>0.0);
        assert(r_out.value>0.0);
        assert(n_bin.value>=0);
        assert(n_glb.value>0);
        assert(n_group_limit.value>0);
        assert(n_interrupt_limit.value>0);
        assert(n_leaf_limit.value>0);
        assert(n_smp_ave.value>0.0);
        assert(theta.value>=0.0);
        assert(eta.value>0.0);
        return true;
    }

};


//! main control class
class PeTar {
public:
    IOParamsPeTar input_parameters;
#ifdef USE_QUAD
    typedef PS::TreeForForceLong<ForceSoft, EPISoft, EPJSoft>::QuadrupoleWithSymmetrySearch TreeForce; 
#else
    typedef PS::TreeForForceLong<ForceSoft, EPISoft, EPJSoft>::MonopoleWithSymmetrySearch TreeForce;
#endif
    typedef PS::ParticleSystem<FPSoft> SystemSoft;

    // For neighbor searching
    typedef PS::TreeForForceShort<ForceSoft, EPISoft, EPJSoft>::Symmetry TreeNB;

#ifdef PROFILE
    PS::S32 dn_loop;
    SysProfile profile;
    SysCounts  n_count;
    SysCounts  n_count_sum;
    FDPSProfile  tree_soft_profile;
    FDPSProfile  tree_nb_profile;
    std::ofstream fprofile;
#endif

    Status stat;
    std::ofstream fstatus;
    PS::F64 time_kick;

    // file system
    FileHeader file_header;
    SystemSoft system_soft;

    // particle index map
    std::map<PS::S64, PS::S32> id_adr_map;

    // domain
    PS::S64 n_loop; // count for domain decomposition
    PS::F64 domain_decompose_weight;
    PS::DomainInfo dinfo;
    PS::F64ort * pos_domain;

    // tree time step manager
    KickDriftStep dt_manager;

    // tree
    TreeNB tree_nb;
    TreeForce tree_soft;

    // hard integrator
    HardManager hard_manager;
    SystemHard system_hard_one_cluster;
    SystemHard system_hard_isolated;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    SystemHard system_hard_connected;
#endif
    int n_interrupt_glb;

    // remove list
    PS::ReallocatableArray<PS::S32> remove_list;

    // search cluster
    SearchCluster search_cluster;

    // MPI 
    PS::S32 my_rank;
    PS::S32 n_proc;

    // safety check flag
    bool initial_fdps_flag;
    bool read_parameters_flag;
    bool read_data_flag;
    bool initial_parameters_flag;
    bool initial_step_flag;

    // initialization
    PeTar(): 
#ifdef PROFILE
        // profile
        dn_loop(0), profile(), n_count(), n_count_sum(), tree_soft_profile(), fprofile(), 
#endif
        stat(), fstatus(), time_kick(0.0),
        file_header(), system_soft(), id_adr_map(),
        n_loop(0), domain_decompose_weight(1.0), dinfo(), pos_domain(NULL), 
        dt_manager(),
        tree_nb(), tree_soft(), 
        hard_manager(), system_hard_one_cluster(), system_hard_isolated(), 
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        system_hard_connected(), 
#endif
        n_interrupt_glb(0),
        remove_list(),
        search_cluster(),
        initial_fdps_flag(false), read_parameters_flag(false), read_data_flag(false), initial_parameters_flag(false), initial_step_flag(false) {
        // set print format
        std::cout<<std::setprecision(PRINT_PRECISION);
        std::cerr<<std::setprecision(PRINT_PRECISION);
     }


private:

    //!initializaton of system parameters
    /*! Obtain Radius parameters, consistent with the input help information
      @param[in]     _tsys:   particle system (soft)
      @param[in,out] _r_in:   changeover function inner boundary
      @param[in,out] _r_out:  changeover function outer boundary
      @param[in,out] _r_bin:  arc group radius criterion
      @param[in,out] _r_search_min: minimum searching radius
      @param[in,out] _r_search_max: maximum searching radius
      @param[out] _v_max:  maximum velocity to calculate r_search
      @param[out] _m_average: averaged mass of particles
      @param[out] _m_max: maximum mass of particles
      @param[in,out] _dt_tree: tree time step
      @param[out] _vel_disp: system velocity dispersion 
      @param[in]     _search_vel_factor: coefficient to calculate r_search
      @param[in]     _ratio_r_cut: _r_out/_r_in
      @param[in]     _n_bin: number of binaries
    */
    void getInitPar(const SystemSoft & _tsys,
                    PS::F64 &_r_in,
                    PS::F64 &_r_out,
                    PS::F64 &_r_bin,
                    PS::F64 &_r_search_min,
                    PS::F64 &_r_search_max,
                    PS::F64 &_v_max,
                    PS::F64 &_m_average,
                    PS::F64 &_m_max,
                    PS::F64 &_dt_tree,
                    PS::F64 &_vel_disp,
                    const PS::F64 _search_vel_factor,
                    const PS::F64 _ratio_r_cut,
                    const PS::S64 _n_bin,
                    const PS::F64 _theta) {

        // local particle number
        const PS::S64 n_loc = _tsys.getNumberOfParticleLocal();

        // local c.m velocity
        PS::F64vec vel_cm_loc = 0.0;
        // local c.m. mass
        PS::F64 mass_cm_loc = 0.0;
        // local maximum mass
        PS::F64 mass_max_loc = 0.0;

        for(PS::S64 i=0; i<n_loc; i++){
            PS::F64 mi = _tsys[i].mass;
            PS::F64vec vi = _tsys[i].vel;

#ifdef PETAR_DEBUG
            assert(mi>0);
#endif
            mass_cm_loc += mi;
            vel_cm_loc += mi * vi;
            mass_max_loc = std::max(mi, mass_max_loc);
        }

        // global c.m. parameters
        PS::F64    mass_cm_glb = PS::Comm::getSum(mass_cm_loc);
        _m_max = PS::Comm::getMaxValue(mass_max_loc);
        PS::F64vec vel_cm_glb  = PS::Comm::getSum(vel_cm_loc);
        vel_cm_glb /= mass_cm_glb;

        // local velocity square
        PS::F64 vel_sq_loc = 0.0;
        PS::S64 n_vel_loc_count = 0;

        // single particle starting index
        PS::S64 single_start_index = 0;
        const PS::S64 bin_last_id = 2*_n_bin;
        if (_tsys[0].id<bin_last_id) {
            single_start_index = std::min(bin_last_id - _tsys[0].id + 1,n_loc);
            if(single_start_index%2!=0) single_start_index--;
        }
        // binary particle starting index
        const PS::S64 binary_start_index = (_tsys[0].id%2==0)?1:0;

        // calculate velocity dispersion
        for (PS::S64 i=binary_start_index; i<single_start_index; i+=2) {
            PS::F64 m1 = _tsys[i].mass;
            PS::F64 m2 = _tsys[i+1].mass;
            PS::F64vec dv = (m1*_tsys[i].vel + m2*_tsys[i+1].vel)/(m1+m2) - vel_cm_glb;
            vel_sq_loc += dv * dv;
            n_vel_loc_count++;
        }
    
        for (PS::S64 i=single_start_index; i<n_loc; i++){
            PS::F64vec dv = _tsys[i].vel - vel_cm_glb;
            vel_sq_loc += dv * dv;
            n_vel_loc_count++;
        }

        const PS::S64    n_vel_glb_count= PS::Comm::getSum(n_vel_loc_count);
        const PS::S64    n_glb          = PS::Comm::getSum(n_loc);
        const PS::F64    vel_sq_glb     = PS::Comm::getSum(vel_sq_loc);
        _vel_disp   = sqrt(vel_sq_glb / 3.0 / (PS::F64)n_vel_glb_count);

        PS::F64 average_mass_glb = mass_cm_glb/(PS::F64)n_glb;
        _m_average = average_mass_glb;

        // flag to check whether r_ous is already defined
        bool r_out_flag = (_r_out>0);
    
        // if r_out is already defined, calculate r_in based on _ratio_r_cut
        if (r_out_flag) _r_in = _r_out * _ratio_r_cut;
        // calculate r_in based on velocity dispersion and averaged mass, calculate r_out by _ratio_r_cut
        else {
            _r_in = 0.5*average_mass_glb / (_vel_disp*_vel_disp);
            _r_out = _r_in / _ratio_r_cut;
        }

        // if tree time step is not defined, calculate tree time step by r_out and velocity dispersion
        if (_dt_tree==0.0) {
            PS::F64 dt_origin = 0.1*_r_out / _vel_disp;
            _dt_tree = 1.0;
            if (dt_origin<1) while (_dt_tree>dt_origin) _dt_tree *= 0.5;
            else {
                while (_dt_tree<=dt_origin) _dt_tree *= 2.0;
                _dt_tree *= 0.5;
            }
        }
        else {
            // regularize dt_tree
            PS::F64 dt_origin = _dt_tree;
            _dt_tree = 1.0;
            if (dt_origin<1) while (_dt_tree>dt_origin) _dt_tree *= 0.5;
            else {
                while (_dt_tree<=dt_origin) _dt_tree *= 2.0;
                _dt_tree *= 0.5;
            }
            // if r_out is not defined, adjust r_out to minimum based on tree step
            if (!r_out_flag) {
                _r_out = 10.0*_dt_tree*_vel_disp;
                _r_in = _r_out*_ratio_r_cut;
            }
        }

        // if r_bin is not defined, set to theta * r_in
        if (_r_bin==0.0) _r_bin = _theta*_r_in;

        // if r_search_min is not defined, calculate by search_vel_factor*velocity_dispersion*tree_time_step + r_out
        if (_r_search_min==0.0) _r_search_min = _search_vel_factor*_vel_disp*_dt_tree + _r_out;
        // if r_search_max is not defined, calcualte by 5*r_out
        if (_r_search_max==0.0) _r_search_max = 5*_r_out;
        // calculate v_max based on r_search_max, tree time step and search_vel_factor
        _v_max = (_r_search_max - _r_out) / _dt_tree / _search_vel_factor;
    }

    //! tree for neighbor searching.
    void treeNeighborSearch() {
#ifdef PROFILE
        profile.tree_nb.start();
        tree_nb.clearNumberOfInteraction();
        tree_nb.clearTimeProfile();
#endif
#ifndef USE_SIMD
        tree_nb.calcForceAllAndWriteBack(SearchNeighborEpEpNoSimd(), system_soft, dinfo);
#else
        tree_nb.calcForceAllAndWriteBack(SearchNeighborEpEpSimd(), system_soft, dinfo);
#endif
        
#ifdef PROFILE
        tree_nb_profile += tree_nb.getTimeProfile();
        profile.tree_nb.barrier();
        PS::Comm::barrier();
        profile.tree_nb.end();
#endif
    }

    //! search clusters
    void searchCluster() {
#ifdef PROFILE
        profile.search_cluster.start();
#endif
        // >2.1 search clusters ----------------------------------------
        search_cluster.searchNeighborOMP<SystemSoft, TreeNB, EPJSoft>
            (system_soft, tree_nb, pos_domain, 1.0, input_parameters.search_peri_factor.value);

        search_cluster.searchClusterLocal();
        search_cluster.setIdClusterLocal();

        // >2.2 Send/receive connect cluster
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
        search_cluster.connectNodes(pos_domain,tree_nb);
        search_cluster.setIdClusterGlobalIteration();
        search_cluster.sendAndRecvCluster(system_soft);
#endif

#ifdef PROFILE
        profile.search_cluster.barrier();
        PS::Comm::barrier();
        profile.search_cluster.end();
#endif    
    }

    //! and ground and create artificial particles
    /*! @param[in] _dt_tree: tree time step for calculating r_search and set stablility checker period limit
     */
    void createGroup(const PS::F64 _dt_tree) {
#ifdef PROFILE
        profile.create_group.start();
#endif

        // >2.3 Find ARC groups and create artificial particles
        // Set local ptcl_hard for isolated  clusters
        system_hard_isolated.setPtclForIsolatedMultiClusterOMP(system_soft, search_cluster.adr_sys_multi_cluster_isolated_, search_cluster.n_ptcl_in_multi_cluster_isolated_);

//#ifdef CLUSTER_DEBUG
//        for (PS::S32 i=0; i<n_loc; i++) system_soft[i].status = -1000000;
//#endif

        // Find groups and add artificial particles to global particle system
        system_hard_isolated.findGroupsAndCreateArtificialParticlesOMP<SystemSoft, FPSoft>(system_soft, _dt_tree);

//#ifdef CLUSTER_DEBUG
//        not correct check, isolated clusters are not reset above
//        for (PS::S32 i=0; i<n_loc; i++) {
//            if (system_soft[i].status ==-1000000) {
//                std::cerr<<"Find un initialized status i="<<i<<" status="<<system_soft[i].status<<std::endl;
//            }
//        }
//#endif
        // update n_loc_iso, n_glb_iso for isolated clusters
        //PS::S64 n_loc_iso = system_soft.getNumberOfParticleLocal();
        //PS::S64 n_glb_iso = system_soft.getNumberOfParticleGlobal();

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        // For connected clusters
        system_hard_connected.setPtclForConnectedCluster(system_soft, search_cluster.mediator_sorted_id_cluster_, search_cluster.ptcl_recv_);
        system_hard_connected.findGroupsAndCreateArtificialParticlesOMP<SystemSoft, FPSoft>(system_soft, _dt_tree);
        // send updated particle back to original (set zero mass particle to origin)
        search_cluster.writeAndSendBackPtcl(system_soft, system_hard_connected.getPtcl(), remove_list);
        system_hard_connected.updateTimeWriteBack();
#endif

        // update total particle number including artificial particles
        stat.n_all_loc = system_soft.getNumberOfParticleLocal();
        stat.n_all_glb = system_soft.getNumberOfParticleGlobal();
            
        // >2.4 set adr/rank for artificial particles in GPS
#pragma omp parallel for
        for(PS::S32 i=stat.n_real_loc; i<stat.n_all_loc; i++){
            system_soft[i].rank_org = my_rank;
            system_soft[i].adr = i;
#ifdef HARD_DEBUG
            assert(system_soft[i].group_data.artificial.isArtificial());
#endif
        }
        // >3 Tree for force ----------------------------------------
#ifdef PROFILE
        profile.create_group.barrier();
        PS::Comm::barrier();
        profile.create_group.end();
#endif
    }

    //! calculate tree solf force
    void treeSoftForce() {
#ifdef PROFILE
        profile.tree_soft.start();

        tree_soft.clearNumberOfInteraction();
        tree_soft.clearTimeProfile();
#endif

#ifndef USE_SIMD
        tree_soft.calcForceAllAndWriteBack(CalcForceEpEpWithLinearCutoffNoSimd(),
#ifdef USE_QUAD
                                           CalcForceEpSpQuadNoSimd(),
#else
                                           CalcForceEpSpMonoNoSimd(),
#endif
                                           system_soft,
                                           dinfo);
#elif USE_GPU
        const PS::S32 n_walk_limit = 200;
        const PS::S32 tag_max = 1;
#ifdef PARTICLE_SIMULATOR_GPU_MULIT_WALK_INDEX
        tree_soft.calcForceAllAndWriteBackMultiWalkIndex(DispatchKernelWithSPIndex,
                                                         RetrieveKernel,
                                                         tag_max,
                                                         system_soft,
                                                         dinfo,
                                                         n_walk_limit);
#else // no multi-walk index
        tree_soft.calcForceAllAndWriteBackMultiWalk(DispatchKernelWithSP,
                                                    RetrieveKernel,
                                                    tag_max,
                                                    system_soft,
                                                    dinfo,
                                                    n_walk_limit);
#endif // multi-walk index

#else // end gpu
        tree_soft.calcForceAllAndWriteBack(CalcForceEpEpWithLinearCutoffSimd(),
#ifdef USE_QUAD
                                           CalcForceEpSpQuadSimd(),
#else // no quad
                                           CalcForceEpSpMonoSimd(),
#endif // end quad
                                           system_soft,
                                           dinfo);
#endif // end gpu

#ifdef PROFILE
        n_count.ep_ep_interact     += tree_soft.getNumberOfInteractionEPEPLocal();
        n_count_sum.ep_ep_interact += tree_soft.getNumberOfInteractionEPEPGlobal();
        n_count.ep_sp_interact     += tree_soft.getNumberOfInteractionEPSPLocal();
        n_count_sum.ep_sp_interact += tree_soft.getNumberOfInteractionEPSPGlobal(); 

        tree_soft_profile += tree_soft.getTimeProfile();
        domain_decompose_weight = tree_soft_profile.calc_force;

        profile.tree_soft.barrier();
        PS::Comm::barrier();
        profile.tree_soft.end();
#endif
    }

    // correct force due to change over function
    void treeForceCorrectChangeover() {
#ifdef PROFILE
        profile.force_correct.start();
#endif

#ifdef CORRECT_FORCE_DEBUG
        // backup particle data
        PS::S64 n_loc_all = system_soft.getNumberOfParticleLocal();

        FPSoft psys_bk[n_loc_all];
#pragma omp parallel for
        for (int i=0; i<n_loc_all; i++) 
            psys_bk[i] = system_soft[i];
#endif

        // single 
        system_hard_one_cluster.correctPotWithCutoffOMP(system_soft, search_cluster.getAdrSysOneCluster());

        // Isolated clusters
        system_hard_isolated.correctForceWithCutoffClusterOMP(system_soft);

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
        // Connected clusters
        system_hard_connected.correctForceWithCutoffTreeNeighborAndClusterOMP<SystemSoft, FPSoft, TreeForce, EPJSoft>(system_soft, tree_soft, search_cluster.getAdrSysConnectClusterSend());
#endif

#ifdef CORRECT_FORCE_DEBUG
        // switch backup and sys particle data
#pragma omp parallel for
        for (int i=0; i<n_loc_all; i++) {
            FPSoft ptmp = psys_bk[i];
            psys_bk[i] = system_soft[i];
            system_soft[i] = ptmp;
        }

        // all particles
        system_hard_isolated.correctForceWithCutoffTreeNeighborOMP<SystemSoft, FPSoft, TreeForce, EPJSoft>(system_soft, tree_soft, n_loc);

        // single 
        //system_hard_one_cluster.correctPotWithCutoffOMP(system_soft, search_cluster.getAdrSysOneCluster());
        // Isolated clusters
        //system_hard_isolated.correctForceWithCutoffTreeNeighborOMP<SystemSoft, FPSoft, TreeForce, EPJSoft>(system_soft, tree_soft, n_loc);
//#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
//        // Connected clusters
//        system_hard_connected.correctForceWithCutoffTreeNeighborOMP<SystemSoft, FPSoft, TreeForce, EPJSoft>(system_soft, tree_soft, n_loc_all);
//#endif

        for (int i=0; i<n_loc_all; i++) {
            PS::F64vec dacci=psys_bk[i].acc- system_soft[i].acc;
            PS::F64 dpoti = psys_bk[i].pot_tot- system_soft[i].pot_tot;
            if(dacci*dacci/(system_soft[i].acc*system_soft[i].acc)>1e-8) {
                std::cerr<<"Corrected Acc diff >1e-8: i "<<i<<" acc(tree): "<<system_soft[i].acc<<" acc(cluster): "<<psys_bk[i].acc<<std::endl;
                abort();
            }
            // Notice, the cluster members can include particles are not in the tree neighbor searching. If the extra neighbors in clusters are group members. The correct from clusters will be more accurate than tree neighbor correction. Since the potential is calculated from the real members instead of artificial particles. This can result in different potential
            if(abs(dpoti/system_soft[i].pot_tot)>1e-6) {
                std::cerr<<"Corrected pot diff >1e-6: i "<<i<<" pot(tree): "<<system_soft[i].pot_tot<<" pot(cluster): "<<psys_bk[i].pot_tot<<std::endl;
                abort();
            }
        }
        
#endif
#ifdef PROFILE
        profile.force_correct.barrier();
        PS::Comm::barrier();
        profile.force_correct.end();
#endif
    }

#ifdef KDKDK_4TH
    //! gradient kick for KDKDK_4TH method
    void GradientKick() {
#ifdef PROFILE
        profile.tree_soft.start();
        tree_soft.clearNumberOfInteraction();
        tree_soft.clearTimeProfile();
#endif
        // correction calculation
        //tree_soft.setParticaleLocalTree(system_soft, false);
        
        tree_soft.calcForceAllAndWriteBack(CalcCorrectEpEpWithLinearCutoffNoSimd(),
#ifdef USE_QUAD
                                           CalcForceEpSpQuadNoSimd(),
#else
                                           CalcForceEpSpMonoNoSimd(),
#endif
                                           system_soft,
                                           dinfo);

#ifdef PROFILE
        n_count.ep_ep_interact     += tree_soft.getNumberOfInteractionEPEPLocal();
        n_count_sum.ep_ep_interact += tree_soft.getNumberOfInteractionEPEPGlobal();
        n_count.ep_sp_interact     += tree_soft.getNumberOfInteractionEPSPLocal();
        n_count_sum.ep_sp_interact += tree_soft.getNumberOfInteractionEPSPGlobal(); 

        tree_soft_profile += tree_soft.getTimeProfile();
        domain_decompose_weight += tree_soft_profile.calc_force;

        profile.tree_soft.barrier();
        PS::Comm::barrier();
        profile.tree_soft.end();
        profile.force_correct.start();
#endif 

        // Isolated clusters
        system_hard_isolated.correctForceWithCutoffClusterOMP(system_soft, true);

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
        // Connected clusters
        system_hard_connected.correctForceWithCutoffTreeNeighborAndClusterOMP<SystemSoft, FPSoft, TreeForce, EPJSoft>(system_soft, tree_soft, search_cluster.getAdrSysConnectClusterSend(), true);
#endif

#ifdef PROFILE
        profile.force_correct.barrier();
        PS::Comm::barrier();
        profile.force_correct.end();
#endif
//        if (true) {


// for debug
//            FILE* fdump;
//            if ( (fdump = fopen("acorr.dump","w")) == NULL) {
//                fprintf(stderr,"Error: Cannot open file acorr.dump\n");
//                abort();
//            }
//            for (int i=0; i<n_loc; i++) {
//                fprintf(fdump, "%d %.12g %.12g %.12g %.12g %.12g %.12g\n", i, system_soft[i].acc[0], system_soft[i].acc[1], system_soft[i].acc[2], system_soft[i].acorr[0], system_soft[i].acorr[1], system_soft[i].acorr[2]);
//            }
//            fclose(fdump);
//            abort();
// debug
//            }

    }
#endif

    //!leap frog kick for single----------------------------------------------
    /* modify the velocity of particle in global system
       reset particle type to single
       @param[in,out] _sys: particle system
       @param[in]: _dt: tree step
       @param[in]; _adr: address for single particles
    */
    void kickOne(SystemSoft & _sys, 
                 const PS::F64 _dt, 
                 const PS::ReallocatableArray<PS::S32>& _adr) {
        const PS::S64 n= _adr.size();
#pragma omp parallel for
        for(PS::S32 i=0; i<n; i++){
            const PS::S32 k=_adr[i];
#ifdef KDKDK_4TH
            _sys[k].vel  += _dt*(_sys[k].acc + 9.0/192.0*_dt*_dt*_sys[k].acorr); 
#else
            _sys[k].vel  += _sys[k].acc * _dt;
#endif
            _sys[k].group_data.artificial.setParticleTypeToSingle();
        }
    }

    //!leap frog kick for clusters
    /*! modify the velocity of particle in local, if particle is from remote note and is not group member, do nothing, need MPI receive to update data
       Recover the mass of members for energy calculation
       @param[in,out] _sys: particle system
       @param[in,out] _ptcl: local particle array in system hard
       @param[in]: _dt: tree step
    */
    template<class Tptcl>
    void kickClusterAndRecoverGroupMemberMass(SystemSoft& _sys,
                                              PS::ReallocatableArray<Tptcl>& _ptcl,
                                              const PS::F64 _dt) {
        assert(Ptcl::group_data_mode == GroupDataMode::artificial);
        const PS::S64 n= _ptcl.size();
#pragma omp parallel for
        for(PS::S32 i=0; i<n; i++) {
            auto& pi_artificial = _ptcl[i].group_data.artificial;
            // if is group member, recover mass and kick due to c.m. force
            if (pi_artificial.isMember()) {
                const PS::S64 cm_adr = _ptcl[i].getParticleCMAddress(); 
                if (cm_adr>0) {
#ifdef HARD_DEBUG
                    assert(pi_artificial.getMassBackup()>0); 
                    assert(cm_adr<stat.n_all_loc);
#endif
                    _ptcl[i].mass = pi_artificial.getMassBackup();
#ifdef KDKDK_4TH
                    _ptcl[i].vel  += _dt*(_sys[cm_adr].acc + 9.0/192.0*_dt*_dt*_sys[cm_adr].acorr); 
#else
#ifdef NAN_CHECK_DEBUG
                    assert(!std::isnan(_sys[cm_adr].acc[0]));
                    assert(!std::isnan(_sys[cm_adr].acc[1]));
                    assert(!std::isnan(_sys[cm_adr].acc[2]));
#endif
                    _ptcl[i].vel += _sys[cm_adr].acc * _dt;
#endif
                    continue;
                }
            }
            // non-member particle
            const PS::S64 i_adr =_ptcl[i].adr_org;
            if(i_adr>=0) {
                // not remote particles
#ifdef HARD_DEBUG
                assert(pi_artificial.getStatus()==_sys[i_adr].group_data.artificial.getStatus());
#endif

#ifdef KDKDK_4TH
                _ptcl[i].vel  += _dt*(_sys[i_adr].acc + 9.0/192.0*_dt*_dt*_sys[i_adr].acorr); 
#else
                _ptcl[i].vel += _sys[i_adr].acc * _dt;
#endif
            }
        }
    }

    //!leap frog kick for sending list
    /*! Kick single particles in sending list 
       @param[in,out] _sys: particle system
       @param[in,out] _ptcl: local particle array in system hard
       @param[in]: _dt: tree step
    */
    void kickSend(SystemSoft& _sys,
                  const PS::ReallocatableArray<PS::S32>& _adr_ptcl_send,
                  const PS::F64 _dt) {
        assert(Ptcl::group_data_mode == GroupDataMode::artificial);
        const PS::S64 n= _adr_ptcl_send.size();
#pragma omp parallel for
        for(PS::S32 i=0; i<n; i++) {
            const PS::S64 adr = _adr_ptcl_send[i];
            // if it is group member with artificial particles, should not do kick since the required c.m. forces is on remote nodes;
            // if it is group member without artificial particles, also kick
            if(_sys[adr].group_data.artificial.isSingle() || (_sys[adr].group_data.artificial.isMember() && _sys[adr].getParticleCMAddress()<0)) {
                _sys[adr].vel += _sys[adr].acc * _dt;
#ifdef KDKDK_4TH
                _sys[adr].vel += _dt*_dt* _sys[adr].acorr /48; 
#endif
            }

#ifdef HARD_DEBUG
            if(_sys[adr].group_data.artificial.isSingle()) assert(_sys[adr].mass>0);
            else assert(_sys[adr].group_data.artificial.getMassBackup()>0);
#endif
        }
    }

    //! kick for artifical c.m. particles
    /*! Kick c.m. velocity (both isolated and connected clusters)
      @param[in,out] _sys: particle system
      @param[in] _adr_artificial_start: c.m. particle starting address
      @param[in] _ap_manager: artificial particle manager
      @param[in] _dt: tree step
    */
    void kickCM(SystemSoft& _sys,
                const PS::S32 _adr_artificial_start,
                ArtificialParticleManager& _ap_manager,
                const PS::F64 _dt) {
        const PS::S64 n_tot= _sys.getNumberOfParticleLocal();
        const PS::S32 n_artifical_per_group = _ap_manager.getArtificialParticleN();
#pragma omp parallel for
        for(PS::S32 i=_adr_artificial_start; i<n_tot; i+= n_artifical_per_group) {
            auto* pcm = _ap_manager.getCMParticles(&(_sys[i]));
            pcm->vel += pcm->acc * _dt;
#ifdef KDKDK_4TH
            pcm->vel += _dt*_dt* pcm->acorr /48; 
#endif
#ifdef HARD_DEBUG
            assert(pcm->group_data.artificial.isCM());
#endif
        }
    }

    //! drift for single particles
    void driftSingle(SystemSoft & system,
                     const PS::F64 dt){
        const PS::S32 n = system.getNumberOfParticleLocal();
#pragma omp parallel for
        for(PS::S32 i=0; i<n; i++){
            if(system[i].n_ngb <= 0){
                system[i].pos  += system[i].vel * dt;
            }
        }
    }

    //! kick
    void kick(const PS::F64 _dt_kick) {
#ifdef PROFILE
        profile.kick.start();
#endif
        /// Member mass are recovered
        // single and reset particle type to single (due to binary disruption)
        kickOne(system_soft, _dt_kick, search_cluster.getAdrSysOneCluster());
        // isolated
        kickClusterAndRecoverGroupMemberMass(system_soft, system_hard_isolated.getPtcl(), _dt_kick);
        // c.m. artificial
        kickCM(system_soft, stat.n_real_loc, hard_manager.ap_manager, _dt_kick);
        
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        // connected
        kickClusterAndRecoverGroupMemberMass(system_soft, system_hard_connected.getPtcl(), _dt_kick);
        // sending list for connected clusters
        kickSend(system_soft, search_cluster.getAdrSysConnectClusterSend(), _dt_kick);
        // send kicked particle from sending list, and receive remote single particle
        search_cluster.SendSinglePtcl(system_soft, system_hard_connected.getPtcl());
#endif

        time_kick += _dt_kick;

#ifdef HARD_DEBUG
        PS::S32 kick_regist[stat.n_real_loc];
        for(int i=0; i<stat.n_real_loc; i++) kick_regist[i] = 0;
        for(int i=0; i<search_cluster.getAdrSysOneCluster().size(); i++) {
            kick_regist[search_cluster.getAdrSysOneCluster()[i]]++;
        }
        for(int i=0; i<system_hard_isolated.getPtcl().size(); i++) {
            PS::S64 adr= system_hard_isolated.getPtcl()[i].adr_org;
            assert(adr>=0);
            kick_regist[adr]++;
        }
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        for(int i=0; i<system_hard_connected.getPtcl().size(); i++) {
            PS::S64 adr= system_hard_connected.getPtcl()[i].adr_org;
            if(adr>=0) kick_regist[adr]++;
        }
        for(int i=0; i<search_cluster.getAdrSysConnectClusterSend().size(); i++) {
            PS::S64 adr= search_cluster.getAdrSysConnectClusterSend()[i];
            kick_regist[adr]++;
        }
#endif
        for(int i=0; i<stat.n_real_loc; i++) {
            assert(kick_regist[i]==1);
        }
#endif

#ifdef PROFILE
        profile.kick.barrier();
        PS::Comm::barrier();
        profile.kick.end();
#endif
    }

    //! hard drift
    /*! \return interrupted cluster total number in all MPI processors
     */
    PS::S32 drift(const PS::F64 _dt_drift) {
        ////// set time
        //system_hard_one_cluster.setTimeOrigin(stat.time);
        //system_hard_isolated.setTimeOrigin(stat.time);
        //system_hard_connected.setTimeOrigin(stat.time);
        ////// set time

#ifdef PROFILE
        profile.hard_single.start();
#endif
        ////// integrater one cluster
        system_hard_one_cluster.initializeForOneCluster(search_cluster.getAdrSysOneCluster().size());
        system_hard_one_cluster.setPtclForOneClusterOMP(system_soft, search_cluster.getAdrSysOneCluster());
        system_hard_one_cluster.driveForOneClusterOMP(_dt_drift);
        //system_hard_one_cluster.writeBackPtclForOneClusterOMP(system_soft, search_cluster.getAdrSysOneCluster());
        system_hard_one_cluster.writeBackPtclForOneClusterOMP(system_soft);
        ////// integrater one cluster
#ifdef PROFILE
        profile.hard_single.barrier();
        PS::Comm::barrier();
        profile.hard_single.end();
#endif
        ////////////////

        /////////////
#ifdef PROFILE
        profile.hard_isolated.start();
#endif
        // reset slowdown energy correction
        system_hard_isolated.energy.resetEnergyCorrection();
        // integrate multi cluster A
        system_hard_isolated.driveForMultiClusterOMP(_dt_drift, &(system_soft[0]));
        //system_hard_isolated.writeBackPtclForMultiCluster(system_soft, search_cluster.adr_sys_multi_cluster_isolated_,remove_list);
        PS::S32 n_interrupt_isolated = system_hard_isolated.getNumberOfInterruptClusters();
        if(n_interrupt_isolated==0) system_hard_isolated.writeBackPtclForMultiCluster(system_soft, remove_list);
        // integrate multi cluster A

#ifdef PROFILE
        n_count.hard_interrupt += n_interrupt_isolated;
        profile.hard_isolated.barrier();
#endif

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
        n_interrupt_glb = PS::Comm::getSum(n_interrupt_isolated);
#else 
        n_interrupt_glb = n_interrupt_isolated;
#endif

#ifdef PROFILE
        n_count_sum.hard_interrupt += n_interrupt_glb;
        profile.hard_isolated.end();
#endif
        /////////////

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
        /////////////
#ifdef PROFILE
        profile.hard_connected.start();
#endif
        // reset slowdown energy correction
        system_hard_connected.energy.resetEnergyCorrection();
        // integrate multi cluster B
        system_hard_connected.driveForMultiClusterOMP(_dt_drift, &(system_soft[0]));
        PS::S32 n_interrupt_connected = system_hard_connected.getNumberOfInterruptClusters();

#ifdef PROFILE
        n_count.hard_interrupt += n_interrupt_connected;
        profile.hard_connected.barrier();
#endif

        PS::S32 n_interrupt_connected_glb = PS::Comm::getSum(n_interrupt_connected);

#ifdef PROFILE
        n_count_sum.hard_interrupt += n_interrupt_connected_glb;
        profile.hard_connected.end();
        profile.hard_connected.start();
#endif


        if (n_interrupt_connected_glb==0) {
            search_cluster.writeAndSendBackPtcl(system_soft, system_hard_connected.getPtcl(), remove_list);
            system_hard_connected.updateTimeWriteBack();
        }
        // integrate multi cluster B

        n_interrupt_glb += n_interrupt_connected_glb;
#ifdef PROFILE
        profile.hard_connected.barrier();
        PS::Comm::barrier();
        profile.hard_connected.end();
#endif
#endif
        
        if (n_interrupt_glb==0) Ptcl::group_data_mode = GroupDataMode::cm;
        
        return n_interrupt_glb;
    }

    //! finish interrupted drift
    /*! \return new interrupt number
     */
    PS::S32 finishInterruptDrift() {
#ifdef PROFILE
        profile.hard_interrupt.start();
        profile.hard_isolated.start();
#endif

        // finish interrupt clusters first
        // isolated clusters
        PS::S32 n_interrupt_isolated = 0;
        if (system_hard_isolated.getNumberOfInterruptClusters()>0) {
            system_hard_isolated.finishIntegrateInterruptClustersOMP();
            n_interrupt_isolated = system_hard_isolated.getNumberOfInterruptClusters();
            if(n_interrupt_isolated==0) system_hard_isolated.writeBackPtclForMultiCluster(system_soft, remove_list);
        }

#ifdef PROFILE
        n_count.hard_interrupt += n_interrupt_isolated;
        profile.hard_isolated.barrier();
#endif

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        n_interrupt_glb = PS::Comm::getSum(n_interrupt_isolated);
#else
        n_interrupt_glb = n_interrupt_isolated;
#endif

#ifdef PROFILE
        n_count_sum.hard_interrupt += n_interrupt_glb;
        profile.hard_isolated.end();
#endif

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
#ifdef PROFILE
        profile.hard_connected.start();
#endif

        // connected clusters
        PS::S32 n_interrupt_connected = 0;
        if (system_hard_connected.getNumberOfInterruptClusters()>0) {
            system_hard_connected.finishIntegrateInterruptClustersOMP();
            n_interrupt_connected = system_hard_connected.getNumberOfInterruptClusters();
        }

#ifdef PROFILE
        n_count.hard_interrupt += n_interrupt_connected;
        profile.hard_connected.barrier();
#endif

        PS::S32 n_interrupt_connected_glb = PS::Comm::getSum(n_interrupt_connected);

#ifdef PROFILE
        n_count_sum.hard_interrupt += n_interrupt_connected_glb;
        profile.hard_connected.end();
        profile.hard_connected.start();
#endif

        if (n_interrupt_connected_glb==0) {
            search_cluster.writeAndSendBackPtcl(system_soft, system_hard_connected.getPtcl(), remove_list);
            system_hard_connected.updateTimeWriteBack();
        }

        n_interrupt_glb += n_interrupt_connected_glb;

#ifdef PROFILE
        profile.hard_connected.barrier();
        PS::Comm::barrier();
        profile.hard_connected.end();
#endif
#endif

#ifdef HARD_INTERRUPT_PRINT
        if (n_interrupt_glb>0)  {
            std::cerr<<"Interrupt detected, number: "<<n_interrupt_glb<<std::endl;
        }
#endif

#ifdef PROFILE
        profile.hard_interrupt.barrier();
        PS::Comm::barrier();
        profile.hard_interrupt.end();
#endif

        if (n_interrupt_glb==0) Ptcl::group_data_mode = GroupDataMode::cm;

        // if interrupt cluster still exist, return the number immediately without new integration.
        return n_interrupt_glb;
    }


    //! Calculate the maximum time step limit for next block step
    /*! Get the maximum time step allown for next block step
      Basic algorithm: the integer of time/dt_min is the binary tree for block step, counting from the minimum digital, the last zero indicate the maximum block step level allown for next step
      @param[in] _time: current time
      @param[in] _dt_max: maximum time step allown
      @param[in] _dt_min: minimum time step allown
    */
    inline PS::F64 calcDtLimit(const PS::F64 _time,
                               const PS::F64 _dt_max,
                               const PS::F64 _dt_min){
        // for first step, the maximum time step is OK
        if(_time==0.0) return _dt_max;
        else {
            // the binary tree for current time position in block step 
            PS::U64 bitmap = _time/_dt_min;
//#ifdef __GNUC__ 
//        PS::S64 dts = __builtin_ctz(bitmap) ;
//        PS::U64 c = (1<<dts);
////        std::cerr<<"time = "<<_time<<"  dt_min = "<<_dt_min<<"  bitmap = "<<bitmap<<"  dts = "<<dts<<std::endl;
//#else

            // block step multiply factor 
            PS::U64 c=1;
            // find the last zero in the binary tree to obtain the current block step level
            while((bitmap&1)==0) {
                bitmap = (bitmap>>1);
                c = (c<<1);
            }
//#endif
            // return the maximum step allown
            return std::min(c*_dt_min,_dt_max);
        }
    }
    
    //! check time consistence
    bool checkTimeConsistence() {
        assert(abs(time_kick-stat.time)<1e-13);
        time_kick = stat.time; // escape the problem of round-off error
        assert(stat.time == system_hard_one_cluster.getTimeOrigin());
        assert(stat.time == system_hard_isolated.getTimeOrigin());
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        assert(stat.time == system_hard_connected.getTimeOrigin());
#endif
        return true;
    }
    
    //! write back hard particles to global system
    void writeBackHardParticles() {
#ifdef PROFILE
        profile.output.start();
#endif
        system_hard_isolated.writeBackPtclForMultiCluster(system_soft, remove_list);
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        // update gloabl particle system and send receive remote particles
        search_cluster.writeAndSendBackPtcl(system_soft, system_hard_connected.getPtcl(), remove_list);
        system_hard_connected.updateTimeWriteBack();
#endif        
#ifdef PROFILE
        profile.output.barrier();
        PS::Comm::barrier();
        profile.output.end();
#endif
    }

    // update group_data.cm to pcm data for search cluster after restart
    void setParticleGroupDataToCMData() {
#ifdef PROFILE
        profile.search_cluster.start();
#endif
#ifdef CLUSTER_VELOCITY
        // update status and mass_bk to pcm data for search cluster after restart
        system_hard_one_cluster.resetParticleGroupData(system_soft);
        system_hard_isolated.setParticleGroupDataToCMData(system_soft);
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        system_hard_connected.setParticleGroupDataToCMData(system_soft);
        search_cluster.writeAndSendBackPtcl(system_soft, system_hard_connected.getPtcl(), remove_list);
        system_hard_connected.updateTimeWriteBack();
#endif
#ifdef PROFILE
        profile.search_cluster.barrier();
        PS::Comm::barrier();
        profile.search_cluster.end();
#endif
    }
    
#endif    

    //! correct force due to the change over update
    void correctForceChangeOverUpdate() {
#ifdef PROFILE
        profile.force_correct.start();
#endif
        // correct changeover for first step
        // Isolated clusters
        system_hard_isolated.correctForceForChangeOverUpdateOMP<SystemSoft, TreeForce, EPJSoft>(system_soft, tree_soft);

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
        // Connected clusters
        auto& adr_send = search_cluster.getAdrSysConnectClusterSend();
        system_hard_connected.correctForceForChangeOverUpdateOMP<SystemSoft, TreeForce, EPJSoft>(system_soft, tree_soft, adr_send.getPointer(), adr_send.size());
#endif

#ifdef PROFILE
        profile.force_correct.barrier();
        PS::Comm::barrier();
        profile.force_correct.end();
#endif

    }

    //! domain decomposition
    void domainDecompose() {
#ifdef PROFILE
        // > 6. Domain decomposition
        profile.domain.start();
#endif
        // Domain decomposition, parrticle exchange and force calculation
        if(n_loop % 16 == 0) {
            dinfo.decomposeDomainAll(system_soft,domain_decompose_weight);
            //std::cout<<"rank: "<<my_rank<<" weight: "<<domain_decompose_weight<<std::endl;
        }
#ifdef PROFILE
        profile.domain.barrier();
        PS::Comm::barrier();
        profile.domain.end();
#endif
    }

    //! adjust dt_reduce_factor for tree time step
    /*!
      \return true: if tree time need update
     */
    bool adjustDtTreeReduce(PS::F64& _dt_reduce_factor, const PS::F64 _dt_tree_base, const PS::F64 _dt_tree_now) {
        bool dt_mod_flag = false;
        PS::F64 dt_tree_new = _dt_tree_base / _dt_reduce_factor;
        if(_dt_tree_now!= dt_tree_new) {
            // in increasing case, need to make sure the time is consistent with block step time/dt_tree
            if(_dt_tree_now<dt_tree_new) {
                PS::F64 dt_tree_max = calcDtLimit(stat.time, _dt_tree_base, _dt_tree_now);
                if (dt_tree_max == _dt_tree_now) {
                    // no modificatoin, recover original reduce factor
                    _dt_reduce_factor= _dt_tree_base/_dt_tree_now;
                }
                else {
                    dt_mod_flag = true;
                    // limit dt_reduce_factor to the maximum block step allown
                    _dt_reduce_factor = std::min(_dt_tree_base/dt_tree_max, _dt_reduce_factor);

                }
            }
            else dt_mod_flag = true;
        }

        return dt_mod_flag;
    }

    //! update system status
    void updateStatus(const bool _initial_flag) {
#ifdef PROFILE
        profile.status.start();
#endif
        if (_initial_flag) {
            // calculate initial energy
            stat.energy.clear();
            stat.energy.calc(&system_soft[0], stat.n_real_loc, true);
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            stat.energy.getSumMultiNodes(true);
#endif
#ifdef HARD_CHECK_ENERGY
            stat.energy.ekin_sd = stat.energy.ekin;
            stat.energy.epot_sd = stat.energy.epot;
            stat.energy_hard_diff = 0;
            stat.energy_hard_sd_diff = 0;
#endif
            stat.calcCenterOfMass(&system_soft[0], stat.n_real_loc);

        }
        else {
#ifdef HARD_CHECK_ENERGY
            // reset sd energy reference to no slowdown case
            PS::F64 energy_old = stat.energy.ekin + stat.energy.epot;
            stat.energy.etot_sd_ref -= stat.energy.ekin_sd + stat.energy.epot_sd - energy_old;
#endif

            stat.energy.calc(&system_soft[0], stat.n_real_loc);
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            stat.energy.getSumMultiNodes();
#endif

#ifdef HARD_CHECK_ENERGY
            HardEnergy energy_local = system_hard_isolated.energy;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            energy_local += system_hard_connected.energy;
#endif
            // hard energy error
            stat.energy_hard_diff += PS::Comm::getSum(energy_local.de);
            stat.energy_hard_sd_diff += PS::Comm::getSum(energy_local.de_sd);
            // energy correction due to slowdown 
            PS::F64 ekin_sd_correction = PS::Comm::getSum(energy_local.ekin_sd_correction);
            stat.energy.ekin_sd = stat.energy.ekin + ekin_sd_correction;
            PS::F64 epot_sd_correction = PS::Comm::getSum(energy_local.epot_sd_correction);
            stat.energy.epot_sd = stat.energy.epot + epot_sd_correction;
            PS::F64 etot_sd_correction = ekin_sd_correction + epot_sd_correction;
            // for total energy reference, first add the cumulative change due to slowdown change in the integration (referring to no slowdown case), then add slowdown energy correction from current time
            stat.energy.etot_sd_ref += PS::Comm::getSum(energy_local.de_sd_change_cum) + etot_sd_correction;

            system_hard_isolated.energy.clear();
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            system_hard_connected.energy.clear();
#endif
#endif
            stat.calcCenterOfMass(&system_soft[0], stat.n_real_loc);
        }
#ifdef PROFILE
        profile.status.barrier();
        PS::Comm::barrier();
        profile.status.end();
#endif
    }

    void output() {
#ifdef PROFILE
        profile.output.start();
#endif
        bool print_flag = input_parameters.print_flag;
        int write_style = input_parameters.write_style.value;

        // print status
        if(print_flag) {
            std::cout<<std::endl;
            stat.print(std::cout);
        }
        // write status, output to separate snapshots
        if(write_style==1) {
            // status output
            if(my_rank==0) {
                stat.printColumn(fstatus, WRITE_WIDTH);
                fstatus<<std::endl;
            }

            // data output
            file_header.n_body = stat.n_real_glb;
            file_header.time = stat.time;
            file_header.nfile++;
            std::string fname = input_parameters.fname_snp.value+"."+std::to_string(file_header.nfile);
#ifdef PETAR_DEBUG
            assert(system_soft.getNumberOfParticleLocal()== stat.n_all_loc);
#endif
            system_soft.setNumberOfParticleLocal(stat.n_real_loc);
            if (input_parameters.data_format.value==1||input_parameters.data_format.value==3)
                system_soft.writeParticleAscii(fname.c_str(), file_header);
            else if(input_parameters.data_format.value==0||input_parameters.data_format.value==2)
                system_soft.writeParticleBinary(fname.c_str(), file_header);
            system_soft.setNumberOfParticleLocal(stat.n_all_loc);
        }
        // write all information in to fstatus
        else if(write_style==2&&my_rank==0) {
            // write snapshot with one line
            stat.printColumn(fstatus, WRITE_WIDTH);
            for (int i=0; i<stat.n_real_loc; i++) system_soft[i].printColumn(fstatus, WRITE_WIDTH);
            fstatus<<std::endl;
        }
        // write status only
        else if(write_style==3&&my_rank==0) {
            stat.printColumn(fstatus, WRITE_WIDTH);
            fstatus<<std::endl;
        }

#ifdef PROFILE
        profile.output.barrier();
        PS::Comm::barrier();
        profile.output.end();
#endif
    }

#ifdef PROFILE
    //! measure profiles
    void calcProfile() {
        
        // profile analysis

        PS::S32 n_hard_single     = system_hard_one_cluster.getPtcl().size();
        PS::S32 n_hard_isolated   = system_hard_isolated.getPtcl().size();

        n_count.hard_single      += n_hard_single;
        n_count.hard_isolated    += n_hard_isolated;

        n_count_sum.hard_single      += PS::Comm::getSum(n_hard_single);
        n_count_sum.hard_isolated    += PS::Comm::getSum(n_hard_isolated);

        PS::S64 ARC_substep_sum   = system_hard_isolated.ARC_substep_sum;
        PS::S64 ARC_tsyn_step_sum   = system_hard_isolated.ARC_tsyn_step_sum;
        PS::S64 ARC_n_groups      = system_hard_isolated.ARC_n_groups;
        PS::S64 ARC_n_groups_iso  = system_hard_isolated.ARC_n_groups_iso;
        PS::S64 H4_step_sum       = system_hard_isolated.H4_step_sum;
#ifdef HARD_COUNT_NO_NEIGHBOR
        PS::S64 n_neighbor_zero   = system_hard_isolated.n_neighbor_zero;
#endif

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        PS::S32 n_hard_connected  = system_hard_connected.getPtcl().size();
        n_count.hard_connected   += n_hard_connected;
        n_count_sum.hard_connected += PS::Comm::getSum(n_hard_connected);

        ARC_substep_sum += system_hard_connected.ARC_substep_sum;
        ARC_tsyn_step_sum += system_hard_connected.ARC_tsyn_step_sum;
        ARC_n_groups += system_hard_connected.ARC_n_groups;
        ARC_n_groups_iso += system_hard_connected.ARC_n_groups_iso;
        H4_step_sum +=  system_hard_connected.H4_step_sum;
#ifdef HARD_COUNT_NO_NEIGHBOR
        n_neighbor_zero+= system_hard_connected.n_neighbor_zero;
#endif
#endif
                                           
        n_count.ARC_substep_sum  += ARC_substep_sum;
        n_count.ARC_tsyn_step_sum+= ARC_tsyn_step_sum;
        n_count.ARC_n_groups     += ARC_n_groups;
        n_count.ARC_n_groups_iso += ARC_n_groups_iso;
        n_count.H4_step_sum      += H4_step_sum;
#ifdef HARD_COUNT_NO_NEIGHBOR
        n_count.n_neighbor_zero  += n_neighbor_zero;
#endif

        n_count_sum.ARC_substep_sum  += PS::Comm::getSum(ARC_substep_sum);
        n_count_sum.ARC_tsyn_step_sum+= PS::Comm::getSum(ARC_tsyn_step_sum);
        n_count_sum.ARC_n_groups     += PS::Comm::getSum(ARC_n_groups);
        n_count_sum.ARC_n_groups_iso     += PS::Comm::getSum(ARC_n_groups_iso);
        n_count_sum.H4_step_sum      += PS::Comm::getSum(H4_step_sum);
#ifdef HARD_COUNT_NO_NEIGHBOR
        n_count_sum.n_neighbor_zero  += PS::Comm::getSum(n_neighbor_zero);
#endif

        system_hard_isolated.ARC_substep_sum = 0;
        system_hard_isolated.ARC_tsyn_step_sum=0;
        system_hard_isolated.ARC_n_groups = 0;
        system_hard_isolated.ARC_n_groups_iso = 0;
        system_hard_isolated.H4_step_sum = 0;
#ifdef HARD_COUNT_NO_NEIGHBOR
        system_hard_isolated.n_neighbor_zero = 0;
#endif

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        system_hard_connected.ARC_substep_sum = 0;
        system_hard_connected.ARC_tsyn_step_sum=0;
        system_hard_connected.ARC_n_groups = 0;
        system_hard_connected.ARC_n_groups_iso = 0;
        system_hard_connected.H4_step_sum = 0;
#ifdef HARD_COUNT_NO_NEIGHBOR
        system_hard_connected.n_neighbor_zero = 0;
#endif
#endif
                                           
        n_count.cluster_count(1, n_hard_single);

        const PS::S32  n_isolated_cluster = system_hard_isolated.getNumberOfClusters();
        n_count.cluster_isolated += n_isolated_cluster;
        n_count_sum.cluster_isolated += PS::Comm::getSum(n_isolated_cluster);

        const PS::S32* isolated_cluster_n_list = system_hard_isolated.getClusterNumberOfMemberList();
        for (PS::S32 i=0; i<n_isolated_cluster; i++) n_count.cluster_count(isolated_cluster_n_list[i]);

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        const PS::S32  n_connected_cluster = system_hard_connected.getNumberOfClusters();
        n_count.cluster_connected += n_connected_cluster;
        n_count_sum.cluster_connected += PS::Comm::getSum(n_connected_cluster);
        const PS::S32* connected_cluster_n_list = system_hard_connected.getClusterNumberOfMemberList();
        for (PS::S32 i=0; i<n_connected_cluster; i++) n_count.cluster_count(connected_cluster_n_list[i]);
#endif
        dn_loop++;

    }

    //! clear profile
    void clearProfile() {
        profile.clear();
        tree_soft_profile.clear();
        tree_nb_profile.clear();
#if defined(USE_GPU) && defined(GPU_PROFILE)
        gpu_profile.clear();
        gpu_counter.clear();
#endif
        n_count.clear();
        n_count_sum.clear();
        dn_loop=0;
    }

    //! output profile data
    void printProfile() {

        const SysProfile& profile_min = profile.getMin();
        
        if(input_parameters.print_flag) {
            std::cout<<std::setprecision(5);
            std::cout<<"Tree step number: "<<dn_loop<<std::endl;
                

            std::cout<<"**** Wallclock time per step (local): [Min/Max]\n";
            //std::cout<<std::setw(PRINT_WIDTH)<<"Rank";
            profile.dumpName(std::cout,PRINT_WIDTH);
            std::cout<<std::endl;

            //std::cout<<std::setw(PRINT_WIDTH)<<my_rank;
            profile_min.dump(std::cout,PRINT_WIDTH,dn_loop);
            std::cout<<std::endl;
            profile.dump(std::cout,PRINT_WIDTH,dn_loop);
            std::cout<<std::endl;

            std::cout<<"**** FDPS tree soft force time profile (local):\n";
            tree_soft_profile.dumpName(std::cout,PRINT_WIDTH);
            std::cout<<std::endl;
            tree_soft_profile.dump(std::cout,PRINT_WIDTH,dn_loop);
            std::cout<<std::endl;

            std::cout<<"**** Tree neighbor time profile (local):\n";
            tree_nb_profile.dumpName(std::cout,PRINT_WIDTH);
            std::cout<<std::endl;
            tree_nb_profile.dump(std::cout,PRINT_WIDTH,dn_loop);
            std::cout<<std::endl;

#if defined(USE_GPU) && defined(GPU_PROFILE)
            std::cout<<"**** GPU time profile (local):\n";
            gpu_profile.dumpName(std::cout,PRINT_WIDTH);
            gpu_counter.dumpName(std::cout,PRINT_WIDTH);
            std::cout<<std::endl;
            gpu_profile.dump(std::cout,PRINT_WIDTH,dn_loop);
            gpu_counter.dump(std::cout,PRINT_WIDTH,dn_loop);
            std::cout<<std::endl;
#endif

            std::cout<<"**** Number per step (global):\n";
            n_count_sum.dumpName(std::cout,PRINT_WIDTH);
            std::cout<<std::endl;
            n_count_sum.dump(std::cout,PRINT_WIDTH,dn_loop);
            std::cout<<std::endl;
                
            std::cout<<"**** Number of members in clusters (local):\n";
            n_count.printHist(std::cout,PRINT_WIDTH,dn_loop);
        }

        if(input_parameters.write_style.value>0) {
            fprofile<<std::setprecision(WRITE_PRECISION);
            fprofile<<std::setw(WRITE_WIDTH)<<my_rank;
            fprofile<<std::setw(WRITE_WIDTH)<<stat.time
                    <<std::setw(WRITE_WIDTH)<<dn_loop
                    <<std::setw(WRITE_WIDTH)<<stat.n_real_loc;
            profile.dump(fprofile, WRITE_WIDTH, dn_loop);
            profile.dumpBarrier(fprofile, WRITE_WIDTH, dn_loop);
            tree_soft_profile.dump(fprofile, WRITE_WIDTH, dn_loop);
            tree_nb_profile.dump(fprofile, WRITE_WIDTH, dn_loop);
#if defined(USE_GPU) && defined(GPU_PROFILE)
            gpu_profile.dump(fprofile, WRITE_WIDTH, dn_loop);
            gpu_counter.dump(fprofile, WRITE_WIDTH, dn_loop);
#endif
            n_count.dump(fprofile, WRITE_WIDTH, dn_loop);
            fprofile<<std::endl;
        }
    }
#endif

public:

    //! get address of particle from an id, if not found, return -1
    PS::S32 getParticleAdrFromID(const PS::S64 _id) {
        auto item = id_adr_map.find(_id);
        if (item==id_adr_map.end()) return -1;
        else return item->second;
    }

    //! regist a particle 
    void addParticleInIdAdrMap(FPSoft& _ptcl) {
        id_adr_map[_ptcl.id] = _ptcl.adr;
    }

    //! reconstruct ID-address map
    void reconstructIdAdrMap() {
        id_adr_map.clear();
        for (PS::S32 i=0; i<stat.n_real_loc; i++) id_adr_map[system_soft[i].id] = i;
    }

    //! remove particle with id from map, return particle index, if not found return -1
    PS::S32 removeParticleFromIdAdrMap(const PS::S64 _id) {
        auto item = id_adr_map.find(_id);
        if (item==id_adr_map.end()) return -1;
        else {
            PS::S32 adr = item->second;
            id_adr_map.erase(item);
            return adr;
        }
    }

    //! remove artificial and unused particles
    void removeParticles() {
        /////////////
        assert(n_interrupt_glb==0);
#ifdef PROFILE
        profile.other.start();
#endif

        // Remove ghost particles
        system_soft.removeParticle(remove_list.getPointer(), remove_list.size());
        // reset particle number
        stat.n_real_loc = stat.n_real_loc-remove_list.size();
        system_soft.setNumberOfParticleLocal(stat.n_real_loc);
        stat.n_real_glb = system_soft.getNumberOfParticleGlobal();
        remove_list.resizeNoInitialize(0);
#ifdef PROFILE
        profile.other.barrier();
        PS::Comm::barrier();
        profile.other.end();
#endif
    }


    //! exchange particles
    void exchangeParticle() {
        assert(n_interrupt_glb==0);

#ifdef PROFILE
        profile.exchange.start();
#endif
        system_soft.exchangeParticle(dinfo);

        const PS::S32 n_loc = system_soft.getNumberOfParticleLocal();

#pragma omp parallel for
        for(PS::S32 i=0; i<n_loc; i++){
            system_soft[i].rank_org = my_rank;
            system_soft[i].adr = i;
        }

        // record real particle n_loc/glb
        stat.n_real_loc = n_loc;
#ifdef PETAR_DEBUG
        assert(stat.n_real_glb==system_soft.getNumberOfParticleGlobal());
#endif
        
#ifdef PROFILE
        profile.exchange.barrier();
        PS::Comm::barrier();
        profile.exchange.end();
#endif
    }

    //! initial FDPS (print LOGO) and MPI
    void initialFDPS(int argc, char *argv[]) {
        PS::Initialize(argc, argv);
        my_rank = PS::Comm::getRank();
        n_proc = PS::Comm::getNumberOfProc();
        initial_fdps_flag = true;
    }

    //! reading input parameters using getopt method
    /*! 
      @param[in] argc: number of options
      @param[in] argv: string of options
      \return -1 if help is used
     */
    int readParameters(int argc, char *argv[]) {
        //assert(initial_fdps_flag);
        // reading parameters
        read_parameters_flag = true;
        if (my_rank==0) input_parameters.print_flag=true;
        else input_parameters.print_flag=false;
        int read_flag = input_parameters.read(argc,argv);
        return read_flag;
    }
    
    //! reading data set from file, filename is given by readParameters
    void readDataFromFile() {
        assert(read_parameters_flag);

        // particle system
        if (!read_data_flag) system_soft.initialize();

        PS::S32 data_format = input_parameters.data_format.value;
        auto* data_filename = input_parameters.fname_inp.value.c_str();
        if(data_format==1||data_format==2||data_format==4)
            system_soft.readParticleAscii(data_filename, file_header);
        else
            system_soft.readParticleBinary(data_filename, file_header);
        PS::Comm::broadcast(&file_header, 1, 0);
        PS::S64 n_glb = system_soft.getNumberOfParticleGlobal();
        PS::S64 n_loc = system_soft.getNumberOfParticleLocal();

        if(input_parameters.print_flag)
            std::cout<<"Reading file "<<data_filename<<std::endl
                     <<"N_tot = "<<n_glb<<"\nN_loc = "<<n_loc<<std::endl;

        stat.time = file_header.time;
        stat.n_real_glb = stat.n_all_glb = n_glb;
        stat.n_real_loc = stat.n_all_loc = n_loc;

        input_parameters.n_glb.value = n_glb;

        read_data_flag = true;
    }

    //! reading data from particle array
    /*!
      @param[in] _n_partcle: number of particles
      @param[in] _mass: mass array
      @param[in] _x: position x array
      @param[in] _y: position y array
      @param[in] _z: position z array
      @param[in] _vx: position vx array
      @param[in] _vy: position vy array
      @param[in] _vz: position vz array
      @param[in] _id: id of particles (if provided, must >0)
     */
    void readDataFromArray(PS::S64 _n_partcle, PS::F64* _mass, PS::F64* _x, PS::F64* _y, PS::F64* _z, PS::F64* _vx, PS::F64* _vy, PS::F64* _vz, PS::S64* _id=NULL) {
        assert(read_parameters_flag);

        // particle system
        if (!read_data_flag) system_soft.initialize();

        PS::S64 n_glb = _n_partcle;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
        PS::S64 n_loc = n_glb / n_proc; 
        if( n_glb % n_proc > my_rank) n_loc++;
        PS::S64 i_h = n_glb/n_proc*my_rank;
        if( n_glb % n_proc  > my_rank) i_h += my_rank;
        else i_h += n_glb % n_proc;
#else
        PS::S64 n_loc = n_glb;
        PS::S64 i_h = 0;
#endif

        system_soft.setNumberOfParticleLocal(n_loc);
        if (_id!=NULL) {
            for(PS::S32 i=0; i<n_loc; i++){
                system_soft[i].mass = _mass[i_h+i];
                system_soft[i].pos.x = _x[i_h+i];
                system_soft[i].pos.y = _y[i_h+i];
                system_soft[i].pos.z = _z[i_h+i];
                system_soft[i].vel.x = _vx[i_h+i];
                system_soft[i].vel.y = _vy[i_h+i];
                system_soft[i].vel.z = _vz[i_h+i];
                system_soft[i].id = _id[i_h+i];
                system_soft[i].group_data.artificial.setParticleTypeToSingle();
                assert(_id[i_h+i]>0);
            }
        }
        else {
            for(PS::S32 i=0; i<n_loc; i++){
                system_soft[i].mass = _mass[i_h+i];
                system_soft[i].pos.x = _x[i_h+i];
                system_soft[i].pos.y = _y[i_h+i];
                system_soft[i].pos.z = _z[i_h+i];
                system_soft[i].vel.x = _vx[i_h+i];
                system_soft[i].vel.y = _vy[i_h+i];
                system_soft[i].vel.z = _vz[i_h+i];
                system_soft[i].id = i_h+i+1;
                system_soft[i].group_data.artificial.setParticleTypeToSingle();
            }
        }
        file_header.nfile = 0;
        file_header.n_body = n_glb;
        file_header.time = 0.0;

        stat.time = 0.0;
        stat.n_real_glb = stat.n_all_glb = n_glb;
        stat.n_real_loc = stat.n_all_loc = n_loc;

        input_parameters.n_glb.value = n_glb;

        read_data_flag = true;
    }

    //! generate data from plummer model
    void generatePlummer() {
        // ensure parameters are used
        assert(read_parameters_flag);

        // particle system
        if (!read_data_flag) system_soft.initialize();

        PS::S64 n_glb = input_parameters.n_glb.value;
        assert(n_glb>0);

        PS::S64 n_loc = n_glb / n_proc; 
        if( n_glb % n_proc > my_rank) n_loc++;
        system_soft.setNumberOfParticleLocal(n_loc);

        PS::F64 * mass;
        PS::F64vec * pos;
        PS::F64vec * vel;

        const PS::F64 m_tot = 1.0;
        const PS::F64 eng = -0.25;

        ParticleDistributionGenerator::makePlummerModel(m_tot, n_glb, n_loc, mass, pos, vel, eng);

        PS::S64 i_h = n_glb/n_proc*my_rank;
        if( n_glb % n_proc  > my_rank) i_h += my_rank;
        else i_h += n_glb % n_proc;
        for(PS::S32 i=0; i<n_loc; i++){
            system_soft[i].mass = mass[i];
            system_soft[i].pos = pos[i];
            system_soft[i].vel = vel[i];
            system_soft[i].id = i_h + i + 1;
            system_soft[i].group_data.artificial.setParticleTypeToSingle();
        }

        file_header.nfile = 0;
        file_header.n_body = n_glb;
        file_header.time = 0.0;

        stat.time = 0.0;
        stat.n_real_glb = stat.n_all_glb = n_glb;
        stat.n_real_loc = stat.n_all_loc = n_loc;

        read_data_flag = true;

        delete [] mass;
        delete [] pos;
        delete [] vel;
    }

    //! create kepler disk
    /*! No center sun is set, suppressed 
     */
    void generateKeplerDisk(const PS::F64 ax_in, // [AU]
                            const PS::F64 ax_out, // [AU]
                            const PS::F64 ecc_rms, // normalized
                            const PS::F64 inc_rms, // normalized
                            const PS::F64 dens = 10.0, // [g/cm^2]
                            const PS::F64 mass_sun = 1.0, //[m_sun]
                            const double a_ice = 0.0,
                            const double f_ice = 1.0,
                            const double power = -1.5,
                            const PS::S32 seed = 0) {
        
        // ensure parameters are used
        assert(read_parameters_flag);

        // particle system
        if (!read_data_flag) system_soft.initialize();

        PS::S64 n_glb = input_parameters.n_glb.value;
        assert(n_glb>0);

        PS::S64 n_loc = n_glb / n_proc; 
        if( n_glb % n_proc > my_rank) n_loc++;
        system_soft.setNumberOfParticleLocal(n_loc);

        PS::F64 * mass;
        PS::F64vec * pos;
        PS::F64vec * vel;

        PS::F64 mass_planet_glb;
        ParticleDistributionGenerator::makeKeplerDisk(mass_planet_glb, mass, pos, vel, n_glb, n_loc,
                                                      ax_in, ax_out, ecc_rms, inc_rms, dens, mass_sun, a_ice, f_ice, power, seed);

        PS::S64 i_h = n_glb/n_proc*my_rank;
        if( n_glb % n_proc  > my_rank) i_h += my_rank;
        else i_h += n_glb % n_proc;
        for(PS::S32 i=0; i<n_loc; i++){
            system_soft[i].mass = mass[i];
            system_soft[i].pos = pos[i];
            system_soft[i].vel = vel[i];
            system_soft[i].id = i_h + i + 1;
            system_soft[i].group_data.artificial.setParticleTypeToSingle();
        }

        file_header.nfile = 0;
        stat.time = file_header.time = 0.0;
        file_header.n_body = n_glb;

        stat.n_real_glb = stat.n_all_glb = n_glb;
        stat.n_real_loc = stat.n_all_loc = n_loc;

        read_data_flag = true;

        delete [] mass;
        delete [] pos;
        delete [] vel;
    }

    //! initial the system parameters
    void initialParameters() {
        // ensure data is read
        assert(read_data_flag);
        assert(stat.n_real_glb>0);

        bool print_flag = input_parameters.print_flag;
        int write_style = input_parameters.write_style.value;

#ifdef PROFILE
        if(print_flag) {
            std::cout<<"----- Parallelization information -----\n";
            std::cout<<"MPI processors: "<<n_proc<<std::endl;
            std::cout<<"OMP threads:    "<<PS::Comm::getNumberOfThread()<<std::endl;
        }
    
        // open profile file
        if(write_style>0) {
            std::string rank_str;
            std::stringstream atmp;
            atmp<<my_rank;
            atmp>>rank_str;
            std::string fproname=input_parameters.fname_snp.value+".prof.rank."+rank_str;
            if(input_parameters.app_flag) fprofile.open(fproname.c_str(),std::ofstream::out|std::ofstream::app);
            else  {
                fprofile.open(fproname.c_str(),std::ofstream::out);

                fprofile<<std::setprecision(WRITE_PRECISION);
                fprofile<<std::setw(WRITE_WIDTH)<<"my_rank"
                        <<std::setw(WRITE_WIDTH)<<"Time"
                        <<std::setw(WRITE_WIDTH)<<"N_steps"
                        <<std::setw(WRITE_WIDTH)<<"N_real_loc";
                profile.dumpName(fprofile, WRITE_WIDTH);
                profile.dumpBarrierName(fprofile, WRITE_WIDTH);
                tree_soft_profile.dumpName(fprofile, WRITE_WIDTH);
                tree_nb_profile.dumpName(fprofile, WRITE_WIDTH);
#if defined(USE_GPU) && defined(GPU_PROFILE)
                gpu_profile.dumpName(fprofile, WRITE_WIDTH);
                gpu_counter.dumpName(fprofile, WRITE_WIDTH);
#endif
                n_count.dumpName(fprofile, WRITE_WIDTH);
                fprofile<<std::endl;
            }
        }
#endif    


        // calculate system parameters
        PS::F64 r_in, m_average, m_max, v_disp, v_max;
        PS::F64& r_out = input_parameters.r_out.value;
        PS::F64& r_bin = input_parameters.r_bin.value;
        PS::F64& r_search_min = input_parameters.r_search_min.value;
        PS::F64& r_search_max = input_parameters.r_search_max.value;
        PS::F64& dt_soft = input_parameters.dt_soft.value;
        PS::F64& search_vel_factor =  input_parameters.search_vel_factor.value;
        PS::F64& ratio_r_cut   =  input_parameters.ratio_r_cut.value;
        PS::S64& n_bin         =  input_parameters.n_bin.value;
        PS::F64& theta         =  input_parameters.theta.value;

        getInitPar(system_soft, r_in, r_out, r_bin, r_search_min, r_search_max, v_max, m_average, m_max, dt_soft, v_disp, search_vel_factor, ratio_r_cut, n_bin, theta);

        EPISoft::eps   = input_parameters.eps.value;
        EPISoft::r_out = r_out;
        ForceSoft::grav_const = input_parameters.gravitational_constant.value;
        Ptcl::search_factor = search_vel_factor;
        Ptcl::r_search_min = r_search_min;
        Ptcl::mean_mass_inv = 1.0/m_average;
        Ptcl::r_group_crit_ratio = r_bin/r_in;

        if(print_flag) {
            std::cout<<"----- Parameter list: -----\n";
            std::cout<<" m_average    = "<<m_average      <<std::endl
                     <<" r_in         = "<<r_in           <<std::endl
                     <<" r_out        = "<<r_out          <<std::endl
                     <<" r_bin        = "<<r_bin          <<std::endl
                     <<" r_search_min = "<<r_search_min   <<std::endl
                     <<" vel_disp     = "<<v_disp         <<std::endl
                     <<" dt_soft      = "<<dt_soft        <<std::endl;
        }

        // check restart
        bool restart_flag = file_header.nfile; // nfile = 0 is assumed as initial data file

        // open output files
        // status information output
        std::string& fname_snp = input_parameters.fname_snp.value;
        if(write_style>0&&my_rank==0) {
            if(input_parameters.app_flag) 
                fstatus.open((fname_snp+".status").c_str(),std::ofstream::out|std::ofstream::app);
            else 
                fstatus.open((fname_snp+".status").c_str(),std::ofstream::out);
            fstatus<<std::setprecision(WRITE_PRECISION);

            // write titles of columns
            if(!restart_flag&&input_parameters.app_flag==false)  {
                stat.printColumnTitle(fstatus,WRITE_WIDTH);
                if (write_style==2) {
                    for (int i=0; i<stat.n_real_loc; i++) system_soft[0].printColumnTitle(fstatus, WRITE_WIDTH);
                }
                fstatus<<std::endl;
            }
        }

        // set id_offset
#ifdef PETAR_DEBUG
        assert(stat.n_real_glb == input_parameters.n_glb.value);
        assert(stat.n_real_glb == system_soft.getNumberOfParticleGlobal());
#endif
        PS::S64& id_offset = input_parameters.id_offset.value;
        id_offset = id_offset==-1 ? stat.n_real_glb+1 : id_offset;

        // initial particles paramters
        if (!restart_flag) {
#pragma omp parallel for
            for (PS::S32 i=0; i<stat.n_real_loc; i++) {
                // ID safety check 
                PS::S64 id = system_soft[i].id;
                assert(id_offset>id);
                assert(id>0);

                // set changeover 
                PS::F64 m_fac = system_soft[i].mass*Ptcl::mean_mass_inv;
                system_soft[i].changeover.setR(m_fac, r_in, r_out);

                // calculate r_search for particles, for binary, r_search depend on v_disp
                if(id<=2*n_bin) system_soft[i].r_search = std::max(r_search_min,v_disp*dt_soft*search_vel_factor + system_soft[i].changeover.getRout());
                else system_soft[i].calcRSearch(dt_soft);
            }
        }
        else {
            // clear up group_data.cm to avoid issue in search neighbor
#pragma omp parallel for
            for (PS::S32 i=0; i<stat.n_real_loc; i++) {
                auto& pi_cm = system_soft[i].group_data.cm;
                pi_cm.mass = pi_cm.vel.x = pi_cm.vel.y = pi_cm.vel.z = 0.0;
            }
        }
    
        // set system hard paramters
        hard_manager.setDtRange(input_parameters.dt_soft.value/input_parameters.dt_limit_hard_factor.value, input_parameters.dt_min_hermite_index.value);
        hard_manager.setEpsSq(input_parameters.eps.value);
        hard_manager.setGravitationalConstant(input_parameters.gravitational_constant.value);
        hard_manager.r_in_base = r_in;
        hard_manager.r_out_base = r_out;
#ifdef HARD_CHECK_ENERGY
        hard_manager.energy_error_max = input_parameters.e_err_hard.value;
#else
        hard_manager.energy_error_max = PS::LARGE_FLOAT;
#endif
        hard_manager.n_step_per_orbit = input_parameters.n_step_per_orbit.value;
        hard_manager.ap_manager.r_tidal_tensor = r_bin;
        hard_manager.ap_manager.id_offset = id_offset;
        hard_manager.ap_manager.setParticleSplitN(input_parameters.n_split.value);
        hard_manager.h4_manager.step.eta_4th = input_parameters.eta.value;
        hard_manager.h4_manager.step.eta_2nd = 0.01*input_parameters.eta.value;
        hard_manager.h4_manager.step.calcAcc0OffsetSq(m_average, r_out);
        hard_manager.ar_manager.energy_error_relative_max = input_parameters.e_err_arc.value;
        hard_manager.ar_manager.step_count_max = input_parameters.step_limit_arc.value;
        hard_manager.ar_manager.step.initialSymplecticCofficients(-6);
        hard_manager.ar_manager.slowdown_pert_ratio_ref = input_parameters.sd_factor.value;
        hard_manager.ar_manager.slowdown_timescale_max = dt_soft*input_parameters.n_step_per_orbit.value;
        //hard_manager.ar_manager.slowdown_timescale_max = dt_soft;
#ifdef SLOWDOWN_MASSRATIO
        hard_manager.ar_manager.slowdown_mass_ref = m_average;
#endif
        hard_manager.ar_manager.interrupt_detection_option = input_parameters.interrupt_detection_option.value;

        // check consistence of paramters
        input_parameters.checkParams();
        hard_manager.checkParams();

        // initial hard class and parameters
        system_hard_one_cluster.manager = &hard_manager;
        system_hard_one_cluster.setTimeOrigin(stat.time);

        system_hard_isolated.allocateHardIntegrator(input_parameters.n_interrupt_limit.value);
        system_hard_isolated.manager = &hard_manager;
        system_hard_isolated.setTimeOrigin(stat.time);

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        system_hard_connected.allocateHardIntegrator(input_parameters.n_interrupt_limit.value);
        system_hard_connected.manager = &hard_manager;
        system_hard_connected.setTimeOrigin(stat.time);
#endif

        time_kick = stat.time;

        if(write_style>0&&my_rank==0) {
            // save initial parameters
            std::string& fname_par = input_parameters.fname_par.value;
            if (print_flag) std::cout<<"Save input parameters to file "<<fname_par<<std::endl;
            FILE* fpar_out;
            if( (fpar_out = fopen(fname_par.c_str(),"w")) == NULL) {
                fprintf(stderr,"Error: Cannot open file %s.\n", fname_par.c_str());
                abort();
            }
            input_parameters.input_par_store.writeAscii(fpar_out);
            fclose(fpar_out);

            // save hard paramters 
            std::string fhard_par = input_parameters.fname_par.value + ".hard";
            if (print_flag) std::cout<<"Save hard_manager parameters to file "<<fhard_par<<std::endl;
            if( (fpar_out = fopen(fhard_par.c_str(),"w")) == NULL) {
                fprintf(stderr,"Error: Cannot open file %s.\n", fhard_par.c_str());
                abort();
            }
            hard_manager.writeBinary(fpar_out);
            fclose(fpar_out);
        }

#ifdef HARD_DUMP
        // initial hard_dump 
        const PS::S32 num_thread = PS::Comm::getNumberOfThread();
        hard_dump.initial(num_thread);
#endif

        // domain decomposition
        system_soft.setAverageTargetNumberOfSampleParticlePerProcess(input_parameters.n_smp_ave.value);
        const PS::F32 coef_ema = 0.2;
        domain_decompose_weight=1.0;
        dinfo.initialize(coef_ema);
        dinfo.decomposeDomainAll(system_soft,domain_decompose_weight);

        if(pos_domain!=NULL) {
            pos_domain = new PS::F64ort[n_proc];
            for(PS::S32 i=0; i<n_proc; i++) pos_domain[i] = dinfo.getPosDomain(i);
        }

        // exchange particles
        exchangeParticle();

        // tree for neighbor search
        tree_nb.initialize(stat.n_real_glb, input_parameters.theta.value, input_parameters.n_leaf_limit.value, input_parameters.n_group_limit.value);

        // tree for force
        PS::S64 n_tree_init = stat.n_real_glb + input_parameters.n_bin.value;
        tree_soft.initialize(n_tree_init, input_parameters.theta.value, input_parameters.n_leaf_limit.value, input_parameters.n_group_limit.value);

        // initial search cluster
        search_cluster.initialize();

        // initial tree step manager
        dt_manager.setKDKMode();

        initial_parameters_flag = true;
    }

    //! the first initial step to find groups and energy calculation
    void initialStep() {
        assert(initial_parameters_flag);
        if (initial_step_flag) return;

        assert(checkTimeConsistence());

        PS::F64 dt_tree = input_parameters.dt_soft.value;
        dt_manager.setStep(dt_tree);

        // >1. Tree for neighbor searching 
        /// get neighbor list to tree_nb
        treeNeighborSearch();

        // >2. search clusters
        /// gether clusters information to search_cluster, using tree_nb and velocity criterion (particles status/mass_bk)
        Ptcl::group_data_mode = GroupDataMode::cm;
        searchCluster();

        // >3. find group and create artificial particles
        /// find group and create artificial particles, using search_cluster, save to system_hard and system_soft (particle status/mass_bk updated)
        createGroup(dt_tree);

        // >4 tree soft force
        /// calculate tree force with linear cutoff, save to system_soft.acc
        treeSoftForce() ;

        // >5 correct change over
        /// correct system_soft.acc with changeover, using system_hard and system_soft particles
        treeForceCorrectChangeover();

        // correct force due to the change over update
        correctForceChangeOverUpdate();

        // recover mass
        kickClusterAndRecoverGroupMemberMass(system_soft, system_hard_isolated.getPtcl(), 0);

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        // connected
        kickClusterAndRecoverGroupMemberMass(system_soft, system_hard_connected.getPtcl(), 0);
        // sending list for connected clusters
        kickSend(system_soft, search_cluster.getAdrSysConnectClusterSend(), 0);
        // send kicked particle from sending list, and receive remote single particle
        search_cluster.SendSinglePtcl(system_soft, system_hard_connected.getPtcl());
#endif

        // update global particle system due to kick
        writeBackHardParticles();

        // initial status and energy
        updateStatus(true);

        // output initial data
        output();

        // remove artificial particles
        system_soft.setNumberOfParticleLocal(stat.n_real_loc);

#ifdef CLUSTER_VELOCITY
        setParticleGroupDataToCMData();
#endif

#ifdef PROFILE
        clearProfile();
#endif
        initial_step_flag = true;
    }

    
    //! integrate the system
    /*! @param[in] _time_break: additional breaking time to interrupt the integration, in default (0.0) the system integrate to time_end
      \return interrupted cluster number
     */
    PS::S32 evolveToTime(const PS::F64 _time_break=0.0) {

        // ensure it is initialized
        assert(initial_step_flag);

        // finish interrupted integrations
        if (n_interrupt_glb>0) finishInterruptDrift();
        // if interrupt still exist, do not continue
        if (n_interrupt_glb>0) return n_interrupt_glb;

        // check time break
        PS::F64 time_break = _time_break==0.0? input_parameters.time_end.value: std::min(_time_break,input_parameters.time_end.value);
        if (stat.time>=time_break) return 0;

        PS::F64 dt_output = input_parameters.dt_snp.value;
        PS::F64 dt_tree = dt_manager.getStep();

        /// Main loop
        while(stat.time <= time_break) {
#ifdef PROFILE
            profile.total.start();
#endif

            // remove deletec particles in system_soft.
            removeParticles();

            // >9. Domain decomposition
            domainDecompose();

            // >10. exchange particles
            exchangeParticle();

            // >1. Tree for neighbor searching 
            /// get neighbor list to tree_nb
            treeNeighborSearch();

            // >2. search clusters
            /// gether clusters information to search_cluster, using tree_nb and velocity criterion (particles status/mass_bk)
            searchCluster();

            // >3. find group and create artificial particles
            /// find group and create artificial particles, using search_cluster, save to system_hard and system_soft (particle status/mass_bk updated)
            createGroup(dt_tree);

            // >4 tree soft force
            /// calculate tree force with linear cutoff, save to system_soft.acc
            treeSoftForce() ;

            // >5 correct change over
            /// correct system_soft.acc with changeover, using system_hard and system_soft particles
            treeForceCorrectChangeover();

#ifdef KDKDK_4TH
            // only do correction at middle step
            if (dt_manager.getCountContinue() == 1) GradientKick();
#endif

            bool interrupt_flag = false;  // for interrupt integration when time reach end
            bool output_flag = false;    // for output snapshot and information
            //bool dt_mod_flag = false;    // for check whether tree time step need update
            bool changeover_flag = false; // for check whether changeover need update
            PS::F64 dt_kick, dt_drift;

            // for initial the system
            if (dt_manager.isNextStart()) {

                // update changeover if last time it is modified.
                correctForceChangeOverUpdate();

                // set step to the begining step
                dt_kick = dt_manager.getDtStartContinue();
            }
            else {
#ifdef PROFILE
                profile.other.start();
#endif
                // increase loop counter
                n_loop++;

                // for next kick-drift pair
                dt_manager.nextContinue();

                // update stat time 
                stat.time = system_hard_one_cluster.getTimeOrigin();

                // check whether output or changeover change are needed (only at the ending step)
                if (dt_manager.isNextEndPossible()) {

                    // adjust tree step
                    //dt_mod_flag = adjustDtTreeReduce(dt_reduce_factor, dt_tree, dt_manager.getStep());
                    //if(dt_mod_flag) dt_kick = dt_manager.getDtEndContinue();

                    // output step, get last kick step
                    output_flag = (fmod(stat.time, dt_output) == 0.0);

                    // check changeover change
                    changeover_flag = (system_hard_isolated.getNClusterChangeOverUpdate()>0);

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
                    PS::S32 n_changeover_modify_local  = system_hard_connected.getNClusterChangeOverUpdate() + system_hard_isolated.getNClusterChangeOverUpdate();
                    PS::S32 n_changeover_modify_global = PS::Comm::getSum(n_changeover_modify_local);
                    if (n_changeover_modify_global>0) changeover_flag = true;
#endif

                    // check interruption
                    interrupt_flag = (stat.time>=time_break);

                    // set next step to be last
                    if (output_flag||changeover_flag||interrupt_flag) dt_kick = dt_manager.getDtEndContinue();
                    else dt_kick = dt_manager.getDtKickContinue();
                }
                else dt_kick = dt_manager.getDtKickContinue();
#ifdef PROFILE
                profile.other.barrier();
                PS::Comm::barrier();
                profile.other.end();
#endif
            }


            // >6. kick 
            kick(dt_kick);

            // >7. write back data
            if(output_flag||interrupt_flag) {
                // update global particle system due to kick
                writeBackHardParticles();
            }

            // output information
            if(output_flag) {
                // update status
                updateStatus(false);
                output();

#ifdef PROFILE
                // profile
                profile.total.barrier();
                PS::Comm::barrier();
                profile.total.end();

                printProfile();
                clearProfile();

                PS::Comm::barrier();
                profile.total.start();
#endif
            }

            // modify the tree step
            //if(dt_mod_flag) {
            //    dt_manager.setStep(dt_tree/dt_reduce_factor);
            //    if (input_parameters.print_flag) {
            //        std::cout<<"Tree time step change, time = "<<stat.time
            //                 <<"  dt_soft = "<<dt_tree
            //                 <<"  reduce factor = "<<dt_reduce_factor
            //                 <<std::endl;
            //    }
            //}

            // interrupt
            if(interrupt_flag) {
#ifdef CLUSTER_VELOCITY
                setParticleGroupDataToCMData();
#endif
                // correct force due to the change over update
                correctForceChangeOverUpdate();

                // need remove artificial particles
                system_soft.setNumberOfParticleLocal(stat.n_real_loc);

                assert(checkTimeConsistence());

#ifdef PROFILE
                profile.total.barrier();
                PS::Comm::barrier();
                profile.total.end();
#endif
                return 0;
            }

            // second kick if output exists or changeover is modified
            //if(dt_mod_flag||output_flag||changeover_flag) {
            if(output_flag||changeover_flag) {

                assert(checkTimeConsistence());

                // determine the next tree time step
                //PS::F64 dt_reduce_factor_org = 1.0;
                //dt_reduce_factor_org = PS::Comm::getMaxValue(dt_reduce_factor_org);
                //dt_reduce_factor = 1.0;
                //while(dt_reduce_factor<dt_reduce_factor_org) dt_reduce_factor *=2.0;

                //update new tree step if reduce factor is changed
                dt_kick = dt_manager.getDtStartContinue();

                correctForceChangeOverUpdate();

                kick(dt_kick);
            }


            // >8. Hard integration 
            // get drift step
            dt_drift = dt_manager.getDtDriftContinue();

            drift(dt_drift);

#ifdef PROFILE
            // calculate profile
            profile.total.barrier();
            PS::Comm::barrier();
            profile.total.end();

            calcProfile();
#endif
            
            // when interrupt exist, quit the loop
            if (n_interrupt_glb>0) {
#ifdef HARD_INTERRUPT_PRINT
                std::cerr<<"Interrupt detected, number: "<<n_interrupt_glb<<std::endl;
#endif
                return n_interrupt_glb;
            }

        }

        return 0;
    }

    void clear() {

        if (fstatus.is_open()) fstatus.close();
#ifdef PROFILE
        if (fprofile.is_open()) fprofile.close();
#endif
        if (pos_domain) {
            delete[] pos_domain;
            pos_domain=NULL;
        }

        if (initial_fdps_flag) PS::Finalize();
        remove_list.resizeNoInitialize(0);
        n_interrupt_glb = 0;
        initial_fdps_flag = false;
        read_parameters_flag = false;
        read_data_flag = false;
        initial_parameters_flag = false;
        initial_step_flag = false;
    }

    ~PeTar() { clear();}
};

