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

#include<get_version.hpp>
#include<particle_simulator.hpp>
#include"hard_assert.hpp"
#include"soft_ptcl.hpp"
#include"soft_force.hpp"
#include"astro_units.hpp"
#ifdef USE_GPU
#include"force_gpu_cuda.hpp"
#endif
#ifdef USE_FUGAKU
#include "force_fugaku.hpp"
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
#include"escaper.hpp"
#ifdef GALPY
#include"galpy_interface.h"
#endif
#ifdef BSE_BASE
#include"rand_interface.hpp"
#endif

//! IO parameters for Petar
class IOParamsPeTar{
public:
    // IO parameters
    IOParamsContainer input_par_store;
    IOParams<PS::F64> ratio_r_cut;
    IOParams<PS::F64> theta;
    IOParams<PS::S64> n_leaf_limit;
    IOParams<PS::S64> n_group_limit;
//    IOParams<PS::S64> n_interrupt_limit;
    IOParams<PS::S64> n_smp_ave;
#ifdef ORBIT_SAMPLING
    IOParams<PS::S64> n_split;
#endif
    IOParams<PS::S64> n_bin;
    IOParams<PS::S64> n_step_per_orbit;
    IOParams<PS::F64> time_end;
    IOParams<PS::F64> eta;
    IOParams<PS::F64> gravitational_constant;
    IOParams<PS::S64> unit_set;
    IOParams<PS::S64> n_glb;
    IOParams<PS::S64> id_offset;
    IOParams<PS::F64> dt_soft;
    IOParams<PS::F64> dt_snap;
    IOParams<PS::F64> nstep_dt_soft_kepler;
    IOParams<PS::F64> search_vel_factor;
    IOParams<PS::F64> search_peri_factor;
    IOParams<PS::F64> dt_limit_hard_factor;
    IOParams<PS::S64> dt_min_hermite_index;
    //IOParams<PS::S64> dt_min_ar_index;
    //IOParams<PS::F64> dt_err_pert;
    //IOParams<PS::F64> dt_err_soft;
    IOParams<PS::F64> e_err_ar;
#ifdef HARD_CHECK_ENERGY
    IOParams<PS::F64> e_err_hard;
#endif
    IOParams<PS::S64> step_limit_ar;
    IOParams<PS::F64> eps;
    IOParams<PS::F64> r_out;
    IOParams<PS::F64> r_bin;
//    IOParams<PS::F64> r_search_max;
    IOParams<PS::F64> r_search_min;
    IOParams<PS::F64> r_escape;
    IOParams<PS::F64> sd_factor;
    IOParams<PS::S64> data_format;
    IOParams<PS::S64> write_style;
#ifdef STELLAR_EVOLUTION
#ifdef BSE_BASE
    IOParams<PS::S64> stellar_evolution_option;
#endif
    IOParams<PS::S64> interrupt_detection_option;
#endif
#ifdef ADJUST_GROUP_PRINT
    IOParams<PS::S64> adjust_group_write_option;
#endif
    IOParams<PS::S64> append_switcher;
    IOParams<PS::S64> record_id_start_one;
    IOParams<PS::S64> record_id_end_one;
    IOParams<PS::S64> record_id_start_two;
    IOParams<PS::S64> record_id_end_two;
    IOParams<std::string> fname_snp;
    IOParams<std::string> fname_par;
    IOParams<std::string> fname_inp;

    // flag
    bool print_flag; 
    bool update_changeover_flag;
    bool update_rsearch_flag;

    IOParamsPeTar(): input_par_store(), 
                     ratio_r_cut      (input_par_store, 0.1,  "r-ratio", "r_in / r_out"),
                     theta            (input_par_store, 0.3,  "T",  "Particle-tree opening angle theta"),
                     n_leaf_limit     (input_par_store, 20,   "number-leaf-limit", "Particle-tree leaf number limit", "Optimal value should be slightly >= 11 + N_bin_sample (20)"),
#ifdef USE__AVX512
                     n_group_limit    (input_par_store, 1024, "number-group-limit", "Particle-tree group number limit", "Optimized for x86-AVX512 (1024)"),    
#else
                     n_group_limit    (input_par_store, 512,  "number-group-limit", "Particle-tree group number limit", "Optimized for x86-AVX2 (512)"),
#endif
//                     n_interrupt_limit(input_par_store, 128,  "number-interrupt-limit", "Interrupted hard integrator limit"),
                     n_smp_ave        (input_par_store, 100,  "number-sample-average", "Average target number of sample particles per process"),
#ifdef ORBIT_SAMPLING
                     n_split          (input_par_store, 4,    "number-split", "Number of binary sample points for tree perturbation force"),
#endif
                     n_bin            (input_par_store, 0,    "b", "Number of primordial binaries for initialization (assuming the binaries' IDs are 1,2*n_bin)"),
                     n_step_per_orbit (input_par_store, 8,    "number-step-tt", "Number of steps per slow-down binary orbits (binary period/tree timestep) for isolated binaries; also the maximum criterion for activating tidal tensor method"),
                     time_end         (input_par_store, 10.0, "t", "End time of simulation"),
                     eta              (input_par_store, 0.1,  "hermite-eta", "Hermite timestep coefficient eta"),
                     gravitational_constant(input_par_store, 1.0, "G", "Gravitational constant"),
                     unit_set         (input_par_store, 0,    "u", "Input data unit, 0: unknown, referring to G; 1: mass:Msun, length:pc, time:Myr, velocity:pc/Myr"),
                     n_glb            (input_par_store, 100000, "n", "Total number of particles, used only for a test with the internal equal-mass Plummer model generator (assuming G=1 and the input data filename is __Plummer)"),
                     id_offset        (input_par_store, -1,   "id-offset", "Starting ID for artificial particles, total number of real particles must always be smaller than this","n_glb+1"),
                     dt_soft          (input_par_store, 0.0,  "s", "Tree timestep (dt_soft), if the value is zero (default) and --nstep-dt-soft-kepler is not used, then dt_soft = 0.1*r_out/sigma_1D"),
                     dt_snap          (input_par_store, 1.0,  "o", "Output time interval for particle dataset snapshots"),
                     nstep_dt_soft_kepler (input_par_store, 0.0, "nstep-dt-soft-kepler", "Determines the tree timestep by P(r_in)/nstep, where P(r_in) is the binary period with the semi-major axis of r_in, nstep is the argument of this option (e.g., 32.0)", "not used"),
                     search_vel_factor(input_par_store, 3.0,  "search-vel-factor", "Neighbor search coefficient for velocity check (v*dt)"),
                     search_peri_factor  (input_par_store, 1.5, "search-peri-factor", "Neighbor search coefficient for periapsis check"),
                     dt_limit_hard_factor(input_par_store, 4.0, "dt-max-factor", "Limit of tree timestep/hard timestep"),
                     dt_min_hermite_index(input_par_store, 40,  "dt-min-hermite",  "Power index n for the smallest timestep (0.5^n) allowed in the Hermite integrator"),
                     //dt_min_ar_index     (input_par_store, 64,  "dt-min-ar",  "Power index n for the smallest timestep (0.5^n) allowed in the ARC integrator, suppressed"),
                     //dt_err_pert  (input_par_store, 1e-6, "dt-error-pert", "Maximum time synchronization error (relative) for perturbed ARC integrator, suppressed"),
                     //dt_err_soft  (input_par_store, 1e-3, "dt-error-iso", "Maximum time synchronization error (relative) for no-perturber (only soft perturbation) ARC integrator, suppressed"),
                     e_err_ar     (input_par_store, 1e-8, "energy-err-ar", "Maximum energy error allowed for the ARC integrator"),
#ifdef HARD_CHECK_ENERGY
                     e_err_hard   (input_par_store, 1e-4, "energy-err-hard", "Maximum energy error allowed for the hard integrator"),
#endif
                     step_limit_ar(input_par_store, 1000000, "step-limit-ar", "Maximum step allowed for the ARC sym integrator"),
                     eps          (input_par_store, 0.0,  "soft-eps", "Softening epsilon"),
                     r_out        (input_par_store, 0.0,  "r", "Outer boundary radius for the changeover function (r_out), if value is zero and -s is not used, use 0.1 GM/[N^(1/3) sigma_3D^2]; if -s is given, calculate r_out from dt_soft"),
                     r_bin        (input_par_store, 0.0,  "r-bin", "Tidal tensor box size and the radial criterion for detecting multiple systems (binaries, triples, etc.), if value is zero, use 0.8*r_in"),
//                     r_search_max (input_par_store, 0.0,  "Maximum search radius criterion", "5*r_out"),
                     r_search_min (input_par_store, 0.0,  "r-search-min", "Minimum neighbor search radius for hard clusters","auto"),
                     r_escape     (input_par_store, PS::LARGE_FLOAT,  "r-escape", "Escape radius criterion, 0: no escaper removal; <0: remove particles when r>-r_escape; >0: remove particles when r>r_escape and energy>0"),
                     sd_factor    (input_par_store, 1e-4, "slowdown-factor", "Slowdown perturbation criterion"),
                     data_format  (input_par_store, 1,    "i", "Data read(r)/write(w) format BINARY(B)/ASCII(A): r-B/w-A (3), r-A/w-B (2), rw-A (1), rw-B (0)"),
                     write_style  (input_par_store, 1,    "w", "File writing style: 0, no output; 1. write snapshots, status, and profile separately; 2. write snapshot and status in one line per step (no MPI support); 3. write only status and profile"),
#ifdef STELLAR_EVOLUTION
#ifdef BSE_BASE
                     stellar_evolution_option  (input_par_store, 1, "stellar-evolution", "Stellar evolution of stars in Hermite+SDAR: 0: off; >=1: using SSE/BSE based codes; ==2: activate dynamical tide and hyperbolic gravitational wave radiation"),
                     interrupt_detection_option(input_par_store, 1, "detect-interrupt", "Stellar evolution of binaries in SDAR: 0: off; 1: using BSE based code (if '--stellar-evolution != 0)"),
#else // NO BSE_BASE
                     interrupt_detection_option(input_par_store, 0, "detect-interrupt", "Interrupt integration in SDAR: 0: turn off; 1: merge two particles if their surfaces overlap; 2. record two particle information if their surfaces overlap without merger"),
#endif // END BSE_BASE
#endif // END STELLAR_EVOLUTION

#ifdef ADJUST_GROUP_PRINT
                     adjust_group_write_option(input_par_store, 1, "write-group-info", "Print new and end of groups: 0: no print; 1: print to file [data filename prefix].group.[MPI rank] if -w >0"),
#endif
                     append_switcher(input_par_store, 1, "a", "Data output style: 0 - create new output files and overwrite existing ones except snapshots; 1 - append new data to existing files"),
                     record_id_start_one(input_par_store, 0, "record-id-start-one", "Starting of the first id range for hard dump recording every tree step"),
                     record_id_end_one(input_par_store, 0, "record-id-end-one", "Ending of the first id range for hard dump; notice that the ending id is not included in hard dump"),
                     record_id_start_two(input_par_store, 0, "record-id-start-two", "Starting of the 2nd id range for hard dump recording every tree step"),
                     record_id_end_two(input_par_store, 0, "record-id-end-two", "Ending of the 2nd id range for hard dump; notice that the ending id is not included in hard dump"),
                     fname_snp(input_par_store, "data", "f", "Prefix of filenames for output data: [prefix].**"),
                     fname_par(input_par_store, "input.par", "p", "Input parameter file (this option should be used first before any other options)"),
                     fname_inp(input_par_store, "__NONE__", "snap-filename", "Input data file", NULL, false),
                     print_flag(false), update_changeover_flag(false), update_rsearch_flag(false) {}

    
    //! reading parameters from GNU option API
    /*!
      @param[in] argc: number of options
      @param[in] argv: string of options
      @param[in] opt_used_pre: already used option number from previous reading, use to correctly count the remaining argument number
      \return -1 if help is used; else the used number of argv
     */
    int read(int argc, char *argv[], const int opt_used_pre=0) {
        static int petar_flag=-1;
        static struct option long_options[] = {
#ifdef ORBIT_SAMPLING
            {n_split.key,              required_argument, &petar_flag, 0},        
#endif
            {search_vel_factor.key,    required_argument, &petar_flag, 1},  
            {dt_limit_hard_factor.key, required_argument, &petar_flag, 2},  
            {dt_min_hermite_index.key, required_argument, &petar_flag, 3}, 
            {n_group_limit.key,        required_argument, &petar_flag, 4},
            {n_leaf_limit.key,         required_argument, &petar_flag, 5},
            {n_smp_ave.key,            required_argument, &petar_flag, 6},
            {e_err_ar.key,             required_argument, &petar_flag, 7}, 
            {eps.key,                  required_argument, &petar_flag, 8},       
            {sd_factor.key,            required_argument, &petar_flag, 9},
            {ratio_r_cut.key,          required_argument, &petar_flag, 10},
            {r_bin.key,                required_argument, &petar_flag, 11},       
            {search_peri_factor.key,   required_argument, &petar_flag, 12}, 
            {eta.key,                  required_argument, &petar_flag, 13}, 
#ifdef HARD_CHECK_ENERGY
            {e_err_hard.key,           required_argument, &petar_flag, 14},  
#endif
            {step_limit_ar.key,        required_argument, &petar_flag, 15},   
            {"disable-print-info",     no_argument,       &petar_flag, 16},
            {n_step_per_orbit.key,     required_argument, &petar_flag, 17},
#ifdef STELLAR_EVOLUTION
//            {n_interrupt_limit.key,    required_argument, &petar_flag, 18},
#ifdef BSE_BASE
            {stellar_evolution_option.key,    required_argument, &petar_flag, 19},
#endif
            {interrupt_detection_option.key,  required_argument, &petar_flag, 20},
#endif
            {r_escape.key,             required_argument, &petar_flag, 21},
            {r_search_min.key,         required_argument, &petar_flag, 22},
            {id_offset.key,            required_argument, &petar_flag, 23},
#ifdef ADJUST_GROUP_PRINT
            {adjust_group_write_option.key,   required_argument, &petar_flag, 24},
#endif            
            {nstep_dt_soft_kepler.key,  required_argument, &petar_flag, 25},
            {record_id_start_one.key,  required_argument, &petar_flag, 26},
            {record_id_end_one.key,    required_argument, &petar_flag, 27},
            {record_id_start_two.key,  required_argument, &petar_flag, 28},
            {record_id_end_two.key,    required_argument, &petar_flag, 29},
            {"help",                  no_argument, 0, 'h'},        
            {0,0,0,0}
        };

        int opt_used = opt_used_pre;
        int copt;
        int option_index;
        optind = 0; // reset getopt
        while ((copt = getopt_long(argc, argv, "-i:a:t:s:o:r:b:n:G:u:T:f:p:w:h", long_options, &option_index)) != -1) 
            switch (copt) {
            case 0:
                switch (petar_flag) {
#ifdef ORBIT_SAMPLING
                case 0:
                    n_split.value = atoi(optarg);
                    if(print_flag) n_split.print(std::cout);
                    opt_used += 2;
                    assert(n_split.value>=0);
                    break;
#endif
                case 1:
                    search_vel_factor.value = atof(optarg);
                    if(print_flag) search_vel_factor.print(std::cout);
                    opt_used += 2;
                    update_rsearch_flag = true;
                    assert(search_vel_factor.value>0.0);
                    break;
                case 2:
                    dt_limit_hard_factor.value = atof(optarg);
                    if(print_flag) dt_limit_hard_factor.print(std::cout);
                    opt_used += 2;
                    assert(dt_limit_hard_factor.value > 0.0);
                    break;
                case 3:
                    dt_min_hermite_index.value = atoi(optarg);
                    if(print_flag) dt_min_hermite_index.print(std::cout);
                    opt_used += 2;
                    assert(dt_min_hermite_index.value > 0);
                    break;
                case 4:
                    n_group_limit.value = atoi(optarg);
                    if(print_flag) n_group_limit.print(std::cout);
                    opt_used += 2;
                    assert(n_group_limit.value>0);
                    break;
                case 5:
                    n_leaf_limit.value = atoi(optarg);
                    if(print_flag) n_leaf_limit.print(std::cout);
                    opt_used += 2;
                    assert(n_leaf_limit.value>0);
                    break;
                case 6:
                    n_smp_ave.value = atoi(optarg);
                    if(print_flag) n_smp_ave.print(std::cout);
                    opt_used += 2;
                    assert(n_smp_ave.value>0.0);
                    break;
                case 7:
                    e_err_ar.value = atof(optarg);
                    if(print_flag) e_err_ar.print(std::cout);
                    opt_used += 2;
                    assert(e_err_ar.value > 0.0);
                    break;
                case 8:
                    eps.value = atof(optarg);
                    if(print_flag) eps.print(std::cout);
                    opt_used += 2;
                    assert(eps.value>=0.0);
                    break;
                case 9:
                    sd_factor.value = atof(optarg);
                    if(print_flag) sd_factor.print(std::cout);
                    opt_used += 2;
                    assert(sd_factor.value>0.0);
                    break;
                case 10:
                    ratio_r_cut.value = atof(optarg);
                    if(print_flag) ratio_r_cut.print(std::cout);
                    update_changeover_flag = true;
                    opt_used += 2;
                    assert(ratio_r_cut.value>0.0);
                    assert(ratio_r_cut.value<1.0);
                    break;
                case 11:
                    r_bin.value = atof(optarg);
                    if(print_flag) r_bin.print(std::cout);
                    opt_used += 2;
                    assert(r_bin.value>=0.0);
                    break;
                case 12:
                    search_peri_factor.value = atof(optarg);
                    if(print_flag) search_peri_factor.print(std::cout);
                    opt_used += 2;
                    assert(search_peri_factor.value>=1.0);
                    break;
                case 13:
                    eta.value = atof(optarg);
                    if(print_flag) eta.print(std::cout);
                    opt_used += 2;
                    assert(eta.value>0.0);
                    break;
#ifdef HARD_CHECK_ENERGY
                case 14:
                    e_err_hard.value = atof(optarg);
                    if(print_flag) e_err_hard.print(std::cout);
                    opt_used += 2;
                    break;
#endif
                case 15:
                    step_limit_ar.value = atoi(optarg);
                    if(print_flag) step_limit_ar.print(std::cout);
                    opt_used += 2;
                    break;
                case 16:
                    print_flag = false;
                    opt_used ++;
                    break;
                case 17:
                    n_step_per_orbit.value = atof(optarg);
                    if(print_flag) n_step_per_orbit.print(std::cout);
                    opt_used += 2;
                    assert(n_step_per_orbit.value>=1.0 || n_step_per_orbit.value==0.0);
                    break;
#ifdef STELLAR_EVOLUTION
//                case 18:
//                    n_interrupt_limit.value = atoi(optarg);
//                    if(print_flag) n_interrupt_limit.print(std::cout);
//                    opt_used += 2;
//                    assert(n_interrupt_limit.value>0);
//                    break;
#ifdef BSE_BASE
                case 19:
                    stellar_evolution_option.value = atoi(optarg);
                    if(print_flag) stellar_evolution_option.print(std::cout);
                    opt_used += 2;
                    break;
#endif
                case 20:
                    interrupt_detection_option.value = atoi(optarg);
                    if(print_flag) interrupt_detection_option.print(std::cout);
                    opt_used += 2;
                    break;
#endif
                case 21:
                    r_escape.value = atof(optarg);
                    if(print_flag) r_escape.print(std::cout);
                    opt_used += 2;
                    break;
                case 22:
                    r_search_min.value = atof(optarg);
                    if(print_flag) r_search_min.print(std::cout);
                    update_rsearch_flag = true;
                    opt_used += 2;
                    break;
                case 23:
                    id_offset.value = atoi(optarg);
                    if(print_flag) id_offset.print(std::cout);
                    opt_used += 2;
                    break;
#ifdef ADJUST_GROUP_PRINT
                case 24:
                    adjust_group_write_option.value = atoi(optarg);
                    if(print_flag) adjust_group_write_option.print(std::cout);
                    opt_used += 2;
                    break;
#endif
                case 25:
                    nstep_dt_soft_kepler.value = atof(optarg);
                    if(print_flag) nstep_dt_soft_kepler.print(std::cout);
                    opt_used += 2;
                    break;
                case 26:
                    record_id_start_one.value = atoi(optarg);
                    if(print_flag) record_id_start_one.print(std::cout);
                    opt_used += 2;
                    break;
                case 27:
                    record_id_end_one.value = atoi(optarg);
                    if(print_flag) record_id_end_one.print(std::cout);
                    opt_used += 2;
                    break;
                case 28:
                    record_id_start_two.value = atoi(optarg);
                    if(print_flag) record_id_start_two.print(std::cout);
                    opt_used += 2;
                    break;
                case 29:
                    record_id_end_two.value = atoi(optarg);
                    if(print_flag) record_id_end_two.print(std::cout);
                    opt_used += 2;
                    break;
                default:
                    break;
                }
                break;
            case 'i':
                data_format.value = atoi(optarg);
                if(print_flag) data_format.print(std::cout);
                assert(data_format.value>=0&&data_format.value<=3);
                opt_used += 2;
                break;
            case 'a':
                append_switcher.value=atoi(optarg);
                if(print_flag) append_switcher.print(std::cout);
                assert(append_switcher.value>=0&&append_switcher.value<=1);
                opt_used += 2;
                break;
            case 't':
                time_end.value = atof(optarg);
                if(print_flag) time_end.print(std::cout);
                opt_used += 2;
                assert(time_end.value>=0.0);
                break;
            case 's':
                dt_soft.value = atof(optarg);
                if(print_flag) dt_soft.print(std::cout);
                update_rsearch_flag = true;
                opt_used += 2;
                assert(dt_soft.value>=0.0);
                break;
            case 'o':
                dt_snap.value = atof(optarg);
                if(print_flag) dt_snap.print(std::cout);
                opt_used += 2;
                assert(dt_snap.value>0.0);
                break;
            case 'r':
                r_out.value = atof(optarg);
                if(print_flag) r_out.print(std::cout);
                update_rsearch_flag = true;
                update_changeover_flag = true;
                opt_used += 2;
                assert(r_out.value>=0.0);
                break;
            case 'b':
                n_bin.value = atoi(optarg);
                if(print_flag) n_bin.print(std::cout);
                opt_used += 2;
                assert(n_bin.value>=0);
                break;
            case 'n':
                n_glb.value = atol(optarg);
                if(print_flag) n_glb.print(std::cout);
                opt_used += 2;
                assert(n_glb.value>0);
                break;
            case 'G':
                gravitational_constant.value = atof(optarg);
                if(print_flag) gravitational_constant.print(std::cout);
                opt_used += 2;
                assert(gravitational_constant.value>0.0);
                break;
            case 'u':
                unit_set.value = atoi(optarg);
                if(print_flag) unit_set.print(std::cout);
                opt_used += 2;
                assert(unit_set.value>=0);
                break;
            case 'T':
                theta.value = atof(optarg);
                if(print_flag) theta.print(std::cout);
                opt_used += 2;
                assert(theta.value>=0.0);
                break;
            case 'f':
                fname_snp.value = optarg;
                if(print_flag) fname_snp.print(std::cout);
                opt_used += 2;
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
                opt_used += 2;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
                input_par_store.mpi_broadcast();
                PS::Comm::barrier();
#endif
                break;
            case 'w':
                write_style.value = atoi(optarg);
                if(print_flag) write_style.print(std::cout);
                opt_used += 2;
                break;
            case 'h':
                if(print_flag){
                    std::cout<<"Usage: petar [option] filename"<<std::endl;
                    std::cout<<"       filename: initial or restart data (snapshot) filename\n"
                             <<"       Data file content:\n"
                             <<"            First line: \n"
                             <<"             1. File_ID: 0 for initialization, else for restarting\n"
                             <<"             2. N_particle \n"
                             <<"             3. Time\n"
#ifdef RECORD_CM_IN_HEADER
                             <<"             4-6: offset of system centeral position (assuming galactic center is the coordiante origin in the case of Galpy)\n"
                             <<"             7-9: offset of system centeral velocity\n"
#endif
                             <<"            Following lines:\n";
                    FPSoft::printTitleWithMeaning(std::cout,0,13);
                    std::cout<<"          PS: (*) show initialization values which should be used together with FILE_ID = 0"<<std::endl;
                    std::cout<<"              [formatted] indicates that the value is only for save, cannot be directly read"<<std::endl;
                    std::cout<<"Options:\n";
                    input_par_store.printHelp(std::cout, 2, 10, 23);
                    std::cout<<"        --disable-print-info:  "<<"Do not print information"<<std::endl;
                    std::cout<<"        --disable-write-info:  "<<"Do not write information"<<std::endl;
                    std::cout<<"  -h(--help):               print help"<<std::endl;
                    std::cout<<"*** PS: dt_soft: tree time step\n"
                             <<"        r_in : transit function inner boundary radius\n"
                             <<"        r_out: transit function outer boundary radius\n"
                             <<"        sigma: half-mass radius velocity dispersion\n"
                             <<"        n_bin: number of primordial binaries\n"
                             <<"        <m>  : averaged mass"<<std::endl;
                }
                return -1;
            case '?':
                opt_used +=2;
                break;
            default:
                break;
            }

        // count used options
        opt_used ++;
        //std::cout<<"Opt used:"<<opt_used<<std::endl;
        if (opt_used<argc) {
            fname_inp.value =argv[argc-1];
            if(print_flag) std::cout<<"Reading data file name: "<<fname_inp.value<<std::endl;
        }

        if(print_flag) std::cout<<"----- Finish reading input options of PeTar -----\n";

        return opt_used-1;
    }

    //! check paramters
    bool checkParams() {
#ifdef ORBIT_SAMPLING
        assert(n_split.value>=0);
#endif
        assert(search_vel_factor.value>0.0);
        assert(dt_limit_hard_factor.value > 0.0);
        assert(dt_min_hermite_index.value > 0);
        assert(e_err_ar.value > 0.0);
        assert(eps.value>=0.0 && eps.value<=ratio_r_cut.value); // avoid incorrect self-potential correction after tree force, when eps>r_out, self-potential is G m /r_eps instead of G m/r_cut;
        assert(sd_factor.value>0.0);
        assert(ratio_r_cut.value>0.0);
        assert(ratio_r_cut.value<1.0);
        assert(r_bin.value>=0.0);
        assert(search_peri_factor.value>=1.0);
        assert(data_format.value>=0||data_format.value<=3);
        assert(time_end.value>=0.0);
        assert(dt_soft.value>=0.0);
        assert(dt_snap.value>0.0);
        assert(r_out.value>=0.0);
        assert(n_bin.value>=0);
        assert(n_glb.value>0);
        assert(n_group_limit.value>0);
//        assert(n_interrupt_limit.value>0);
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
#ifdef USE_QUAD
    typedef PS::TreeForForceLong<ForceSoft, EPISoft, EPJSoft>::QuadrupoleWithSymmetrySearch TreeForce; 
#else
    typedef PS::TreeForForceLong<ForceSoft, EPISoft, EPJSoft>::MonopoleWithSymmetrySearch TreeForce;
#endif
    typedef PS::ParticleSystem<FPSoft> SystemSoft;

    // For neighbor searching
    typedef PS::TreeForForceShort<ForceSoft, EPISoft, EPJSoft>::Symmetry TreeNB;

    // IO
    IOParamsPeTar input_parameters;
#ifdef BSE_BASE
    IOParamsBSE bse_parameters;
    std::string fbse_par_suffix = BSEManager::getBSEOutputFilenameSuffix();
    std::string fsse_par_suffix = BSEManager::getSSEOutputFilenameSuffix();
    IOParamsRand rand_parameters;
#endif // BSE_BASE
#ifdef GALPY
    IOParamsGalpy galpy_parameters;
#endif
#ifdef EXTERNAL_HARD
    IOParamsExternalHard external_hard_parameters;
#endif

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

    // escaper
    Escaper escaper;
    std::ofstream fesc;

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
#ifdef FDPS_COMM
    PS::CommInfo comm_info; // MPI communicator
#endif

    // tree time step manager
    KickDriftStep dt_manager;

    // tree
    TreeNB tree_nb;
    TreeForce tree_soft;

    // random manager
#ifdef BSE_BASE
    RandomManager rand_manager;
#endif

#ifdef GALPY
    GalpyManager galpy_manager;
#endif

    // hard integrator
    HardManager hard_manager;
    SystemHard system_hard_one_cluster;
    SystemHard system_hard_isolated;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    SystemHard system_hard_connected;
#endif
    int n_interrupt_glb;

    // mass change particle list
    PS::ReallocatableArray<PS::S32> mass_modify_list;

    // remove list
    PS::ReallocatableArray<PS::S32> remove_list;
    PS::ReallocatableArray<PS::S64> remove_id_record;
    PS::S32 n_remove;

    // search cluster
    SearchCluster search_cluster;

    // safety check flag
    bool read_parameters_flag;
    bool read_data_flag;
    bool initial_parameters_flag;
    bool initial_step_flag;

    // MPI 
    PS::S32 my_rank;
    PS::S32 n_proc;

    // FPDS flag
    static bool initial_fdps_flag;

    //! initialization
    PeTar(): 
        input_parameters(),
#ifdef BSE_BASE
        bse_parameters(),
        rand_parameters(),
#endif
#ifdef GALPY
        galpy_parameters(),
#endif
#ifdef EXTERNAL_HARD
        external_hard_parameters(),
#endif
#ifdef PROFILE
        // profile
        dn_loop(0), profile(), n_count(), n_count_sum(), tree_soft_profile(), fprofile(), 
#endif
        stat(), fstatus(), time_kick(0.0),
        escaper(), fesc(),
        file_header(), system_soft(), id_adr_map(),
        n_loop(0), domain_decompose_weight(1.0), dinfo(), pos_domain(NULL), 
        dt_manager(),
        tree_nb(), tree_soft(), 
#ifdef BSE_BASE
        rand_manager(),
#endif
#ifdef GALPY
        galpy_manager(),
#endif
        hard_manager(), system_hard_one_cluster(), system_hard_isolated(), 
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        system_hard_connected(), 
#endif
        n_interrupt_glb(0),
        mass_modify_list(), remove_list(), remove_id_record(),
        search_cluster(),
        read_parameters_flag(false), read_data_flag(false), initial_parameters_flag(false), initial_step_flag(false) {
        assert(initial_fdps_flag);
        my_rank = PS::Comm::getRank();
        n_proc = PS::Comm::getNumberOfProc();
     }


    //! regular block time step
    PS::F64 regularTimeStep(const PS::F64 _dt) {
        // regularize dt_tree
        PS::F64 dt = 1.0;
        if (_dt<1) while (dt>_dt) dt *= 0.5;
        else {
            while (dt<=_dt) dt *= 2.0;
            dt *= 0.5;
        }
        return dt;
    }

    //! tree for neighbor searching.
    void treeNeighborSearch() {
#ifdef PROFILE
        profile.tree_nb.start();
        tree_nb.clearNumberOfInteraction();
        tree_nb.clearTimeProfile();
#endif
#ifdef USE_SIMD
        tree_nb.calcForceAllAndWriteBack(SearchNeighborEpEpSimd(), system_soft, dinfo);
#elif USE_FUGAKU
        tree_nb.calcForceAllAndWriteBack(SearchNeighborEpEpFugaku(), system_soft, dinfo);
#else
        tree_nb.calcForceAllAndWriteBack(SearchNeighborEpEpNoSimd(), system_soft, dinfo);
#endif
        
#ifdef PROFILE
        tree_nb_profile += tree_nb.getTimeProfile();
        //profile.tree_nb.barrier();
        //PS::Comm::barrier();
        profile.tree_nb.end();
        profile.tree_nb.tbar = profile.tree_nb.time - tree_nb_profile.getTotalTime();
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
        search_cluster.writeAndSendBackPtcl(system_soft, system_hard_connected.getPtcl(), mass_modify_list);
        mass_modify_list.resizeNoInitialize(0);
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

#ifdef USE_GPU
        const PS::S32 n_walk_limit = 200;
        const PS::S32 tag_max = 1;
        PS::F64 eps2 = EPISoft::eps*EPISoft::eps;
        PS::F64 rout2 = EPISoft::r_out*EPISoft::r_out;
        PS::F64 G= ForceSoft::grav_const;
#ifdef PARTICLE_SIMULATOR_GPU_MULIT_WALK_INDEX
        tree_soft.calcForceAllAndWriteBackMultiWalkIndex(CalcForceWithLinearCutoffCUDAMultiWalk(my_rank, eps2, rout2, G),
                                                         RetrieveForceCUDA,
                                                         tag_max,
                                                         system_soft,
                                                         dinfo,
                                                         n_walk_limit);
#else // no multi-walk index
        tree_soft.calcForceAllAndWriteBackMultiWalk(CalcForceWithLinearCutoffCUDA(my_rank, eps2, rout2, G),
                                                    RetrieveForceCUDA,
                                                    tag_max,
                                                    system_soft,
                                                    dinfo,
                                                    n_walk_limit);
#endif // multi-walk index

#elif USE_FUGAKU
        PS::F64 eps2 = EPISoft::eps*EPISoft::eps;
        PS::F64 rout2 = EPISoft::r_out*EPISoft::r_out;
        PS::F64 G= ForceSoft::grav_const;
        tree_soft.calcForceAllAndWriteBack(CalcForceEpEpWithLinearCutoffFugaku(eps2, rout2, G),
#ifdef USE_QUAD
                                           CalcForceEpSpQuadFugaku(eps2, G),
#else // no quad
                                           CalcForceEpSpMonoFugaku(eps2, G),
#endif // end quad
                                           system_soft,
                                           dinfo);
        
#elif USE_SIMD // end use_gpu
        tree_soft.calcForceAllAndWriteBack(CalcForceEpEpWithLinearCutoffSimd(),
#ifdef USE_QUAD
                                           CalcForceEpSpQuadSimd(),
#else // no quad
                                           CalcForceEpSpMonoSimd(),
#endif // end quad
                                           system_soft,
                                           dinfo);
#else // end use_simd
        tree_soft.calcForceAllAndWriteBack(CalcForceEpEpWithLinearCutoffNoSimd(),
#ifdef USE_QUAD
                                           CalcForceEpSpQuadNoSimd(),
#else
                                           CalcForceEpSpMonoNoSimd(),
#endif
                                           system_soft,
                                           dinfo);
#endif // end else

#ifdef PROFILE
        n_count.ep_ep_interact     += tree_soft.getNumberOfInteractionEPEPLocal();
        n_count_sum.ep_ep_interact += tree_soft.getNumberOfInteractionEPEPGlobal();
        n_count.ep_sp_interact     += tree_soft.getNumberOfInteractionEPSPLocal();
        n_count_sum.ep_sp_interact += tree_soft.getNumberOfInteractionEPSPGlobal(); 

        tree_soft_profile += tree_soft.getTimeProfile();
        domain_decompose_weight = tree_soft_profile.calc_force;

        //profile.tree_soft.barrier();
        //PS::Comm::barrier();
        profile.tree_soft.end();
        profile.tree_soft.tbar = profile.tree_soft.time - tree_soft_profile.getTotalTime();
#endif
    }

    //! correct force due to change over function by using particle tree neighbor search
    void treeForceCorrectChangeoverTreeNeighbor() {
        // all particles
        SystemHard::correctForceWithCutoffTreeNeighborOMP<SystemSoft, FPSoft, TreeForce, EPJSoft>(system_soft, tree_soft, system_soft.getNumberOfParticleLocal(), hard_manager.ap_manager);        
    }

    //! correct force due to change over function
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
        SystemHard::correctForceWithCutoffTreeNeighborOMP<SystemSoft, FPSoft, TreeForce, EPJSoft>(system_soft, tree_soft, n_loc, hard_manager.ap_manager);

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

    //! calculate external force
    void externalForce() {
#ifdef PROFILE
        profile.other.start();
#endif

#ifdef GALPY
        // external force and potential
        // update the types and arguments 
        galpy_manager.updatePotential(stat.time, true);

        galpy_manager.resetPotAcc();
        galpy_manager.calcMovePotAccFromPot(stat.time, &stat.pcm.pos[0]);

        PS::S64 n_loc_all = system_soft.getNumberOfParticleLocal();
#pragma omp parallel for
        for (int i=0; i<n_loc_all; i++) {
            auto& pi = system_soft[i];
            double acc[3], pot;
#ifdef RECORD_CM_IN_HEADER
            PS::F64vec pos_correct=pi.pos + stat.pcm.pos;
            galpy_manager.calcAccPot(acc, pot, stat.time, input_parameters.gravitational_constant.value*pi.mass, &pos_correct[0], &pi.pos[0]);
#else
            PS::F64vec pos_center=pi.pos - stat.pcm.pos;
            galpy_manager.calcAccPot(acc, pot, stat.time, input_parameters.gravitational_constant.value*pi.mass, &pi.pos[0], &pos_center[0]);
#endif
            assert(!std::isinf(acc[0]));
            assert(!std::isnan(acc[0]));
            assert(!std::isinf(pot));
            assert(!std::isnan(pot));
            pi.acc[0] += acc[0]; 
            pi.acc[1] += acc[1]; 
            pi.acc[2] += acc[2]; 
            pi.pot_tot += pot;
            pi.pot_soft += pot;
#ifdef EXTERNAL_POT_IN_PTCL
            pi.pot_ext = pot;
#endif
        }
#endif //GALPY

#ifdef EXTERNAL_HARD
#ifndef GALPY
        // update time and gas density
        hard_manager.h4_manager.interaction.ext_force.updateTime(stat.time);
#endif
#endif
        
#ifdef PROFILE
        profile.other.barrier();
        PS::Comm::barrier();
        profile.other.end();
#endif
    }

    // correct c.m. vel due to external potential to avoid large rsearch 
    void correctPtclVelCM(const PS::F64& _dt) {
        PS::F64vec dv=PS::F64vec(0.0);
 
#ifdef GALPY
        PS::F64vec acc;
        PS::F64 pot;
        // evaluate center of mass acceleration
        PS::F64vec pos_zero=PS::F64vec(0.0);
        // set zero mass to avoid duplicate anti force to potential set
        galpy_manager.calcAccPot(&acc[0], pot, stat.time, 0, &stat.pcm.pos[0], &pos_zero[0]);
        dv = acc*_dt;
#endif        

        stat.pcm.vel += dv;
        for (int i=0; i<stat.n_all_loc; i++) system_soft[i].vel -= dv;

        //auto& adr = search_cluster.getAdrSysOneCluster();
        //for (int i=0; i<adr.size(); i++) system_soft[adr[i]].vel -= dv;

        auto& ptcl_iso=system_hard_isolated.getPtcl();
        for (int i=0; i<ptcl_iso.size(); i++) ptcl_iso[i].vel -= dv;

        //const PS::S64 n_tot= system_soft.getNumberOfParticleLocal();
        //const PS::S32 n_artifical_per_group = hard_manager.ap_manager.getArtificialParticleN();
        //for(PS::S32 i=stat.n_real_loc; i<n_tot; i+= n_artifical_per_group) {
        //    auto* pcm = hard_manager.ap_manager.getCMParticles(&(system_soft[i]));
        //    pcm->vel -= dv;
        //}        

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        auto& ptcl_con=system_hard_connected.getPtcl();
        for (int i=0; i<ptcl_con.size(); i++) ptcl_con[i].vel -= dv;
#endif

      //Ptcl::vel_cm += dv;
//#ifdef PETAR_DEBUG
//        if (my_rank==0) std::cerr<<"Ptcl::vel_cm: "<<Ptcl::vel_cm<<std::endl;
//#endif
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
//#ifdef STELLAR_EVOLUTION
                // shift time interrupt in order to get consistent time for stellar evolution in the next drift
                //system[i].time_record    -= dt;
                //system[i].time_interrupt -= dt;
//#endif
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

#ifdef GALPY
        galpy_manager.kickMovePot(_dt_kick);
#endif

#ifdef RECORD_CM_IN_HEADER
        // correct Ptcl:vel_cm
        correctPtclVelCM(_dt_kick);
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
    /*! 
     */
    void drift(const PS::F64 _dt_drift) {
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
        system_hard_one_cluster.writeBackPtclForOneClusterOMP(system_soft, mass_modify_list);
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
#ifdef HARD_CHECK_ENERGY
        // reset slowdown energy correction
        system_hard_isolated.energy.resetEnergyCorrection();
#endif
        // integrate multi cluster A
        system_hard_isolated.driveForMultiClusterOMP(_dt_drift, &(system_soft[0]));
        //system_hard_isolated.writeBackPtclForMultiCluster(system_soft, search_cluster.adr_sys_multi_cluster_isolated_,remove_list);
        system_hard_isolated.writeBackPtclForMultiCluster(system_soft, mass_modify_list);
        // integrate multi cluster A

#ifdef PROFILE
        profile.hard_isolated.barrier();
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
        
        search_cluster.writeAndSendBackPtcl(system_soft, system_hard_connected.getPtcl(), mass_modify_list);
        system_hard_connected.updateTimeWriteBack();

        // integrate multi cluster B
#ifdef PROFILE
        profile.hard_connected.barrier();
        PS::Comm::barrier();
        profile.hard_connected.end();
#endif
#endif
        // drift cm
        stat.pcm.pos += stat.pcm.vel*_dt_drift;

#ifdef GALPY
        galpy_manager.driftMovePot(_dt_drift);
#endif
        
        Ptcl::group_data_mode = GroupDataMode::cm;
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
        if (abs(time_kick-stat.time)>1e-13) {
            std::cerr<<"Error: kick time ("<<time_kick<<") and system time ("<<stat.time<<") are inconsistent!"<<std::endl;
            abort();
        }
        time_kick = stat.time; // escape the problem of round-off error
        if (stat.n_real_glb>1) {
            assert(stat.time == system_hard_one_cluster.getTimeOrigin());
            assert(stat.time == system_hard_isolated.getTimeOrigin());
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            assert(stat.time == system_hard_connected.getTimeOrigin());
#endif
        }
        return true;
    }
    
    //! write back hard particles to global system
    void writeBackHardParticles() {
#ifdef PROFILE
        profile.output.start();
#endif
        system_hard_isolated.writeBackPtclForMultiCluster(system_soft, mass_modify_list);
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        // update gloabl particle system and send receive remote particles
        search_cluster.writeAndSendBackPtcl(system_soft, system_hard_connected.getPtcl(), mass_modify_list);
        system_hard_connected.updateTimeWriteBack();
#endif        
        mass_modify_list.resizeNoInitialize(0);

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
        search_cluster.writeAndSendBackPtcl(system_soft, system_hard_connected.getPtcl(), mass_modify_list);
        system_hard_connected.updateTimeWriteBack();
        mass_modify_list.resizeNoInitialize(0);
#endif
#endif    
#ifdef PROFILE
        profile.search_cluster.barrier();
        PS::Comm::barrier();
        profile.search_cluster.end();
#endif
    }
    
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
    /*!
      @param[in] _enforce: do domain decompose without check n_loop (false)
     */
    void domainDecompose(const bool _enforce=false) {
#ifdef PROFILE
        // > 6. Domain decomposition
        profile.domain.start();
#endif
        // Domain decomposition, parrticle exchange and force calculation
        if(n_loop % 16 == 0 || _enforce) {
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
//        stat.shiftToOriginFrame(&system_soft[0], stat.n_real_loc);
        if (_initial_flag) {
            // calculate initial energy
            stat.energy.clear();
#ifdef RECORD_CM_IN_HEADER
            stat.energy.calc(&system_soft[0], stat.n_real_loc, true, &(stat.pcm.pos), &(stat.pcm.vel));
#else
            stat.energy.calc(&system_soft[0], stat.n_real_loc, true);
#endif
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            stat.energy.getSumMultiNodes(true);
#endif
#ifdef HARD_CHECK_ENERGY
            stat.energy.ekin_sd = stat.energy.ekin;
            stat.energy.epot_sd = stat.energy.epot;
#endif
#ifndef RECORD_CM_IN_HEADER
            stat.calcCenterOfMass(&system_soft[0], stat.n_real_loc);
#endif
            //Ptcl::vel_cm = stat.pcm.vel;

        }
        else {
#ifdef HARD_CHECK_ENERGY
            // reset sd energy reference to no slowdown case
            PS::F64 energy_old = stat.energy.ekin + stat.energy.epot;
            stat.energy.etot_sd_ref -= stat.energy.ekin_sd + stat.energy.epot_sd - energy_old;
#endif
#ifdef RECORD_CM_IN_HEADER
            stat.energy.calc(&system_soft[0], stat.n_real_loc,false, &(stat.pcm.pos), &(stat.pcm.vel));
#else
            stat.energy.calc(&system_soft[0], stat.n_real_loc);
#endif
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            stat.energy.getSumMultiNodes();
#endif

#ifdef HARD_CHECK_ENERGY
            HardEnergy energy_local = system_hard_one_cluster.energy;
            energy_local += system_hard_isolated.energy;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            energy_local += system_hard_connected.energy;
#endif
            // hard energy error
            stat.energy.error_hard_cum += PS::Comm::getSum(energy_local.de);
            stat.energy.error_hard_sd_cum += PS::Comm::getSum(energy_local.de_sd);
            // energy correction due to slowdown 
            PS::F64 ekin_sd_correction = PS::Comm::getSum(energy_local.ekin_sd_correction);
            stat.energy.ekin_sd = stat.energy.ekin + ekin_sd_correction;
            PS::F64 epot_sd_correction = PS::Comm::getSum(energy_local.epot_sd_correction);
            stat.energy.epot_sd = stat.energy.epot + epot_sd_correction;
            PS::F64 etot_sd_correction = ekin_sd_correction + epot_sd_correction;
            PS::F64 de_change_cum = PS::Comm::getSum(energy_local.de_change_cum);
            PS::F64 de_change_binary_interrupt = PS::Comm::getSum(energy_local.de_change_binary_interrupt);
            PS::F64 de_change_modify_single = PS::Comm::getSum(energy_local.de_change_modify_single);
            stat.energy.etot_ref += de_change_cum;
            stat.energy.de_change_cum += de_change_cum;
            stat.energy.de_change_binary_interrupt += de_change_binary_interrupt;
            stat.energy.de_change_modify_single += de_change_modify_single;
            // for total energy reference, first add the cumulative change due to slowdown change in the integration (referring to no slowdown case), then add slowdown energy correction from current time
            PS::F64 de_sd_change_cum  = PS::Comm::getSum(energy_local.de_sd_change_cum) + etot_sd_correction;
            PS::F64 de_sd_change_binary_interrupt = PS::Comm::getSum(energy_local.de_sd_change_binary_interrupt);
            PS::F64 de_sd_change_modify_single = PS::Comm::getSum(energy_local.de_sd_change_modify_single);
            stat.energy.etot_sd_ref += de_sd_change_cum; 
            stat.energy.de_sd_change_cum += de_sd_change_cum;
            stat.energy.de_sd_change_binary_interrupt += de_sd_change_binary_interrupt;
            stat.energy.de_sd_change_modify_single += de_sd_change_modify_single;

            system_hard_one_cluster.energy.clear();
            system_hard_isolated.energy.clear();
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            system_hard_connected.energy.clear();
#endif
#endif
#ifndef RECORD_CM_IN_HEADER
            stat.calcCenterOfMass(&system_soft[0], stat.n_real_loc);
#endif
            //Ptcl::vel_cm = stat.pcm.vel;
        }
        //stat.shiftToCenterOfMassFrame(&system_soft[0], stat.n_real_loc);

#ifdef PROFILE
        profile.status.barrier();
        PS::Comm::barrier();
        profile.status.end();
#endif
    }

    //! output data
    void output() {
#ifdef PROFILE
        profile.output.start();
#endif
        bool print_flag = input_parameters.print_flag;
        int write_style = input_parameters.write_style.value;
        std::cout<<std::setprecision(PRINT_PRECISION);

        // print status
        if(print_flag) {
            std::cout<<std::endl;
            stat.print(std::cout);
        }
#ifdef GALPY
        if (print_flag) galpy_manager.printData(std::cout);
#endif
        // write status, output to separate snapshots
        if(write_style==1) {

            // data output
            file_header.n_body = stat.n_real_glb;
            file_header.time = stat.time;
            file_header.nfile++;
#ifdef RECORD_CM_IN_HEADER
            file_header.pos_offset = stat.pcm.pos;
            file_header.vel_offset = stat.pcm.vel;
#endif

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

            if(my_rank==0) {
                // status output
                stat.printColumn(fstatus, WRITE_WIDTH);
                fstatus<<std::endl;

#ifdef GALPY
                // for External potential
                galpy_manager.writePotentialPars(fname+".galpy", stat.time);
#endif
            }
#ifdef BSE_BASE
            std::string fname_seed = fname+".randseeds";
            rand_manager.writeRandSeeds(fname_seed.c_str());
#endif
        }
        // write all information in to fstatus
        else if(write_style==2&&my_rank==0) {
            // write snapshot with one line
            stat.printColumn(fstatus, WRITE_WIDTH);
            for (int i=0; i<stat.n_real_loc; i++) {
#ifdef RECORD_CM_IN_HEADER
                system_soft[i].printColumnWithOffset(stat.pcm, fstatus, WRITE_WIDTH);
#else
                system_soft[i].printColumn(fstatus, WRITE_WIDTH);
#endif
            }
            fstatus<<std::endl;
        }
        // write status only
        else if(write_style==3&&my_rank==0) {
            stat.printColumn(fstatus, WRITE_WIDTH);
            fstatus<<std::endl;
        }

        // save current error
        stat.energy.saveEnergyError();


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
        
        n_count.clearClusterCount();
        n_count.clusterCount(1, n_hard_single);

        const PS::S32  n_isolated_cluster = system_hard_isolated.getNumberOfClusters();
        n_count.cluster_isolated += n_isolated_cluster;
        n_count_sum.cluster_isolated += PS::Comm::getSum(n_isolated_cluster);

        const PS::S32* isolated_cluster_n_list = system_hard_isolated.getClusterNumberOfMemberList();
        for (PS::S32 i=0; i<n_isolated_cluster; i++) n_count.clusterCount(isolated_cluster_n_list[i]);

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        const PS::S32  n_connected_cluster = system_hard_connected.getNumberOfClusters();
        n_count.cluster_connected += n_connected_cluster;
        n_count_sum.cluster_connected += PS::Comm::getSum(n_connected_cluster);
        const PS::S32* connected_cluster_n_list = system_hard_connected.getClusterNumberOfMemberList();
        for (PS::S32 i=0; i<n_connected_cluster; i++) n_count.clusterCount(connected_cluster_n_list[i]);
#endif

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        PS::S32 cluster_hist_size_loc = n_count.n_cluster.size();
        PS::S32 n_member_in_cluster_loc[cluster_hist_size_loc];
        PS::S32 n_cluster_loc[cluster_hist_size_loc];
        n_count.getherClusterCount(n_member_in_cluster_loc, n_cluster_loc, cluster_hist_size_loc);
        PS::S32 cluster_hist_size_recv[n_proc];
        PS::S32 cluster_hist_size_disp[n_proc+1];
        PS::Comm::allGather(&cluster_hist_size_loc, 1, cluster_hist_size_recv);
        cluster_hist_size_disp[0] = 0;
        for (PS::S32 i=0; i<n_proc; i++){
            cluster_hist_size_disp[i+1] = cluster_hist_size_recv[i] + cluster_hist_size_disp[i];
        }
        PS::S32 cluster_hist_size_total = cluster_hist_size_disp[n_proc];
        PS::S32 n_member_in_cluster_recv[cluster_hist_size_total];
        PS::S32 n_cluster_recv[cluster_hist_size_total];
        PS::Comm::allGatherV(n_member_in_cluster_loc, cluster_hist_size_loc, 
                             n_member_in_cluster_recv, cluster_hist_size_recv, cluster_hist_size_disp);
        PS::Comm::allGatherV(n_cluster_loc, cluster_hist_size_loc, 
                             n_cluster_recv, cluster_hist_size_recv, cluster_hist_size_disp);
        for (PS::S32 i=0; i<cluster_hist_size_total; i++) {
            n_count_sum.clusterCount(n_member_in_cluster_recv[i], n_cluster_recv[i]);
        }
#else
        n_count_sum.addClusterCount(n_count);
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
            profile.dumpName(std::cout);
            std::cout<<std::endl;

            //std::cout<<std::setw(PRINT_WIDTH)<<my_rank;
            profile_min.dump(std::cout,dn_loop);
            std::cout<<std::endl;
            profile.dump(std::cout,dn_loop);
            std::cout<<std::endl;

            std::cout<<"**** FDPS tree soft force time profile (local):\n";
            tree_soft_profile.dumpName(std::cout);
            std::cout<<std::endl;
            tree_soft_profile.dump(std::cout,dn_loop);
            std::cout<<std::endl;

            std::cout<<"**** Tree neighbor time profile (local):\n";
            tree_nb_profile.dumpName(std::cout);
            std::cout<<std::endl;
            tree_nb_profile.dump(std::cout,dn_loop);
            std::cout<<std::endl;

#if defined(USE_GPU) && defined(GPU_PROFILE)
            std::cout<<"**** GPU time profile (local):\n";
            gpu_profile.dumpName(std::cout);
            gpu_counter.dumpName(std::cout);
            std::cout<<std::endl;
            gpu_profile.dump(std::cout,dn_loop);
            gpu_counter.dump(std::cout,dn_loop);
            std::cout<<std::endl;
#endif

            std::cout<<"**** Number per step (global):\n";
            n_count_sum.dumpName(std::cout);
            std::cout<<std::endl;
            n_count_sum.dump(std::cout,dn_loop);
            std::cout<<std::endl;
                
            std::cout<<"**** Number of members in clusters (global):\n";
            n_count_sum.printHist(std::cout,dn_loop);
        }

        if(input_parameters.write_style.value>0) {
            fprofile<<std::setprecision(WRITE_PRECISION);
            fprofile<<std::setw(WRITE_WIDTH)<<my_rank;
            fprofile<<std::setw(WRITE_WIDTH)<<stat.time
                    <<std::setw(WRITE_WIDTH)<<dn_loop
                    <<std::setw(WRITE_WIDTH)<<stat.n_real_loc;
            profile.dump(fprofile, dn_loop, WRITE_WIDTH);
            profile.dumpBarrier(fprofile, dn_loop, WRITE_WIDTH);
            tree_soft_profile.dump(fprofile, dn_loop, WRITE_WIDTH);
            tree_nb_profile.dump(fprofile, dn_loop, WRITE_WIDTH);
#if defined(USE_GPU) && defined(GPU_PROFILE)
            gpu_profile.dump(fprofile, dn_loop, WRITE_WIDTH);
            gpu_counter.dump(fprofile, dn_loop, WRITE_WIDTH);
#endif
            n_count.dump(fprofile, dn_loop, WRITE_WIDTH);
            fprofile<<std::endl;
        }
    }
#endif

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

#ifdef STELLAR_EVOLUTION
    //! Correct potential energy due to modificaiton of particle mass
    void correctSoftPotMassChange() {
        // correct soft potential energy due to mass change
#pragma omp parallel for
        for (int k=0; k<mass_modify_list.size(); k++)  {
            PS::S32 i = mass_modify_list[k];
            auto& pi = system_soft[i];
            PS::F64 dpot = pi.dm*pi.pot_soft;
            stat.energy.etot_ref += dpot;
            stat.energy.de_change_cum += dpot;
            stat.energy.etot_sd_ref += dpot;
            stat.energy.de_sd_change_cum += dpot;
            pi.dm = 0.0;
            // ghost particle case, check in remove_particle instead
            //if(pi.mass==0.0&&pi.group_data.artificial.isUnused()) {
            //    remove_list.push_back(i);
            //}
        }
        mass_modify_list.resizeNoInitialize(0);
    }
#endif

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

        // first remove artificial particles
        system_soft.setNumberOfParticleLocal(stat.n_real_loc);

        // set flag for particles in remove_list
        PS::S32 n_remove = remove_list.size();
        for (PS::S32 i=0; i<n_remove; i++) {
            auto& pi = system_soft[remove_list[i]];
            pi.group_data.artificial.setParticleTypeToUnused(); // sign for removing
            // if mass is not zero, correct energy
            if (pi.mass>0) {
                PS::F64 dpot = pi.mass*pi.pot_tot;
                PS::F64 dkin = 0.5*pi.mass*(pi.vel*pi.vel);
                PS::F64 eloss = dpot + dkin;
                stat.energy.etot_ref -= eloss;
#ifdef HARD_CHECK_ENERGY
                stat.energy.de_change_cum -= eloss;
                stat.energy.etot_sd_ref -= eloss;
                stat.energy.de_sd_change_cum -= eloss;
#endif
                pi.mass = 0.0;
            }
        }
        remove_list.resizeNoInitialize(0);

        // check all particles escaper status and remove
        const PS::S32 num_thread = PS::Comm::getNumberOfThread();
        PS::ReallocatableArray<PS::S32> remove_list_thx[num_thread];
        for (PS::S32 i=0; i<num_thread; i++) remove_list_thx[i].resizeNoInitialize(0);

#pragma omp parallel
        { 
            const PS::S32 ith = PS::Comm::getThreadNum();
#pragma omp for 
            for (PS::S32 i=0; i<stat.n_real_loc; i++) {
                auto& pi = system_soft[i];
                if (escaper.isEscaper(pi,stat.pcm)) {
                    remove_list_thx[ith].push_back(i);
                    PS::F64 dpot = pi.mass*pi.pot_tot;
                    PS::F64 dkin = 0.5*pi.mass*(pi.vel*pi.vel);
                    PS::F64 eloss = dpot + dkin;
                    stat.energy.etot_ref -= eloss;
#ifdef HARD_CHECK_ENERGY
                    stat.energy.de_change_cum -= eloss;
                    stat.energy.etot_sd_ref -= eloss;
                    stat.energy.de_sd_change_cum -= eloss;
#endif                    
                }
                // Registered removed particles have already done energy correction
                else if (pi.mass==0.0&&pi.group_data.artificial.isUnused()) 
                    remove_list_thx[ith].push_back(i);
            }
        }

        // gether remove indices
        int n_esc=0;
        for (PS::S32 i=0; i<num_thread; i++) 
            for (PS::S32 k=0; k<remove_list_thx[i].size(); k++) {
                PS::S32 index=remove_list_thx[i][k];
                remove_list.push_back(index);
                remove_id_record.push_back(system_soft[index].id);
                if (system_soft[index].mass>0) {
                    if (input_parameters.write_style.value>0) {
                        fesc<<std::setw(WRITE_WIDTH)<<stat.time;
                        system_soft[index].printColumn(fesc,WRITE_WIDTH);
                        fesc<<std::endl;
                    }
                    n_esc++;
                }
            }

        // Remove particles
        n_remove = remove_list.size();
        system_soft.removeParticle(remove_list.getPointer(), remove_list.size());

        stat.n_escape_glb += PS::Comm::getSum(n_esc);
        stat.n_remove_glb += PS::Comm::getSum(n_remove);

        // reset particle number
        stat.n_real_loc = stat.n_real_loc-remove_list.size();
        system_soft.setNumberOfParticleLocal(stat.n_real_loc);
        stat.n_real_glb = system_soft.getNumberOfParticleGlobal();
        remove_list.resizeNoInitialize(0);

#ifdef PETAR_DEBUG
#pragma omp parallel for
        for(PS::S32 i=0; i<stat.n_real_loc; i++){
            assert(system_soft[i].mass!=0.0);
            assert(!std::isinf(system_soft[i].pos[0]));
            assert(!std::isnan(system_soft[i].pos[0]));
            assert(!std::isinf(system_soft[i].vel[0]));
            assert(!std::isnan(system_soft[i].vel[0]));
        }
#endif

#ifdef PROFILE
        profile.other.barrier();
        PS::Comm::barrier();
        profile.other.end();
#endif
    }


    //! get number of removed particles in local
    PS::S32 getNumberOfRemovedParticlesLocal() const {
        return remove_id_record.size();
    }

    //! get number of removed particles in global
    PS::S32 getNumberOfRemovedParticlesGlobal() const {
        PS::S32 n_remove_loc = remove_id_record.size();
        return PS::Comm::getSum(n_remove_loc);
    }

    //! get removed particle ID list local
    PS::S64* getRemovedIDListLocal() const {
        return remove_id_record.getPointer();
    }

    //! clear removed particle ID list
    void clearRemovedIDList() {
        remove_id_record.resizeNoInitialize(0);
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
    static void initialFDPS(int argc, char *argv[]) {
        if (!initial_fdps_flag) {
            PS::Initialize(argc, argv);
            initial_fdps_flag = true;
        }
    }

    //! finalize FDPS (close MPI) 
    static void finalizeFDPS() {
        if (initial_fdps_flag) {
            PS::Finalize();
            initial_fdps_flag = false;
        }
    }

#ifdef FDPS_COMM
    // create a MPI communicator
    /*! create a MPI communicator with selected rank list
      Notice this cannot be done during the integration, only once at the initialization
      @param[in] n: number of ranks
      @param[in] rank: rank list
     */
    void createCommunicator(const int n, const int rank[]) {
        assert(initial_fdps_flag);
        assert(!initial_parameters_flag);
        comm_info.create(n,rank);
        system_soft.setCommInfo(comm_info);
        dinfo.setCommInfo(comm_info);
        tree_nb.setCommInfo(comm_info);
        tree_soft.setCommInfo(comm_info);
        my_rank = comm_info.getRank();
        n_proc = comm_info.getNumberOfProc();
    }

    // set a MPI communicator
    /*! set a MPI communicator from input
      Notice this cannot be done during the integration, only once at the initialization
      @param[in] _comm_info: FDPS CommInfo class
     */
    void setCommunicator(const PS::CommInfo &_comm_info) {
        assert(initial_fdps_flag);
        assert(!initial_parameters_flag);
        comm_info = _comm_info;
        system_soft.setCommInfo(comm_info);
        dinfo.setCommInfo(comm_info);
        tree_nb.setCommInfo(comm_info);
        tree_soft.setCommInfo(comm_info);
        my_rank = comm_info.getRank();
        n_proc = comm_info.getNumberOfProc();
    }
#endif

    //! print terminal Logo
    void printLogo(std::ostream & fout) const {
        fout<<"\n ******************************************\n"
            <<"  ██████╗ ███████╗████████╗ █████╗ ██████╗\n"
            <<"  ██╔══██╗██╔════╝╚══██╔══╝██╔══██╗██╔══██╗\n"
            <<"  ██████╔╝█████╗     ██║   ███████║██████╔╝\n"
            <<"  ██╔═══╝ ██╔══╝     ██║   ██╔══██║██╔══██╗\n"
            <<"  ██║     ███████╗   ██║   ██║  ██║██║  ██║\n"
            <<"  ╚═╝     ╚══════╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝\n"
            <<" ******************************************\n"
            <<"Version: "<<GetVersion()
            <<std::endl;
    }

    //! print reference to cite
    void printReference(std::ostream & fout) const{
        fout<<"============== PeTar ================\n"
            <<" Online document: (preparing)\n"
            <<" GitHub page: https://github.com/lwang-astro/PeTar\n"
            <<" License: MIT\n"
            <<" Please cite the following papers when you use PeTar for any publications\n"
            <<"    FDPS: see above FPDS Logo message\n";
        AR::printReference(fout);
        fout<<"    PeTar: Wang L., Iwasawa M., Nitadori K., Makino J., 2020, MNRAS, 497, 536\n";
#ifdef BSE_BASE
        BSEManager::printReference(fout);
#endif
#ifdef GALPY
        GalpyManager::printReference(fout);
#endif
        fout<<" Copyright (C) 2017\n"
            <<"    Long Wang, Masaki Iwasawa, Keigo Nitadori, Junichiro Makino and many others\n";
        fout<<"====================================="
            <<std::endl;
    }

    //! print features
    void printFeatures(std::ostream & fout) const {
#ifdef USE_QUAD
        fout<<"Use quadrupole moment in tree force calculation\n";
#endif        

#ifdef SOFT_PERT
#ifdef TIDAL_TENSOR_3RD
        fout<<"Use 3rd order tidal tensor method\n";
#else
        fout<<"Use 2nd order tidal tensor method\n";
#endif
#endif

#ifdef ORBIT_SAMPLING
        fout<<"Use orbit-sampling method\n";
#else
        fout<<"Use Pseudoparticle multipole method\n";
#endif

#ifdef STELLAR_EVOLUTION
        fout<<"Use stellar evolution method: ";
#ifdef BSE_BASE
        fout<<BSEManager::getBSEName()<<std::endl;
#else
        fout<<"Base\n";
#endif
#endif

#ifdef GALPY
        fout<<"Use external potential: Galpy\n";
#endif 

#ifdef KDKDK_2ND
        fout<<"Use 2nd-order KDKDK mode for tree step\n";
#endif
#ifdef KDKDK_4TH
        fout<<"Use 4th-order KDKDK mode for tree step\n";
#endif

#ifdef CLUSTER_VELOCITY
        fout<<"Use orbit-dependent neighbor criterion\n";
#endif

#ifdef HARD_CHECK_ENERGY
        fout<<"Check energy of short-range (hard) integration\n";
#endif

#ifdef HARD_COUNT_NO_NEIGHBOR
        fout<<"Count isolated particles in hard clusters\n";
#endif

#ifdef USE_SIMD
        fout<<"Use SIMD\n";
#ifdef P3T_64BIT
        fout<<"Use 64 bit SIMD n";
#endif
#endif

#ifdef USE_FUGAKU
        fout<<"Use Fugaku\n";
#endif

#ifdef USE_GPU
        fout<<"Use GPU\n";
#endif

#ifdef GPERF_PROFILE
        fout<<"Use gperftools\n";
#endif

#ifdef PROFILE
        fout<<"Calculate profile\n";
#endif

#ifdef GPU_PROFILE
        fout<<"Calculate GPU profile\n";
#endif

        AR::printFeatures(fout);
        H4::printFeatures(fout);
    }

    //! print features
    void printDebugFeatures(std::ostream & fout) const {
#ifdef PETAR_DEBUG
        fout<<"Debug mode: PETAR\n";
#endif
#ifdef STABLE_CHECK_DEBUG
        fout<<"Debug mode: STABLE_CHECK\n";
#endif
#ifdef CLUSTER_DEBUG
        fout<<"Debug mode: CLUSTER\n";
#endif
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        fout<<"Debug mode: ARTIFICIAL_PARTICLE\n";
#endif
#ifdef HARD_DEBUG
        fout<<"Debug mode: HARD\n";
#endif
#ifdef NAN_CHECK_DEBUG
        fout<<"Debug mode: NAN_CHECK\n";
#endif
        AR::printDebugFeatures(fout);
        H4::printDebugFeatures(fout);
    }

    //! reading input parameters using getopt method
    /*! 
      @param[in] argc: number of options
      @param[in] argv: string of options
      \return -1 if help is used
     */
    int readParameters(int argc, char *argv[]) {
        // print reference
        if (my_rank==0) {
            printLogo(std::cerr);
            printReference(std::cerr);
            printFeatures(std::cerr);
            printDebugFeatures(std::cerr);
            std::cerr<<"====================================="<<std::endl;
        }

        //assert(initial_fdps_flag);
        assert(!read_parameters_flag);
        // reading parameters
        opterr = 0;
        read_parameters_flag = true;
        if (my_rank==0) input_parameters.print_flag=true;
        else input_parameters.print_flag=false;
        int read_flag = input_parameters.read(argc,argv);
#ifdef BSE_BASE
        if (my_rank==0) bse_parameters.print_flag=true;
        else bse_parameters.print_flag=false;
        bse_parameters.read(argc,argv);
        if (my_rank==0) rand_parameters.print_flag=true;
        else rand_parameters.print_flag=false;
        rand_parameters.read(argc,argv);
#endif
#ifdef GALPY
        if (my_rank==0) galpy_parameters.print_flag=true;
        else galpy_parameters.print_flag=false;
        galpy_parameters.read(argc,argv);
#endif
#ifdef EXTERNAL_HARD
        if (my_rank==0) external_hard_parameters.print_flag=true;
        else external_hard_parameters.print_flag=false;
        external_hard_parameters.read(argc,argv);
#endif

        // help case, return directly
        if (read_flag==-1) {
            // avoid segmentation fault due to FDPS clear function bug
            tree_nb.initialize(input_parameters.n_glb.value, input_parameters.theta.value, input_parameters.n_leaf_limit.value, input_parameters.n_group_limit.value);
            tree_soft.initialize(input_parameters.n_glb.value, input_parameters.theta.value, input_parameters.n_leaf_limit.value, input_parameters.n_group_limit.value);

            return read_flag;
        }

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
            std::string fproname=input_parameters.fname_snp.value+".prof.rank."+std::to_string(my_rank);
            if(input_parameters.append_switcher.value==1) fprofile.open(fproname.c_str(),std::ofstream::out|std::ofstream::app);
            else  {
                fprofile.open(fproname.c_str(),std::ofstream::out);

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
            fprofile<<std::setprecision(WRITE_PRECISION);
        }
#endif    

        // open output files
        // status information output
        std::string& fname_snp = input_parameters.fname_snp.value;
        if(write_style>0&&my_rank==0) {
            if(input_parameters.append_switcher.value==1) 
                fstatus.open((fname_snp+".status").c_str(),std::ofstream::out|std::ofstream::app);
            else {
                fstatus.open((fname_snp+".status").c_str(),std::ofstream::out);
                // write titles of columns
                stat.printColumnTitle(fstatus,WRITE_WIDTH);
                if (write_style==2) {
                    for (int i=0; i<stat.n_real_loc; i++) system_soft[0].printColumnTitle(fstatus, WRITE_WIDTH);
                }
                fstatus<<std::endl;
            }
            fstatus<<std::setprecision(WRITE_PRECISION);
        }

        if(write_style>0) {
            // open escaper file
            std::string my_rank_str = std::to_string(my_rank);
            std::string fname_esc = fname_snp + ".esc." + my_rank_str;
            if(input_parameters.append_switcher.value==1) 
                fesc.open(fname_esc.c_str(), std::ofstream::out|std::ofstream::app);
            else {
                fesc.open(fname_esc.c_str(), std::ofstream::out);
                // write titles of columns, will cause issue when gether different MPI ranks
                // fesc<<std::setw(WRITE_WIDTH)<<"Time";
                // FPSoft::printColumnTitle(fesc,WRITE_WIDTH);
                // fesc<<std::endl;
            }
            fesc<<std::setprecision(WRITE_PRECISION);

#ifdef STELLAR_EVOLUTION
#ifdef BSE_BASE
            // open SSE/BSE file
            std::string fsse_name = fname_snp + fsse_par_suffix + "." + my_rank_str;
            std::string fbse_name = fname_snp + fbse_par_suffix + "." + my_rank_str;
            if(input_parameters.append_switcher.value==1) {
                hard_manager.ar_manager.interaction.fout_sse.open(fsse_name.c_str(), std::ofstream::out|std::ofstream::app);
                hard_manager.ar_manager.interaction.fout_bse.open(fbse_name.c_str(), std::ofstream::out|std::ofstream::app);
            }
            else {
                hard_manager.ar_manager.interaction.fout_sse.open(fsse_name.c_str(), std::ofstream::out);
                hard_manager.ar_manager.interaction.fout_bse.open(fbse_name.c_str(), std::ofstream::out);
            }
            hard_manager.ar_manager.interaction.fout_sse<<std::setprecision(WRITE_PRECISION);
            hard_manager.ar_manager.interaction.fout_bse<<std::setprecision(WRITE_PRECISION);
#else
            // open interrupt file
            std::string finterrupt_name = fname_snp + ".interrupt." + my_rank_str;
            if(input_parameters.append_switcher.value==1) 
                hard_manager.ar_manager.interaction.fout_interrupt.open(finterrupt_name.c_str(), std::ofstream::out|std::ofstream::app);
            else 
                hard_manager.ar_manager.interaction.fout_interrupt.open(finterrupt_name.c_str(), std::ofstream::out);
            hard_manager.ar_manager.interaction.fout_interrupt<<std::setprecision(WRITE_PRECISION);
#endif 
#endif

#ifdef ADJUST_GROUP_PRINT
            // open file for new/end group information
            if (input_parameters.adjust_group_write_option.value==1) {
                std::string fgroup_name = fname_snp + ".group." + my_rank_str;
                if(input_parameters.append_switcher.value==1) 
                    hard_manager.h4_manager.fgroup.open(fgroup_name.c_str(), std::ofstream::out|std::ofstream::app);
                else 
                    hard_manager.h4_manager.fgroup.open(fgroup_name.c_str(), std::ofstream::out);
                hard_manager.h4_manager.fgroup<<std::setprecision(WRITE_PRECISION);
            }
#endif
        }

#ifdef HARD_DUMP
        // initial hard_dump 
        const PS::S32 num_thread = PS::Comm::getNumberOfThread();
        hard_dump.initial(num_thread, my_rank);
#endif

        // particle system
        system_soft.initialize();

        // domain decomposition
        system_soft.setAverageTargetNumberOfSampleParticlePerProcess(input_parameters.n_smp_ave.value);
        const PS::F32 coef_ema = 0.2;
        domain_decompose_weight=1.0;
        dinfo.initialize(coef_ema);

        if(pos_domain==NULL) {
            pos_domain = new PS::F64ort[n_proc];
            for(PS::S32 i=0; i<n_proc; i++) pos_domain[i] = dinfo.getPosDomain(i);
        }

        // tree for neighbor search
        tree_nb.initialize(input_parameters.n_glb.value, input_parameters.theta.value, input_parameters.n_leaf_limit.value, input_parameters.n_group_limit.value);

        // tree for force
        PS::S64 n_tree_init = input_parameters.n_glb.value + input_parameters.n_bin.value;
        tree_soft.initialize(n_tree_init, input_parameters.theta.value, input_parameters.n_leaf_limit.value, input_parameters.n_group_limit.value);

        // initial search cluster
        search_cluster.initialize();

        return read_flag;
    }
    
    //! reading data set from file, filename is given by readParameters
    void readDataFromFile() {
        assert(read_parameters_flag);
        if(input_parameters.fname_inp.value=="__NONE__") {
            std::cerr<<"Error: the input data filename is not provided, make sure your options have the correct format."<<std::endl;
            abort();
        }

        PS::S32 data_format = input_parameters.data_format.value;
        auto* data_filename = input_parameters.fname_inp.value.c_str();
                
        if(data_format==1||data_format==2||data_format==4)
            system_soft.readParticleAscii(data_filename, file_header);
        else
            system_soft.readParticleBinary(data_filename, file_header);
        PS::Comm::broadcast(&file_header, 1, 0);
        PS::S64 n_glb = system_soft.getNumberOfParticleGlobal();
        PS::S64 n_loc = system_soft.getNumberOfParticleLocal();

        if(input_parameters.print_flag) {
            std::cout<<std::setprecision(WRITE_PRECISION);
            std::cout<<"----- Reading file: "<<data_filename<<" -----"<<std::endl
                     <<"Number of particles = "<<n_glb<<std::endl
                     <<"Time = "<<file_header.time<<std::endl;
        }

        stat.time = file_header.time;
        stat.n_real_glb = stat.n_all_glb = n_glb;
        stat.n_real_loc = stat.n_all_loc = n_loc;
#ifdef RECORD_CM_IN_HEADER
        PS::F64 mass = 0.0;
        for (int i=0; i<n_loc; i++) mass += system_soft[i].mass;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        stat.pcm.mass = PS::Comm::getSum(mass);
#else        
        stat.pcm.mass = mass;
#endif
        stat.pcm.pos = file_header.pos_offset;
        stat.pcm.vel = file_header.vel_offset;
        stat.pcm.is_center_shift_flag = true;
#endif        

        input_parameters.n_glb.value = n_glb;

        read_data_flag = true;
    }

    //! reading data from particle array
    /*!
      @param[in] _n_partcle: number of particles
      @param[in] _particle_array
     */
    template <class Tptcl>
    void readDataFromArray(PS::S64 _n_partcle, Tptcl* _particle) {
        assert(read_parameters_flag);

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
        for(PS::S32 i=0; i<n_loc; i++){
            int k = i_h+i;
            system_soft[i].mass  = _particle[k].getMass();
            system_soft[i].pos.x = _particle[k].getPos()[0];
            system_soft[i].pos.y = _particle[k].getPos()[1];
            system_soft[i].pos.z = _particle[k].getPos()[2];
            system_soft[i].vel.x = _particle[k].getVel()[0];
            system_soft[i].vel.y = _particle[k].getVel()[1];
            system_soft[i].vel.z = _particle[k].getVel()[2];
            system_soft[i].id = i_h+i+1;
            system_soft[i].group_data.artificial.setParticleTypeToSingle();
        }

        file_header.nfile = 0;
        file_header.n_body = n_glb;
        file_header.time = 0.0;

        stat.time = 0.0;
        stat.n_real_glb = stat.n_all_glb = n_glb;
        stat.n_real_loc = stat.n_all_loc = n_loc;
#ifdef RECORD_CM_IN_HEADER
        stat.calcAndShiftCenterOfMass(&system_soft[0], n_loc, 3, true);
        file_header.pos_offset = stat.pcm.pos;
        file_header.vel_offset = stat.pcm.vel;
#endif
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
#ifdef RECORD_CM_IN_HEADER
        stat.calcAndShiftCenterOfMass(&system_soft[0], n_loc, 3, true);

        file_header.pos_offset = stat.pcm.pos;
        file_header.vel_offset = stat.pcm.vel;
#endif

        input_parameters.n_glb.value = n_glb;

        read_data_flag = true;
    }

    //! generate data from plummer model
    void generatePlummer() {
        // ensure parameters are used
        assert(read_parameters_flag);

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
#ifdef STELLAR_EVOLUTION
            system_soft[i].radius = 0.0;
            system_soft[i].dm = 0.0;
            system_soft[i].time_record = 0.0;
            system_soft[i].time_interrupt = 0.0;
            system_soft[i].binary_state = 0;
#ifdef BSE_BASE
            system_soft[i].star.initial(mass[i]*bse_parameters.mscale.value);
#endif
#endif
            system_soft[i].group_data.artificial.setParticleTypeToSingle();
        }

        file_header.nfile = 0;
        file_header.n_body = n_glb;
        file_header.time = 0.0;

        stat.time = 0.0;
        stat.n_real_glb = stat.n_all_glb = n_glb;
        stat.n_real_loc = stat.n_all_loc = n_loc;

#ifdef RECORD_CM_IN_HEADER
        stat.calcAndShiftCenterOfMass(&system_soft[0], n_loc, 1, true);
        file_header.pos_offset = stat.pcm.pos;
        file_header.vel_offset = stat.pcm.vel;
#endif

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

#ifdef RECORD_CM_IN_HEADER
        stat.calcAndShiftCenterOfMass(&system_soft[0], n_loc, 1, true);
        file_header.pos_offset = stat.pcm.pos;
        file_header.vel_offset = stat.pcm.vel;
#endif

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
        std::cout<<std::setprecision(WRITE_PRECISION);

        // units
        if (input_parameters.unit_set.value==1) {
            input_parameters.gravitational_constant.value = G_ASTRO;
#ifdef BSE_BASE
            bse_parameters.tscale.value = 1.0; // Myr
            bse_parameters.rscale.value = PC_TO_RSUN;
            bse_parameters.mscale.value = 1.0; // Msun
            bse_parameters.vscale.value = PCMYR_TO_KMS;
#endif
            if(print_flag) {
                std::cout<<"----- Unit set 1: Msun, pc, Myr -----\n"
                         <<"gravitational_constant = "<<input_parameters.gravitational_constant.value<<" pc^3/(Msun*Myr^2)\n";
#ifdef BSE_BASE
                std::cout<<"----- Unit conversion for BSE ----- \n"
                         <<" tscale = "<<bse_parameters.tscale.value<<"  Myr / Myr\n"
                         <<" mscale = "<<bse_parameters.mscale.value<<"  Msun / Msun\n"
                         <<" rscale = "<<bse_parameters.rscale.value<<"  Rsun / pc\n"
                         <<" vscale = "<<bse_parameters.vscale.value<<"  [km/s] / [pc/Myr]\n";
#endif

            }
//#ifdef GALPY
//            galpy_parameters.setStdUnit(print_flag);
//#endif

        }

        // calculate system parameters
        PS::F64 r_in, mass_average, vel_disp;// mass_max, vel_max;
        PS::F64& r_out = input_parameters.r_out.value;
        PS::F64& r_bin = input_parameters.r_bin.value;
        PS::F64& r_search_min = input_parameters.r_search_min.value;
//        PS::F64& r_search_max = input_parameters.r_search_max.value;
        PS::F64& dt_soft = input_parameters.dt_soft.value;
        PS::F64& dt_snap = input_parameters.dt_snap.value;
        PS::F64& search_vel_factor =  input_parameters.search_vel_factor.value;
        PS::F64& ratio_r_cut   =  input_parameters.ratio_r_cut.value;
        PS::S64& n_bin         =  input_parameters.n_bin.value;
        //PS::F64& theta         =  input_parameters.theta.value;
        PS::F64& G             =  input_parameters.gravitational_constant.value;
        PS::F64& nstep_dt_soft_kepler = input_parameters.nstep_dt_soft_kepler.value;

        // local particle number
        const PS::S64 n_loc = system_soft.getNumberOfParticleLocal();

        // local c.m velocity
        PS::F64vec vel_cm_loc = 0.0;
        // local c.m. pos
        PS::F64vec pos_cm_loc = 0.0;
        // local c.m. mass
        PS::F64 mass_cm_loc = 0.0;
        // local maximum mass
        //PS::F64 mass_max_loc = 0.0;
        // box size
        PS::F64 rmax=0.0;

        for(PS::S64 i=0; i<n_loc; i++){
            PS::F64 mi = system_soft[i].mass;
            PS::F64vec vi = system_soft[i].vel;
            PS::F64vec ri = system_soft[i].pos;

#ifdef PETAR_DEBUG
            assert(mi>0);
#endif
            mass_cm_loc += mi;
            vel_cm_loc += mi * vi;
            pos_cm_loc += mi * ri;
            PS::F64 r2 = ri*ri;
            rmax = std::max(r2,rmax);
            //mass_max_loc = std::max(mi, mass_max_loc);
        }
        rmax = std::sqrt(rmax);
        PS::F64 rmax_glb = PS::Comm::getMaxValue(rmax);

        // global c.m. parameters
        PS::F64    mass_cm_glb = PS::Comm::getSum(mass_cm_loc);
        //mass_max = PS::Comm::getMaxValue(mass_max_loc);
        PS::F64vec pos_cm_glb  = PS::Comm::getSum(pos_cm_loc);
        PS::F64vec vel_cm_glb  = PS::Comm::getSum(vel_cm_loc);
        pos_cm_glb /= mass_cm_glb;
        vel_cm_glb /= mass_cm_glb;

        PS::F64 rmin_glb = std::sqrt(pos_cm_glb*pos_cm_glb);

        // local velocity square
        PS::F64 vel_sq_loc = 0.0;
        PS::S64 n_vel_loc_count = 0;

        // single particle starting index
        PS::S64 single_start_index = 0;
        const PS::S64 bin_last_id = 2*n_bin;
        if (system_soft[0].id<bin_last_id) {
            single_start_index = std::min(bin_last_id - system_soft[0].id + 1,n_loc);
            if(single_start_index%2!=0) single_start_index--;
        }
        // binary particle starting index
        const PS::S64 binary_start_index = (system_soft[0].id%2==0)?1:0;

        // calculate velocity dispersion
        for (PS::S64 i=binary_start_index; i<single_start_index; i+=2) {
            PS::F64 m1 = system_soft[i].mass;
            PS::F64 m2 = system_soft[i+1].mass;
            PS::F64vec dv = (m1*system_soft[i].vel + m2*system_soft[i+1].vel)/(m1+m2) - vel_cm_glb;
            vel_sq_loc += dv * dv;
            n_vel_loc_count++;
        }
    
        if (single_start_index <n_loc) 
            for (PS::S64 i=single_start_index; i<n_loc; i++){
                PS::F64vec dv = system_soft[i].vel - vel_cm_glb;
                vel_sq_loc += dv * dv;
                n_vel_loc_count++;
            }

        const PS::S64    n_vel_glb_count= PS::Comm::getSum(n_vel_loc_count);
        const PS::S64    n_glb          = PS::Comm::getSum(n_loc);
        const PS::F64    vel_sq_glb     = PS::Comm::getSum(vel_sq_loc);
        vel_disp   = sqrt(vel_sq_glb / 3.0 / (PS::F64)n_vel_glb_count);

        PS::F64 mass_average_glb = mass_cm_glb/(PS::F64)n_glb;
        mass_average = mass_average_glb;

        // flag to check whether r_ous is already defined
        bool r_out_flag = (r_out>0);
    
        // if r_out is already defined, calculate r_in based on ratio_r_cut
        if (r_out_flag) r_in = r_out * ratio_r_cut;
        // calculate r_out based on virial radius scaled with (N)^(1/3), calculate r_in by ratio_r_cut
        else {
            if (n_glb>1) {
                r_out = std::min(0.1*G*mass_cm_glb/(std::pow(n_glb,1.0/3.0)) / (3*vel_disp*vel_disp), 3.0*(rmax_glb-rmin_glb));
                r_in = r_out * ratio_r_cut;
            }
            else {
                // give two small values, no meaning at all
                r_out = 1e-16;
                r_in = r_out*ratio_r_cut;
                if (print_flag) std::cout<<"In one particle case, changeover radius is set to a small value\n";
            }
        }

        // if tree time step is not defined, calculate tree time step by r_out and velocity dispersion
        if (dt_soft==0.0) {
            if (n_glb==1) {
                if (print_flag) std::cout<<"In one particle case, tree time step is finishing - starting time\n";
                dt_soft = input_parameters.time_end.value - stat.time;
            }
            else {
                // 1/nstep of a binary period with semi-major axis = r_int.
                if (nstep_dt_soft_kepler>0)  
                    dt_soft = regularTimeStep(COMM::Binary::semiToPeriod(r_in, mass_average, G)/nstep_dt_soft_kepler);
                else 
                    dt_soft = regularTimeStep(0.1*r_out / vel_disp);
            }
        }
        else {
            dt_soft = regularTimeStep(dt_soft);
            // if r_out is not defined, adjust r_out to minimum based on tree step
            if (!r_out_flag) {
                if (n_glb>1) {
                    if (nstep_dt_soft_kepler>0) {
                            r_in = COMM::Binary::periodToSemi(dt_soft*nstep_dt_soft_kepler, mass_average, G);
                            r_out = r_in / ratio_r_cut;
                    }
                    else {
                        r_out = 10.0*dt_soft*vel_disp;
                        r_in = r_out * ratio_r_cut;
                    }
                }
                else {
                    r_out = 1e-16;
                    r_in = r_out*ratio_r_cut;
                    if (print_flag) std::cout<<"In one particle case, changeover radius is set to a small value\n";
                }
            }
        }

        // if r_bin is not defined, set to theta * r_in
        if (r_bin==0.0) r_bin = 0.8*r_in;

        // if r_search_min is not defined, calculate by search_vel_factor*velocity_dispersion*tree_time_step + r_out
        if (r_search_min==0.0) r_search_min = search_vel_factor*vel_disp*dt_soft + r_out;
        // if r_search_max is not defined, calcualte by 5*r_out
//        if (r_search_max==0.0) r_search_max = 5*r_out;
        // calculate v_max based on r_search_max, tree time step and search_vel_factor
        //vel_max = (r_search_max - r_out) / dt_soft / search_vel_factor;

        // regularize output time to be integer times of dt_soft
        if (dt_snap<dt_soft) 
            dt_snap = dt_soft;
        else
            dt_snap = int(dt_snap/dt_soft)*dt_soft;

        EPISoft::eps   = input_parameters.eps.value;
        EPISoft::r_out = r_out;
        ForceSoft::grav_const = input_parameters.gravitational_constant.value;
        Ptcl::search_factor = search_vel_factor;
        Ptcl::r_search_min = r_search_min;
        Ptcl::mean_mass_inv = 1.0/mass_average;
        Ptcl::r_group_crit_ratio = r_bin/r_in;
        //Ptcl::vel_cm = vel_cm_glb;
        escaper.r_escape_sq = input_parameters.r_escape.value*input_parameters.r_escape.value;
        escaper.check_energy_flag = (input_parameters.r_escape.value>=0);

        if(print_flag) {
        // set print format
            std::cout<<"----- Parameter list: -----\n";
            std::cout<<" Average mass                      = "<<mass_average   <<std::endl
                     <<" Mean inner changeover radius      = "<<r_in           <<std::endl
                     <<" Mean outer changeover radius      = "<<r_out          <<std::endl
                     <<" Mean SDAR group detection radius  = "<<r_bin          <<std::endl
                     <<" Minimum neighbor searching radius = "<<r_search_min   <<std::endl
                     <<" Velocity dispersion               = "<<vel_disp       <<std::endl
                     <<" Tree time step                    = "<<dt_soft        <<std::endl
                     <<" Output time step                  = "<<dt_snap        <<std::endl;
        }

        // check restart
        bool restart_flag = file_header.nfile; // nfile = 0 is assumed as initial data file

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

                // calculate r_search for particles, for binary, r_search depend on vel_disp
                if(id<=2*n_bin) system_soft[i].r_search = std::max(r_search_min,vel_disp*dt_soft*search_vel_factor + system_soft[i].changeover.getRout());
                else system_soft[i].calcRSearch(dt_soft);

                auto& pi_cm = system_soft[i].group_data.cm;
                pi_cm.mass = pi_cm.vel.x = pi_cm.vel.y = pi_cm.vel.z = 0.0;
            }
        }
        else {
            // clear up group_data.cm to avoid issue in search neighbor; update changeover and r_search
#pragma omp parallel for
            for (PS::S32 i=0; i<stat.n_real_loc; i++) {
                auto& pi_cm = system_soft[i].group_data.cm;

                // update changeover
                if (input_parameters.update_changeover_flag) {
                    PS::F64 m_fac = system_soft[i].mass*Ptcl::mean_mass_inv;
                    system_soft[i].changeover.setR(m_fac, r_in, r_out);
                }

                // update research
                if (input_parameters.update_rsearch_flag) {
                    if (pi_cm.mass!=0.0) {
//#ifdef GROUP_DATA_WRITE_ARTIFICIAL
                        system_soft[i].r_search = std::max(r_search_min, vel_disp*dt_soft*search_vel_factor + system_soft[i].changeover.getRout());
// not correct, when data is dumped, it is not cm information but artificial particle data
//#else
//                        PS::F64 vcm = std::sqrt(pi_cm.vel*pi_cm.vel);
//                        system_soft[i].r_search = std::max(r_search_min,vcm*dt_soft*search_vel_factor + system_soft[i].changeover.getRout());
//#endif
                    }
                    else system_soft[i].calcRSearch(dt_soft);
                }

                pi_cm.mass = pi_cm.vel.x = pi_cm.vel.y = pi_cm.vel.z = 0.0;
            }
#ifdef RECORD_CM_IN_HEADER
            stat.calcAndShiftCenterOfMass(&system_soft[0], stat.n_real_loc);
#else
            stat.calcCenterOfMass(&system_soft[0], stat.n_real_loc);
#endif
            //Ptcl::vel_cm = stat.pcm.vel;

#ifdef BSE_BASE
            // set randseeds filename to read later
            std::string &data_filename = input_parameters.fname_inp.value;
            rand_parameters.seedfile.value = data_filename + ".randseeds";
#endif
        }

#ifdef GALPY
        std::string galpy_conf_filename = input_parameters.fname_inp.value+".galpy";
        galpy_manager.initial(galpy_parameters, stat.time, galpy_conf_filename, restart_flag, print_flag);
#endif
    
        // set system hard paramters
        hard_manager.setDtRange(input_parameters.dt_soft.value/input_parameters.dt_limit_hard_factor.value, input_parameters.dt_min_hermite_index.value);
        hard_manager.setEpsSq(input_parameters.eps.value*input_parameters.eps.value);
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
#ifdef ORBIT_SAMPLING
        hard_manager.ap_manager.orbit_manager.setParticleSplitN(input_parameters.n_split.value);
#endif
        hard_manager.h4_manager.step.eta_4th = input_parameters.eta.value;
        hard_manager.h4_manager.step.eta_2nd = 0.01*input_parameters.eta.value;
        hard_manager.h4_manager.step.calcAcc0OffsetSq(mass_average, r_out, input_parameters.gravitational_constant.value);
        hard_manager.ar_manager.energy_error_relative_max = input_parameters.e_err_ar.value;
        hard_manager.ar_manager.step_count_max = input_parameters.step_limit_ar.value;
        hard_manager.ar_manager.step.initialSymplecticCofficients(-6);
        hard_manager.ar_manager.slowdown_pert_ratio_ref = input_parameters.sd_factor.value;
        hard_manager.ar_manager.slowdown_timescale_max = dt_soft*input_parameters.n_step_per_orbit.value;
        //hard_manager.ar_manager.slowdown_timescale_max = dt_soft;
#ifdef SLOWDOWN_MASSRATIO
        hard_manager.ar_manager.slowdown_mass_ref = mass_average;
#endif
#ifdef STELLAR_EVOLUTION
        hard_manager.ar_manager.interaction.interrupt_detection_option = input_parameters.interrupt_detection_option.value;
#ifdef BSE_BASE
        hard_manager.ar_manager.interaction.stellar_evolution_option = input_parameters.stellar_evolution_option.value;
        if (write_style) hard_manager.ar_manager.interaction.stellar_evolution_write_flag = true;
        else hard_manager.ar_manager.interaction.stellar_evolution_write_flag = false;
        if (input_parameters.stellar_evolution_option.value>0) {
            hard_manager.ar_manager.interaction.bse_manager.initial(bse_parameters, print_flag);
            hard_manager.ar_manager.interaction.tide.speed_of_light = hard_manager.ar_manager.interaction.bse_manager.getSpeedOfLight();
        }
        rand_manager.initialAll(rand_parameters);
        rand_manager.printRandSeeds(std::cout);


        // initial stellar evolution for each star
        if (!restart_flag) {
#pragma omp parallel for
            for (PS::S32 i=0; i<stat.n_real_loc; i++) {
                auto& pi = system_soft[i];
                hard_manager.ar_manager.interaction.modifyOneParticle(pi, stat.time, stat.time);
            }
        }

#endif
#endif
#ifdef ADJUST_GROUP_PRINT
        // group information
        if (write_style&&input_parameters.adjust_group_write_option.value==1) 
            hard_manager.h4_manager.adjust_group_write_flag=true;
        else 
            hard_manager.h4_manager.adjust_group_write_flag=false;
#endif

#ifdef EXTERNAL_HARD
#ifdef GALPY
        hard_manager.h4_manager.interaction.ext_force.initial(external_hard_parameters, galpy_manager, stat, print_flag);
#else
        hard_manager.h4_manager.interaction.ext_force.initial(external_hard_parameters, stat.time, print_flag);
#endif
        hard_manager.ar_manager.interaction.ext_force = &hard_manager.h4_manager.interaction.ext_force;
#endif        

        // record id range
        hard_manager.record_id_range.id_start_one = input_parameters.record_id_start_one.value;
        hard_manager.record_id_range.id_end_one = input_parameters.record_id_end_one.value;
        hard_manager.record_id_range.id_start_two = input_parameters.record_id_start_two.value;
        hard_manager.record_id_range.id_end_two = input_parameters.record_id_end_two.value;

        // check consistence of paramters
        input_parameters.checkParams();
        hard_manager.checkParams();

        // initial hard class and parameters
        system_hard_one_cluster.manager = &hard_manager;
        system_hard_one_cluster.setTimeOrigin(stat.time);

        system_hard_isolated.manager = &hard_manager;
        system_hard_isolated.setTimeOrigin(stat.time);

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        system_hard_connected.manager = &hard_manager;
        system_hard_connected.setTimeOrigin(stat.time);
#endif

        time_kick = stat.time;

        if(write_style>0&&my_rank==0) {
            if (print_flag) std::cout<<"-----  Dump parameter files -----"<<std::endl;
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

#ifdef BSE_BASE
            // save bse parameters
            std::string fbse_par = input_parameters.fname_par.value + fbse_par_suffix;
            if (print_flag) std::cout<<"Save bse_parameters to file "<<fbse_par<<std::endl;
            if( (fpar_out = fopen(fbse_par.c_str(),"w")) == NULL) {
                fprintf(stderr,"Error: Cannot open file %s.\n", fbse_par.c_str());
                abort();
            }
            bse_parameters.input_par_store.writeAscii(fpar_out);
            fclose(fpar_out);

            // save random parameters
            std::string frand_par = input_parameters.fname_par.value + ".rand";
            if (print_flag) std::cout<<"Save rand_parameters to file "<<frand_par<<std::endl;
            if( (fpar_out = fopen(frand_par.c_str(),"w")) == NULL) {
                fprintf(stderr,"Error: Cannot open file %s.\n", frand_par.c_str());
                abort();
            }
            rand_parameters.input_par_store.writeAscii(fpar_out);
            fclose(fpar_out);
#endif

#ifdef GALPY
            // save galpy parameters
            std::string fgalpy_par = input_parameters.fname_par.value + ".galpy";
            if (print_flag) std::cout<<"Save galpy_parameters to file "<<fgalpy_par<<std::endl;
            if( (fpar_out = fopen(fgalpy_par.c_str(),"w")) == NULL) {
                fprintf(stderr,"Error: Cannot open file %s.\n", fgalpy_par.c_str());
                abort();
            }
            galpy_parameters.input_par_store.writeAscii(fpar_out);
            fclose(fpar_out);
#endif

#ifdef EXTERNAL_HARD
            // save galpy parameters
            std::string fexthard_par = input_parameters.fname_par.value + ".exthard";
            if (print_flag) std::cout<<"Save external_hard_parameters to file "<<fexthard_par<<std::endl;
            if( (fpar_out = fopen(fexthard_par.c_str(),"w")) == NULL) {
                fprintf(stderr,"Error: Cannot open file %s.\n", fexthard_par.c_str());
                abort();
            }
            external_hard_parameters.input_par_store.writeAscii(fpar_out);
            fclose(fpar_out);
#endif
        }

        // initial tree step manager
        dt_manager.setKDKMode();

        if (print_flag) std::cout<<"-----  Finish parameter initialization -----"<<std::endl;

        initial_parameters_flag = true;
    }

    //! the first initial step to find groups and energy calculation
    void initialStep() {
        assert(initial_parameters_flag);
        if (initial_step_flag) return;

        assert(checkTimeConsistence());

        // one particle case
        if (stat.n_real_glb==1) {
            Ptcl::group_data_mode = GroupDataMode::artificial;

            if (stat.n_real_loc==1) system_soft[0].clearForce();

            /// force from external potential
            externalForce();

            // initial status and energy
            updateStatus(true);

            PS::F64 dt_tree = input_parameters.dt_soft.value;
            dt_manager.setStep(dt_tree);

            // output initial data
            file_header.nfile--; // avoid repeating files
            output();

#ifdef PROFILE
            clearProfile();
#endif
            initial_step_flag = true;

            return;
        }

        // domain decomposition
        domainDecompose();

        // exchange particles
        exchangeParticle();

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
        treeSoftForce();

        /// force from external potential
        externalForce();

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
        file_header.nfile--; // avoid repeating files
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

    //! integrate one particle with modification function
    PS::S32 integrateOneToTime(const PS::F64 _time_break=0.0) {
        // ensure it is initialized
        assert(initial_step_flag);
        assert(stat.n_real_glb==1);

        // check time break
        PS::F64 time_break = _time_break==0.0? input_parameters.time_end.value: std::min(_time_break,input_parameters.time_end.value);
        if (stat.time>=time_break) return 0;

        //PS::F64 dt = std::min(input_parameters.dt_soft.value, time_break-stat.time);
        PS::F64 dt_output = input_parameters.dt_snap.value;
        //bool start_flag=true;

        if (stat.n_real_loc==1) { 
            assert(system_soft.getNumberOfParticleLocal()==1);
            auto& p = system_soft[0];

            while (true) {
#ifdef PROFILE
                profile.total.start();
#endif
                
                // reset total potential
                p.clearForce();

                /// force from external potential and kick
                externalForce();

#ifdef EXTERNAL_HARD
                /// force from external hard
                if (hard_manager.h4_manager.interaction.ext_force.mode>0) {
                    H4::ForceH4 f;
                    hard_manager.h4_manager.interaction.ext_force.calcAccJerkExternal(f, p);
                    p.acc[0] += f.acc0[0];
                    p.acc[1] += f.acc0[1];
                    p.acc[2] += f.acc0[2];
                }
#endif

#ifdef RECORD_CM_IN_HEADER
                stat.calcAndShiftCenterOfMass(&p, stat.n_real_loc);
#endif

                bool interrupt_flag = false;  // for interrupt integration when time reach end
                bool output_flag = false;    // for output snapshot and information

                PS::F64 dt_kick, dt_drift;
                // for initial the system
                if (dt_manager.isNextStart()) {

                    // set step to the begining step
                    dt_kick = dt_manager.getDtStartContinue();
                }
                else {
                    // increase loop counter
                    n_loop++;

                    // for next kick-drift pair
                    dt_manager.nextContinue();
                    
                    if (dt_manager.isNextEndPossible()) {

                        // output step, get last kick step
                        output_flag = (fmod(stat.time, dt_output) == 0.0);

                        // check interruption
                        interrupt_flag = (stat.time>=time_break);

                        if (output_flag||interrupt_flag) dt_kick = dt_manager.getDtEndContinue();
                        else dt_kick = dt_manager.getDtKickContinue();
                    }
                    else dt_kick = dt_manager.getDtKickContinue();
                }
                
                //kick 
                p.vel += p.acc * dt_kick;
#ifdef GALPY
                galpy_manager.kickMovePot(dt_kick);
#endif
                time_kick += dt_kick;

                // output information
                if(output_flag) {

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

                // interrupt
                if(interrupt_flag) {
                    assert(checkTimeConsistence());
#ifdef PROFILE
                    profile.total.barrier();
                    PS::Comm::barrier();
                    profile.total.end();
#endif
                    return 0;
                }

                // second kick if output exists 
                if(output_flag) {
                    dt_kick = dt_manager.getDtStartContinue();
                    p.vel += p.acc * dt_kick;
#ifdef GALPY
                    galpy_manager.kickMovePot(dt_kick);
#endif
#ifdef RECORD_CM_IN_HEADER
                    // correct Ptcl:vel_cm
                    correctPtclVelCM(dt_kick);
#endif
                    time_kick += dt_kick;
                }

                // get drift step
                dt_drift = dt_manager.getDtDriftContinue();

#ifdef PETAR_USE_MPFRC
                mprealVec pos_mp(p.pos, p.pos_high);
                pos_mp += p.vel * dt_drift;
                pos_mp.split(p.pos, p.pos_high);
#else
                p.pos += p.vel * dt_drift;
#endif
                // drift cm
                stat.pcm.pos += stat.pcm.vel*dt_drift;

#ifdef GALPY
                galpy_manager.driftMovePot(dt_drift);
#endif

#ifdef STELLAR_EVOLUTION
                PS::F64 mbk = p.mass;
                PS::F64vec vbk = p.vel; //back up velocity in case of change
                int modify_flag = hard_manager.ar_manager.interaction.modifyOneParticle(p, stat.time, stat.time + dt_drift);
                if (modify_flag) {
                    auto& v = p.vel;
                    PS::F64 de_kin = 0.5*(p.mass*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]) - mbk*(vbk[0]*vbk[0]+vbk[1]*vbk[1]+vbk[2]*vbk[2]));
                    stat.energy.de_sd_change_cum += de_kin;
                    stat.energy.de_sd_change_modify_single += de_kin;
                    stat.energy.de_change_cum += de_kin;
                    stat.energy.de_change_modify_single += de_kin;
                    stat.energy.etot_ref += de_kin;
                    stat.energy.etot_sd_ref += de_kin;
                }
                // shift time interrupt in order to get consistent time for stellar evolution in the next drift
                //p.time_record    -= dt;
                //p.time_interrupt -= dt;
#endif
                stat.time += dt_drift;

#ifdef PROFILE
                profile.total.barrier();
                PS::Comm::barrier();
                profile.total.end();

                calcProfile();
#endif
            }
        }
        else {
            PS::F64 dt_tree = dt_manager.getStep();
            while (stat.time  < time_break) {
                stat.time += dt_tree;
                time_kick = stat.time;
            }
        }

        return 0;
    }

    
    //! integrate the system
    /*! @param[in] _time_break: additional breaking time to interrupt the integration, in default (0.0) the system integrate to time_end 
      \return interrupted cluster number
     */
    PS::S32 integrateToTime(const PS::F64 _time_break=0.0) {

        // ensure it is initialized
        assert(initial_step_flag);

        // for one particle case
        if (stat.n_real_glb==1) return integrateOneToTime(_time_break);

        // check time break
        PS::F64 time_break = _time_break==0.0? input_parameters.time_end.value: std::min(_time_break,input_parameters.time_end.value);
        if (stat.time>=time_break) return 0;

        PS::F64 dt_output = input_parameters.dt_snap.value;
        PS::F64 dt_tree = dt_manager.getStep();

        /// Main loop
        while(true) {
#ifdef PROFILE
            profile.total.start();
#endif

#ifdef STELLAR_EVOLUTION
            // correct soft potential energy due to mass change
            correctSoftPotMassChange();
#endif

            // remove artificial and ununsed particles in system_soft.
            removeParticles();

#ifdef RECORD_CM_IN_HEADER
            // update center
            stat.calcAndShiftCenterOfMass(&system_soft[0], stat.n_real_loc);
#endif

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

            /// force from external potential
            externalForce();

            // >5 correct change over and potential energy due to mass change
            /// correct system_soft.acc with changeover, using system_hard and system_soft particles
            /// substract tidal tensor measure point force
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
            
#ifdef STELLAR_EVOLUTION
            hard_manager.ar_manager.interaction.time_interrupt_max = stat.time + dt_drift;
#endif            
            
            drift(dt_drift);
            // update stat time 
            stat.time = system_hard_one_cluster.getTimeOrigin();

#ifdef PROFILE
            // calculate profile
            profile.total.barrier();
            PS::Comm::barrier();
            profile.total.end();

            calcProfile();
#endif
            
        }

        return 0;
    }

    void clear() {

        if (fstatus.is_open()) fstatus.close();
        if (fesc.is_open()) fesc.close();
#ifdef PROFILE
        if (fprofile.is_open()) fprofile.close();
#endif

#ifdef STELLAR_EVOLUTION
#ifdef BSE_BASE
        auto& interaction = hard_manager.ar_manager.interaction;
        if (interaction.fout_sse.is_open()) interaction.fout_sse.close();
        if (interaction.fout_bse.is_open()) interaction.fout_bse.close();
#else
        auto& interaction = hard_manager.ar_manager.interaction;
        if (interaction.fout_interrupt.is_open()) interaction.fout_interrupt.close();
#endif
#endif
#ifdef ADJUST_GROUP_PRINT
        if (hard_manager.h4_manager.fgroup.is_open()) hard_manager.h4_manager.fgroup.close();
#endif
        if (pos_domain) {
            delete[] pos_domain;
            pos_domain=NULL;
        }

        //if (initial_fdps_flag) PS::Finalize();
        remove_list.resizeNoInitialize(0);
        n_interrupt_glb = 0;
        //initial_fdps_flag = false;
        read_parameters_flag = false;
        read_data_flag = false;
        initial_parameters_flag = false;
        initial_step_flag = false;
    }

    ~PeTar() { clear();}
};

bool PeTar::initial_fdps_flag = false;
