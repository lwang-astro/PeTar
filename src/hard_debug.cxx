#include <iostream>
#include <iomanip>
#include <string>
#include<getopt.h>

#include <particle_simulator.hpp>
#define HARD_DEBUG_PRINT_FEQ 1024

#include "io.hpp"
#include "hard_assert.hpp"
#include "cluster_list.hpp"
#include "hard.hpp"
#include "soft_ptcl.hpp"
#include "static_variables.hpp"

#ifdef BSE_BASE
#include "../parallel-random/rand_interface.hpp"
#endif

int main(int argc, char **argv){
  int mode=0; // 0: integrate to time; 1: times stability
  PS::F64 tstart = -1;
  PS::F64 tend = -1;
  PS::F64 slowdown_factor=0;
  PS::F64 eta_4th=0;
  PS::F64 eta_2nd=0;
  PS::F64 e_err_ar = -1;
  PS::F64 e_err_hard = 1e-4;
  PS::S32 dt_min_power = -1;
  PS::F64 dt_max = -1;
  PS::F64 ds_scale = -1.0;
  PS::S32 step_arc_limit = 100000;
  PS::S32 n_crit = 0;
  PS::S32 istart = -1;
  PS::S32 iend = -1;
  std::string filename="hard_dump";
  std::string fhardpar="input.par.hard";
#ifdef BSE_BASE
  uint64_t seed = 0;
  std::string bse_name = BSEManager::getBSEName();
  std::string fsse_suffix = BSEManager::getSSEOutputFilenameSuffix();
  std::string fbse_suffix = BSEManager::getBSEOutputFilenameSuffix();

  std::string fbsepar = "input.par" + fbse_suffix;
#endif
#ifdef STELLAR_EVOLUTION
  int stellar_evolution_option = -1;
#endif
#ifdef SOFT_PERT
  bool soft_pert_flag=true;
#endif

  int copt;
  int option_index;
  static int opt_flag = -1;
  static struct option long_options[] = {
      {"hermite-dt-max",        required_argument, &opt_flag, 0},
      {"hermite-dt-min-power",  required_argument, &opt_flag, 1},
      {"energy-err-ar",         required_argument, &opt_flag, 2},
#ifdef HARD_CHECK_ENERGY
      {"energy-err-hard",       required_argument, &opt_flag, 3},
#endif
      {"slowdown-factor",   required_argument, &opt_flag, 4},
      {"step-limit-ar",     required_argument, &opt_flag, 5},
      {"step-scale-ar",     required_argument, &opt_flag, 6},
      {"hermite-eta-4th",   required_argument, &opt_flag, 7},
      {"hermite-eta-2nd",   required_argument, &opt_flag, 8},
#ifdef STELLAR_EVOLUTION
      {"stellar-evolution", required_argument, &opt_flag, 9},
#endif
#ifdef BSE_BASE
      {"rand-seed",         required_argument, &opt_flag, 10},
#endif
      {"tstart",            required_argument, &opt_flag, 11},
      {"tend",              required_argument, &opt_flag, 12},
      {"istart",            required_argument, &opt_flag, 13},
      {"iend",              required_argument, &opt_flag, 14},
      {"help",        no_argument, 0, 'h'},        
      {0,0,0,0}
  };

  while ((copt = getopt_long(argc, argv, "m:n:b:p:h", long_options, &option_index)) != -1)
    switch (copt) {
    case 0:
        switch (opt_flag) {
        case 0:
            dt_max = atof(optarg);
            break;
        case 1:
            dt_min_power = atoi(optarg);
            break;
        case 2:
            e_err_ar = atof(optarg);
            break;
#ifdef HARD_CHECK_ENERGY
        case 3:
            e_err_hard = atof(optarg);
            break;
#endif
        case 4:
            slowdown_factor = atof(optarg);
            break;
        case 5:
            step_arc_limit = atoi(optarg);
            break;
        case 6:
            ds_scale = atof(optarg);
            break;
        case 7:
            eta_4th = atof(optarg);
            break;
        case 8:
            eta_2nd = atof(optarg);
            break;
#ifdef STELLAR_EVOLUTION
        case 9:
            stellar_evolution_option = atoi(optarg);
            break;
#endif
#ifdef BSE_BASE
        case 10:
            seed = atoi(optarg);
            break;
#endif
        case 11:
            tstart = atof(optarg);
            break;
        case 12:
            tend = atof(optarg);
            break;
        case 13:
            istart = atoi(optarg);
            break;
        case 14:
            iend = atoi(optarg);
            break;
        default:
            break;
        }
        break;
    case 'm':
        mode = atoi(optarg);
        break;
    case 'n':
        n_crit = atoi(optarg);
        break;
    case 'p':
        fhardpar = optarg;
        break;
#ifdef SOFT_PERT
    case 'S':
        soft_pert_flag=false;
        break;
#endif
#ifdef BSE_BASE
    case 'b':
        fbsepar = optarg;
        break;
#endif
    case 'h':
        std::cout<<"A tool to integrate a dumped cluster of neighbor particles using particle-particle method (Hermite/SDAR)\n"
                 <<"Usage: petar.hard.debug [options] [hard parameter filename (defaulted: input.par.hard)] [dumped data filename (defaulted: hard_dump)]\n"
                 <<"   dumped data file: Hard dump file from petar simulation\n"
                 <<"options (default values shown at the end):\n"
                 <<"    -m [int]:     running mode: 0: evolve system to time_end; 1: stability check: "<<mode<<std::endl
                 <<"    -n [int]:     if >0, only do integration when particle number matches the given value: "<<n_crit<<std::endl
                 <<"    -p [string]:  hard parameter file name: "<<fhardpar<<std::endl
#ifdef BSE_BASE
                 <<"    -b [string]:  bse parameter file name: "<<fbsepar<<std::endl
#endif
#ifdef SOFT_PERT
                 <<"    -S:           suppress soft perturbation (tidal tensor)\n"
#endif
                 <<"    -h (--help):  help"<<std::endl
                 <<"long options (if no default values, use values from input.par.hard):\n"
                 <<"        --tstart            [double]:  if >0 only do integration when physical time >= tstart: "<<tstart<<std::endl
                 <<"        --tend              [double]:  if >0 only do integration when physical time < tend: "<<tend<<std::endl
                 <<"        --istart            [int]:     if >0 only do integration when dump index >= istart: "<<istart<<std::endl
                 <<"        --iend              [int]:     if >0 only do integration when dump index < iend: "<<iend<<std::endl
#ifdef HARD_CHECK_ENERGY
                 <<"        --energy-err-hard   [double]:  hard energy limit\n"
#endif
                 <<"        --energy-err-ar     [double]:  AR energy limit \n"
                 <<"        --hermite-dt-max    [double]:  hard time step max (should use together with -d)\n"
                 <<"        --hermite-dt-min-power [int]:  hard time step min power (should use together with -D)\n"
                 <<"        --hermite-eta-4th   [double]:  Eta 4th for hermite \n"
                 <<"        --hermite-eta-2nd   [double]:  Eta 2nd for hermite \n"
#ifdef BSE_BASE
                 <<"        --rand-seed         [int]:     random seed to generate kick velocity\n"
#endif
#ifdef STELLAR_EVOLUTION
                 <<"        --stellar-evolution [int]:     Stellar evolution option: \n"
#endif
                 <<"        --slowdown-factor   [double]:  change slowdown factor reference\n"
                 <<"        --step-limit-ar     [int]:     AR step count limit\n"
                 <<"        --step-scale-ar     [double]:  AR step scaling factor\n";
        return 0;
    default:
        std::cerr<<"Unknown argument. check '-h' for help.\n";
        abort();
    }

  if (optind<argc) {
      filename=argv[argc-1];
  }

  std::cerr<<"Reading dump file:"<<filename<<std::endl;
  std::cerr<<"Hard manager parameter file:"<<fhardpar<<std::endl;

  std::cout<<std::setprecision(WRITE_PRECISION);

  HardManager hard_manager;
  FILE* fpar_in;
  if( (fpar_in = fopen(fhardpar.c_str(),"r")) == NULL) {
      fprintf(stderr,"Error: Cannot open file %s.\n", fhardpar.c_str());
      abort();
  }
  hard_manager.readBinary(fpar_in);
  fclose(fpar_in);

#ifdef STELLAR_EVOLUTION
#ifdef BSE_BASE
  if (stellar_evolution_option>=0) {
      hard_manager.ar_manager.interaction.stellar_evolution_option = stellar_evolution_option;
      if (stellar_evolution_option==0) 
          hard_manager.ar_manager.interaction.stellar_evolution_write_flag = false;
  }
  IOParamsBSE bse_io;
  std::cerr<<bse_name<<" parameter file:"<<fbsepar<<std::endl;
  if( (fpar_in = fopen(fbsepar.c_str(),"r")) == NULL) {
      fprintf(stderr,"Error: Cannot open file %s.\n", fbsepar.c_str());
      abort();
  }
  bse_io.input_par_store.readAscii(fpar_in);
  fclose(fpar_in);
  hard_manager.ar_manager.interaction.bse_manager.initial(bse_io);
  hard_manager.ar_manager.interaction.tide.gravitational_constant = hard_manager.ar_manager.interaction.gravitational_constant;
  hard_manager.ar_manager.interaction.tide.speed_of_light = hard_manager.ar_manager.interaction.bse_manager.getSpeedOfLight();

  if (hard_manager.ar_manager.interaction.stellar_evolution_write_flag) {
      hard_manager.ar_manager.interaction.fout_sse.open((filename+fsse_suffix).c_str(), std::ofstream::out);
      hard_manager.ar_manager.interaction.fout_bse.open((filename+fbse_suffix).c_str(), std::ofstream::out);
      hard_manager.ar_manager.interaction.fout_sse<<std::setprecision(WRITE_PRECISION);
      hard_manager.ar_manager.interaction.fout_bse<<std::setprecision(WRITE_PRECISION);      
  }
#endif // BSE_BASE
#endif //STELLAR_EVOLUTION

#ifdef ADJUST_GROUP_PRINT
  if (hard_manager.h4_manager.adjust_group_write_flag) {
      hard_manager.h4_manager.fgroup.open((filename+".group").c_str(), std::ofstream::out);
      hard_manager.h4_manager.fgroup<<std::setprecision(WRITE_PRECISION);
  }
#endif

#ifdef HARD_CHECK_ENERGY
  // Set hard energy limit
  if (e_err_hard>0 ) {
      std::cerr<<"New hard energy error max: "<<e_err_hard<<std::endl;
      hard_manager.energy_error_max = e_err_hard;
  }
#endif

  // Set step limit for ARC sym
  if (step_arc_limit>0) {
      std::cerr<<"New AR step count max: "<<step_arc_limit<<std::endl;
      hard_manager.ar_manager.step_count_max = step_arc_limit;
  }

  // set slowdown factor
  if(slowdown_factor>0) {
      std::cerr<<"New slowdown perturbation ratio reference: "<<slowdown_factor<<std::endl;
      hard_manager.ar_manager.slowdown_pert_ratio_ref = slowdown_factor;
  }

  // set eta
  if(eta_4th>0) {
      std::cerr<<"New eta_4th: "<<eta_4th<<std::endl;
      hard_manager.h4_manager.step.eta_4th = eta_4th;
  }

  if(eta_2nd>0) {
      std::cerr<<"New eta_2nd: "<<eta_2nd<<std::endl;
      hard_manager.h4_manager.step.eta_2nd = eta_2nd;
  }

  // time step
  if(dt_min_power>0&&dt_max>0) {
      std::cerr<<"New time step max: "<<dt_max<<"   min power index: "<<dt_min_power<<std::endl;
      hard_manager.setDtRange(dt_max, dt_min_power);
  }

  if(ds_scale>0) {
      std::cerr<<"New ds scale: "<<ds_scale<<std::endl;
      hard_manager.ar_manager.ds_scale = ds_scale;
  }

  if(e_err_ar>0) {
      std::cerr<<"New AR relative energy error maximum: "<<e_err_ar<<std::endl;
      hard_manager.ar_manager.energy_error_relative_max = e_err_ar;
  }

  hard_manager.checkParams();
  hard_manager.print(std::cerr);

  std::FILE* fp = std::fopen(filename.c_str(),"r");
  if (fp==NULL) {
      std::cerr<<"Error: filename "<<filename<<" cannot be open!\n";
      abort();
  }

  HardDump hard_dump;
  int ncount = 0;
  
  while (true) {
      hard_dump.readOneClusterBinary(fp);
      if (feof(fp)) break;

#ifdef BSE_BASE
      hard_dump.rand_manager.readRandSeedLocalBinary(fp);
      if (seed!=0) hard_dump.rand_manager.initialFromSeed(seed);
#endif

      ncount++;

      // skip if particle number not match n_crit when n_crit>0
      if (n_crit>0 && hard_dump.n_ptcl != n_crit) continue;
      if (tstart>0 && hard_dump.time_offset < tstart) continue;
      if (tend>0 && hard_dump.time_offset >= tend) continue;

      std::cerr<<"Dump "<<ncount<<"\nTime: "<<hard_dump.time_offset<<std::endl;
#ifdef BSE_BASE
      hard_dump.rand_manager.printRandSeeds(std::cerr);
#endif

#ifdef SOFT_PERT
      if (!soft_pert_flag) {
          std::cerr<<"Suppress soft perturbation\n";
          if (hard_dump.n_group>0) {
              // if no artificial particles, continue
              if (hard_dump.ptcl_arti_bk.getPointer()!=NULL) {
                  // set all tidal tensor force to zero
                  for (int i=0; i<hard_dump.n_group; i++) {
                      int offset= i*hard_manager.ap_manager.getArtificialParticleN();
                      auto* pi = &(hard_dump.ptcl_arti_bk[offset]);
                      //auto* pcm = ap_manager.getCMParticles(pi);
                      auto* ptt = hard_manager.ap_manager.getTidalTensorParticles(pi);
                      for (int j=0; j<hard_manager.ap_manager.getTidalTensorParticleN(); j++) {
                          ptt[j].acc = PS::F64vec(0.0);
                      }
                  }
              }
          }
      }
#endif

      // running mode
      if (mode==0) {
          //SystemHard sys;
          //sys.manager = &hard_manager;
          //sys.allocateHardIntegrator();

          // change ARC parameters
          //sys.driveForMultiClusterImpl(hard_dump.ptcl_bk.getPointer(), hard_dump.n_ptcl, hard_dump.ptcl_arti_bk.getPointer(), hard_dump.n_group, hard_dump.time_end, 0);
          HardIntegrator hard_int;
          hard_int.output_filename_prefix = filename;
          auto* ptcl_artificial_ptr =  hard_dump.ptcl_arti_bk.getPointer();
          if (hard_dump.n_arti == 0) ptcl_artificial_ptr = NULL; // if no artificial particle, avoid reading artificial data from last hard_dump
          hard_int.initial(hard_dump.ptcl_bk.getPointer(), hard_dump.n_ptcl, ptcl_artificial_ptr, hard_dump.n_group, hard_dump.n_member_in_group.getPointer(), &hard_manager, hard_dump.time_offset);

          hard_int.integrateToTime(hard_dump.time_end);
          hard_int.driftClusterCMRecordGroupCMDataAndWriteBack(hard_dump.time_end);

      }
      // test stability
      else if (mode==1) {
          typedef H4::ParticleH4<PtclHard> PtclH4;

          SearchGroupCandidate<PtclH4> group_candidate;
          auto* ptcl = hard_dump.ptcl_bk.getPointer();
          PS::S32 n_ptcl = hard_dump.n_ptcl;

          group_candidate.searchAndMerge(ptcl, n_ptcl);

          PS::ReallocatableArray<PtclH4> ptcl_new;
          PS::S32 n_group_in_cluster;

          SystemHard sys;
          sys.manager = &hard_manager;

          PS::ReallocatableArray<COMM::BinaryTree<PtclH4,COMM::Binary>> binary_table;
          PS::ReallocatableArray<SystemHard::GroupIndexInfo> n_member_in_group;
          PS::ReallocatableArray<PS::S32> i_cluster_changeover_update;
          // generate artificial particles, stability test is included
          sys.findGroupsAndCreateArtificialParticlesOneCluster(0, ptcl, n_ptcl, ptcl_new, binary_table, n_group_in_cluster, n_member_in_group, i_cluster_changeover_update, group_candidate, hard_dump.time_end);
      }
  } 

  fclose(fp);

#ifdef STELLAR_EVOLUTION
#ifdef BSE_BASE
  auto& interaction = hard_manager.ar_manager.interaction;
  if (interaction.fout_sse.is_open()) interaction.fout_sse.close();
  if (interaction.fout_bse.is_open()) interaction.fout_bse.close();
#endif
#endif
  return 0;
}
