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

int main(int argc, char **argv){
  int arg_label;
  int mode=0; // 0: integrate to time; 1: times stability
  PS::F64 slowdown_factor=0;
  PS::F64 eta_4th=0;
  PS::F64 eta_2nd=0;
  PS::F64 e_err_ar = -1;
  PS::F64 e_err_hard = 1e-4;
  PS::S32 dt_min_power = -1;
  PS::F64 dt_max = -1;
  PS::F64 ds_scale = -1.0;
  PS::S32 par_version = -1;
  PS::S32 step_arc_limit = 100000;
  std::string filename="hard_dump";
  std::string fhardpar="input.par.hard";
#ifdef BSE_BASE
  int idum=0;
  std::string bse_name = BSEManager::getBSEName();
  std::string fsse_suffix = BSEManager::getSSEOutputFilenameSuffix();
  std::string fbse_suffix = BSEManager::getBSEOutputFilenameSuffix();

  std::string fbsepar = "input.par" + fbse_suffix;
  std::string fbserandpar = "bse.rand.par";
#endif
#ifdef STELLAR_EVOLUTION
  int stellar_evolution_option = -1;
#endif
#ifdef SOFT_PERT
  bool soft_pert_flag=true;
#endif

  while ((arg_label = getopt(argc, argv, "k:E:A:a:D:d:e:s:c:m:b:B:p:I:v:i:Sh")) != -1)
    switch (arg_label) {
    case 'k':
        slowdown_factor = atof(optarg);
        break;
    case 'E':
        eta_4th = atof(optarg);
        break;
    case 'A':
        eta_2nd = atof(optarg);
        break;
    case 'a':
        e_err_ar = atof(optarg);
        break;
    case 'D':
        dt_max = atof(optarg);
        break;
    case 'd':
        dt_min_power = atoi(optarg);
        break;
#ifdef HARD_CHECK_ENERGY
    case 'e':
        e_err_hard = atof(optarg);
        break;
#endif
    case 's':
        step_arc_limit = atoi(optarg);
        break;
    case 'c':
        ds_scale = atof(optarg);
        break;
    case 'm':
        mode = atoi(optarg);
        break;
    case 'p':
        fhardpar = optarg;
        break;
#ifdef SOFT_PERT
    case 'S':
        soft_pert_flag=false;
        break;
#endif
#ifdef STELLAR_EVOLUTION
    case 'I':
        stellar_evolution_option = atoi(optarg);
        break;
#endif
#ifdef BSE_BASE
    case 'i':
        idum = atoi(optarg);
        break;
    case 'b':
        fbsepar = optarg;
        break;
    case 'B':
        fbserandpar = optarg;
        break;
    case 'v':
        par_version = atoi(optarg);
        break;
#endif
    case 'h':
        std::cout<<"petar.hard.debug [options] [hard_manager (defaulted: input.par.hard)] [cluster_data] (defaulted: hard_dump)\n"
                 <<"options:\n"
                 <<"    -k [double]:  change slowdown factor reference\n"
#ifdef HARD_CHECK_ENERGY
                 <<"    -e [double]:  hard energy limit\n"
#endif
                 <<"    -s [int]:     AR step count limit\n"
                 <<"    -c [double]:  AR step scaling factor\n"
                 <<"    -E [double]:  Eta 4th for hermite \n"
                 <<"    -A [double]:  Eta 2nd for hermite \n"
                 <<"    -a [double]:  AR energy limit \n"
                 <<"    -D [double]:  hard time step max (should use together with -d)\n"
                 <<"    -d [int]:     hard time step min power (should use together with -D)\n"
                 <<"    -m [int]:     running mode: 0: evolve system to time_end; 1: stability check: "<<mode<<std::endl
                 <<"    -p [string]:  hard parameter file name: "<<fhardpar<<std::endl
#ifdef STELLAR_EVOLUTION
                 <<"    -I [int]:     Stellar evolution option: \n"
#endif
#ifdef BSE_BASE
                 <<"    -i [int]      random seed to generate kick velocity\n"
                 <<"    -B [string]:  read bse random parameter dump file with filename: "<<fbserandpar<<"\n"
                 <<"    -b [string]:  bse parameter file name: "<<fbsepar<<std::endl
#endif
#ifdef SOFT_PERT
                 <<"    -S:           Suppress soft perturbation (tidal tensor)\n"
#endif
                 <<"    -v [int]:     version of hard parameters: 0: default, 1: mssing ds_scale in ar_manager: 0\n"
                 <<"    -h:           help\n";
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
  if (par_version<0) par_version = 0;
  hard_manager.readBinary(fpar_in, par_version);
  fclose(fpar_in);

#ifdef STELLAR_EVOLUTION
  if (stellar_evolution_option>=0) {
      hard_manager.ar_manager.interaction.stellar_evolution_option = stellar_evolution_option;
      if (stellar_evolution_option==0) 
          hard_manager.ar_manager.interaction.stellar_evolution_write_flag = false;
  }
#ifdef BSE_BASE
  IOParamsBSE bse_io;
  std::cerr<<bse_name<<" parameter file:"<<fbsepar<<std::endl;
  if( (fpar_in = fopen(fbsepar.c_str(),"r")) == NULL) {
      fprintf(stderr,"Error: Cannot open file %s.\n", fbsepar.c_str());
      abort();
  }
  bse_io.input_par_store.readAscii(fpar_in);
  fclose(fpar_in);
  if (idum!=0) bse_io.idum.value = idum;
  hard_manager.ar_manager.interaction.bse_manager.initial(bse_io);

  std::cerr<<"Check "<<bse_name<<" rand parameter file: "<<fbserandpar<<std::endl;
  hard_manager.ar_manager.interaction.bse_manager.readRandConstant(fbserandpar.c_str());
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

  HardDump hard_dump;
  hard_dump.readOneCluster(filename.c_str());
  std::cerr<<"Time end: "<<hard_dump.time_end<<std::endl;

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
      hard_int.initial(hard_dump.ptcl_bk.getPointer(), hard_dump.n_ptcl, hard_dump.ptcl_arti_bk.getPointer(), hard_dump.n_group, hard_dump.n_member_in_group.getPointer(), &hard_manager, hard_dump.time_offset);

      auto& interrupt_binary = hard_int.integrateToTime(hard_dump.time_end);
      if (interrupt_binary.status!=AR::InterruptStatus::none) {
          hard_int.printInterruptBinaryInfo(std::cerr);
      }
      else {
          hard_int.driftClusterCMRecordGroupCMDataAndWriteBack(hard_dump.time_end);
      }

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

#ifdef STELLAR_EVOLUTION
#ifdef BSE_BASE
  auto& interaction = hard_manager.ar_manager.interaction;
  if (interaction.fout_sse.is_open()) interaction.fout_sse.close();
  if (interaction.fout_bse.is_open()) interaction.fout_bse.close();
#endif
#endif
  return 0;
}
