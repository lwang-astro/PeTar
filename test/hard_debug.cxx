#include <iostream>
#include <cstdio>
#include <iomanip>
#include <string>
#include <unistd.h>
#include <unordered_map>
#include <particle_simulator.hpp>
#define HARD_DEBUG_PRINT_FEQ 1024

#include "io.hpp"
#include "hard_assert.hpp"
#include "kepler.hpp"
#include "cluster_list.hpp"
#include "hard.hpp"
#include "soft_ptcl.hpp"


int main(int argc, char **argv){
  int n_opt=0;
  int arg_label;
  PS::F64 slowdown_factor=0;
  PS::F64 eta=0;
  PS::F64 e_err_ar = -1;
  PS::F64 e_err_hard = 1e-4;
  PS::S32 dt_min_power = -1;
  PS::F64 dt_max = -1;
  PS::S32 step_arc_limit = 100000;

  while ((arg_label = getopt(argc, argv, "k:E:a:D:d:e:s:h")) != -1)
    switch (arg_label) {
    case 'k':
        slowdown_factor = atof(optarg);
        n_opt+=2;
        break;
    case 'E':
        eta = atof(optarg);
        n_opt+=2;
        break;
    case 'a':
        e_err_ar = atof(optarg);
        n_opt+=2;
        break;
    case 'D':
        dt_max = atof(optarg);
        n_opt+=2;
        break;
    case 'd':
        dt_min_power = atoi(optarg);
        n_opt+=2;
        break;
#ifdef HARD_CHECK_ENERGY
    case 'e':
        e_err_hard = atof(optarg);
        n_opt+=2;
        break;
#endif
#ifdef AR_SYM
    case 's':
        step_arc_limit = atoi(optarg);
        n_opt+=2;
        break;
#endif
    case 'h':
        std::cout<<"hard_debug.out [options] [hard_manager (defaulted: input.par.hard)] [cluster_data] (defaulted: hard_dump)\n"
                 <<"options:\n"
                 <<"    -k [double]:  change slowdown factor\n"
#ifdef HARD_CHECK_ENERGY
                 <<"    -e [double]:  hard energy limit ("<<e_err_hard<<")\n"
#endif
#ifdef AR_SYM
                 <<"    -s [int]:     AR step count limit ("<<step_arc_limit<<")\n"
#endif
                 <<"    -E [double]:  Eta for hermite \n"
                 <<"    -a [double]:  AR energy limit \n"
                 <<"    -D [double]:  hard time step max \n"
                 <<"    -d [int]:     hard time step min power \n"
                 <<"    -h         :  help\n";
        return 0;
    default:
        std::cerr<<"Unknown argument. check '-h' for help.\n";
        abort();
    }

  std::string filename="hard_dump";
  std::string fhardpar="input.par.hard";
  if (argc-n_opt>1) {
      filename=argv[argc-1];
      if (argc-n_opt>2) 
          fhardpar=argv[argc-2];
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

#ifdef HARD_CHECK_ENERGY
  // Set hard energy limit
  if (e_err_hard>0 )
      hard_manager.energy_error_max = e_err_hard;
#endif
#ifdef AR_SYM
  // Set step limit for ARC sym
  if (step_arc_limit>0) 
      hard_manager.ar_manager.step_count_max = step_arc_limit;
#endif
  // set slowdown factor
  if(slowdown_factor>0)    
      hard_manager.ar_manager.slowdown_pert_ratio_ref = slowdown_factor;

  // set eta
  if(eta>0) {
      hard_manager.h4_manager.step.eta_4th = eta;
      hard_manager.h4_manager.step.eta_2nd = 0.1*eta;
  }
  // time step
  if(dt_min_power>0&&dt_max>0) 
      hard_manager.setDtRange(dt_max, dt_min_power);

  if(e_err_ar>0) 
      hard_manager.ar_manager.energy_error_relative_max = e_err_ar;

  hard_manager.checkParams();
  hard_manager.print(std::cerr);

  HardDump hard_dump;
  hard_dump.readOneCluster(filename.c_str());
  std::cerr<<"Dt: "<<hard_dump.time_end<<std::endl;

  SystemHard sys;
  sys.manager = &hard_manager;

  // change ARC parameters
  sys.driveForMultiClusterImpl(hard_dump.ptcl_bk.getPointer(), hard_dump.n_ptcl, hard_dump.ptcl_arti_bk.getPointer(), hard_dump.n_group, hard_dump.time_end, 0);

  return 0;
}
