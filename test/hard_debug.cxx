#include <iostream>
#include <cstdio>
#include <iomanip>
#include <string>
#include <unistd.h>
#include <unordered_map>
#include <particle_simulator.hpp>
#include "Newtonian_acceleration.h"
#include "ptree.h"
#include "kepler.hpp"
#include "cluster_list.hpp"
#include "hard.hpp"
#include "soft.hpp"

int main(int argc, char **argv){
  std::string filename="hard_dump";
  int n_opt=0;
  int arg_label;
  PS::F64 slowdown_factor=0;
  PS::F64 eta=0;
  PS::F64 dE_arc_limit = -1;
  PS::F64 dt_err_pert = -1;
  PS::F64 dt_err_soft = -1;
  PS::F64 dt_arc_min = -1;
  PS::F64 dE_hard_limit = 1e-4;
  PS::S32 step_arc_limit = 100000;
  while ((arg_label = getopt(argc, argv, "k:E:a:p:f:d:e:s:h")) != -1)
    switch (arg_label) {
    case 'k':
        slowdown_factor = atof(optarg);
        n_opt++;
        break;
    case 'E':
        eta = atof(optarg);
        n_opt++;
        break;
    case 'a':
        dE_arc_limit = atof(optarg);
        n_opt++;
        break;
    case 'p':
        dt_err_pert = atof(optarg);
        n_opt++;
        break;
    case 'f':
        dt_err_soft = atof(optarg);
        n_opt++;
        break;
    case 'd':
        dt_arc_min = atof(optarg);
        n_opt++;
        break;
#ifdef HARD_CHECK_ENERGY
    case 'e':
        dE_hard_limit = atof(optarg);
        n_opt++;
        break;
#endif
#ifdef ARC_SYM
    case 's':
        step_arc_limit = atoi(optarg);
        n_opt++;
        break;
#endif
    case 'h':
        std::cout<<"options:\n"
                 <<"    -k [double]:  change slowdown factor\n"
#ifdef HARD_CHECK_ENERGY
                 <<"    -e [double]:  hard energy limit ("<<dE_hard_limit<<")\n"
#endif
#ifdef ARC_SYM
                 <<"    -s [int]:     ARC SYM step count limit ("<<step_arc_limit<<")\n"
#endif
                 <<"    -E [double]:  Eta for hermite \n"
                 <<"    -a [double]:  ARC energy limit \n"
                 <<"    -p [double]:  ARC time error perturbation case \n"
                 <<"    -f [double]:  ARC time error no perturber case \n"
                 <<"    -d [double]:  ARC time step limit \n"
                 <<"    -h         :  help\n";
        return 0;
    default:
        std::cerr<<"Unknown argument. check '-h' for help.\n";
        abort();
    }

  if (argc-n_opt*2>1) filename=argv[argc-1];

  std::cout<<"Reading dump file:"<<filename<<std::endl;

  SystemHard sys;
  PS::ReallocatableArray<FPSoft> ptcl_artifical;

  PS::F64 time_end;
  sys.readOneCluster(filename.c_str(), time_end, ptcl_artifical);

  std::cout<<"Time_end: "<<time_end<<std::endl;

#ifdef HARD_CHECK_ENERGY
    // Set hard energy limit
  sys.hard_dE_limit = dE_hard_limit;
#endif
#ifdef ARC_SYM
  // Set step limit for ARC sym
  sys.arc_step_count_limit = step_arc_limit;
#endif

  // set slowdown factor
  if(slowdown_factor>0) sys.setSlowdownFactor(slowdown_factor);

  // set eta
  if(eta>0) sys.setEta(eta);

  // change ARC parameters
  if(dt_arc_min>0||dE_arc_limit>0||dt_err_pert>0||dt_err_soft>0) 
      sys.setARCParam(dE_arc_limit, dt_err_pert, dt_err_soft, dt_arc_min);
  
  sys.driveForMultiCluster(time_end, ptcl_artifical);

  return 0;
}
