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
  PS::F64 dE_limit = 1e-4;
  PS::S32 step_limit = 100000;
  while ((arg_label = getopt(argc, argv, "k:e:s:h")) != -1)
    switch (arg_label) {
    case 'k':
        slowdown_factor = atof(optarg);
        n_opt++;
        break;
#ifdef HARD_CHECK_ENERGY
    case 'e':
        dE_limit = atof(optarg);
        n_opt++;
        break;
#endif
#ifdef ARC_SYM
    case 's':
        step_limit = atoi(optarg);
        n_opt++;
        break;
#endif
    case 'h':
        std::cout<<"options:\n"
                 <<"    -k [double]:  change slowdown factor\n"
#ifdef HARD_CHECK_ENERGY
                 <<"    -e [double]:  hard energy limit ("<<dE_limit<<")\n"
#endif
#ifdef ARC_SYM
                 <<"    -s [int]:     ARC SYM step count limit ("<<step_limit<<")\n"
#endif
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
  sys.hard_dE_limit = dE_limit;
#endif
#ifdef ARC_SYM
  // Set step limit for ARC sym
  sys.arc_step_count_limit = step_limit;
#endif

  // set slowdown factor
  if(slowdown_factor>0) sys.set_slowdown_factor(slowdown_factor);
  
  sys.driveForMultiCluster(time_end, ptcl_artifical);

  return 0;
}
