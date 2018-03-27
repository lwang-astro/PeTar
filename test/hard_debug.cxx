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
  while ((arg_label = getopt(argc, argv, "s:h")) != -1)
    switch (arg_label) {
    case 's':
        slowdown_factor = atof(optarg);
        n_opt++;
        break;
    case 'h':
        std::cout<<"options:\n"
                 <<"    -s [double]:  change slowdown factor\n"
                 <<"    -h         :  help\n";
        return 0;
    default:
        std::cerr<<"Unknown argument. check '-h' for help.\n";
        abort();
    }

  if (argc-n_opt*2>1) filename=argv[argc-1];
  
  std::FILE* fp = std::fopen(filename.c_str(),"r");
  if (fp==NULL) {
    std::cerr<<"Error: Filename "<<filename<<" not found\n";
    abort();
  }

  PS::F64 time_end;
  fread(&time_end,sizeof(PS::F64),1,fp);
  std::cout<<"Time_end: "<<time_end<<std::endl;
  PS::S32 first_step_flag;
  fread(&first_step_flag,sizeof(PS::S32),1,fp);
  
  PS::ReallocatableArray<PtclHard> ptcl;
  PtclHardRead(fp,ptcl);
  Ptcl::r_search_min = ptcl[0].r_search;

  std::cout<<"n: "<<ptcl.size();

  SystemHard sys;
  sys.parread(fp);
  fclose(fp);

  if(slowdown_factor>0) sys.set_slowdown_factor(slowdown_factor);
  
  sys.driveForMultiClusterOneDebug(ptcl.getPointer(), ptcl.size(), time_end, first_step_flag);

  return 0;
}
