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
  
  std::FILE* fp = std::fopen(filename.c_str(),"r");
  if (fp==NULL) {
    std::cerr<<"Error: Filename "<<filename<<" not found\n";
    abort();
  }
  std::cout<<"Reading dump file:"<<filename<<std::endl;

  PS::F64 time_end;
  fread(&time_end,sizeof(PS::F64),1,fp);
  std::cout<<"Time_end: "<<time_end<<std::endl;
  
  PS::ReallocatableArray<PtclHard> ptcl;
  PtclHardRead(fp,ptcl);

  PS::S32 n_artifical, n_group;
  fread(&n_artifical, sizeof(PS::S32),1,fp);
  fread(&n_group, sizeof(PS::S32),1,fp);
  if(n_artifical>0) assert(n_artifical%n_group==0);

  PS::ReallocatableArray<FPSoft> ptcl_artifical;
  ptcl_artifical.resizeNoInitialize(n_artifical);
  
  for(int i=0; i<n_artifical; i++) ptcl_artifical[i].readBinary(fp);

  std::cout<<"n: "<<ptcl.size()<<std::endl;
  std::cout<<std::setprecision(14);
  std::cerr<<std::setprecision(12);
  PS::F64 r_search_max=0.0;
  PS::S32 i_r_search_max=-1;
  if(n_group>0) std::cerr<<"member :\n";
  for(int i=0; i<ptcl.size(); i++) {
      if(ptcl[i].status<0) {
          ptcl[i].print(std::cerr);
          std::cerr<<std::endl;
          //std::cerr<<ptcl[i].mass_bk<<" "<<ptcl[i].pos<<" "<<ptcl[i].vel<<std::endl;
      }
  }
  std::cerr<<"single:\n";
  PS::S32 n_single=0;
  PS::ReallocatableArray<PS::S32> single_adr;
  for(int i=0; i<ptcl.size(); i++) {
      if(ptcl[i].status==0) {
          ptcl[i].print(std::cerr);
          std::cerr<<std::endl;
          //std::cerr<<ptcl[i].mass<<" "<<ptcl[i].pos<<" "<<ptcl[i].vel<<std::endl;
          n_single++;
          single_adr.push_back(i);
      }
      if(r_search_max<ptcl[i].r_search) {
          r_search_max =ptcl[i].r_search;
          i_r_search_max = i;
      }
  }

  std::cout<<"n_group: "<<n_group<<std::endl;

  if(n_single==2) {
      std::cout<<"Kepler parameters for two single:\n";
      Binary bin;
      PosVel2OrbParam(bin, ptcl[single_adr[0]], ptcl[single_adr[1]]);
      bin.print(std::cout,14,true);
  }
  
  std::cout<<"R_search_max = "<<r_search_max;
  ptcl[i_r_search_max].print(std::cout);
  std::cout<<std::endl;

  SystemHard sys;
  sys.parread(fp);

  fclose(fp);

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
  
  sys.driveForMultiClusterOneDebug(ptcl.getPointer(), ptcl.size(), ptcl_artifical.getPointer(), n_group, 1000.0, time_end);

  return 0;
}
