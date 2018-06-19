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
  std::cerr<<std::setprecision(20);
  PS::F64 r_search_max=0.0;
  PS::S32 i_r_search_max=-1;
  std::cerr<<"member :\n";
  for(int i=0; i<ptcl.size(); i++) {
      if(ptcl[i].status<0) {
          ptcl[i].print(std::cerr);
          std::cerr<<std::endl;
          //std::cerr<<ptcl[i].mass_bk<<" "<<ptcl[i].pos<<" "<<ptcl[i].vel<<std::endl;
      }
  }
  std::cerr<<"single:\n";
  for(int i=0; i<ptcl.size(); i++) {
      if(ptcl[i].status==0) {
          ptcl[i].print(std::cerr);
          std::cerr<<std::endl;
          //std::cerr<<ptcl[i].mass<<" "<<ptcl[i].pos<<" "<<ptcl[i].vel<<std::endl;
      }
      if(r_search_max<ptcl[i].r_search) {
          r_search_max =ptcl[i].r_search;
          i_r_search_max = i;
      }
  }

  std::cout<<"n_group: "<<n_group<<std::endl;
  
  std::cout<<"R_search_max = "<<r_search_max;
  ptcl[i_r_search_max].print(std::cout);
  std::cout<<std::endl;

  SystemHard sys;
  sys.parread(fp);
  fclose(fp);

  if(slowdown_factor>0) sys.set_slowdown_factor(slowdown_factor);
  
  sys.driveForMultiClusterOneDebug(ptcl.getPointer(), ptcl.size(), ptcl_artifical.getPointer(), n_group, time_end);

  return 0;
}
