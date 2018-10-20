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

struct HardPars{
    PS::F64 dt_limit_hard;
    PS::F64 dt_min_hard;
    PS::F64 eta_s;
    PS::F64 time_origin;
    PS::F64 r_bin;
    PS::F64 sdfactor;
    PS::S64 id_offset;
    PS::S32 n_split;

    void read(FILE* fp) {
        size_t rcount = fread(this, sizeof(PS::S32),15,fp);
        if (rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            abort();
        }
    }
};

int main(int argc, char **argv){
  std::string filename="hard_dump";
  int n_opt=0;
  int arg_label;

  while ((arg_label = getopt(argc, argv, "h")) != -1)
    switch (arg_label) {
    case 'h':
        std::cout<<"options:\n"
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

  PS::S32 n_ptcl = ptcl.size();
  std::cout<<"n: "<<n_ptcl<<std::endl;
  //std::cerr<<std::setprecision(20);

  PS::S32 n_artifical, n_group;
  fread(&n_artifical, sizeof(PS::S32),1,fp);
  fread(&n_group, sizeof(PS::S32),1,fp);
  if(n_artifical>0) assert(n_artifical%n_group==0);

  PS::ReallocatableArray<FPSoft> ptcl_artifical;
  ptcl_artifical.resizeNoInitialize(n_artifical);
  
  for(int i=0; i<n_artifical; i++) ptcl_artifical[i].readBinary(fp);

  ARC_int_pars int_pars;
  HardPars hard_pars;

  hard_pars.read(fp);
  int_pars.read(fp);
  std::cout<<"rin="<<int_pars.rin<<std::endl
           <<"rout="<<int_pars.rout<<std::endl
           <<"rbin="<<hard_pars.r_bin<<std::endl;
  
  fclose(fp);

  for(int i=0;i<ptcl.size();i++) {
      std::cerr<<"i="<<i;
      ptcl[i].print(std::cerr);
      std::cerr<<std::endl;
  }

  SearchGroup<PtclHard> group;
  //  group.findGroups(ptcl.getPointer(), ptcl.size(), hard_pars.n_split);
  if (n_ptcl==2) group.searchAndMerge(ptcl.getPointer(), n_ptcl, int_pars.rout);
  else group.searchAndMerge(ptcl.getPointer(), n_ptcl, int_pars.rin);
  std::cout<<"SearchAndMerge\n";
    
  for(int i=0; i<group.getNumOfGroups(); i++) {
      std::cout<<"group["<<i<<"]: ";
      for(int j=0; j<group.getGroupN(i); j++) {
          std::cout<<std::setw(10)<<group.getGroup(i)[j];
      }
      std::cout<<std::endl;
  }

//  std::cout<<"Ptcl List:";
//  for(int i=0; i<group.getPtclN(); i++) {
//      std::cout<<std::setw(10)<<group.getPtclList()[i];
//  }
//  std::cout<<std::endl;

  PS::ReallocatableArray<PtclHard> ptcl_new;
  PS::S32 n_group_in_cluster;

  group.generateList(0, ptcl.getPointer(), n_ptcl, ptcl_new, n_group_in_cluster, hard_pars.r_bin, int_pars.rin, int_pars.rout, time_end, hard_pars.id_offset, hard_pars.n_split);
  
  std::cout<<"GenerateList\n";
  for (int i=0; i<ptcl.size(); i++) {
      ptcl[i].print(std::cout);
      std::cout<<std::endl;
  }
  
  std::cout<<"New ptcl\n";
  for (int i=0; i<ptcl_new.size(); i++) {
      ptcl_new[i].print(std::cout);
      std::cout<<std::endl;
  }
  
  return 0;
}
  
