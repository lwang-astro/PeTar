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
#include "hard_ptcl.hpp"
#include "cluster_list.hpp"
#include "hard.hpp"


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

  hard_manager.checkParams();
  hard_manager.print(std::cerr);

  HardDump hard_dump;
  hard_dump.readOneCluster(filename.c_str());
  std::cerr<<"Time_end: "<<hard_dump.time_end<<std::endl;

  typedef H4::ParticleH4<PtclHard> PtclH4;

  SearchGroup<PtclH4> group;
  auto* ptcl = hard_dump.ptcl_bk.getPointer();
  PS::S32 n_ptcl = hard_dump.n_ptcl;

  group.searchAndMerge(ptcl, n_ptcl);

  // generate artifical particles,
  //std::cout<<"SearchAndMerge\n";
  //  
  //for(int i=0; i<group.getNumOfGroups(); i++) {
  //    std::cout<<"group["<<i<<"]: ";
  //    for(int j=0; j<group.getGroupN(i); j++) {
  //        std::cout<<std::setw(10)<<group.getGroup(i)[j];
  //    }
  //    std::cout<<std::endl;
  //}
  // 
////  std::cout<<"Ptcl List:";
////  for(int i=0; i<group.getPtclN(); i++) {
////      std::cout<<std::setw(10)<<group.getPtclList()[i];
////  }
////  std::cout<<std::endl;

  PS::ReallocatableArray<PtclH4> ptcl_new;
  PS::S32 n_group_in_cluster;
  group.generateList(0, ptcl, n_ptcl, ptcl_new, n_group_in_cluster, hard_manager.r_tidal_tensor, hard_manager.r_in_base, hard_manager.r_out_base, hard_dump.time_end, hard_manager.id_offset, hard_manager.n_split);

  //std::cout<<"GenerateList\n";
  //for (int i=0; i<ptcl.size(); i++) {
  //    ptcl[i].print(std::cout);
  //    std::cout<<std::endl;
  //}
  // 
  //std::cout<<"New ptcl\n";
  //for (int i=0; i<ptcl_new.size(); i++) {
  //    ptcl_new[i].print(std::cout);
  //    std::cout<<std::endl;
  //}
  
  return 0;
}
  
