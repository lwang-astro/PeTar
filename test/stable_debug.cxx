#include <iostream>
#include <cstdio>
#include <iomanip>
#include <string>
#include <unistd.h>
#include <unordered_map>
#include <particle_simulator.hpp>
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

  HardDump hard_dump;
  hard_dump.readOneCluster(filename.c_str());
  std::cerr<<"Time_end: "<<hard_dump.time_end<<std::endl;

  //SystemHard sys;
  //sys.manage = &hard_manager; 

  //PS::ParticleSystem<FPSoft> sys_soft;
  //sys_soft.initialize();
  //sys_soft.createParticle(1000);

  //sys.findGroupsAndCreateArtificalParticlesOMP<PS::ParticleSystem<FPSoft>, FPSoft>(sys_soft, time_end);

  //std::FILE* fp = std::fopen(filename.c_str(),"r");
  //if (fp==NULL) {
  //  std::cerr<<"Error: Filename "<<filename<<" not found\n";
  //  abort();
  //}
  // 
  //PS::F64 time_end;
  //fread(&time_end,sizeof(PS::F64),1,fp);
  //std::cout<<"Time_end: "<<time_end<<std::endl;
  // 
  //PS::ReallocatableArray<PtclHard> ptcl;
  //PtclHardRead(fp,ptcl);
  // 
  //PS::S32 n_ptcl = ptcl.size();
  //std::cout<<"n: "<<n_ptcl<<std::endl;
  ////std::cerr<<std::setprecision(20);
  // 
  //PS::S32 n_artifical, n_group;
  //fread(&n_artifical, sizeof(PS::S32),1,fp);
  //fread(&n_group, sizeof(PS::S32),1,fp);
  //if(n_artifical>0) assert(n_artifical%n_group==0);

  //PS::ReallocatableArray<FPSoft> ptcl_artifical;
  //ptcl_artifical.resizeNoInitialize(n_artifical);
  // 
  //for(int i=0; i<n_artifical; i++) ptcl_artifical[i].readBinary(fp);
  // 
  //ARC_int_pars int_pars;
  //HardPars hard_pars;
  // 
  //hard_pars.read(fp);
  //int_pars.read(fp);
  //std::cout<<"rin="<<int_pars.rin<<std::endl
  //         <<"rout="<<int_pars.rout<<std::endl
  //         <<"rbin="<<hard_pars.r_bin<<std::endl;
  // 
  //fclose(fp);
  // 
  //for(int i=0;i<ptcl.size();i++) {
  //    std::cerr<<"i="<<i;
  //    ptcl[i].print(std::cerr);
  //    std::cerr<<std::endl;
  //}

  typedef H4::ParticleH4<PtclHard> PtclH4;

  SearchGroup<PtclH4> group;
  auto* ptcl = hard_dump.ptcl_bk.getPointer();
  PS::S32 n_ptcl = hard_dump.n_ptcl;
  PS::F64 rout = hard_manager.changeover.getRout();
  PS::F64 rin = hard_manager.changeover.getRin(); 

  if (n_ptcl==2) group.searchAndMerge(ptcl, n_ptcl, rout);
  else group.searchAndMerge(ptcl, n_ptcl, rin);

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
  group.generateList(0, ptcl, n_ptcl, ptcl_new, n_group_in_cluster, hard_manager.h4_manager.r_break_crit, rin, rout, hard_dump.time_end, hard_manager.id_offset, hard_manager.n_split);

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
  
