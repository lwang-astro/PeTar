#include <iostream>
#include <fstream>
#include <iomanip>
#include <unordered_map>
#include <particle_simulator.hpp>
#include "Newtonian_acceleration.h"
#include "ptree.h"
#include "kepler.hpp"
#include "hard.hpp"
#include "soft.hpp"


#ifndef NAN_CHECK
#define NAN_CHECK(val) assert((val) == (val));
#endif


//struct params{
//    double rin,rsearch,dt_tree;
//    int n_split;
//};

PtclHard pshift(const PtclHard& a, const PtclHard& ref) {
    PtclHard p(a);
    p.pos += ref.pos;
    p.vel += ref.vel;
    return p;
}

PtclHard kepler_print(const std::size_t id, const std::size_t ib, PtclHard* c[2], PS::F64& rsearch){
    PS::F64vec x[2],v[2];
    double m[2];
    for (std::size_t i=0; i<2; i++) {
      x[i]=c[i]->pos;
      v[i]=c[i]->vel;
      m[i]=c[i]->mass;
    }
        
    double ax,per,ecc,angle[3],true_anomaly,ecc_anomaly,mean_anomaly; 
    double dx[3] = {x[1][0]-x[0][0], x[1][1]-x[0][1], x[1][2]-x[0][2]};
    double dv[3] = {v[1][0]-v[0][0], v[1][1]-v[0][1], v[1][2]-v[0][2]};
    double mt = m[0]+m[1];
    
    NTA::calc_kepler_orbit_par(ax,per,ecc,angle,true_anomaly,ecc_anomaly,mean_anomaly,mt,dx,dv);
    std::cout<<std::setw(20)<<id
             <<std::setw(20)<<ib
             <<std::setw(20)<<m[0]
             <<std::setw(20)<<m[1]
             <<std::setw(20)<<ax
             <<std::setw(20)<<ecc
             <<std::setw(20)<<per
             <<std::setw(20)<<angle[0]
             <<std::setw(20)<<angle[1]
             <<std::setw(20)<<angle[2]
             <<std::setw(20)<<ecc_anomaly
             <<std::setw(20)<<true_anomaly
             <<std::setw(20)<<mean_anomaly
             <<std::endl;

    PS::F64vec xcm((x[0][0]*m[0]+x[1][0]*m[1])/mt, 
                   (x[0][1]*m[0]+x[1][1]*m[1])/mt, 
                   (x[0][2]*m[0]+x[1][2]*m[1])/mt);

    PS::F64vec vcm((v[0][0]*m[0]+v[1][0]*m[1])/mt, 
                   (v[0][1]*m[0]+v[1][1]*m[1])/mt, 
                   (v[0][2]*m[0]+v[1][2]*m[1])/mt);

    return PtclHard(Ptcl(ParticleBase(mt, xcm, vcm), rsearch, 0, ib, 0));
}

void print_p(PtclHard* p, const int n) {
    std::cout<<std::setw(20)<<"mass"
             <<std::setw(20)<<"x1"
             <<std::setw(20)<<"x2"
             <<std::setw(20)<<"x3"
             <<std::setw(20)<<"v1"
             <<std::setw(20)<<"v2"
             <<std::setw(20)<<"v3"
             <<std::setw(20)<<"rsearch"
             <<std::setw(20)<<"mass_bk"
             <<std::setw(20)<<"status"
             <<std::setw(20)<<"id"
             <<std::setw(20)<<"id_cluster"
             <<std::setw(20)<<"adr"
             <<std::endl;
    for (int i=0; i<n; i++) {
        std::cout<<std::setw(20)<<p[i].mass
                 <<std::setw(20)<<p[i].pos[0]
                 <<std::setw(20)<<p[i].pos[1]
                 <<std::setw(20)<<p[i].pos[2]
                 <<std::setw(20)<<p[i].vel[0]
                 <<std::setw(20)<<p[i].vel[1]
                 <<std::setw(20)<<p[i].vel[2]
                 <<std::setw(20)<<p[i].r_search
                 <<std::setw(20)<<p[i].mass_bk
                 <<std::setw(20)<<p[i].status
                 <<std::setw(20)<<p[i].id
                 <<std::setw(20)<<p[i].id_cluster
                 <<std::setw(20)<<p[i].adr_org
                 <<std::endl;
    }
}

void print_p(FPSoft* p, const int n, const int adr_ref) {
    std::cout<<std::setw(20)<<"mass"
             <<std::setw(20)<<"x1"
             <<std::setw(20)<<"x2"
             <<std::setw(20)<<"x3"
             <<std::setw(20)<<"v1"
             <<std::setw(20)<<"v2"
             <<std::setw(20)<<"v3"
             <<std::setw(20)<<"rsearch"
             <<std::setw(20)<<"mass_bk"
             <<std::setw(20)<<"status"
             <<std::setw(20)<<"id"
             <<std::setw(20)<<"adr"
             <<std::endl;
    for (int i=0; i<n; i++) {
        std::cout<<std::setw(20)<<p[i].mass
                 <<std::setw(20)<<p[i].pos[0]
                 <<std::setw(20)<<p[i].pos[1]
                 <<std::setw(20)<<p[i].pos[2]
                 <<std::setw(20)<<p[i].vel[0]
                 <<std::setw(20)<<p[i].vel[1]
                 <<std::setw(20)<<p[i].vel[2]
                 <<std::setw(20)<<p[i].r_search
                 <<std::setw(20)<<p[i].mass_bk
                 <<std::setw(20)<<p[i].status
                 <<std::setw(20)<<p[i].id
                 <<std::setw(20)<<i+adr_ref
                 <<std::endl;
    }
}

int main(int argc, char** argv)
{
  // data file name
  char* filename = argv[argc-1];
  // open data file
  std::fstream fs;
  fs.open(filename,std::fstream::in);
  if(!fs.is_open()) {
    std::cerr<<"Error: Filename "<<filename<<" not found\n";
    abort();
  }

  PS::ReallocatableArray<PtclHard> p;

  PS::S32 n_cluster, n_split;
  PS::F64 rbin, rin, rout, rsearch, dt_tree;
  fs>>n_cluster>>rbin>>rin>>rout>>rsearch>>dt_tree>>n_split;;
  Ptcl::r_search_min = rsearch;

  PS::S32 n_offset[n_cluster+1];
  PS::ReallocatableArray<PS::S32> n_ptcl_cluster;
  n_ptcl_cluster.resizeNoInitialize(n_cluster);
  PS::ReallocatableArray<PS::S32> n_group;
  
  for (int i=0; i<n_cluster;i++) {

      PS::S32 N;
      fs>>N;
      n_offset[i]=p.size();
      n_ptcl_cluster[i]=N;
      p.increaseSize(N);
  
      ptree<PtclHard,PS::F64> plist;
      int idc=1;
      for(int i=0;i<N-1;i++) {
          int id,ib;
          double m1,m2,ax,ecc,angle[3],ecc_anomaly;
          fs>>id>>ib>>m1>>m2>>ax>>ecc>>angle[0]>>angle[1]>>angle[2]>>ecc_anomaly;
          if (fs.eof()) {
              std::cerr<<"Error: data file reach end when reading pairs (current loaded pair number is "<<i<<"; required pair number "<<N-1<<std::endl;
              abort();
          }
    
          double x1[3],x2[3],v1[3],v2[3];
          NTA::kepler_orbit_generator(x1,x2,v1,v2,m1,m2,ax,ecc,angle,ecc_anomaly);
    
          PS::F64vec xx1(x1[0],x1[1],x1[2]);
          PS::F64vec xx2(x2[0],x2[1],x2[2]);
          PS::F64vec vv1(v1[0],v1[1],v1[2]);
          PS::F64vec vv2(v2[0],v2[1],v2[2]);    

          PtclHard a(Ptcl(ParticleBase(m1,xx1,vv1),rsearch,0,idc++,0));
          PtclHard b(Ptcl(ParticleBase(m2,xx2,vv2),rsearch,0,idc++,0));

          bool flag=plist.link(id,ib,a,b,pshift);
          if (!flag) {
              std::cerr<<"Error: particle id "<<id<<", ib "<<ib<<" are inconsistent with global particle tree structure, cannot created pairs!\n";
              abort();
          }
      }

      int count=plist.collect_and_store(p.getPointer(n_offset[i]),N);
      if (count<0) {
          std::cerr<<"Error: particle number mismatched particle tree!\n";
          abort();
      }
      std::cout<<std::setprecision(13);
      //  plist.kepler_print(0,0,18,10);
      // header
      std::cout<<"        " 
               <<std::setw(20)<<"id          " 
               <<std::setw(20)<<"ib          " 
               <<std::setw(20)<<"m[0]        " 
               <<std::setw(20)<<"m[1]        " 
               <<std::setw(20)<<"ax          " 
               <<std::setw(20)<<"ecc         " 
               <<std::setw(20)<<"per         " 
               <<std::setw(20)<<"angle[0]    " 
               <<std::setw(20)<<"angle[1]    " 
               <<std::setw(20)<<"angle[2]    " 
               <<std::setw(20)<<"ecc_anomaly " 
               <<std::setw(20)<<"true_anomaly" 
               <<std::setw(20)<<"mean_anomaly" 
               <<std::endl;                 
      
      std::cout<<std::setw(20)<<"Group"<<std::setw(20)<<i<<std::endl;
      plist.pair_process(0,0,kepler_print,rsearch);

      for (int j=0; j<n_ptcl_cluster[i]; j++) p[j+n_offset[i]].adr_org = j+n_offset[i];
      print_p(p.getPointer(),N);
  }

  PS::ParticleSystem<FPSoft> sys;
  for (int i=0; i<p.size(); i++) sys.addOneParticle(FPSoft(p[i],0,i));
  const PS::S32 n_sys = sys.getNumberOfParticleLocal();

  PS::ReallocatableArray<PS::S32> ptcl_list;
  ptcl_list.resizeNoInitialize(p.size());
  for (int i=0; i<n_cluster; i++) {
      std::cout<<"index list group "<<i<<std::endl;
      for (int j=0; j<n_ptcl_cluster[i]; j++) {
          ptcl_list[n_offset[i]+j] = n_offset[i]+n_ptcl_cluster[i]-j-1;
          std::cout<<std::setw(10)<<ptcl_list[n_offset[i]+j];
      }
      std::cout<<std::endl;
  }
  PS::S64 id_offset=p.size();
  
  SystemHard sys_hard;
  sys_hard.setParam(rbin, rout, rin, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, id_offset, n_split);

  sys_hard.setPtclForIsolatedMultiCluster(sys, ptcl_list, n_ptcl_cluster);
  sys_hard.findGroupsAndCreateArtificalParticlesOMP<PS::ParticleSystem<FPSoft>, FPSoft>(sys, dt_tree);

  for(int i=0; i<n_cluster; i++) {
      std::cout<<"new index list group "<<i<<std::endl;
      for (int j=0; j<n_ptcl_cluster[i]; j++) {
          std::cout<<std::setw(10)<<sys_hard.getPtcl()[n_offset[i]+j].adr_org;
      }
      std::cout<<"\nPtcl:\n";
      print_p(&sys[n_offset[i]], n_ptcl_cluster[i], n_offset[i]);
      std::cout<<std::endl;
  }
  std::cout<<"New artifical particles: N="<<sys.getNumberOfParticleLocal()-n_sys<<"\n";
  print_p(&sys[n_sys], sys.getNumberOfParticleLocal()-n_sys, n_sys);

  std::cout<<"N group, offset and first adr of artifical:\n";
  for(int i=0; i<n_cluster; i++) {
      std::cout<<std::setw(10)<<*sys_hard.getGroupNList(i);
      std::cout<<" offset: "<<std::setw(10)<<*sys_hard.getGroupNOffset(i+1)<<" adr: ";
      for(int j=*sys_hard.getGroupNOffset(i); j<*sys_hard.getGroupNOffset(i+1); j++) 
          std::cout<<std::setw(10)<<(*sys_hard.getAdrPtclArtFirstList(j));
      std::cout<<std::endl;
  }
  
  
  //group.findGroups(p.getPointer(),p.size(), par.n_split);
  // 
  //PS::ReallocatableArray<PtclHard> pn;
  //PS::ReallocatableArray<PtclHard> pg;
  //PS::ReallocatableArray<PtclHard> pp;
  //std::cout<<"find groups\n";
  //for(int i=0;i<group.getPtclN();i++) pn.push_back(p[group.getPtclList()[i]]);
  //print_p(pn.getPointer(),pn.size());
  // 
  //for(int i=0;i<group.getNumOfGroups();i++)  {
  //    pg.clearSize();
  //    std::cout<<"Group["<<i<<"]: ";
  //    for(int j=0;j<group.getGroupN(i);j++) 
  //        std::cout<<std::setw(10)<<group.getGroup(i)[j];
  //    std::cout<<std::endl;
  //    for(int j=0;j<group.getGroupN(i);j++) 
  //        pg.push_back(p[group.getGroup(i)[j]]);
  //    print_p(pg.getPointer(),pg.size());
  //    //for(int j=0;j<group.getGroupN(i);j++) {
  //    //    std::cout<<"Member "<<j<<": ";
  //    //    for(int k=0;k<8;k++) std::cout<<std::setw(10)<<group.getGroupPertAdr(i,j,k);
  //    //    std::cout<<std::endl;
  //    //}
  //    pp.clearSize();
  //    for(int j=0;j<16;j++) {
  //        pp.push_back(p[group.getGroupPertList(i,par.n_split)[j]]);
  //    }
  //    std::cout<<"Pert force:"<<std::endl;
  //    print_p(pp.getPointer(),pp.size());
  // 
  //}
  //group.resolveGroups();
  //std::cout<<"resolvegroups\n";
  //pn.clearSize();
  //for(int i=0;i<group.getPtclN();i++) pn.push_back(p[group.getPtclList()[i]]);
  //print_p(pn.getPointer(),pn.size());
  // 
  //// second try
  //group.searchAndMerge(p.getPointer(), par.rin);
  //std::cout<<"SearchAndMerge 2nd\n";
  //print_p(p.getPointer(),N);
  // 
  //for(int i=0; i<group.getNumOfGroups(); i++) {
  //    std::cout<<"group["<<i<<"]: ";
  //    for(int j=0; j<group.getGroupN(i); j++) {
  //        std::cout<<std::setw(10)<<group.getGroup(i)[j];
  //    }
  //    std::cout<<std::endl;
  //}
  // 
  //std::cout<<"Ptcl List 2nd:";
  //for(int i=0; i<group.getPtclN(); i++) {
  //    std::cout<<std::setw(10)<<group.getPtclList()[i];
  //}
  //std::cout<<std::endl;
  // 
  ////PS::ReallocatableArray<PS::S32> adr_group_glb;
  ////PS::ReallocatableArray<std::vector<Ptcl>> group_ptcl_glb; 
  ////PS::ReallocatableArray<PS::S32> group_ptcl_glb_empty_list;
  //ptcl_new.clearSize();
  // 
  ////group.generateList(p, N, adr_group_glb, group_ptcl_glb, group_ptcl_glb_empty_list);
  //group.generateList(p.getPointer(), ptcl_new, par.rin, par.rin, par.rin, par.dt_tree,N,par.n_split);
  //std::cout<<"GenerateList 2nd\n";
  //print_p(p.getPointer(),p.size());
  // 
  ////std::cout<<"adr_group_glb: ";
  ////for(int i=0; i<adr_group_glb.size(); i++) std::cout<<std::setw(6)<<adr_group_glb[i];
  ////std::cout<<std::endl;
  // 
  //std::cout<<"new ptcl 2nd: "<<ptcl_new.size()<<"\n ";
  //print_p(ptcl_new.getPointer(),ptcl_new.size());
  ////for(int i=0; i<ptcl_new.size(); i++) {
  ////    print_p(group_ptcl_glb[i].data(),group_ptcl_glb[i].size());
  ////}
  ////p.reserveEmptyAreaAtLeast(ptcl_new.size());
  ////for(int i=0;i<ptcl_new.size();i++) p.pushBackNoCheck(ptcl_new[i]);
  // 
  //group.findGroups(p.getPointer(),p.size(),par.n_split);
  // 
  //pn.clearSize();
  //pg.clearSize();
  //pp.clearSize();
  //std::cout<<"find groups 2nd\n";
  //for(int i=0;i<group.getPtclN();i++) pn.push_back(p[group.getPtclList()[i]]);
  //print_p(pn.getPointer(),pn.size());
  // 
  //for(int i=0;i<group.getNumOfGroups();i++)  {
  //    pg.clearSize();
  //    std::cout<<"Group["<<i<<"]: ";
  //    for(int j=0;j<group.getGroupN(i);j++) 
  //        std::cout<<std::setw(10)<<group.getGroup(i)[j];
  //    std::cout<<std::endl;
  //    for(int j=0;j<group.getGroupN(i);j++) 
  //        pg.push_back(p[group.getGroup(i)[j]]);
  //    print_p(pg.getPointer(),pg.size());
  //    //for(int j=0;j<group.getGroupN(i);j++) {
  //    //    std::cout<<"Member "<<j<<": ";
  //    //    for(int k=0;k<8;k++) std::cout<<std::setw(10)<<group.getGroupPertAdr(i,j,k);
  //    //    std::cout<<std::endl;
  //    //}
  //    pp.clearSize();
  //    for(int j=0;j<16;j++) {
  //        pp.push_back(p[group.getGroupPertList(i,par.n_split)[j]]);
  //    }
  //    std::cout<<"Pert force:"<<std::endl;
  //    print_p(pp.getPointer(),pp.size());
  // 
  //}
  //group.resolveGroups();
  //std::cout<<"resolvegroups 2nd\n";
  //pn.clearSize();
  //for(int i=0;i<group.getPtclN();i++) pn.push_back(p[group.getPtclList()[i]]);
  //print_p(pn.getPointer(),pn.size());
      
  
  return 0;
}
