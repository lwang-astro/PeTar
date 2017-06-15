#include <iostream>
#include <fstream>
#include <iomanip>
#include <unordered_map>
#include <particle_simulator.hpp>
#include "Newtonian_acceleration.h"
#include "ptree.h"
#include "kepler.hpp"
#include "rsearch.hpp"
#include "cluster_list.hpp"
#include "hard.hpp"


#ifndef NAN_CHECK
#define NAN_CHECK(val) assert((val) == (val));
#endif


struct params{
    double rin,rout,gmin;
};

PtclHard pshift(const PtclHard& a, const PtclHard& ref) {
    PtclHard p(a);
    p.pos += ref.pos;
    p.vel += ref.vel;
    return p;
}

PtclHard kepler_print(const std::size_t id, const std::size_t ib, PtclHard* c[2], params& ppars){
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
    std::cout<<std::setw(12)<<id
             <<std::setw(12)<<ib
             <<std::setw(12)<<m[0]
             <<std::setw(12)<<m[1]
             <<std::setw(12)<<ax
             <<std::setw(12)<<ecc
             <<std::setw(12)<<per
             <<std::setw(12)<<angle[0]
             <<std::setw(12)<<angle[1]
             <<std::setw(12)<<angle[2]
             <<std::setw(12)<<ecc_anomaly
             <<std::setw(12)<<true_anomaly
             <<std::setw(12)<<mean_anomaly
             <<std::endl;
    if(id==0&&ib==0) 
        std::cout<<"      " 
                 <<std::setw(12)<<"id          " 
                 <<std::setw(12)<<"ib          " 
                 <<std::setw(12)<<"m[0]        " 
                 <<std::setw(12)<<"m[1]        " 
                 <<std::setw(12)<<"ax          " 
                 <<std::setw(12)<<"ecc         " 
                 <<std::setw(12)<<"per         " 
                 <<std::setw(12)<<"angle[0]    " 
                 <<std::setw(12)<<"angle[1]    " 
                 <<std::setw(12)<<"angle[2]    " 
                 <<std::setw(12)<<"ecc_anomaly " 
                 <<std::setw(12)<<"true_anomaly" 
                 <<std::setw(12)<<"mean_anomaly" 
                 <<std::endl;                 
            

    PS::F64vec xcm((x[0][0]*m[0]+x[1][0]*m[1])/mt, 
                   (x[0][1]*m[0]+x[1][1]*m[1])/mt, 
                   (x[0][2]*m[0]+x[1][2]*m[1])/mt);

    PS::F64vec vcm((v[0][0]*m[0]+v[1][0]*m[1])/mt, 
                   (v[0][1]*m[0]+v[1][1]*m[1])/mt, 
                   (v[0][2]*m[0]+v[1][2]*m[1])/mt);

    return PtclHard(ib,mt,xcm,vcm,ppars.rout,0,0);
}


void print_p(PtclHard* p, const int n) {
    std::cout<<std::setw(12)<<"mass"
             <<std::setw(12)<<"x1"
             <<std::setw(12)<<"x2"
             <<std::setw(12)<<"x3"
             <<std::setw(12)<<"v1"
             <<std::setw(12)<<"v2"
             <<std::setw(12)<<"v3"
             <<std::setw(12)<<"rout"
             <<std::setw(12)<<"status"
             <<std::setw(12)<<"id"
             <<std::endl;
    for (int i=0; i<n; i++) {
        std::cout<<std::setw(12)<<p[i].mass
                 <<std::setw(12)<<p[i].pos[0]
                 <<std::setw(12)<<p[i].pos[1]
                 <<std::setw(12)<<p[i].pos[2]
                 <<std::setw(12)<<p[i].vel[0]
                 <<std::setw(12)<<p[i].vel[1]
                 <<std::setw(12)<<p[i].vel[2]
                 <<std::setw(12)<<p[i].r_out
                 <<std::setw(12)<<p[i].status
                 <<std::setw(12)<<p[i].id
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

  int N;
  params par;
  fs>>N>>par.rin>>par.rout>>par.gmin;
  PS::ReallocatableArray<PtclHard> p;
  p.resizeNoInitialize(N);
  
  ptree<PtclHard,params> plist;
  int idc=1;
  double m_average=0.0;
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

    PtclHard a(idc++,m1,xx1,vv1,par.rout,0,0,0);
    PtclHard b(idc++,m2,xx2,vv2,par.rout,0,0,0);

    bool flag=plist.link(id,ib,a,b,pshift);
    if (!flag) {
      std::cerr<<"Error: particle id "<<id<<", ib "<<ib<<" are inconsistent with global particle tree structure, cannot created pairs!\n";
      abort();
    }
  }

  int count=plist.collect_and_store(p.getPointer(),N);

  for(int i=0;i<N;i++) m_average += p[i].mass;
  m_average /=N;
  
  if (count<0) {
    std::cerr<<"Error: particle number mismatched particle tree!\n";
    abort();
  }
  
  //  plist.kepler_print(0,0,18,10);
  plist.pair_process(0,0,kepler_print,par);
  SearchGroup<PtclHard> group;
  print_p(p.getPointer(),N);

  group.findGroups(p.getPointer(), N);
  group.resolveGroups();

  group.searchAndMerge(p.getPointer(), N, par.rin);
  std::cout<<"SearchAndMerge\n";
  print_p(p.getPointer(),N);

  for(int i=0; i<group.getNGroups(); i++) {
      std::cout<<"group["<<i<<"]: ";
      for(int j=0; j<group.getGroupN(i); j++) {
          std::cout<<std::setw(10)<<group.getGroup(i)[j];
      }
      std::cout<<std::endl;
  }

  std::cout<<"Ptcl List:";
  for(int i=0; i<group.getNPtcl(); i++) {
      std::cout<<std::setw(10)<<group.getPtclList()[i];
  }
  std::cout<<std::endl;

  //PS::ReallocatableArray<PS::S32> adr_group_glb;
  //PS::ReallocatableArray<std::vector<PtclHard>> group_ptcl_glb; 
  //PS::ReallocatableArray<PS::S32> group_ptcl_glb_empty_list;
  PS::ReallocatableArray<PtclHard> ptcl_new;

  //group.generateList(p, N, adr_group_glb, group_ptcl_glb, group_ptcl_glb_empty_list);
  group.generateList<ptclTree>(p.getPointer(), N, ptcl_new, par.rin);
  std::cout<<"GenerateList\n";
  print_p(p.getPointer(),N);
  
  //std::cout<<"adr_group_glb: ";
  //for(int i=0; i<adr_group_glb.size(); i++) std::cout<<std::setw(6)<<adr_group_glb[i];
  //std::cout<<std::endl;
  
  std::cout<<"new ptcl:\n ";
  print_p(ptcl_new.getPointer(),ptcl_new.size());
  //for(int i=0; i<ptcl_new.size(); i++) {
  //    print_p(group_ptcl_glb[i].data(),group_ptcl_glb[i].size());
  //}
  p.reserveEmptyAreaAtLeast(ptcl_new.size());
  for(int i=0;i<ptcl_new.size();i++) p.pushBackNoCheck(ptcl_new[i]);

  group.findGroups(p.getPointer(),p.size());

  PS::ReallocatableArray<PtclHard> pn;
  PS::ReallocatableArray<PtclHard> pg;
  PS::ReallocatableArray<PtclHard> pp;
  std::cout<<"find groups\n";
  for(int i=0;i<group.getNPtcl();i++) pn.push_back(p[group.getPtclList()[i]]);
  print_p(pn.getPointer(),pn.size());
  
  for(int i=0;i<group.getNGroups();i++)  {
      pg.clearSize();
      std::cout<<"Group["<<i<<"]: ";
      for(int j=0;j<group.getGroupN(i);j++) 
          std::cout<<std::setw(10)<<group.getGroup(i)[j];
      std::cout<<std::endl;
      for(int j=0;j<group.getGroupN(i);j++) 
          pg.push_back(p[group.getGroup(i)[j]]);
      print_p(pg.getPointer(),pg.size());
      for(int j=0;j<group.getGroupN(i);j++) {
          std::cout<<"Member "<<j<<": ";
          for(int k=0;k<8;k++) std::cout<<std::setw(10)<<group.getGroupPertAdr(i,j,k);
          std::cout<<std::endl;
      }
      for(int j=0;j<16;j++) {
          pp.push_back(p[group.getGroupPertList(i)[j]]);
      }
      std::cout<<"Pert force:"<<std::endl;
      print_p(pp.getPointer(),pp.size());

  }
  group.resolveGroups();
  std::cout<<"resolvegroups\n";
  pn.clearSize();
  for(int i=0;i<group.getNPtcl();i++) pn.push_back(p[group.getPtclList()[i]]);
  print_p(pn.getPointer(),pn.size());
      
  
  return 0;
}
