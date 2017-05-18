#include <iostream>
#include <fstream>
#include <particle_simulator.hpp>
#include "Newtonian_acceleration.h"
#include "ptree.h"
#include "kepler.hpp"
#include "rsearch.hpp"
#include "hard.hpp"

struct params{
    double rin,rout,gmin;
};

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
  PtclHard* p[N];
  
  ptree<PtclHard,params> plist;
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
    PtclHard a(i,m1,xx1,vv1,par.rout,0,0);
    PtclHard b(-i,m2,xx2,vv2,par.rout,0,0);
    bool flag=plist.link(id,ib,a,b);
    if (!flag) {
      std::cerr<<"Error: particle id "<<id<<", ib "<<ib<<" are inconsistent with global particle tree structure, cannot created pairs!\n";
      abort();
    }
  }

  int count=plist.collect(p,N);

  for(int i=0;i<N;i++) m_average += p[i]->mass;
  m_average /=N;
  
  if (count<0) {
    std::cerr<<"Error: particle number mismatched particle tree!\n";
    abort();
  }
  
  //  plist.kepler_print(0,0,18,10);

  updateRout(p, N, par.rin, par.rout, par.gmin, m_average);
  for (int i=0; i<N; i++) {
    std::cout<<i<<" "<<p[i]->r_out<<std::endl;
  }

  return 0;
}
