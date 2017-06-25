#include <iostream>
#include <cstdio>
#include <iomanip>
#include <unordered_map>
#include <particle_simulator.hpp>
#include "Newtonian_acceleration.h"
#include "ptree.h"
#include "kepler.hpp"
#include "rsearch.hpp"
#include "cluster_list.hpp"
#include "hard.hpp"

const PS::F64 SAFTY_FACTOR_FOR_SEARCH = 1.05;
const PS::F64 SAFTY_FACTOR_FOR_SEARCH_SQ = SAFTY_FACTOR_FOR_SEARCH * SAFTY_FACTOR_FOR_SEARCH;
const PS::F64 SAFTY_OFFSET_FOR_SEARCH = 1e-7;

#ifndef NAN_CHECK
#define NAN_CHECK(val) assert((val) == (val));
#endif

void print_p(PtclHard* p, const int n) {
    std::cout<<std::setw(12)<<"mass"
             <<std::setw(12)<<"x1"
             <<std::setw(12)<<"x2"
             <<std::setw(12)<<"x3"
             <<std::setw(12)<<"v1"
             <<std::setw(12)<<"v2"
             <<std::setw(12)<<"v3"
             <<std::setw(12)<<"rsearch"
             <<std::setw(12)<<"mass_bk"
             <<std::setw(12)<<"status"
             <<std::setw(12)<<"id"
             <<std::setw(12)<<"id_cluster"
             <<std::setw(12)<<"adr"
             <<std::endl;
    for (int i=0; i<n; i++) {
        std::cout<<std::setw(12)<<p[i].mass
                 <<std::setw(12)<<p[i].pos[0]
                 <<std::setw(12)<<p[i].pos[1]
                 <<std::setw(12)<<p[i].pos[2]
                 <<std::setw(12)<<p[i].vel[0]
                 <<std::setw(12)<<p[i].vel[1]
                 <<std::setw(12)<<p[i].vel[2]
                 <<std::setw(12)<<p[i].r_search
                 <<std::setw(12)<<p[i].mass_bk
                 <<std::setw(12)<<p[i].status
                 <<std::setw(12)<<p[i].id
                 <<std::setw(12)<<p[i].id_cluster
                 <<std::setw(12)<<p[i].adr_org
                 <<std::endl;
    }
}

int main(int argc, char** argv)
{
  // data file name
  char* filename = argv[argc-1];
  // open data file

  FILE* fin;
  if ( (fin = fopen(filename,"r")) == NULL) {
    fprintf(stderr,"Error: Cannot open input file %s.\n",filename);
    abort();
  }

  int N;
  PS::F64 rin, rout, rsearch, gmin, eps, eta, dt_limit, time, m_average=0;
  PS::S32 rcount = fscanf(fin, "%lf %d %lf %lf %lf %lf %lf %lf\n", 
                          &time, &N, &rin, &rout, &dt_limit, &eta, &eps, &gmin);
  if (rcount<8) {
      std::cerr<<"Error: parameter reading fail!\n";
      abort();
  }
  rsearch = rout;

  fprintf(stderr,"t_end = %e\nN = %d\nr_in = %e\nr_out = %e\neta = %e\ndt_limit = %e\neps = %e\n",time,N,rin,rout,eta,dt_limit,eps);

  PS::ReallocatableArray<ParticleBase> pin;
  PS::ReallocatableArray<PtclHard> p;
  PS::ReallocatableArray<PS::S32> adr;
  PS::ReallocatableArray<PS::S32> np;
  pin.resizeNoInitialize(N);
  for (int i=0; i<N; i++) {
      pin[i].readAscii(fin);
      p.push_back(PtclHard(pin[i]));
      p.back().r_search = rsearch;
      //p.back().mass_bk = 0.0;
      p.back().id = i+1;
      p.back().status = 0;
      m_average += p.back().mass;
      adr.push_back(i);
  }
  m_average /= N;
  np.push_back(N);

  print_p(p.getPointer(),N);
  
  SystemHard sys;
  PS::F64 time_sys = 0.0;
  sys.setParam(rout, rin, eps, dt_limit, eta, time_sys, gmin, m_average);
  sys.setARCParam();
  
  sys.setPtclForIsolatedMultiCluster(p,adr,np);
  while(time_sys < time){
      sys.driveForMultiClusterOMP(dt_limit);
      time_sys += dt_limit;
  }
  
  return 0;
}

