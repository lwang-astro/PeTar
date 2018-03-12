#include <iostream>
#include <cstdio>
#include <iomanip>
#include <string>
#include <particle_simulator.hpp>
#include "ptree.h"
#include "kepler.hpp"
#include "hard.hpp"

void print_tree(PtclTree<PtclHard> &bins) {
    if(bins.member[0]->status!=0) print_tree(*(PtclTree<PtclHard>*)bins.member[0]);
    if(bins.member[1]->status!=0) print_tree(*(PtclTree<PtclHard>*)bins.member[1]);

    std::cout<<std::setw(18)<<bins.ax
             <<std::setw(18)<<bins.ecc
             <<std::setw(18)<<bins.peri
             <<std::setw(18)<<bins.m1
             <<std::setw(18)<<bins.m2
             <<std::setw(18)<<bins.tstep
             <<std::endl;
}

int main(int argc, char** argv)
{
  // data file name
  char* filename = argv[1];
  // open data file

  FILE* fin;
  if ( (fin = fopen(filename,"r")) == NULL) {
    fprintf(stderr,"Error: Cannot open input file %s.\n",filename);
    abort();
  }

  int N;
  PS::F64 rin, rout, rbin, rsearch, eps, eta, dt_limit, time;
  PS::S32 rcount = fscanf(fin, "%lf %d %lf %lf %lf %lf %lf %lf %lf\n", 
                          &time, &N, &rin, &rout, &rsearch, &rbin, &dt_limit, &eta, &eps);
  if (rcount<8) {
      std::cerr<<"Error: parameter reading fail!\n";
      abort();
  }

  PS::ReallocatableArray<ParticleBase> pin;
  PS::ReallocatableArray<PtclHard> p;

  pin.resizeNoInitialize(N);
  for (int i=0; i<N; i++) {
      pin[i].readAscii(fin);
      p.push_back(PtclHard(pin[i]));
      p.back().r_search = rsearch;
      //p.back().mass_bk = 0.0;
      p.back().id = i+1;
      p.back().status = 0;
  }

  PS::ReallocatableArray<PtclTree<PtclHard>> bins;
  PS::ReallocatableArray<PtclTree<PtclHard>*> stab_bins;
  bins.reserve(1);
  stab_bins.reserve(1);
  bins.resizeNoInitialize(N-1);
  
  PS::S32 list[N];
  for (int i=0; i<N; i++) list[i]=i;
  
  keplerTreeGenerator(bins.getPointer(), list, N, p.getPointer(), dt_limit);

  std::cout<<std::setw(18)<<"A"
           <<std::setw(18)<<"ecc"
           <<std::setw(18)<<"peri"
           <<std::setw(18)<<"m1"
           <<std::setw(18)<<"m2"
           <<std::setw(18)<<"tstep"
           <<std::endl;
  print_tree(bins.back());

  PS::F64 fstab = stabilityCheck<PtclHard>(stab_bins,bins.back(),rbin,rin,rout);
  
  return 0;
}

