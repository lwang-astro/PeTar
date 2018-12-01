#include <iostream>
#include <cstdio>
#include <iomanip>
#include <string>
#include <particle_simulator.hpp>
#include "ptree.h"
#include "kepler.hpp"
#include "hard.hpp"

void print_ptcl(Ptcl* p, const int n) {
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
                 <<std::endl;
    }
}

void print_tree(PtclTree<Ptcl> &bins) {
    if(bins.member[0]->status!=0) print_tree(*(PtclTree<Ptcl>*)bins.member[0]);
    if(bins.member[1]->status!=0) print_tree(*(PtclTree<Ptcl>*)bins.member[1]);

    std::cout<<std::setw(18)<<bins.semi
             <<std::setw(18)<<bins.ecc
             <<std::setw(18)<<bins.peri
             <<std::setw(18)<<bins.m1
             <<std::setw(18)<<bins.m2
             <<std::setw(18)<<bins.tstep
             <<std::setw(18)<<bins.fpert
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
  // time as tree step used to stability check
  PS::F64 rin, rout, rsearch, rbin, eps, eta, dt_limit, time;
  PS::S32 rcount = fscanf(fin, "%lf %d %lf %lf %lf %lf %lf %lf %lf\n", 
                          &time, &N, &rin, &rout, &rsearch, &rbin, &dt_limit, &eta, &eps);
  if (rcount<9) {
      std::cerr<<"Error: parameter reading fail! rcount = "<<rcount<<", required 9 \n";
      abort();
  }

  Ptcl::r_search_min = rout;

  fprintf(stderr,"N = %d\nr_in = %e\nr_out = %e\nr_bin = %e\nr_search = %e\ndt_tree = %e\n",N,rin,rout,rbin,rsearch,time);
  ParticleBase ptcl_in;
  PS::ReallocatableArray<Ptcl> ptcl;
  PS::ReallocatableArray<PS::S32> ptcl_list;

  ptcl.resizeNoInitialize(N);
  ptcl_list.resizeNoInitialize(N);

  for (int i=0; i<N; i++) {
      ptcl_in.readAscii(fin);
      ptcl[i]=Ptcl(ptcl_in, rsearch, 0.0, i+1, 0);
      ptcl[i].status=0;
      ptcl_list[i]=i;
  }

  print_ptcl(ptcl.getPointer(),N);

  PtclTree<Ptcl> ptcl_tree[N-1];
  keplerTreeGenerator<Ptcl>(ptcl_tree, ptcl_list.getPointer(), N, ptcl.getPointer(), dt_limit, 100.0);

  for (PS::S32 i=0; i<N; i++) ptcl[i].status=0;

  std::cout<<std::setw(18)<<"semi"
           <<std::setw(18)<<"ecc"
           <<std::setw(18)<<"peri"
           <<std::setw(18)<<"m1"
           <<std::setw(18)<<"m2"
           <<std::setw(18)<<"tstep"
           <<std::setw(18)<<"fpert"
           <<std::endl;
  print_tree(ptcl_tree[N-2]);

  PS::ReallocatableArray<PtclTree<Ptcl>*> stab_bins; // stable checked binary tree
  stab_bins.reserve(N);

  stabilityCheck<Ptcl>(stab_bins, ptcl_tree[N-2], rbin, rin, rout, time);
  
  return 0;
}

