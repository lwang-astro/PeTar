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

//const PS::F64 SAFTY_FACTOR_FOR_SEARCH = 1.05;
//const PS::F64 SAFTY_FACTOR_FOR_SEARCH_SQ = SAFTY_FACTOR_FOR_SEARCH * SAFTY_FACTOR_FOR_SEARCH;
//const PS::F64 SAFTY_OFFSET_FOR_SEARCH = 1e-7;

template <class Tptcl>
void write_p(FILE* fout, const PS::F64 time, const HardEnergy &E, const PS::F64 err, const Tptcl* p, const int n) {
    fprintf(fout,"%e %e %e %e %e ",time, err, E.kin, E.pot, E.tot);
    //PS::ReallocatableArray<PtclHard> pp;
    //for (int i=0; i<n; i++) {
    //    if(flag==2&&(p[i].status>0||p[i].id<0)) continue;
    //    if(flag==1&&(p[i].id>=0||p[i].status<0)) continue;
    //    pp.push_back(p[i]);
    //    if((flag==2&&p[i].status!=0)||flag==1) pp.back().mass = pp.back().mass_bk;
    //}
    //for (int i=0; i<pp.size(); i++) {
    for (int i=0; i<n; i++) {
        fprintf(fout,"%e %e %e %e %e %e %e ", 
                p[i].mass, p[i].pos[0], p[i].pos[1], p[i].pos[2], 
                p[i].vel[0], p[i].vel[1], p[i].vel[2]);
    }
    fprintf(fout,"\n");
}

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
  PS::F64 rin, rout, rsearch, eps, eta, dt_limit, time, m_average=0;
  PS::S32 rcount = fscanf(fin, "%lf %d %lf %lf %lf %lf %lf\n", 
                          &time, &N, &rin, &rout, &dt_limit, &eta, &eps);
  if (rcount<7) {
      std::cerr<<"Error: parameter reading fail!\n";
      abort();
  }
  rsearch = 10*rout;

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

  print_p(p.getPointer(),N);

  const PS::S32 n_split = 8;

  SearchGroup<PtclHard> group;
  group.findGroups(p.getPointer(), N, n_split);
  group.searchAndMerge(p.getPointer(), N, rsearch);
  //std::cout<<"SearchAndMerge\n";
  // 
  //for(int i=0; i<group.getNGroups(); i++) {
  //    std::cout<<"group["<<i<<"]: ";
  //    for(int j=0; j<group.getGroupN(i); j++) {
  //        std::cout<<std::setw(10)<<group.getGroup(i)[j];
  //    }
  //    std::cout<<std::endl;
  //}

  //std::cout<<"Ptcl List:";
  //for(int i=0; i<group.getNPtcl(); i++) {
  //    std::cout<<std::setw(10)<<group.getPtclList()[i];
  //}
  //std::cout<<std::endl;
  
  PS::ReallocatableArray<PtclHard> ptcl_new;

  group.generateList(p.getPointer(), N, ptcl_new, rsearch, n_split);
  //std::cout<<"GenerateList\n";
  //print_p(p.getPointer(),p.size());

  //std::cout<<"new ptcl: "<<ptcl_new.size()<<"\n ";
  //print_p(ptcl_new.getPointer(),ptcl_new.size());

  p.reserveEmptyAreaAtLeast(ptcl_new.size());
  for(int i=0;i<ptcl_new.size();i++) {
      p.pushBackNoCheck(ptcl_new[i]);
      adr.push_back(i+N);
  }

  np.push_back(p.size());

  std::cout<<"new p: "<<p.size()<<"\n ";
  print_p(p.getPointer(),np[0]);

  group.findGroups(p.getPointer(),p.size(),n_split);

  ARC::chainpars ARC_control;
  ARC_control.setA(Newtonian_cut_AW<PtclHard,ARC_pert_pars>,Newtonian_extA<PtclHard,PtclHard*,PtclForce*,ARC_pert_pars>,Newtonian_timescale<ARC_pert_pars>);
  ARC_control.setabg(0,1,0);
  ARC_control.setErr(1e-10,1e-24,1e-6);
  ARC_control.setIterSeq(20,20);
  ARC_control.setIntp(1);
  ARC_control.setIterConst(0);
  ARC_control.setAutoStep(3);

  ARC_int_pars Int_pars;
  Int_pars.rin = rin;
  Int_pars.rout = rout;
  Int_pars.eps2 = eps*eps;

  ARCIntegrator<PtclHard, PtclHard, PtclForce, ARC_int_pars, ARC_pert_pars> Aint(ARC_control, Int_pars);
  PtclHard pcm = p[group.getPtclList()[0]];
  Aint.reserveARMem(1);
  Aint.reservePertMem(10);
  Aint.addOneGroup(p.getPointer(),group.getGroup(0), group.getGroupN(0),group.getGroupPertList(0,n_split), n_split, &pcm, NULL, NULL, 0);

  std::cerr<<"Add group, N = "<<group.getGroupN(0)<<std::endl;
  
  Aint.initial();
  
  HardEnergy e0,e1;

  FILE* fout;
  if ( (fout = fopen("arc.dat","w")) == NULL) {
    fprintf(stderr,"Error: Cannot open file hard.dat.\n");
    abort();
  }
  
  PS::S32 nstep = time/dt_limit;
  PS::F64 time_i = 0.0;
  PS::F64 err;

  Aint.EnergyRecord(e0);
  write_p(fout,0.0,e0,0.0,p.getPointer(),N);
  for (int i=0; i<nstep; i++) {
      time_i += dt_limit;
      Aint.integrateOneStep(0, time_i, dt_limit);
      Aint.resolve();
      Aint.EnergyRecord(e1);
      err = (e1.kin + e1.pot - (e0.kin + e0.pot))/(e0.kin + e0.pot);
      write_p(fout,time_i,e1,err,p.getPointer(),N);
      //std::cerr<<"Time = "<<time_i<<std::endl;
      //print_p(p.getPointer(),p.size());
      Aint.shift();
  }

  printf("Energy error: %e, kin: %e, pot: %e, init: %e, end: %e\n",err, e1.kin, e1.pot, e0.tot, e1.tot);

  fclose(fout);
  
  return 0;

}
