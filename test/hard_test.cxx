#include <iostream>
#include <cstdio>
#include <iomanip>
#include <string>
#include <unordered_map>
#include <particle_simulator.hpp>
#define HARD_DEBUG_PRINT_FEQ 1

#include "io.hpp"
#include "integrate.hpp"
#include "hard_assert.hpp"
#include "cluster_list.hpp"
#include "hard.hpp"
#include "soft_ptcl.hpp"

#ifndef NAN_CHECK
#define NAN_CHECK(val) assert((val) == (val));
#endif

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
  PS::F64 rin, rout, rbin, rsearch, eps, eta, dt_limit, time;
  PS::S32 rcount = fscanf(fin, "%lf %d %lf %lf %lf %lf %lf %lf %lf\n", 
                          &time, &N, &rin, &rout, &rsearch, &rbin, &dt_limit, &eta, &eps);
  if (rcount<9) {
      std::cerr<<"Error: parameter reading fail!\n";
      abort();
  }
  Ptcl::r_search_min = rsearch;
  const PS::F64 n_split = 8;

  fprintf(stderr,"t_end = %e\nN = %d\nr_in = %e\nr_out = %e\neta = %e\ndt_limit = %e\neps = %e\n",time,N,rin,rout,eta,dt_limit,eps);

  ParticleBase pin;
  PS::ReallocatableArray<PS::S32> p_list;
  PS::ReallocatableArray<PS::S32> n_cluster;
  PS::ReallocatableArray<PS::S32> np;
  p_list.resizeNoInitialize(N);
  n_cluster.resizeNoInitialize(1);
  n_cluster[0] = N;

  PS::ParticleSystem<FPSoft> sys;

  PS::F64 m_average =0;
  for (int i=0; i<N; i++) {
      pin.readAscii(fin);
      ChangeOver co;
      sys.addOneParticle(FPSoft(Ptcl(pin, rsearch, 0.0, i+1, 0, co), 0, i));
      p_list[i]=i;
      m_average = pin.mass;
  }
  m_average = pin.mass/N;
  Ptcl::mean_mass_inv = 1.0/m_average;
  Ptcl::search_factor = 3;

  for (int i=0; i<N; i++) {
      sys[i].changeover.setR(sys[i].mass*Ptcl::mean_mass_inv, rin, rout);
  }

  PS::F64 time_sys = 0.0;
  PS::F64 dt_min_hard = 1.0;
  for (PS::S32 i=0;i<40;i++) dt_min_hard *= 0.5;

  // system hard paramters
  HardManager hard_manager;
  hard_manager.setDtRange(dt_limit, 40);
  hard_manager.setEpsSq(eps);
  hard_manager.setG(1.0);
#ifdef HARD_CHECK_ENERGY
  hard_manager.energy_error_max = 1e-4;
#else
  hard_manager.energy_error_max = NUMERIC_FLOAT_MAX;
#endif
  hard_manager.fp_manager.r_tidal_tensor = rbin;
  hard_manager.fp_manager.r_in_base = sys[0].changeover.getRin();
  hard_manager.fp_manager.r_out_base = sys[0].changeover.getRout();
  hard_manager.fp_manager.id_offset = N;
  hard_manager.fp_manager.n_split = n_split;
  hard_manager.h4_manager.step.eta_4th = eta*eta;
  hard_manager.h4_manager.step.eta_2nd = 0.01*eta*eta;
  hard_manager.h4_manager.step.calcAcc0OffsetSq(m_average, sys[0].changeover.getRout());
  hard_manager.ar_manager.energy_error_relative_max = 1e-8;
#ifdef AR_SYM
  hard_manager.ar_manager.step_count_max = 1e6;
#endif
  // set symplectic order
  hard_manager.ar_manager.step.initialSymplecticCofficients(-6);
  hard_manager.ar_manager.slowdown_pert_ratio_ref = 1e-4;
  hard_manager.ar_manager.slowdown_timescale_max = dt_limit;
  hard_manager.ar_manager.slowdown_mass_ref = m_average;

  // check consistence of paramters
  hard_manager.checkParams();
   
  SystemHard sys_hard;
  sys_hard.manager = &hard_manager;
  sys_hard.setPtclForIsolatedMultiCluster(sys, p_list, n_cluster);

  PS::S32 n_sys = sys.getNumberOfParticleLocal();
  sys_hard.findGroupsAndCreateArtificalParticlesOMP<PS::ParticleSystem<FPSoft>, FPSoft>(sys, dt_limit);
  
  PS::ReallocatableArray<PS::S32> remove_list;

  // correct change over
  auto& hard_ptcl = sys_hard.getPtcl();
  int n =hard_ptcl.size();
  for (int i=0; i<n; i++) hard_ptcl[i].changeover.updateWithRScale();

  kickCluster(sys, sys_hard.getPtcl(), 0.0);

  std::cout<<std::setprecision(WRITE_PRECISION);

  sys_hard.setTimeOrigin(time_sys);
  while(time_sys < time){
      fprintf(stderr,"Time = %e\n", time_sys+dt_limit);
      sys_hard.driveForMultiCluster(dt_limit, &sys[0]);
      sys_hard.writeBackPtclForMultiCluster(sys, remove_list);
      time_sys += dt_limit;
      sys.setNumberOfParticleLocal(n_sys);
      sys_hard.setTimeOrigin(time_sys);
      sys_hard.ARC_substep_sum = 0;
      sys_hard.findGroupsAndCreateArtificalParticlesOMP<PS::ParticleSystem<FPSoft>, FPSoft>(sys, dt_limit);
      for (int i=0; i<n; i++) hard_ptcl[i].changeover.updateWithRScale();
      kickCluster(sys, sys_hard.getPtcl(), 0.0);
  }
  

  fclose(fin);
  
  return 0;
}

