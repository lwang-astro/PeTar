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
#include "kepler.hpp"
#include "cluster_list.hpp"
#include "hard.hpp"
#include "soft_ptcl.hpp"

//template<class Teng>
//void CalcEnergyHard(const PtclHard ptcl[], const PS::S32 n_tot, Teng & eng, 
//                    const PS::F64 r_in, const PS::F64 r_out, const PS::F64 eps_sq = 0.0){
//    eng.kin = eng.pot = eng.tot = 0.0;
//#ifndef INTEGRATED_CUTOFF_FUNCTION
//    PS::F64 r_oi_inv = 1.0/(r_out-r_in);
//    PS::F64 r_A = (r_out-r_in)/(r_out+r_in);
//    PS::F64 pot_off = cutoff_pot(1.0,r_oi_inv,r_A,r_in)/r_out;
//#endif
//    for(PS::S32 i=0; i<n_tot; i++){
//        eng.kin += 0.5 * ptcl[i].mass * ptcl[i].vel * ptcl[i].vel;
// 
//        for(PS::S32 j=i+1; j<n_tot; j++){
//            //PS::F64 r_out = std::max(ptcl[i].r_out,ptcl[j].r_out);
//            PS::F64vec rij = ptcl[i].pos - ptcl[j].pos;
//            PS::F64 dr = sqrt(rij*rij + eps_sq);
//#ifdef INTEGRATED_CUTOFF_FUNCTION
//            eng.pot -= ptcl[j].mass*ptcl[i].mass/dr*(1.0 - CalcW(dr/r_out, r_in/r_out));
//#else
//            if(dr<r_out) eng.pot -= ptcl[j].mass*ptcl[i].mass*(1.0/dr*cutoff_pot(dr, r_oi_inv, r_A, r_in) - pot_off);
//#endif
//        }
//    }
//    eng.tot = eng.kin + eng.pot;
//}

//const PS::F64 SAFTY_FACTOR_FOR_SEARCH = 1.05;
//const PS::F64 SAFTY_FACTOR_FOR_SEARCH_SQ = SAFTY_FACTOR_FOR_SEARCH * SAFTY_FACTOR_FOR_SEARCH;
//const PS::F64 SAFTY_OFFSET_FOR_SEARCH = 1e-7;

#ifndef NAN_CHECK
#define NAN_CHECK(val) assert((val) == (val));
#endif

//template<class Tptcl>
//void print_p(Tptcl* p, const int n, const int width= 20, const int precision = 14) {
//    std::cout<<std::setprecision(precision);
//    std::cout<<std::setw(width)<<"mass"
//             <<std::setw(width)<<"x1"
//             <<std::setw(width)<<"x2"
//             <<std::setw(width)<<"x3"
//             <<std::setw(width)<<"v1"
//             <<std::setw(width)<<"v2"
//             <<std::setw(width)<<"v3"
//             <<std::setw(width)<<"rsearch"
//             <<std::setw(width)<<"mass_bk"
//             <<std::setw(width)<<"status"
//             <<std::setw(width)<<"id"
//             <<std::setw(width)<<"id_cluster"
//             <<std::setw(width)<<"adr"
//             <<std::endl;
//    for (int i=0; i<n; i++) {
//        std::cout<<std::setw(width)<<p[i].mass
//                 <<std::setw(width)<<p[i].pos[0]
//                 <<std::setw(width)<<p[i].pos[1]
//                 <<std::setw(width)<<p[i].pos[2]
//                 <<std::setw(width)<<p[i].vel[0]
//                 <<std::setw(width)<<p[i].vel[1]
//                 <<std::setw(width)<<p[i].vel[2]
//                 <<std::setw(width)<<p[i].r_search
//                 <<std::setw(width)<<p[i].mass_bk
//                 <<std::setw(width)<<p[i].status
//                 <<std::setw(width)<<p[i].id
//                 <<std::setw(width)<<p[i].id_cluster
//                 <<std::setw(width)<<p[i].adr_org
//                 <<std::endl;
//    }
//}
// 
// 
//template<class Teng>
//void write_p(FILE* fout, const PS::F64 time, PtclHard* p, const int n, PtclHard& pcm, Teng &et, const PS::F64 rin, const PS::F64 rout, const PS::F64 eps2, const PS::F64 et0=0) {
//    fprintf(fout,"%20.14e ",time);
//    //PS::ReallocatableArray<PtclHard> pp;
//    //for (int i=0; i<n; i++) {
//    //    if(flag==2&&(p[i].status>0||p[i].id<0)) continue;
//    //    if(flag==1&&((p[i].id>0&&p[i].status!=0)||(p[i].id<=0&&p[i].status<0))) continue;
//    //    pp.push_back(p[i]);
//    //    if((flag==2||flag==1)&&p[i].status!=0) {
//    //        pp.back().mass = pp.back().mass_bk;
//    //    }
//    //    
//    //}
//    //for (int i=0; i<n; i++) p[i].mass = p[i].mass_bk;
//    calc_center_of_mass(pcm, p, n);
//    //printf("n = %d\n",pp.size());
//    Teng etb;
//    CalcEnergyHard(p, n, et,rin,rout,eps2);
//    PS::F64 err = et0==0?0:(et.tot-et0)/et0;
//    fprintf(fout,"%20.14e %20.14e %20.14e %20.14e ",err,et.kin,et.pot,et.tot);
//    for (int i=0; i<n; i++) {
//        fprintf(fout,"%20.14e %20.14e %20.14e %20.14e %20.14e %20.14e %20.14e ", 
//                p[i].mass, p[i].pos[0], p[i].pos[1], p[i].pos[2], 
//                p[i].vel[0], p[i].vel[1], p[i].vel[2]);
//    }
//    fprintf(fout,"\n");
//}
// 
//// flag: 1: c.m; 2: individual; 
//template<class Teng>
//void write_p(FILE* fout, const PS::F64 time, const PtclHard* p, const int n, const Teng *E, const Teng *E0=NULL) {
//    fprintf(fout,"%20.14e ",time);
//    //PS::ReallocatableArray<PtclHard> pp;
//    //for (int i=0; i<n; i++) {
//    //    if(p[i].status>0||p[i].id<0) continue;
//    //    pp.push_back(p[i]);
//    //    if(p[i].status!=0) pp.back().mass = pp.back().mass_bk;
//    //}
//    PS::F64 err = E0==NULL?0:(E->tot-E0->tot)/E0->tot;
//    fprintf(fout,"%20.14e %20.14e %20.14e %20.14e ",
//            err,E->kin,E->pot,E->tot);
//    for (int i=0; i<n; i++) {
//        fprintf(fout,"%20.14e %20.14e %20.14e %20.14e %20.14e %20.14e %20.14e ", 
//                p[i].mass, p[i].pos[0], p[i].pos[1], p[i].pos[2], 
//                p[i].vel[0], p[i].vel[1], p[i].vel[2]);
//    }
//    fprintf(fout,"\n");
//}

int main(int argc, char** argv)
{
  // data file name
  char* filename = argv[argc-1];
  //char* foutname = argv[argc-1];
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
//      sys[i].changeover.setR(1.5, rin, rout);
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
  hard_manager.r_tidal_tensor = rbin;
  hard_manager.r_in_base = sys[0].changeover.getRin();
  hard_manager.r_out_base = sys[0].changeover.getRout();
  hard_manager.id_offset = N;
  hard_manager.n_split = n_split;
  hard_manager.h4_manager.r_break_crit = rbin;
  hard_manager.h4_manager.r_neighbor_crit = rsearch;
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

//  EnergyAndMomemtum et0,et,etcm0,etcm;
//  PS::F64 eps2 = eps*eps;
//  PtclHard pcm0,pcm1;//,ppcm0,ppcm1;
//  FILE* fout;
//  std::string fname="hard.dat.";
//  if ( (fout = fopen((fname+foutname).c_str(),"w")) == NULL) {
//    fprintf(stderr,"Error: Cannot open file hard.dat.\n");
//    abort();
//  }
//  fprintf(stderr,"Time = %e\n", time_sys);
//  print_p(sys_hard.getPtcl().getPointer(),N);
//  write_p(fout,time_sys,sys_hard.getPtcl().getPointer(),sys_hard.getPtcl().size(),pcm0,et0,rin,rout,eps2,0.0);

  PS::S32 n_sys = sys.getNumberOfParticleLocal();
  sys_hard.findGroupsAndCreateArtificalParticlesOMP<PS::ParticleSystem<FPSoft>, FPSoft>(sys, dt_limit);
  
  kickCluster(sys, sys_hard.getPtcl(), 0.0);
  
//  fprintf(stderr,"After findgroup\n");
//  print_p(sys_hard.getPtcl().getPointer(),N);
 // FILE* fout2;
 // fname = "hardc.dat.";
 // if ( (fout2 = fopen((fname+foutname).c_str(),"w")) == NULL) {
 //   fprintf(stderr,"Error: Cannot open file hardcm.dat.\n");
 //   abort();
 // }

  std::cout<<std::setprecision(WRITE_PRECISION);

  //print_p(sys_hard.getPtcl().getPointer(),sys_hard.getPtcl().size());
  //write_p(fout,time_sys,sys_hard.getPtcl().getPointer(),sys_hard.getPtcl().size(), &sys_hard.ESD0);
  //write_p(fout2,time_sys,sys_hard.getPtcl().getPointer(),sys_hard.getPtcl().size(),ppcm0,etcm0,rin,rout,eps2,0.0,1);
  sys_hard.setTimeOrigin(time_sys);
  while(time_sys < time){
      fprintf(stderr,"Time = %e\n", time_sys+dt_limit);
      sys_hard.driveForMultiCluster(dt_limit, &sys[0]);
      time_sys += dt_limit;
      sys.setNumberOfParticleLocal(n_sys);
      sys_hard.setTimeOrigin(time_sys);
//      print_p(sys_hard.getPtcl().getPointer(),sys_hard.getPtcl().size());
      //write_p(fout,time_sys,sys_hard.getPtcl().getPointer(),sys_hard.getPtcl().size(), &sys_hard.ESD1, &sys_hard.ESD0);
//      write_p(fout,time_sys,sys_hard.getPtcl().getPointer(),sys_hard.getPtcl().size(),pcm1,et,rin,rout,eps2,et0.tot);
      //write_p(fout2,time_sys,sys_hard.getPtcl().getPointer(),sys_hard.getPtcl().size(),ppcm1,etcm,rin,rout,eps2,etcm0.tot,1);
//      std::cerr<<"CM: pos="<<pcm1.pos<<" vel="<<pcm1.vel<<" shift pos="<<pcm1.pos-pcm0.pos<<" shift vel="<<pcm1.vel-pcm0.vel<<std::endl;
      //std::cerr<<"CMHint: pos="<<ppcm1.pos<<" vel="<<ppcm1.vel<<" shift pos="<<ppcm1.pos-ppcm0.pos<<" shift vel="<<ppcm1.vel-ppcm0.vel<<std::endl;
//      std::cerr<<"Profile: ARC steps: "<<sys_hard.ARC_substep_sum<<std::endl;
      sys_hard.ARC_substep_sum = 0;
      sys_hard.findGroupsAndCreateArtificalParticlesOMP<PS::ParticleSystem<FPSoft>, FPSoft>(sys, dt_limit);
      kickCluster(sys, sys_hard.getPtcl(), 0.0);
  }
  

  fclose(fin);
//  fclose(fout);
//  fclose(fout2);
  
  return 0;
}

