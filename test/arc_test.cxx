//arctest [input] [(arc.dat.)outname]

#include <iostream>
#include <cstdio>
#include <string>
#include <iomanip>
#include <unordered_map>
#include <particle_simulator.hpp>
//#include "Newtonian_acceleration.h"
#include "hard.hpp"
#include "ptree.h"
#include "kepler.hpp"
#include "soft.hpp"
//#include "rsearch.hpp"
#include "cluster_list.hpp"

#ifndef FIX_STEP_DEBUG
#define STEP_DIVIDER 1
#endif  

//const PS::F64 SAFTY_FACTOR_FOR_SEARCH = 1.05;
//const PS::F64 SAFTY_FACTOR_FOR_SEARCH_SQ = SAFTY_FACTOR_FOR_SEARCH * SAFTY_FACTOR_FOR_SEARCH;
//const PS::F64 SAFTY_OFFSET_FOR_SEARCH = 1e-7;

//class EnergyAndMomemtum{
//public:
//    PS::F64 kin;
//    PS::F64 pot;
//    PS::F64 tot;
//    PS::F64vec L; // angular momentum
//    PS::F64 Lt; // total angular momemtum
// 
//    EnergyAndMomemtum() {
//        clear();
//    }
// 
//    void clear(){
//        kin = pot = tot = Lt = 0.0;
//        L = PS::F64vec(0.0);
//    }
// 
//    EnergyAndMomemtum operator -(const EnergyAndMomemtum& eng){
//        EnergyAndMomemtum diff;
//        diff.kin = this->kin - eng.kin;
//        diff.pot = this->pot - eng.pot;
//        diff.tot = this->tot - eng.tot;
//        diff.L   = this->L   - eng.L;
//        diff.Lt  = std::sqrt(diff.L*diff.L);
//        return diff;
//    }
//};

template<class Tptcl, class Teng>
void CalcEnergyHard(const Tptcl ptcl[], const PS::S32 n_tot, Teng & eng, 
                    const PS::F64 r_in, const PS::F64 r_out, const PS::F64 eps_sq = 0.0){
    eng.clear();
#ifndef INTEGRATED_CUTOFF_FUNCTION
    PS::F64 r_oi_inv = 1.0/(r_out-r_in);
    PS::F64 r_A = (r_out-r_in)/(r_out+r_in);
    PS::F64 pot_off = (1.0 + r_A)/r_out;
#endif
    for(PS::S32 i=0; i<n_tot; i++){
        eng.kin += 0.5 * ptcl[i].mass * ptcl[i].vel * ptcl[i].vel;
        eng.L   += ptcl[i].pos ^ (ptcl[i].mass*ptcl[i].vel);

        for(PS::S32 j=i+1; j<n_tot; j++){
            //PS::F64 r_out = std::max(ptcl[i].r_out,ptcl[j].r_out);
            PS::F64vec rij = ptcl[i].pos - ptcl[j].pos;
            PS::F64 dr = sqrt(rij*rij + eps_sq);
#ifdef INTEGRATED_CUTOFF_FUNCTION
            eng.pot -= ptcl[j].mass*ptcl[i].mass/dr*(1.0 - CalcW(dr/r_out, r_in/r_out));
#else
            if(dr<r_out) eng.pot -= ptcl[j].mass*ptcl[i].mass*(1.0/dr*cutoff_pot(dr, r_oi_inv, r_A, r_in) - pot_off);
#endif
        }
    }
    eng.tot = eng.kin + eng.pot;
    eng.Lt  = std::sqrt(eng.L*eng.L);
}

template <class Tptcl, class Teng>
void write_ptcl(FILE* fout, const PS::F64 time, const Teng &E, const Teng &Ediff, const Tptcl* p, const int n) {
    fprintf(fout,"%20.14e %20.14e %20.14e %20.14e %20.14e %20.14e %20.14e %20.14e %20.14e %20.14e %20.14e %20.14e %20.14e ",
            time, Ediff.tot/E.tot, E.kin, E.pot, E.tot,
            Ediff.Lt/E.Lt, Ediff.L[0]/E.Lt, Ediff.L[1]/E.Lt, Ediff.L[2]/E.Lt,
            E.Lt,    E.L[0],    E.L[1],    E.L[2]);
    //PS::ReallocatableArray<PtclH4> pp;
    //for (int i=0; i<n; i++) {
    //    if(flag==2&&(p[i].status>0||p[i].id<0)) continue;
    //    if(flag==1&&(p[i].id>=0||p[i].status<0)) continue;
    //    pp.push_back(p[i]);
    //    if((flag==2&&p[i].status!=0)||flag==1) pp.back().mass = pp.back().mass_bk;
    //}
    //for (int i=0; i<pp.size(); i++) {
    for (int i=0; i<n; i++) {
        fprintf(fout,"%20.14e %20.14e %20.14e %20.14e %20.14e %20.14e %20.14e ", 
                p[i].mass, p[i].pos[0], p[i].pos[1], p[i].pos[2], 
                p[i].vel[0], p[i].vel[1], p[i].vel[2]);
    }
    fprintf(fout,"\n");
}

#ifndef NAN_CHECK
#define NAN_CHECK(val) assert((val) == (val));
#endif

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

int main(int argc, char** argv)
{
  // data file name
  char* filename = argv[argc-2];
  char* foutname = argv[argc-1];
  // open data file

  FILE* fin;
  if ( (fin = fopen(filename,"r")) == NULL) {
    fprintf(stderr,"Error: Cannot open input file %s.\n",filename);
    abort();
  }

  int N;
  PS::F64 rin, rout, rsearch, rbin, eps, eta, dt_limit, time;
  PS::S32 rcount = fscanf(fin, "%lf %d %lf %lf %lf %lf %lf %lf %lf\n", 
                          &time, &N, &rin, &rout, &rsearch, &rbin, &dt_limit, &eta, &eps);
  if (rcount<9) {
      std::cerr<<"Error: parameter reading fail! rcount = "<<rcount<<", required 9 \n";
      abort();
  }
  Ptcl::r_search_min = rout;

  fprintf(stderr,"t_end = %e\nN = %d\nr_in = %e\nr_out = %e\neta = %e\ndt_limit = %e\neps = %e\n",time,N,rin,rout,eta,dt_limit,eps);

  ParticleBase ptcl_in;
  PS::ReallocatableArray<Ptcl> ptcl;
  PS::ReallocatableArray<PS::S32> ptcl_list;
  //PS::ReallocatableArray<PS::S32> n_cluster;
  //PS::ReallocatableArray<PS::S32> np;
  ptcl.resizeNoInitialize(N);
  ptcl_list.resizeNoInitialize(N);
  //n_cluster.resizeNoInitialize(1);
  //n_cluster[0] = N;
  for (int i=0; i<N; i++) {
      ptcl_in.readAscii(fin);
      ptcl[i]=Ptcl(ptcl_in, rsearch, 0.0, i+1, 0);
      ptcl_list[i]=i;
  }

  //PS::ParticleSystem<FPSoft> sys;
  //for (int i=0; i<p.size(); i++) sys.addOneParticle(FPSoft(p[i],0,i));
  
  const PS::S32 n_split = 8;

  //SystemHard sys_hard;
  //sys_hard.setParam(rbin, rout, rin, eps, dt_limit, dt_limit, eta, 0.0, 0.0, id_offset, n_split);
  // 
  //sys_hard.setPtclForIsolatedMultiCluster(sys, p_list, n_cluster);
  //sys_hard.findGroupsAndCreateArtificalParticlesOMP<PS::ParticleSystem<FPSoft>, FPSoft>(sys, dt_limit);

  print_ptcl(ptcl.getPointer(),N);

  PtclTree<Ptcl> ptcl_tree[N-1];
  keplerTreeGenerator<Ptcl>(ptcl_tree, ptcl_list.getPointer(), N, ptcl.getPointer(), dt_limit, 100.0);
  PS::F64 arc_step=PS::LARGE_FLOAT;
  for (PS::S32 i=0; i<N-1; i++) {
      const PS::F64 semi = ptcl_tree[i].semi;
      const PS::F64 m1 = ptcl_tree[i].m1;
      const PS::F64 m2 = ptcl_tree[i].m2;
      const PS::F64 mcm = m1+m2;
      ptcl_tree[i].tstep=0.78539816339*std::sqrt(semi/mcm)*m1*m2;  
      arc_step = std::min(arc_step, ptcl_tree[i].tstep);
  }
  //SearchGroup<Ptcl> group;
  
  //group.findGroups(p.getPointer(), N, n_split);
  //group.searchAndMerge(p.getPointer(), N, rsearch);
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
  
  //PS::ReallocatableArray<Ptcl> ptcl_new;
  //PS::S32 n_group;
  //group.generateList(0, p.getPointer(), N, ptcl_new, n_group, rbin, rin, rout, dt_limit, N, n_split);
  //std::cout<<"GenerateList\n";
  //print_p(p.getPointer(),p.size());

  //std::cout<<"new ptcl: "<<ptcl_new.size()<<"\n ";
  //print_p(ptcl_new.getPointer(),ptcl_new.size());


  //PS::ReallocatableArray<FPSoft> psys;
  //psys.reserveEmptyAreaAtLeast(ptcl_new.size());
  //for(int i=0;i<ptcl_new.size();i++) {
  //psys.pushBackNoCheck(FPSoft(ptcl_new[i],0,i));
  ////    adr.push_back(i+N);
//}

  //np.push_back(p.size());

  //std::cout<<"create artifical particles: N= "<<ptcl_new.size()<<"\n ";
  //print_p(ptcl_new.getPointer(),ptcl_new.size());

  //group.findGroups(p.getPointer(),p.size(),n_split);

  ARC::chainpars ARC_control;
  ARC_control.setA(Newtonian_cut_AW<Ptcl,ARC_pert_pars>,Newtonian_extA_soft<Ptcl,PtclH4*,PtclForce*,ARC_pert_pars>,Newtonian_timescale<ARC_pert_pars>);
  ARC_control.setabg(0,1,0);
  ARC_control.setErr(1e-10,1e-24,1e-6);
#ifdef ARC_SYM
  ARC_control.setSymOrder(-6);
#else
  ARC_control.setIterSeq(20,3,20);
#endif
  ARC_control.setIntp(1);
  ARC_control.setIterConst(0);
  ARC_control.setAutoStep(3);

  ARC_control.print(std::cerr);

  ARC_int_pars Int_pars;
  Int_pars.rin = rin;
  Int_pars.rout = rout;
  Int_pars.r_oi_inv = 1.0/(rout-rin);
  Int_pars.r_A = (rout-rin)/(rout+rin);
  Int_pars.pot_off = (1.0 + Int_pars.r_A)/rout;
  Int_pars.eps2 = eps*eps;

  ARCIntegrator<Ptcl, PtclH4, PtclForce> Aint(ARC_control, Int_pars);
  //GroupPars gpars(n_split);
  //gpars.getGroupIndex(ptcl_new.getPointer());
  //Ptcl* pcm = &ptcl_new[gpars.offset_cm];
#ifdef ARC_SYM
  Aint.step_count_limit = 10000;
#endif
  Aint.reserveARMem(1);
  Aint.reservePertMem(1,1);
  //gpars.getBinPars(Aint.bininfo[0],ptcl_new.getPointer());
  Aint.bininfo[0].tstep = arc_step;
  Aint.addOneGroup(ptcl.getPointer(), NULL, N, (FPSoft*)NULL, n_split);
  //PS::S32 iact=0;
  //Aint.updateCM(pcm, &iact, 1);

  //Aint.initialSlowDown(_time_end, sdfactor_);
  Aint.initial();
  Ptcl ptcl_cm = *Aint.getCM(0);
  ptcl_cm.print(std::cout);
  std::cout<<std::endl;
  
  EnergyAndMomemtum e0,e1,ediff;

  FILE* fout;
  std::string fname="arc.dat.";
  if ( (fout = fopen((fname+foutname).c_str(),"w")) == NULL) {
      fprintf(stderr,"Error: Cannot open file arc.dat.\n");
      abort();
  }

  PS::S32 nstep = time/dt_limit*STEP_DIVIDER;
  PS::F64 time_i = 0.0;
  // PS::F64 err;

  std::cerr<<"Time = "<<time_i<<std::endl;
  print_ptcl(ptcl.getPointer(),N);
  CalcEnergyHard(ptcl.getPointer(),N,e0,rin,rout,eps*eps);
  //Aint.EnergyRecord(e0);
  write_ptcl(fout,0.0,e0,ediff,ptcl.getPointer(),N);
  PS::S64 stepcount = 0;
  for (int i=0; i<nstep; i++) {
      time_i += dt_limit/STEP_DIVIDER;
      PS::S64 stepcount_i;
#ifdef ARC_SYM
      stepcount_i =Aint.integrateOneStepSym(0, time_i, dt_limit/STEP_DIVIDER);
#else
      stepcount_i =Aint.integrateOneStepExt(0, time_i, dt_limit/STEP_DIVIDER);
#endif
      if(stepcount_i<0) {
          std::cout<<"Error happen!"<<std::endl;
          abort();
      }
      else stepcount += stepcount_i;
      if((i+1)%(int)STEP_DIVIDER==0) {
          printf("step_d=%d, i=%d\n",(int)STEP_DIVIDER,i);
          Aint.resolve();
      //Aint.EnergyRecord(e1);
          CalcEnergyHard(ptcl.getPointer(),N,e1,rin,rout,eps*eps);
          //err = (e1.kin + e1.pot - (e0.kin + e0.pot))/(e0.kin + e0.pot);
          ediff = e1 - e0;
          write_ptcl(fout,time_i,e1,ediff,ptcl.getPointer(),N);
          std::cerr<<"Time = "<<time_i<<" Energy Error = "<<ediff.tot/e1.tot<<std::endl;
          print_ptcl(ptcl.getPointer(),N);
          Aint.shift2CM();
      }
  }

  printf("Energy error: %e, kin: %e, pot: %e, init: %e, end: %e, nstep: %lld\n",ediff.tot/e1.tot, e1.kin, e1.pot, e0.tot, e1.tot, stepcount);

  fclose(fout);
  
  return 0;

}
