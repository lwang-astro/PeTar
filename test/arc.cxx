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
//#include "rsearch.hpp"
#include "cluster_list.hpp"

#ifndef FIX_STEP_DEBUG
#define STEP_DIVIDER 1
#endif  

//const PS::F64 SAFTY_FACTOR_FOR_SEARCH = 1.05;
//const PS::F64 SAFTY_FACTOR_FOR_SEARCH_SQ = SAFTY_FACTOR_FOR_SEARCH * SAFTY_FACTOR_FOR_SEARCH;
//const PS::F64 SAFTY_OFFSET_FOR_SEARCH = 1e-7;

class EnergyAndMomemtum{
public:
    PS::F64 kin;
    PS::F64 pot;
    PS::F64 tot;
    PS::F64vec L; // angular momentum
    PS::F64 Lt; // total angular momemtum

    EnergyAndMomemtum() {
        clear();
    }

    void clear(){
        kin = pot = tot = Lt = 0.0;
        L = PS::F64vec(0.0);
    }

    EnergyAndMomemtum operator -(const EnergyAndMomemtum& eng){
        EnergyAndMomemtum diff;
        diff.kin = this->kin - eng.kin;
        diff.pot = this->pot - eng.pot;
        diff.tot = this->tot - eng.tot;
        diff.L   = this->L   - eng.L;
        diff.Lt  = std::sqrt(diff.L*diff.L);
        return diff;
    }
};

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
void write_p(FILE* fout, const PS::F64 time, const Teng &E, const Teng &Ediff, const Tptcl* p, const int n) {
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

void print_p(PtclH4* p, const int n) {
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
  PS::F64 rin, rout, rsearch, rbin, eps, eta, dt_limit, time, m_average=0;
  PS::S32 rcount = fscanf(fin, "%lf %d %lf %lf %lf %lf %lf %lf %lf\n", 
                          &time, &N, &rin, &rout, &rsearch, &rbin, &dt_limit, &eta, &eps);
  if (rcount<9) {
      std::cerr<<"Error: parameter reading fail! rcount = "<<rcount<<", required 9 \n";
      abort();
  }
  Ptcl::r_search_min = rout;

  fprintf(stderr,"t_end = %e\nN = %d\nr_in = %e\nr_out = %e\neta = %e\ndt_limit = %e\neps = %e\n",time,N,rin,rout,eta,dt_limit,eps);

  PS::ReallocatableArray<ParticleBase> pin;
  PS::ReallocatableArray<PtclH4> p;
  PS::ReallocatableArray<PS::S32> adr;
  PS::ReallocatableArray<PS::S32> np;
  pin.resizeNoInitialize(N);
  for (int i=0; i<N; i++) {
      pin[i].readAscii(fin);
      p.push_back(PtclH4(pin[i]));
      p.back().r_search = rsearch;
      //p.back().mass_bk = 0.0;
      p.back().id = i+1;
      p.back().status = 0;
      m_average += p.back().mass;
      adr.push_back(i);
  }
  m_average /= N;

  PtclH4 pcm;
  calc_center_of_mass(pcm, p.getPointer(), N);
  center_of_mass_shift(pcm, p.getPointer(), N);

  print_p(p.getPointer(),N);

  const PS::S32 n_split = 8;

  SearchGroup<PtclH4> group;
  group.findGroups(p.getPointer(), N, n_split);
  group.searchAndMerge(p.getPointer(), rsearch);
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
  
  PS::ReallocatableArray<PtclH4> ptcl_new;

  group.generateList(p.getPointer(), ptcl_new, rsearch, rsearch, dt_limit, N, n_split);
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
  ARC_control.setA(Newtonian_cut_AW<PtclH4,ARC_pert_pars>,Newtonian_extA<PtclH4,PtclH4*,PtclForce*,ARC_pert_pars>,Newtonian_timescale<ARC_pert_pars>);
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

  ARCIntegrator<PtclH4, PtclH4, PtclForce> Aint(ARC_control, Int_pars);
  pcm = p[group.getPtclList()[0]];
  Aint.reserveARMem(1);
  Aint.reservePertMem(10);
  group.getBinPars(Aint.bininfo[0],p.getPointer(),0,n_split);
  Aint.addOneGroup(p.getPointer(),group.getGroup(0), group.getGroupN(0),group.getGroupPertList(0,n_split), n_split, &pcm, NULL, NULL, 0);

  std::cerr<<"Add group, N = "<<group.getGroupN(0)<<std::endl;
  
  Aint.initial();
  
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
  print_p(p.getPointer(),p.size());
  CalcEnergyHard(p.getPointer(),N,e0,rin,rout,eps*eps);
  //Aint.EnergyRecord(e0);
  write_p(fout,0.0,e0,ediff,p.getPointer(),N);
  PS::S64 stepcount = 0;
  for (int i=0; i<nstep; i++) {
      time_i += dt_limit/STEP_DIVIDER;
#ifdef ARC_SYM
      stepcount +=Aint.integrateOneStepSym(0, time_i, dt_limit/STEP_DIVIDER);
#else
      stepcount +=Aint.integrateOneStepExt(0, time_i, dt_limit/STEP_DIVIDER);
#endif
      if((i+1)%(int)STEP_DIVIDER==0) {
          printf("step_d=%d, i=%d\n",(int)STEP_DIVIDER,i);
          Aint.resolve();
      //Aint.EnergyRecord(e1);
          CalcEnergyHard(p.getPointer(),N,e1,rin,rout,eps*eps);
          //err = (e1.kin + e1.pot - (e0.kin + e0.pot))/(e0.kin + e0.pot);
          ediff = e1 - e0;
          write_p(fout,time_i,e1,ediff,p.getPointer(),N);
          std::cerr<<"Time = "<<time_i<<std::endl;
          print_p(p.getPointer(),p.size());
          Aint.shift();
      }
  }

  printf("Energy error: %e, kin: %e, pot: %e, init: %e, end: %e, nstep: %lld\n",ediff.tot/e1.tot, e1.kin, e1.pot, e0.tot, e1.tot, stepcount);

  fclose(fout);
  
  return 0;

}
