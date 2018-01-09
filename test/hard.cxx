#include <iostream>
#include <cstdio>
#include <iomanip>
#include <string>
#include <unordered_map>
#include <particle_simulator.hpp>
#include "Newtonian_acceleration.h"
#include "ptree.h"
#include "kepler.hpp"
#include "cluster_list.hpp"
#include "hard.hpp"
#include "soft.hpp"


template<class Teng>
void CalcEnergyHard(const PtclHard ptcl[], const PS::S32 n_tot, Teng & eng, 
                    const PS::F64 r_in, const PS::F64 r_out, const PS::F64 eps_sq = 0.0){
    eng.kin = eng.pot = eng.tot = 0.0;
#ifndef INTEGRATED_CUTOFF_FUNCTION
    PS::F64 r_oi_inv = 1.0/(r_out-r_in);
    PS::F64 r_A = (r_out-r_in)/(r_out+r_in);
    PS::F64 pot_off = cutoff_pot(1.0,r_oi_inv,r_A,r_in)/r_out;
#endif
    for(PS::S32 i=0; i<n_tot; i++){
        eng.kin += 0.5 * ptcl[i].mass * ptcl[i].vel * ptcl[i].vel;

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
}

//const PS::F64 SAFTY_FACTOR_FOR_SEARCH = 1.05;
//const PS::F64 SAFTY_FACTOR_FOR_SEARCH_SQ = SAFTY_FACTOR_FOR_SEARCH * SAFTY_FACTOR_FOR_SEARCH;
//const PS::F64 SAFTY_OFFSET_FOR_SEARCH = 1e-7;

#ifndef NAN_CHECK
#define NAN_CHECK(val) assert((val) == (val));
#endif

void print_p(PtclHard* p, const int n, const int width = 20, const int precision = 14) {
    std::cout<<std::setprecision(precision);
    std::cout<<std::setw(width)<<"mass"
             <<std::setw(width)<<"x1"
             <<std::setw(width)<<"x2"
             <<std::setw(width)<<"x3"
             <<std::setw(width)<<"v1"
             <<std::setw(width)<<"v2"
             <<std::setw(width)<<"v3"
             <<std::setw(width)<<"rsearch"
             <<std::setw(width)<<"mass_bk"
             <<std::setw(width)<<"status"
             <<std::setw(width)<<"id"
             <<std::setw(width)<<"id_cluster"
             <<std::setw(width)<<"adr"
             <<std::endl;
    for (int i=0; i<n; i++) {
        std::cout<<std::setw(width)<<p[i].mass
                 <<std::setw(width)<<p[i].pos[0]
                 <<std::setw(width)<<p[i].pos[1]
                 <<std::setw(width)<<p[i].pos[2]
                 <<std::setw(width)<<p[i].vel[0]
                 <<std::setw(width)<<p[i].vel[1]
                 <<std::setw(width)<<p[i].vel[2]
                 <<std::setw(width)<<p[i].r_search
                 <<std::setw(width)<<p[i].mass_bk
                 <<std::setw(width)<<p[i].status
                 <<std::setw(width)<<p[i].id
                 <<std::setw(width)<<p[i].id_cluster
                 <<std::setw(width)<<p[i].adr_org
                 <<std::endl;
    }
}

// flag: 1: c.m; 2: individual; 
template<class Teng>
void write_p(FILE* fout, const PS::F64 time, const PtclHard* p, const int n, PtclHard & pcm, Teng &et, const PS::F64 rin, const PS::F64 rout, const PS::F64 eps2, const PS::F64 et0=0, const int flag=2) {
    fprintf(fout,"%20.14e ",time);
    PS::ReallocatableArray<PtclHard> pp;
    for (int i=0; i<n; i++) {
        if(flag==2&&(p[i].status>0||p[i].id<0)) continue;
        if(flag==1&&((p[i].id>0&&p[i].status!=0)||(p[i].id<=0&&p[i].status<0))) continue;
        pp.push_back(p[i]);
        if((flag==2||flag==1)&&p[i].status!=0) {
            pp.back().mass = pp.back().mass_bk;
        }
        
    }
    calc_center_of_mass(pcm, pp.getPointer(), pp.size());
    //printf("n = %d\n",pp.size());
    Teng etb;
    CalcEnergyHard(pp.getPointer(),pp.size(),et,rin,rout,eps2);
    PS::F64 err = et0==0?0:(et.tot-et0)/et0;
    fprintf(fout,"%20.14e %20.14e %20.14e %20.14e ",err,et.kin,et.pot,et.tot);
    for (int i=0; i<pp.size(); i++) {
        fprintf(fout,"%20.14e %20.14e %20.14e %20.14e %20.14e %20.14e %20.14e ", 
                pp[i].mass, pp[i].pos[0], pp[i].pos[1], pp[i].pos[2], 
                pp[i].vel[0], pp[i].vel[1], pp[i].vel[2]);
    }
    fprintf(fout,"\n");
}

// flag: 1: c.m; 2: individual; 
template<class Teng>
void write_p(FILE* fout, const PS::F64 time, const PtclHard* p, const int n, const Teng *E, const Teng *E0=NULL) {
    fprintf(fout,"%20.14e ",time);
    PS::ReallocatableArray<PtclHard> pp;
    for (int i=0; i<n; i++) {
        if(p[i].status>0||p[i].id<0) continue;
        pp.push_back(p[i]);
        if(p[i].status!=0) pp.back().mass = pp.back().mass_bk;
    }
    PS::F64 err = E0==NULL?0:(E->tot-E0->tot)/E0->tot;
    fprintf(fout,"%20.14e %20.14e %20.14e %20.14e ",err,E->kin,E->pot,E->tot);
    for (int i=0; i<pp.size(); i++) {
        fprintf(fout,"%20.14e %20.14e %20.14e %20.14e %20.14e %20.14e %20.14e ", 
                pp[i].mass, pp[i].pos[0], pp[i].pos[1], pp[i].pos[2], 
                pp[i].vel[0], pp[i].vel[1], pp[i].vel[2]);
    }
    fprintf(fout,"\n");
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
  PS::F64 rin, rout, rbin, rsearch, eps, eta, dt_limit, time, m_average=0;
  PS::S32 rcount = fscanf(fin, "%lf %d %lf %lf %lf %lf %lf %lf %lf\n", 
                          &time, &N, &rin, &rout, &rsearch, &rbin, &dt_limit, &eta, &eps);
  if (rcount<8) {
      std::cerr<<"Error: parameter reading fail!\n";
      abort();
  }
  Ptcl::r_search_min = rsearch;
  const PS::F64 n_split = 8;

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

  PtclHard pcm;
  calc_center_of_mass(pcm, p.getPointer(), p.size());
  center_of_mass_shift(pcm, p.getPointer(), p.size());

  print_p(p.getPointer(),N);

  SearchGroup<PtclHard> group;
  group.findGroups(p.getPointer(), N, n_split);
  group.searchAndMerge(p.getPointer(), rbin);
  std::cerr<<"SearchAndMerge\n";
  //print_p(p.getPointer(),N);
  
  for(int i=0; i<group.getNumOfGroups(); i++) {
      std::cerr<<"group["<<i<<"]: ";
      for(int j=0; j<group.getGroupN(i); j++) {
          std::cerr<<std::setw(10)<<group.getGroup(i)[j];
      }
      std::cerr<<std::endl;
  }

  std::cerr<<"Ptcl List:";
  for(int i=0; i<group.getPtclN(); i++) {
      std::cerr<<std::setw(10)<<group.getPtclList()[i];
  }
  std::cerr<<std::endl;
  
  PS::ReallocatableArray<PtclHard> ptcl_new;

  group.generateList(p.getPointer(), ptcl_new, rbin, dt_limit, N, n_split);
  std::cerr<<"GenerateList\n";
  print_p(p.getPointer(),p.size());

  std::cerr<<"new ptcl: "<<ptcl_new.size()<<"\n ";
  print_p(ptcl_new.getPointer(),ptcl_new.size());

  p.reserveEmptyAreaAtLeast(ptcl_new.size());
  for(int i=0;i<ptcl_new.size();i++) {
      p.pushBackNoCheck(ptcl_new[i]);
      adr.push_back(i+N);
  }

  np.push_back(p.size());

  std::cerr<<"new p: "<<p.size()<<"\n ";
  print_p(p.getPointer(),np[0]);
    
  for(PS::S32 i=0; i<p.size(); i++) {
      if (p[i].id>0&&p[i].status>0) {
          p[i].vel = PS::F64vec(0.0);
      }
  }


  SystemHard sys;
  PS::ParticleSystem<FPSoft> fp;
  PS::F64 time_sys = 0.0;
  PS::F64 dt_min_hard = 1.0;
  for (PS::S32 i=0;i<40;i++) dt_min_hard *= 0.5;
  sys.setParam(rbin, rout, rin, eps, dt_limit, dt_min_hard, eta, time_sys, 1e-8, N, n_split);
  sys.setARCParam();
  
  sys.setPtclForIsolatedMultiCluster(p,adr,np);
  //sys.initialMultiCluserOMP(fp);

  FILE* fout;
  std::string fname="hard.dat.";
  if ( (fout = fopen((fname+foutname).c_str(),"w")) == NULL) {
    fprintf(stderr,"Error: Cannot open file hard.dat.\n");
    abort();
  }

  FILE* fout2;
  fname = "hardc.dat.";
  if ( (fout2 = fopen((fname+foutname).c_str(),"w")) == NULL) {
    fprintf(stderr,"Error: Cannot open file hardcm.dat.\n");
    abort();
  }

  EnergyAndMomemtum et0,et,etcm0,etcm;
  PS::F64 eps2 = eps*eps;
  PtclHard pcm0,pcm1,ppcm0,ppcm1;
  fprintf(stderr,"Time = %e\n", time_sys);
  print_p(sys.ptcl_hard_.getPointer(),sys.ptcl_hard_.size());
  write_p(fout,time_sys,sys.ptcl_hard_.getPointer(),sys.ptcl_hard_.size(),pcm0,et0,rin,rout,eps2,0.0);
  //write_p(fout,time_sys,sys.ptcl_hard_.getPointer(),sys.ptcl_hard_.size(), &sys.ESD0);
  write_p(fout2,time_sys,sys.ptcl_hard_.getPointer(),sys.ptcl_hard_.size(),ppcm0,etcm0,rin,rout,eps2,0.0,1);
  while(time_sys < time){
      fprintf(stderr,"Time = %e\n", time_sys+dt_limit);
      sys.driveForMultiCluster<PS::ParticleSystem<FPSoft>,FPSoft>(dt_limit, fp);
      time_sys += dt_limit;
      for(PS::S32 i=0; i<sys.ptcl_hard_.size(); i++) {
          if (sys.ptcl_hard_[i].id>0&&sys.ptcl_hard_[i].status>0) {
              sys.ptcl_hard_[i].vel = PS::F64vec(0.0);
          }
      }
      print_p(sys.ptcl_hard_.getPointer(),sys.ptcl_hard_.size());
      //write_p(fout,time_sys,sys.ptcl_hard_.getPointer(),sys.ptcl_hard_.size(), &sys.ESD1, &sys.ESD0);
      write_p(fout,time_sys,sys.ptcl_hard_.getPointer(),sys.ptcl_hard_.size(),pcm1,et,rin,rout,eps2,et0.tot);
      write_p(fout2,time_sys,sys.ptcl_hard_.getPointer(),sys.ptcl_hard_.size(),ppcm1,etcm,rin,rout,eps2,etcm0.tot,1);
      std::cerr<<"CM: pos="<<pcm1.pos<<" vel="<<pcm1.vel<<" shift pos="<<pcm1.pos-pcm0.pos<<" shift vel="<<pcm1.vel-pcm0.vel<<std::endl;
      std::cerr<<"CMHint: pos="<<ppcm1.pos<<" vel="<<ppcm1.vel<<" shift pos="<<ppcm1.pos-ppcm0.pos<<" shift vel="<<ppcm1.vel-ppcm0.vel<<std::endl;
  }
  

  fclose(fin);
  fclose(fout);
  fclose(fout2);
  
  return 0;
}

