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

PtclHard kepler_print(const std::size_t id, const std::size_t ib, PtclHard* c[2], params& ppars){
    const double* x[2];
    const double* v[2];
    double m[2];
    for (std::size_t i=0; i<2; i++) {
      x[i]=c[i]->getPos();
      v[i]=c[i]->getVel();
      m[i]=c[i]->getMass();
    }
        
    double ax,per,ecc,angle[3],true_anomaly,ecc_anomaly,mean_anomaly; 
    double dx[3] = {x[1][0]-x[0][0], x[1][1]-x[0][1], x[1][2]-x[0][2]};
    double dv[3] = {v[1][0]-v[0][0], v[1][1]-v[0][1], v[1][2]-v[0][2]};
    double mt = m[0]+m[1];
    
    NTA::calc_kepler_orbit_par(ax,per,ecc,angle,true_anomaly,ecc_anomaly,mean_anomaly,mt,dx,dv);
    std::cout<<std::setw(12)<<id
             <<std::setw(12)<<ib
             <<std::setw(12)<<ax
             <<std::setw(12)<<ecc
             <<std::setw(12)<<per
             <<std::setw(12)<<angle[0]
             <<std::setw(12)<<angle[1]
             <<std::setw(12)<<angle[2]
             <<std::setw(12)<<ecc_anomaly
             <<std::setw(12)<<true_anomaly
             <<std::setw(12)<<mean_anomaly
             <<std::endl;
    if(id==0&&ib==0) 
        std::cout<<"      " 
                 <<std::setw(12)<<"id          " 
                 <<std::setw(12)<<"ib          " 
                 <<std::setw(12)<<"ax          " 
                 <<std::setw(12)<<"ecc         " 
                 <<std::setw(12)<<"per         " 
                 <<std::setw(12)<<"angle[0]    " 
                 <<std::setw(12)<<"angle[1]    " 
                 <<std::setw(12)<<"angle[2]    " 
                 <<std::setw(12)<<"ecc_anomaly " 
                 <<std::setw(12)<<"true_anomaly" 
                 <<std::setw(12)<<"mean_anomaly" 
                 <<std::endl;                 
            

    PS::F64vec xcm((x[0][0]*m[0]+x[1][0]*m[1])/mt, 
                   (x[0][1]*m[0]+x[1][1]*m[1])/mt, 
                   (x[0][2]*m[0]+x[1][2]*m[1])/mt);

    PS::F64vec vcm((v[0][0]*m[0]+v[1][0]*m[1])/mt, 
                   (v[0][1]*m[0]+v[1][1]*m[1])/mt, 
                   (v[0][2]*m[0]+v[1][2]*m[1])/mt);

    return PtclHard(ib,mt,xcm,vcm,ppars.rout,0,0);
}


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
  plist.pair_process(0,0,kepler_print,par);

  updateRout(p, N, par.rin, par.rout, par.gmin, m_average);
  for (int i=0; i<N; i++) {
    std::cout<<i<<" "<<p[i]->r_out<<std::endl;
  }

  return 0;
}
