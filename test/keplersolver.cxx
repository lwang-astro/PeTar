#include <iostream>
#include <fstream>
#include <iomanip>
#include <particle_simulator.hpp>
#include "Newtonian_acceleration.h"
#include "usr_define.hpp"


#include "kepler.hpp"


// center of mass shift
template <class particle>
void center_of_mass_shift(particle &cm, particle p[], const int num) {
  for (int i=0;i<num;i++) {
    p[i].pos -= cm.pos;
    p[i].vel -= cm.vel;
  }
}

void print_p(ParticleBase* p, const int n) {
    std::cout<<std::setw(12)<<"mass"
             <<std::setw(12)<<"x1"
             <<std::setw(12)<<"x2"
             <<std::setw(12)<<"x3"
             <<std::setw(12)<<"v1"
             <<std::setw(12)<<"v2"
             <<std::setw(12)<<"v3"
             <<std::endl;
    for (int i=0; i<n; i++) {
        std::cout<<std::setw(12)<<p[i].mass
                 <<std::setw(12)<<p[i].pos[0]
                 <<std::setw(12)<<p[i].pos[1]
                 <<std::setw(12)<<p[i].pos[2]
                 <<std::setw(12)<<p[i].vel[0]
                 <<std::setw(12)<<p[i].vel[1]
                 <<std::setw(12)<<p[i].vel[2]
                 <<std::endl;
    }
}

ParticleBase kepler_print(const std::size_t id, const std::size_t ib, ParticleBase c[2]){
    PS::F64vec x[2],v[2];
    double m[2];
    for (std::size_t i=0; i<2; i++) {
      x[i]=c[i].pos;
      v[i]=c[i].vel;
      m[i]=c[i].mass;
    }
        
    double ax,per,ecc,angle[3],true_anomaly,ecc_anomaly,mean_anomaly; 
    double dx[3] = {x[1][0]-x[0][0], x[1][1]-x[0][1], x[1][2]-x[0][2]};
    double dv[3] = {v[1][0]-v[0][0], v[1][1]-v[0][1], v[1][2]-v[0][2]};
    double mt = m[0]+m[1];
    
    NTA::calc_kepler_orbit_par(ax,per,ecc,angle,true_anomaly,ecc_anomaly,mean_anomaly,mt,dx,dv);
    std::cout<<std::setw(12)<<id
             <<std::setw(12)<<ib
             <<std::setw(12)<<m[0]
             <<std::setw(12)<<m[1]
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
                 <<std::setw(12)<<"m[0]        " 
                 <<std::setw(12)<<"m[1]        " 
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

    return ParticleBase(mt, xcm, vcm);
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
  ParticleBase p[2];

  for(int i=0;i<2;i++) {
      double m,x,y,z,vx,vy,vz;
      fs>>m>>x>>y>>z>>vx>>vy>>vz;
      if (fs.eof()) {
          std::cerr<<"Error: data file reach end when reading pairs (current loaded pair number is "<<i<<"; required pair number "<<2<<std::endl;
          abort();
      }
      PS::F64vec r(x,y,z);
      PS::F64vec v(vx,vy,vz);
      p[i]=ParticleBase(m,r,v);
  }
  
  kepler_print(0,0,p);
  ParticleBase pcm;
  calc_center_of_mass(pcm,p,2);

  double ax, ecc, inc, OMG, omg, tperi, ecca, peri;
  ecca=PosVel2OrbParam(ax,ecc,inc,OMG,omg,tperi,peri,p[0].pos,p[1].pos,p[0].vel,p[1].vel,p[0].mass,p[1].mass);
  //peri = 8.0*std::atan(1.0)*std::abs(ax)*std::sqrt(std::abs(ax)/pcm.mass);

  printf("From P2P fun: ax=%e, ecc=%e, inc=%e, OMG=%e, omg=%e, tperi=%e, peri=%e, ecca=%e\n",ax, ecc, inc, OMG, omg, tperi, peri, ecca);
  center_of_mass_shift(pcm,p,2);
  
  printf("original:\n");
  print_p(p,2);

  double dt = peri/8.0;
  for(int i=0; i<8; i++) {
      DriveKeplerOrbParam(p[0].pos, p[1].pos, p[0].vel, p[1].vel, p[0].mass, p[1].mass, (i+1)*dt, ax, ecc, inc, OMG, omg, peri, ecca);      
      printf("peri = %e, t = %e\n",peri,(i+1)*dt);
      print_p(p,2);
      kepler_print(0,0,p);
  }

  return 0;
}
