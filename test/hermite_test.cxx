#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <particle_simulator.hpp>
#include "Newtonian_acceleration.h"
#include "hard.hpp"
#include "soft.hpp"
#include "integrate.hpp"
#include "cluster_list.hpp"


#ifndef NAN_CHECK
#define NAN_CHECK(val) assert((val) == (val));
#endif

template <class Tptcl, class TEnergy>
void write_p(FILE* fout, const PS::F64 time, const TEnergy &E, const TEnergy &Ediff, const Tptcl* p, const int n) {
    fprintf(fout,"%e %e %e %e %e ",time, Ediff.tot, E.kin, E.pot, E.tot);
    for (int i=0; i<n; i++) {
        fprintf(fout,"%e %e %e %e %e %e %e ", 
                p[i].mass, p[i].pos[0], p[i].pos[1], p[i].pos[2], 
                p[i].vel[0], p[i].vel[1], p[i].vel[2]);
    }
    fprintf(fout,"\n");
}


void print_p(PtclHard* p, const int n) {
    std::cout<<std::setw(12)<<"mass"
             <<std::setw(12)<<"x1"
             <<std::setw(12)<<"x2"
             <<std::setw(12)<<"x3"
             <<std::setw(12)<<"v1"
             <<std::setw(12)<<"v2"
             <<std::setw(12)<<"v3"
             <<std::setw(12)<<"rsearch"
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
                 <<std::setw(12)<<p[i].status
                 <<std::setw(12)<<p[i].id
                 <<std::endl;
    }
}

struct params{
    double rin,rout,r_oi_inv,rsearch,rbin,dt_limit_hard,eta,eps;
};

int main(int argc, char** argv)
{
    // data file name
    char* filename = argv[argc-2];
    char* foutname = argv[argc-1];

    // open data file
    std::fstream fs;
    fs.open(filename,std::fstream::in);
    if(!fs.is_open()) {
        std::cerr<<"Error: Filename "<<filename<<" not found\n";
        abort();
    }

    int N;
    PS::F64 time_end;
    params par;

    EnergyAndMomemtum E0,E1;

    fs>>time_end>>N>>par.rin>>par.rout>>par.rsearch>>par.rbin>>par.dt_limit_hard>>par.eta>>par.eps;
    par.r_oi_inv = 1.0/(par.rout-par.rin);

    PS::F64 dt_min_hard = 1.0;
    for(int i=0;i<40;i++) dt_min_hard *= 0.5;
    
    fprintf(stderr,"t_end = %e\nN = %d\nr_in = %e\nr_out = %e\nr_search = %e\neta = %e\ndt_limit = %e\neps = %e\n",time_end,N,par.rin,par.rout,par.rsearch,par.eta,par.dt_limit_hard,par.eps);

    PS::ReallocatableArray<PtclHard> p;
    p.resizeNoInitialize(N);

    for(int i=0;i<N;i++) {
        fs>>p[i].mass
          >>p[i].pos[0]>>p[i].pos[1]>>p[i].pos[2]
          >>p[i].vel[0]>>p[i].vel[1]>>p[i].vel[2];
        p[i].id = i;
        p[i].status = 0;
        p[i].r_search = par.rsearch;
      
        if (fs.eof()) {
            std::cerr<<"Error: data file reach end when reading pairs (current loaded pair number is "<<i<<"; required pair number "<<N<<std::endl;
            abort();
        }
    }    

    PtclHard pcm;
    calc_center_of_mass(pcm, p.getPointer(), p.size());
    center_of_mass_shift(pcm, p.getPointer(), p.size());
    
    HermiteIntegrator<PtclHard> Hint;
    Hint.setParams(par.eta, par.rin, par.rout, par.eps*par.eps, p.size());
    Hint.reserveMem(p.size());
    Hint.addPtclList(p.getPointer(), NULL, p.size(), 0, 0.0, false);
    Hint.searchPerturber(p.size());
            
    Hint.CalcEnergy(E0);
    
    std::cerr<<"Energy: init: ";
    E0.print(std::cerr);

    //PS::S32 group_act_n = 0;
    //PS::ReallocatableArray<PS::S32> group_act_list; //active group_list act adr

    //group_act_list.resizeNoInitialize(p.size());
            
    PS::F64 time_sys=0.0;
#ifdef FIX_STEP_DEBUG
    PS::F64 dt_limit = par.dt_limit_hard;
#else
    PS::F64 dt_limit = calcDtLimit(time_sys, par.dt_limit_hard, dt_min_hard);
#endif
    ARCIntegrator<PtclHard, PtclH4, PtclForce>* arcint = NULL;
    Hint.calcA0offset();
    Hint.initial(NULL, p.size(), 0.0, dt_limit, dt_min_hard, arcint, false);
    Hint.SortAndSelectIp();

    FILE* fout;
    std::string fname="hermite.dat.";
    if ( (fout = fopen((fname+foutname).c_str(),"w")) == NULL) {
        fprintf(stderr,"Error: Cannot open file hard.dat.\n");
        abort();
    }
    
    PS::F64 time_pre=0.0;
    EnergyAndMomemtum Edif;
    PtclH4 pcm0,pcm1;
    calc_center_of_mass(pcm0, Hint.getPtcl(), Hint.getPtclN());
    write_p(fout,time_sys,E0,Edif,Hint.getPtcl(),Hint.getPtclN());
    while(time_sys<time_end) {
        time_pre = time_sys;
        time_sys = Hint.getNextTime();
        assert(time_sys>time_pre);
#ifdef FIX_STEP_DEBUG
        dt_limit = par.dt_limit_hard;
#else
        dt_limit = calcDtLimit(time_sys, par.dt_limit_hard, dt_min_hard);
#endif
        Hint.integrateOneStepAct(time_sys,dt_limit,dt_min_hard,arcint);
        Hint.SortAndSelectIp();
        if(fmod(time_sys,par.dt_limit_hard)==0) {
            //std::cout<<"Time = "<<time_sys<<std::endl;
            Hint.CalcEnergy(E1);
            Edif = E1-E0;
            calc_center_of_mass(pcm1, Hint.getPtcl(), Hint.getPtclN());
            write_p(fout,time_sys,E1,Edif,Hint.getPtcl(),Hint.getPtclN());
            std::cerr<<"CM: pos="<<pcm1.pos<<" vel="<<pcm1.vel<<" shift pos="<<pcm1.pos-pcm0.pos<<" shift vel="<<pcm1.vel-pcm0.vel<<std::endl;
            
        }
    }
    Hint.writeBackPtcl(0);
    
    Hint.CalcEnergy(E1);

    EnergyAndMomemtum Ediff = E1-E0;
    
    std::cerr<<"final: ";
    E1.print(std::cerr);
    std::cerr<<"error: "<<Ediff.tot/E0.tot<<std::endl;

#ifdef HARD_DEBUG
    Hint.printStepHist();
#endif
    
    return 0;

}
