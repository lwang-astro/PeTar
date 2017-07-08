#include <iostream>
#include <fstream>
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


PS::F64 calcDtLimit(const PS::F64 time_sys,
                    const PS::F64 dt_limit_org,
                    const PS::F64 time_offset = 0.0){
    PS::F64 dt_limit_ret = dt_limit_org;
    PS::F64 s = (time_sys-time_offset) / dt_limit_ret;
    PS::F64 time_head = ((PS::S64)(s)) * dt_limit_org;
    PS::F64 time_shifted = time_sys - time_head;
    while( fmod(time_shifted, dt_limit_ret) != 0.0) dt_limit_ret *= 0.5;
    return dt_limit_ret;
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
    double rin,rout,dt_limit_hard,eta,eps;
};

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
    PS::F64 time_end;
    params par;

    Energy E0,E1;

    fs>>time_end>>N>>par.rin>>par.rout>>par.dt_limit_hard>>par.eta>>par.eps;

    fprintf(stderr,"t_end = %e\nN = %d\nr_in = %e\nr_out = %e\neta = %e\ndt_limit = %e\neps = %e\n",time_end,N,par.rin,par.rout,par.eta,par.dt_limit_hard,par.eps);

    PS::ReallocatableArray<PtclHard> p;
    p.resizeNoInitialize(N);

    for(int i=0;i<N;i++) {
        fs>>p[i].mass
          >>p[i].pos[0]>>p[i].pos[1]>>p[i].pos[2]
          >>p[i].vel[0]>>p[i].vel[1]>>p[i].vel[2];
        p[i].id = i;
        p[i].status = 0;
        p[i].r_search = par.rout*10;
      
        if (fs.eof()) {
            std::cerr<<"Error: data file reach end when reading pairs (current loaded pair number is "<<i<<"; required pair number "<<N-1<<std::endl;
            abort();
        }
    }    

    SearchGroup<PtclHard> group;
            
    group.findGroups(p.getPointer(), N);
    
    HermiteIntegrator<PtclHard> Hint;
    Hint.setParams(par.dt_limit_hard, par.eta, par.rin, par.rout, par.eps*par.eps);
    Hint.setPtcl(p.getPointer(),p.size(),group.getPtclList(),group.getNPtcl());
    Hint.searchPerturber();
            
    Hint.CalcEnergyHard(E0);
    
    std::cerr<<"Energy: init: ";
    E0.dump(std::cerr);

    PS::S32 group_act_n = 0;
    PS::ReallocatableArray<PS::S32> group_act_list; //active group_list act adr

    group_act_list.resizeNoInitialize(group.getNPtcl());
            
    PS::S32 n_groups = 0;
    PS::F64 time_sys=0.0;
    PS::F64 dt_limit = calcDtLimit(time_sys, par.dt_limit_hard);
    ARCIntegrator<PtclHard, PtclH4, PtclForce, ARC_int_pars, ARC_pert_pars>* arcint = NULL;
    Hint.initialize(dt_limit, group_act_list.getPointer(), group_act_n, n_groups, arcint);

    FILE* fout;
    if ( (fout = fopen("hermite.dat","w")) == NULL) {
        fprintf(stderr,"Error: Cannot open file hard.dat.\n");
        abort();
    }
    
    PS::F64 time_pre=0.0;
    Energy Edif = E0.calcDiff(E0);
    PtclH4 pcm0,pcm1;
    calc_center_of_mass(pcm0, Hint.getPtcl(), Hint.getPtclN());
    write_p(fout,time_sys,E0,Edif,Hint.getPtcl(),Hint.getPtclN());
    while(time_sys<time_end) {
        time_pre = time_sys;
        time_sys = Hint.getNextTime();
        assert(time_sys>time_pre);
        dt_limit = calcDtLimit(time_sys, par.dt_limit_hard);
        Hint.integrateOneStep(time_sys,dt_limit,true,arcint);
        Hint.SortAndSelectIp(group_act_list.getPointer(), group_act_n, n_groups);
        if(fmod(time_sys,par.dt_limit_hard)==0) {
            //std::cout<<"Time = "<<time_sys<<std::endl;
            Hint.CalcEnergyHard(E1);
            Edif = E1.calcDiff(E0);
            calc_center_of_mass(pcm1, Hint.getPtcl(), Hint.getPtclN());
            write_p(fout,time_sys,E1,Edif,Hint.getPtcl(),Hint.getPtclN());
            std::cerr<<"CM: pos="<<pcm1.pos<<" vel="<<pcm1.vel<<" shift pos="<<pcm1.pos-pcm0.pos<<" shift vel="<<pcm1.vel-pcm0.vel<<std::endl;
            
        }
    }
    Hint.writeBackPtcl(p.getPointer(),p.size(),group.getPtclList(),group.getNPtcl());
    
    Hint.CalcEnergyHard(E1);

    Energy Ediff = E1.calcDiff(E0);
    
    std::cerr<<"final: ";
    E1.dump(std::cerr);
    std::cerr<<"error: "<<Ediff.tot/E0.tot<<std::endl;

#ifdef HARD_DEBUG
    Hint.printStepHist();
#endif
    
    return 0;

}
