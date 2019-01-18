#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <particle_simulator.hpp>
#include "integrate.hpp"
#include "hard.hpp"

void print_ptcl(PtclHard* p, const int n) {
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

struct Psoft{
    PS::F64vec acc;
    PS::F64 mass_bk;
};

int main(int argc, char** argv)
{
    // data file name
    char* filename = argv[argc-1];

    FILE* fin;
    if ( (fin = fopen(filename,"r")) == NULL) {
        fprintf(stderr,"Error: Cannot open input file %s.\n",filename);
        abort();
    }

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

    int check_flag; // 1: break arc; 2: break hermite
    fs>>time_end>>N>>par.rin>>par.rout>>par.rsearch>>par.rbin>>par.dt_limit_hard>>par.eta>>par.eps>>check_flag;
    par.r_oi_inv = 1.0/(par.rout-par.rin);

    fprintf(stderr,"t_end = %e\nN = %d\nr_in = %e\nr_out = %e\nr_search = %e\neta = %e\ndt_limit = %e\neps = %e\n",time_end,N,par.rin,par.rout,par.rsearch,par.eta,par.dt_limit_hard,par.eps);

    PS::ReallocatableArray<PtclHard> p;
    p.resizeNoInitialize(N);

    PS::F64 mave=0.0;
    for(int i=0;i<N;i++) {
        fs>>p[i].mass
          >>p[i].pos[0]>>p[i].pos[1]>>p[i].pos[2]
          >>p[i].vel[0]>>p[i].vel[1]>>p[i].vel[2];
        p[i].id = i;
        p[i].status = 0;
        p[i].r_search = par.rsearch;
        mave += p[i].mass;
      
        if (fs.eof()) {
            std::cerr<<"Error: data file reach end when reading particle "<<i<<"; required particle number "<<N<<std::endl;
            abort();
        }
    }    
    mave /= N;

    Ptcl::search_factor = 3.0;
    Ptcl::r_search_min = par.rsearch;
    Ptcl::mean_mass_inv = 1.0/mave;

    print_ptcl(p.getPointer(),N);

    ARC_int_pars Int_pars;
    Int_pars.rin = par.rin;
    Int_pars.rout = par.rout;
    Int_pars.r_oi_inv = 1.0/(par.rout-par.rin);
    Int_pars.r_A = (par.rout-par.rin)/(par.rout+par.rin);
    Int_pars.pot_off = (1.0 + Int_pars.r_A)/par.rout;
    Int_pars.eps2 = par.eps*par.eps;

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

    HermiteIntegrator <PtclHard> Hint;
    ARCIntegrator<Ptcl, PtclH4, PtclForce> Aint(ARC_control,Int_pars);

    Hint.setParams(par.eta, par.rin, par.rout, par.eps*par.eps, p.size());
    Hint.reserveMem(p.size());

    Aint.reserveARMem(p.size());
    Aint.reservePertMem(p.size(),p.size());

    //PS::ReallocatableArray<PS::S32> group_member_index[p.size()];  // Group member index in _ptcl_local array
    //PS::S32 single_index[p.size()]; 
    //PS::S32 n_group=0; // number of groups, including blend groups

    PS::F64 dt_min_hard = std::pow(0.5,32);

    if(check_flag==2) {
        Hint.addPtclList(p.getPointer(), NULL, p.size(), 0, 0.0, false);
        Hint.searchPerturber(p.size());
    
    
        Hint.calcA0offset();
        Hint.initial(NULL, p.size(), 0.0, par.dt_limit_hard, dt_min_hard, &Aint, false);
        Hint.SortAndSelectIp();
    }
    else {
        Aint.addOneGroup(p.getPointer(), NULL, p.size(), (Psoft*)NULL, 0, Hint.getPtcl(), Hint.getForce());
        Aint.initialOneChain(0);
        Aint.initialOneSys(0,0.0);
        //group_member_index[0].resizeNoInitialize(p.size());
        //for(int i=0; i<p.size(); i++) group_member_index[0][i]=i;
//        n_group = Aint.getNGroups();

        Hint.addPtclList((PtclHard*)Aint.getCM(0), NULL, 1, 0, 0.0, false);
        //n_single = 0;
    }

    SystemHard sys;
    sys.setParam(par.rbin, par.rout, par.rin, par.eps, par.dt_limit_hard, dt_min_hard, par.eta, 0.0, 1e-8, N+10, 8);
    sys.setARCParam(1e-8, 1e-6, 1e-6, dt_min_hard);

    sys.adjustGroup<Psoft>(Hint, Aint, p.getPointer(), 0.0, par.dt_limit_hard, time_end, par.rbin);
    
    return 0;
}
