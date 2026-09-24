// tsv110 kernel benchmark driver: PeTar NoSimd vs NEON prototypes
#include <particle_simulator.hpp>
#include "soft_ptcl.hpp"
#include "soft_force.hpp"
#include "static_variables.hpp"
#ifdef USE_NEON_KERNEL
#include "force_tsv110.hpp"
#endif
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <algorithm>
#include <cmath>

using SPJQuad = PS::SPJQuadrupoleInAndOut;
using SJPMono = PS::SPJMonopoleInAndOut;

static double now_s(){
    using namespace std::chrono;
    return duration_cast<duration<double>>(steady_clock::now().time_since_epoch()).count();
}

struct Dataset {
    std::vector<EPISoft> epi;
    std::vector<EPJSoft> epj;
    std::vector<SPJQuad>  spq;
    std::vector<SJPMono>  spm;
};

static void gen(Dataset& D, int ni, int nj, int ns){
    unsigned s=12345u;
    auto rnd=[&](){ s=s*1664525u+1013904223u; return (double)s/4294967296.0; };
    D.epi.assign(ni, EPISoft());
    D.epj.assign(nj, EPJSoft());
    D.spq.assign(ns, SPJQuad());
    D.spm.assign(ns, SJPMono());
    for(int i=0;i<ni;i++){
        D.epi[i].id=i;
        D.epi[i].pos=PS::F64vec(0.1+0.8*rnd(), 0.1+0.8*rnd(), 0.1+0.8*rnd());
        D.epi[i].r_search=0.01+0.04*rnd();
        D.epi[i].type=1;
        D.epi[i].rank_org=0;
    }
    for(int j=0;j<nj;j++){
        D.epj[j].id=j;
        D.epj[j].mass=1e-3*(0.5+rnd());
        D.epj[j].pos=PS::F64vec(0.1+0.8*rnd(), 0.1+0.8*rnd(), 0.1+0.8*rnd());
        D.epj[j].r_search=0.01+0.04*rnd();
    }
    for(int j=0;j<ns;j++){
        double m=1e-3*(0.5+rnd());
        PS::F64vec p(0.1+0.8*rnd(), 0.1+0.8*rnd(), 0.1+0.8*rnd());
        D.spq[j].mass=m; D.spq[j].pos=p;
        D.spq[j].quad.xx=1e-7*rnd(); D.spq[j].quad.yy=1e-7*rnd(); D.spq[j].quad.zz=1e-7*rnd();
        D.spq[j].quad.xy=1e-7*rnd(); D.spq[j].quad.xz=1e-7*rnd(); D.spq[j].quad.yz=1e-7*rnd();
        D.spm[j].mass=m; D.spm[j].pos=p;
    }
}

static std::vector<ForceSoft> clear_force(int n){
    std::vector<ForceSoft> f(n);
    for(int i=0;i<n;i++) f[i].clear();
    return f;
}

template<class Call>
static double bench_call(Call c, double target, long* reps_used){
    c(); // warm up
    double t1s=now_s(); c(); double t1=now_s()-t1s;
    long reps = (long)(target/std::max(t1,1e-12));
    if(reps<3) reps=3;
    if(reps>2000000) reps=2000000;
    double best=1e300;
    for(int r=0;r<3;r++){
        double a=now_s();
        for(long k=0;k<reps;k++) c();
        double b=now_s()-a;
        if(b<best) best=b;
    }
    if(reps_used) *reps_used=reps;
    return best/reps; // seconds per call
}

/* ------------- correctness ------------- */
static void compare(const std::vector<ForceSoft>& ref, const std::vector<ForceSoft>& got,
                    const char* name, double* maxrel, double* maxpot, int* nnmis){
    double mr=0, mp=0; int nm=0;
    for(size_t i=0;i<ref.size();i++){
        double a0=std::sqrt(ref[i].acc.x*ref[i].acc.x+ref[i].acc.y*ref[i].acc.y+ref[i].acc.z*ref[i].acc.z);
        double a1=std::sqrt(got[i].acc.x*got[i].acc.x+got[i].acc.y*got[i].acc.y+got[i].acc.z*got[i].acc.z);
        double denom=std::max(a0,1e-300);
        double d=std::fabs(a1-a0)/denom;
        if(d>mr) mr=d;
        double p0=ref[i].pot, p1=got[i].pot;
        double dp=std::fabs(p1-p0)/std::max(std::fabs(p0),1e-300);
        if(dp>mp) mp=dp;
        if(ref[i].n_ngb!=got[i].n_ngb) nm++;
    }
    if(maxrel) *maxrel=mr;
    if(maxpot) *maxpot=mp;
    if(nnmis) *nnmis=nm;
    printf("# %-28s max_rel_acc=%.3e max_rel_pot=%.3e nnb_mismatch=%d\n", name, mr, mp, nm);
}

int main(int argc, char** argv){
    EPISoft::eps=1e-4;
    EPISoft::r_out=0.01;
    ForceSoft::grav_const=1.0;

    const char* mode = argc>1?argv[1]:"check";

    if(!strcmp(mode,"check")){
        const int ni=500, nj=2000, ns=1000;
        Dataset D; gen(D,ni,nj,ns);
        auto f0=clear_force(ni); auto f1=clear_force(ni); auto f2=clear_force(ni);
        SearchNeighborEpEpNoSimd nb0;
        nb0(D.epi.data(),ni,D.epj.data(),nj,f0.data());
        CalcForceEpEpWithLinearCutoffNoSimd ee0;
        ee0(D.epi.data(),ni,D.epj.data(),nj,f0.data());
        CalcForceEpSpQuadNoSimd q0;
        q0(D.epi.data(),ni,D.spq.data(),ns,f0.data());
#ifdef USE_NEON_KERNEL
        tsv110::SearchNeighborEpEpNeon nb1;
        nb1(D.epi.data(),ni,D.epj.data(),nj,f1.data());
        tsv110::CalcForceEpEpWithLinearCutoffNeon ee1(EPISoft::eps*EPISoft::eps, EPISoft::r_out*EPISoft::r_out, ForceSoft::grav_const);
        ee1(D.epi.data(),ni,D.epj.data(),nj,f1.data());
        tsv110::CalcForceEpSpQuadNeon<SPJQuad> q1(EPISoft::eps*EPISoft::eps, ForceSoft::grav_const);
        q1(D.epi.data(),ni,D.spq.data(),ns,f1.data());
        double a,p; int nm;
        compare(f0,f1,"NEON-vs-NoSimd EPEP+SPquad+NB", &a,&p,&nm);
        // isolated kernels
        auto g0=clear_force(ni); auto g1=clear_force(ni);
        nb0(D.epi.data(),ni,D.epj.data(),nj,g0.data());
        nb1(D.epi.data(),ni,D.epj.data(),nj,g1.data());
        compare(g0,g1,"NEON-vs-NoSimd NB only", &a,&p,&nm);
        g0=clear_force(ni); g1=clear_force(ni);
        ee0(D.epi.data(),ni,D.epj.data(),nj,g0.data());
        ee1(D.epi.data(),ni,D.epj.data(),nj,g1.data());
        compare(g0,g1,"NEON-vs-NoSimd EP-EP only", &a,&p,&nm);
        g0=clear_force(ni); g1=clear_force(ni);
        q0(D.epi.data(),ni,D.spq.data(),ns,g0.data());
        q1(D.epi.data(),ni,D.spq.data(),ns,g1.data());
        compare(g0,g1,"NEON-vs-NoSimd SPquad only", &a,&p,&nm);
#else
        printf("# built without USE_NEON_KERNEL\n");
#endif
        return 0;
    }

    if(!strcmp(mode,"scan")){
        int ni_list[] = {4,16,64,256,1024};
        int nj_list[] = {8,32,128,512,2048};
        const char* ksel = argc>2?argv[2]:"all";
        printf("# kernel,variant,ni,nj,time_ns_per_call,interactions_per_s\n");
        for(int a=0;a<5;a++) for(int b=0;b<5;b++){
            int ni=ni_list[a], nj=nj_list[b];
            Dataset D; gen(D,ni,std::max(nj,ni),std::max(nj,ni));
            const double target=0.15;
            auto run=[&](const char* kname, const char* vname, auto call){
                long reps;
                double t=bench_call(call,target,&reps);
                printf("%s,%s,%d,%d,%.1f,%.3e\n",kname,vname,ni,nj,t*1e9,(double)ni*nj/t);
                fflush(stdout);
            };
            bool all = !strcmp(ksel,"all");
            if(all || !strcmp(ksel,"nb")){
                SearchNeighborEpEpNoSimd nb; auto f=clear_force(ni);
                run("nb","nosimd",[&](){ nb(D.epi.data(),ni,D.epj.data(),nj,f.data()); });
#ifdef USE_NEON_KERNEL
                tsv110::SearchNeighborEpEpNeon nbN;
                run("nb","neon4",[&](){ nbN.Kernel_I4_J1(D.epi.data(),ni,D.epj.data(),nj,f.data()); });
                run("nb","neon1",[&](){ nbN.Kernel_I1_J4(D.epi.data(),ni,D.epj.data(),nj,f.data()); });
#endif
            }
            if(all || !strcmp(ksel,"epep")){
                CalcForceEpEpWithLinearCutoffNoSimd ee; auto f=clear_force(ni);
                run("epep","nosimd",[&](){ ee(D.epi.data(),ni,D.epj.data(),nj,f.data()); });
#ifdef USE_NEON_KERNEL
                tsv110::CalcForceEpEpWithLinearCutoffNeon eeN(EPISoft::eps*EPISoft::eps, EPISoft::r_out*EPISoft::r_out, ForceSoft::grav_const);
                run("epep","neon4",[&](){ eeN.Kernel_I4_J1(D.epi.data(),ni,D.epj.data(),nj,f.data()); });
                run("epep","neon1",[&](){ eeN.Kernel_I1_J4(D.epi.data(),ni,D.epj.data(),nj,f.data()); });
#endif
            }
            if(all || !strcmp(ksel,"mono")){
                CalcForceEpSpMonoNoSimd m; auto f=clear_force(ni);
                run("mono","nosimd",[&](){ m(D.epi.data(),ni,D.spm.data(),nj,f.data()); });
#ifdef USE_NEON_KERNEL
                tsv110::CalcForceEpSpMonoNeon<SJPMono> mN(EPISoft::eps*EPISoft::eps, ForceSoft::grav_const);
                run("mono","neon4",[&](){ mN.Kernel_I4_J1(D.epi.data(),ni,D.spm.data(),nj,f.data()); });
                run("mono","neon1",[&](){ mN.Kernel_I1_J4(D.epi.data(),ni,D.spm.data(),nj,f.data()); });
#endif
            }
            if(all || !strcmp(ksel,"quad")){
                CalcForceEpSpQuadNoSimd q; auto f=clear_force(ni);
                run("quad","nosimd",[&](){ q(D.epi.data(),ni,D.spq.data(),nj,f.data()); });
#ifdef USE_NEON_KERNEL
                tsv110::CalcForceEpSpQuadNeon<SPJQuad> qN(EPISoft::eps*EPISoft::eps, ForceSoft::grav_const);
                run("quad","neon4",[&](){ qN.Kernel_I4_J1(D.epi.data(),ni,D.spq.data(),nj,f.data()); });
                run("quad","neon1",[&](){ qN.Kernel_I1_J4(D.epi.data(),ni,D.spq.data(),nj,f.data()); });
#endif
            }
        }
        return 0;
    }

    if(!strcmp(mode,"one")){
        // one <kernel> <variant> <ni> <nj>
        const char* k = argv[2]; const char* v = argv[3];
        int ni=atoi(argv[4]), nj=atoi(argv[5]);
        Dataset D; gen(D,ni,std::max(nj,ni),std::max(nj,ni));
        auto f=clear_force(ni);
        long reps; double t=0;
        if(!strcmp(k,"epep") && !strcmp(v,"nosimd")){
            CalcForceEpEpWithLinearCutoffNoSimd ee; t=bench_call([&](){ee(D.epi.data(),ni,D.epj.data(),nj,f.data());},0.2,&reps);
        }
        printf("%s,%s,%d,%d,%.1f,%.3e\n",k,v,ni,nj,t*1e9,(double)ni*nj/t);
        return 0;
    }

    fprintf(stderr,"usage: %s check|scan [kernel]\n",argv[0]);
    return 1;
}
