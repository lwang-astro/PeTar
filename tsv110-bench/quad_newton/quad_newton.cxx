/* quad_newton.cxx - kernel-level comparison of the NEON quadrupole kernel with and
 * without the extra half Newton step (NEON_QUAD_NEWTON), against the NoSimd F64
 * reference.  Also checks the EP-EP mass>0 filtering introduced for parity with the
 * x86 SIMD / Fugaku kernels.
 *
 * usage:
 *   ./quad_newton err     <ni> <ns> <scale> <offset> <seed>
 *   ./quad_newton errdump <ni> <ns> <scale> <offset> <seed>
 *   ./quad_newton time    <ni> <ns>
 *   ./quad_newton masszero <ni> <nj> <seed>
 */
#include <particle_simulator.hpp>
#include "soft_ptcl.hpp"
#include "soft_force.hpp"
#include "static_variables.hpp"
#include "force_tsv110.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>
#include <chrono>

using SPJQuad = PS::SPJQuadrupoleInAndOut;

#if NEON_QUAD_NEWTON
static const char* NEON_VARIANT = "neon_newton";
#else
static const char* NEON_VARIANT = "neon_cubic";
#endif

static double now_s(){
    using namespace std::chrono;
    return duration_cast<duration<double>>(steady_clock::now().time_since_epoch()).count();
}
static double norm3(const PS::F64vec& v){ return std::sqrt(v.x*v.x+v.y*v.y+v.z*v.z); }

struct Dataset {
    std::vector<EPISoft> epi;
    std::vector<SPJQuad> sp;
};

static void gen(Dataset& D, int ni, int ns, double scale, double offset, unsigned seed){
    D.epi.assign(ni, EPISoft());
    D.sp.assign(ns, SPJQuad());
    auto rnd=[&](){ seed=seed*1664525u+1013904223u; return (double)seed/4294967296.0; };
    for(int i=0;i<ni;i++){
        D.epi[i].id=i;
        D.epi[i].pos=PS::F64vec(offset+scale*rnd(), offset+scale*rnd(), offset+scale*rnd());
        D.epi[i].r_search=0.05*scale;
        D.epi[i].type=1;
        D.epi[i].rank_org=0;
    }
    const double qsc = 0.3*scale;
    for(int j=0;j<ns;j++){
        D.sp[j].mass=1e-3*(0.5+rnd());
        D.sp[j].pos=PS::F64vec(offset+scale*rnd(), offset+scale*rnd(), offset+scale*rnd());
        D.sp[j].quad.xx=D.sp[j].mass*qsc*qsc*(rnd()-0.5);
        D.sp[j].quad.yy=D.sp[j].mass*qsc*qsc*(rnd()-0.5);
        D.sp[j].quad.zz=D.sp[j].mass*qsc*qsc*(rnd()-0.5);
        D.sp[j].quad.xy=D.sp[j].mass*qsc*qsc*(rnd()-0.5);
        D.sp[j].quad.xz=D.sp[j].mass*qsc*qsc*(rnd()-0.5);
        D.sp[j].quad.yz=D.sp[j].mass*qsc*qsc*(rnd()-0.5);
    }
}

static std::vector<ForceSoft> clear_force(int n){
    std::vector<ForceSoft> f(n);
    for(int i=0;i<n;i++) f[i].clear();
    return f;
}

struct Stats { double mx, mean, rms, p50, p90, p99; };
static Stats stat_of(std::vector<double>& v){
    Stats s; s.mx=0; s.mean=0; s.rms=0; s.p50=s.p90=s.p99=0;
    if(v.empty()) return s;
    double sum=0,sumsq=0;
    for(double x:v){ sum+=x; sumsq+=x*x; if(x>s.mx) s.mx=x; }
    std::sort(v.begin(),v.end());
    s.mean=sum/v.size(); s.rms=std::sqrt(sumsq/v.size());
    s.p50=v[(size_t)(v.size()*0.50)];
    s.p90=v[(size_t)(std::min(v.size()-1,(size_t)(v.size()*0.90)))];
    s.p99=v[(size_t)(std::min(v.size()-1,(size_t)(v.size()*0.99)))];
    return s;
}

static void error_vectors(const std::vector<ForceSoft>& ref, const std::vector<ForceSoft>& got,
                          std::vector<double>& ea, std::vector<double>& ep){
    ea.clear(); ep.clear();
    for(size_t i=0;i<ref.size();i++){
        double a0=norm3(ref[i].acc), a1=norm3(got[i].acc);
        ea.push_back(std::fabs(a1-a0)/std::max(a0,1e-300));
        ep.push_back(std::fabs(got[i].pot-ref[i].pot)/std::max(std::fabs(ref[i].pot),1e-300));
    }
}

template<class F>
static double bench_call(F f){
    const double target=0.15;
    f(); double t1=now_s(); f(); double t0=now_s(); double dt=t0-t1;
    long reps=(long)(target/std::max(dt,1e-12));
    if(reps<3) reps=3;
    if(reps>2000000) reps=2000000;
    double best=1e300;
    for(int r=0;r<3;r++){
        double a=now_s();
        for(long k=0;k<reps;k++) f();
        double b=now_s()-a;
        if(b<best) best=b;
    }
    return best/reps;
}

int main(int argc, char** argv){
    if(argc<2){ fprintf(stderr,"usage: %s err|errdump|time|masszero ...\n",argv[0]); return 1; }
    EPISoft::eps=1e-4;
    EPISoft::r_out=0.01;
    ForceSoft::grav_const=1.0;

    if(!strcmp(argv[1],"err") || !strcmp(argv[1],"errdump")){
        if(argc<7) return 1;
        int ni=atoi(argv[2]), ns=atoi(argv[3]);
        double scale=atof(argv[4]), offset=atof(argv[5]);
        unsigned seed=(unsigned)strtoul(argv[6],NULL,10);
        Dataset D; gen(D,ni,ns,scale,offset,seed);
        auto fr=clear_force(ni), fg=clear_force(ni);
        CalcForceEpSpQuadNoSimd ref;
        ref(D.epi.data(),ni,D.sp.data(),ns,fr.data());
        tsv110::CalcForceEpSpQuadNeon<SPJQuad> got(EPISoft::eps*EPISoft::eps, ForceSoft::grav_const);
        got(D.epi.data(),ni,D.sp.data(),ns,fg.data());
        std::vector<double> ea,ep; error_vectors(fr,fg,ea,ep);
        if(!strcmp(argv[1],"errdump")){
            printf("i,rel_acc,rel_pot\n");
            for(int i=0;i<ni;i++) printf("%d,%.6e,%.6e\n",i,ea[i],ep[i]);
            return 0;
        }
        Stats sa=stat_of(ea), sp=stat_of(ep);
        printf("ERR,%s,%d,%d,%.3f,%.1f,%u,"
               "%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,"
               "%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
               NEON_VARIANT,ni,ns,scale,offset,seed,
               sa.mx,sa.mean,sa.rms,sa.p50,sa.p90,sa.p99,
               sp.mx,sp.mean,sp.rms,sp.p50,sp.p90,sp.p99);
        return 0;
    }

    if(!strcmp(argv[1],"time")){
        if(argc<4) return 1;
        int ni=atoi(argv[2]), ns=atoi(argv[3]);
        Dataset D; gen(D,ni,ns,1.0,0.0,12345u);
        auto f1=clear_force(ni), f2=clear_force(ni);
        CalcForceEpSpQuadNoSimd ref;
        double t_ref=bench_call([&](){ ref(D.epi.data(),ni,D.sp.data(),ns,f1.data()); });
        tsv110::CalcForceEpSpQuadNeon<SPJQuad> got(EPISoft::eps*EPISoft::eps, ForceSoft::grav_const);
        double t_ne=bench_call([&](){ got(D.epi.data(),ni,D.sp.data(),ns,f2.data()); });
        printf("TIME,nosimd,%d,%d,%.1f\n",ni,ns,t_ref*1e9);
        printf("TIME,%s,%d,%d,%.1f\n",NEON_VARIANT,ni,ns,t_ne*1e9);
        return 0;
    }

    if(!strcmp(argv[1],"masszero")){
        if(argc<5) return 1;
        int ni=atoi(argv[2]), nj=atoi(argv[3]);
        unsigned seed=(unsigned)strtoul(argv[4],NULL,10);
        auto rnd=[&](){ seed=seed*1664525u+1013904223u; return (double)seed/4294967296.0; };
        std::vector<EPISoft> epi(ni);
        for(int i=0;i<ni;i++){
            epi[i].id=i;
            epi[i].pos=PS::F64vec(rnd(),rnd(),rnd());
            epi[i].r_search=0.05+0.05*rnd();
            epi[i].type=1;
        }
        std::vector<EPJSoft> epj(nj), epjf;
        for(int j=0;j<nj;j++){
            epj[j].id=j;
            epj[j].mass = (j%3==0) ? 0.0 : 1e-3*(0.5+rnd());   /* 1/3 zero-mass entries */
            epj[j].pos=PS::F64vec(rnd(),rnd(),rnd());
            epj[j].r_search=0.05+0.05*rnd();
            if(epj[j].mass>0) epjf.push_back(epj[j]);
        }
        auto f_all=clear_force(ni), f_filt=clear_force(ni), f_ne=clear_force(ni);
        CalcForceEpEpWithLinearCutoffNoSimd ref;
        ref(epi.data(),ni,epj.data(),nj,f_all.data());
        ref(epi.data(),ni,epjf.data(),(int)epjf.size(),f_filt.data());
        tsv110::CalcForceEpEpWithLinearCutoffNeon got(EPISoft::eps*EPISoft::eps, EPISoft::r_out*EPISoft::r_out, ForceSoft::grav_const);
        got(epi.data(),ni,epj.data(),nj,f_ne.data());
        std::vector<double> ea,ep; error_vectors(f_filt,f_ne,ea,ep);
        Stats sa=stat_of(ea), sp=stat_of(ep);
        int mism_ne_filt=0, mism_all_filt=0;
        for(int i=0;i<ni;i++){
            if(f_ne[i].n_ngb!=f_filt[i].n_ngb) mism_ne_filt++;
            if(f_all[i].n_ngb!=f_filt[i].n_ngb) mism_all_filt++;
        }
        printf("MASSZERO,%s,%d,%d,%u,acc_max=%.6e,pot_max=%.6e,acc_mean=%.6e,pot_mean=%.6e,"
               "nnb_mismatch_vs_filtered=%d,nnb_mismatch_unfiltered_nosimd_vs_filtered=%d\n",
               NEON_VARIANT,ni,nj,seed,
               sa.mx,sp.mx,sa.mean,sp.mean,mism_ne_filt,mism_all_filt);
        return 0;
    }

    fprintf(stderr,"unknown mode\n");
    return 1;
}
