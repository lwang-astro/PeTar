#include <iostream>
#include <cstdio>
#include <iomanip>
#include <string>
#include <unistd.h>
#include <particle_simulator.hpp>
#include "soft.hpp"
#include "io.hpp"
#include <cstdlib>

typedef PS::ParticleSystem<FPSoft> SystemSoft;

class SPJSoft: public EPJSoft{
public:
    PS::F64mat quad;

    void setQuad(const PS::F64 N) {
        mass =  10.0/N*rand()/(float)RAND_MAX;
        pos.x = 10.0*rand()/(float)RAND_MAX;
        pos.y = 10.0*rand()/(float)RAND_MAX;
        pos.z = 10.0*rand()/(float)RAND_MAX;
        quad.xx = 10.0*rand()/(float)RAND_MAX;
        quad.yy = 10.0*rand()/(float)RAND_MAX;
        quad.zz = 10.0*rand()/(float)RAND_MAX;
        quad.xy = 10.0*rand()/(float)RAND_MAX;
        quad.yz = 10.0*rand()/(float)RAND_MAX;
        quad.xz = 10.0*rand()/(float)RAND_MAX;
    };
};

int main(int argc, char **argv){
    const int Nepi = 10000;
    const int Nepj = 50000;
    const int Nspj = 50000;
    const int N = Nepi+Nepj;
    EPISoft::r_out = FPSoft::r_out = EPJSoft::r_out = 0.01;
    EPISoft::r_in = 0.001;
    EPISoft::eps = 0.0;
    const PS::F64 DF_MAX=1e-4;
    
    SystemSoft ptcl;
    PS::S32 n=0;
    SetParticlePlummer(ptcl, N, n);
    EPISoft epi[Nepi];
    EPJSoft epj[Nepj];
    SPJSoft spj[Nspj];
    ForceSoft force[Nepi], force_simd[Nepi];
    ForceSoft force_sp[Nepi], force_sp_simd[Nepi];
    ForceSoft force_sp_quad[Nepi], force_sp_quad_simd[Nepi];
    for (int i=0; i<N; i++) ptcl[i].calcRSearch(1.0/2048.0);
    for (int i=0; i<Nepi; i++) {
        epi[i].copyFromFP(ptcl[i]);
        epj[i].copyFromFP(ptcl[i+Nepi]);
        force[i].clear();
        force_simd[i].clear();
        force_sp_simd[i].clear();
        force_sp[i].clear();
        force_sp_quad_simd[i].clear();
        force_sp_quad[i].clear();
    }
    //for (int i=0; i<Nepj; i++) epj[i].copyFromFP(ptcl[i]);
    for (int i=0; i<Nspj; i++) spj[i].setQuad((PS::F64)N);

    CalcForceEpEpWithLinearCutoffSimd f_ep_ep_simd;
    CalcForceEpEpWithLinearCutoffNoSimd f_ep_ep;
    CalcForceEpSpMonoSimd f_ep_sp_simd;
    CalcForceEpSpMonoNoSimd f_ep_sp;
    CalcForceEpSpQuadSimd f_ep_sp_quad_simd;
    CalcForceEpSpQuadNoSimd  f_ep_sp_quad;
        
    PS::F64 t_ep_simd=0, t_ep_no=0, t_sp_simd=0, t_sp_no=0, t_sp_quad_simd=0, t_sp_quad_no=0;
    t_ep_simd -= PS::GetWtime();
    f_ep_ep_simd(epi, Nepi, epj, Nepj, force_simd);
    t_ep_simd += PS::GetWtime();

    t_ep_no -= PS::GetWtime();
    f_ep_ep(epi, Nepi, epj, Nepj, force);
    t_ep_no += PS::GetWtime();

    t_sp_simd -= PS::GetWtime();
    f_ep_sp_simd(epi, Nepi, epj, Nepj, force_sp_simd);
    t_sp_simd += PS::GetWtime();

    t_sp_no -= PS::GetWtime();
    f_ep_sp(epi, Nepi, epj, Nepj, force_sp);
    t_sp_no += PS::GetWtime();

    t_sp_quad_simd -= PS::GetWtime();
    f_ep_sp_quad_simd(epi, Nepi, spj, Nspj, force_sp_quad_simd);
    t_sp_quad_simd += PS::GetWtime();

    t_sp_quad_no -= PS::GetWtime();
    f_ep_sp_quad(epi, Nepi, spj, Nspj, force_sp_quad);
    t_sp_quad_no += PS::GetWtime();

    PS::F64 dfmax=0,dsmax=0,dqmax=0,dpmax=0;
    for(int i=0; i<Nepi; i++) {
        for (int j=0; j<3; j++) {
            PS::F64 df=(force[i].acc[j]-force_simd[i].acc[j])/force[i].acc[j];
            dfmax = std::max(dfmax, df);
            if(df>DF_MAX) std::cerr<<"Force diff: i="<<i<<" nosimd["<<j<<"] "<<force[i].acc[j]<<" simd["<<j<<"] "<<force_simd[i].acc[j]<<std::endl;
            df = (force_sp[i].acc[j]-force_sp_simd[i].acc[j])/force_sp[i].acc[j];
            dsmax = std::max(dsmax, df);
            if(df>DF_MAX) std::cerr<<"Force sp diff: i="<<i<<" nosimd["<<j<<"] "<<force_sp[i].acc[j]<<" simd["<<j<<"] "<<force_sp_simd[i].acc[j]<<std::endl;
            df = (force_sp_quad[i].acc[j]-force_sp_quad_simd[i].acc[j])/force_sp_quad[i].acc[j];
            dqmax = std::max(dqmax, df);
            if(df>DF_MAX) std::cerr<<"Force sp_quad diff: i="<<i<<" nosimd["<<j<<"] "<<force_sp_quad[i].acc[j]<<" simd["<<j<<"] "<<force_sp_quad_simd[i].acc[j]<<std::endl;
            
        }
        dpmax = std::max(dpmax, (force[i].pot-force_simd[i].pot)/force[i].pot);
        if(force[i].n_ngb!=force_simd[i].n_ngb) {
            std::cerr<<"Neighbor diff: i="<<i<<" nosimd "<<force[i].n_ngb<<" simd "<<force_simd[i].n_ngb<<std::endl;
        }
    }
    std::cout<<"Force diff max: "<<dfmax<<"\nPot diff max: "<<dpmax<<"\nSp diff max: "<<dsmax<<"\nQuad diff max: "<<dqmax<<std::endl;
    
    std::cout<<"Time: epj  simd="<<t_ep_simd<<" no="<<t_ep_no<<" ratio="<<t_ep_no/t_ep_simd<<std::endl;
    std::cout<<"Time: spj  simd="<<t_sp_simd<<" no="<<t_sp_no<<" ratio="<<t_sp_no/t_sp_simd<<std::endl;
    std::cout<<"Time: quad simd="<<t_sp_quad_simd<<" no="<<t_sp_quad_no<<" ratio="<<t_sp_quad_no/t_sp_quad_simd<<std::endl;
    
    return 0;
}
