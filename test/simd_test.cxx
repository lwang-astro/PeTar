#include <iostream>
#include <cstdio>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <cassert>
#define ASSERT assert
#include <unistd.h>
#include <particle_simulator.hpp>
#include "soft_ptcl.hpp"
#include "soft_force.hpp"
#include "io.hpp"
#include "particle_distribution_generator.hpp"
#include "force_gpu_cuda.hpp"
#include "static_variables.hpp"

void setQuad(const PS::F64 N, SPJSoft& sp) {
    sp.mass =  1.0/N+0.001/N*rand()/(float)RAND_MAX;
    sp.pos.x = 1.0+10.0*rand()/(float)RAND_MAX;
    sp.pos.y = 1.0+10.0*rand()/(float)RAND_MAX;
    sp.pos.z = 1.0+10.0*rand()/(float)RAND_MAX;
    sp.quad.xx = 10.0*rand()/(float)RAND_MAX;
    sp.quad.yy = 10.0*rand()/(float)RAND_MAX;
    sp.quad.zz = 10.0*rand()/(float)RAND_MAX;
    sp.quad.xy = 10.0*rand()/(float)RAND_MAX;
    sp.quad.yz = 10.0*rand()/(float)RAND_MAX;
    sp.quad.xz = 10.0*rand()/(float)RAND_MAX;
}

int main(int argc, char **argv){
    const int Nepi = 10000;
    const int Nepj = 50000;
    const int Nspj = 50000;
    const int N = Nepi+Nepj;
    EPISoft::r_out = 0.01;
    EPISoft::eps = 0.0;
    ForceSoft::grav_const = 1.0;
    const PS::F64 DF_MAX=7e-3;
    
    std::cout<<"Make model N="<<N<<"\n";
    FPSoft ptcl[N];
    PS::S32 n=0;
    PS::F64 * mass;
    PS::F64vec * pos;
    PS::F64vec * vel;
    ParticleDistributionGenerator::makePlummerModel(1.0, N, N, mass, pos, vel, -0.25, n);
    std::cout<<"Initial data"<<std::endl;
    for(PS::S32 i=0; i<N; i++){
        ptcl[i].mass = mass[i];
        ptcl[i].pos = pos[i];
        ptcl[i].vel = vel[i];
        ptcl[i].id =  i + 1;
        ptcl[i].group_data.artificial.setParticleTypeToSingle();
        ptcl[i].changeover.setR(1.0, 0.001, 0.01);
    }    

    EPISoft epi[Nepi];
    EPJSoft epj[Nepj];
    SPJSoft spj[Nspj];
    ForceSoft force[Nepi], force_simd[Nepi], force_gpu[Nepi];
    ForceSoft force_sp[Nepi], force_sp_simd[Nepi];
    ForceSoft force_sp_quad[Nepi], force_sp_quad_simd[Nepi];
#pragma omp parallel for
    for (int i=0; i<N; i++) ptcl[i].calcRSearch(1.0/2048.0);
#pragma omp parallel for    
    for (int i=0; i<Nepi; i++) {
        epi[i].copyFromFP(ptcl[i]);
        epj[i].copyFromFP(ptcl[i+Nepi]);
        force[i].clear();
        force_simd[i].clear();
        force_gpu[i].clear();
        force_sp_simd[i].clear();
        force_sp[i].clear();
        force_sp_quad_simd[i].clear();
        force_sp_quad[i].clear();
    }
    //for (int i=0; i<Nepj; i++) epj[i].copyFromFP(ptcl[i]);
#ifdef USE_QUAD    
    for (int i=0; i<Nspj; i++) setQuad((PS::F64)N, spj[i]);
#endif

    CalcForceEpEpWithLinearCutoffSimd f_ep_ep_simd;
    CalcForceEpEpWithLinearCutoffNoSimd f_ep_ep;
    CalcForceEpSpMonoSimd f_ep_sp_simd;
    CalcForceEpSpMonoNoSimd f_ep_sp;
    CalcForceEpSpQuadSimd f_ep_sp_quad_simd;
    CalcForceEpSpQuadNoSimd  f_ep_sp_quad;

    PS::F64 t_ep_simd=0, t_ep_no=0, t_sp_simd=0, t_sp_no=0, t_sp_quad_simd=0, t_sp_quad_no=0, t_gpu=0;
    std::cout<<"calc Ep Ep simd\n";
    t_ep_simd -= PS::GetWtime();
    f_ep_ep_simd(epi, Nepi, epj, Nepj, force_simd);
    t_ep_simd += PS::GetWtime();

    std::cout<<"calc Ep Ep\n";
    t_ep_no -= PS::GetWtime();
    f_ep_ep(epi, Nepi, epj, Nepj, force);
    t_ep_no += PS::GetWtime();

    std::cout<<"calc Ep Sp mono simd\n";
    t_sp_simd -= PS::GetWtime();
    f_ep_sp_simd(epi, Nepi, epj, Nepj, force_sp_simd);
    t_sp_simd += PS::GetWtime();

    std::cout<<"calc Ep Sp mono\n";
    t_sp_no -= PS::GetWtime();
    f_ep_sp(epi, Nepi, epj, Nepj, force_sp);
    t_sp_no += PS::GetWtime();

#ifdef USE_QUAD
    std::cout<<"calc Ep Sp quad simd\n";
    t_sp_quad_simd -= PS::GetWtime();
    f_ep_sp_quad_simd(epi, Nepi, spj, Nspj, force_sp_quad_simd);
    t_sp_quad_simd += PS::GetWtime();

    std::cout<<"calc Ep Sp quad\n";
    t_sp_quad_no -= PS::GetWtime();
    f_ep_sp_quad(epi, Nepi, spj, Nspj, force_sp_quad);
    t_sp_quad_no += PS::GetWtime();
#endif

    std::cout<<"calc GPU\n";
    t_gpu -= PS::GetWtime();
    const EPISoft *epi_ptr = epi;
    const EPJSoft *epj_ptr = epj;
    const SPJSoft *spj_ptr = spj;
    ForceSoft *force_gpu_ptr = force_gpu;
    DispatchKernelWithSP(1, 1, &epi_ptr, &Nepi, &epj_ptr, &Nepj, &spj_ptr, &Nspj);
    RetrieveKernel(1, 1, &Nepi, &force_gpu_ptr);
    t_gpu += PS::GetWtime();

    std::cout<<"compare results\n";
    PS::F64 dfmax=0,dsmax=0,dqmax=0,dgmax=0,dfpmax=0,dspmax=0,dqpmax=0,dgpmax=0;
    for(int i=0; i<Nepi; i++) {
        for (int j=0; j<3; j++) {
            PS::F64 df=(force[i].acc[j]-force_simd[i].acc[j])/force[i].acc[j];
            dfmax = std::max(dfmax, df);
            if(df>DF_MAX) std::cerr<<"Force diff: i="<<i<<" nosimd["<<j<<"] "<<force[i].acc[j]<<" simd["<<j<<"] "<<force_simd[i].acc[j]<<std::endl;
            df = (force_sp[i].acc[j]-force_sp_simd[i].acc[j])/force_sp[i].acc[j];
            dsmax = std::max(dsmax, df);
            if(df>DF_MAX) std::cerr<<"Force sp diff: i="<<i<<" nosimd["<<j<<"] "<<force_sp[i].acc[j]<<" simd["<<j<<"] "<<force_sp_simd[i].acc[j]<<std::endl;
#ifdef USE_QUAD
            df = (force_sp_quad[i].acc[j]-force_sp_quad_simd[i].acc[j])/force_sp_quad[i].acc[j];
            dqmax = std::max(dqmax, df);
            if(df>DF_MAX) std::cerr<<"Force sp_quad diff: i="<<i<<" nosimd["<<j<<"] "<<force_sp_quad[i].acc[j]<<" simd["<<j<<"] "<<force_sp_quad_simd[i].acc[j]<<std::endl;

            df = (force[i].acc[j]+force_sp_quad[i].acc[j] - force_gpu[i].acc[j])/force_gpu[i].acc[j];
            dgmax = std::max(dgmax, df);
            if(df>DF_MAX) std::cerr<<"Force diff: i="<<i<<" nosimd["<<j<<"] "<<force[i].acc[j]+force_sp_quad[i].acc[j]<<" gpu["<<j<<"] "<<force_gpu[i].acc[j]<<std::endl;
#else
            df = (force[i].acc[j]+force_sp[i].acc[j] - force_gpu[i].acc[j])/force_gpu[i].acc[j];
            dgmax = std::max(dgmax, df);
            if(df>DF_MAX) std::cerr<<"Force diff: i="<<i<<" nosimd["<<j<<"] "<<force[i].acc[j]+force_sp[i].acc[j]<<" gpu["<<j<<"] "<<force_gpu[i].acc[j]<<std::endl;
#endif
        }
        dfpmax = std::max(dfpmax, (force[i].pot-force_simd[i].pot)/force[i].pot);
        dspmax = std::max(dspmax, (force_sp[i].pot-force_sp_simd[i].pot)/force_sp[i].pot);
#ifdef USE_QUAD
        dqpmax = std::max(dqpmax, (force_sp_quad[i].pot-force_sp_quad_simd[i].pot)/force_sp_quad[i].pot);
        dgpmax = std::max(dqpmax, (force_sp_quad[i].pot+force[i].pot - force_gpu[i].pot)/force_gpu[i].pot);
#else
        dgpmax = std::max(dqpmax, (force_sp[i].pot+force[i].pot - force_gpu[i].pot)/force_gpu[i].pot);
#endif
        if(force[i].n_ngb!=force_simd[i].n_ngb) {
            std::cerr<<"Neighbor diff: i="<<i<<" nosimd "<<force[i].n_ngb<<" simd "<<force_simd[i].n_ngb<<std::endl;
        }
        if(force[i].n_ngb!=force_gpu[i].n_ngb) {
            std::cerr<<"Neighbor diff: i="<<i<<" nosimd "<<force[i].n_ngb<<" gpu "<<force_gpu[i].n_ngb<<std::endl;
        }
    }
#ifdef USE__AVX512
    std::cout<<"Use AVX512";
#elif defined(__AVX2__)
    std::cout<<"Use AVX2";
#elif defined(__AVX__)
    std::cout<<"Use AVX";
#endif
#ifdef USE_QUAD
    std::cout<<" GPU_quad";
#else
    std::cout<<" GPU_mono";
#endif
#ifdef CALC_EP_64bit
    std::cout<<" EP_64bit";
#else
    std::cout<<" EP_32bit";
#endif
#ifdef CALC_SP_64bit
    std::cout<<" SP_64bit";
#else
    std::cout<<" SP_32bit";
#endif
#ifdef RSQRT_NR_EPJ_X4
    std::cout<<" EP_RSQRT_X4";
#elif defined (RSQRT_NR_EPJ_X2)
    std::cout<<" EP_RSQRT_X2";
#endif
#ifdef RSQRT_NR_SPJ_X4
    std::cout<<" SP_RSQRT_X4";
#elif defined (RSQRT_NR_SPJ_X2)
    std::cout<<" SP_RSQRT_X2";
#endif
#ifdef AVX_PRELOAD
    std::cout<<" PRELOAD";
#endif
    std::cout<<std::endl;
    
    std::cout<<"Force diff max: "<<dfmax<<" Pot diff max: "<<dfpmax
             <<"\nSp diff max: "<<dsmax<<" Pot diff max: "<<dspmax
             <<"\nQuad diff max: "<<dqmax<<" Pot diff max: "<<dqpmax
             <<"\nGPU diff max: "<<dgmax<<" Pot diff max: "<<dgpmax
             <<std::endl;
    
    std::cout<<"Time: epj  simd="<<t_ep_simd<<" no="<<t_ep_no<<" ratio="<<t_ep_no/t_ep_simd<<std::endl;
    std::cout<<"Time: spj  simd="<<t_sp_simd<<" no="<<t_sp_no<<" ratio="<<t_sp_no/t_sp_simd<<std::endl;
    std::cout<<"Time: quad simd="<<t_sp_quad_simd<<" no="<<t_sp_quad_no<<" ratio="<<t_sp_quad_no/t_sp_quad_simd<<std::endl;
    std::cout<<"Time: gpu ="<<t_gpu<<" no="<<t_ep_no+t_sp_no<<" ratio="<<t_gpu/(t_ep_no+t_sp_no)<<std::endl;

//    for (int i=0; i<Nepi; i++) {
//        force_sp_quad_simd[i].clear();
//        force_sp_quad[i].clear();
//    }
//    for (int i=0; i<Nspj; i++) {
//        int id=4783;
//        f_ep_sp_quad_simd(&epi[id], 16, &spj[i], 1, &force_sp_quad_simd[id]);
//        f_ep_sp_quad(&epi[id], 16, &spj[i], 1, &force_sp_quad[id]);
// 
//        for(int j=0; j<3; j++) {
//            double df = (force_sp_quad[id].acc[j]-force_sp_quad_simd[id].acc[j])/force_sp_quad[id].acc[j];
//            if(df>DF_MAX) {
//                std::cerr<<"Force sp_quad diff: i="<<i<<" nosimd["<<j<<"] "<<force_sp_quad[id].acc[j]<<" simd["<<j<<"] "<<force_sp_quad_simd[id].acc[j]<<std::endl;
//            }
//            
//        }
//        force_sp_quad_simd[id].clear();
//        force_sp_quad[id].clear();
//    }
    return 0;
}
