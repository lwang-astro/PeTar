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
#include "static_variables.hpp"
#ifdef USE_GPU
#include "force_gpu_cuda.hpp"
#else
#ifdef USE_QUAD
#define SPJSoft PS::SPJQuadrupoleInAndOut
#else
#define SPJSoft PS::SPJMonopoleInAndOut
#endif
#ifdef USE_FUGAKU
#include "force_fugaku.hpp"
#endif
#endif

#ifdef USE_QUAD
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
#endif

int main(int argc, char **argv){
    const int Nepi = 10;
    const int Nepj = 20;
#ifdef USE_QUAD
    const int Nspj = 10;
#endif
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

    ForceSoft force[Nepi];
    ForceSoft force_sp[Nepi];
#ifdef USE_QUAD
    SPJSoft spj[Nspj];
    ForceSoft force_sp_quad[Nepi];
#endif
#ifdef USE_GPU
    ForceSoft force_gpu[Nepi];
#endif
#ifdef USE_SIMD
    ForceSoft force_simd[Nepi];
    ForceSoft force_sp_simd[Nepi];
#ifdef USE_QUAD
    ForceSoft force_sp_quad_simd[Nepi];
#endif
#endif
    #ifdef USE_FUGAKU
    ForceSoft force_fgk[Nepi];
    //ForceSoft force_sp_fgk[Nepi];
#endif

#pragma omp parallel for
    for (int i=0; i<N; i++) ptcl[i].calcRSearch(1.0/2048.0);
#pragma omp parallel for    
    for (int i=0; i<Nepi; i++) {
        epi[i].copyFromFP(ptcl[i]);
        epj[i].copyFromFP(ptcl[i+Nepi]);
        force[i].clear();
#ifdef USE_GPU
        force_gpu[i].clear();
#endif
        force_sp[i].clear();
#ifdef USE_QUAD
        force_sp_quad[i].clear();
#endif
#ifdef USE_SIMD
        force_simd[i].clear();
        force_sp_simd[i].clear();
#ifdef USE_QUAD
        force_sp_quad_simd[i].clear();
#endif
#endif
#ifdef USE_FUGAKU
        force_fgk[i].clear();
        //force_sp_fgk[i].clear();
#endif
    }
    for (int i=Nepi; i<Nepj; i++) 
      epj[i].copyFromFP(ptcl[i+Nepi]);
    //for (int i=0; i<Nepj; i++) epj[i].copyFromFP(ptcl[i]);
#ifdef USE_QUAD    
    for (int i=0; i<Nspj; i++) setQuad((PS::F64)N, spj[i]);
#endif

#ifdef USE_GPU
    PS::F64 t_gpu=0;

    std::cout<<"calc GPU\n";
    t_gpu -= PS::GetWtime();
    const EPISoft *epi_ptr = epi;
    const EPJSoft *epj_ptr = epj;
    const SPJSoft *spj_ptr = spj;
    ForceSoft *force_gpu_ptr = force_gpu;
#ifdef PARTICLE_SIMULATOR_GPU_MULIT_WALK_INDEX
    PS::S32 id_epj[Nepj], id_spj[Nspj];
    const PS::S32* id_epj_ptr=id_epj;
    const PS::S32* id_spj_ptr=id_spj;
    for (int i=0; i<Nepj; i++) id_epj[i]=i;
    for (int i=0; i<Nspj; i++) id_spj[i]=i;
    DispatchKernelWithSPIndex(1, 1, &epi_ptr, &Nepi, &id_epj_ptr, &Nepj, &id_spj_ptr, &Nspj, epj_ptr, Nepj, spj_ptr, Nspj, true);
    DispatchKernelWithSPIndex(1, 1, &epi_ptr, &Nepi, &id_epj_ptr, &Nepj, &id_spj_ptr, &Nspj, epj_ptr, Nepj, spj_ptr, Nspj, false);
#else
    DispatchKernelWithSP(1, 1, &epi_ptr, &Nepi, &epj_ptr, &Nepj, &spj_ptr, &Nspj);
#endif
    RetrieveKernel(1, 1, &Nepi, &force_gpu_ptr);
    t_gpu += PS::GetWtime();
#endif

#ifdef USE_SIMD
    std::cout<<"calc Ep Ep simd\n";
    CalcForceEpEpWithLinearCutoffSimd f_ep_ep_simd;
    PS::F64 t_ep_simd=0;
    t_ep_simd -= PS::GetWtime();
    f_ep_ep_simd(epi, Nepi, epj, Nepj, force_simd);
    t_ep_simd += PS::GetWtime();

    std::cout<<"calc Ep Sp mono simd\n";
    CalcForceEpSpMonoSimd f_ep_sp_simd;
    PS::F64 t_sp_simd=0;
    t_sp_simd -= PS::GetWtime();
    f_ep_sp_simd(epi, Nepi, epj, Nepj, force_sp_simd);
    t_sp_simd += PS::GetWtime();

#ifdef USE_QUAD
    std::cout<<"calc Ep Sp quad simd\n";
    CalcForceEpSpQuadSimd f_ep_sp_quad_simd;
    PS::F64 t_sp_quad_simd=0;
    t_sp_quad_simd -= PS::GetWtime();
    f_ep_sp_quad_simd(epi, Nepi, spj, Nspj, force_sp_quad_simd);
    t_sp_quad_simd += PS::GetWtime();
#endif
#endif

#ifdef USE_FUGAKU
    std::cout<<"calc Ep Ep fugaku\n";
    CalcForceEpEpWithLinearCutoffFugaku f_ep_ep_fgk(EPISoft::eps*EPISoft::eps, EPISoft::r_out*EPISoft::r_out, ForceSoft::grav_const);
    PS::F64 t_ep_fgk=0;
    t_ep_fgk -= PS::GetWtime();
    f_ep_ep_fgk(epi, Nepi, epj, Nepj, force_fgk);
    t_ep_fgk += PS::GetWtime();
#endif

    std::cout<<"calc Ep Ep\n";
    CalcForceEpEpWithLinearCutoffNoSimd f_ep_ep;
    PS::F64 t_ep_no=0;
    t_ep_no -= PS::GetWtime();
    f_ep_ep(epi, Nepi, epj, Nepj, force);
    t_ep_no += PS::GetWtime();

    std::cout<<"calc Ep Sp mono\n";
    CalcForceEpSpMonoNoSimd f_ep_sp;
    PS::F64 t_sp_no=0;
    t_sp_no -= PS::GetWtime();
    f_ep_sp(epi, Nepi, epj, Nepj, force_sp);
    t_sp_no += PS::GetWtime();

#ifdef USE_QUAD
    std::cout<<"calc Ep Sp quad\n";
    CalcForceEpSpQuadNoSimd  f_ep_sp_quad;
    PS::F64 t_sp_quad_no=0;
    t_sp_quad_no -= PS::GetWtime();
    f_ep_sp_quad(epi, Nepi, spj, Nspj, force_sp_quad);
    t_sp_quad_no += PS::GetWtime();
#endif

    std::cout<<"compare results\n";
#ifdef USE_SIMD
    PS::F64 dfmax=0,dsmax=0;
    PS::F64 dfpmax=0,dspmax=0;
#ifdef USE_QUAD
    PS::F64 dqmax=0,dqpmax=0;
#endif
#endif
#ifdef USE_GPU
    PS::F64 dgmax=0,dgpmax=0;
#endif
#ifdef USE_FUGAKU
    PS::F64 dfgkmax=0,dfgkpmax=0;
#endif
    PS::F64 df;
    for(int i=0; i<Nepi; i++) {
        for (int j=0; j<3; j++) {
#ifdef USE_SIMD
            df=(force[i].acc[j]-force_simd[i].acc[j])/force[i].acc[j];
            dfmax = std::max(dfmax, df);
            if(df>DF_MAX) std::cerr<<"Force diff: i="<<i<<" nosimd["<<j<<"] "<<force[i].acc[j]<<" simd["<<j<<"] "<<force_simd[i].acc[j]<<std::endl;
            df = (force_sp[i].acc[j]-force_sp_simd[i].acc[j])/force_sp[i].acc[j];
            dsmax = std::max(dsmax, df);
            if(df>DF_MAX) std::cerr<<"Force sp diff: i="<<i<<" nosimd["<<j<<"] "<<force_sp[i].acc[j]<<" simd["<<j<<"] "<<force_sp_simd[i].acc[j]<<std::endl;
#ifdef USE_QUAD
            df = (force_sp_quad[i].acc[j]-force_sp_quad_simd[i].acc[j])/force_sp_quad[i].acc[j];
            dqmax = std::max(dqmax, df);
            if(df>DF_MAX) std::cerr<<"Force sp_quad diff: i="<<i<<" nosimd["<<j<<"] "<<force_sp_quad[i].acc[j]<<" simd["<<j<<"] "<<force_sp_quad_simd[i].acc[j]<<std::endl;
#endif
#endif
#ifdef USE_GPU
#ifdef USE_QUAD
            df = (force[i].acc[j]+force_sp_quad[i].acc[j] - force_gpu[i].acc[j])/force_gpu[i].acc[j];
            dgmax = std::max(dgmax, df);
            if(df>DF_MAX) std::cerr<<"Force diff: i="<<i<<" nosimd["<<j<<"] "<<force[i].acc[j]+force_sp_quad[i].acc[j]<<" gpu["<<j<<"] "<<force_gpu[i].acc[j]<<std::endl;
            df = (force[i].acc[j]+force_sp[i].acc[j] - force_gpu[i].acc[j])/force_gpu[i].acc[j];
#else
            dgmax = std::max(dgmax, df);
            if(df>DF_MAX) std::cerr<<"Force diff: i="<<i<<" nosimd["<<j<<"] "<<force[i].acc[j]+force_sp[i].acc[j]<<" gpu["<<j<<"] "<<force_gpu[i].acc[j]<<std::endl;
#endif
#endif
#ifdef USE_FUGAKU
            df=(force[i].acc[j]-force_fgk[i].acc[j])/force[i].acc[j];
            dfgkmax = std::max(dfgkmax, df);
            if(df>DF_MAX) std::cerr<<"Force diff: i="<<i<<" nosimd["<<j<<"] "<<force[i].acc[j]<<" fugaku["<<j<<"] "<<force_fgk[i].acc[j]<<std::endl;
#endif
        }
#ifdef USE_SIMD
        dfpmax = std::max(dfpmax, (force[i].pot-force_simd[i].pot)/force[i].pot);
        dspmax = std::max(dspmax, (force_sp[i].pot-force_sp_simd[i].pot)/force_sp[i].pot);
#ifdef USE_QUAD
        dqpmax = std::max(dqpmax, (force_sp_quad[i].pot-force_sp_quad_simd[i].pot)/force_sp_quad[i].pot);
#endif
        if(force[i].n_ngb!=force_simd[i].n_ngb) {
            std::cerr<<"Neighbor diff: i="<<i<<" nosimd "<<force[i].n_ngb<<" simd "<<force_simd[i].n_ngb<<std::endl;
        }
#endif
#ifdef USE_GPU
#ifdef USE_QUAD
        dgpmax = std::max(dqpmax, (force_sp_quad[i].pot+force[i].pot - force_gpu[i].pot)/force_gpu[i].pot);
#else
        dgpmax = std::max(dqpmax, (force_sp[i].pot+force[i].pot - force_gpu[i].pot)/force_gpu[i].pot);
#endif
        if(force[i].n_ngb!=force_gpu[i].n_ngb) {
            std::cerr<<"Neighbor diff: i="<<i<<" nosimd "<<force[i].n_ngb<<" gpu "<<force_gpu[i].n_ngb<<std::endl;
        }
#endif
#ifdef USE_FUGAKU
	dfgkpmax = std::max(dfgkpmax, (force[i].pot-force_fgk[i].pot)/force[i].pot);
        if(force[i].n_ngb!=force_fgk[i].n_ngb) {
            std::cerr<<"Neighbor diff: i="<<i<<" nosimd "<<force[i].n_ngb<<" fugaku "<<force_fgk[i].n_ngb<<std::endl;
        }
#endif
    }
#ifdef USE__AVX512
    std::cout<<"Use AVX512";
#elif defined(__AVX2__)
    std::cout<<"Use AVX2";
#elif defined(__AVX__)
    std::cout<<"Use AVX";
#endif
#ifdef USE_GPU
#ifdef USE_QUAD
    std::cout<<" GPU_quad";
#else
    std::cout<<" GPU_mono";
#endif
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

#ifdef USE_SIMD    
    std::cout<<"Force diff max: "<<dfmax<<" Pot diff max: "<<dfpmax
             <<"\nSp diff max: "<<dsmax<<" Pot diff max: "<<dspmax;
#ifdef USE_QUAD       
    std::cout<<"\nQuad diff max: "<<dqmax<<" Pot diff max: "<<dqpmax;
#endif
#endif
#ifdef USE_GPU
    std::cout<<"\nGPU diff max: "<<dgmax<<" Pot diff max: "<<dgpmax;
#endif
#ifdef USE_FUGAKU
    std::cout<<"\nFugaku diff max: "<<dfgkmax<<" Pot diff max: "<<dfgkpmax;
#endif
    std::cout<<std::endl;

#ifdef USE_SIMD
    std::cout<<"Time: epj  simd="<<t_ep_simd<<" no="<<t_ep_no<<" ratio="<<t_ep_no/t_ep_simd<<std::endl;
    std::cout<<"Time: spj  simd="<<t_sp_simd<<" no="<<t_sp_no<<" ratio="<<t_sp_no/t_sp_simd<<std::endl;
#ifdef USE_QUAD
    std::cout<<"Time: quad simd="<<t_sp_quad_simd<<" no="<<t_sp_quad_no<<" ratio="<<t_sp_quad_no/t_sp_quad_simd<<std::endl;
#endif
#endif
#ifdef USE_GPU
    std::cout<<"Time: gpu ="<<t_gpu<<" no="<<t_ep_no+t_sp_no<<" ratio="<<(t_ep_no+t_sp_no)/t_gpu<<std::endl;
#endif
#ifdef USE_FUGAKU
    std::cout<<"Time: fugaku ="<<t_ep_fgk<<" no="<<t_ep_no<<" ratio="<<t_ep_no/t_ep_fgk<<std::endl;
#endif

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
