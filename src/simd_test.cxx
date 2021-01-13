#include <iostream>
#include <cstdio>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <cassert>
#define ASSERT assert
#include <unistd.h>
#include <particle_simulator.hpp>

#ifdef P3T_64BIT
#define CALC_EP_64bit
#define CALC_SP_64bit
#define RSQRT_NR_EPJ_X4
#define RSQRT_NR_SPJ_X4

#elif P3T_MIXBIT
#define CALC_EP_64bit
#define RSQRT_NR_EPJ_X4

#else
#define RSQRT_NR_EPJ_X2
//#define RSQRT_NR_SPJ_X2
#endif 

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

void setSpj(const PS::F64 N, SPJSoft& sp) {
    sp.mass =  1.0/N+0.001/N*rand()/(float)RAND_MAX;
    sp.pos.x = 1.0+10.0*rand()/(float)RAND_MAX;
    sp.pos.y = 1.0+10.0*rand()/(float)RAND_MAX;
    sp.pos.z = 1.0+10.0*rand()/(float)RAND_MAX;
#ifdef USE_QUAD
    sp.quad.xx = 10.0*rand()/(float)RAND_MAX;
    sp.quad.yy = 10.0*rand()/(float)RAND_MAX;
    sp.quad.zz = 10.0*rand()/(float)RAND_MAX;
    sp.quad.xy = 10.0*rand()/(float)RAND_MAX;
    sp.quad.yz = 10.0*rand()/(float)RAND_MAX;
    sp.quad.xz = 10.0*rand()/(float)RAND_MAX;
#endif
}

int main(int argc, char **argv){
    const int Nepi = 1000;
    const int Nepj = 2000;
    const int Nspj = 1000;

    const int N = std::max(Nepi,Nepj);
    EPISoft::r_out = 0.01;
    EPISoft::eps = 1e-4;
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

    ForceSoft force[Nepi];
    ForceSoft force_sp[Nepi];
    ForceSoft force_nb[Nepi];
#ifdef USE_GPU
    ForceSoft force_gpu[Nepi];
#endif
#ifdef USE_SIMD
    ForceSoft force_simd[Nepi];
    ForceSoft force_sp_simd[Nepi];
    ForceSoft force_nb_simd[Nepi];
#endif
#ifdef USE_FUGAKU
    ForceSoft force_fgk[Nepi];
    ForceSoft force_sp_fgk[Nepi];
    ForceSoft force_nb_fgk[Nepi];
#endif

    for (int i=0; i<N; i++) ptcl[i].calcRSearch(1.0/2048.0);

    for (int i=0; i<Nepi; i++) {
        epi[i].copyFromFP(ptcl[i]);
        force[i].clear();
#ifdef USE_GPU
        force_gpu[i].clear();
#endif
        force_sp[i].clear();
#ifdef USE_SIMD
        force_simd[i].clear();
        force_sp_simd[i].clear();
        force_nb_simd[i].clear();
#endif
#ifdef USE_FUGAKU
        force_fgk[i].clear();
        force_sp_fgk[i].clear();
        force_nb_fgk[i].clear();
#endif
    }
    for (int i=0; i<Nepj; i++) 
        epj[i].copyFromFP(ptcl[i]);

    for (int i=0; i<Nspj; i++) setSpj((PS::F64)N, spj[i]);

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
    CalcForceWithLinearCutoffCUDAMultiWalk f_ep_ep_gpu(0, EPISoft::eps*EPISoft::eps, EPISoft::r_out*EPISoft::r_out, ForceSoft::grav_const);
    f_ep_ep_gpu(1, 1, &epi_ptr, &Nepi, &id_epj_ptr, &Nepj, &id_spj_ptr, &Nspj, epj_ptr, Nepj, spj_ptr, Nspj, true);
    f_ep_ep_gpu(1, 1, &epi_ptr, &Nepi, &id_epj_ptr, &Nepj, &id_spj_ptr, &Nspj, epj_ptr, Nepj, spj_ptr, Nspj, false);
#else
    CalcForceWithLinearCutoffCUDA f_ep_ep_gpu(0, EPISoft::eps*EPISoft::eps, EPISoft::r_out*EPISoft::r_out, ForceSoft::grav_const);
    f_ep_ep_gpu(1, 1, &epi_ptr, &Nepi, &epj_ptr, &Nepj, &spj_ptr, &Nspj);
#endif
    RetrieveForceCUDA(1, 1, &Nepi, &force_gpu_ptr);
    t_gpu += PS::GetWtime();
#endif

#ifdef USE_SIMD
    std::cout<<"calc Ep Ep simd\n";
    CalcForceEpEpWithLinearCutoffSimd f_ep_ep_simd;
    PS::F64 t_ep_simd=0;
    t_ep_simd -= PS::GetWtime();
    f_ep_ep_simd(epi, Nepi, epj, Nepj, force_simd);
    t_ep_simd += PS::GetWtime();

#ifdef USE_QUAD
    std::cout<<"calc Ep Sp quad simd\n";
    CalcForceEpSpQuadSimd f_ep_sp_simd;
#else
    std::cout<<"calc Ep Sp mono simd\n";
    CalcForceEpSpMonoSimd f_ep_sp_simd;
#endif
    PS::F64 t_sp_simd=0;
    t_sp_simd -= PS::GetWtime();
    f_ep_sp_simd(epi, Nepi, spj, Nspj, force_sp_simd);
    t_sp_simd += PS::GetWtime();

    std::cout<<"neighbor search simd\n";
    SearchNeighborEpEpSimd f_nb_simd;
    PS::F64 t_nb_simd=0;
    t_nb_simd -= PS::GetWtime();
    f_nb_simd(epi, Nepi, epj, Nepj, force_nb_simd);
    t_nb_simd += PS::GetWtime();
#endif

#ifdef USE_FUGAKU
    std::cout<<"calc Ep Ep fugaku\n";
    CalcForceEpEpWithLinearCutoffFugaku f_ep_ep_fgk(EPISoft::eps*EPISoft::eps, EPISoft::r_out*EPISoft::r_out, ForceSoft::grav_const);
    PS::F64 t_ep_fgk=0;
    t_ep_fgk -= PS::GetWtime();
    f_ep_ep_fgk(epi, Nepi, epj, Nepj, force_fgk);
    t_ep_fgk += PS::GetWtime();

#ifdef USE_QUAD
    std::cout<<"calc Ep Sp quad fugaku\n";
    CalcForceEpSpQuadFugaku f_ep_sp_fgk(EPISoft::eps*EPISoft::eps, ForceSoft::grav_const);
#else
    std::cout<<"calc Ep Sp mono fugaku\n";
    CalcForceEpSpMonoFugaku f_ep_sp_fgk(EPISoft::eps*EPISoft::eps, ForceSoft::grav_const);
#endif
    PS::F64 t_sp_fgk=0;
    t_sp_fgk -= PS::GetWtime();
    f_ep_sp_fgk(epi, Nepi, spj, Nspj, force_sp_fgk);
    t_sp_fgk += PS::GetWtime();

    std::cout<<"neighbor search fugaku\n";
    SearchNeighborEpEpFugaku f_nb_fgk;
    PS::F64 t_nb_fgk=0;
    t_nb_fgk -= PS::GetWtime();
    f_nb_fgk(epi, Nepi, epj, Nepj, force_nb_fgk);
    t_nb_fgk += PS::GetWtime();
#endif

    std::cout<<"calc Ep Ep\n";
    CalcForceEpEpWithLinearCutoffNoSimd f_ep_ep;
    PS::F64 t_ep_no=0;
    t_ep_no -= PS::GetWtime();
    f_ep_ep(epi, Nepi, epj, Nepj, force);
    t_ep_no += PS::GetWtime();

#ifdef USE_QUAD
    std::cout<<"calc Ep Sp quad\n";
    CalcForceEpSpQuadNoSimd f_ep_sp;
#else
    std::cout<<"calc Ep Sp mono\n";
    CalcForceEpSpMonoNoSimd f_ep_sp;
#endif
    PS::F64 t_sp_no=0;
    t_sp_no -= PS::GetWtime();
    f_ep_sp(epi, Nepi, spj, Nspj, force_sp);
    t_sp_no += PS::GetWtime();

    std::cout<<"neighbor search\n";
    SearchNeighborEpEpNoSimd f_nb;
    PS::F64 t_nb=0;
    t_nb -= PS::GetWtime();
    f_nb(epi, Nepi, epj, Nepj, force_nb);
    t_nb += PS::GetWtime();

    std::cout<<"compare results\n";
    PS::S32 nbcount[20];
    for(int i=0; i<20; i++) nbcount[i]=0;
    PS::F64 nbcount_ave=0;

#ifdef USE_SIMD
    PS::F64 dfmax_simd=0, dfpmax_simd=0;
    PS::F64 dsmax_simd=0, dspmax_simd=0;
    PS::F64 nbcount_ave_simd=0;
#endif
#ifdef USE_GPU
    PS::F64 dfmax_gpu=0, dfpmax_gpu=0;
    PS::F64 nbcount_ave_gpu=0;
#endif
#ifdef USE_FUGAKU
    PS::F64 dfmax_fgk=0,dfpmax_fgk=0;
    PS::F64 dsmax_fgk=0,dspmax_fgk=0;
    PS::F64 nbcount_ave_fgk=0;
#endif
    PS::F64 df;

    for(int i=0; i<Nepi; i++) {
        for (int j=0; j<3; j++) {
#ifdef USE_SIMD
            df=(force[i].acc[j]-force_simd[i].acc[j])/force[i].acc[j];
            dfmax_simd = std::max(dfmax_simd, df);
            if(df>DF_MAX) std::cerr<<"Force diff: i="<<i<<" nosimd["<<j<<"] "<<force[i].acc[j]<<" simd["<<j<<"] "<<force_simd[i].acc[j]<<std::endl;

            df = (force_sp[i].acc[j]-force_sp_simd[i].acc[j])/force_sp[i].acc[j];
            dsmax_simd = std::max(dsmax_simd, df);
            if(df>DF_MAX) std::cerr<<"Force sp diff: i="<<i<<" nosimd["<<j<<"] "<<force_sp[i].acc[j]<<" simd["<<j<<"] "<<force_sp_simd[i].acc[j]<<std::endl;
#endif
#ifdef USE_GPU
            dfmax_gpu = std::max(dfmax_gpu, df);
            if(df>DF_MAX) std::cerr<<"Force diff: i="<<i<<" nosimd["<<j<<"] "<<force[i].acc[j]+force_sp[i].acc[j]<<" gpu["<<j<<"] "<<force_gpu[i].acc[j]<<std::endl;
#endif
#ifdef USE_FUGAKU
            df=(force[i].acc[j]-force_fgk[i].acc[j])/force[i].acc[j];
            dfmax_fgk = std::max(dfmax_fgk, df);
            if(df>DF_MAX) std::cerr<<"Force diff: i="<<i<<" nosimd["<<j<<"] "<<force[i].acc[j]<<" fugaku["<<j<<"] "<<force_fgk[i].acc[j]<<std::endl;

            df=(force_sp[i].acc[j]-force_sp_fgk[i].acc[j])/force_sp[i].acc[j];
            dsmax_fgk = std::max(dsmax_fgk, df);
            if(df>DF_MAX) std::cerr<<"Force sp diff: i="<<i<<" nosimd["<<j<<"] "<<force_sp[i].acc[j]<<" fugaku["<<j<<"] "<<force_sp_fgk[i].acc[j]<<std::endl;
#endif
        }
#ifdef USE_SIMD
        dfpmax_simd = std::max(dfpmax_simd, (force[i].pot-force_simd[i].pot)/force[i].pot);
        dspmax_simd = std::max(dspmax_simd, (force_sp[i].pot-force_sp_simd[i].pot)/force_sp[i].pot);

        if(force[i].n_ngb!=force_simd[i].n_ngb) {
            std::cerr<<"Neighbor diff: i="<<i<<" nosimd "<<force[i].n_ngb<<" simd "<<force_simd[i].n_ngb<<std::endl;
        }
        if(force_nb[i].n_ngb!=force_nb_simd[i].n_ngb) {
            std::cerr<<"NB search diff: i="<<i<<" nosimd "<<force[i].n_ngb<<" simd "<<force_nb_simd[i].n_ngb<<std::endl;
        }
        nbcount_ave_simd += force_simd[i].n_ngb;
#endif
#ifdef USE_GPU
        dfpmax_gpu = std::max(dfpmax_gpu, (force_sp[i].pot+force[i].pot - force_gpu[i].pot)/force_gpu[i].pot);

        if(force[i].n_ngb!=force_gpu[i].n_ngb) {
            std::cerr<<"Neighbor diff: i="<<i<<" nosimd "<<force[i].n_ngb<<" gpu "<<force_gpu[i].n_ngb<<std::endl;
        }
        nbcount_ave_gpu += force_gpu[i].n_ngb;
#endif
#ifdef USE_FUGAKU
        dfpmax_fgk = std::max(dfpmax_fgk, (force[i].pot-force_fgk[i].pot)/force[i].pot);
        dspmax_fgk = std::max(dspmax_fgk, (force_sp[i].pot-force_sp_fgk[i].pot)/force_sp[i].pot);

        if(force[i].n_ngb!=force_fgk[i].n_ngb) {
            std::cerr<<"Neighbor diff: i="<<i<<" nosimd "<<force[i].n_ngb<<" fugaku "<<force_fgk[i].n_ngb<<std::endl;
        }
        if(force_nb[i].n_ngb!=force_nb_fgk[i].n_ngb) {
            std::cerr<<"NB search diff: i="<<i<<" nosimd "<<force_nb[i].n_ngb<<" fugaku "<<force_nb_fgk[i].n_ngb<<std::endl;
        }
        nbcount_ave_fgk += force_fgk[i].n_ngb;
#endif
        nbcount_ave += force[i].n_ngb;
        if (force[i].n_ngb<20) nbcount[force[i].n_ngb]++;
    }

#ifdef USE__AVX512
    std::cout<<"Use AVX512";
#elif defined(__AVX2__)
    std::cout<<"Use AVX2";
#elif defined(__AVX__)
    std::cout<<"Use AVX";
#endif
#ifdef USE_FUGAKU
#ifdef USE_QUAD
    std::cout<<" FUGAKU_quad";
#else
    std::cout<<" FUGAKU_mono";
#endif
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
    std::cout<<"SIMD EP-EP force diff max: "<<dfmax_simd<<" Pot diff max: "<<dfpmax_simd<<std::endl
             <<"SIMD EP-Sp force diff max: "<<dsmax_simd<<" Pot diff max: "<<dspmax_simd<<std::endl;
#endif
#ifdef USE_GPU
    std::cout<<"GPU EP+SP force diff max: "<<dfmax_gpu<<" Pot diff max: "<<dfpmax_gpu<<std::endl;
#endif
#ifdef USE_FUGAKU
    std::cout<<"Fugaku EP-EP diff max: "<<dfmax_fgk<<" Pot diff max: "<<dfpmax_fgk<<std::endl;
    std::cout<<"Fugaku EP-SP diff max: "<<dsmax_fgk<<" Pot diff max: "<<dspmax_fgk<<std::endl;
#endif

    for (int i=0; i<20; i++)
      if (nbcount[i]>0) std::cout<<"NNB: "<<i<<" "<<nbcount[i]<<std::endl;
    std::cout<<"<NNB>:";
    std::cout<<" no_simd: "<<nbcount_ave;
#ifdef USE_SIMD
    std::cout<<" simd: "<<nbcount_ave_simd;
#endif
#ifdef USE_GPU
    std::cout<<" gpu: "<<nbcount_ave_gpu;
#endif
#ifdef USE_FUGAKU
    std::cout<<" fugaku: "<<nbcount_ave_fgk;
#endif
    std::cout<<std::endl;
    
#ifdef USE_SIMD
    std::cout<<"Time: epj  simd="<<t_ep_simd<<" no="<<t_ep_no<<" ratio="<<t_ep_no/t_ep_simd<<std::endl;
    std::cout<<"Time: spj  simd="<<t_sp_simd<<" no="<<t_sp_no<<" ratio="<<t_sp_no/t_sp_simd<<std::endl;
#endif
#ifdef USE_GPU
    std::cout<<"Time: gpu ="<<t_gpu<<" no="<<t_ep_no+t_sp_no<<" ratio="<<(t_ep_no+t_sp_no)/t_gpu<<std::endl;
#endif
#ifdef USE_FUGAKU
    std::cout<<"Time: fugaku ="<<t_ep_fgk<<" no="<<t_ep_no<<" ratio="<<t_ep_no/t_ep_fgk<<std::endl;
    std::cout<<"Time: fugaku ="<<t_sp_fgk<<" no="<<t_sp_no<<" ratio="<<t_sp_no/t_sp_fgk<<std::endl;
#endif

    return 0;
}
