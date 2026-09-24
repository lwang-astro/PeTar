// PeTar NEON tree-force kernels for ARMv8 (Kunpeng-920 / TSV110), 128-bit NEON.
// API-compatible with the reference kernels in soft_force.hpp; enabled by
// ./configure --with-arch=tsv110 (-D USE_NEON_KERNEL). Validation and benchmark
// results: tsv110-bench/REPORT.md.

#pragma once
#ifdef USE_NEON_KERNEL

#include <arm_neon.h>
#include <vector>
#include"soft_ptcl.hpp"
#include"soft_force.hpp"

namespace tsv110{
// F32 must be the native float for the NEON intrinsics below (note that
// PS::F32 can be double under PARTICLE_SIMULATOR_ALL_64BIT_PRECISION)
using F32 = float;
using F64 = PS::F64;
using S32 = PS::S32;
using S64 = PS::S64;

//! fast reciprocal sqrt (4 x F32): estimate + cubic correction
//! r = r0 [1 + h (1/2 + 3h/8)], h = 1 - x r0^2 (relative error ~1e-7)
static inline float32x4_t rsqrt4(float32x4_t x){
    float32x4_t r = vrsqrteq_f32(x);
    float32x4_t h = vmulq_f32(x, r);
    h = vfmsq_f32(vdupq_n_f32(1.0f), h, r);
    float32x4_t p = vfmaq_n_f32(vdupq_n_f32(0.5f), h, 0.375f);
    p = vmulq_f32(p, h);
    return vfmaq_f32(r, r, p);
}
// Optional extra half Newton step for the quadrupole kernels, as used by the
// Fugaku implementation:  0 = cubic correction only (default)
//                         1 = cubic correction + one half Newton step
// Cost/accuracy study: tsv110-bench/quad_newton/REPORT.md
#ifndef NEON_QUAD_NEWTON
#define NEON_QUAD_NEWTON 0
#endif
static inline float32x4_t rsqrt4_quad(float32x4_t x){
    float32x4_t r = rsqrt4(x);
#if NEON_QUAD_NEWTON
    float32x4_t h = vmulq_f32(r, r);
    h = vfmsq_f32(vdupq_n_f32(3.0f), x, h);            // 3 - x r^2
    r = vmulq_f32(r, vmulq_f32(h, vdupq_n_f32(0.5f))); // r *= 0.5 h
#endif
    return r;
}

//! compact EPI layout (position + search radius, 16 B)
struct EPI32{
    F32 x, y, z, rs;
};
//! compact EPJ layout (position + mass + search radius, 20 B)
struct EPJ32{
    F32 x, y, z, m, rs;
};

//! per-thread scratch buffers, reused across kernel calls
struct Scratch{
    std::vector<S32>   lst;                            // original indices of active i
    std::vector<EPI32> iloc;                           // compact i (padded to multiple of 4)
    std::vector<EPJ32> jloc;                           // compact j (EP-EP force, mass>0 only)
    std::vector<F32>   jx, jy, jz, jm, jrs;            // compact j for I1_J4
    std::vector<F32>   qxx, qyy, qzz, qxy, qxz, qyz;   // compact super-particle quadrupole
};
//! return the thread-local scratch buffer
static inline Scratch& scratch(){
    static thread_local Scratch s;
    return s;
}

//! compact the active EPI list (type==1) into EPI32, padded to a multiple of 4
static inline void compact_i(const EPISoft* epi, const S32 ni, std::vector<S32>& lst,
                             std::vector<EPI32>& loc, const bool shift,
                             const F64 cx, const F64 cy, const F64 cz){
    lst.clear(); loc.clear();
    for(S32 i=0;i<ni;i++){
        if(epi[i].type!=1) continue;                   // same filter as the x86 SIMD/Fugaku kernels
        EPI32 p;
        p.x = (F32)(epi[i].pos.x - (shift ? cx : 0.0));
        p.y = (F32)(epi[i].pos.y - (shift ? cy : 0.0));
        p.z = (F32)(epi[i].pos.z - (shift ? cz : 0.0));
        p.rs= (F32)epi[i].r_search;
        lst.push_back(i); loc.push_back(p);
    }
    while(loc.size()%4) loc.push_back(EPI32{0.f,0.f,0.f,0.f});
}

//! compact the EPJ list for the EP-EP force, skipping zero-mass entries
//! (same behaviour as the x86 SIMD and Fugaku kernels)
static inline S32 compact_j_ep(const EPJSoft* epj, const S32 nj, std::vector<EPJ32>& loc){
    loc.clear();
    for(S32 j=0;j<nj;j++){
        if(epj[j].mass<=0.0) continue;
        EPJ32 p;
        p.x=(F32)epj[j].pos.x; p.y=(F32)epj[j].pos.y; p.z=(F32)epj[j].pos.z;
        p.m=(F32)epj[j].mass;  p.rs=(F32)epj[j].r_search;
        loc.push_back(p);
    }
    return (S32)loc.size();
}

//! kernel orientation selector: vectorize over i (I4_J1) or over j (I1_J4)
static inline bool use_i4(const S32 ni, const S32 nj){
    return (nj<=8) || (ni<=4) || ((S64)ni*(S64)nj<=512);
}
//! below this interaction count the scalar reference kernel is faster
static inline bool use_scalar(const S32 ni, const S32 nj){
    return ((S64)ni*(S64)nj < 256);
}

//! tree neighbor search (EP-EP)
struct SearchNeighborEpEpNeon{
    void operator()(const EPISoft* epi, const S32 ni, const EPJSoft* epj, const S32 nj, ForceSoft* force) const {
        if(use_scalar(ni,nj)){ SearchNeighborEpEpNoSimd()(epi,ni,epj,nj,force); return; }
        if(use_i4(ni,nj)) Kernel_I4_J1(epi,ni,epj,nj,force);
        else              Kernel_I1_J4(epi,ni,epj,nj,force);
    }
    void Kernel_I4_J1(const EPISoft* epi, const S32 ni, const EPJSoft* epj, const S32 nj, ForceSoft* force) const {
        Scratch& S = scratch();
        compact_i(epi,ni,S.lst,S.iloc,false,0,0,0);
        const S32 na=(S32)S.lst.size();
        for(S32 ib=0; ib<na; ib+=4){
            float32x4x4_t t = vld4q_f32((const F32*)&S.iloc[ib]);
            float32x4_t xi=t.val[0], yi=t.val[1], zi=t.val[2], rsi=t.val[3];
            float32x4_t rsi2=vmulq_f32(rsi,rsi);
            int32x4_t nn=vdupq_n_s32(0);
            for(S32 j=0;j<nj;j++){
                float32x4_t dx=vsubq_f32(xi,vdupq_n_f32((F32)epj[j].pos.x));
                float32x4_t dy=vsubq_f32(yi,vdupq_n_f32((F32)epj[j].pos.y));
                float32x4_t dz=vsubq_f32(zi,vdupq_n_f32((F32)epj[j].pos.z));
                float32x4_t r2=vfmaq_f32(vfmaq_f32(vmulq_f32(dx,dx),dy,dy),dz,dz);
                F32 rs2=(F32)(epj[j].r_search*epj[j].r_search);
                uint32x4_t cmp=vcltq_f32(r2,vmaxq_f32(rsi2,vdupq_n_f32(rs2)));
                nn=vsubq_s32(nn,vreinterpretq_s32_u32(cmp));
            }
            S32 tmp[4]; vst1q_s32(tmp,nn);
            for(S32 k=0;k<4;k++){
                S32 ii=ib+k; if(ii>=na) break;
                force[S.lst[ii]].n_ngb = tmp[k];
            }
        }
    }
    void Kernel_I1_J4(const EPISoft* epi, const S32 ni, const EPJSoft* epj, const S32 nj, ForceSoft* force) const {
        Scratch& S = scratch();
        compact_i(epi,ni,S.lst,S.iloc,false,0,0,0);
        const S32 na=(S32)S.lst.size();
        const S32 n4=(nj+3)/4*4;
        S.jx.resize(n4); S.jy.resize(n4); S.jz.resize(n4); S.jrs.resize(n4);
        for(S32 j=0;j<nj;j++){
            S.jx[j]=(F32)epj[j].pos.x; S.jy[j]=(F32)epj[j].pos.y; S.jz[j]=(F32)epj[j].pos.z;
            S.jrs[j]=(F32)epj[j].r_search;
        }
        for(S32 j=nj;j<n4;j++){ S.jx[j]=1e15f; S.jy[j]=1e15f; S.jz[j]=1e15f; S.jrs[j]=0.f; }
        for(S32 ii=0; ii<na; ii++){
            const F32 xi=S.iloc[ii].x, yi=S.iloc[ii].y, zi=S.iloc[ii].z, rsi=S.iloc[ii].rs;
            float32x4_t xv=vdupq_n_f32(xi), yv=vdupq_n_f32(yi), zv=vdupq_n_f32(zi);
            float32x4_t rsi2=vdupq_n_f32(rsi*rsi);
            int32x4_t nn=vdupq_n_s32(0);
            for(S32 j=0;j<n4;j+=4){
                float32x4_t dx=vsubq_f32(xv,vld1q_f32(&S.jx[j]));
                float32x4_t dy=vsubq_f32(yv,vld1q_f32(&S.jy[j]));
                float32x4_t dz=vsubq_f32(zv,vld1q_f32(&S.jz[j]));
                float32x4_t r2=vfmaq_f32(vfmaq_f32(vmulq_f32(dx,dx),dy,dy),dz,dz);
                float32x4_t rs=vld1q_f32(&S.jrs[j]);
                uint32x4_t cmp=vcltq_f32(r2,vmaxq_f32(rsi2,vmulq_f32(rs,rs)));
                nn=vsubq_s32(nn,vreinterpretq_s32_u32(cmp));
            }
            force[S.lst[ii]].n_ngb = vaddvq_s32(nn);
        }
    }
};

//! EP-EP force with linear cutoff and neighbor counting
struct CalcForceEpEpWithLinearCutoffNeon{
    F32 eps2, rcut2; F64 G;
    CalcForceEpEpWithLinearCutoffNeon(){}
    CalcForceEpEpWithLinearCutoffNeon(F64 e2, F64 rc2, F64 g):eps2((F32)e2),rcut2((F32)rc2),G(g){}
    void initialize(F64 e2, F64 rc2, F64 g){ eps2=(F32)e2; rcut2=(F32)rc2; G=g; }
    void operator()(const EPISoft* epi, const S32 ni, const EPJSoft* epj, const S32 nj, ForceSoft* force) const {
        if(use_scalar(ni,nj)){ CalcForceEpEpWithLinearCutoffNoSimd()(epi,ni,epj,nj,force); return; }
        if(use_i4(ni,nj)) Kernel_I4_J1(epi,ni,epj,nj,force);
        else              Kernel_I1_J4(epi,ni,epj,nj,force);
    }
    void Kernel_I4_J1(const EPISoft* epi, const S32 ni, const EPJSoft* epj, const S32 nj, ForceSoft* force) const {
        Scratch& S = scratch();
        compact_i(epi,ni,S.lst,S.iloc,false,0,0,0);
        const S32 na=(S32)S.lst.size();
        const S32 nja=compact_j_ep(epj,nj,S.jloc);
        const float32x4_t e2v=vdupq_n_f32(eps2), rcv=vdupq_n_f32(rcut2), one=vdupq_n_f32(1.f);
        for(S32 ib=0; ib<na; ib+=4){
            float32x4x4_t t = vld4q_f32((const F32*)&S.iloc[ib]);
            float32x4_t xi=t.val[0], yi=t.val[1], zi=t.val[2], rsi=t.val[3];
            float32x4_t rsi2=vmulq_f32(rsi,rsi);
            float32x4_t ax=vdupq_n_f32(0.f), ay=ax, az=ax, ap=ax;
            int32x4_t nn=vdupq_n_s32(0);
            for(S32 j=0;j<nja;j++){
                const EPJ32& q=S.jloc[j];
                F32 mj=q.m;
                float32x4_t dx=vsubq_f32(xi,vdupq_n_f32(q.x));
                float32x4_t dy=vsubq_f32(yi,vdupq_n_f32(q.y));
                float32x4_t dz=vsubq_f32(zi,vdupq_n_f32(q.z));
                float32x4_t r2=vfmaq_f32(vfmaq_f32(vmulq_f32(dx,dx),dy,dy),dz,dz);
                float32x4_t rc=vmaxq_f32(vaddq_f32(r2,e2v),rcv);
                float32x4_t ri=rsqrt4(rc);
                float32x4_t ri2=vmulq_f32(ri,ri);
                float32x4_t mr=vmulq_f32(vdupq_n_f32(mj),ri);
                float32x4_t mr3=vmulq_f32(ri2,mr);
                ax=vfmsq_f32(ax,mr3,dx);
                ay=vfmsq_f32(ay,mr3,dy);
                az=vfmsq_f32(az,mr3,dz);
                ap=vfmsq_f32(ap,mr,one);
                F32 rs2=q.rs*q.rs;
                uint32x4_t cmp=vcltq_f32(r2,vmaxq_f32(rsi2,vdupq_n_f32(rs2)));
                nn=vsubq_s32(nn,vreinterpretq_s32_u32(cmp));
            }
            F32 fx[4],fy[4],fz[4],fp[4]; S32 fn[4];
            vst1q_f32(fx,ax); vst1q_f32(fy,ay); vst1q_f32(fz,az); vst1q_f32(fp,ap); vst1q_s32(fn,nn);
            for(S32 k=0;k<4;k++){
                S32 ii=ib+k; if(ii>=na) break;
                S32 jj=S.lst[ii];
                force[jj].acc.x += G*fx[k];
                force[jj].acc.y += G*fy[k];
                force[jj].acc.z += G*fz[k];
                force[jj].pot   += G*fp[k];
                force[jj].n_ngb  = fn[k];
            }
        }
    }
    void Kernel_I1_J4(const EPISoft* epi, const S32 ni, const EPJSoft* epj, const S32 nj, ForceSoft* force) const {
        Scratch& S = scratch();
        compact_i(epi,ni,S.lst,S.iloc,false,0,0,0);
        const S32 na=(S32)S.lst.size();
        // single pass: filter mass>0 and pack SoA
        S.jx.resize(nj+3); S.jy.resize(nj+3); S.jz.resize(nj+3); S.jm.resize(nj+3); S.jrs.resize(nj+3);
        S32 nja=0;
        for(S32 j=0;j<nj;j++){
            if(epj[j].mass<=0.0) continue;
            S.jx[nja]=(F32)epj[j].pos.x; S.jy[nja]=(F32)epj[j].pos.y; S.jz[nja]=(F32)epj[j].pos.z;
            S.jm[nja]=(F32)epj[j].mass;  S.jrs[nja]=(F32)epj[j].r_search;
            nja++;
        }
        const S32 n4=(nja+3)/4*4;
        for(S32 j=nja;j<n4;j++){ S.jx[j]=1e15f; S.jy[j]=1e15f; S.jz[j]=1e15f; S.jm[j]=0.f; S.jrs[j]=0.f; }
        const float32x4_t e2v=vdupq_n_f32(eps2), rcv=vdupq_n_f32(rcut2), one=vdupq_n_f32(1.f);
        for(S32 ii=0; ii<na; ii++){
            const F32 xi=S.iloc[ii].x, yi=S.iloc[ii].y, zi=S.iloc[ii].z, rsi=S.iloc[ii].rs;
            float32x4_t xv=vdupq_n_f32(xi), yv=vdupq_n_f32(yi), zv=vdupq_n_f32(zi);
            float32x4_t rsi2=vdupq_n_f32(rsi*rsi);
            float32x4_t ax=vdupq_n_f32(0.f), ay=ax, az=ax, ap=ax;
            int32x4_t nn=vdupq_n_s32(0);
            for(S32 j=0;j<n4;j+=4){
                float32x4_t dx=vsubq_f32(xv,vld1q_f32(&S.jx[j]));
                float32x4_t dy=vsubq_f32(yv,vld1q_f32(&S.jy[j]));
                float32x4_t dz=vsubq_f32(zv,vld1q_f32(&S.jz[j]));
                float32x4_t r2=vfmaq_f32(vfmaq_f32(vmulq_f32(dx,dx),dy,dy),dz,dz);
                float32x4_t rc=vmaxq_f32(vaddq_f32(r2,e2v),rcv);
                float32x4_t ri=rsqrt4(rc);
                float32x4_t ri2=vmulq_f32(ri,ri);
                float32x4_t mr=vmulq_f32(vld1q_f32(&S.jm[j]),ri);
                float32x4_t mr3=vmulq_f32(ri2,mr);
                ax=vfmsq_f32(ax,mr3,dx);
                ay=vfmsq_f32(ay,mr3,dy);
                az=vfmsq_f32(az,mr3,dz);
                ap=vfmsq_f32(ap,mr,one);
                float32x4_t rs=vld1q_f32(&S.jrs[j]);
                uint32x4_t cmp=vcltq_f32(r2,vmaxq_f32(rsi2,vmulq_f32(rs,rs)));
                nn=vsubq_s32(nn,vreinterpretq_s32_u32(cmp));
            }
            S32 jj=S.lst[ii];
            force[jj].acc.x += G*vaddvq_f32(ax);
            force[jj].acc.y += G*vaddvq_f32(ay);
            force[jj].acc.z += G*vaddvq_f32(az);
            force[jj].pot   += G*vaddvq_f32(ap);
            force[jj].n_ngb  = vaddvq_s32(nn);
        }
    }
};

//! EP-SP monopole force
template<class Tsp>
struct CalcForceEpSpMonoNeon{
    F32 eps2; F64 G;
    CalcForceEpSpMonoNeon(){}
    CalcForceEpSpMonoNeon(F64 e2, F64 g):eps2((F32)e2),G(g){}
    void initialize(F64 e2, F64 g){ eps2=(F32)e2; G=g; }
    void operator()(const EPISoft* epi,const S32 ni,const Tsp* epj,const S32 nj,ForceSoft* force) const {
        if(use_scalar(ni,nj)){ CalcForceEpSpMonoNoSimd()(epi,ni,epj,nj,force); return; }
        if(use_i4(ni,nj)) Kernel_I4_J1(epi,ni,epj,nj,force);
        else              Kernel_I1_J4(epi,ni,epj,nj,force);
    }
    void Kernel_I4_J1(const EPISoft* epi,const S32 ni,const Tsp* epj,const S32 nj,ForceSoft* force) const {
        Scratch& S = scratch();
        const F64 cx=epi[0].pos.x, cy=epi[0].pos.y, cz=epi[0].pos.z;
        compact_i(epi,ni,S.lst,S.iloc,true,cx,cy,cz);
        const S32 na=(S32)S.lst.size();
        const float32x4_t e2v=vdupq_n_f32(eps2), one=vdupq_n_f32(1.f);
        for(S32 ib=0; ib<na; ib+=4){
            float32x4x4_t t = vld4q_f32((const F32*)&S.iloc[ib]);
            float32x4_t xi=t.val[0], yi=t.val[1], zi=t.val[2];
            float32x4_t ax=vdupq_n_f32(0.f), ay=ax, az=ax, ap=ax;
            for(S32 j=0;j<nj;j++){
                F32 mj=(F32)epj[j].mass;
                float32x4_t dx=vsubq_f32(xi,vdupq_n_f32((F32)(epj[j].pos.x-cx)));
                float32x4_t dy=vsubq_f32(yi,vdupq_n_f32((F32)(epj[j].pos.y-cy)));
                float32x4_t dz=vsubq_f32(zi,vdupq_n_f32((F32)(epj[j].pos.z-cz)));
                float32x4_t r2=vaddq_f32(vfmaq_f32(vfmaq_f32(vmulq_f32(dx,dx),dy,dy),dz,dz),e2v);
                float32x4_t ri=rsqrt4(r2);
                float32x4_t ri2=vmulq_f32(ri,ri);
                float32x4_t mr=vmulq_f32(vdupq_n_f32(mj),ri);
                float32x4_t mr3=vmulq_f32(ri2,mr);
                ax=vfmsq_f32(ax,mr3,dx);
                ay=vfmsq_f32(ay,mr3,dy);
                az=vfmsq_f32(az,mr3,dz);
                ap=vfmsq_f32(ap,mr,one);
            }
            F32 fx[4],fy[4],fz[4],fp[4];
            vst1q_f32(fx,ax); vst1q_f32(fy,ay); vst1q_f32(fz,az); vst1q_f32(fp,ap);
            for(S32 k=0;k<4;k++){
                S32 ii=ib+k; if(ii>=na) break;
                S32 jj=S.lst[ii];
                force[jj].acc.x += G*fx[k];
                force[jj].acc.y += G*fy[k];
                force[jj].acc.z += G*fz[k];
                force[jj].pot   += G*fp[k];
            }
        }
    }
    void Kernel_I1_J4(const EPISoft* epi,const S32 ni,const Tsp* epj,const S32 nj,ForceSoft* force) const {
        Scratch& S = scratch();
        const F64 cx=epi[0].pos.x, cy=epi[0].pos.y, cz=epi[0].pos.z;
        compact_i(epi,ni,S.lst,S.iloc,true,cx,cy,cz);
        const S32 na=(S32)S.lst.size();
        const S32 n4=(nj+3)/4*4;
        S.jx.resize(n4); S.jy.resize(n4); S.jz.resize(n4); S.jm.resize(n4);
        for(S32 j=0;j<nj;j++){
            S.jx[j]=(F32)(epj[j].pos.x-cx); S.jy[j]=(F32)(epj[j].pos.y-cy); S.jz[j]=(F32)(epj[j].pos.z-cz);
            S.jm[j]=(F32)epj[j].mass;
        }
        for(S32 j=nj;j<n4;j++){ S.jx[j]=1e15f; S.jy[j]=1e15f; S.jz[j]=1e15f; S.jm[j]=0.f; }
        const float32x4_t e2v=vdupq_n_f32(eps2), one=vdupq_n_f32(1.f);
        for(S32 ii=0; ii<na; ii++){
            const F32 xi=S.iloc[ii].x, yi=S.iloc[ii].y, zi=S.iloc[ii].z;
            float32x4_t xv=vdupq_n_f32(xi), yv=vdupq_n_f32(yi), zv=vdupq_n_f32(zi);
            float32x4_t ax=vdupq_n_f32(0.f), ay=ax, az=ax, ap=ax;
            for(S32 j=0;j<n4;j+=4){
                float32x4_t dx=vsubq_f32(xv,vld1q_f32(&S.jx[j]));
                float32x4_t dy=vsubq_f32(yv,vld1q_f32(&S.jy[j]));
                float32x4_t dz=vsubq_f32(zv,vld1q_f32(&S.jz[j]));
                float32x4_t r2=vaddq_f32(vfmaq_f32(vfmaq_f32(vmulq_f32(dx,dx),dy,dy),dz,dz),e2v);
                float32x4_t ri=rsqrt4(r2);
                float32x4_t ri2=vmulq_f32(ri,ri);
                float32x4_t mr=vmulq_f32(vld1q_f32(&S.jm[j]),ri);
                float32x4_t mr3=vmulq_f32(ri2,mr);
                ax=vfmsq_f32(ax,mr3,dx);
                ay=vfmsq_f32(ay,mr3,dy);
                az=vfmsq_f32(az,mr3,dz);
                ap=vfmsq_f32(ap,mr,one);
            }
            S32 jj=S.lst[ii];
            force[jj].acc.x += G*vaddvq_f32(ax);
            force[jj].acc.y += G*vaddvq_f32(ay);
            force[jj].acc.z += G*vaddvq_f32(az);
            force[jj].pot   += G*vaddvq_f32(ap);
        }
    }
};

//! EP-SP quadrupole force
template<class Tsp>
struct CalcForceEpSpQuadNeon{
    F32 eps2; F64 G;
    CalcForceEpSpQuadNeon(){}
    CalcForceEpSpQuadNeon(F64 e2, F64 g):eps2((F32)e2),G(g){}
    void initialize(F64 e2, F64 g){ eps2=(F32)e2; G=g; }
    void operator()(const EPISoft* epi,const S32 ni,const Tsp* epj,const S32 nj,ForceSoft* force) const {
        if(use_scalar(ni,nj)){ CalcForceEpSpQuadNoSimd()(epi,ni,epj,nj,force); return; }
        if(use_i4(ni,nj)) Kernel_I4_J1(epi,ni,epj,nj,force);
        else              Kernel_I1_J4(epi,ni,epj,nj,force);
    }
    void Kernel_I4_J1(const EPISoft* epi,const S32 ni,const Tsp* epj,const S32 nj,ForceSoft* force) const {
        Scratch& S = scratch();
        const F64 cx=epi[0].pos.x, cy=epi[0].pos.y, cz=epi[0].pos.z;
        compact_i(epi,ni,S.lst,S.iloc,true,cx,cy,cz);
        const S32 na=(S32)S.lst.size();
        const float32x4_t e2v=vdupq_n_f32(eps2), one=vdupq_n_f32(1.f), half=vdupq_n_f32(0.5f);
        for(S32 ib=0; ib<na; ib+=4){
            float32x4x4_t t = vld4q_f32((const F32*)&S.iloc[ib]);
            float32x4_t xi=t.val[0], yi=t.val[1], zi=t.val[2];
            float32x4_t ax=vdupq_n_f32(0.f), ay=ax, az=ax, ap=ax;
            for(S32 j=0;j<nj;j++){
                const Tsp& q = epj[j];
                F32 mj=(F32)q.mass;
                float32x4_t dx=vsubq_f32(xi,vdupq_n_f32((F32)(q.pos.x-cx)));
                float32x4_t dy=vsubq_f32(yi,vdupq_n_f32((F32)(q.pos.y-cy)));
                float32x4_t dz=vsubq_f32(zi,vdupq_n_f32((F32)(q.pos.z-cz)));
                float32x4_t r2=vaddq_f32(vfmaq_f32(vfmaq_f32(vmulq_f32(dx,dx),dy,dy),dz,dz),e2v);
                float32x4_t ri=rsqrt4_quad(r2);
                float32x4_t ri2=vmulq_f32(ri,ri);
                float32x4_t ri3=vmulq_f32(ri2,ri);
                float32x4_t r5=vmulq_f32(vmulq_f32(ri2,ri3),vdupq_n_f32(1.5f));
                float32x4_t qxx=vdupq_n_f32((F32)q.quad.xx), qyy=vdupq_n_f32((F32)q.quad.yy), qzz=vdupq_n_f32((F32)q.quad.zz);
                float32x4_t qxy=vdupq_n_f32((F32)q.quad.xy), qyz=vdupq_n_f32((F32)q.quad.yz), qxz=vdupq_n_f32((F32)q.quad.xz);
                float32x4_t qrx=vfmaq_f32(vfmaq_f32(vmulq_f32(qxx,dx),qxy,dy),qxz,dz);
                float32x4_t qry=vfmaq_f32(vfmaq_f32(vmulq_f32(qyy,dy),qyz,dz),qxy,dx);
                float32x4_t qrz=vfmaq_f32(vfmaq_f32(vmulq_f32(qzz,dz),qxz,dx),qyz,dy);
                float32x4_t qrr=vfmaq_f32(vfmaq_f32(vmulq_f32(qrx,dx),qry,dy),qrz,dz);
                float32x4_t tr=vaddq_f32(vaddq_f32(qxx,qyy),qzz);
                float32x4_t qrr5=vmulq_f32(r5,qrr);
                float32x4_t qrr7=vmulq_f32(ri2,qrr5);
                float32x4_t A=vfmsq_f32(vfmaq_f32(vmulq_f32(vdupq_n_f32(mj),ri3),qrr7,vdupq_n_f32(5.f)),tr,r5);
                float32x4_t m2r5=vmulq_f32(vdupq_n_f32(-2.f),r5);
                ax=vfmsq_f32(ax,A,dx); ay=vfmsq_f32(ay,A,dy); az=vfmsq_f32(az,A,dz);
                ax=vfmsq_f32(ax,m2r5,qrx); ay=vfmsq_f32(ay,m2r5,qry); az=vfmsq_f32(az,m2r5,qrz);
                ap=vfmsq_f32(ap,vmulq_f32(vdupq_n_f32(mj),ri),one);
                ap=vfmaq_f32(ap,vmulq_f32(half,vmulq_f32(tr,ri3)),one);
                ap=vfmsq_f32(ap,qrr5,one);
            }
            F32 fx[4],fy[4],fz[4],fp[4];
            vst1q_f32(fx,ax); vst1q_f32(fy,ay); vst1q_f32(fz,az); vst1q_f32(fp,ap);
            for(S32 k=0;k<4;k++){
                S32 ii=ib+k; if(ii>=na) break;
                S32 jj=S.lst[ii];
                force[jj].acc.x += G*fx[k];
                force[jj].acc.y += G*fy[k];
                force[jj].acc.z += G*fz[k];
                force[jj].pot   += G*fp[k];
            }
        }
    }
    void Kernel_I1_J4(const EPISoft* epi,const S32 ni,const Tsp* epj,const S32 nj,ForceSoft* force) const {
        Scratch& S = scratch();
        const F64 cx=epi[0].pos.x, cy=epi[0].pos.y, cz=epi[0].pos.z;
        compact_i(epi,ni,S.lst,S.iloc,true,cx,cy,cz);
        const S32 na=(S32)S.lst.size();
        const S32 n4=(nj+3)/4*4;
        S.jx.resize(n4); S.jy.resize(n4); S.jz.resize(n4); S.jm.resize(n4);
        S.qxx.resize(n4); S.qyy.resize(n4); S.qzz.resize(n4);
        S.qxy.resize(n4); S.qxz.resize(n4); S.qyz.resize(n4);
        for(S32 j=0;j<nj;j++){
            const Tsp& q=epj[j];
            S.jx[j]=(F32)(q.pos.x-cx); S.jy[j]=(F32)(q.pos.y-cy); S.jz[j]=(F32)(q.pos.z-cz);
            S.jm[j]=(F32)q.mass;
            S.qxx[j]=(F32)q.quad.xx; S.qyy[j]=(F32)q.quad.yy; S.qzz[j]=(F32)q.quad.zz;
            S.qxy[j]=(F32)q.quad.xy; S.qxz[j]=(F32)q.quad.xz; S.qyz[j]=(F32)q.quad.yz;
        }
        for(S32 j=nj;j<n4;j++){
            S.jx[j]=1e15f; S.jy[j]=1e15f; S.jz[j]=1e15f; S.jm[j]=0.f;
            S.qxx[j]=S.qyy[j]=S.qzz[j]=S.qxy[j]=S.qxz[j]=S.qyz[j]=0.f;
        }
        const float32x4_t e2v=vdupq_n_f32(eps2), one=vdupq_n_f32(1.f), half=vdupq_n_f32(0.5f);
        for(S32 ii=0; ii<na; ii++){
            const F32 xi=S.iloc[ii].x, yi=S.iloc[ii].y, zi=S.iloc[ii].z;
            float32x4_t xv=vdupq_n_f32(xi), yv=vdupq_n_f32(yi), zv=vdupq_n_f32(zi);
            float32x4_t ax=vdupq_n_f32(0.f), ay=ax, az=ax, ap=ax;
            for(S32 j=0;j<n4;j+=4){
                float32x4_t dx=vsubq_f32(xv,vld1q_f32(&S.jx[j]));
                float32x4_t dy=vsubq_f32(yv,vld1q_f32(&S.jy[j]));
                float32x4_t dz=vsubq_f32(zv,vld1q_f32(&S.jz[j]));
                float32x4_t r2=vaddq_f32(vfmaq_f32(vfmaq_f32(vmulq_f32(dx,dx),dy,dy),dz,dz),e2v);
                float32x4_t ri=rsqrt4_quad(r2);
                float32x4_t ri2=vmulq_f32(ri,ri);
                float32x4_t ri3=vmulq_f32(ri2,ri);
                float32x4_t r5=vmulq_f32(vmulq_f32(ri2,ri3),vdupq_n_f32(1.5f));
                float32x4_t Qxx=vld1q_f32(&S.qxx[j]), Qyy=vld1q_f32(&S.qyy[j]), Qzz=vld1q_f32(&S.qzz[j]);
                float32x4_t Qxy=vld1q_f32(&S.qxy[j]), Qxz=vld1q_f32(&S.qxz[j]), Qyz=vld1q_f32(&S.qyz[j]);
                float32x4_t qrx=vfmaq_f32(vfmaq_f32(vmulq_f32(Qxx,dx),Qxy,dy),Qxz,dz);
                float32x4_t qry=vfmaq_f32(vfmaq_f32(vmulq_f32(Qyy,dy),Qyz,dz),Qxy,dx);
                float32x4_t qrz=vfmaq_f32(vfmaq_f32(vmulq_f32(Qzz,dz),Qxz,dx),Qyz,dy);
                float32x4_t qrr=vfmaq_f32(vfmaq_f32(vmulq_f32(qrx,dx),qry,dy),qrz,dz);
                float32x4_t tr=vaddq_f32(vaddq_f32(Qxx,Qyy),Qzz);
                float32x4_t qrr5=vmulq_f32(r5,qrr);
                float32x4_t qrr7=vmulq_f32(ri2,qrr5);
                float32x4_t A=vfmsq_f32(vfmaq_f32(vmulq_f32(vld1q_f32(&S.jm[j]),ri3),qrr7,vdupq_n_f32(5.f)),tr,r5);
                float32x4_t m2r5=vmulq_f32(vdupq_n_f32(-2.f),r5);
                ax=vfmsq_f32(ax,A,dx); ay=vfmsq_f32(ay,A,dy); az=vfmsq_f32(az,A,dz);
                ax=vfmsq_f32(ax,m2r5,qrx); ay=vfmsq_f32(ay,m2r5,qry); az=vfmsq_f32(az,m2r5,qrz);
                ap=vfmsq_f32(ap,vmulq_f32(vld1q_f32(&S.jm[j]),ri),one);
                ap=vfmaq_f32(ap,vmulq_f32(half,vmulq_f32(tr,ri3)),one);
                ap=vfmsq_f32(ap,qrr5,one);
            }
            S32 jj=S.lst[ii];
            force[jj].acc.x += G*vaddvq_f32(ax);
            force[jj].acc.y += G*vaddvq_f32(ay);
            force[jj].acc.z += G*vaddvq_f32(az);
            force[jj].pot   += G*vaddvq_f32(ap);
        }
    }
};

} // namespace tsv110

#endif // USE_NEON_KERNEL
