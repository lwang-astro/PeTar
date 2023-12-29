#pragma once
#ifdef INTRINSIC_K
#include"phantomquad_for_p3t_k.hpp"
#endif
#ifdef INTRINSIC_X86
#include"phantomquad_for_p3t_x86.hpp"
#endif


// Neighbor search function
struct SearchNeighborEpEpNoSimd{
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
                      const EPJSoft * ep_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F64vec xi = ep_i[i].pos;
            PS::S32 n_ngb_i = 0;
            for(PS::S32 j=0; j<n_jp; j++){
                const PS::F64vec rij = xi - ep_j[j].pos;
                const PS::F64 r2 = rij * rij;
                const PS::F64 r_search = std::max(ep_i[i].r_search,ep_j[j].r_search);
                if(r2 < r_search*r_search){
#ifdef SAVE_NEIGHBOR_ID_IN_FORCE_KERNEL
                    force[i].id_ngb[n_ngb_i & 0x3] = ep_j[j].id;
#endif
                    n_ngb_i++;
                }
            }
            force[i].n_ngb = n_ngb_i;
        }
    }    
};

////////////////////
/// FORCE FUNCTOR
struct CalcForceEpEpWithLinearCutoffNoSimd{
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
                      const EPJSoft * ep_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
        const PS::F64 r_out2 = EPISoft::r_out*EPISoft::r_out;
        const PS::F64 G = ForceSoft::grav_const;
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F64vec xi = ep_i[i].pos;
            //PS::S64 id_i = ep_i[i].id;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            PS::S32 n_ngb_i = 0;
            for(PS::S32 j=0; j<n_jp; j++){
                //if(id_i == ep_j[j].id){
                //    n_ngb_i++;
                //    continue;
                //}
                const PS::F64vec rij = xi - ep_j[j].pos;
                const PS::F64 r2 = rij * rij;
                const PS::F64 r2_eps = r2 + eps2;
                const PS::F64 r_search = std::max(ep_i[i].r_search,ep_j[j].r_search);
                if(r2 < r_search*r_search){
                    n_ngb_i++;
                }
                const PS::F64 r2_tmp = (r2_eps > r_out2) ? r2_eps : r_out2;
                const PS::F64 r_inv = 1.0/sqrt(r2_tmp);
                const PS::F64 m_r = ep_j[j].mass * r_inv;
                const PS::F64 m_r3 = m_r * r_inv * r_inv;
                ai -= m_r3 * rij;
                poti -= m_r;
            }
            //std::cerr<<"poti= "<<poti<<std::endl;
            force[i].acc += G*ai;
#ifdef KDKDK_4TH
            force[i].acorr = 0.0;
#endif
            force[i].pot += G*poti;
#ifdef NAN_CHECK_DEBUG
            assert(!std::isnan(ai[0]));
            assert(!std::isnan(ai[1]));
            assert(!std::isnan(ai[2]));
            assert(!std::isnan(poti));
#endif
            force[i].n_ngb = n_ngb_i;
        }
    }
};

#ifdef KDKDK_4TH
struct CalcCorrectEpEpWithLinearCutoffNoSimd{
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
                      const EPJSoft * ep_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
        const PS::F64 r_out2 = EPISoft::r_out*EPISoft::r_out;

        for(PS::S32 i=0; i<n_ip; i++){
            PS::F64vec acorr = 0.0;
            const PS::F64vec posi = ep_i[i].pos;
            const PS::F64vec acci = ep_i[i].acc;
            for(PS::S32 j=0; j<n_jp; j++){
                const PS::F64vec dr = posi - ep_j[j].pos;
                const PS::F64vec da = acci - ep_j[j].acc; 
                const PS::F64 r2    = dr * dr + eps2;
                const PS::F64 drda  = dr * da;
                const PS::F64 r2_tmp = (r2 > r_out2) ? r2 : r_out2;
                const PS::F64 r_inv = 1.0/sqrt(r2_tmp);
                const PS::F64 r2_inv = r_inv*r_inv;
                const PS::F64 m_r = ep_j[j].mass * r_inv;
                const PS::F64 m_r3 = m_r * r2_inv;

                const PS::F64 alpha = 3.0 * drda * r2_inv;
                acorr -= m_r3 * (da - alpha * dr); 
            }
            //std::cerr<<"poti= "<<poti<<std::endl;
            force[i].acorr += 2.0 * acorr;
            force[i].acc = acci;
        }
    }
};
#endif

struct CalcForceEpSpMonoNoSimd {
    template<class Tsp>
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
                      const Tsp * sp_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
        const PS::F64 G = ForceSoft::grav_const;
        for(PS::S32 i=0; i<n_ip; i++){
            PS::F64vec xi = ep_i[i].pos;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            for(PS::S32 j=0; j<n_jp; j++){
                PS::F64vec rij = xi - sp_j[j].getPos();
                PS::F64 r3_inv = rij * rij + eps2;
                PS::F64 r_inv = 1.0/sqrt(r3_inv);
                r3_inv = r_inv * r_inv;
                r_inv *= sp_j[j].getCharge();
                r3_inv *= r_inv;
                ai -= r3_inv * rij;
                poti -= r_inv;
            }
            force[i].acc += G*ai;
            force[i].pot += G*poti;
#ifdef NAN_CHECK_DEBUG
            assert(!std::isnan(ai[0]));
            assert(!std::isnan(ai[1]));
            assert(!std::isnan(ai[2]));
            assert(!std::isnan(poti));
#endif
        }
    }
};

struct CalcForceEpSpQuadNoSimd{
    template<class Tsp>
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
                      const Tsp * sp_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
        const PS::F64 G = ForceSoft::grav_const;
//        assert(n_jp==0);
        for(PS::S32 ip=0; ip<n_ip; ip++){
            PS::F64vec xi = ep_i[ip].pos;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            for(PS::S32 jp=0; jp<n_jp; jp++){
                PS::F64 mj = sp_j[jp].mass;
                PS::F64vec xj= sp_j[jp].pos;
                PS::F64vec rij= xi - xj;
                PS::F64 r2 = rij * rij + eps2;
                PS::F64mat qj = sp_j[jp].quad;
                PS::F64 tr = qj.getTrace();
                PS::F64vec qr( (qj.xx*rij.x + qj.xy*rij.y + qj.xz*rij.z),
                               (qj.yy*rij.y + qj.yz*rij.z + qj.xy*rij.x),
                               (qj.zz*rij.z + qj.xz*rij.x + qj.yz*rij.y) );
                PS::F64 qrr = qr * rij;
                PS::F64 r_inv = 1.0f/sqrt(r2);
                PS::F64 r2_inv = r_inv * r_inv;
                PS::F64 r3_inv = r2_inv * r_inv;
                PS::F64 r5_inv = r2_inv * r3_inv * 1.5;
                PS::F64 qrr_r5 = r5_inv * qrr;
                PS::F64 qrr_r7 = r2_inv * qrr_r5;
                PS::F64 A = mj*r3_inv - tr*r5_inv + 5*qrr_r7;
                PS::F64 B = -2.0*r5_inv;
                ai -= A*rij + B*qr;
                poti -= mj*r_inv - 0.5*tr*r3_inv + qrr_r5;
            }
            force[ip].acc += G*ai;
            force[ip].pot += G*poti;
        }
    }
};

template<class Tpi, class Tpj>
struct CalcForcePPNoSimd {
    void operator () (const Tpi * ep_i,
                      const PS::S32 n_ip,
                      const Tpj * ep_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
      //const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
        const PS::F64 eps2 = 0;
        const PS::F64 G = ForceSoft::grav_const;
        for(PS::S32 i=0; i<n_ip; i++){
            PS::F64vec xi = ep_i[i].pos;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            for(PS::S32 j=0; j<n_jp; j++){
                PS::F64vec rij = xi - ep_j[j].pos;
                PS::F64 r3_inv = rij * rij + eps2;
                PS::F64 r_inv = 1.0/sqrt(r3_inv);
                r3_inv = r_inv * r_inv;
                r_inv *= ep_j[j].mass;
                r3_inv *= r_inv;
                ai -= r3_inv * rij;
                poti -= r_inv;
            }
            force[i].acc += G*ai;
            force[i].pot += G*poti;
#ifdef NAN_CHECK_DEBUG
            assert(!std::isnan(ai[0]));
            assert(!std::isnan(ai[1]));
            assert(!std::isnan(ai[2]));
            assert(!std::isnan(poti));
#endif
        }
    }
};

#ifdef USE_SIMD
struct SearchNeighborEpEpSimd{
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
                      const EPJSoft * ep_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
    #ifdef __HPC_ACE__
        PhantomGrapeQuad pg;
    #else
        #if defined(CALC_EP_64bit) || defined(CALC_EP_MIX)
        static thread_local PhantomGrapeQuad64Bit pg;
        #else
        static thread_local PhantomGrapeQuad pg;
        #endif
    #endif
        if(n_ip > pg.NIMAX || n_jp > pg.NJMAX){
            std::cout<<"ni= "<<n_ip<<" NIMAX= "<<pg.NIMAX<<" nj= "<<n_jp<<" NJMAX= "<<pg.NJMAX<<std::endl;
        }
        assert(n_ip<=pg.NIMAX);
        assert(n_jp<=pg.NJMAX);
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F64vec pos_i = ep_i[i].getPos();
            pg.set_xi_one(i, pos_i.x, pos_i.y, pos_i.z, ep_i[i].r_search);
        }
        PS::S32 loop_max = (n_jp-1) / PhantomGrapeQuad::NJMAX + 1;
        for(PS::S32 loop=0; loop<loop_max; loop++){
            const PS::S32 ih = PhantomGrapeQuad::NJMAX*loop;
            const PS::S32 n_jp_tmp = ( (n_jp - ih) < PhantomGrapeQuad::NJMAX) ? (n_jp - ih) : PhantomGrapeQuad::NJMAX;
            const PS::S32 it =ih + n_jp_tmp;
            PS::S32 i_tmp = 0;
            for(PS::S32 i=ih; i<it; i++, i_tmp++){
                const PS::F64 m_j = ep_j[i].getCharge();
                const PS::F64vec pos_j = ep_j[i].getPos();
                pg.set_epj_one(i_tmp, pos_j.x, pos_j.y, pos_j.z, m_j, ep_j[i].r_search);

            }
            pg.run_epj_for_neighbor_count(n_ip, n_jp_tmp);
            for(PS::S32 i=0; i<n_ip; i++){
                PS::F64 n_ngb = 0;
                pg.accum_accp_one(i, n_ngb);
                force[i].n_ngb += (PS::S32)(n_ngb*1.00001);
            }
        }
    }
};

template <class Tpi, class Tpj>
struct CalcForcePPSimd{
    void operator () (const Tpi * ep_i,
                      const PS::S32 n_ip,
                      const Tpj * ep_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 G = ForceSoft::grav_const;
        PS::S32 ep_j_list[n_jp], n_jp_local=0;
        for (PS::S32 i=0; i<n_jp; i++){
            if(ep_j[i].mass>0) ep_j_list[n_jp_local++] = i;
        }
    #ifdef __HPC_ACE__
        PhantomGrapeQuad pg;
    #else
        #if defined(CALC_EP_64bit) || defined(CALC_EP_MIX)
        static thread_local PhantomGrapeQuad64Bit pg;
        #else
        static thread_local PhantomGrapeQuad pg;
        #endif
    #endif
        if(n_ip > pg.NIMAX || n_jp > pg.NJMAX){
            std::cout<<"ni= "<<n_ip<<" NIMAX= "<<pg.NIMAX<<" nj= "<<n_jp<<" NJMAX= "<<pg.NJMAX<<std::endl;
        }
        assert(n_ip<=pg.NIMAX);
        assert(n_jp<=pg.NJMAX);
        pg.set_eps2(0.0);
        pg.set_r_crit2(0.0);
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F64vec pos_i = ep_i[i].pos;
            pg.set_xi_one(i, pos_i.x, pos_i.y, pos_i.z, 0.0);
        }
        PS::S32 loop_max = (n_jp_local-1) / PhantomGrapeQuad::NJMAX + 1;
        for(PS::S32 loop=0; loop<loop_max; loop++){
            const PS::S32 ih = PhantomGrapeQuad::NJMAX*loop;
            const PS::S32 n_jp_tmp = ( (n_jp_local - ih) < PhantomGrapeQuad::NJMAX) ? (n_jp_local - ih) : PhantomGrapeQuad::NJMAX;
            const PS::S32 it =ih + n_jp_tmp;
            PS::S32 i_tmp = 0;
            for(PS::S32 i=ih; i<it; i++, i_tmp++){
                const PS::S32 ij = ep_j_list[i];
                const PS::F64 m_j = ep_j[ij].mass;
                const PS::F64vec pos_j = ep_j[ij].pos;
                pg.set_epj_one(i_tmp, pos_j.x, pos_j.y, pos_j.z, m_j, 0.0);

            }
            pg.run_epj_for_p3t_with_linear_cutoff(n_ip, n_jp_tmp);
            for(PS::S32 i=0; i<n_ip; i++){
                PS::F64 p = 0;
                PS::F64 a[3]= {0,0,0};
                PS::F64 n_ngb = 0;
                pg.accum_accp_one(i, a[0], a[1], a[2], p, n_ngb);
                force[i].acc[0] += G*a[0];
                force[i].acc[1] += G*a[1];
                force[i].acc[2] += G*a[2];
                force[i].pot += G*p;
                force[i].n_ngb += (PS::S32)(n_ngb*1.00001);
            }
        }
    }
};

//! force calculation kernel for EP EP
/*! Notice this function cannot be used for neighbor searching because type of EPISoft is not correct (the group_data is used as cm). 
    Member particles will be excluded in the I particle list
 */
struct CalcForceEpEpWithLinearCutoffSimd{
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
                      const EPJSoft * ep_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
        const PS::F64 G = ForceSoft::grav_const;
        PS::S32 ep_j_list[n_jp], n_jp_local=0;
        PS::S32 ep_i_list[n_ip], n_ip_local=0;
        for (PS::S32 i=0; i<n_jp; i++){
            if(ep_j[i].mass>0) ep_j_list[n_jp_local++] = i;
        }
//        std::cerr<<"n_jp="<<n_jp<<" reduced n_jp="<<n_jp_local<<std::endl;
//        const PS::F64 r_crit2 = EPJSoft::r_search * EPJSoft::r_search;
    #ifdef __HPC_ACE__
        PhantomGrapeQuad pg;
    #else
        #if defined(CALC_EP_64bit) || defined(CALC_EP_MIX)
        static thread_local PhantomGrapeQuad64Bit pg;
        #else
        static thread_local PhantomGrapeQuad pg;
        #endif
    #endif
        if(n_ip > pg.NIMAX || n_jp > pg.NJMAX){
            std::cout<<"ni= "<<n_ip<<" NIMAX= "<<pg.NIMAX<<" nj= "<<n_jp<<" NJMAX= "<<pg.NJMAX<<std::endl;
        }
        assert(n_ip<=pg.NIMAX);
        assert(n_jp<=pg.NJMAX);
        pg.set_eps2(eps2);
        pg.set_r_crit2(EPISoft::r_out*EPISoft::r_out);
        for(PS::S32 i=0; i<n_ip; i++){
            // remove the orbital sample for the force calculation
            if (ep_i[i].type==1) {
                ep_i_list[n_ip_local] = i;
                const PS::F64vec pos_i = ep_i[i].getPos();
                pg.set_xi_one(n_ip_local, pos_i.x, pos_i.y, pos_i.z, ep_i[i].r_search);
                n_ip_local++;
            }
        }
//        std::cerr<<"n_ip="<<n_ip<<" reduced n_ip="<<n_ip_local<<std::endl;
        PS::S32 loop_max = (n_jp_local-1) / PhantomGrapeQuad::NJMAX + 1;
        for(PS::S32 loop=0; loop<loop_max; loop++){
            const PS::S32 ih = PhantomGrapeQuad::NJMAX*loop;
            const PS::S32 n_jp_tmp = ( (n_jp_local - ih) < PhantomGrapeQuad::NJMAX) ? (n_jp_local - ih) : PhantomGrapeQuad::NJMAX;
            const PS::S32 it =ih + n_jp_tmp;
            PS::S32 i_tmp = 0;
            for(PS::S32 i=ih; i<it; i++, i_tmp++){
                const PS::S32 ij = ep_j_list[i];
                const PS::F64 m_j = ep_j[ij].getCharge();
                const PS::F64vec pos_j = ep_j[ij].getPos();
                pg.set_epj_one(i_tmp, pos_j.x, pos_j.y, pos_j.z, m_j, ep_j[ij].r_search);

            }
            pg.run_epj_for_p3t_with_linear_cutoff(n_ip, n_jp_tmp);
            for(PS::S32 k=0; k<n_ip_local; k++){
                PS::S32 i=ep_i_list[k];
                PS::F64 p = 0;
                PS::F64 a[3]= {0,0,0};
                PS::F64 n_ngb = 0;
                pg.accum_accp_one(k, a[0], a[1], a[2], p, n_ngb);
#ifdef NAN_CHECK_DEBUG
                assert(!std::isnan(a[0]));
                assert(!std::isnan(a[1]));
                assert(!std::isnan(a[2]));
                assert(!std::isnan(p));
#endif
                force[i].acc[0] += G*a[0];
                force[i].acc[1] += G*a[1];
                force[i].acc[2] += G*a[2];
                force[i].pot += G*p;
                force[i].n_ngb += (PS::S32)(n_ngb*1.00001);
            }
        }
    }
};

struct CalcForceEpSpMonoSimd{
    template<class Tsp>
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
                      //const PS::SPJMonopoleScatter * sp_j,
                      const Tsp * sp_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
        const PS::F64 G = ForceSoft::grav_const;
        PS::S32 ep_i_list[n_ip], n_ip_local=0;
#ifdef __HPC_ACE__
        PhantomGrapeQuad pg;
#else
    #if defined(CALC_EP_64bit)
        static thread_local PhantomGrapeQuad64Bit pg;
    #else
        static thread_local PhantomGrapeQuad pg;
    #endif
#endif
        assert(n_ip<=pg.NIMAX);
        assert(n_jp<=pg.NJMAX);
        pg.set_eps2(eps2);
        for(PS::S32 i=0; i<n_ip; i++){
            // remove the orbital sample for the force calculation
            if (ep_i[i].type==1) {
                ep_i_list[n_ip_local] = i;
                const PS::F64vec pos_i = ep_i[i].getPos();
                pg.set_xi_one(n_ip_local, pos_i.x, pos_i.y, pos_i.z, 0.0);
                n_ip_local++;
            }                
        }
        PS::S32 loop_max = (n_jp-1) / PhantomGrapeQuad::NJMAX + 1;
        for(PS::S32 loop=0; loop<loop_max; loop++){
            const PS::S32 ih = PhantomGrapeQuad::NJMAX*loop;
            const PS::S32 n_jp_tmp = ( (n_jp - ih) < PhantomGrapeQuad::NJMAX) ? (n_jp - ih) : PhantomGrapeQuad::NJMAX;
            const PS::S32 it = ih + n_jp_tmp;
            PS::S32 i_tmp = 0;
            for(PS::S32 i=ih; i<it; i++, i_tmp++){
                const PS::F64 m_j = sp_j[i].getCharge();
                const PS::F64vec pos_j = sp_j[i].getPos();
                pg.set_epj_one(i_tmp, pos_j.x, pos_j.y, pos_j.z, m_j, 0.0);
            }
            pg.run_epj(n_ip, n_jp_tmp);
            for(PS::S32 k=0; k<n_ip_local; k++){
                PS::S32 i=ep_i_list[k];
                PS::F64 p = 0;
                PS::F64 a[3]= {0,0,0};
                pg.accum_accp_one(k, a[0], a[1], a[2], p);
                force[i].acc[0] += G*a[0];
                force[i].acc[1] += G*a[1];
                force[i].acc[2] += G*a[2];
                force[i].pot += G*p;
#ifdef NAN_CHECK_DEBUG
                assert(!std::isnan(a[0]));
                assert(!std::isnan(a[1]));
                assert(!std::isnan(a[2]));
                assert(!std::isnan(p));
#endif
            }
        }
    }
};

struct CalcForceEpSpQuadSimd{
    template<class Tsp>
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
                      const Tsp * sp_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
        const PS::F64 G = ForceSoft::grav_const;
        PS::S32 ep_i_list[n_ip], n_ip_local=0;
    #ifdef __HPC_ACE__
        PhantomGrapeQuad pg;
    #else
        #if defined(CALC_EP_64bit)
        static thread_local PhantomGrapeQuad64Bit pg;
        #else
        static thread_local PhantomGrapeQuad pg;
        #endif
    #endif
        assert(n_ip<=pg.NIMAX);
        assert(n_jp<=pg.NJMAX);
        pg.set_eps2(eps2);
        for(PS::S32 i=0; i<n_ip; i++){
            // remove the orbital sample for the force calculation
            if (ep_i[i].type==1) {
                ep_i_list[n_ip_local] = i;
                const PS::F64vec pos_i = ep_i[i].getPos();
                pg.set_xi_one(n_ip_local, pos_i.x, pos_i.y, pos_i.z, 0.0);
                n_ip_local++;
            }                
        }
        PS::S32 loop_max = (n_jp-1) / PhantomGrapeQuad::NJMAX + 1;
        for(PS::S32 loop=0; loop<loop_max; loop++){
            const PS::S32 ih = PhantomGrapeQuad::NJMAX*loop;
            const PS::S32 n_jp_tmp = ( (n_jp - ih) < PhantomGrapeQuad::NJMAX) ? (n_jp - ih) : PhantomGrapeQuad::NJMAX;
            const PS::S32 it = ih + n_jp_tmp;
            PS::S32 i_tmp = 0;
            for(PS::S32 i=ih; i<it; i++, i_tmp++){
                const PS::F64 m_j = sp_j[i].getCharge();
                const PS::F64vec pos_j = sp_j[i].getPos();
                const PS::F64mat q = sp_j[i].quad;
                pg.set_spj_one(i, pos_j.x, pos_j.y, pos_j.z, m_j,
                               q.xx, q.yy, q.zz, q.xy, q.yz, q.xz);
            }
            pg.run_spj(n_ip, n_jp_tmp);
            for(PS::S32 k=0; k<n_ip_local; k++){
                PS::S32 i=ep_i_list[k];
                PS::F64 p = 0;
                PS::F64 a[3]= {0,0,0};
                pg.accum_accp_one(k, a[0], a[1], a[2], p);
#ifdef NAN_CHECK_DEBUG
                assert(!std::isnan(a[0]));
                assert(!std::isnan(a[1]));
                assert(!std::isnan(a[2]));
                assert(!std::isnan(p));
#endif
                force[i].acc[0] += G*a[0];
                force[i].acc[1] += G*a[1];
                force[i].acc[2] += G*a[2];
                force[i].pot += G*p;
            }
        }
    }
};
#endif
