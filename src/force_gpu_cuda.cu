//#include "class.hpp"
//#include "force.hpp"
#include<particle_simulator.hpp>
#include "cuda_pointer.h"
#include "force_gpu_cuda.hpp"

#ifdef GPU_PROFILE
GPUProfile gpu_profile;
GPUCounter gpu_counter;
#endif

enum{
	N_THREAD_GPU = 32,
	N_WALK_LIMIT = 1000,
	NI_LIMIT     = N_WALK_LIMIT*1000,
	NJ_LIMIT     = N_WALK_LIMIT*10000,
};

struct EpiGPU{
	float3 pos;
    float  r_search;
	int    id_walk;
};

struct EpiDev{
	float3 pos;
    float  r_search;
};

struct EpjGPU{
	float3 pos;
    float  m;
    float  r_search;
};

struct SpjGPU{
    float3 pos;
    float  m;
#ifdef USE_QUAD
    float  qxx, qyy, qzz, qxy, qxz, qyz;
#endif
};

struct ForceGPU{
	float4 accp;
    int    nnb;
};

//! device pair force of Epi and Epi with linear cutoff
inline __device__ ForceGPU dev_gravity_ep_ep(
		float  eps2,
        float  rcut2,
        float  G,
		EpiDev epii,
		EpjGPU epjj,
        ForceGPU forcei)
{
	float dx = epjj.pos.x - epii.pos.x;
	float dy = epjj.pos.y - epii.pos.y;
	float dz = epjj.pos.z - epii.pos.z;

	float r2   = eps2 + dx*dx + dy*dy + dz*dz;
    float rsmin = max(epii.r_search, epjj.r_search);
    if (r2 < rsmin*rsmin) forcei.nnb ++;

    float r2_cut = (r2 > rcut2)? r2 : rcut2;
	float rinv = rsqrtf(r2_cut);
	float pij  = epjj.m * rinv;
	float mri3 = G*rinv*rinv * pij;

	forcei.accp.x += mri3 * dx;
	forcei.accp.y += mri3 * dy;
	forcei.accp.z += mri3 * dz;
	forcei.accp.w -= G*pij;

    return forcei;
}

__device__ ForceGPU force_kernel_ep_ep_1walk(
		EpjGPU       *jpsh,
		const EpiDev  epii,
		const int     id_walk,
		const int3   *ij_disp,
		const EpjGPU *epj, 
#ifdef PARTICLE_SIMULATOR_GPU_MULIT_WALK_INDEX
        const int*   *id_epj,
#endif
		ForceGPU      forcei,
		const float   eps2,
        const float   rcut2,
        const float   G) {

    const int tid = threadIdx.x;
    const int j_head = ij_disp[id_walk  ].y;
    const int j_tail = ij_disp[id_walk+1].y;

	for(int j=j_head; j<j_tail; j+=N_THREAD_GPU){
		// __syncthreads();
#ifdef PARTICLE_SIMULATOR_GPU_MULIT_WALK_INDEX
		jpsh[tid] = ((EpjGPU *)(epj + id_epj[j])) [tid];
#else
		jpsh[tid] = ((EpjGPU *)(epj + j)) [tid];
#endif
		// __syncthreads();

		if(j_tail-j < N_THREAD_GPU){
			for(int jj=0; jj<j_tail-j; jj++){
				forcei = dev_gravity_ep_ep(eps2, rcut2, G, epii, jpsh[jj], forcei);
			}
		}else{
#pragma unroll
			for(int jj=0; jj<N_THREAD_GPU; jj++){
				forcei = dev_gravity_ep_ep(eps2, rcut2, G, epii, jpsh[jj], forcei);
			}
		}
	}
	
	return forcei;
}

__device__ ForceGPU force_kernel_ep_ep_2walk(
		EpjGPU        jpsh[2][N_THREAD_GPU],
		const EpiDev  epii,
		const int     id_walk,
		const int     iwalk0,
		const int     iwalk1,
		const int3   *ij_disp,
		const EpjGPU *epj, 
#ifdef PARTICLE_SIMULATOR_GPU_MULIT_WALK_INDEX
        const int*   *id_epj,
#endif
		ForceGPU      forcei,
		const float   eps2,
        const float   rcut2,
        const float   G) {

	const int jbeg0 = ij_disp[iwalk0].y;
	const int jbeg1 = ij_disp[iwalk1].y;
	const int jend0 = ij_disp[iwalk0 + 1].y;
	const int jend1 = ij_disp[iwalk1 + 1].y;
	const int nj0   = jend0 - jbeg0;
	const int nj1   = jend1 - jbeg1;

	const int nj_longer  = nj0 > nj1 ? nj0 : nj1;
	const int nj_shorter = nj0 > nj1 ? nj1 : nj0;
	const int walk_longer= nj0 > nj1 ? 0 : 1;
	const int jbeg_longer = nj0 > nj1 ? jbeg0 : jbeg1;

	const int mywalk = id_walk==iwalk0 ? 0 : 1;

    const int tid = threadIdx.x;
	for(int j=0; j<nj_shorter; j+=N_THREAD_GPU){
#ifdef PARTICLE_SIMULATOR_GPU_MULIT_WALK_INDEX
		jpsh[0][tid] = ((EpjGPU *)(epj + id_epj[jbeg0 + j])) [tid];
		jpsh[1][tid] = ((EpjGPU *)(epj + id_epj[jbeg1 + j])) [tid];
#else
		jpsh[0][tid] = ((EpjGPU *)(epj + jbeg0 + j)) [tid];
		jpsh[1][tid] = ((EpjGPU *)(epj + jbeg1 + j)) [tid];
#endif
		if(nj_shorter-j < N_THREAD_GPU){
			for(int jj=0; jj<nj_shorter-j; jj++){
				forcei = dev_gravity_ep_ep(eps2, rcut2, G, epii, jpsh[mywalk][jj], forcei);
			}
		}else{
#pragma unroll
			for(int jj=0; jj<N_THREAD_GPU; jj++){
				forcei = dev_gravity_ep_ep(eps2, rcut2, G, epii, jpsh[mywalk][jj], forcei);
			}
		}
	}
	for(int j=nj_shorter; j<nj_longer; j+=N_THREAD_GPU){
		jpsh[0][tid] = ((EpjGPU *)(epj + jbeg_longer +  j)) [tid];
		int jrem = nj_longer - j;
		if(jrem < N_THREAD_GPU){
			for(int jj=0; jj<jrem; jj++){
				if(mywalk == walk_longer)
                    forcei = dev_gravity_ep_ep(eps2, rcut2, G, epii, jpsh[0][jj], forcei);
			}
		}else{
#pragma unroll
			for(int jj=0; jj<N_THREAD_GPU; jj++){
				if(mywalk == walk_longer)
                    forcei = dev_gravity_ep_ep(eps2, rcut2, G, epii, jpsh[0][jj], forcei);
			}
		}
	}

	return forcei;
}

__device__ ForceGPU force_kernel_ep_ep_multiwalk(
		const EpiDev  epii,
		const int     id_walk,
		const int3   *ij_disp,
		const EpjGPU *epj, 
#ifdef PARTICLE_SIMULATOR_GPU_MULIT_WALK_INDEX
        const int*   *id_epj,
#endif
		ForceGPU      forcei,
		const float   eps2,
        const float   rcut2,
        const float   G) {

    const int j_head = ij_disp[id_walk  ].y;
    const int j_tail = ij_disp[id_walk+1].y;

    for(int j=j_head; j<j_tail; j++){
#ifdef PARTICLE_SIMULATOR_GPU_MULIT_WALK_INDEX
		EpjGPU epjj = epj[id_epj[j]];
#else
		EpjGPU epjj = epj[j];
#endif
		forcei = dev_gravity_ep_ep(eps2, rcut2, G, epii, epjj, forcei);
	}
	return forcei;
}

__global__ void force_kernel_ep_ep(
		const int3   * ij_disp,
		const EpiGPU * epi,
		const EpjGPU * epj, 
#ifdef PARTICLE_SIMULATOR_GPU_MULIT_WALK_INDEX
        const int*   *id_epj,
#endif
		ForceGPU     * force,
		const float    eps2,
        const float    rcut2,
        const float    G) {

    int tid = blockDim.x * blockIdx.x + threadIdx.x;
	EpiDev epii;
    epii.pos       = epi[tid].pos;
    epii.r_search  = epi[tid].r_search;
	int    id_walk = epi[tid].id_walk;
	ForceGPU forcei;
    forcei.accp = make_float4(0.f, 0.f, 0.f, 0.f);
    forcei.nnb  = 0;

	int t_head = blockDim.x * blockIdx.x;
	int t_tail = t_head + N_THREAD_GPU - 1;
	int nwalk_in_block = 1 + (epi[t_tail].id_walk - epi[t_head].id_walk);

	__shared__ EpjGPU jpsh[2][N_THREAD_GPU];

#ifdef PARTICLE_SIMULATOR_GPU_MULIT_WALK_INDEX
	if(1 == nwalk_in_block){
		forcei = force_kernel_ep_ep_1walk(jpsh[0], epii, id_walk, ij_disp, epj, id_epj, forcei, eps2, rcut2, G);
	} else if(2 == nwalk_in_block){
		// accp = force_kernel_ep_ep_multiwalk(epii, id_walk, ij_disp, epj, accp, eps2);
		int iwalk0 = epi[t_head].id_walk;
		int iwalk1 = epi[t_tail].id_walk;
		forcei = force_kernel_ep_ep_2walk(jpsh, epii, id_walk, iwalk0, iwalk1, ij_disp, epj, id_epj, forcei, eps2, rcut2, G);
	} else{
		forcei = force_kernel_ep_ep_multiwalk(epii, id_walk, ij_disp, epj, id_epj, forcei, eps2, rcut2, G);
	}
#else
	if(1 == nwalk_in_block){
		forcei = force_kernel_ep_ep_1walk(jpsh[0], epii, id_walk, ij_disp, epj, forcei, eps2, rcut2, G);
	} else if(2 == nwalk_in_block){
		// accp = force_kernel_ep_ep_multiwalk(epii, id_walk, ij_disp, epj, accp, eps2);
		int iwalk0 = epi[t_head].id_walk;
		int iwalk1 = epi[t_tail].id_walk;
		forcei = force_kernel_ep_ep_2walk(jpsh, epii, id_walk, iwalk0, iwalk1, ij_disp, epj, forcei, eps2, rcut2, G);
	} else{
		forcei = force_kernel_ep_ep_multiwalk(epii, id_walk, ij_disp, epj, forcei, eps2, rcut2, G);
	}
#endif
	force[tid] = forcei;
}

//! device pair force of Epi and Spi 
inline __device__ float4 dev_gravity_ep_sp(
		float  eps2,
        float  G,
		float3 posi,
		SpjGPU spjj,
        float4 accpi) {

	float dx = posi.x - spjj.pos.x;
	float dy = posi.y - spjj.pos.y;
	float dz = posi.z - spjj.pos.z;

	float r2   = eps2 + dx*dx + dy*dy + dz*dz;
	float rinv = rsqrtf(r2);

#ifdef USE_QUAD
    float qrx = spjj.qxx*dx + spjj.qxy*dy + spjj.qxz*dz;
    float qry = spjj.qxy*dx + spjj.qyy*dy + spjj.qyz*dz;
    float qrz = spjj.qxz*dx + spjj.qyz*dy + spjj.qzz*dz;
    float tr = spjj.qxx + spjj.qyy + spjj.qzz;
    
    float qrr = qrx*dx + qry*dy + qrz*dz;
    float rinv2 = rinv*rinv;
    float rinv3 = rinv2*rinv;
    float rinv5 = rinv2*rinv3*1.5f;
    float qrr_r5 = rinv5*qrr;
    float qrr_r7 = rinv2*qrr_r5;
    float A = G*(spjj.m*rinv3 - tr*rinv5 + 5.0f*qrr_r7);
    float B = -2.0f*G*rinv5;
    
    accpi.x -= A*dx + B*qrx;
    accpi.y -= A*dy + B*qry;
    accpi.z -= A*dz + B*qrz;
    accpi.w -= G*(spjj.m*rinv - 0.5f*tr*rinv3 + qrr_r5);

#else
    
	float pij  = spjj.m * rinv;
	float mri3 = G*rinv*rinv * pij;

	accpi.x += mri3 * dx;
	accpi.y += mri3 * dy;
	accpi.z += mri3 * dz;
	accpi.w -= G*pij;

#endif

    return accpi;
}

__device__ float4 force_kernel_ep_sp_1walk(
		SpjGPU   *jpsh,
		const float3  posi,
		const int     id_walk,
		const int3   *ij_disp,
		const SpjGPU *spj, 
#ifdef PARTICLE_SIMULATOR_GPU_MULIT_WALK_INDEX
        const int*   *id_spj,
#endif
		float4        accpi,
		const float   eps2,
        const float   G) {

    const int tid = threadIdx.x;
    const int j_head = ij_disp[id_walk  ].z;
    const int j_tail = ij_disp[id_walk+1].z;

	for(int j=j_head; j<j_tail; j+=N_THREAD_GPU){
		// __syncthreads();
#ifdef PARTICLE_SIMULATOR_GPU_MULIT_WALK_INDEX
		jpsh[tid] = ((SpjGPU *)(spj + id_spj[j])) [tid];
#else
		jpsh[tid] = ((SpjGPU *)(spj + j)) [tid];
#endif
		// __syncthreads();

		if(j_tail-j < N_THREAD_GPU){
			for(int jj=0; jj<j_tail-j; jj++){
				accpi = dev_gravity_ep_sp(eps2, G, posi, jpsh[jj], accpi);
			}
		}else{
#pragma unroll
			for(int jj=0; jj<N_THREAD_GPU; jj++){
				accpi = dev_gravity_ep_sp(eps2, G, posi, jpsh[jj], accpi);
			}
		}
	}
	
	return accpi;
}

__device__ float4 force_kernel_ep_sp_2walk(
		SpjGPU        jpsh[2][N_THREAD_GPU],
		const float3  posi,
		const int     id_walk,
		const int     iwalk0,
		const int     iwalk1,
		const int3   *ij_disp,
		const SpjGPU *spj, 
#ifdef PARTICLE_SIMULATOR_GPU_MULIT_WALK_INDEX
        const int*   *id_spj,
#endif
		float4        accpi,
		const float   eps2,
        const float   G)
{
	const int jbeg0 = ij_disp[iwalk0].z;
	const int jbeg1 = ij_disp[iwalk1].z;
	const int jend0 = ij_disp[iwalk0 + 1].z;
	const int jend1 = ij_disp[iwalk1 + 1].z;
	const int nj0   = jend0 - jbeg0;
	const int nj1   = jend1 - jbeg1;

	const int nj_longer  = nj0 > nj1 ? nj0 : nj1;
	const int nj_shorter = nj0 > nj1 ? nj1 : nj0;
	const int walk_longer= nj0 > nj1 ? 0 : 1;
	const int jbeg_longer = nj0 > nj1 ? jbeg0 : jbeg1;

	const int mywalk = id_walk==iwalk0 ? 0 : 1;

    const int tid = threadIdx.x;
	for(int j=0; j<nj_shorter; j+=N_THREAD_GPU){
#ifdef PARTICLE_SIMULATOR_GPU_MULIT_WALK_INDEX
		jpsh[0][tid] = ((SpjGPU *)(spj + id_spj[jbeg0 + j])) [tid];
		jpsh[1][tid] = ((SpjGPU *)(spj + id_spj[jbeg1 + j])) [tid];
#else
		jpsh[0][tid] = ((SpjGPU *)(spj + jbeg0 + j)) [tid];
		jpsh[1][tid] = ((SpjGPU *)(spj + jbeg1 + j)) [tid];
#endif
		if(nj_shorter-j < N_THREAD_GPU){
			for(int jj=0; jj<nj_shorter-j; jj++){
				accpi = dev_gravity_ep_sp(eps2, G, posi, jpsh[mywalk][jj], accpi);
			}
		}else{
#pragma unroll
			for(int jj=0; jj<N_THREAD_GPU; jj++){
				accpi = dev_gravity_ep_sp(eps2, G, posi, jpsh[mywalk][jj], accpi);
			}
		}
	}
	for(int j=nj_shorter; j<nj_longer; j+=N_THREAD_GPU){
		jpsh[0][tid] = ((SpjGPU *)(spj + jbeg_longer +  j)) [tid];
		int jrem = nj_longer - j;
		if(jrem < N_THREAD_GPU){
			for(int jj=0; jj<jrem; jj++){
				if(mywalk == walk_longer)
                    accpi = dev_gravity_ep_sp(eps2, G, posi, jpsh[0][jj], accpi);
			}
		}else{
#pragma unroll
			for(int jj=0; jj<N_THREAD_GPU; jj++){
				if(mywalk == walk_longer)
                    accpi = dev_gravity_ep_sp(eps2, G, posi, jpsh[0][jj], accpi);
			}
		}
	}

	return accpi;
}

__device__ float4 force_kernel_ep_sp_multiwalk(
		const float3  posi,
		const int     id_walk,
		const int3   *ij_disp,
		const SpjGPU *spj, 
		float4        accpi,
		const float   eps2,
        const float   G)
{
    const int j_head = ij_disp[id_walk  ].z;
    const int j_tail = ij_disp[id_walk+1].z;

    for(int j=j_head; j<j_tail; j++){
		SpjGPU spjj = spj[j];
		accpi = dev_gravity_ep_sp(eps2, G, posi, spjj, accpi);
	}
	return accpi;
}

__global__ void force_kernel_ep_sp(
		const int3   * ij_disp,
		const EpiGPU * epi,
		const SpjGPU * spj, 
#ifdef PARTICLE_SIMULATOR_GPU_MULIT_WALK_INDEX
        const int*   *id_spj,
#endif
		ForceGPU     * force,
		const float    eps2,
        const float    G) {

    int tid = blockDim.x * blockIdx.x + threadIdx.x;
	float3 posi    = epi[tid].pos;
	int    id_walk = epi[tid].id_walk;
	float4 accpi   = force[tid].accp;

	int t_head = blockDim.x * blockIdx.x;
	int t_tail = t_head + N_THREAD_GPU - 1;
	int nwalk_in_block = 1 + (epi[t_tail].id_walk - epi[t_head].id_walk);

	__shared__ SpjGPU jpsh[2][N_THREAD_GPU];

#ifdef PARTICLE_SIMULATOR_GPU_MULIT_WALK_INDEX
	if(1 == nwalk_in_block){
		accpi = force_kernel_ep_sp_1walk(jpsh[0], posi, id_walk, ij_disp, spj, id_spj, accpi, eps2, G);
	} else if(2 == nwalk_in_block){
		// accpi = force_kernel_ep_sp_multiwalk(posi, id_walk, ij_disp, spj, accpi, eps2);
		int iwalk0 = epi[t_head].id_walk;
		int iwalk1 = epi[t_tail].id_walk;
		accpi = force_kernel_ep_sp_2walk(jpsh, posi, id_walk, iwalk0, iwalk1, ij_disp, spj, id_spj, accpi, eps2, G);
	} else{
		accpi = force_kernel_ep_sp_multiwalk(posi, id_walk, ij_disp, spj, id_spj, accpi, eps2, G);
	}
#else
	if(1 == nwalk_in_block){
		accpi = force_kernel_ep_sp_1walk(jpsh[0], posi, id_walk, ij_disp, spj, accpi, eps2, G);
	} else if(2 == nwalk_in_block){
		// accpi = force_kernel_ep_sp_multiwalk(posi, id_walk, ij_disp, spj, accpi, eps2);
		int iwalk0 = epi[t_head].id_walk;
		int iwalk1 = epi[t_tail].id_walk;
		accpi = force_kernel_ep_sp_2walk(jpsh, posi, id_walk, iwalk0, iwalk1, ij_disp, spj, accpi, eps2, G);
	} else{
		accpi = force_kernel_ep_sp_multiwalk(posi, id_walk, ij_disp, spj, accpi, eps2, G);
    }
#endif
	force[tid].accp = accpi;
}

static cudaPointer<EpiGPU>    dev_epi;
static cudaPointer<EpjGPU>    dev_epj;
static cudaPointer<SpjGPU>    dev_spj;
static cudaPointer<ForceGPU>  dev_force;
static cudaPointer<int3>      ij_disp;

static bool init_call = true;
#ifdef GPU_PROFILE
static cudaEvent_t cu_event_disp;
static cudaEvent_t cu_event_htod;
static cudaEvent_t cu_event_calc;
static cudaEvent_t cu_event_retr;
static cudaEvent_t cu_event_dtoh;
#endif

#ifdef PARTICLE_SIMULATOR_GPU_MULIT_WALK_INDEX

static cudaPointer<int>      dev_id_epj;
static cudaPointer<int>      dev_id_spj;

PS::S32 DispatchKernelWithSPIndex(const PS::S32 tag,
                                  const PS::S32 n_walk,
                                  const EPISoft ** epi,
                                  const PS::S32 *  n_epi,
                                  const PS::S32 ** id_epj,
                                  const PS::S32 *  n_epj,
                                  const PS::S32 ** id_spj,
                                  const PS::S32 *  n_spj,
                                  const EPJSoft * epj,
                                  const PS::S32 n_epj_tot,
                                  const SPJSoft * spj,
                                  const PS::S32 n_spj_tot,
                                  const bool send_flag) {
    assert(n_walk <= N_WALK_LIMIT);
    if(init_call){
		dev_epi  .allocate(NI_LIMIT);
		dev_epj  .allocate(NJ_LIMIT);
        dev_spj  .allocate(NJ_LIMIT);
		dev_force.allocate(NI_LIMIT);
		ij_disp  .allocate(N_WALK_LIMIT+2);
        dev_id_epj .allocate(NJ_LIMIT);
        dev_id_spj .allocate(NJ_LIMIT);
        cudaEventCreate(&cu_event_disp);
        cudaEventCreate(&cu_event_htod);
        cudaEventCreate(&cu_event_calc);
        cudaEventCreate(&cu_event_retr);
        cudaEventCreate(&cu_event_dtoh);
		init_call = false;
    }

    const float eps2 = EPISoft::eps * EPISoft::eps;
    const PS::F64 rcut2 = EPISoft::r_out*EPISoft::r_out;
    const PS::F64 G = ForceSoft::grav_const;

    if(send_flag==true){
#ifdef GPU_PROFILE
        gpu_profile.copy.start();
#endif
        /*
        if(dev_epj.size < n_epj_tot+n_spj_tot){
            dev_epj.free();
            dev_epj.allocate(n_epj_tot+n_spj_tot);
        }
        */
#pragma omp parallel for
        for(PS::S32 i=0; i<n_epj_tot; i++){
            dev_epj[i].pos.x  = epj[i].pos.x;
            dev_epj[i].pos.y  = epj[i].pos.y;
            dev_epj[i].pos.z  = epj[i].pos.z;
            dev_epj[i].m      = epj[i].mass;
            dev_epj[i].r_search = epj[i].r_search;
        }
#pragma omp parallel for
        for(PS::S32 i=0; i<n_spj_tot; i++){
            dev_spj[i].pos.x  = spj[i].pos.x;
            dev_spj[i].pos.y  = spj[i].pos.y;
            dev_spj[i].pos.z  = spj[i].pos.z;
            dev_spj[i].m      = spj[i].mass;
#ifdef USE_QUAD
            dev_spj[i].qxx    = spj[i].quad.xx;
            dev_spj[i].qyy    = spj[i].quad.yy;
            dev_spj[i].qzz    = spj[i].quad.zz;
            dev_spj[i].qxy    = spj[i].quad.xy;
            dev_spj[i].qxz    = spj[i].quad.xz;
            dev_spj[i].qyz    = spj[i].quad.yz;
#endif
        }
#ifdef GPU_PROFILE
        gpu_profile.copy.end();
        cudaEventRecord(cu_event_disp);
#endif
        dev_epj.htod(n_epj_tot);
        dev_spj.htod(n_spj_tot);
#ifdef GPU_PROFILE
        cudaEventRecord(cu_event_htod);
#endif
        return 0;
    }
    else{
#ifdef GPU_PROFILE
        gpu_profile.copy.start();
#endif
        ij_disp[0].x = 0;
        ij_disp[0].y = 0;
        ij_disp[0].z = 0;
        for(int k=0; k<n_walk; k++){
            ij_disp[k+1].x = ij_disp[k].x + n_epi[k];
            ij_disp[k+1].y = ij_disp[k].y + n_epj[k];
            ij_disp[k+1].z = ij_disp[k].z + n_spj[k];
        }
        ij_disp[n_walk+1] = ij_disp[n_walk];
        assert(ij_disp[n_walk].x < NI_LIMIT);
        assert(ij_disp[n_walk].y < NJ_LIMIT);
        assert(ij_disp[n_walk].z < NJ_LIMIT);

        int ni_tot_reg = ij_disp[n_walk].x;
        if(ni_tot_reg % N_THREAD_GPU){
            ni_tot_reg /= N_THREAD_GPU;
            ni_tot_reg++;
            ni_tot_reg *= N_THREAD_GPU;
        }

        int ni_tot = ij_disp[n_walk].x;
        int nej_tot = ij_disp[n_walk].y;
        int nsj_tot = ij_disp[n_walk].z;

#pragma omp parallel for schedule(dynamic)
        for(int iw=0; iw<n_walk; iw++){
            for(int i=0; i<n_epi[iw]; i++){
                int ik = i+ij_disp[iw].x;
                dev_epi[ik].pos.x = epi[iw][i].pos.x;
                dev_epi[ik].pos.y = epi[iw][i].pos.y;
                dev_epi[ik].pos.z = epi[iw][i].pos.z;
                dev_epi[ik].r_search = epi[iw][i].r_search;
                dev_epi[ik].id_walk = iw;
            }
            for(int j=0; j<n_epj[iw]; j++){
                int jk = j+ij_disp[iw].y;
                dev_id_epj[jk] = id_epj[iw][j];
            }
            for(int j=0; j<n_spj[iw]; j++){
                int jk = j+ij_disp[iw].z;
                dev_id_spj[jk] = id_spj[iw][j];
            }
        }
        for(int i=ni_tot; i<ni_tot_reg; i++){
            dev_epi[i].id_walk = n_walk;
        }

#ifdef GPU_PROFILE
        gpu_profile.copy.end();
        gpu_counter.n_walk+= n_walk;
        gpu_counter.n_epi += ni_tot;
        gpu_counter.n_epj += nej_tot;
        gpu_counter.n_spj += nsj_tot;
        gpu_counter.n_call+= 1;
        cudaEventRecord(cu_event_disp);
#endif
        ij_disp.htod(n_walk + 2);
        dev_epi.htod(ni_tot_reg);

#ifdef GPU_PROFILE
        cudaEventRecord(cu_event_htod);
#endif

        int nblocks  = ni_tot_reg / N_THREAD_GPU;
        int nthreads = N_THREAD_GPU;
        force_kernel_ep_ep <<<nblocks, nthreads>>> (ij_disp, dev_epi, dev_epj, dev_id_epj, dev_force, eps2, rcut2, G);
        force_kernel_ep_sp <<<nblocks, nthreads>>> (ij_disp, dev_epi, dev_spj, dev_id_spj, dev_force, eps2, G);

#ifdef GPU_PROFILE
        cudaEventRecord(cu_event_calc);
#endif
        return 0;
    }
}

#else
PS::S32 DispatchKernelWithSP(const PS::S32  tag,
                             const PS::S32  n_walk,
                             const EPISoft *epi[],
                             const PS::S32  n_epi[],
                             const EPJSoft  *epj[],
                             const PS::S32  n_epj[],
                             const SPJSoft *spj[],
                             const PS::S32  n_spj[]){
    assert(n_walk <= N_WALK_LIMIT);
    if(init_call){
		dev_epi  .allocate(NI_LIMIT);
		dev_epj  .allocate(NJ_LIMIT);
        dev_spj  .allocate(NJ_LIMIT);
		dev_force.allocate(NI_LIMIT);
		ij_disp  .allocate(N_WALK_LIMIT+2);
        cudaEventCreate(&cu_event_disp);
        cudaEventCreate(&cu_event_htod);
        cudaEventCreate(&cu_event_calc);
        cudaEventCreate(&cu_event_retr);
        cudaEventCreate(&cu_event_dtoh);
		init_call = false;
    }
#ifdef GPU_PROFILE
    gpu_profile.copy.start();
#endif
    const float eps2 = EPISoft::eps * EPISoft::eps;
    const PS::F64 rcut2 = EPISoft::r_out*EPISoft::r_out;
    const PS::F64 G = ForceSoft::grav_const;
    ij_disp[0].x = 0;
    ij_disp[0].y = 0;
    ij_disp[0].z = 0;
    for(int k=0; k<n_walk; k++){
        ij_disp[k+1].x = ij_disp[k].x + n_epi[k];
        ij_disp[k+1].y = ij_disp[k].y + n_epj[k];
        ij_disp[k+1].z = ij_disp[k].z + n_spj[k];
    }
    ij_disp[n_walk+1] = ij_disp[n_walk];

    assert(ij_disp[n_walk].x < NI_LIMIT);
    assert(ij_disp[n_walk].y < NJ_LIMIT);
    assert(ij_disp[n_walk].z < NJ_LIMIT);

    int ni_tot_reg = ij_disp[n_walk].x;
    if(ni_tot_reg % N_THREAD_GPU){
        ni_tot_reg /= N_THREAD_GPU;
        ni_tot_reg++;
        ni_tot_reg *= N_THREAD_GPU;
    }

    int ni_tot = ij_disp[n_walk].x;
    int nej_tot = ij_disp[n_walk].y;
    int nsj_tot = ij_disp[n_walk].z;

#pragma omp parallel for schedule(dynamic)
    for(int iw=0; iw<n_walk; iw++){
        for(int i=0; i<n_epi[iw]; i++){
            int ik = i+ij_disp[iw].x;
            dev_epi[ik].pos.x = epi[iw][i].pos.x;
            dev_epi[ik].pos.y = epi[iw][i].pos.y;
            dev_epi[ik].pos.z = epi[iw][i].pos.z;
            dev_epi[ik].r_search = epi[iw][i].r_search;
            dev_epi[ik].id_walk = iw;
        }
        for(int j=0; j<n_epj[iw]; j++){
            int jk = j+ij_disp[iw].y;
            dev_epj[jk].pos.x  = epj[iw][j].pos.x;
            dev_epj[jk].pos.y  = epj[iw][j].pos.y;
            dev_epj[jk].pos.z  = epj[iw][j].pos.z;
            dev_epj[jk].m      = epj[iw][j].mass;
            dev_epj[jk].r_search = epj[iw][j].r_search;
        }
        for(int j=0; j<n_spj[iw]; j++){
            int jk = j+ij_disp[iw].z;
            dev_spj[jk].pos.x  = spj[iw][j].pos.x;
            dev_spj[jk].pos.y  = spj[iw][j].pos.y;
            dev_spj[jk].pos.z  = spj[iw][j].pos.z;
            dev_spj[jk].m      = spj[iw][j].getCharge();
#ifdef USE_QUAD
            dev_spj[jk].qxx    = spj[iw][j].quad.xx;
            dev_spj[jk].qyy    = spj[iw][j].quad.yy;
            dev_spj[jk].qzz    = spj[iw][j].quad.zz;
            dev_spj[jk].qxy    = spj[iw][j].quad.xy;
            dev_spj[jk].qxz    = spj[iw][j].quad.xz;
            dev_spj[jk].qyz    = spj[iw][j].quad.yz;
#endif
        }
    }
    for(int i=ni_tot; i<ni_tot_reg; i++){
        dev_epi[i].id_walk = n_walk;
    }

#ifdef GPU_PROFILE
    gpu_profile.copy.end();
    gpu_counter.n_walk+= n_walk;
    gpu_counter.n_epi += ni_tot;
    gpu_counter.n_epj += nej_tot;
    gpu_counter.n_spj += nsj_tot;
    gpu_counter.n_call+= 1;
    cudaEventRecord(cu_event_disp);
#endif
    ij_disp.htod(n_walk + 2);
    dev_epi.htod(ni_tot_reg);
    dev_epj.htod(nej_tot);
    dev_spj.htod(nsj_tot);

#ifdef GPU_PROFILE
    cudaEventRecord(cu_event_htod);
#endif
    int nblocks  = ni_tot_reg / N_THREAD_GPU;
    int nthreads = N_THREAD_GPU;
    force_kernel_ep_ep <<<nblocks, nthreads>>> (ij_disp, dev_epi, dev_epj, dev_force, eps2, rcut2, G);
    force_kernel_ep_sp <<<nblocks, nthreads>>> (ij_disp, dev_epi, dev_spj, dev_force, eps2, G);

#ifdef GPU_PROFILE
    cudaEventRecord(cu_event_calc);
#endif
    return 0;
}

#endif

PS::S32 RetrieveKernel(const PS::S32 tag,
                       const PS::S32 n_walk,
                       const PS::S32 ni[],
                       ForceSoft    *force[]) {

#ifdef GPU_PROFILE
    cudaEventRecord(cu_event_retr);
#endif
    int ni_tot = 0;
    for(int k=0; k<n_walk; k++){
        ni_tot += ni[k];
    }
    dev_force.dtoh(ni_tot);

#ifdef GPU_PROFILE
    cudaEventRecord(cu_event_dtoh);
    cudaEventSynchronize(cu_event_dtoh);
    float send_time=0, calc_time=0, recv_time=0;
    cudaEventElapsedTime(&send_time, cu_event_disp, cu_event_htod);
    cudaEventElapsedTime(&calc_time, cu_event_htod, cu_event_calc);
    cudaEventElapsedTime(&recv_time, cu_event_retr, cu_event_dtoh);
    gpu_profile.send.time += 0.001f*send_time;
    gpu_profile.calc.time += 0.001f*calc_time;
    gpu_profile.recv.time += 0.001f*recv_time;
    gpu_profile.copy.start();
#endif

    int n_cnt = 0;
    for(int iw=0; iw<n_walk; iw++){
        for(int i=0; i<ni[iw]; i++){
            force[iw][i].acc.x = dev_force[n_cnt].accp.x;
            force[iw][i].acc.y = dev_force[n_cnt].accp.y;
            force[iw][i].acc.z = dev_force[n_cnt].accp.z;
            force[iw][i].pot   = dev_force[n_cnt].accp.w;
            force[iw][i].n_ngb = dev_force[n_cnt].nnb;
            n_cnt++;
        }
    }
#ifdef GPU_PROFILE
    gpu_profile.copy.end();
#endif
    return 0;
}

