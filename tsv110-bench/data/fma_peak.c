/* tsv110 microbenchmark: NEON FMA peak / latency, div, sqrt, rsqrt-est, recpe */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <signal.h>
#include <setjmp.h>
#include <arm_neon.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static double now(void){
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9*ts.tv_nsec;
}

static volatile float  sinkf;
static volatile double sinkd;

/* ---- PMCCNTR_EL0 probe (user-space cycles counter) ---- */
static sigjmp_buf jb;
static void sigill_handler(int sig){ (void)sig; siglongjmp(jb,1); }
static inline uint64_t rd_pmccntr(void){
    uint64_t v; asm volatile("mrs %0, pmccntr_el0" : "=r"(v)); return v;
}
static inline uint64_t rd_cntvct(void){
    uint64_t v; asm volatile("mrs %0, cntvct_el0" : "=r"(v)); return v;
}
static int pmccntr_ok(void){
    struct sigaction sa, old;
    memset(&sa,0,sizeof sa); sa.sa_handler = sigill_handler; sigaction(SIGILL,&sa,&old);
    int ok = 0;
    if(sigsetjmp(jb,1)==0){ volatile uint64_t v = rd_pmccntr(); (void)v; ok = 1; }
    sigaction(SIGILL,&old,NULL);
    return ok;
}

enum MMODE { M_FMA32, M_LAT32, M_FMA64, M_LAT64,
             M_RSQRTE32, M_RECPE32, M_SQRT32, M_DIV32,
             M_RSQRTE64, M_RECPE64, M_SQRT64, M_DIV64 };

int main(int argc, char**argv){
    const char* mode = argc>1 ? argv[1] : "fma32";
    long N = argc>2 ? atol(argv[2]) : 100000000L;
    int nth = 1;
#ifdef _OPENMP
    nth = omp_get_max_threads();
#endif

    if(strcmp(mode,"info")==0){
        printf("threads %d pmccntr_ok %d\n", nth, pmccntr_ok());
        uint64_t c0 = rd_cntvct(); double t0 = now();
        volatile double x = 0; for(long i=0;i<200000000L;i++) x += i*1e-9;
        double t1 = now(); uint64_t c1 = rd_cntvct();
        printf("cntfrq=%.3f MHz (sink %.3f)\n", (c1-c0)/(t1-t0)/1e6, (double)x);
        return 0;
    }

    int m = -1;
    if(!strcmp(mode,"fma32")) m=M_FMA32; else if(!strcmp(mode,"lat32")) m=M_LAT32;
    else if(!strcmp(mode,"fma64")) m=M_FMA64; else if(!strcmp(mode,"lat64")) m=M_LAT64;
    else if(!strcmp(mode,"rsqrte32")) m=M_RSQRTE32; else if(!strcmp(mode,"recpe32")) m=M_RECPE32;
    else if(!strcmp(mode,"sqrt32")) m=M_SQRT32; else if(!strcmp(mode,"div32")) m=M_DIV32;
    else if(!strcmp(mode,"rsqrte64")) m=M_RSQRTE64; else if(!strcmp(mode,"recpe64")) m=M_RECPE64;
    else if(!strcmp(mode,"sqrt64")) m=M_SQRT64; else if(!strcmp(mode,"div64")) m=M_DIV64;
    if(m<0){ fprintf(stderr,"unknown mode %s\n",mode); return 1; }

    double t0=now(), t1=now(), ops=0, gflop=0;

    if(m==M_FMA32 || m==M_FMA64){
        t0 = now();
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            float sf=0; double sd=0;
            if(m==M_FMA32){
                float32x4_t a[10];
                for(int u=0;u<10;u++) a[u]=vdupq_n_f32(1.0f+u*1e-3f);
                const float32x4_t b=vdupq_n_f32(0.9999f), c=vdupq_n_f32(1e-7f);
#ifdef _OPENMP
#pragma omp for
#endif
                for(long i=0;i<N;i++)
                    for(int u=0;u<10;u++) a[u]=vfmaq_f32(a[u],b,c);
                for(int u=0;u<10;u++) sf+=vaddvq_f32(a[u]);
            }else{
                float64x2_t a[10];
                for(int u=0;u<10;u++) a[u]=vdupq_n_f64(1.0+u*1e-3);
                const float64x2_t b=vdupq_n_f64(0.9999), c=vdupq_n_f64(1e-7);
#ifdef _OPENMP
#pragma omp for
#endif
                for(long i=0;i<N;i++)
                    for(int u=0;u<10;u++) a[u]=vfmaq_f64(a[u],b,c);
                for(int u=0;u<10;u++) sd+=vaddvq_f64(a[u]);
            }
            sinkf=sf; sinkd=sd;
        }
        t1 = now();
        ops = (double)N*10.0;
        gflop = (m==M_FMA32 ? 8.0 : 4.0) * ops;
    } else if(m==M_LAT32 || m==M_LAT64){
        t0 = now();
        if(m==M_LAT32){
            float32x4_t a=vdupq_n_f32(1.0f);
            const float32x4_t b=vdupq_n_f32(0.9999f), c=vdupq_n_f32(1e-7f);
            for(long i=0;i<N;i++) a=vfmaq_f32(a,b,c);
            sinkf=vaddvq_f32(a);
        }else{
            float64x2_t a=vdupq_n_f64(1.0);
            const float64x2_t b=vdupq_n_f64(0.9999), c=vdupq_n_f64(1e-7);
            for(long i=0;i<N;i++) a=vfmaq_f64(a,b,c);
            sinkd=vaddvq_f64(a);
        }
        t1 = now();
        ops = (double)N;
    } else {
        /* independent-op throughput: x is perturbed by a cheap integer op,
           op result goes to a separate accumulator so latency doesn't limit */
        t0 = now();
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            float sf=0; double sd=0;
            float32x4_t xf[8]; float64x2_t xd[8];
            float32x4_t yf[8]; float64x2_t yd[8];
            for(int u=0;u<8;u++){
                xf[u]=vdupq_n_f32(1.0f+0.01f*u); xd[u]=vdupq_n_f64(1.0+0.01*u);
                yf[u]=vdupq_n_f32(0.0f); yd[u]=vdupq_n_f64(0.0);
            }
#ifdef _OPENMP
#pragma omp for
#endif
            for(long i=0;i<N;i++){
                int u = (int)(i&7);
                if(m>=M_RSQRTE32 && m<=M_DIV32){
                    xf[u]=vreinterpretq_f32_u32(vaddq_u32(vreinterpretq_u32_f32(xf[u]),vdupq_n_u32(0x00800000u)));
                    switch(m){
                        case M_RSQRTE32: yf[u]=vrsqrteq_f32(xf[u]); break;
                        case M_RECPE32:  yf[u]=vrecpeq_f32(xf[u]); break;
                        case M_SQRT32:   yf[u]=vsqrtq_f32(xf[u]); break;
                        default:         yf[u]=vdivq_f32(xf[u],vdupq_n_f32(1.3f)); break;
                    }
                }else{
                    xd[u]=vreinterpretq_f64_u64(vaddq_u64(vreinterpretq_u64_f64(xd[u]),vdupq_n_u64(0x0008000000000000ull)));
                    switch(m){
                        case M_RSQRTE64: yd[u]=vrsqrteq_f64(xd[u]); break;
                        case M_RECPE64:  yd[u]=vrecpeq_f64(xd[u]); break;
                        case M_SQRT64:   yd[u]=vsqrtq_f64(xd[u]); break;
                        default:         yd[u]=vdivq_f64(xd[u],vdupq_n_f64(1.3)); break;
                    }
                }
            }
            for(int u=0;u<8;u++){ sf+=vaddvq_f32(yf[u]); sd+=vaddvq_f64(yd[u]); }
            sinkf=sf; sinkd=sd;
        }
        t1 = now();
        ops = (double)N;
    }

    double dt = t1-t0;
    if(gflop>0)
        printf("%-9s threads=%2d N=%ld time=%.4f s  %.3f ns/vector-instr  -> %.2f GFLOP/s (%.2f /thread)\n",
               mode, nth, N, dt, dt/ops*1e9, gflop/dt/1e9, gflop/dt/1e9/nth);
    else
        printf("%-9s threads=%2d N=%ld time=%.4f s  %.3f ns/vector-instr  -> %.3f Gvector-instr/s\n",
               mode, nth, N, dt, dt/ops*1e9, ops/dt/1e9);
    return 0;
}
