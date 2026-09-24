/* rsinv_test.c - accuracy and throughput of the NEON rsqrt used by the quadrupole kernel
 *
 *   cubic : r = rsqrt4(x)                      (3rd order corrected estimate)
 *   newton: cubic + one half Newton step        (Fugaku quadrupole variant)
 *
 * usage: ./rsinv_test <cubic|newton> [n] [nthr]
 * output: CSV lines STAT/HIST/THR
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <arm_neon.h>

static double now(void){
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC,&ts);
    return ts.tv_sec+1e-9*ts.tv_nsec;
}

static inline float32x4_t rsqrt4_base(float32x4_t x){
    float32x4_t r = vrsqrteq_f32(x);
    float32x4_t h = vmulq_f32(x, r);
    h = vfmsq_f32(vdupq_n_f32(1.0f), h, r);
    float32x4_t p = vfmaq_n_f32(vdupq_n_f32(0.5f), h, 0.375f);
    p = vmulq_f32(p, h);
    return vfmaq_f32(r, r, p);
}
static inline float32x4_t rsqrt4_newton(float32x4_t x){
    float32x4_t r = rsqrt4_base(x);
    float32x4_t h = vmulq_f32(r, r);
    h = vfmsq_f32(vdupq_n_f32(3.0f), x, h);
    return vmulq_f32(r, vmulq_f32(h, vdupq_n_f32(0.5f)));
}

static int cmp_double(const void*a, const void*b){
    double x=*(const double*)a, y=*(const double*)b;
    return (x>y)-(x<y);
}

int main(int argc, char** argv){
    const char* method = argc>1 ? argv[1] : "cubic";
    long n = argc>2 ? atol(argv[2]) : 2000000L;
    int newton = strcmp(method,"newton")==0;

    double* err = malloc(n*sizeof(double));
    if(!err) return 1;
    unsigned seed = 20260924u;
    double sum=0, sumsq=0, emax=0;
    const int NB=300;
    long hist[NB+2]; for(int i=0;i<NB+2;i++) hist[i]=0;
    for(long i=0;i<n;i++){
        /* log-uniform x in [1e-4, 1e4] (covers the r^2 range of the tree) */
        double e10 = -4.0 + 8.0*((double)rand_r(&seed)/RAND_MAX);
        float x = (float)pow(10.0, e10);
        float32x4_t xv = vdupq_n_f32(x);
        float32x4_t r = newton ? rsqrt4_newton(xv) : rsqrt4_base(xv);
        float r0 = vgetq_lane_f32(r,0);
        double ref = 1.0/sqrt((double)x);
        double e = fabs((double)r0-ref)/ref;
        err[i]=e; sum+=e; sumsq+=e*e; if(e>emax) emax=e;
        int b;
        if(e<=0) b=0;
        else {
            double lg = log10(e);
            b = (int)floor((lg+10.0)/10.0*NB);
            if(b<0) b=0; if(b>NB+1) b=NB+1;
        }
        hist[b]++;
    }
    qsort(err,n,sizeof(double),cmp_double);
    printf("STAT,%s,n=%ld,max=%.6e,mean=%.6e,rms=%.6e,p50=%.6e,p90=%.6e,p99=%.6e\n",
           method,(long)n,emax,sum/n,sqrt(sumsq/n),err[n/2],err[(long)(n*0.9)],err[(long)(n*0.99)]);
    for(int b=0;b<NB+2;b++){
        if(hist[b]==0) continue;
        double lg = b==0 ? -10.0 : (b==NB+1 ? 0.0 : -10.0+10.0*(b+0.5)/NB);
        printf("HIST,%s,%.5f,%ld\n", method, lg, hist[b]);
    }
    free(err);

    /* throughput: independent ops */
    const long N=500000000L;
    float32x4_t x[8]; float32x4_t y[8];
    for(int u=0;u<8;u++){ x[u]=vdupq_n_f32(1.0f+0.01f*u); y[u]=vdupq_n_f32(0.0f); }
    volatile float sink;
    double t0=now();
    for(long i=0;i<N;i++){
        int u=(int)(i&7);
        x[u]=vreinterpretq_f32_u32(vaddq_u32(vreinterpretq_u32_f32(x[u]),vdupq_n_u32(0x00800000u)));
        y[u] = newton ? rsqrt4_newton(x[u]) : rsqrt4_base(x[u]);
    }
    double t1=now();
    float s=0; for(int u=0;u<8;u++) s+=vaddvq_f32(y[u]);
    sink=s; (void)sink;
    printf("THR,%s,nops=%ld,ns_per_op=%.5f\n", method, N, (t1-t0)/N*1e9);
    return 0;
}
