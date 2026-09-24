/* tsv110: memory bandwidth (single & multi thread) */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>

static double now(void){
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC,&ts);
    return ts.tv_sec+1e-9*ts.tv_nsec;
}
static volatile double sink;

static double bench_triad(double*a,double*b,double*c,size_t n,double s,long reps,int nth,int do_read){
    double t0=now();
    for(long r=0;r<reps;r++){
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for(size_t i=0;i<n;i++){
            if(do_read) c[i]=a[i]+s*b[i];
            else c[i]=a[i];
        }
    }
    double dt=now()-t0;
    sink=c[n-1];
    return dt;
}

int main(int argc, char**argv){
    double mb[] = {16e3, 512e3, 4e6, 32e6, 256e6};
    const char* names[] = {"L1 16KB","L2 512KB","L3 4MB","L3 32MB","DRAM 256MB"};
    int ns = sizeof(mb)/sizeof(mb[0]);
    int nth = omp_get_max_threads();
    printf("threads=%d\n", nth);
    printf("%-12s %10s %10s %10s %10s\n","size","read GB/s","triad GB/s","dt(ms)","reps");
    for(int k=0;k<ns;k++){
        size_t S=(size_t)mb[k]; size_t n=S/sizeof(double);
        double *a,*b,*c;
        if(posix_memalign((void**)&a,64,S)||posix_memalign((void**)&b,64,S)||posix_memalign((void**)&c,64,S)) return 1;
        for(size_t i=0;i<n;i++){ a[i]=1.0+1e-9*i; b[i]=2.0; c[i]=0.0; }
        long reps = (long)(2e9/S); if(reps<3) reps=3; if(reps>200000) reps=200000;
        double dtR=1e9, dtT=1e9;
        for(int t=0;t<3;t++){
            double dt = bench_triad(a,b,c,n,0.5,reps,nth,0);
            if(dt<dtR) dtR=dt;
            dt = bench_triad(a,b,c,n,0.5,reps,nth,1);
            if(dt<dtT) dtT=dt;
        }
        double bwR = (double)n*8*reps/dtR/1e9;
        double bwT = (double)n*8*3*reps/dtT/1e9;
        printf("%-12s %10.1f %10.1f %10.2f %10ld\n", names[k], bwR, bwT, dtT*1e3, reps);
        free(a);free(b);free(c);
    }
    return 0;
}
