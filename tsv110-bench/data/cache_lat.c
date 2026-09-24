/* tsv110: cache/DRAM latency via random pointer chase */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

static double now(void){
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC,&ts);
    return ts.tv_sec+1e-9*ts.tv_nsec;
}

int main(void){
    size_t sizes[] = {16*1024UL, 256*1024UL, 8*1024UL*1024UL, 512UL*1024UL*1024UL};
    const char* names[] = {"L1 16KB","L2 256KB","L3 8MB","DRAM 512MB"};
    long steps[] = {200000000L, 200000000L, 50000000L, 8000000L};
    for(int k=0;k<4;k++){
        size_t S = sizes[k];
        size_t n = S / sizeof(uint32_t);
        uint32_t *p = malloc(S);
        uint32_t *perm = malloc(n*sizeof(uint32_t));
        if(!p||!perm){ fprintf(stderr,"alloc fail\n"); return 1; }
        for(size_t i=0;i<n;i++) perm[i]=(uint32_t)i;
        for(size_t i=n-1;i>0;i--){ size_t j=(size_t)rand()% (i+1); uint32_t t=perm[i];perm[i]=perm[j];perm[j]=t; }
        for(size_t i=0;i<n-1;i++) p[perm[i]]=perm[i+1];
        p[perm[n-1]]=perm[0];
        uint32_t idx=0; uint64_t sum=0;
        double dt=1e9;
        for(int rep=0;rep<3;rep++){
            double t0=now();
            for(long s=0;s<steps[k];s++){ idx=p[idx]; sum+=idx; }
            double t=now()-t0;
            if(t<dt) dt=t;
        }
        printf("%-12s n=%9zu  %8.2f ns/access  (checksum %llu)\n", names[k], n, dt/steps[k]*1e9, (unsigned long long)sum);
        free(p); free(perm);
    }
    return 0;
}
