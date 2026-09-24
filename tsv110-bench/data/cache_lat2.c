/* tsv110: clean cache/DRAM latency: pure 8-byte pointer chase, pinned, optional THP */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include <sys/mman.h>
#include <sched.h>

static double now(void){
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC,&ts);
    return ts.tv_sec+1e-9*ts.tv_nsec;
}

static uint64_t rng_state=88172645463325252ull;
static inline uint64_t xorshift64(void){
    uint64_t x=rng_state;
    x^=x<<13; x^=x>>7; x^=x<<17;
    return rng_state=x;
}

static void print_thp(void){
    FILE* f=fopen("/proc/self/smaps_rollup","r");
    char line[256]; if(!f) return;
    while(fgets(line,sizeof line,f)) if(strstr(line,"AnonHugePages")){ printf("  %s",line); break; }
    fclose(f);
}

int main(int argc, char** argv){
    cpu_set_t set; CPU_ZERO(&set); CPU_SET(0,&set);
    if(sched_setaffinity(0,sizeof(set),&set)!=0) perror("setaffinity");
    int use_thp = argc>1 ? atoi(argv[1]) : 0;

    size_t sizes[] = {16UL*1024, 256UL*1024, 8UL*1024*1024, 512UL*1024*1024};
    const char* names[] = {"L1 16KB","L2 256KB","L3 8MB","DRAM 512MB"};
    long steps[] = {100000000L, 100000000L, 30000000L, 5000000L};

    for(int k=0;k<4;k++){
        size_t S = sizes[k];
        size_t n = S/8;                      /* number of 8-byte nodes */
        uint64_t* arr = mmap(NULL, S, PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);
        uint32_t* perm = malloc(n*sizeof(uint32_t));
        if(arr==MAP_FAILED || !perm){ fprintf(stderr,"alloc fail\n"); return 1; }
        if(use_thp) madvise(arr, S, MADV_HUGEPAGE);
        for(size_t i=0;i<n;i++) perm[i]=(uint32_t)i;
        for(size_t i=n-1;i>0;i--){
            size_t j=(size_t)(xorshift64()%(i+1));
            uint32_t t=perm[i]; perm[i]=perm[j]; perm[j]=t;
        }
        for(size_t i=0;i+1<n;i++) arr[perm[i]] = (uint64_t)(uintptr_t)&arr[perm[i+1]];
        arr[perm[n-1]] = (uint64_t)(uintptr_t)&arr[perm[0]];

        volatile uint64_t* q = (volatile uint64_t*)&arr[perm[0]];
        for(long s=0;s<1000000;s++) q=(volatile uint64_t*)*q;   /* warmup */

        double best=1e300;
        for(int r=0;r<3;r++){
            double t0=now();
            for(long s=0;s<steps[k];s++) q=(volatile uint64_t*)*q;
            double dt=now()-t0;
            if(dt<best) best=dt;
        }
        printf("%-12s n=%10zu  %8.2f ns/access  thp=%d\n", names[k], n, best/steps[k]*1e9, use_thp);
        if(use_thp) print_thp();
        fflush(stdout);
        munmap(arr,S); free(perm);
    }
    return 0;
}
