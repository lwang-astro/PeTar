/* tsv110: independent clock estimate via dependent-ALU chains (ns/op) */
#define _GNU_SOURCE
#include <stdio.h>
#include <time.h>
#include <sched.h>

static double now(void){
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC,&ts);
    return ts.tv_sec+1e-9*ts.tv_nsec;
}

int main(void){
    cpu_set_t set; CPU_ZERO(&set); CPU_SET(0,&set);
    sched_setaffinity(0,sizeof(set),&set);
    const long N=500000000L;
    double t0,t1;

    unsigned long x=0;
    t0=now();
    for(long i=0;i<N;i++) __asm__ __volatile__("add %0,%0,#1":"+r"(x));
    t1=now();
    double ns_add=(t1-t0)/N*1e9;
    printf("dependent int add : %.4f ns/op  -> clock=%.3f GHz (if latency=1 cyc)\n", ns_add, 1.0/ns_add);

    double d=1.0, one=1.0, small=1e-18;
    t0=now();
    for(long i=0;i<N;i++) __asm__ __volatile__("fadd %d0,%d0,%d1":"+w"(d):"w"(small));
    t1=now();
    double ns_fadd=(t1-t0)/N*1e9;
    printf("dependent fadd    : %.4f ns/op  (%.2f cyc @2.6GHz)\n", ns_fadd, ns_fadd*2.6);

    double a=1.0, b=1e-9, c=1e-9;
    t0=now();
    for(long i=0;i<N;i++) __asm__ __volatile__("fmadd %d0,%d1,%d2,%d0":"+w"(a):"w"(b),"w"(c));
    t1=now();
    double ns_fma=(t1-t0)/N*1e9;
    printf("dependent scalar fmadd: %.4f ns/op  (%.2f cyc @2.6GHz)\n", ns_fma, ns_fma*2.6);
    (void)x;(void)d;(void)one;(void)a;
    return 0;
}
