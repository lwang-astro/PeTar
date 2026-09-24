/* tsv110: accuracy of vrsqrte/vrecpe + Newton iterations */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <arm_neon.h>

static double relerr(double got, double ref){ return fabs(got-ref)/fabs(ref); }

int main(void){
    const int N = 2000000;
    double e0=0,e1=0,e2=0,ec=0, e0d=0,e1d=0,e2d=0;
    double m0=0,m1=0,m2=0,mc=0, m0d=0,m1d=0,m2d=0;
    unsigned seed = 12345;
    for(int i=0;i<N;i++){
        float x = (float)pow(10.0,-3.0+6.0*((double)rand_r(&seed)/RAND_MAX));
        double ref = 1.0/sqrt((double)x);

        float32x2_t xv = vdup_n_f32(x);
        float32x2_t r0 = vrsqrte_f32(xv);
        double v0 = vget_lane_f32(r0,0);
        e0 = fmax(e0, relerr(v0,ref)); m0 += relerr(v0,ref);

        float32x2_t h = vmul_f32(vmul_f32(xv,r0),r0);              /* x*r^2 */
        float32x2_t t = vsub_f32(vdup_n_f32(1.5f), vmul_f32(vdup_n_f32(0.5f),h));
        float32x2_t r1 = vmul_f32(r0,t);
        double v1 = vget_lane_f32(r1,0);
        e1 = fmax(e1, relerr(v1,ref)); m1 += relerr(v1,ref);

        h = vmul_f32(vmul_f32(xv,r1),r1);
        t = vsub_f32(vdup_n_f32(1.5f), vmul_f32(vdup_n_f32(0.5f),h));
        float32x2_t r2 = vmul_f32(r1,t);
        double v2 = vget_lane_f32(r2,0);
        e2 = fmax(e2, relerr(v2,ref)); m2 += relerr(v2,ref);

        /* Fugaku-style cubic correction: h=1-x*r^2; r=r+r*h*(0.5+0.375h) */
        float32x2_t hh = vsub_f32(vdup_n_f32(1.0f), vmul_f32(vmul_f32(xv,r0),r0));
        float32x2_t poly = vmul_f32(vadd_f32(vmul_f32(vdup_n_f32(0.375f),hh), vdup_n_f32(0.5f)), hh);
        float32x2_t rc = vadd_f32(r0, vmul_f32(r0,poly));
        double vc = vget_lane_f32(rc,0);
        ec = fmax(ec, relerr(vc,ref)); mc += relerr(vc,ref);

        /* f64 rsqrt estimate */
        double xd = pow(10.0,-3.0+6.0*((double)rand_r(&seed)/RAND_MAX));
        double refd = 1.0/sqrt(xd);
        float64x1_t xv64 = vdup_n_f64(xd);
        float64x1_t q0 = vrsqrte_f64(xv64);
        double w0 = vget_lane_f64(q0,0);
        e0d = fmax(e0d, relerr(w0,refd)); m0d += relerr(w0,refd);
        float64x1_t h64 = vmul_f64(vmul_f64(xv64,q0),q0);
        float64x1_t t64 = vsub_f64(vdup_n_f64(1.5), vmul_f64(vdup_n_f64(0.5),h64));
        float64x1_t q1 = vmul_f64(q0,t64);
        double w1 = vget_lane_f64(q1,0);
        e1d = fmax(e1d, relerr(w1,refd)); m1d += relerr(w1,refd);
        h64 = vmul_f64(vmul_f64(xv64,q1),q1);
        t64 = vsub_f64(vdup_n_f64(1.5), vmul_f64(vdup_n_f64(0.5),h64));
        float64x1_t q2 = vmul_f64(q1,t64);
        double w2 = vget_lane_f64(q2,0);
        e2d = fmax(e2d, relerr(w2,refd)); m2d += relerr(w2,refd);
    }
    printf("=== rsqrt f32 (2-lane NEON est) vs exact, N=%d ===\n", N);
    printf("vrsqrte          : max_rel=%.3e mean_rel=%.3e\n", e0, m0/N);
    printf("+1 NR            : max_rel=%.3e mean_rel=%.3e\n", e1, m1/N);
    printf("+2 NR            : max_rel=%.3e mean_rel=%.3e\n", e2, m2/N);
    printf("+1 cubic(Fugaku) : max_rel=%.3e mean_rel=%.3e\n", ec, mc/N);
    printf("=== rsqrt f64 (1-lane NEON est) ===\n");
    printf("vrsqrte          : max_rel=%.3e mean_rel=%.3e\n", e0d, m0d/N);
    printf("+1 NR            : max_rel=%.3e mean_rel=%.3e\n", e1d, m1d/N);
    printf("+2 NR            : max_rel=%.3e mean_rel=%.3e\n", e2d, m2d/N);
    return 0;
}
