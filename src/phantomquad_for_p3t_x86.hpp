// gcc -O2 -march=core-avx2
#include <cassert>
#include<immintrin.h>
class PhantomGrapeQuad{
public:
    enum{
	NIMAX = 32768,
	NJMAX = 131072,
    };
    
private:
#if 1
    float xibuf  [NIMAX/8]  [3][8];   // x, y, z
    double xibufd  [NIMAX/8]  [3][8];   // x, y, z
    float accpbuf[NIMAX/8]  [5][8];   // ax, ay, az, pot, nngb
    double accpbufd[NIMAX/8]  [5][8];   // ax, ay, az, pot, nngb
    float epjbuf [NJMAX]    [4];      // x, y, z, m
    double epjbufd [NJMAX]    [4];      // x, y, z, m
    float spjbuf [NJMAX]    [3][4];   // x, y, z, m, | xx, yy, zz, pad, | xy, yz, zx, tr
    double spjbufd [NJMAX]    [3][4];   // x, y, z, m, | xx, yy, zz, pad, | xy, yz, zx, tr
#else
    float *** xibuf;
    float *** accpbuf;
    float ** epjbuf;
    float *** spjbuf;
#endif
    double eps2;
    static double get_a_NaN(){
	union{ long   l; double d; } m;
	m.l = -1;
	return m.d;
    }
    double r_crit2;
    double r_out;
    double r_in;
    double denominator; // for cut off
public:

    PhantomGrapeQuad() : eps2(get_a_NaN()) {} // default NaN

    void set_cutoff(const double _r_out, const double _r_in){
        r_out = _r_out;
        r_in = _r_in;
        denominator = 1.0 / (r_out - r_in);
    }
    
    void set_eps2(const double _eps2){
        this->eps2 = _eps2;
    }
    void set_r_crit2(const double _r_crit2){
        this->r_crit2 = _r_crit2;
        //this->r_crit2 = _r_crit2 * 1.01;
    }
    void set_epj_one(const int addr, const double x, const double y, const double z, const double m) {
        epjbuf[addr][0] = x;
        epjbuf[addr][1] = y;
        epjbuf[addr][2] = z;
        epjbuf[addr][3] = m;
    }

    void set_epj_one_d(const int addr, const double x, const double y, const double z, const double m){
        epjbufd[addr][0] = x;
        epjbufd[addr][1] = y;
        epjbufd[addr][2] = z;
        epjbufd[addr][3] = m;
    }

    // please specialize the following
    template <typename EPJ_t>
    void set_epj(const int nj, const EPJ_t epj[]);
    
    void set_spj_one(
		     const int addr, 
		     const double x,   const double y,   const double z,   const double m,
		     const double qxx, const double qyy, const double qzz,
		     const double qxy, const double qyz, const double qzx)
    {
	const double tr = qxx + qyy + qzz;
	spjbuf[addr][0][0] = x;
	spjbuf[addr][0][1] = y;
	spjbuf[addr][0][2] = z;
	spjbuf[addr][0][3] = m;

	spjbuf[addr][1][0] = 3.0 * qxx - tr;
	spjbuf[addr][1][1] = 3.0 * qyy - tr;
	spjbuf[addr][1][2] = 3.0 * qzz - tr;
	spjbuf[addr][1][3] = m;

	spjbuf[addr][2][0] = 3.0 * qxy;
	spjbuf[addr][2][1] = 3.0 * qyz;
	spjbuf[addr][2][2] = 3.0 * qzx;
	spjbuf[addr][2][3] = -(eps2 * tr);
    }

    void set_spj_one_d(
		       const int addr, 
		       const double x,   const double y,   const double z,   const double m,
		       const double qxx, const double qyy, const double qzz,
		       const double qxy, const double qyz, const double qzx)
    {
	const double tr = qxx + qyy + qzz;
	spjbufd[addr][0][0] = x;
	spjbufd[addr][0][1] = y;
	spjbufd[addr][0][2] = z;
	spjbufd[addr][0][3] = m;

	spjbufd[addr][1][0] = 3.0 * qxx - tr;
	spjbufd[addr][1][1] = 3.0 * qyy - tr;
	spjbufd[addr][1][2] = 3.0 * qzz - tr;
	spjbufd[addr][1][3] = m;

	spjbufd[addr][2][0] = 3.0 * qxy;
	spjbufd[addr][2][1] = 3.0 * qyz;
	spjbufd[addr][2][2] = 3.0 * qzx;
	spjbufd[addr][2][3] = -(eps2 * tr);
    }

    // please specialize the following
    template <typename SPJ_t>
    void set_spj(const int nj, const SPJ_t spj[]);

    void set_xi_one(const int addr, const double x, const double y, const double z){
        const int ah = addr / 8;
        const int al = addr % 8;
        xibuf[ah][0][al] = x;
        xibuf[ah][1][al] = y;
        xibuf[ah][2][al] = z;
    }
    void set_xi_one_d(const int addr, const double x, const double y, const double z){
        const int ah = addr / 8;
        const int al = addr % 8;
        xibufd[ah][0][al] = x;
        xibufd[ah][1][al] = y;
        xibufd[ah][2][al] = z;
    }

    // please specialize the following
    template <typename EPI_t>
    void set_xi(const int ni, const EPI_t spj[]);

    template <typename real_t>
    void get_accp_one(const int addr, real_t &ax, real_t &ay, real_t &az, real_t &pot){
	const int ah = addr / 8;
	const int al = addr % 8;
	ax  = accpbuf[ah][0][al];
	ay  = accpbuf[ah][1][al];
	az  = accpbuf[ah][2][al];
	pot = accpbuf[ah][3][al];
    }
    template <typename real_t>
    void get_accp_one(const int addr, real_t &ax, real_t &ay, real_t &az, real_t &pot, real_t &nngb){
	const int ah = addr / 8;
	const int al = addr % 8;
	ax  = accpbuf[ah][0][al];
	ay  = accpbuf[ah][1][al];
	az  = accpbuf[ah][2][al];
	pot = accpbuf[ah][3][al];
	nngb = accpbuf[ah][4][al];
    }

    template <typename real_t>
    void get_accp_one_d(const int addr, real_t &ax, real_t &ay, real_t &az, real_t &pot, real_t &nngb){
	const int ah = addr / 8;
	const int al = addr % 8;
	ax  = accpbufd[ah][0][al];
	ay  = accpbufd[ah][1][al];
	az  = accpbufd[ah][2][al];
	pot = accpbufd[ah][3][al];
	nngb = accpbufd[ah][4][al];
    }

    template <typename real_t>
    void accum_accp_one(const int addr, real_t &ax, real_t &ay, real_t &az, real_t &pot){
	const int ah = addr / 8;
	const int al = addr % 8;
	ax  += accpbuf[ah][0][al];
	ay  += accpbuf[ah][1][al];
	az  += accpbuf[ah][2][al];
	pot += accpbuf[ah][3][al];
    }

    template <typename real_t>
    void accum_accp_one_d(const int addr, real_t &ax, real_t &ay, real_t &az, real_t &pot){
	const int ah = addr / 8;
	const int al = addr % 8;
	ax  += accpbufd[ah][0][al];
	ay  += accpbufd[ah][1][al];
	az  += accpbufd[ah][2][al];
	pot += accpbufd[ah][3][al];
    }

    template <typename real_t>
    void accum_accp_one(const int addr, real_t &ax, real_t &ay, real_t &az, real_t &pot, real_t &nngb){
	const int ah = addr / 8;
	const int al = addr % 8;
	ax  += accpbuf[ah][0][al];
	ay  += accpbuf[ah][1][al];
	az  += accpbuf[ah][2][al];
	pot += accpbuf[ah][3][al];
	nngb += accpbuf[ah][4][al];
    }

    template <typename real_t>
    void accum_accp_one_d(const int addr, real_t &ax, real_t &ay, real_t &az, real_t &pot, real_t &nngb){
	const int ah = addr / 8;
	const int al = addr % 8;
	ax  += accpbufd[ah][0][al];
	ay  += accpbufd[ah][1][al];
	az  += accpbufd[ah][2][al];
	pot += accpbufd[ah][3][al];
	nngb += accpbufd[ah][4][al];
    }

    void run_epj(const int ni, const int nj){
	if(ni > NIMAX || nj > NJMAX){
	    std::cout<<"ni= "<<ni<<" NIMAX= "<<NIMAX<<" nj= "<<nj<<" NJMAX= "<<NJMAX<<std::endl;
	    for(PS::S32 i=0; i<(ni-1)/8+1; i++){
		for(PS::S32 j=0; j<8; j++){
		    for(PS::S32 k=0; k<3; k++){
			std::cout<<"i,j,k="<<i<<" "<<j<<" "<<k<<std::endl;
			std::cout<<"xibuf[i][k][j]="<<xibuf[i][k][j]<<std::endl;
		    }
		    std::cout<<std::endl;
		}
	    }
	}
	assert(ni <= NIMAX);
	assert(nj <= NJMAX);

	kernel_epj_nounroll(ni, nj);
    }

    void run_epj_d(const int ni, const int nj){
	if(ni > NIMAX || nj > NJMAX){
	    std::cout<<"ni= "<<ni<<" NIMAX= "<<NIMAX<<" nj= "<<nj<<" NJMAX= "<<NJMAX<<std::endl;
	    for(PS::S32 i=0; i<(ni-1)/8+1; i++){
		for(PS::S32 j=0; j<8; j++){
		    for(PS::S32 k=0; k<3; k++){
			std::cout<<"i,j,k="<<i<<" "<<j<<" "<<k<<std::endl;
			std::cout<<"xibuf[i][k][j]="<<xibuf[i][k][j]<<std::endl;
		    }
		    std::cout<<std::endl;
		}
	    }
	}
	assert(ni <= NIMAX);
	assert(nj <= NJMAX);
	kernel_epj_nounroll_64bit(ni, nj);
    }

    //////////
    // include linear cutoff
    void run_epj_for_p3t_with_linear_cutoff(const int ni, const int nj){
        if(ni > NIMAX || nj > NJMAX){
            std::cout<<"ni= "<<ni<<" NIMAX= "<<NIMAX<<" nj= "<<nj<<" NJMAX= "<<NJMAX<<std::endl;
            for(PS::S32 i=0; i<(ni-1)/8+1; i++){
                for(PS::S32 j=0; i<8; j++){
                    for(PS::S32 k=0; k<3; k++){
                        std::cout<<"xibufd[i][k][j]="<<xibufd[i][k][j]<<std::endl;
                    }
                    std::cout<<std::endl;
                }
            }
        }
        assert(ni <= NIMAX);
        assert(nj <= NJMAX);
        kernel_epj_nounroll_for_p3t_with_linear_cutoff(ni, nj);
    }

    void run_epj_for_p3t_with_linear_cutoff_d(const int ni, const int nj){
	if(ni > NIMAX || nj > NJMAX){
	    std::cout<<"ni= "<<ni<<" NIMAX= "<<NIMAX<<" nj= "<<nj<<" NJMAX= "<<NJMAX<<std::endl;
	    for(PS::S32 i=0; i<(ni-1)/8+1; i++){
		for(PS::S32 j=0; i<8; j++){
		    for(PS::S32 k=0; k<3; k++){
			std::cout<<"xibufd[i][k][j]="<<xibufd[i][k][j]<<std::endl;
		    }
		    std::cout<<std::endl;
		}
	    }
	}
	assert(ni <= NIMAX);
	assert(nj <= NJMAX);
	kernel_epj_64bit_nounroll_for_p3t_with_linear_cutoff(ni, nj);
    }
    // include linear cutoff
    //////////

    void run_epj_for_p3t(const int ni, const int nj){
	if(ni > NIMAX || nj > NJMAX){
	    std::cout<<"ni= "<<ni<<" NIMAX= "<<NIMAX<<" nj= "<<nj<<" NJMAX= "<<NJMAX<<std::endl;
	    for(PS::S32 i=0; i<(ni-1)/8+1; i++){
		for(PS::S32 j=0; i<8; j++){
		    for(PS::S32 k=0; k<3; k++){
			std::cout<<"xibufd[i][k][j]="<<xibufd[i][k][j]<<std::endl;
		    }
		    std::cout<<std::endl;
		}
	    }
	}
	assert(ni <= NIMAX);
	assert(nj <= NJMAX);
	kernel_epj_nounroll_for_p3t(ni, nj);
    }

    void run_epj_for_p3t_d(const int ni, const int nj){
	if(ni > NIMAX || nj > NJMAX){
	    std::cout<<"ni= "<<ni<<" NIMAX= "<<NIMAX<<" nj= "<<nj<<" NJMAX= "<<NJMAX<<std::endl;
	    for(PS::S32 i=0; i<(ni-1)/8+1; i++){
		for(PS::S32 j=0; i<8; j++){
		    for(PS::S32 k=0; k<3; k++){
			std::cout<<"xibufd[i][k][j]="<<xibufd[i][k][j]<<std::endl;
		    }
		    std::cout<<std::endl;
		}
	    }
	}
	assert(ni <= NIMAX);
	assert(nj <= NJMAX);
	kernel_epj_64bit_nounroll_for_p3t(ni, nj);
    }

    void run_spj(const int ni, const int nj){
	if(ni > NIMAX || nj > NJMAX){
	    std::cout<<"ni= "<<ni<<" NIMAX= "<<NIMAX<<" nj= "<<nj<<" NJMAX= "<<NJMAX<<std::endl;
	}
	assert(ni <= NIMAX);
	assert(nj <= NJMAX);

	kernel_spj_nounroll(ni, nj);
	// kernel_spj_unroll2(ni, nj);
    }

    void run_spj_d(const int ni, const int nj){
	if(ni > NIMAX || nj > NJMAX){
	    std::cout<<"ni= "<<ni<<" NIMAX= "<<NIMAX<<" nj= "<<nj<<" NJMAX= "<<NJMAX<<std::endl;
	}
	assert(ni <= NIMAX);
	assert(nj <= NJMAX);

	kernel_spj_64bit_nounroll(ni, nj);
	// kernel_spj_unroll2(ni, nj);
    }

private:
	typedef float v4sf __attribute__((vector_size(16)));
	typedef float v8sf __attribute__((vector_size(32)));
	typedef double v4df __attribute__((vector_size(32)));
	
	__attribute__ ((noinline))
	void kernel_epj_nounroll(const int ni, const int nj){

	    const v8sf veps2 = {(float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2};
		for(int i=0; i<ni; i+=8){
		        const v8sf xi = *(v8sf *)(xibuf[i/8][0]);
			const v8sf yi = *(v8sf *)(xibuf[i/8][1]);
			const v8sf zi = *(v8sf *)(xibuf[i/8][2]);

			v8sf ax, ay, az, pot;
			ax = ay = az = pot = (v8sf){0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
			v8sf jbuf = __builtin_ia32_vbroadcastf128_ps256((v4sf *)epjbuf);
			v8sf xj =  __builtin_ia32_shufps256(jbuf, jbuf, 0x00);
			v8sf yj =  __builtin_ia32_shufps256(jbuf, jbuf, 0x55);
			v8sf zj =  __builtin_ia32_shufps256(jbuf, jbuf, 0xaa);
			v8sf mj =  __builtin_ia32_shufps256(jbuf, jbuf, 0xff);
			for(int j=0; j<nj; j++){
				jbuf = __builtin_ia32_vbroadcastf128_ps256((v4sf *)(epjbuf + j+1));

				v8sf dx = xj - xi;
				v8sf dy = yj - yi;
				v8sf dz = zj - zi;

				v8sf r2   = ((veps2 + dx*dx) + dy*dy) + dz*dz;
				v8sf ri1  = __builtin_ia32_rsqrtps256(r2);
#ifdef RSQRT_NR_EPJ_X2
				v8sf v3p0 = {3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f}; 
				ri1 *= (v3p0 - r2*(ri1*ri1));
#endif
				v8sf mri1 = mj * ri1;
				v8sf ri2  = ri1 * ri1;
				v8sf mri3 = mri1 * ri2;

				xj =  __builtin_ia32_shufps256(jbuf, jbuf, 0x00);
				yj =  __builtin_ia32_shufps256(jbuf, jbuf, 0x55);
				zj =  __builtin_ia32_shufps256(jbuf, jbuf, 0xaa);
				mj =  __builtin_ia32_shufps256(jbuf, jbuf, 0xff);

				pot -= mri1;
				ax += mri3 * dx;
				ay += mri3 * dy;
				az += mri3 * dz;
			}
#ifdef RSQRT_NR_EPJ_X2
			v8sf v0p5 = {0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f}; 
			v8sf v0p125 = {0.125f, 0.125f, 0.125f, 0.125f, 0.125f, 0.125f, 0.125f, 0.125f}; 
			pot *= v0p5;
			ax  *= v0p125;
			ay  *= v0p125;
			az  *= v0p125;
#endif

			*(v8sf *)(accpbuf[i/8][0]) = ax;
			*(v8sf *)(accpbuf[i/8][1]) = ay;
			*(v8sf *)(accpbuf[i/8][2]) = az;
			*(v8sf *)(accpbuf[i/8][3]) = pot;
		}
	}

    __attribute__ ((noinline))
    void kernel_epj_nounroll_64bit(const int ni, const int nj){
	//const v4df vzero = {0.0, 0.0, 0.0, 0.0};

	const v4df veps2 = {eps2, eps2, eps2, eps2};
	for(int i=0; i<ni; i+=4){
	    const int il = i%8;
	    const v4df xi = *(v4df *)(&xibufd[i/8][0][il]);
	    const v4df yi = *(v4df *)(&xibufd[i/8][1][il]);
	    const v4df zi = *(v4df *)(&xibufd[i/8][2][il]);

	    v4df ax, ay, az, pot;
	    ax = ay = az = pot = (v4df){0.0, 0.0, 0.0, 0.0};

	    v4df jbuf = *((v4df*)epjbufd);
	    v4df xj =  _mm256_permute4x64_pd(jbuf, 0x00);
	    v4df yj =  _mm256_permute4x64_pd(jbuf, 0x55);
	    v4df zj =  _mm256_permute4x64_pd(jbuf, 0xaa);
	    v4df mj =  _mm256_permute4x64_pd(jbuf, 0xff);

	    for(int j=0; j<nj; j++){
            jbuf = *((v4df*)(epjbufd+j+1));
            v4df dx = xj - xi;
            v4df dy = yj - yi;
            v4df dz = zj - zi;
            v4df r2   = ((veps2 + dx*dx) + dy*dy) + dz*dz;
            //v4df mask = _mm256_cmp_pd(vrcrit2, r2, 0x01); // vrcrit2 < r2
            v4df mask = _mm256_cmp_pd(veps2, r2, 0x4); // veps2 != r2
            v4df ri1  = _mm256_and_pd( __builtin_ia32_cvtps2pd256( __builtin_ia32_rsqrtps( __builtin_ia32_cvtpd2ps256(r2))), mask );
#ifdef RSQRT_NR_EPJ_X2
            //x2
            v4df v0p5 = {0.5, 0.5, 0.5, 0.5};
            v4df v3p0 = {3.0, 3.0, 3.0, 3.0};
            ri1 *= (v3p0 - r2*(ri1*ri1))*v0p5;
#elif defined(RSQRT_NR_EPJ_X4)
            // x4
            v4df vone = {1.0, 1.0, 1.0, 1.0};
            v4df v8p0 = {8.0, 8.0, 8.0, 8.0};
            v4df v6p0 = {6.0, 6.0, 6.0, 6.0};
            v4df v5p0 = {5.0, 5.0, 5.0, 5.0};
            v4df v0p0625 = {1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0};
            v4df h = vone - r2*(ri1*ri1);
            ri1 *= vone + h*(v8p0+h*(v6p0+v5p0*h))*v0p0625;
#endif
            v4df mri1 = mj * ri1;
            v4df ri2  = ri1 * ri1;
            v4df mri3 = mri1 * ri2;
		
            xj =  _mm256_permute4x64_pd(jbuf, 0x00);
            yj =  _mm256_permute4x64_pd(jbuf, 0x55);
            zj =  _mm256_permute4x64_pd(jbuf, 0xaa);
            mj =  _mm256_permute4x64_pd(jbuf, 0xff);

            pot -= mri1;
            ax += mri3 * dx;
            ay += mri3 * dy;
            az += mri3 * dz;
	    }
	    *(v4df *)(&accpbufd[i/8][0][il]) = ax;
	    *(v4df *)(&accpbufd[i/8][1][il]) = ay;
	    *(v4df *)(&accpbufd[i/8][2][il]) = az;
	    *(v4df *)(&accpbufd[i/8][3][il]) = pot;
	}
    }

    ///////////////////////
    // with linear cutoff
    __attribute__ ((noinline))
    void kernel_epj_nounroll_for_p3t_with_linear_cutoff(const int ni, const int nj){
        const v8sf vone = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
        const v8sf veps2 = {(float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2};
        const v8sf vr_out  = {(float)r_out,  (float)r_out,  (float)r_out,  (float)r_out,  (float)r_out,  (float)r_out,  (float)r_out,  (float)r_out};
        const v8sf vr_out2  = vr_out * vr_out;
        const v8sf vrcrit2 = {(float)r_crit2, (float)r_crit2, (float)r_crit2, (float)r_crit2, (float)r_crit2, (float)r_crit2, (float)r_crit2, (float)r_crit2};
        const v8sf allbits = _mm256_cmp_ps(vone, vone, 0x00);
        for(int i=0; i<ni; i+=8){
            const v8sf xi = *(v8sf *)(xibuf[i/8][0]);
            const v8sf yi = *(v8sf *)(xibuf[i/8][1]);
            const v8sf zi = *(v8sf *)(xibuf[i/8][2]);
            v8sf ax, ay, az, pot, nngb;
            ax = ay = az = pot = nngb = (v8sf){0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
            v8sf jbuf = __builtin_ia32_vbroadcastf128_ps256((v4sf *)epjbuf);
            v8sf xj =  __builtin_ia32_shufps256(jbuf, jbuf, 0x00);
            v8sf yj =  __builtin_ia32_shufps256(jbuf, jbuf, 0x55);
            v8sf zj =  __builtin_ia32_shufps256(jbuf, jbuf, 0xaa);
            v8sf mj =  __builtin_ia32_shufps256(jbuf, jbuf, 0xff);
            for(int j=0; j<nj; j++){
                jbuf = __builtin_ia32_vbroadcastf128_ps256((v4sf *)(epjbuf + j+1));
                v8sf dx = xj - xi;
                v8sf dy = yj - yi;
                v8sf dz = zj - zi;
                v8sf r2_real   = ((veps2 + dx*dx) + dy*dy) + dz*dz;
		v8sf r2 = _mm256_max_ps(r2_real, vr_out2);
		v8sf ri1  = __builtin_ia32_rsqrtps256(r2);
#ifdef RSQRT_NR_EPJ_X2
                v8sf v3p0 = {3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f}; 
                v8sf v0p5 = {0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f};
                ri1 *= (v3p0 - r2*(ri1*ri1))*v0p5;
#elif defined(RSQRT_NR_EPJ_X4)
                // x4
                v8sf v8p0 = {8.0f, 8.0f, 8.0f, 8.0f, 8.0f, 8.0f, 8.0f, 8.0f};
                v8sf v6p0 = {6.0f, 6.0f, 6.0f, 6.0f, 6.0f, 6.0f, 6.0f, 6.0f};
                v8sf v5p0 = {5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f};
                v8sf v0p0625 = {(float)1.0/16.0, (float)1.0/16.0, (float)1.0/16.0, (float)1.0/16.0};
                v8sf h = vone - r2*(ri1*ri1);
                ri1 *= vone + h*(v8p0+h*(v6p0+v5p0*h))*v0p0625;
#endif
                v8sf ri2 = ri1*ri1;
                v8sf mri1 = mj*ri1;
                v8sf mri3 = mri1 * ri2;
                xj =  __builtin_ia32_shufps256(jbuf, jbuf, 0x00);
                yj =  __builtin_ia32_shufps256(jbuf, jbuf, 0x55);
                zj =  __builtin_ia32_shufps256(jbuf, jbuf, 0xaa);
                mj =  __builtin_ia32_shufps256(jbuf, jbuf, 0xff);
                pot -= mri1;
                ax += mri3 * dx;
                ay += mri3 * dy;
                az += mri3 * dz;
                v8sf mask = _mm256_cmp_ps(vrcrit2, r2_real, 0x01); // for neighbour search
                nngb += _mm256_and_ps( vone, _mm256_xor_ps(mask, allbits) ); // can remove
            }
            *(v8sf *)(accpbuf[i/8][0]) = ax;
            *(v8sf *)(accpbuf[i/8][1]) = ay;
            *(v8sf *)(accpbuf[i/8][2]) = az;
            *(v8sf *)(accpbuf[i/8][3]) = pot;
            *(v8sf *)(accpbuf[i/8][4]) = nngb;
        }
    }

    __attribute__ ((noinline))
    void kernel_epj_64bit_nounroll_for_p3t_with_linear_cutoff(const int ni, const int nj){
        const v4df vone = {1.0, 1.0, 1.0, 1.0};
        const v4df veps2 = {eps2, eps2, eps2, eps2};
        const v4df vr_out  = {r_out,  r_out,  r_out,  r_out};
        const v4df vr_out2  = vr_out * vr_out;
        const v4df vrcrit2 = {r_crit2, r_crit2, r_crit2, r_crit2};
        const v4df allbits = _mm256_cmp_pd(vone, vone, 0x00);
        for(int i=0; i<ni; i+=4){
            const int il = i%8;
            const v4df xi = *(v4df *)(&xibufd[i/8][0][il]);
            const v4df yi = *(v4df *)(&xibufd[i/8][1][il]);
            const v4df zi = *(v4df *)(&xibufd[i/8][2][il]);
            v4df ax, ay, az, pot, nngb;
            ax = ay = az = pot = nngb = (v4df){0.0, 0.0, 0.0, 0.0};
            v4df jbuf = *((v4df*)epjbufd);
            v4df xj =  _mm256_permute4x64_pd(jbuf, 0x00);
            v4df yj =  _mm256_permute4x64_pd(jbuf, 0x55);
            v4df zj =  _mm256_permute4x64_pd(jbuf, 0xaa);
            v4df mj =  _mm256_permute4x64_pd(jbuf, 0xff);
            for(int j=0; j<nj; j++){
                jbuf = *((v4df*)(epjbufd+j+1));
                v4df dx = xj - xi;
                v4df dy = yj - yi;
                v4df dz = zj - zi;
                v4df r2_real   = ((veps2 + dx*dx) + dy*dy) + dz*dz;
                v4df r2 = _mm256_max_pd( r2_real, vr_out2);
                v4df ri1  = __builtin_ia32_cvtps2pd256(__builtin_ia32_rsqrtps( __builtin_ia32_cvtpd2ps256(r2)));
#ifdef RSQRT_NR_EPJ_X2
                //x2
                v4df v0p5 = {0.5, 0.5, 0.5, 0.5};
                v4df v3p0 = {3.0, 3.0, 3.0, 3.0};
                ri1 *= (v3p0 - r2*(ri1*ri1))*v0p5;
#elif defined(RSQRT_NR_EPJ_X4)
                // x4
                v4df v8p0 = {8.0, 8.0, 8.0, 8.0};
                v4df v6p0 = {6.0, 6.0, 6.0, 6.0};
                v4df v5p0 = {5.0, 5.0, 5.0, 5.0};
                v4df v0p0625 = {1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0};
                v4df h = vone - r2*(ri1*ri1);
                ri1 *= vone + h*(v8p0+h*(v6p0+v5p0*h))*v0p0625;
#endif
                v4df ri2 = ri1*ri1;
                v4df mri1 = mj*ri1;
                v4df mri3 = mri1 * ri2;
                xj =  _mm256_permute4x64_pd(jbuf, 0x00);
                yj =  _mm256_permute4x64_pd(jbuf, 0x55);
                zj =  _mm256_permute4x64_pd(jbuf, 0xaa);
                mj =  _mm256_permute4x64_pd(jbuf, 0xff);
                pot -= mri1;
                ax += mri3 * dx;
                ay += mri3 * dy;
                az += mri3 * dz;
                v4df mask = _mm256_cmp_pd(vrcrit2, r2_real, 0x01); // for neighbour search
                nngb += _mm256_and_pd( vone, _mm256_xor_pd(mask, allbits) ); // can remove
            }
            *(v4df *)(&accpbufd[i/8][0][il]) = ax;
            *(v4df *)(&accpbufd[i/8][1][il]) = ay;
            *(v4df *)(&accpbufd[i/8][2][il]) = az;
            *(v4df *)(&accpbufd[i/8][3][il]) = pot;
            *(v4df *)(&accpbufd[i/8][4][il]) = nngb;
        }
    }
    // linear cutoff
    //////////////

    //////////////
    // includ cutoff
    __attribute__ ((noinline))
    void kernel_epj_nounroll_for_p3t(const int ni, const int nj){
        const v8sf vzero = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
        const v8sf vone = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
        const v8sf veps2 = {(float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2};
        const v8sf vr_in  = {(float)r_in,  (float)r_in,  (float)r_in,  (float)r_in,  (float)r_in,  (float)r_in,  (float)r_in,  (float)r_in};
        const v8sf vr_in2  = vr_in * vr_in;
        const v8sf vrcrit2 = {(float)r_crit2, (float)r_crit2, (float)r_crit2, (float)r_crit2, (float)r_crit2, (float)r_crit2, (float)r_crit2, (float)r_crit2};
        const v8sf vdenominator = {(float)denominator, (float)denominator, (float)denominator, (float)denominator, (float)denominator, (float)denominator, (float)denominator, (float)denominator}; 
        const v8sf vn20 = {-20.0f, -20.0f, -20.0f, -20.0f, -20.0f, -20.0f, -20.0f, -20.0f};
        const v8sf v70  = { 70.0f,  70.0f,  70.0f,  70.0f,  70.0f,  70.0f,  70.0f,  70.0f};
        const v8sf vn84 = {-84.0f, -84.0f, -84.0f, -84.0f, -84.0f, -84.0f, -84.0f, -84.0f}; 
        const v8sf v35  = { 35.0f,  35.0f,  35.0f,  35.0f,  35.0f,  35.0f,  35.0f,  35.0f}; 
        const v8sf allbits = _mm256_cmp_ps(vone, vone, 0x00);
        for(int i=0; i<ni; i+=8){
            const v8sf xi = *(v8sf *)(xibuf[i/8][0]);
            const v8sf yi = *(v8sf *)(xibuf[i/8][1]);
            const v8sf zi = *(v8sf *)(xibuf[i/8][2]);
            v8sf ax, ay, az, pot, nngb;
            ax = ay = az = pot = nngb = (v8sf){0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
            v8sf jbuf = __builtin_ia32_vbroadcastf128_ps256((v4sf *)epjbuf);
            v8sf xj =  __builtin_ia32_shufps256(jbuf, jbuf, 0x00);
            v8sf yj =  __builtin_ia32_shufps256(jbuf, jbuf, 0x55);
            v8sf zj =  __builtin_ia32_shufps256(jbuf, jbuf, 0xaa);
            v8sf mj =  __builtin_ia32_shufps256(jbuf, jbuf, 0xff);
            for(int j=0; j<nj; j++){
                jbuf = __builtin_ia32_vbroadcastf128_ps256((v4sf *)(epjbuf + j+1));
                v8sf dx = xj - xi;
                v8sf dy = yj - yi;
                v8sf dz = zj - zi;
                v8sf r2   = ((veps2 + dx*dx) + dy*dy) + dz*dz;
                v8sf mask = _mm256_cmp_ps(vr_in2, r2, 0x01); // 1op, can remove?
                v8sf ri1  = __builtin_ia32_rsqrtps256(r2);
                ri1  = _mm256_and_ps(ri1, mask); // 1op
#ifdef RSQRT_NR_EPJ_X2
                v8sf v3p0 = {3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f}; 
                v8sf v0p5 = {0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f};
                ri1 *= (v3p0 - r2*(ri1*ri1))*v0p5;
#elif defined(RSQRT_NR_EPJ_X4)
                // x4
                v8sf v8p0 = {8.0f, 8.0f, 8.0f, 8.0f, 8.0f, 8.0f, 8.0f, 8.0f};
                v8sf v6p0 = {6.0f, 6.0f, 6.0f, 6.0f, 6.0f, 6.0f, 6.0f, 6.0f};
                v8sf v5p0 = {5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f};
                v8sf v0p0625 = {(float)1.0/16.0, (float)1.0/16.0, (float)1.0/16.0, (float)1.0/16.0};
                v8sf h = vone - r2*(ri1*ri1);
                ri1 *= vone + h*(v8p0+h*(v6p0+v5p0*h))*v0p0625;
#endif
                v8sf r1 = r2*ri1; // 1op
                v8sf x = _mm256_min_ps( _mm256_max_ps( ((r1-vr_in)*vdenominator), vzero), vone); // 4op
                v8sf x2 = x * x; // 1op
                v8sf x4 = x2 * x2; // 1op
                v8sf k = _mm256_fmadd_ps( _mm256_fmadd_ps( _mm256_fmadd_ps(vn20, x, v70), x, vn84), x, v35)*x4; // 4op
                v8sf mkri1 = mj * k * ri1 ;
                v8sf ri2  = ri1 * ri1;
                v8sf mkri3 = mkri1 * ri2;
                xj =  __builtin_ia32_shufps256(jbuf, jbuf, 0x00);
                yj =  __builtin_ia32_shufps256(jbuf, jbuf, 0x55);
                zj =  __builtin_ia32_shufps256(jbuf, jbuf, 0xaa);
                mj =  __builtin_ia32_shufps256(jbuf, jbuf, 0xff);
                pot -= mkri1;
                ax += mkri3 * dx;
                ay += mkri3 * dy;
                az += mkri3 * dz;
                mask = _mm256_cmp_ps(vrcrit2, r2, 0x01); // for neighbour search
                nngb += _mm256_and_ps( vone, _mm256_xor_ps(mask, allbits) ); // can remove
            }
            *(v8sf *)(accpbuf[i/8][0]) = ax;
            *(v8sf *)(accpbuf[i/8][1]) = ay;
            *(v8sf *)(accpbuf[i/8][2]) = az;
            *(v8sf *)(accpbuf[i/8][3]) = pot;
            *(v8sf *)(accpbuf[i/8][4]) = nngb;
        }
    }

    __attribute__ ((noinline))
    void kernel_epj_64bit_nounroll_for_p3t(const int ni, const int nj){
	const v4df vzero = {0.0, 0.0, 0.0, 0.0};
	const v4df vone = {1.0, 1.0, 1.0, 1.0};
	const v4df veps2 = {eps2, eps2, eps2, eps2};
	const v4df vr_in  = {r_in,  r_in,  r_in,  r_in};
	const v4df vr_in2  = vr_in * vr_in;
	const v4df vrcrit2 = {r_crit2, r_crit2, r_crit2, r_crit2};
	const v4df vdenominator = {denominator, denominator, denominator, denominator};
	const v4df vn20 = {-20.0, -20.0, -20.0, -20.0};
	const v4df v70  = { 70.0,  70.0,  70.0,  70.0};
	const v4df vn84 = {-84.0, -84.0, -84.0, -84.0};
	const v4df v35  = { 35.0,  35.0,  35.0,  35.0};
	const v4df allbits = _mm256_cmp_pd(vone, vone, 0x00);
	for(int i=0; i<ni; i+=4){
	    const int il = i%8;
	    const v4df xi = *(v4df *)(&xibufd[i/8][0][il]);
	    const v4df yi = *(v4df *)(&xibufd[i/8][1][il]);
	    const v4df zi = *(v4df *)(&xibufd[i/8][2][il]);
	    v4df ax, ay, az, pot, nngb;
	    ax = ay = az = pot = nngb = (v4df){0.0, 0.0, 0.0, 0.0};
	    v4df jbuf = *((v4df*)epjbufd);
	    v4df xj =  _mm256_permute4x64_pd(jbuf, 0x00);
	    v4df yj =  _mm256_permute4x64_pd(jbuf, 0x55);
	    v4df zj =  _mm256_permute4x64_pd(jbuf, 0xaa);
	    v4df mj =  _mm256_permute4x64_pd(jbuf, 0xff);
	    for(int j=0; j<nj; j++){
		jbuf = *((v4df*)(epjbufd+j+1));
		v4df dx = xj - xi;
		v4df dy = yj - yi;
		v4df dz = zj - zi;
		v4df r2   = ((veps2 + dx*dx) + dy*dy) + dz*dz;
		v4df mask = _mm256_cmp_pd(vr_in2, r2, 0x01);
		v4df ri1  = _mm256_and_pd( __builtin_ia32_cvtps2pd256( __builtin_ia32_rsqrtps( __builtin_ia32_cvtpd2ps256(r2))), mask );
		ri1  = _mm256_and_pd(ri1, mask); // 1op
#ifdef RSQRT_NR_EPJ_X2
		//x2
		v4df v0p5 = {0.5, 0.5, 0.5, 0.5};
		v4df v3p0 = {3.0, 3.0, 3.0, 3.0};
		ri1 *= (v3p0 - r2*(ri1*ri1))*v0p5;
#elif defined(RSQRT_NR_EPJ_X4)
		// x4
		v4df v8p0 = {8.0, 8.0, 8.0, 8.0};
		v4df v6p0 = {6.0, 6.0, 6.0, 6.0};
		v4df v5p0 = {5.0, 5.0, 5.0, 5.0};
		v4df v0p0625 = {1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0};
		v4df h = vone - r2*(ri1*ri1);
		ri1 *= vone + h*(v8p0+h*(v6p0+v5p0*h))*v0p0625;
#endif
		v4df r1 = r2*ri1; // 1op
		v4df x = _mm256_min_pd( _mm256_max_pd( ((r1-vr_in)*vdenominator), vzero), vone); // 4op
		v4df x2 = x * x; // 1op
		v4df x4 = x2 * x2; // 1op
		v4df k = _mm256_fmadd_pd( _mm256_fmadd_pd( _mm256_fmadd_pd(vn20, x, v70), x, vn84), x, v35)*x4; // 4op
		v4df mkri1 = mj * k * ri1 ;
		v4df ri2  = ri1 * ri1;
		v4df mkri3 = mkri1 * ri2;
		
		xj =  _mm256_permute4x64_pd(jbuf, 0x00);
		yj =  _mm256_permute4x64_pd(jbuf, 0x55);
		zj =  _mm256_permute4x64_pd(jbuf, 0xaa);
		mj =  _mm256_permute4x64_pd(jbuf, 0xff);
		pot -= mkri1;
		ax += mkri3 * dx;
		ay += mkri3 * dy;
		az += mkri3 * dz;
		mask = _mm256_cmp_pd(vrcrit2, r2, 0x01); // for neighbour search
		nngb += _mm256_and_pd( vone, _mm256_xor_pd(mask, allbits) ); // can remove
	    }
	    *(v4df *)(&accpbufd[i/8][0][il]) = ax;
	    *(v4df *)(&accpbufd[i/8][1][il]) = ay;
	    *(v4df *)(&accpbufd[i/8][2][il]) = az;
	    *(v4df *)(&accpbufd[i/8][3][il]) = pot;
	    *(v4df *)(&accpbufd[i/8][4][il]) = nngb;
	}
    }

    /*
    __attribute__ ((noinline))
    void kernel_epj_mix_nounroll_for_p3t(const int ni, const int nj){
	const v8sf vzero = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
	const v8sf vone = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
	const v8sf veps2 = {(float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2};
	const v8sf vr_in  = {(float)r_in,  (float)r_in,  (float)r_in,  (float)r_in,  (float)r_in,  (float)r_in,  (float)r_in,  (float)r_in};
	const v8sf vr_in2  = vr_in * vr_in;
	const v8sf vrcrit2 = {(float)r_crit2, (float)r_crit2, (float)r_crit2, (float)r_crit2, (float)r_crit2, (float)r_crit2, (float)r_crit2, (float)r_crit2};
	const v8sf vdenominator = {(float)denominator, (float)denominator, (float)denominator, (float)denominator, (float)denominator, (float)denominator, (float)denominator, (float)denominator}; 
	const v8sf vn20 = {-20.0f, -20.0f, -20.0f, -20.0f, -20.0f, -20.0f, -20.0f, -20.0f};
	const v8sf v70  = { 70.0f,  70.0f,  70.0f,  70.0f,  70.0f,  70.0f,  70.0f,  70.0f};
	const v8sf vn84 = {-84.0f, -84.0f, -84.0f, -84.0f, -84.0f, -84.0f, -84.0f, -84.0f}; 
	const v8sf v35  = { 35.0f,  35.0f,  35.0f,  35.0f,  35.0f,  35.0f,  35.0f,  35.0f}; 
	const v8sf allbits = _mm256_cmp_ps(vone, vone, 0x00);
	for(int i=0; i<ni; i+=8){
	    const v4df xi0 = *(v4df *)(&xibufd[i/8][0][0]);
	    const v4df yi0 = *(v4df *)(&xibufd[i/8][1][0]);
	    const v4df zi0 = *(v4df *)(&xibufd[i/8][2][0]);

	    const v4df xi1 = *(v4df *)(&xibufd[i/8][0][1]);
	    const v4df yi1 = *(v4df *)(&xibufd[i/8][1][1]);
	    const v4df zi1 = *(v4df *)(&xibufd[i/8][2][1]);

	    //const v8sf xi = *(v8sf *)(xibuf[i/8][0]);
	    //const v8sf yi = *(v8sf *)(xibuf[i/8][1]);
	    //const v8sf zi = *(v8sf *)(xibuf[i/8][2]);
	    //v8sf ax, ay, az, pot, nngb;
	    //ax = ay = az = pot = nngb = (v8sf){0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};

	    v4df ax0, ay0, az0, pot0; 
	    ax = ay = az = pot = (v4df){0.0, 0.0, 0.0, 0.0};
	    v8sf nngb;
	    nngb = (v8sf){0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};

	    //v8sf jbuf = __builtin_ia32_vbroadcastf128_ps256((v4sf *)epjbuf);
	    //v8sf xj =  __builtin_ia32_shufps256(jbuf, jbuf, 0x00);
	    //v8sf yj =  __builtin_ia32_shufps256(jbuf, jbuf, 0x55);
	    //v8sf zj =  __builtin_ia32_shufps256(jbuf, jbuf, 0xaa);
	    //v8sf mj =  __builtin_ia32_shufps256(jbuf, jbuf, 0xff);
	    v4df jbuf = *((v4df*)epjbufd);
	    v4df xj =  _mm256_permute4x64_pd(jbuf, 0x00);
	    v4df yj =  _mm256_permute4x64_pd(jbuf, 0x55);
	    v4df zj =  _mm256_permute4x64_pd(jbuf, 0xaa);
	    v4df mj =  _mm256_permute4x64_pd(jbuf, 0xff);
	    for(int j=0; j<nj; j++){
		//jbuf = __builtin_ia32_vbroadcastf128_ps256((v4sf *)(epjbuf + j+1));
		jbuf = *((v4df*)(epjbufd+j+1));
		//v8sf dx = xj - xi;
		//v8sf dy = yj - yi;
		//v8sf dz = zj - zi;
		v4df dx0 = xj - xi0;
		v4df dy0 = yj - yi0;
		v4df dz0 = zj - zi0;
		v4df dx1 = xj - xi1;
		v4df dy1 = yj - yi1;
		v4df dz1 = zj - zi1;
		v8sf dx0s = _mm256_cvtpd_ps(dx0);
		v8sf dx1s = _mm256_cvtpd_ps(dx1);
		v8sf dx = __mm256_permute2f128_ps(dx0s, dx1s, 0x20);
		v8sf dy0s = _mm256_cvtpd_ps(dy0);
		v8sf dy1s = _mm256_cvtpd_ps(dy1);
		v8sf dy = __mm256_permute2f128_ps(dy0s, dy1s, 0x20);
		v8sf dz0s = _mm256_cvtpd_ps(dz0);
		v8sf dz1s = _mm256_cvtpd_ps(dz1);
		v8sf dz = __mm256_permute2f128_ps(dz0s, dz1s, 0x20);
		v8sf r2   = ((veps2 + dx*dx) + dy*dy) + dz*dz;
		v8sf mask = _mm256_cmp_ps(vr_in2, r2, 0x01); // 1op, can remove?
		v8sf ri1  = __builtin_ia32_rsqrtps256(r2);
		ri1  = _mm256_and_ps(ri1, mask); // 1op
#ifdef RSQRT_NR_EPJ_X2
		v8sf v3p0 = {3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f}; 
		v8sf v0p5 = {0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f};
		ri1 *= (v3p0 - r2*(ri1*ri1))*v0p5;
#elif defined(RSQRT_NR_EPJ_X4)
		// x4
		v8sf v8p0 = {8.0f, 8.0f, 8.0f, 8.0f, 8.0f, 8.0f, 8.0f, 8.0f};
		v8sf v6p0 = {6.0f, 6.0f, 6.0f, 6.0f, 6.0f, 6.0f, 6.0f, 6.0f};
		v8sf v5p0 = {5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f};
		v8sf v0p0625 = {(float)1.0/16.0, (float)1.0/16.0, (float)1.0/16.0, (float)1.0/16.0};
		v8sf h = vone - r2*(ri1*ri1);
		ri1 *= vone + h*(v8p0+h*(v6p0+v5p0*h))*v0p0625;
#endif
		v8sf r1 = r2*ri1; // 1op
		v8sf x = _mm256_min_ps( _mm256_max_ps( ((r1-vr_in)*vdenominator), vzero), vone); // 4op
		v8sf x2 = x * x; // 1op
		v8sf x4 = x2 * x2; // 1op
		v8sf k = _mm256_fmadd_ps( _mm256_fmadd_ps( _mm256_fmadd_ps(vn20, x, v70), x, vn84), x, v35)*x4; // 4op
		v8sf mkri1 = mj * k * ri1 ;
		v8sf ri2  = ri1 * ri1;
		v8sf mkri3 = mkri1 * ri2;
		xj =  _mm256_permute4x64_pd(jbuf, 0x00);
		yj =  _mm256_permute4x64_pd(jbuf, 0x55);
		zj =  _mm256_permute4x64_pd(jbuf, 0xaa);
		mj =  _mm256_permute4x64_pd(jbuf, 0xff);
		//xj =  __builtin_ia32_shufps256(jbuf, jbuf, 0x00);
		//yj =  __builtin_ia32_shufps256(jbuf, jbuf, 0x55);
		//zj =  __builtin_ia32_shufps256(jbuf, jbuf, 0xaa);
		//mj =  __builtin_ia32_shufps256(jbuf, jbuf, 0xff);
		pot -= mkri1;
		ax += mkri3 * dx;
		ay += mkri3 * dy;
		az += mkri3 * dz;
		mask = _mm256_cmp_ps(vrcrit2, r2, 0x01); // for neighbour search
		nngb += _mm256_and_ps( vone, _mm256_xor_ps(mask, allbits) ); // can remove
	    }
	    *(v8sf *)(accpbuf[i/8][0]) = ax;
	    *(v8sf *)(accpbuf[i/8][1]) = ay;
	    *(v8sf *)(accpbuf[i/8][2]) = az;
	    *(v8sf *)(accpbuf[i/8][3]) = pot;
	    *(v8sf *)(accpbuf[i/8][4]) = nngb;
	}
    }
    */

	__attribute__ ((noinline))
	void kernel_spj_nounroll(const int ni, const int nj){
	    const v8sf veps2 = {(float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2};
		for(int i=0; i<ni; i+=8){
			const v8sf xi = *(v8sf *)(xibuf[i/8][0]);
			const v8sf yi = *(v8sf *)(xibuf[i/8][1]);
			const v8sf zi = *(v8sf *)(xibuf[i/8][2]);

			v8sf ax, ay, az, pot;
			ax = ay = az = pot = (v8sf){0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};

#define PRELOAD_SPJ

#ifdef PRELOAD_SPJ
			v8sf jbuf0 = __builtin_ia32_vbroadcastf128_ps256((v4sf *)&spjbuf[0][0]);
			v8sf jbuf1 = __builtin_ia32_vbroadcastf128_ps256((v4sf *)&spjbuf[0][1]);
			v8sf jbuf2 = __builtin_ia32_vbroadcastf128_ps256((v4sf *)&spjbuf[0][2]);
#else
			v8sf jbuf0, jbuf1, jbuf2;
#endif
			for(int j=0; j<nj; j++){
#ifndef PRELOAD_SPJ
				jbuf0 = __builtin_ia32_vbroadcastf128_ps256((v4sf *)&spjbuf[j+0][0]);
#endif
				v8sf xj  = __builtin_ia32_shufps256(jbuf0, jbuf0, 0x00);
				v8sf yj  = __builtin_ia32_shufps256(jbuf0, jbuf0, 0x55);
				v8sf zj  = __builtin_ia32_shufps256(jbuf0, jbuf0, 0xaa);
#ifdef PRELOAD_SPJ
				jbuf0 = __builtin_ia32_vbroadcastf128_ps256((v4sf *)&spjbuf[j+1][0]);
#endif

#ifndef PRELOAD_SPJ
				jbuf1 = __builtin_ia32_vbroadcastf128_ps256((v4sf *)&spjbuf[j+0][1]);
#endif
				v8sf qxx = __builtin_ia32_shufps256(jbuf1, jbuf1, 0x00);
				v8sf qyy = __builtin_ia32_shufps256(jbuf1, jbuf1, 0x55);
				v8sf qzz = __builtin_ia32_shufps256(jbuf1, jbuf1, 0xaa);
				v8sf mj  = __builtin_ia32_shufps256(jbuf1, jbuf1, 0xff);
#ifdef PRELOAD_SPJ
				jbuf1 = __builtin_ia32_vbroadcastf128_ps256((v4sf *)&spjbuf[j+1][1]);
#endif

#ifndef PRELOAD_SPJ
				jbuf2 = __builtin_ia32_vbroadcastf128_ps256((v4sf *)&spjbuf[j+0][2]);
#endif
				v8sf qxy = __builtin_ia32_shufps256(jbuf2, jbuf2, 0x00);
				v8sf qyz = __builtin_ia32_shufps256(jbuf2, jbuf2, 0x55);
				v8sf qzx = __builtin_ia32_shufps256(jbuf2, jbuf2, 0xaa);
				v8sf mtr = __builtin_ia32_shufps256(jbuf2, jbuf2, 0xff);
#ifdef PRELOAD_SPJ
				jbuf2 = __builtin_ia32_vbroadcastf128_ps256((v4sf *)&spjbuf[j+1][2]);
#endif

				v8sf dx = xj - xi;
				v8sf dy = yj - yi;
				v8sf dz = zj - zi;

				v8sf r2  = ((veps2 + dx*dx) + dy*dy) + dz*dz;
				v8sf ri1 = __builtin_ia32_rsqrtps256(r2);
				v8sf v0p5 = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
#ifdef RSQRT_NR_SPJ_X2
                                v8sf v3p0 = {3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0};
                                ri1 *= (v3p0 - r2*(ri1*ri1))*v0p5;
#endif
				v8sf ri2 = ri1 * ri1;
				v8sf ri3 = ri1 * ri2;
				v8sf ri4 = ri2 * ri2;
				v8sf ri5 = ri2 * ri3;

				v8sf qr_x = (qxx*dx + qxy*dy) + qzx*dz;
				v8sf qr_y = (qyy*dy + qxy*dx) + qyz*dz;
				v8sf qr_z = (qzz*dz + qzx*dx) + qyz*dy;

				v8sf rqr = ((mtr + qr_x*dx) + qr_y*dy) + qr_z*dz;
				v8sf rqr_ri4 = rqr * ri4;

				//v8sf v0p5 = {0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f};
				v8sf v2p5 = {2.5f, 2.5f, 2.5f, 2.5f, 2.5f, 2.5f, 2.5f, 2.5f}; 

				v8sf meff  =  mj + v0p5 * rqr_ri4;
				v8sf meff3 = (mj + v2p5 * rqr_ri4) * ri3;

				pot -= meff * ri1;

				ax = (ax - ri5*qr_x) + meff3*dx;
				ay = (ay - ri5*qr_y) + meff3*dy;
				az = (az - ri5*qr_z) + meff3*dz;
			}
			*(v8sf *)(accpbuf[i/8][0]) = ax;
			*(v8sf *)(accpbuf[i/8][1]) = ay;
			*(v8sf *)(accpbuf[i/8][2]) = az;
			*(v8sf *)(accpbuf[i/8][3]) = pot;
		}
	}

    __attribute__ ((noinline))
    void kernel_spj_64bit_nounroll(const int ni, const int nj){
	const v4df veps2 = {eps2, eps2, eps2, eps2};
	for(int i=0; i<ni; i+=4){
	    const int il = i%8;
	    const v4df xi = *(v4df *)(&xibufd[i/8][0][il]);
	    const v4df yi = *(v4df *)(&xibufd[i/8][1][il]);
	    const v4df zi = *(v4df *)(&xibufd[i/8][2][il]);

	    v4df ax, ay, az, pot;
	    ax = ay = az = pot = (v4df){0.0, 0.0, 0.0, 0.0};

	    v4df jbuf0 = *((v4df*)spjbufd[0][0]);
	    v4df jbuf1 = *((v4df*)spjbufd[0][1]);
	    v4df jbuf2 = *((v4df*)spjbufd[0][2]);

	    for(int j=0; j<nj; j++){
		v4df xj  = _mm256_permute4x64_pd(jbuf0, 0x00);
		v4df yj  = _mm256_permute4x64_pd(jbuf0, 0x55);
		v4df zj  = _mm256_permute4x64_pd(jbuf0, 0xaa);
		jbuf0 = *((v4df*)spjbufd[j+1][0]);

		v4df qxx = _mm256_permute4x64_pd(jbuf1, 0x00);
		v4df qyy = _mm256_permute4x64_pd(jbuf1, 0x55);
		v4df qzz = _mm256_permute4x64_pd(jbuf1, 0xaa);
		v4df mj  = _mm256_permute4x64_pd(jbuf1, 0xff);
		jbuf1 = *((v4df*)spjbufd[j+1][1]);
		
		v4df qxy = _mm256_permute4x64_pd(jbuf2, 0x00);
		v4df qyz = _mm256_permute4x64_pd(jbuf2, 0x55);
		v4df qzx = _mm256_permute4x64_pd(jbuf2, 0xaa);
		v4df mtr = _mm256_permute4x64_pd(jbuf2, 0xff);
		jbuf2 = *((v4df*)spjbufd[j+1][2]);
		
		v4df dx = xj - xi;
		v4df dy = yj - yi;
		v4df dz = zj - zi;
		
		v4df r2  = ((veps2 + dx*dx) + dy*dy) + dz*dz;
		//v4df ri1 = __builtin_ia32_rsqrtps256(r2);
		v4df ri1 = __builtin_ia32_cvtps2pd256( __builtin_ia32_rsqrtps( __builtin_ia32_cvtpd2ps256(r2)));
		
#ifdef RSQRT_NR_SPJ_X2
		//x2
		v4df v3p0 = {3.0, 3.0, 3.0, 3.0};
		ri1 *= (v3p0 - r2*(ri1*ri1));
#elif defined(RSQRT_NR_SPJ_X4)
		// x4
		v4df v8p0 = {8.0, 8.0, 8.0, 8.0};
		v4df v6p0 = {6.0, 6.0, 6.0, 6.0};
		v4df v5p0 = {5.0, 5.0, 5.0, 5.0};
		v4df v0p0625 = {1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0};
		v4df v1p0 = {1.0, 1.0, 1.0, 1.0};
		v4df h = v1p0 - r2*(ri1*ri1);
		ri1 *= v1p0 + h*(v8p0+h*(v6p0+v5p0*h))*v0p0625;
#endif
		v4df ri2 = ri1 * ri1;
		v4df ri3 = ri1 * ri2;
		v4df ri4 = ri2 * ri2;
		v4df ri5 = ri2 * ri3;
		
		v4df qr_x = (qxx*dx + qxy*dy) + qzx*dz;
		v4df qr_y = (qyy*dy + qxy*dx) + qyz*dz;
		v4df qr_z = (qzz*dz + qzx*dx) + qyz*dy;
		
		v4df rqr = ((mtr + qr_x*dx) + qr_y*dy) + qr_z*dz;
		v4df rqr_ri4 = rqr * ri4;
		
		v4df v0p5 = {0.5, 0.5, 0.5, 0.5};
		v4df v2p5 = {2.5, 2.5, 2.5, 2.5};
		
		v4df meff  =  mj + v0p5 * rqr_ri4;
		v4df meff3 = (mj + v2p5 * rqr_ri4) * ri3;

		pot -= meff * ri1;

		ax = (ax - ri5*qr_x) + meff3*dx;
		ay = (ay - ri5*qr_y) + meff3*dy;
		az = (az - ri5*qr_z) + meff3*dz;
	    }

#ifdef RSQRT_NR_SPJ_X2
	    //x2
	    v4df v0p5 = {0.5, 0.5, 0.5, 0.5};
	    v4df v0p125 = {0.125, 0.125, 0.125, 0.125};
	    pot *= v0p5;
	    ax  *= v0p125;
	    ay  *= v0p125;
	    az  *= v0p125;
#endif
	    
	    *(v4df *)(&accpbufd[i/8][0][il]) = ax;
	    *(v4df *)(&accpbufd[i/8][1][il]) = ay;
	    *(v4df *)(&accpbufd[i/8][2][il]) = az;
	    *(v4df *)(&accpbufd[i/8][3][il]) = pot;
	}
    }


} __attribute__ ((aligned(128)));







