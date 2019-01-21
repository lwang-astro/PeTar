#pragma once
#ifdef INTRINSIC_K
#include"phantomquad_for_p3t_k.hpp"
#endif
#ifdef INTRINSIC_X86
#include"phantomquad_for_p3t_x86.hpp"
#endif

#include"ptcl.hpp"


class ForceSoft{
public:
    PS::F64vec acc; ///> soft acceleration (c.m.: averaged force from orbital particles; tensor: c.m. is substracted)
#ifdef KDKDK_4TH
    PS::F64vec acorr; ///> soft gradient correction for 4th order KDKDK method
#endif
    PS::F64 pot; ///> full potential
    PS::S32 n_ngb; // neighbor number+1
    void clear(){
        acc = 0.0;
        pot = 0.0;
        n_ngb = 0;
    }
};

class FPSoft: public Ptcl{
public:
#ifdef KDKDK_4TH
    PS::F64vec acorr;
#else
    PS::F64vec acc; // soft
#endif
    PS::F64 pot_tot; // soft + hard
    PS::S32 n_ngb;
    PS::S32 rank_org;
    PS::S32 adr;
    static PS::F64 r_out;

    FPSoft() {}

    //! Get position (required for \ref ARC::chain)
    /*! \return position vector (PS::F64[3])
     */
    PS::F64vec getPos() const{
        return pos;
    }

    //! Get velocity (required for \ref ARC::chain)
    /*! \return velocity vector (PS::F64[3])
     */
    PS::F64vec getVel() const{
        return vel;
    }
    
    //template<class Tptcl>
    //FPSoft& operator = (const FPSoft& p) {
    //    Ptcl::DataCopy(p);
    //    acc = p.acc;
    //    pot_tot = p.pot_tot;
    //    n_ngb = p.n_ngb;
    //    return *this;
    //}

    template<class Tptcl>
    FPSoft(const Tptcl& p, const PS::S32 rank_, const PS::S32 adr_) {
        Ptcl::DataCopy(p);
        rank_org = rank_;
        adr = adr_;
        acc = 0;
        pot_tot = 0;
        n_ngb = 0;
    }

    void copyFromForce(const ForceSoft & force){
        acc = force.acc;
#ifdef KDKDK_4TH
        acorr = force.acorr;
#endif
        pot_tot = force.pot;
        n_ngb = force.n_ngb;
    }

    PS::F64 getRSearch() const {
//        return std::max(r_search * std::pow(2.0*mass/m_average,0.3333),r_search_min) * SAFTY_FACTOR_FOR_SEARCH + SAFTY_OFFSET_FOR_SEARCH;
//        return this->r_search * SAFTY_FACTOR_FOR_SEARCH + SAFTY_OFFSET_FOR_SEARCH;
        return this->r_search*SAFTY_FACTOR_FOR_SEARCH;
    }

    void writeAscii(FILE* fp) const{
        Ptcl::writeAscii(fp);
        fprintf(fp, "%26.17e %26.17e %26.17e %26.17e %d\n", 
                this->acc.x, this->acc.y, this->acc.z,  // 9-11
                this->pot_tot, this->n_ngb);
    }

    void writeBinary(FILE* fp) const{
        Ptcl::writeBinary(fp);
        fwrite(&(this->acc), sizeof(PS::F64), 4, fp);
        fwrite(&(this->n_ngb), sizeof(PS::S32), 1, fp);
    }

    void readAscii(FILE* fp) {
        Ptcl::readAscii(fp);
        PS::S64 rcount=fscanf(fp, "%lf %lf %lf %lf %d\n",
                              &this->acc.x, &this->acc.y, &this->acc.z,  // 9-11
                              &this->pot_tot, &this->n_ngb);
        if (rcount<5) {
            std::cerr<<"Error: Data reading fails! requiring data number is 5, only obtain "<<rcount<<".\n";
            abort();
        }
    }

    void readBinary(FILE* fp) {
        Ptcl::readBinary(fp);
        size_t rcount = fread(&(this->acc), sizeof(PS::F64), 4, fp);
        rcount += fread(&(this->n_ngb), sizeof(PS::S32), 1, fp);
        if (rcount<5) {
            std::cerr<<"Error: Data reading fails! requiring data number is 5, only obtain "<<rcount<<".\n";
            abort();
        }
    }

    void print(std::ostream & fout){
        Ptcl::print(fout);
        fout<<" adr= "<<adr
            <<" acc= "<<acc
            <<" pot_tot= "<<pot_tot;
    }
    
};

class EPISoft{
public:
    PS::S64 id;
    PS::F64vec pos;
    PS::F64 r_search;
    PS::S32 rank_org;
#ifdef KDKDK_4TH
    PS::F64vec acc;
#endif
    static PS::F64 eps;
    static PS::F64 r_out;
    static PS::F64 r_in;
//    static PS::F64 r_search;
//    static PS::F64 r_in;
//    static PS::F64 m_average;
//    static PS::F64 r_search_min;
    PS::F64vec getPos() const { return pos;}
    void copyFromFP(const FPSoft & fp){ 
        id = fp.id;
        pos = fp.pos;
#ifdef KDKDK_4TH
        acc = fp.acc;
#endif
        r_search = fp.r_search;
        rank_org = fp.rank_org;
    }
    void print(std::ostream & fout=std::cout) const {
        fout<<" id="<<id
            <<" rank_org="<<rank_org
            <<" pos="<<pos
            <<" eps="<<eps;
    }
    PS::F64 getRSearch() const {
//        return std::max(r_search * std::pow(2.0*mass/m_average,0.3333),r_search_min) * SAFTY_FACTOR_FOR_SEARCH + SAFTY_OFFSET_FOR_SEARCH;
//        return r_search; * SAFTY_FACTOR_FOR_SEARCH + SAFTY_OFFSET_FOR_SEARCH;
        return r_search*SAFTY_FACTOR_FOR_SEARCH;
    }

};


class EPJSoft{
public:
    PS::S64 id;
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64vec vel;
#ifdef KDKDK_4TH
    PS::F64vec acc;
#endif
    PS::F64 r_search;
    PS::F64 mass_bk;
    PS::S64 status;
    PS::S32 rank_org;
    PS::S32 adr_org;
    static PS::F64 r_out;
//    static PS::F64 m_average;
//    static PS::F64 r_search_min;
    void copyFromFP(const FPSoft & fp){
        id = fp.id;
        mass = fp.mass;
        pos = fp.pos;
        vel = fp.vel;
#ifdef KDKDK_4TH
        acc = fp.acc;
#endif
        r_search = fp.r_search;
        mass_bk  = fp.mass_bk;
        status   = fp.status;
        rank_org = fp.rank_org;
        adr_org = fp.adr;
    }
    PS::F64vec getPos() const { return pos; }
    void setPos(const PS::F64vec & pos_new){ pos = pos_new;}
    PS::F64 getCharge() const { return mass; }
    PS::F64 getRSearch() const {
//        return std::max(r_search * std::pow(2.0*mass/m_average,0.3333),r_search_min) * SAFTY_FACTOR_FOR_SEARCH + SAFTY_OFFSET_FOR_SEARCH;
//        return r_search * SAFTY_FACTOR_FOR_SEARCH + SAFTY_OFFSET_FOR_SEARCH;
        return r_search*SAFTY_FACTOR_FOR_SEARCH;
    }
    // FORDEBUG
    void print(std::ostream & fout=std::cout) const {
        fout<<" id="<<id
            <<" rank_org="<<rank_org
            <<" mass="<<mass
            <<" pos="<<pos
            <<" vel="<<vel
            <<" r_search="<<r_search;
    }
    void clear(){
        mass = 0.0;
        pos = vel = 0.0;
        r_search = 0.0;
        id = rank_org = adr_org = -1;
    }
};

PS::F64 EPISoft::eps = 0.0;
PS::F64 EPISoft::r_out = 0.0;
PS::F64 EPISoft::r_in  = 0.0;
PS::F64 EPJSoft::r_out = 0.0;
PS::F64  FPSoft::r_out = 0.0;

//PS::F64 EPISoft::r_in; = 0.0;
//PS::F64 EPJSoft::m_average;
//PS::F64 EPJSoft::r_search_min;
//PS::F64 EPISoft::m_average;
//PS::F64 EPISoft::r_search_min;

class EnergyAndMomemtum{
public:
    PS::F64 kin;
    PS::F64 pot;
    PS::F64 tot;
    PS::F64vec L; // angular momentum
    PS::F64 Lt; // total angular momemtum

    EnergyAndMomemtum() {
        clear();
    }

    void clear(){
        kin = pot = tot = Lt = 0.0;
        L = PS::F64vec(0.0);
    }
    void print(std::ostream & fout=std::cout){
        fout<<"Energy: total = "<<tot<<" kin = "<<kin<<" pot = "<<pot
            <<"\nAngular Momentum: L = "<<L<<" |L| = "<<Lt
            <<std::endl;
    }
    void dumpName(std::ofstream & fout, const PS::S32 width=20) {
        fout<<std::setw(width)<<"Etot"
            <<std::setw(width)<<"Ekin"
            <<std::setw(width)<<"Epot"
            <<std::setw(width)<<"Lx"
            <<std::setw(width)<<"Ly"
            <<std::setw(width)<<"Lz"
            <<std::setw(width)<<"|L|";
    }
    void dump(std::ofstream & fout, const PS::S32 width=20){
        fout<<std::setw(width)<<tot
            <<std::setw(width)<<kin
            <<std::setw(width)<<pot
            <<std::setw(width)<<L.x
            <<std::setw(width)<<L.y
            <<std::setw(width)<<L.z
            <<std::setw(width)<<Lt;
    }
    void writeAscii(FILE* fout) {
        fprintf(fout, "%26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e ",
                this->kin, this->pot, this->tot, this->L[0], this->L[1], this->L[2], this->Lt);
    }
    void writeBinary(FILE* fout) {
        fwrite(&this->kin, sizeof(EnergyAndMomemtum), 1, fout);
    }

    template<class Tptcl>
    void calc(const Tptcl* sys,
              const PS::S32 n) {
        // PS::S32 n = sys.getNumberOfParticleLocal();
        PS::F64 pot_loc = 0.0;
        PS::F64 kin_loc = 0.0;
        PS::F64vec L_loc = PS::F64vec(0.0);
        //#pragma omp parallel for reduction(+:pot_d) reduction (+:pot_loc) reduction(+:kin_loc)
        //PS::ReallocatableArray<PS::S32> plist;
        //plist.reserve(n);
        for(PS::S32 i=0; i<n; i++){
            PS::F64 mi = sys[i].mass;
            if(sys[i].status<0) mi = sys[i].mass_bk;
#ifdef HARD_DEBUG
            assert(sys[i].id>0&&sys[i].status<=0);
            assert(mi>0);
#endif
            PS::F64vec vi = sys[i].vel;
            pot_loc += 0.5 * mi * sys[i].pot_tot;
            kin_loc += 0.5 * mi * vi * vi;
            L_loc += sys[i].pos ^ (mi*vi);
        }
        this->kin += kin_loc;
        this->pot += pot_loc;
        this->L   += L_loc;
        this->Lt   = std::sqrt(L*L);
        this->tot = this->kin + this->pot;
    }

    template<class Tptcl> 
    void calc(const Tptcl* sys,
              const PS::S32* p_list,
              const PS::S32 n) {
        PS::F64 pot_loc = 0.0;
        PS::F64 kin_loc = 0.0;
        PS::F64vec L_loc = PS::F64vec(0.0);
        for(PS::S32 k=0; k<n; k++){
            PS::S32 i = p_list[k];
            PS::F64 mi = sys[i].mass;
            PS::F64vec vi = sys[i].vel;
            pot_loc += 0.5 * mi * sys[i].pot_tot;
            kin_loc += 0.5 * mi * vi * vi;
            L_loc += sys[i].pos ^ (mi*vi);
        }
        this->kin += kin_loc;
        this->pot += pot_loc;
        this->L   += L_loc;
        this->Lt   = std::sqrt(L*L);
        this->tot  = this->kin + this->pot;
    }
              

    void getSumMultiNodes() {
        this->kin = PS::Comm::getSum(this->kin);
        this->pot = PS::Comm::getSum(this->pot);
        this->L   = PS::Comm::getSum(this->L);
        this->Lt  = std::sqrt(L*L);
        this->tot = this->kin + this->pot;
    }

    EnergyAndMomemtum operator -(const EnergyAndMomemtum& eng){
        EnergyAndMomemtum diff;
        diff.kin = this->kin - eng.kin;
        diff.pot = this->pot - eng.pot;
        diff.tot = this->tot - eng.tot;
        diff.L   = this->L   - eng.L;
        diff.Lt  = std::sqrt(diff.L*diff.L);
        return diff;
    }

    void relative(const EnergyAndMomemtum& ref) {
        this->kin /= ref.kin;
        this->pot /= ref.pot;
        this->tot /= ref.tot;
        this->L   /= ref.Lt;
        this->Lt  /= ref.Lt;
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
        //        const PS::F64 r_crit2 = EPJSoft::r_search * EPJSoft::r_search * SAFTY_FACTOR_FOR_SEARCH_SQ;
        // const PS::F64 r_out = EPISoft::r_out; 
        // const PS::F64 r_in = EPISoft::r_in;
        //std::cerr<<"r_out= "<<r_out<<" r_in= "<<r_in<<" eps2= "<<eps2<<" r_crit2= "<<r_crit2<<std::endl;
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
            force[i].acc += ai;
#ifdef KDKDK_4TH
            force[i].acorr = 0.0;
#endif
            force[i].pot += poti;
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
        assert(n_jp==0);
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
            force[i].acc += ai;
            force[i].pot += poti;
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
            force[ip].acc += ai;
            force[ip].pot += poti;
        }
    }
};

#ifdef USE_SIMD
struct CalcForceEpEpWithLinearCutoffSimd{
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
                      const EPJSoft * ep_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
        PS::S32 ep_j_list[n_jp], n_jp_local=0;
        for (PS::S32 i=0; i<n_jp; i++){
            if(ep_j[i].mass>0) ep_j_list[n_jp_local++] = i;
        }
        //std::cerr<<"n_jp="<<n_jp<<" reduced n_jp="<<n_jp_local<<std::endl;
//        const PS::F64 r_crit2 = EPJSoft::r_search * EPJSoft::r_search;
    #ifdef __HPC_ACE__
        PhantomGrapeQuad pg;
    #else
        #if defined(CALC_EP_64bit) || defined(CALC_EP_MIX)
        static __thread PhantomGrapeQuad64Bit pg;
        #else
        static __thread PhantomGrapeQuad pg;
        #endif
    #endif
        if(n_ip > pg.NIMAX || n_jp > pg.NJMAX){
            std::cout<<"ni= "<<n_ip<<" NIMAX= "<<pg.NIMAX<<" nj= "<<n_jp<<" NJMAX= "<<pg.NJMAX<<std::endl;
        }
        assert(n_ip<=pg.NIMAX);
        assert(n_jp<=pg.NJMAX);
        pg.set_eps2(eps2);
        pg.set_r_crit2(EPISoft::r_out*EPISoft::r_out);
        //pg.set_cutoff(EPISoft::r_out, EPISoft::r_in);
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F64vec pos_i = ep_i[i].getPos();
            pg.set_xi_one(i, pos_i.x, pos_i.y, pos_i.z, ep_i[i].r_search);
        }
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
            for(PS::S32 i=0; i<n_ip; i++){
                PS::F64 * p = &(force[i].pot);
                PS::F64 * a = (PS::F64 * )(&force[i].acc[0]);
                PS::F64 n_ngb = 0;
                pg.accum_accp_one(i, a[0], a[1], a[2], *p, n_ngb);
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
#ifdef __HPC_ACE__
        PhantomGrapeQuad pg;
#else
    #if defined(CALC_EP_64bit)
        static __thread PhantomGrapeQuad64Bit pg;
    #else
        static __thread PhantomGrapeQuad pg;
    #endif
#endif
        assert(n_ip<=pg.NIMAX);
        assert(n_jp<=pg.NJMAX);
        pg.set_eps2(eps2);
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F64vec pos_i = ep_i[i].getPos();
            pg.set_xi_one(i, pos_i.x, pos_i.y, pos_i.z, 0.0);
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
            for(PS::S32 i=0; i<n_ip; i++){
                PS::F64 * p = &(force[i].pot);
                PS::F64 * a = (PS::F64 * )(&force[i].acc[0]);
                pg.accum_accp_one(i, a[0], a[1], a[2], *p);
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
    #ifdef __HPC_ACE__
        PhantomGrapeQuad pg;
    #else
        #if defined(CALC_EP_64bit)
        static __thread PhantomGrapeQuad64Bit pg;
        #else
        static __thread PhantomGrapeQuad pg;
        #endif
    #endif
        assert(n_ip<=pg.NIMAX);
        assert(n_jp<=pg.NJMAX);
        pg.set_eps2(eps2);
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F64vec pos_i = ep_i[i].getPos();
            pg.set_xi_one(i, pos_i.x, pos_i.y, pos_i.z, 0.0);
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
            for(PS::S32 i=0; i<n_ip; i++){
                PS::F64 * p = &(force[i].pot);
                PS::F64 * a = (PS::F64 * )(&force[i].acc[0]);
                pg.accum_accp_one(i, a[0], a[1], a[2], *p);
            }
        }
    }
};
#endif
