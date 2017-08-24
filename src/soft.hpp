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
    PS::F64vec acc; // soft
    PS::F64 pot; // soft
    PS::S32 n_ngb;
    void clear(){
        acc = 0.0;
        pot = 0.0;
        n_ngb = 0;
    }
};

class FPSoft: public Ptcl{
public:
    PS::F64vec acc; // soft
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
        pot_tot = force.pot;
        n_ngb = force.n_ngb;
    }

    PS::F64 getRSearch() const {
//        return std::max(r_search * std::pow(2.0*mass/m_average,0.3333),r_search_min) * SAFTY_FACTOR_FOR_SEARCH + SAFTY_OFFSET_FOR_SEARCH;
        return r_out + this->r_search * SAFTY_FACTOR_FOR_SEARCH + SAFTY_OFFSET_FOR_SEARCH;
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

    void dump(std::ofstream & fout){
        Ptcl::dump(fout);
        fout<<"adr= "<<adr<<std::endl;
        fout<<"acc= "<<acc<<std::endl;
        fout<<"pot_tot= "<<pot_tot<<std::endl;
    }

};

class EPISoft{
public:
    PS::S64 id;
    PS::F64vec pos;
    PS::F64 r_search;
    PS::S32 rank_org;
    static PS::F64 eps;
    static PS::F64 r_out;
//    static PS::F64 r_search;
//    static PS::F64 r_in;
//    static PS::F64 m_average;
//    static PS::F64 r_search_min;
    PS::F64vec getPos() const { return pos;}
    void copyFromFP(const FPSoft & fp){ 
        id = fp.id;
        pos = fp.pos;
        r_search = fp.r_search;
        rank_org = fp.rank_org;
    }
    void dump(std::ostream & fout=std::cout) const {
        fout<<"id="<<id<<std::endl;
        fout<<"rank_org="<<rank_org<<std::endl;
        fout<<"pos="<<pos<<std::endl;
        fout<<"eps="<<eps<<std::endl;
    }
    PS::F64 getRSearch() const {
//        return std::max(r_search * std::pow(2.0*mass/m_average,0.3333),r_search_min) * SAFTY_FACTOR_FOR_SEARCH + SAFTY_OFFSET_FOR_SEARCH;
        return r_search * SAFTY_FACTOR_FOR_SEARCH + SAFTY_OFFSET_FOR_SEARCH;
    }

};


class EPJSoft{
public:
    PS::S64 id;
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64vec vel;
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
        return r_search * SAFTY_FACTOR_FOR_SEARCH + SAFTY_OFFSET_FOR_SEARCH;
    }
    // FORDEBUG
    void dump(std::ostream & fout=std::cout) const {
        fout<<"id="<<id<<std::endl;
        fout<<"rank_org="<<rank_org<<std::endl;
        fout<<"mass="<<mass<<std::endl;
        fout<<"pos="<<pos<<std::endl;
        fout<<"vel="<<vel<<std::endl;
        fout<<"r_search="<<r_search<<std::endl;
    }
    void clear(){
        mass = 0.0;
        pos = vel = 0.0;
        r_search = 0.0;
        id = rank_org = adr_org = -1;
    }
};

PS::F64 EPISoft::eps = 1.0/1024.0;
PS::F64 EPISoft::r_out = 0.0;
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
    void dump(std::ostream & fout=std::cout){
        fout<<"Energy: total = "<<tot<<" kin = "<<kin<<" pot = "<<pot
            <<"\nAngular Momentum: L = "<<L<<" |L| = "<<Lt
            <<std::endl;
    }
    void  writeAscii(FILE* fout) {
        fprintf(fout, "%26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e ",
                this->kin, this->pot, this->tot, this->L[0], this->L[1], this->L[2], this->Lt);
    }
    void writeBinary(FILE* fout) {
        fwrite(&this->kin, sizeof(EnergyAndMomemtum), 1, fout);
    }

    template<class Tptcl>
    void calc(const Tptcl* sys,
              const PS::S32 n,
              const PS::F64 dt_soft){
        // PS::S32 n = sys.getNumberOfParticleLocal();
        PS::F64 pot_loc = 0.0;
        PS::F64 kin_loc = 0.0;
        PS::F64vec L_loc = PS::F64vec(0.0);
        //#pragma omp parallel for reduction(+:pot_d) reduction (+:pot_loc) reduction(+:kin_loc)
        PS::ReallocatableArray<PS::S32> plist;
        plist.reserve(n);
        for(PS::S32 i=0; i<n; i++){
            if(sys[i].id<0||sys[i].status>0) continue;
            plist.push_back(i);
        }
        for(PS::S32 ki=0; ki<plist.size(); ki++){
            PS::S32 i = plist[ki];
            PS::F64 mi = sys[i].mass;
            PS::F64vec vi = sys[i].vel;
            if(sys[i].status!=0) {
                mi = sys[i].mass_bk;
                vi += sys[i].acc * dt_soft;
            }
            pot_loc += 0.5 * mi * sys[i].pot_tot;
            //for (PS::S32 kj=0; kj<ki; kj++)  {
            //    PS::S32 j = plist[kj];
            //    PS::F64 mj = sys[j].mass;
            //    if(sys[j].status!=0) mj = sys[j].mass_bk;
            //    PS::F64vec dr = sys[i].pos-sys[j].pos;
            //    PS::F64 dr2 = dr*dr;
            //    PS::F64 drm = 1.0/sqrt(dr2 + EPISoft::eps*EPISoft::eps);
            //    pot_loc += - mi * mj * drm;
            //}
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
};


////////////////////
/// FORCE FUNCTOR
struct CalcForceEpEpWithLinearCutoffNoSIMD{
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
                      const EPJSoft * ep_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
        const PS::F64 r_out = EPISoft::r_out;
        //        const PS::F64 r_crit2 = EPJSoft::r_search * EPJSoft::r_search * SAFTY_FACTOR_FOR_SEARCH_SQ;
        // const PS::F64 r_out = EPISoft::r_out; 
        // const PS::F64 r_in = EPISoft::r_in;
        //std::cerr<<"r_out= "<<r_out<<" r_in= "<<r_in<<" eps2= "<<eps2<<" r_crit2= "<<r_crit2<<std::endl;
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F64vec xi = ep_i[i].pos;
            PS::S64 id_i = ep_i[i].id;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            PS::S32 n_ngb_i = 0;
            for(PS::S32 j=0; j<n_jp; j++){
                const PS::F64 r_search = std::max(ep_i[i].r_search,ep_j[j].r_search);
                if(id_i == ep_j[j].id){
                    n_ngb_i++;
                    continue;
                }
                const PS::F64vec rij = xi - ep_j[j].pos;
                const PS::F64 r2 = rij * rij;
                const PS::F64 r2_eps = rij * rij + eps2;
                if(r2 < r_search*r_search*SAFTY_FACTOR_FOR_SEARCH_SQ){
                    n_ngb_i++;
                }
                const PS::F64 r2_tmp = (r2_eps > r_out*r_out) ? r2_eps : r_out*r_out;
                const PS::F64 r_inv = 1.0/sqrt(r2_tmp);
                const PS::F64 m_r = ep_j[j].mass * r_inv;
                const PS::F64 m_r3 = m_r * r_inv * r_inv;
                ai -= m_r3 * rij;
                poti -= m_r;
            }
            //std::cerr<<"poti= "<<poti<<std::endl;
            force[i].acc += ai;
            force[i].pot += poti;
            force[i].n_ngb = n_ngb_i;
        }
    }
};

struct CalcForceEpSpNoSIMD{
    template<class Tsp>
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
                      const Tsp * sp_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
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
