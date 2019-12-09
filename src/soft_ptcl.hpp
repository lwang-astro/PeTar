#pragma once
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
    PS::F64vec acc; // soft
#ifdef KDKDK_4TH
    PS::F64vec acorr;
#endif
    PS::F64 pot_tot; // soft + hard
    PS::S32 n_ngb;
    PS::S32 rank_org;
    PS::S32 adr;
//    static PS::F64 r_out;

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
        fout<<" acc= "<<acc
            <<" pot_tot= "<<pot_tot
            <<" N_b= "<<n_ngb;
    }

    //! print titles of class members using column style
    /*! print titles of class members in one line for column style
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
     */
    void printColumnTitle(std::ostream & _fout, const int _width=20) {
        Ptcl::printColumnTitle(_fout, _width);
        _fout<<std::setw(_width)<<"acc.x"
             <<std::setw(_width)<<"acc.y"
             <<std::setw(_width)<<"acc.z"
             <<std::setw(_width)<<"pot_tot"
             <<std::setw(_width)<<"n_b";
    }

    //! print data of class members using column style
    /*! print data of class members in one line for column style. Notice no newline is printed at the end
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
     */
    void printColumn(std::ostream & _fout, const int _width=20){
        Ptcl::printColumn(_fout, _width);
        _fout<<std::setw(_width)<<acc.x
             <<std::setw(_width)<<acc.y
             <<std::setw(_width)<<acc.z
             <<std::setw(_width)<<pot_tot
             <<std::setw(_width)<<n_ngb;
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
    PS::F64 r_in;
    PS::F64 r_out;
    PS::F64 r_search;
    PS::F64 r_scale_next;
    Gdat mass_bk;
    Gdat status;
    PS::S32 rank_org;
    PS::S32 adr_org;
//    static PS::F64 r_out;
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
        r_in = fp.changeover.getRin();
        r_out = fp.changeover.getRout();
        r_scale_next = fp.changeover.r_scale_next;
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
        r_in = r_out = 0.0;
        r_search = 0.0;
        r_scale_next = 1.0;
        id = rank_org = adr_org = -1;
    }
};

PS::F64 EPISoft::eps = 0.0;
PS::F64 EPISoft::r_out = 0.0;
PS::F64 EPISoft::r_in  = 0.0;

//PS::F64 EPISoft::r_in; = 0.0;
//PS::F64 EPJSoft::m_average;
//PS::F64 EPJSoft::r_search_min;
//PS::F64 EPISoft::m_average;
//PS::F64 EPISoft::r_search_min;

