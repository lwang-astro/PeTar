#pragma once
#include"ptcl.hpp"

class ForceSoft{
public:
    PS::F64vec acc; ///> soft acceleration (c.m.: averaged force from orbital particles; tensor: c.m. is substracted)
    PS::F64 pot; ///> full potential
#ifdef KDKDK_4TH
    PS::F64vec acorr; ///> soft gradient correction for 4th order KDKDK method
#endif
#ifdef SAVE_NEIGHBOR_ID_IN_FORCE_KERNEL
    PS::S64 id_ngb[4]; /// five neighbor id
#endif
    PS::S64 n_ngb; ///> neighbor number+1
    static PS::F64 grav_const; ///> gravitational constant
    void clear(){
        acc = 0.0;
#ifdef KDKDK_4TH
        acorr = 0.0;
#endif        
        pot = 0.0;
        n_ngb = 0;
#ifdef SAVE_NEIGHBOR_ID_IN_FORCE_KERNEL
        id_ngb[0] = id_ngb[1] = id_ngb[2] = id_ngb[3] = 0;
#endif
    }
};

class FPSoft: public Ptcl{
public:
    PS::F64vec acc; // soft
#ifdef KDKDK_4TH
    PS::F64vec acorr;
#endif
    PS::F64 pot_tot; // soft + hard
    PS::F64 pot_soft; // soft only
#ifdef EXTERNAL_POT_IN_PTCL
    PS::F64 pot_ext; // external potential
#endif
#ifdef SAVE_NEIGHBOR_ID_IN_FORCE_KERNEL
    PS::S64 id_ngb[4];
#endif
    PS::S64 n_ngb;
    PS::S32 rank_org;
    PS::S32 adr;
//    static PS::F64 r_out;

    FPSoft() {}

    
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
        pot_soft= 0;
#ifdef EXTERNAL_POT_IN_PTCL
        pot_ext = 0;
#endif
        n_ngb = 0;
    }

    void copyFromForce(const ForceSoft & force){
        acc = force.acc;
        pot_tot = force.pot;
        pot_soft= pot_tot;
#ifdef KDKDK_4TH
        acorr = force.acorr;
#endif
#ifdef SAVE_NEIGHBOR_ID_IN_FORCE_KERNEL
        for (int k=0; k<4; k++) id_ngb[k] = force.id_ngb[k];
#endif
        n_ngb = force.n_ngb;
    }

    PS::F64 getRSearch() const {
//        return std::max(r_search * std::pow(2.0*mass/m_average,0.3333),r_search_min) * SAFTY_FACTOR_FOR_SEARCH + SAFTY_OFFSET_FOR_SEARCH;
//        return this->r_search * SAFTY_FACTOR_FOR_SEARCH + SAFTY_OFFSET_FOR_SEARCH;
        return this->r_search*SAFTY_FACTOR_FOR_SEARCH;
    }

    void writeAscii(FILE* fp) const{
        Ptcl::writeAscii(fp);
        fprintf(fp, "%26.17e %26.17e %26.17e %26.17e %26.17e ",
                this->acc.x, this->acc.y, this->acc.z,  // 9-11
                this->pot_tot, this->pot_soft);
#ifdef EXTERNAL_POT_IN_PTCL
        fprintf(fp, "%26.17e ",this->pot_ext);
#endif        
        fprintf(fp, "%lld\n",this->n_ngb);
    }

    void writeBinary(FILE* fp) const{
        Ptcl::writeBinary(fp);
#ifdef EXTERNAL_POT_IN_PTCL
        fwrite(&(this->acc), sizeof(PS::F64), 7, fp);
#else
        fwrite(&(this->acc), sizeof(PS::F64), 6, fp);
#endif
    }

    void readAscii(FILE* fp) {
        Ptcl::readAscii(fp);
        PS::S64 rcount=fscanf(fp, "%lf %lf %lf %lf %lf ",
                              &this->acc.x, &this->acc.y, &this->acc.z,  // 9-11
                              &this->pot_tot, &this->pot_soft);
        if (rcount<5) {
            std::cerr<<"Error: FPSoft Data reading fails! requiring data number is 6, only obtain "<<rcount<<".\n";
            std::cerr<<"Check your input data, whether the consistent features (interrupt mode and external mode) are used in configuring petar and the data generation\n";
            abort();
        }
#ifdef EXTERNAL_POT_IN_PTCL
        rcount=fscanf(fp, "%lf ", &this->pot_ext);
        if (rcount<1) {
            std::cerr<<"Error: FPSoft Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            std::cerr<<"Check your input data, whether the consistent features (interrupt mode and external mode) are used in configuring petar and the data generation\n";
            abort();
        }
#endif
        rcount=fscanf(fp, "%lld\n", &this->n_ngb);
        if (rcount<1) {
            std::cerr<<"Error: FPSoft Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            std::cerr<<"Check your input data, whether the consistent features (interrupt mode and external mode) are used in configuring petar and the data generation\n";
            abort();
        }
    }

    void readBinary(FILE* fp) {
        Ptcl::readBinary(fp);
#ifdef EXTERNAL_POT_IN_PTCL
        size_t rcount = fread(&(this->acc), sizeof(PS::F64), 7, fp);
        if (rcount<7) {
            std::cerr<<"Error: Data reading fails! requiring data number is 7, only obtain "<<rcount<<".\n";
            std::cerr<<"Check your input data, whether the consistent features (interrupt mode and external mode) are used in configuring petar and the data generation\n";
            abort();
        }
#else
        size_t rcount = fread(&(this->acc), sizeof(PS::F64), 6, fp);
        if (rcount<6) {
            std::cerr<<"Error: Data reading fails! requiring data number is 6, only obtain "<<rcount<<".\n";
            std::cerr<<"Check your input data, whether the consistent features (interrupt mode and external mode) are used in configuring petar and the data generation\n";
            abort();
        }
#endif
    }

    void print(std::ostream & fout){
        Ptcl::print(fout);
        fout<<" acc= "<<acc
            <<" pot_tot= "<<pot_tot
            <<" pot_soft= "<<pot_soft
#ifdef EXTERNAL_POT_IN_PTCL
            <<" pot_ext= "<<pot_ext
#endif
            <<" N_b= "<<n_ngb;
    }

    //! print titles of class members using column style
    /*! print titles of class members in one line for column style
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
     */
    static void printColumnTitle(std::ostream & _fout, const int _width=20) {
        Ptcl::printColumnTitle(_fout, _width);
        _fout<<std::setw(_width)<<"acc_soft.x"
             <<std::setw(_width)<<"acc_soft.y"
             <<std::setw(_width)<<"acc_soft.z"
             <<std::setw(_width)<<"pot_tot"
             <<std::setw(_width)<<"pot_soft"
#ifdef EXTERNAL_POT_IN_PTCL
             <<std::setw(_width)<<"pot_ext"
#endif
             <<std::setw(_width)<<"n_b";
    }

    //! print column title with meaning (each line for one column)
    /*! @param[out] _fout: std::ostream output object
      @param[in] _counter: offset of the number counter for each line to indicate the column index (defaulted 0)
      @param[in] _offset: the printing whitespace offset for each line (defaulted 0)
      \return: the total counter of columns
     */
    static int printTitleWithMeaning(std::ostream & _fout, const int _counter=0, const int _offset=0) {
        int counter = _counter;
        counter = Ptcl::printTitleWithMeaning(_fout, counter, _offset);
        counter++;
        _fout<<std::setw(_offset)<<" "<<counter<<"-"<<counter+2<<". acc_soft.[x/y/z]: 3D soft (long-range) acceleration (0.0)\n";
        counter+=3;
        _fout<<std::setw(_offset)<<" "<<counter<<". pot_tot: total potential (0.0)\n";
        counter++;
        _fout<<std::setw(_offset)<<" "<<counter<<". pot_soft: soft potential (0.0)\n";
        counter++;
#ifdef EXTERNAL_POT_IN_PTCL
        _fout<<std::setw(_offset)<<" "<<counter<<". pot_ext: external potential (0.0)\n";
        counter++;
#endif
        _fout<<std::setw(_offset)<<" "<<counter<<". n_b: number of neighbors (0)\n";
        return counter;
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
             <<std::setw(_width)<<pot_soft
#ifdef EXTERNAL_POT_IN_PTCL
             <<std::setw(_width)<<pot_ext
#endif
             <<std::setw(_width)<<n_ngb;
    }

    //! print data of class members with pos and vel offset using column style
    /*! print data of class members in one line for column style. Notice no newline is printed at the end
      @param[in] _pcm: particle data with position and velocity offset that are added when print data
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
     */
    template <class Tpcm>
    void printColumnWithOffset(Tpcm& _pcm, std::ostream & _fout, const int _width=20){
        Ptcl::printColumnWithOffset(_pcm, _fout, _width);
        _fout<<std::setw(_width)<<acc.x
             <<std::setw(_width)<<acc.y
             <<std::setw(_width)<<acc.z
             <<std::setw(_width)<<pot_tot
             <<std::setw(_width)<<pot_soft
#ifdef EXTERNAL_POT_IN_PTCL
             <<std::setw(_width)<<pot_ext
#endif
             <<std::setw(_width)<<n_ngb;
    }

    //! clear force
    void clearForce() {
        acc = 0.0;
#ifdef KDKDK_4TH
        acorr = 0.0;
#endif
        pot_tot = 0.0;
        pot_soft = 0.0;
#ifdef EXTERNAL_POT_IN_PTCL
        pot_ext = 0.0;
#endif        
        n_ngb = 0;
    }
};

class EPISoft{
public:
    PS::S64 id;
    PS::F64vec pos;
    PS::F64 r_search;
    PS::S32 rank_org;
    PS::S32 type; // 0: orbital artificial particles; 1: others
#ifdef KDKDK_4TH
    PS::F64vec acc;
#endif
    static PS::F64 eps;
    static PS::F64 r_out;
    PS::F64vec getPos() const { return pos;}
    void copyFromFP(const FPSoft & fp){ 
        id = fp.id;
        pos = fp.pos;
        if (fp.group_data.artificial.isArtificial() && !fp.group_data.artificial.isCM() && fp.mass>0 ) type = 0;
        else type = 1;

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
    GroupDataDeliver group_data;
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
        group_data = fp.group_data;
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

    PS::S64 getId() const {
        return id;
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

