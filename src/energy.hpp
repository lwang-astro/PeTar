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

