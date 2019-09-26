#pragma once

//! class for collecting and calculating the energy and angular momemtum of the system
class EnergyAndMomemtum{
public:
    PS::F64 ekin;
    PS::F64 epot;
    PS::F64 etot_ref; // energy reference (initial value) for calculating energy error
    // slowdown energy
    PS::F64 ekin_sd;
    PS::F64 epot_sd;
    PS::F64 etot_sd_ref; // slowdown energy reference (initial value) for calculating energy error
    PS::F64vec L; // angular momentum
    PS::F64 Lt; // total angular momemtum
    PS::F64vec L_ref; // total angular momemtum reference

    EnergyAndMomemtum() {
        clear();
    }

    void clear(){
        ekin = epot = etot_ref = Lt = 0.0;
        ekin_sd = epot_sd = etot_sd_ref = 0.0;
        L = L_ref = PS::F64vec(0.0);
    }

    //! print title and values in one lines
    /*! print titles and values in one lines
      @param[out] _fout: std::ostream output object
    */
    void print(std::ostream & _fout=std::cout) {
        _fout<<"Energy:"
             <<"  dE: "  <<getEnergyError()
             <<"  Etot: "<<ekin + epot  
             <<"  Ekin: "<<ekin         
             <<"  Epot: "<<epot         
             <<"  dE(SD): "  <<getEnergyErrorSlowDown()
             <<"  Etot(SD): "<<ekin_sd + epot_sd
             <<"  Ekin(SD): "<<ekin_sd
             <<"  Epot(SD): "<<epot_sd
             <<std::endl;
        _fout<<"Angular Momemtum:"
             <<"  dL: "<<getMomentumError()
             <<"  L: "<<L
             <<"  |L|: "<<Lt
             <<std::endl;
    }

    //! print titles of class members using column style
    /*! print titles of class members in one line for column style
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
    */
    void printColumnTitle(std::ofstream & _fout, const PS::S32 _width=20) const {
        _fout<<std::setw(_width)<<"dE"
             <<std::setw(_width)<<"Ekin"
             <<std::setw(_width)<<"Epot"
             <<std::setw(_width)<<"Etot"
             <<std::setw(_width)<<"dE_SD"
             <<std::setw(_width)<<"Ekin_SD"
             <<std::setw(_width)<<"Epot_SD"
             <<std::setw(_width)<<"Etot_SD"
             <<std::setw(_width)<<"d|L|"
             <<std::setw(_width)<<"Lx"
             <<std::setw(_width)<<"Ly"
             <<std::setw(_width)<<"Lz"
             <<std::setw(_width)<<"|L|";
    }

    //! print data of class members using column style
    /*! print data of class members in one line for column style. Notice no newline is printed at the end
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
    */
    void printColumn(std::ofstream & _fout, const PS::S32 _width=20) const {
        _fout<<std::setw(_width)<<getEnergyError()
             <<std::setw(_width)<<ekin
             <<std::setw(_width)<<epot
             <<std::setw(_width)<<ekin+epot
             <<std::setw(_width)<<getEnergyErrorSlowDown()
             <<std::setw(_width)<<ekin_sd
             <<std::setw(_width)<<epot_sd
             <<std::setw(_width)<<ekin_sd+epot_sd
             <<std::setw(_width)<<getMomentumError()
             <<std::setw(_width)<<L.x
             <<std::setw(_width)<<L.y
             <<std::setw(_width)<<L.z
             <<std::setw(_width)<<Lt;
    }

    void writeAscii(FILE* _fout) {
        fprintf(_fout, "%26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e ",
                ekin, epot, etot_ref,
                ekin_sd, epot_sd, etot_sd_ref, 
                L[0], L[1], L[2], Lt,
                L_ref[0], L_ref[1], L_ref[2]);
    }

    void writeBinary(FILE* _fout) {
        fwrite(&ekin, sizeof(EnergyAndMomemtum), 1, _fout);
    }

    //! calculate the system kinetic and potential energy of particles
    /*! 
      @param[in] _particles: particle array
      @param[in] _n_particle: number of particles
      @param[in] _init_flag: if true, set etot, etot_sd and L reference
     */
    template<class Tptcl>
    void calc(const Tptcl* _particles,
              const PS::S32 _n_particle, 
              const bool _init_flag=false) {
        ekin = epot = 0.0;
        L = PS::F64vec(0.0);
        for(PS::S32 i=0; i<_n_particle; i++){
            PS::F64 mi = _particles[i].mass;
            if(_particles[i].status.d<0) mi = _particles[i].mass_bk.d;
#ifdef HARD_DEBUG
            assert(_particles[i].id>0&&_particles[i].status.d<=0);
            assert(mi>0);
#endif
            PS::F64vec vi = _particles[i].vel;
            epot += 0.5 * mi * _particles[i].pot_tot;
            ekin += 0.5 * mi * vi * vi;
            L += _particles[i].pos ^ (mi*vi);
        }
        Lt = std::sqrt(L*L);
        if (_init_flag) {
            etot_ref = ekin + epot;
            etot_sd_ref = etot_ref;
            L_ref = L;
        }
    }

    //! calculate the system kinetic and potential energy of particles
    /*! Using particle index array to select particles
      @param[in] _particles: particle array
      @param[in] _particle_index: index array to select particles
      @param[in] _n_particle: number of particles
      @param[in] _init_flag: if true, set etot, etot_sd and L reference
     */
    template<class Tptcl> 
    void calc(const Tptcl* _particles,
              const PS::S32* _particle_index,
              const PS::S32 _n_particle, 
              const bool _init_flag=false) {
        ekin = epot = 0.0;
        L = PS::F64vec(0.0);
        for(PS::S32 k=0; k<_n_particle; k++){
            PS::S32 i = _particle_index[k];
            PS::F64 mi = _particles[i].mass;
            PS::F64vec vi = _particles[i].vel;
            epot += 0.5 * mi * _particles[i].pot_tot;
            ekin += 0.5 * mi * vi * vi;
            L += _particles[i].pos ^ (mi*vi);
        }
        Lt   = std::sqrt(L*L);
        if (_init_flag) {
            etot_ref = ekin + epot;
            etot_sd_ref = etot_ref;
            L_ref = L;
        }
    }

    //! get summation of kinetic, potential energy and angular momemtum of all MPI processes
    /*!
      @param[in] _init_flag: if true, set etot, etot_sd and L reference
     */
    void getSumMultiNodes(const bool _init_flag=false) {
        ekin = PS::Comm::getSum(ekin);
        epot = PS::Comm::getSum(epot);
        L   = PS::Comm::getSum(L);
        Lt  = std::sqrt(L*L);
        if (_init_flag) {
            etot_ref = ekin + epot;
            etot_sd_ref = etot_ref;
            L_ref = L;
        }
    }

    //! get energy error
    PS::F64 getEnergyError() const {
        return ekin + epot - etot_ref;
    }

    //! get slowdown energy error
    PS::F64 getEnergyErrorSlowDown() const {
        return ekin_sd + epot_sd - etot_sd_ref;
    }

    //! get angular momemtum (value) error
    PS::F64 getMomentumError() const {
        PS::F64vec dL= L- L_ref;
        return std::sqrt(dL*dL);
    }

    /*
    EnergyAndMomemtum operator -(const EnergyAndMomemtum& eng){
        EnergyAndMomemtum diff;
        diff.ekin = ekin - eng.ekin;
        diff.epot = epot - eng.epot;
        diff.etot = etot - eng.etot;
        diff.L   = L   - eng.L;
        diff.Lt  = std::sqrt(diff.L*diff.L);
        return diff;
    }

    void relative(const EnergyAndMomemtum& ref) {
        ekin /= ref.ekin;
        epot /= ref.epot;
        etot /= ref.etot;
        L   /= ref.Lt;
        Lt  /= ref.Lt;
    }
    */
};

