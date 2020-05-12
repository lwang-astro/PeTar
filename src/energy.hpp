#pragma once

//! class for collecting and calculating the energy and angular momemtum of the system
class EnergyAndMomemtum{
public:
    PS::F64 error_cum_pre;  // previous cumulative error
    PS::F64 ekin;
    PS::F64 epot;
    PS::F64 etot_ref; // energy reference (initial value) for calculating energy error
#ifdef HARD_CHECK_ENERGY
    PS::F64 de_change_cum; // cumulative energy change 
    PS::F64 de_change_binary_interrupt; // cumulative energy change due to interruption (only hard, without cluster cm energy change)
    PS::F64 de_change_modify_single; // cumulative energy change due to modification of particle (only hard, without cluster cm energy change)
    PS::F64 error_hard_cum_pre;
    PS::F64 error_hard_cum;
    // slowdown energy
    PS::F64 error_sd_cum_pre;  // previous cumulative slowdown error
    PS::F64 ekin_sd;
    PS::F64 epot_sd;
    PS::F64 etot_sd_ref; // slowdown energy reference (initial value) for calculating energy error
    PS::F64 de_sd_change_cum; // cumulative slowdown energy change 
    PS::F64 de_sd_change_binary_interrupt; // cumulative slowdown energy change du to interruption (only hard, without cluster cm energy change)
    PS::F64 de_sd_change_modify_single; // cumulative slowdown energy change due to modification of particle (only hard, without cluster cm energy change)
    PS::F64 error_hard_sd_cum_pre;
    PS::F64 error_hard_sd_cum;
#endif
    PS::F64 error_Lt_cum_pre; // previous angular momentum error
    PS::F64 Lt; // total angular momemtum
    PS::F64vec L; // angular mommentum;
    PS::F64vec L_ref; // total angular momemtum reference

    EnergyAndMomemtum() {
        clear();
    }

    void clear(){
        error_cum_pre = ekin = epot = etot_ref = 0.0;
#ifdef HARD_CHECK_ENERGY
        de_change_cum = de_change_binary_interrupt = de_change_modify_single = 0.0;
        error_hard_cum = error_hard_cum_pre = 0.0;
        error_sd_cum_pre = ekin_sd = epot_sd = etot_sd_ref = de_sd_change_cum = de_sd_change_binary_interrupt = de_sd_change_modify_single = 0.0;
        error_hard_sd_cum = error_hard_sd_cum_pre = 0.0;
#endif
        error_Lt_cum_pre = Lt = 0.0;
        L = L_ref = PS::F64vec(0.0);
    }

    //! print title and values in one lines
    /*! print titles and values in one lines
      @param[out] _fout: std::ostream output object
    */
    void print(std::ostream & _fout=std::cout, const PS::S32 _width=20) {
        _fout<<"Energy:  "
             <<std::setw(_width)<<"Error/Total"
             <<std::setw(_width)<<"Error"
             <<std::setw(_width)<<"Error_cum"
             <<std::setw(_width)<<"Total"
             <<std::setw(_width)<<"Kinetic"
             <<std::setw(_width)<<"Potential"
#ifdef HARD_CHECK_ENERGY
             <<std::setw(_width)<<"Modify"
             <<std::setw(_width)<<"Modify_group"
             <<std::setw(_width)<<"Modify_single"
             <<std::setw(_width)<<"Error_hard"
             <<std::setw(_width)<<"Error_hard_cum"
#endif
             <<std::endl;
        _fout<<"Physic:  "
             <<std::setw(_width)<<(getEnergyError() - error_cum_pre)/(ekin+epot)
             <<std::setw(_width)<<getEnergyError() - error_cum_pre
             <<std::setw(_width)<<getEnergyError()
             <<std::setw(_width)<<ekin + epot  
             <<std::setw(_width)<<ekin         
             <<std::setw(_width)<<epot         
#ifdef HARD_CHECK_ENERGY
             <<std::setw(_width)<<de_change_cum
             <<std::setw(_width)<<de_change_binary_interrupt
             <<std::setw(_width)<<de_change_modify_single
             <<std::setw(_width)<<error_hard_cum - error_hard_cum_pre
             <<std::setw(_width)<<error_hard_cum
#endif
             <<std::endl;
#ifdef HARD_CHECK_ENERGY
        _fout<<"Slowdown:"
             <<std::setw(_width)<<(getEnergyErrorSlowDown() - error_sd_cum_pre)/(ekin_sd+epot_sd)
             <<std::setw(_width)<<getEnergyErrorSlowDown() - error_sd_cum_pre
             <<std::setw(_width)<<getEnergyErrorSlowDown()
             <<std::setw(_width)<<ekin_sd + epot_sd
             <<std::setw(_width)<<ekin_sd
             <<std::setw(_width)<<epot_sd
             <<std::setw(_width)<<de_sd_change_cum
             <<std::setw(_width)<<de_sd_change_binary_interrupt
             <<std::setw(_width)<<de_sd_change_modify_single
             <<std::setw(_width)<<error_hard_sd_cum - error_hard_sd_cum_pre
             <<std::setw(_width)<<error_hard_sd_cum
             <<std::endl;
#endif
        _fout<<"Angular Momemtum:"
             <<"  |L|err: "<<getMomentumError() - error_Lt_cum_pre
             <<"  |L|err_cum: "<<getMomentumError()
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
        _fout<<std::setw(_width)<<"Error"
             <<std::setw(_width)<<"Error_cum"
             <<std::setw(_width)<<"Ekin"
             <<std::setw(_width)<<"Epot"
             <<std::setw(_width)<<"Etot"
#ifdef HARD_CHECK_ENERGY
             <<std::setw(_width)<<"dE_modify"
             <<std::setw(_width)<<"dE_interrupt"
             <<std::setw(_width)<<"Error_hard"
             <<std::setw(_width)<<"Error_hard_cum"
             <<std::setw(_width)<<"Error_SD"
             <<std::setw(_width)<<"Error_SD_cum"
             <<std::setw(_width)<<"Ekin_SD"
             <<std::setw(_width)<<"Epot_SD"
             <<std::setw(_width)<<"Etot_SD"
             <<std::setw(_width)<<"dE_modify_SD"
             <<std::setw(_width)<<"dE_interrupt_SD"
             <<std::setw(_width)<<"Error_hard_SD"
             <<std::setw(_width)<<"Error_hard_SD_cum"
#endif
             <<std::setw(_width)<<"|L|error"
             <<std::setw(_width)<<"|L|error_cum"
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
        _fout<<std::setw(_width)<<getEnergyError() - error_cum_pre
             <<std::setw(_width)<<getEnergyError()
             <<std::setw(_width)<<ekin
             <<std::setw(_width)<<epot
             <<std::setw(_width)<<ekin+epot
#ifdef HARD_CHECK_ENERGY
             <<std::setw(_width)<<de_change_cum
             <<std::setw(_width)<<de_change_binary_interrupt
             <<std::setw(_width)<<error_hard_cum - error_hard_cum_pre
             <<std::setw(_width)<<error_hard_cum
             <<std::setw(_width)<<getEnergyErrorSlowDown() - error_sd_cum_pre
             <<std::setw(_width)<<getEnergyErrorSlowDown()
             <<std::setw(_width)<<ekin_sd
             <<std::setw(_width)<<epot_sd
             <<std::setw(_width)<<ekin_sd+epot_sd
             <<std::setw(_width)<<de_sd_change_cum
             <<std::setw(_width)<<de_sd_change_binary_interrupt
             <<std::setw(_width)<<error_hard_sd_cum - error_hard_sd_cum_pre
             <<std::setw(_width)<<error_hard_sd_cum
#endif
             <<std::setw(_width)<<getMomentumError() - error_Lt_cum_pre
             <<std::setw(_width)<<getMomentumError()
             <<std::setw(_width)<<L.x
             <<std::setw(_width)<<L.y
             <<std::setw(_width)<<L.z
             <<std::setw(_width)<<Lt;
    }

#ifdef HARD_CHECK_ENERGY
    void writeAscii(FILE* _fout) {
        fprintf(_fout, "%26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e ",
                error_cum_pre, ekin, epot, etot_ref, de_change_cum, de_change_binary_interrupt,
                error_sd_cum_pre, ekin_sd, epot_sd, etot_sd_ref, de_sd_change_cum, de_sd_change_binary_interrupt,
                error_Lt_cum_pre, Lt, L[0], L[1], L[2], 
                L_ref[0], L_ref[1], L_ref[2]);
    }
#else
    void writeAscii(FILE* _fout) {
        fprintf(_fout, "%26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e ",
                error_cum_pre, ekin, epot, etot_ref,
                error_Lt_cum_pre, Lt, L[0], L[1], L[2], 
                L_ref[0], L_ref[1], L_ref[2]);
    }
#endif

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
        assert(Ptcl::group_data_mode == GroupDataMode::artificial);
        ekin = epot = 0.0;
        L = PS::F64vec(0.0);
//#pragma omp declare reduction(+:PS::F64vec:omp_out += omp_in) initializer (omp_priv=PS::F64vec(0.0))
//#pragma omp parallel for reduction(+:epot,ekin,L) 
        for(PS::S32 i=0; i<_n_particle; i++){
            PS::F64 mi = _particles[i].mass;
            auto pi_artificial = _particles[i].group_data.artificial;
            if(pi_artificial.isMember()) mi = pi_artificial.getMassBackup();
#ifdef HARD_DEBUG
            assert(_particles[i].id>0&&(pi_artificial.isMember()||pi_artificial.isSingle()));
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
#ifdef HARD_CHECK_ENERGY
            etot_sd_ref = etot_ref;
#endif
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
//#pragma omp declare reduction(+:PS::F64vec:omp_out += omp_in) initializer (omp_priv=PS::F64vec(0.0))
//#pragma omp parallel for reduction(+:epot,ekin,L) 
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
#ifdef HARD_CHECK_ENERGY
            etot_sd_ref = etot_ref;
#endif
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
#ifdef HARD_CHECK_ENERGY
            etot_sd_ref = etot_ref;
#endif
            L_ref = L;
        }
    }

    //! save current energy error
    void saveEnergyError() {
        error_cum_pre = getEnergyError();
#ifdef HARD_CHECK_ENERGY
        error_sd_cum_pre = getEnergyErrorSlowDown();
        error_hard_cum_pre  = error_hard_cum;
        error_hard_sd_cum_pre = error_hard_sd_cum;
#endif
        error_Lt_cum_pre  = getMomentumError();
    }

    //! get energy error
    PS::F64 getEnergyError() const {
        return ekin + epot - etot_ref;
    }

#ifdef HARD_CHECK_ENERGY
    //! get slowdown energy error
    PS::F64 getEnergyErrorSlowDown() const {
        return ekin_sd + epot_sd - etot_sd_ref;
    }
#endif

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

