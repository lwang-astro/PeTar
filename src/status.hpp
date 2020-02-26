#pragma once

#include "usr_define.hpp"
#include "energy.hpp"

//! class for measure the status of the system
class Status {
public:
    PS::F64 time;
    PS::S64 n_real_loc;
    PS::S64 n_real_glb;
    PS::S64 n_all_loc;
    PS::S64 n_all_glb;
    PS::F64 half_mass_radius;
    ParticleBase pcm;
    EnergyAndMomemtum energy;
    PS::F64 energy_hard_diff;
    PS::F64 energy_hard_sd_diff;

    Status(): time(0.0), n_real_loc(0), n_real_glb(0), n_all_loc(0), n_all_glb(0), half_mass_radius(0), pcm() {}

    template <class Tsoft>
    void calcCenterOfMass(Tsoft* _tsys, const PS::S64 _n) {
        PS::F64 mass = 0.0;
        PS::F64vec pos_cm = PS::F64vec(0.0);
        PS::F64vec vel_cm = PS::F64vec(0.0);

#pragma omp declare reduction(+:PS::F64vec:omp_out += omp_in) initializer (omp_priv=PS::F64vec(0.0))
#pragma omp parallel for reduction(+:mass,pos_cm,vel_cm)
        for (int i=0; i<_n; i++) {
            auto& pi = _tsys[i];
            PS::F64 mi = pi.mass;
            mass  += mi;
#ifdef NAN_CHECK_DEBUG
            assert(!std::isnan(pi.vel.x));
            assert(!std::isnan(pi.vel.y));
            assert(!std::isnan(pi.vel.z));
#endif
            pos_cm += mi*pi.pos;
            vel_cm += mi*pi.vel;
        }        

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        pcm.mass = PS::Comm::getSum(mass);
        pcm.pos  = PS::Comm::getSum(pos_cm);
        pcm.vel  = PS::Comm::getSum(vel_cm);
#else
        pcm.mass = mass;
        pcm.pos  = pos_cm;
        pcm.vel  = vel_cm;
#endif
        pcm.pos /= pcm.mass;
        pcm.vel /= pcm.mass;
    }

    //! print titles of class members using column style
    /*! print titles of class members in one line for column style
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
    */
    void printColumnTitle(std::ofstream & _fout, const PS::S32 _width=20) {
        _fout<<std::setw(_width)<<"Time"
             <<std::setw(_width)<<"N_real_loc"
             <<std::setw(_width)<<"N_real_glb"
             <<std::setw(_width)<<"N_all_loc"
             <<std::setw(_width)<<"N_all_glb";
        energy.printColumnTitle(_fout, _width);
#ifdef HARD_CHECK_ENERGY
        _fout<<std::setw(_width)<<"dE_hard"
             <<std::setw(_width)<<"dE_SD_hard";
#endif
        _fout<<std::setw(_width)<<"CM.mass"
             <<std::setw(_width)<<"CM.pos.x"
             <<std::setw(_width)<<"CM.pos.y"
             <<std::setw(_width)<<"CM.pos.z"
             <<std::setw(_width)<<"CM.vel.x"
             <<std::setw(_width)<<"CM.vel.y"
             <<std::setw(_width)<<"CM.vel.z";
    }

    //! print data of class members using column style
    /*! print data of class members in one line for column style. Notice no newline is printed at the end
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
    */
    void printColumn(std::ofstream & _fout, const PS::S32 _width=20) {
        _fout<<std::setw(_width)<<time
             <<std::setw(_width)<<n_real_loc
             <<std::setw(_width)<<n_real_glb
             <<std::setw(_width)<<n_all_loc
             <<std::setw(_width)<<n_all_glb;
        energy.printColumn(_fout, _width);
#ifdef HARD_CHECK_ENERGY
        _fout<<std::setw(_width)<<energy_hard_diff
             <<std::setw(_width)<<energy_hard_sd_diff;
#endif
        _fout<<std::setw(_width)<<pcm.mass
             <<std::setw(_width)<<pcm.pos.x
             <<std::setw(_width)<<pcm.pos.y
             <<std::setw(_width)<<pcm.pos.z
             <<std::setw(_width)<<pcm.vel.x
             <<std::setw(_width)<<pcm.vel.y
             <<std::setw(_width)<<pcm.vel.z;
    }

    //! print title and values in one lines
    /*! print titles and values in one lines
      @param[out] _fout: std::ostream output object
      @param[in] _precision: floating data precision except time (15 digital)
    */
    void print(std::ostream & _fout, const PS::S32 _precision=7) {
        _fout<<"Time: "<<std::setprecision(15)<<time
             <<std::setprecision(_precision);
        _fout<<"  N_real(loc): "<<n_real_loc
             <<"  N_real(glb): "<<n_real_glb
             <<"  N_all(loc): "<<n_all_loc
             <<"  N_all(glb): "<<n_all_glb;
#ifdef HARD_CHECK_ENERGY
        _fout<<"  Hard Energy: dE: "<<energy_hard_diff
             <<"  dE_SD: "<<energy_hard_sd_diff
             <<std::endl;
#endif
        _fout<<"C.M.: mass: "<<pcm.mass
             <<" pos: "<<pcm.pos
             <<" vel: "<<pcm.vel
             <<std::endl;
        energy.print(_fout);
    }
};
