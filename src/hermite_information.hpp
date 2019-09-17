#pragma once

#include "Common/Float.h"
#include "hermite_interaction.hpp"

class HermiteInformation{
public:
    Float time_origin; // time of origin
    Float time; // current time in integrator
    Float de;    // energy difference
    Float etot_init; // initial total energy;
    Float etot; // total energy
    Float ekin; // kinetic energy
    Float epot; // potential energy
    Float ett;  // tidal tensor energy
    Float de_sd; // slowdown energy difference
    Float etot_sd_ref; // reference slowdown energy for calculating de
    Float etot_sd; // slowdown energy record
    Float ekin_sd; // slowdown kinetic energy record
    Float epot_sd; // slowdown potential energy record

    //! check whether parameters values are correct
    /*! \return true: all correct
     */
    bool checkParams() {
        return true;
    }        

    //! calculate energy of particle group
    template <class Tptcl, class Tgroup>
    void calcEnergySlowDown(Tptcl& _particles, Tgroup& _groups, HermiteInteraction& _interaction, const bool _initial_flag) {
        ekin = epot = etot = ett = 0.0;
        const int n = _particles.getSize();
        // inner
        for (int i=0; i<n; i++) {
            auto& pi = _particles[i];
            ekin += pi.mass* (pi.vel[0]*pi.vel[0] + pi.vel[1]*pi.vel[1] + pi.vel[2]*pi.vel[2]);
            Float poti = 0.0;
            for (int j=0; j<i; j++) {
                poti +=_interaction.calcPotPair(pi, _particles[j]);
            }
            epot += poti*pi.mass;
        }
#ifdef SOFT_PERT
        // tidal
        const int n_group =_groups.getSize();
        for (int i=0; i<n_group; i++) {
            const int n_member = _groups[i].particles.getSize();
            if (_groups[i].perturber.soft_pert!=NULL) {
                auto* pert = _groups[i].perturber.soft_pert;
                for (int j=0; j<n_member; j++) {
                    ett += _groups[i].particles[j].mass*pert->evalPot(_groups[i].particles[j].pos);
                }
            }
        }
#endif
        ekin *= 0.5;
        ett *= 0.5;
        etot = ekin + epot + ett;

        // slowdown energy
        ekin_sd = ekin;
        epot_sd = epot;
        for (int i=0; i<_groups.getSize(); i++) {
            auto& gi = _groups[i];
            Float kappa = gi.slowdown.getSlowDownFactor();
            Float kappa_inv = 1.0/kappa;
            ekin_sd -= gi.getEkin();
            epot_sd -= gi.getEpot(); 
#ifdef AR_TTL_SLOWDOWN_INNER
            ekin_sd += kappa_inv * gi.getEkinSlowDownInner();
            epot_sd += kappa_inv * gi.getEpotSlowDownInner();
#else
            ekin_sd += kappa_inv * gi.getEkin();
            epot_sd += kappa_inv * gi.getEpot(); 
#endif
        }
        etot_sd = ekin_sd + epot_sd;

        // initial 
        if (_initial_flag) {
            etot_init = etot;
            etot_sd_ref = etot_sd;
            de = 0.0;
            de_sd = 0.0;
        }
        else {
            de = etot - etot_init;
            de_sd = etot_sd - etot_sd_ref;
        }
    }

    //! correct Etot slowdown reference due to the groups change
    template <class TGroupList>
    void correctEtotSlowDownRef(TGroupList& _groups) {
        for (int i=0; i<_groups.getSize(); i++) {
            auto& gi = _groups[i];
            etot_sd_ref += gi.getDESlowDownCum();
            gi.resetDESlowDownCum();
        } 
    }

    //! print titles of class members using column style
    /*! print titles of class members in one line for column style
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
    */
    void printColumnTitle(std::ostream & _fout, const int _width=20) {
        _fout<<std::setw(_width)<<"Time_org"
             <<std::setw(_width)<<"Time_int"
             <<std::setw(_width)<<"dE"
             <<std::setw(_width)<<"Etot"
             <<std::setw(_width)<<"Ekin"
             <<std::setw(_width)<<"Epot"
             <<std::setw(_width)<<"Ett"
             <<std::setw(_width)<<"dE_SD"
             <<std::setw(_width)<<"Etot_SD_ref"
             <<std::setw(_width)<<"Etot_SD"
             <<std::setw(_width)<<"Ekin_SD"
             <<std::setw(_width)<<"Epot_SD";
    }

    //! print data of class members using column style
    /*! print data of class members in one line for column style. Notice no newline is printed at the end
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
    */
    void printColumn(std::ostream & _fout, const int _width=20){
        _fout<<std::setw(_width)<<time_origin
             <<std::setw(_width)<<time
             <<std::setw(_width)<<de
             <<std::setw(_width)<<etot
             <<std::setw(_width)<<ekin
             <<std::setw(_width)<<epot
             <<std::setw(_width)<<ett
             <<std::setw(_width)<<de_sd
             <<std::setw(_width)<<etot_sd_ref
             <<std::setw(_width)<<etot_sd
             <<std::setw(_width)<<ekin_sd
             <<std::setw(_width)<<epot_sd;
    }
};
