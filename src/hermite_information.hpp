#pragma once

#include "Common/Float.h"
#include "hermite_interaction.hpp"

class HermiteInformation{
public:
    Float time; // current time
    Float etot0; // initial total energy;
    Float de;    // energy difference
    Float etot; // total energy
    Float ekin; // kinetic energy
    Float epot; // potential energy

    //! check whether parameters values are correct
    /*! \return true: all correct
     */
    bool checkParams() {
        return true;
    }        

    //! calculate energy of particle group
    template <class TList>
    void calcEnergy(TList& _particles, HermiteInteraction& _interaction, const bool _initial_flag) {
        ekin = epot = etot = 0.0;
        const int n = _particles.getSize();
        for (int i=0; i<n; i++) {
            auto& pi = _particles[i];
            ekin += pi.mass* (pi.vel[0]*pi.vel[0] + pi.vel[1]*pi.vel[1] + pi.vel[2]*pi.vel[2]);
            Float poti = 0.0;
            for (int j=0; j<i; j++) {
                poti +=_interaction.calcPotPair(pi, _particles[j]);
            }
            epot += poti*pi.mass;
        }
        ekin *= 0.5;
        etot = ekin + epot;

        if (_initial_flag) {
            etot0 = etot;
            de = 0.0;
        }
        else de = etot - etot0;
    }

    //! print titles of class members using column style
    /*! print titles of class members in one line for column style
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
    */
    void printColumnTitle(std::ostream & _fout, const int _width=20) {
        _fout<<std::setw(_width)<<"Time_int"
             <<std::setw(_width)<<"dE"
             <<std::setw(_width)<<"Etot"
             <<std::setw(_width)<<"Ekin"
             <<std::setw(_width)<<"Epot";
    }

    //! print data of class members using column style
    /*! print data of class members in one line for column style. Notice no newline is printed at the end
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
    */
    void printColumn(std::ostream & _fout, const int _width=20){
        _fout<<std::setw(_width)<<time
             <<std::setw(_width)<<de
             <<std::setw(_width)<<etot
             <<std::setw(_width)<<ekin
             <<std::setw(_width)<<epot;
    }
};
