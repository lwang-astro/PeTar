#pragma once

#include "Common/Float.h"
#include "hermite_interaction.hpp"

class HermiteInformation{
public:
    Float time_origin; // time of origin

    //! check whether parameters values are correct
    /*! \return true: all correct
     */
    bool checkParams() {
        return true;
    }        

    //! print titles of class members using column style
    /*! print titles of class members in one line for column style
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
    */
    void printColumnTitle(std::ostream & _fout, const int _width=20) {
        _fout<<std::setw(_width)<<"Time_org";
    }

    //! print data of class members using column style
    /*! print data of class members in one line for column style. Notice no newline is printed at the end
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
    */
    void printColumn(std::ostream & _fout, const int _width=20){
        _fout<<std::setw(_width)<<time_origin;
    }
};
