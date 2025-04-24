#pragma once
#include <iostream>
#include <iomanip>
#include "Common/Float.h"

//! Changeover function class
class ChangeOver{
private:
    Float r_in_;      ///> r_in
    Float r_out_;     ///> r_out
    //Float norm_;      ///> 1.0/(r_out-r_in)
    //Float coff_;      ///> (r_out-r_in)/(r_out+r_in)
    //Float pot_off_;   ///> (1 + coff_)/r_out = 2/(r_out+r_in)
//#ifdef INTEGRATED_CUTOFF_FUNCTION
//    Float q_;         ///> r_in/r_out
//#endif
  
public:
    Float r_scale_next;   ///> scaling for changeover factor (for next step)

    ChangeOver(): r_in_(-1.0), r_out_(-1.0), r_scale_next(1.0) {}

    //! clear function
    /*! clear all data
     */
    void clear() {
        r_in_     = -1.0;
        r_out_    = -1.0;
        r_scale_next = 1.0;
    }
    
    //! check whether parameters values are correct
    /*! \return true: all correct
     */
    bool checkParams() {
        assert(r_in_>0.0);
        assert(r_out_>0.0);
        return true;
    }

    //! constructor based on inherited class
    template <class Tpars>
    ChangeOver(const Tpars& _par) {
        dataCopy(_par);
    }

    //! set r_in and r_out for changeover function
    /*!
      @param[in] _m_fac: mass factor 
      @param[in] _r_in:   changeover function inner boundary
      @param[in] _r_out:  changeover function outer boundary
    */
    void setR(const Float& _m_fac, const Float& _r_in, const Float& _r_out) {
#ifdef FIX_CHANGEOVER
        Float m_fac3 = 1.0;
#else
        Float m_fac3 = std::max(std::pow(_m_fac,(1.0/3.0)),1.0);
#endif
        r_in_     = m_fac3*_r_in;          
        r_out_    = m_fac3*_r_out;
        //norm_    = 1.0/(r_out_-r_in_);
        //coff_     = (r_out_-r_in_)/(r_out_+r_in_);
        //pot_off_  = (1.0+coff_)/r_out_;
#ifdef CHANGEOVER_DEBUG
        assert(r_in_>0.0);
        assert(r_out_>r_in_);
#endif
    }

    //! set r_in and r_out for changeover function
    /*!
      @param[in] _m_fac: mass factor 
      @param[in] _r_in:   changeover function inner boundary
      @param[in] _r_out:  changeover function outer boundary
    */
    void setR(const Float& _r_in, const Float& _r_out) {
        r_in_     = _r_in;          
        r_out_    = _r_out;
        //norm_    = 1.0/(_r_out-_r_in);
        //coff_     = (_r_out-_r_in)/(_r_out+_r_in);
        //pot_off_  = (1.0+coff_)/_r_out;
#ifdef CHANGEOVER_DEBUG
        assert(_r_in>0.0);
        assert(_r_out>_r_in);
#endif
    }
    //! update radius based on r_scale_next, reset r_scale_next to 1
    void updateWithRScale() {
        setR(r_in_*r_scale_next, r_out_*r_scale_next);
        r_scale_next = 1.0;
    }
     
    //! get r_in
    /*! \return r_in
     */
    const Float &getRin() const {
        return r_in_;
    }

    //! get r_out
    /*! \return r_out
     */
    const Float &getRout() const {
        return r_out_;
    }

    void print(std::ostream & _fout) const{
        _fout<<" r_in="<<r_in_
             <<" r_out="<<r_out_;
    }

    //! print titles of class members using column style
    /*! print titles of class members in one line for column style
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
     */
    static void printColumnTitle(std::ostream & _fout, const int _width=20) {
        _fout<<std::setw(_width)<<"r_in"
             <<std::setw(_width)<<"r_out";
    }

    //! print column title with meaning (each line for one column)
    /*! @param[out] _fout: std::ostream output object
      @param[in] _counter: offset of the number counter for each line to indicate the column index (defaulted 0)
      @param[in] _offset: the printing whitespace offset for each line (defaulted 0)
      \return: the total counter of columns
     */
    static int printTitleWithMeaning(std::ostream & _fout, const int _counter=0, const int _offset=0) {
        int counter = _counter;
        counter++;
        _fout<<std::setw(_offset)<<" "<<counter<<". r_in: inner changeover radius (0.0)\n";
        counter++;
        _fout<<std::setw(_offset)<<" "<<counter<<". r_out: outer changeover radius (0.0)\n";
        return counter;
    }

    //! print data of class members using column style
    /*! print data of class members in one line for column style. Notice no newline is printed at the end
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
     */
    void printColumn(std::ostream & _fout, const int _width=20){
        _fout<<std::setw(_width)<<r_in_
             <<std::setw(_width)<<r_out_;
    }

    //! write class data to file with binary format
    /*! @param[in] _fp: FILE type file for output
     */
    void writeBinary(FILE *_fp) const {
        fwrite(this, sizeof(Float),2,_fp);
    }

    //! read class data to file with binary format
    /*! @param[in] _fp: FILE type file for reading
     */
    void readBinary(FILE *_fin) {
        size_t rcount = fread(this, sizeof(Float),2,_fin);
        if (rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            abort();
        }
        setR(1.0, r_in_, r_out_);
    }

    //! write class data to file with binary format
    /*! @param[in] _fp: FILE type file for output
     */
    void writeAscii(FILE *_fp) const {
        fprintf(_fp, "%26.17e %26.17e ", 
                this->r_in_, this->r_out_);
    }

    //! read class data to file with binary format
    /*! @param[in] _fp: FILE type file for reading
     */
    void readAscii(FILE *_fin) {
        PS::S64 rcount=fscanf(_fin, "%lf %lf ",
                              &this->r_in_, &this->r_out_);
        if (rcount<2) {
            std::cerr<<"Error: Data reading fails! requiring data number is 2, only obtain "<<rcount<<".\n";
            abort();
        }
        setR(1.0, r_in_, r_out_);
    }

    //! copy data from inherited class object
    /*!
      @param[in] _par: object that contain data (should be inherited class)
     */
    template <class Tpars>
    void dataCopy(const Tpars& _par) {
        r_in_   = _par.r_in_   ;
        r_out_  = _par.r_out_  ;
        //norm_   = _par.norm_   ;
        //coff_   = _par.coff_   ;
        //pot_off_= _par.pot_off_;
        r_scale_next= _par.r_scale_next;
//#ifdef INTEGRATED_CUTOFF_FUNCTION
//        q_      = _par.q_;
//#endif        
    }


// Integrated potential changeover function
#ifdef INTEGRATED_CUTOFF_FUNCTION 
    //! changeover function for potential 
    /*! @param[in] _dr: particle separation
      \return \f$ \int{W_0(x) dr} \f$
     */
    inline Float calcPotW(const Float& _dr) const {
#ifdef CHANGEOVER_DEBUG
        assert(r_in_>0.0);
        assert(r_out_>r_in_);
#endif
        Float dr_rout = _dr/r_out_;
        Float q = r_in_/r_out_;
        Float q2 = q*q;
        Float q3 = q2*q;
        Float q4 = q2*q2;
        Float q5 = q3*q2;
        Float q6 = q3*q3;
        Float q7 = q4*q3;
        Float denominator = (q-1.0)*(q-1.0)*(q-1.0)*(q-1.0)*(q-1.0)*(q-1.0)*(q-1.0);
        Float A7 = 20.0/denominator/-6;
        Float A6 = (-70.0*q - 70.0)/denominator/-5;
        Float A5 = (84.0*q2 + 252.0*q + 84.0)/denominator/-4;
        Float A4 = (-35.0*q3 - 315.0*q2 - 315.0*q - 35.0)/denominator/-3;
        Float A3 = (140.0*q3 + 420.0*q2 + 140.0*q)/denominator/-2;
        Float A2 = (-210*q3 - 210.0*q2)/denominator/-1;
        Float A1 = (140*q3)/denominator*-1;
        Float A0 = (-35.0*q4 + 21.0*q5 - 7.0*q6 + q7)/denominator;
        Float x = 1.0; // x=rout/rout
        Float B1 = 1.0 - ( (((((((A7*x + A6)*x + A5)*x + A4)*x + A3)*x + A2)*x + A1*log(x))*x) + A0 ); // to W(r>rout) = 1.0
        Float A1_dash = -7*(60*q3*log(q) - q6 + 9.0*q5 - 45.0*q4 + 45.0*q2 - 9.0*q + 1.0)/(3.0*denominator);
        if(dr_rout <= q) return 1.0 - A1_dash*dr_rout;
        else if(dr_rout >= 1.0) return 0.0;
        else return 1.0 - (((((((A7*dr_rout + A6)*dr_rout + A5)*dr_rout + A4)*dr_rout + A3)*dr_rout + A2)*dr_rout + A1*log(dr_rout) + B1)*dr_rout) - A0;
    }

    //! changeover function (Poly 3rd) for force
    /*! @param[in] _dr: particle separation
      \f$ x = \frac{dr - r_{in}{r_{out} - r_{in}} \f$ \n 
      \f$ W_0(x) = x^4 (-35 + 84 x - 70 x^2 + 20 x^3) + 1 \f$ 
      \return \f$ W_0(x) \f$
     */
    inline Float calcAcc0W(const Float& _dr) const {
#ifdef CHANGEOVER_DEBUG
        assert(r_in_>0.0);
        assert(r_out_>r_in_);
#endif
        Float x = (_dr - r_in_)/(r_out_-r_in_);
        x = (x < 1.0) ? x : 1.0;
        x = (x > 0.0) ? x : 0.0;
        Float x2 = x*x;
        Float x4 = x2*x2;
        Float k = 1-(((-20.0*x+70.0)*x-84.0)*x+35.0)*x4;
        return k;
    }

    //! changeover function for force derivative
    /*!
      @param[in] _dr: particle separation
      @param[in] _drdot: time derivation of _dr
      \f$ x = \frac{dr - r_{in}{r_{out} - r_{in}} \f$ \n 
      \f$ \frac{dx}{dt} = \frac{1}{r_{out} - r_{in}} \frac{dr}{dt} \f$ \n 
      \f$ W_1(x) = x^4 (-35 + 84 x - 70 x^2 + 20 x^3) + 1 \f$ 
      \return \f$ W_1(x) dx/dt \f$
     */
    inline Float calcAcc1W(const Float& _dr, const Float& _drdot) const {
#ifdef CHANGEOVER_DEBUG
        assert(r_in_>0.0);
        assert(r_out_>r_in_);
#endif
        Float norm = 1.0/(r_out_-r_in_);
        Float x = (_dr - r_in_)*norm;
        Float xdot = norm*_drdot;
        Float kdot = 0.0;
        if(x <= 0.0)
            kdot = 0.0;
        else if(1.0 <= x)
            kdot = 0.0;
        else{
            Float x2 = x*x;
            Float x3 = x2*x;
            Float x4 = x2*x2;
            Float x5 = x4*x;
            Float x6 = x4*x2;
            kdot = -(-140.0*x6 + 420.0*x5 - 420.0*x4 + 140.0*x3) * xdot;
        }
        return kdot;
    }

#else
    // changeover function start from potential

    //! changeover function for potential 
    /*! @param[in] _dr: particle separation
      \f$ x = \frac{dr - r_{in}{r_{out} - r_{in}} \f$ \n 
      \f$ R_a = \frac{r_{in}}{r_{out} - r_{in}} \f$ \n 
      \f$ W_pot(x) = \frac{x^5}{2 R_a + 1} (5 x^3 - 20 x^2 + 28 x - 14) \f$ 
      \return \f$ W_pot(x) \f$
     */
    inline Float calcPotW(const Float& _dr) const {
#ifdef CHANGEOVER_DEBUG
        assert(r_in_>0.0);
        assert(r_out_>r_in_);
#endif
        Float x = (_dr - r_in_)/(r_out_-r_in_);
        Float coff = (r_out_-r_in_)/(r_out_+r_in_);
        Float k = coff*(1.0 - 2*x);
        if(x >= 1.0 ) k = 0.0;
        else if(x > 0.0) {
            Float x2 = x*x;
            Float x3 = x2*x;
            Float x5 = x2*x3;
            k -= coff*x5*(5.0*x3 - 20.0*x2 + 28.0*x - 14.0);
        }
        return k;
    }

    //! changeover function for force
    /*! @param[in] _dr: particle separation  
      \f$ x = \frac{dr - r_{in}{r_{out} - r_{in}} \f$ \n 
      \f$ R_a = \frac{r_{in}}{r_{out} - r_{in}} \f$ \n 
      \f$ W_0(x) = (x-1)^4 (1 + 4 x + 10 x^2 + 20 x^3 + \frac{35 x^4}{2 R_a + 1} \f$ 
      \return \f$ W_0(x) \f$
     */
    inline Float calcAcc0W(const Float& _dr) const {
#ifdef CHANGEOVER_DEBUG
        assert(r_in_>0.0);
        assert(r_out_>r_in_);
#endif
        Float x = (_dr - r_in_)/(r_out_-r_in_);
        x = (x < 1.0) ? x : 1.0;
        x = (x > 0.0) ? x : 0.0;
        Float x_1 = x - 1;
        Float x_2 = x_1*x_1;
        Float x_4 = x_2*x_2;
        Float x2 = x*x;
        Float x3 = x2*x;
        Float x4 = x2*x2;
        Float coff = (r_out_-r_in_)/(r_out_+r_in_);
        Float k = x_4*(1.0 + 4.0*x + 10.0*x2 + 20.0*x3 + 35.0*coff*x4);
        return k;
    }

    //! changeover function for force derivative
    /*!
      @param[in] _dr: particle separation
      @param[in] _drdot: time derivation of _dr
      \f$ R_a = \frac{r_{in}}{r_{out} - r_{in}} \f$ \n 
      \f$ x = \frac{dr - r_{in}{r_{out} - r_{in}} \f$ \n 
      \f$ \frac{dx}{dt} = \frac{1}{r_{out} - r_{in}} \frac{dr}{dt} \f$ \n 
      \f$ W_1(x) = \frac{280 x^3 (R_a + x)(x-1)^3}{2 R_a + 1} \f$ 
      \return \f$ W_1(x) dx/dt \f$
     */
    inline Float calcAcc1W(const Float& _dr, const Float& _drdot) const {
#ifdef CHANGEOVER_DEBUG
        assert(r_in_>0.0);
        assert(r_out_>r_in_);
#endif
        Float norm = 1.0/(r_out_-r_in_);
        Float x = (_dr - r_in_)*norm;
        Float xdot = norm*_drdot;
        Float kdot = 0.0;
        if(x > 0.0 && x < 1.0) {
            Float x3 = x*x*x;
            Float x_1 = x - 1;
            Float x_3 = x_1*x_1*x_1;
            Float coff = (r_out_-r_in_)/(r_out_+r_in_);
            kdot = coff*280.0*x3*(r_in_*norm + x)*x_3*xdot;
        }
        return kdot;
    }

#endif

    //! calculate changeover function Pot by selecting maximum rout
    static Float calcPotWTwo(const ChangeOver& _ch1, const ChangeOver& _ch2, const Float& _dr) {
        if (_ch1.getRout()> _ch2.getRout()) return _ch1.calcPotW(_dr);
        else return _ch2.calcPotW(_dr);
    }

    //! calculate changeover function Acc0 by selecting maximum rout
    static Float calcAcc0WTwo(const ChangeOver& _ch1, const ChangeOver& _ch2, const Float& _dr) {
        if (_ch1.getRout()> _ch2.getRout()) return _ch1.calcAcc0W(_dr);
        else return _ch2.calcAcc0W(_dr);
    }

    //! calculate changeover function Acc1 by selecting maximum rout
    static Float calcAcc1WTwo(const ChangeOver& _ch1, const ChangeOver& _ch2, const Float& _dr, const Float& _drdot) {
        if (_ch1.getRout()> _ch2.getRout()) return _ch1.calcAcc1W(_dr, _drdot);
        else return _ch2.calcAcc1W(_dr, _drdot);
    }
    
};

