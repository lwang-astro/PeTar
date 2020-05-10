#pragma once
#include <iostream>
#include <iomanip>
#include <cassert>
#include "../src/io.hpp"

extern "C" {
    extern struct{
        double neta;  ///> the Reimers mass-loss coefficent (neta*4x10^-13; 0.5 normally).
        double bwind; ///> the binary enhanced mass loss parameter (inactive for single).
        double hewind; ///> the helium star mass loss factor (1.0 normally). 
        double mxns;  ///> the maximum NS mass (1.8, nsflag=0; 2.5, nsflag>=1). 
    } value1_;

    extern struct{
        int idum; ///> the random number seed used in the kick routine. 
    } value3_;

    extern struct{
        double sigma; ///> the dispersion in the Maxwellian for the SN kick speed
        int bhflag;   ///> BH kick 
    } value4_;

    extern struct{
        int ceflag; ///> common envelope model
        int tflag;  ///> tidal circularization 
        int ifflag; ///> if > 0 uses WD IFMR of HPE, 1995, MNRAS, 272, 800 (0).  (not work any more)
        int nsflag; ///> NS/BH formation
        int wdflag; ///> WD formation
    } flags_;

    extern struct{
        int psflag; ///> PPSN condition
        int kmech;  ///> kick mechanism
        int ecflag; ///> ECS switcher
    } flags2_;

    extern struct{
        double pts1;  ///> time step of MS (0.05)
        double pts2;  ///> time step of GB, CHeB, AGB, HeGB (0.01)
        double pts3;  ///> time step of HG, HeMS (0.02)
    } points_;


    //! function for initial metallicity parameters
    void zcnsts_(double* z, double* zpars);

    //! SSE function for evolving one star
    void evolv1_(int* kw, double* mass, double* mt, double* r, double* lum, double* mc, double* rc, double* menv, double* renv, double* ospin,
                 double* epoch, double* tm, double* tphys, double* tphysf, double* dtp, double* z, double* zpars, double* vs);

    void star_(int* kw, double* mass, double* mt, double* tm, double* tn, double* tscls, double* lums, double* GB, double* zpars);

    void deltat_(int* kw, double* age, double* tm, double* tn, double* tscls, double* dt, double* dtr);
}

//! SSE/BSE star parameter for saving
struct StarParameter{
    int kw;       ///> stellar type
    double m0;    ///> Initial stellar mass in solar units
    double mt;    ///> Current mass in solar units (used for R)
    double r;     ///> Stellar radius in solar units
    double ospin;  ///> spin of star
    double epoch;  ///> starting time of one evolution phase, age = tphys - epoch
    double tphys;  ///> physical evolve time in Myr

    //! print titles of class members using column style
    /*! print titles of class members in one line for column style
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
     */
    static void printColumnTitle(std::ostream & _fout, const int _width=20) {
        _fout<<std::setw(_width)<<"type"
             <<std::setw(_width)<<"mass0[Msun]"
             <<std::setw(_width)<<"mass[Msun]"
             <<std::setw(_width)<<"radius[Rsun]"
             <<std::setw(_width)<<"spin"
             <<std::setw(_width)<<"epoch[Myr]"
             <<std::setw(_width)<<"time[Myr]";
    }    

    //! print data of class members using column style
    /*! print data of class members in one line for column style. Notice no newline is printed at the end
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
     */
    void printColumn(std::ostream & _fout, const int _width=20){
        _fout<<std::setw(_width)<<kw
             <<std::setw(_width)<<m0
             <<std::setw(_width)<<mt
             <<std::setw(_width)<<r
             <<std::setw(_width)<<ospin
             <<std::setw(_width)<<epoch
             <<std::setw(_width)<<tphys;
    }
};

//! SSE/BSE star parameter for output
struct StarParameterOut{
    double lum;   ///> Landmark luminosities 
    double mc;    ///> core mass in solar units 
    double rc;    ///> core radius in solar units (output)
    double menv;  ///> mass of convective envelope 
    double renv;  ///> radius of convective envelope
    double tm;   ///> Main sequence lifetime
    double vkick[3]; ///> kick velocity for NS/BH formation

    //! print titles of class members using column style
    /*! print titles of class members in one line for column style
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
     */
    static void printColumnTitle(std::ostream & _fout, const int _width=20) {
        _fout<<std::setw(_width)<<"luminosity"
             <<std::setw(_width)<<"mass_core[Msun]"
             <<std::setw(_width)<<"radius_core[Rsun]"
             <<std::setw(_width)<<"mass_CE[Msun]"
             <<std::setw(_width)<<"radius_CE[Rsun]"
             <<std::setw(_width)<<"time_MS[Myr]"
             <<std::setw(_width)<<"vkick.x[km/s]"
             <<std::setw(_width)<<"vkick.y[km/s]"
             <<std::setw(_width)<<"vkick.z[km/s]";
    }
    
    //! print data of class members using column style
    /*! print data of class members in one line for column style. Notice no newline is printed at the end
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
     */
    void printColumn(std::ostream & _fout, const int _width=20){
        _fout<<std::setw(_width)<<lum
             <<std::setw(_width)<<mc
             <<std::setw(_width)<<rc
             <<std::setw(_width)<<menv
             <<std::setw(_width)<<renv
             <<std::setw(_width)<<tm
             <<std::setw(_width)<<vkick[0]
             <<std::setw(_width)<<vkick[1]
             <<std::setw(_width)<<vkick[2];
    }
};

class IOParamsBSE{
public:
    IOParamsContainer input_par_store;
    IOParams<double> neta;
    IOParams<double> bwind;
    IOParams<double> hewind;
    //IOParams<double> mxns;
    IOParams<double> sigma;
    IOParams<int> ceflag;
    IOParams<int> tflag;
    //IOParams<int> ifflag;
    IOParams<int> wdflag;
    IOParams<int> bhflag;
    IOParams<int> nsflag;
    IOParams<int> psflag;
    IOParams<int> kmech;
    IOParams<int> ecflag;
    IOParams<double> pts1;
    IOParams<double> pts2;
    IOParams<double> pts3;
    IOParams<int> idum;

    bool print_flag;

    IOParamsBSE(): input_par_store(),
                   neta  (input_par_store, 0.5, "Reimers mass-loss coefficent [neta*4x10^-13]"),
                   bwind (input_par_store, 0.0, "Binary enhanced mass loss parameter; inactive for single"),
                   hewind(input_par_store, 1.0, "Helium star mass loss factor"),
                   //mxns  (input_par_store, 1.0, "Helium star mass loss factor"),
                   sigma (input_par_store, 265.0, "the dispersion in the Maxwellian for the SN kick speed [km/s]"),
                   ceflag(input_par_store, 0,  "if =3, activates de Kool common-envelope model"),
                   tflag (input_par_store, 1,  "if >0, activates tidal circularisation"),
                   //ifflag(input_par_store, 2,   "if > 0 uses WD IFMR of HPE, 1995, MNRAS, 272, 800"),
                   wdflag(input_par_store, 1,  "if >0, uses WD IFMR of HPE, 1995, MNRAS, 272, 800"),
                   bhflag(input_par_store, 2,  "BH kick option: 0: no kick; 1: same as NS; 2: scaled by fallback"),
                   nsflag(input_par_store, 3,  "NS/BH foramtion options: 0: original SSE; 1: Belczynski (2002); 2: Belczynski (2008); 3: Fryer (2012) rapid SN; 4: Fryer (2012) delayed SN; 5: Eldridge & Tout (2004)"),
                   psflag(input_par_store, 1,  "PPSN condition (Leung 2019): 0: no PPSN; 1: strong; 2: moderate; 3: weak"),
                   kmech (input_par_store, 1,  "Kick mechanism: 1: standard momentum-conserving; 2: convection-asymmetry-driven; 3: collapse-asymmerty-driven; 4: neutrino driven"),
                   ecflag(input_par_store, 1,  "if >0, ECS is switched on"),
                   pts1  (input_par_store, 0.05, "time step of MS"),
                   pts2  (input_par_store, 0.01, "time step of GB, CHeB, AGB, HeGB"),
                   pts3  (input_par_store, 0.02, "time step of HG, HeMS"),
                   idum  (input_par_store, 1234, "random number seed used by the kick routine"),
                   print_flag(false) {}

    int read(int argc, char *argv[]) {
        static int sse_flag=-1;
        static struct option long_options[] = {
            {"neta",   required_argument, &sse_flag, 0},  
            {"bwind",  required_argument, &sse_flag, 1},  
            {"hewind", required_argument, &sse_flag, 2},  
          //{"mxns",   required_argument, &sse_flag, 3}, 
            {"sigma",  required_argument, &sse_flag, 4},
            {"ceflag", required_argument, &sse_flag, 5},
            {"tflag",  required_argument, &sse_flag, 6},
          //{"ifflag", required_argument, &sse_flag, 7},
            {"wdflag", required_argument, &sse_flag, 8},
            {"bhflag", required_argument, &sse_flag, 9}, 
            {"nsflag", required_argument, &sse_flag, 10}, 
            {"psflag", required_argument, &sse_flag, 11},
            {"kmech",  required_argument, &sse_flag, 12},
            {"ecflag", required_argument, &sse_flag, 13},
            {"pts1",   required_argument, &sse_flag, 14},
            {"pts2",   required_argument, &sse_flag, 15},       
            {"pts3",   required_argument, &sse_flag, 16},
            {"idum",   required_argument, &sse_flag, 17}, 
            {"help",   no_argument,       0, 'h'},
            {0,0,0,0}
        };

        int copt;
        int option_index;
        optind = 1;
        while ((copt = getopt_long(argc, argv, "h", long_options, &option_index)) != -1) 
            switch (copt) {
            case 0:
                switch (sse_flag) {
                case 0:
                    neta.value = atof(optarg);
                    if(print_flag) neta.print(std::cout);
                    break;
                case 1:
                    bwind.value = atof(optarg);
                    if(print_flag) bwind.print(std::cout);
                    break;
                case 2:
                    hewind.value = atof(optarg);
                    if(print_flag) hewind.print(std::cout);
                    break;
                //case 3:
                //    mxns.value = atof(optarg);
                //    if(print_flag) mxns.print(std::cout);
                //    break;
                case 4:
                    sigma.value = atof(optarg);
                    if(print_flag) sigma.print(std::cout);
                    break;
                case 5:
                    ceflag.value = atof(optarg);
                    if(print_flag) ceflag.print(std::cout);
                    break;
                case 6:
                    tflag.value = atof(optarg);
                    if(print_flag) tflag.print(std::cout);
                    break;
                case 8:
                    wdflag.value = atof(optarg);
                    if(print_flag) wdflag.print(std::cout);
                    break;
                case 9:
                    bhflag.value = atof(optarg);
                    if(print_flag) bhflag.print(std::cout);
                    break;
                case 10:
                    nsflag.value = atof(optarg);
                    if(print_flag) nsflag.print(std::cout);
                    break;
                case 11:
                    psflag.value = atof(optarg);
                    if(print_flag) psflag.print(std::cout);
                    break;
                case 12:
                    kmech.value = atof(optarg);
                    if(print_flag) kmech.print(std::cout);
                    break;
                case 13:
                    ecflag.value = atof(optarg);
                    if(print_flag) ecflag.print(std::cout);
                    break;
                case 14:
                    pts1.value = atof(optarg);
                    if(print_flag) pts1.print(std::cout);
                    break;
                case 15:
                    pts2.value = atof(optarg);
                    if(print_flag) pts2.print(std::cout);
                    break;
                case 16:
                    pts3.value = atof(optarg);
                    if(print_flag) pts3.print(std::cout);
                    break;
                case 17:
                    idum.value = atof(optarg);
                    if(print_flag) idum.print(std::cout);
                    break;
                default:
                    break;
                }
            case 'h':
                if(print_flag){
                    std::cout<<"SSE/BSE options:"<<std::endl;
                    std::cout<<"       Option defaulted values are shown after ':'"<<std::endl;
                    std::cout<<"        --neta:   [D] "<<neta<<std::endl
                             <<"        --bwind:  [D] "<<bwind<<std::endl
                             <<"        --hewind: [D] "<<hewind<<std::endl
                        //<<"        --mxns:   [D] "<<mxns<<std::endl
                             <<"        --sigma:  [D] "<<sigma<<std::endl
                             <<"        --ceflag: [I] "<<ceflag<<std::endl
                             <<"        --tflag:  [I] "<<tflag<<std::endl
                             <<"        --wdflag: [I] "<<wdflag<<std::endl
                             <<"        --bhflag: [I] "<<bhflag<<std::endl
                             <<"        --nsflag: [I] "<<nsflag<<std::endl
                             <<"        --psflag: [I] "<<psflag<<std::endl
                             <<"        --kmech:  [I] "<<kmech<<std::endl
                             <<"        --ecflag: [I] "<<ecflag<<std::endl
                             <<"        --pts1:   [D] "<<pts1<<std::endl
                             <<"        --pts2:   [D] "<<pts2<<std::endl
                             <<"        --pts3:   [D] "<<pts3<<std::endl
                             <<"        --idum:   [I] "<<idum<<std::endl;
                }
                return -1;
            default:
                break;
            }
        return 0;
    }    
};

//! SSE/BSE interface manager
class BSEManager{
public:
    double z, zpars[20]; ///> metallicity parameters
    double tscale; ///> time scaling factor from NB to Myr (t[Myr]=t[NB]*tscale)

    BSEManager(): z(0.0), zpars{0}, tscale(0.0) {}

    bool checkParams() {
        assert(z>0.0);
        assert(tscale>0.0);
        return true;
    }

    //! initial SSE/BSE global parameters
    void initial(const IOParamsBSE& _input, const double _z) {
        // common block
        value1_.neta  = _input.neta.value;
        value1_.bwind = _input.bwind.value;
        value1_.hewind= _input.hewind.value;
        value1_.mxns  = 1.8;
        if (_input.nsflag.value>0) value1_.mxns = 2.5;
        
        value4_.sigma  = _input.sigma.value;
        value4_.bhflag = _input.bhflag.value;

        flags_.ceflag = _input.ceflag.value;
        flags_.tflag  = _input.tflag.value;
        //flags_.ifflag = _input.ifflag.value;
        flags_.wdflag = _input.wdflag.value;
        flags_.nsflag = _input.nsflag.value;

        flags2_.psflag = _input.psflag.value;
        flags2_.kmech  = _input.kmech.value;
        flags2_.ecflag = _input.ecflag.value;

        points_.pts1 = _input.pts1.value;
        points_.pts2 = _input.pts2.value;
        points_.pts3 = _input.pts3.value;

        // Set parameters which depend on the metallicity 
        z = _z;
        zcnsts_(&z, zpars);
        value3_.idum = (_input.idum.value>0)? -_input.idum.value: _input.idum.value;
    }

    //! initial SSE/BSE global parameters
    void initial(const double _z, 
                 const double _neta, const double _bwind, const double _hewind, 
                 const double _sigma, 
                 const int _ceflag, const int _tflag, const int _wdflag, const int _bhflag, const int _nsflag, 
                 const double _pts1, const double _pts2, const double _pts3, 
                 const int _idum) {

        // common block
        value1_.neta  = _neta;
        value1_.bwind = _bwind;
        value1_.hewind= _hewind;
        value1_.mxns  = 1.8;
        if (_nsflag==1) value1_.mxns = 3.0;
        
        value4_.sigma = _sigma;
        value4_.bhflag = _bhflag;

        flags_.ceflag = _ceflag;
        flags_.tflag = _tflag;
        //flags_.ifflag = _ifflag;
        flags_.wdflag = _wdflag;
        flags_.nsflag = _nsflag;

        points_.pts1 = _pts1;
        points_.pts2 = _pts2;
        points_.pts3 = _pts3;

        // Set parameters which depend on the metallicity 
        z = _z;
        zcnsts_(&z, zpars);
        if (value3_.idum>0) value3_.idum = -_idum;
        else value3_.idum = _idum;
    }
    

    //! call SSE evolve1 for single star
    /*!
      @param[in,out] _star: star parameter
      @param[in] _tphysf: physical time in Myr to evolve
     */
    StarParameterOut evolveStar(StarParameter& _star, const double _time_nb) {
        StarParameterOut out;
        double tphysf = _time_nb*tscale;
        double dtp=tphysf*100.0+1000.0;
        evolv1_(&_star.kw, &_star.m0, &_star.mt, &_star.r, 
                &out.lum, &out.mc, &out.rc, &out.menv, &out.renv, 
                &_star.ospin, &_star.epoch, 
                &out.tm, &_star.tphys, &tphysf, &dtp, &z, zpars, out.vkick);

        return out;
    }

    //! get next time step to check in Myr
    double getTimeStep(StarParameter& _star) {
        double tm, tn, tscls[20], lums[10], gb[10], dtm, dtr;
        
        // obtain star parameters
        star_(&_star.kw, &_star.m0, &_star.mt, &tm, &tn, tscls, lums, gb, zpars);

        // get next step
        double age = _star.tphys-_star.epoch;
        deltat_(&_star.kw, &age, &tm, &tn, tscls, &dtm, &dtr);

        assert(dtr>0.0);
        assert(dtm>0.0);

        return std::min(dtr, dtm);
    }
    
};
