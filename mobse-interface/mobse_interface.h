#pragma once
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cstdio>
#include <string>
#include <getopt.h>
#include "../src/io.hpp"

extern "C" {
    extern struct{
        double neta;  ///> the Reimers mass-loss coefficent (neta*4x10^-13; 0.5 normally).
        double bwind; ///> the binary enhanced mass loss parameter (0.0 normally).
        double hewind; ///> the helium star mass loss factor (1.0 normally). 
    } value1_;

    extern struct{
        double alpha; ///>the common-envelope efficiency parameter (3.0).
        double lambda; ///>the binding energy factor for common envelope evolution (0.1).
    } value2_;

    extern struct{
        int idum; ///> the random number seed used in the kick routine. 
    } value3_;

    extern struct{
        double sigma1; ///> the dispersion in the Maxwellian for the CCSN kick speed
        double sigma2; ///> the dispersion in the Maxwellian for the ECSN kick speed
        double mxns;  ///> the maximum NS mass (1.8, nsflag=0; 3.0, nsflag>=1). 
        int bhflag;   ///> BH kick (3 normally Giacobbo & Mapelli ApJ 2020)
    } value4_;

    extern struct{
        double beta;   ///> wind velocity factor: proportional to vwind**2 (1/8). 
        double xi;     ///> wind accretion efficiency factor (1.0). 
        double bhwacc;   ///> the Bondi-Hoyle wind accretion factor (3/2). 
        double epsnov; ///>the fraction of accreted matter retained in nova eruption (0.001). 
        double eddfac; ///> Eddington limit factor for mass transfer (1.0).
        double gamma;  ///> the angular momentum factor for mass lost during Roche (-1.0). 
    } value5_;

    extern struct{
        int ceflag; ///> common envelope model (0 normally)
        int tflag;  ///> tidal circularization (1 normally)
        int ifflag; ///> if > 0 uses WD IFMR of HPE, 1995, MNRAS, 272, 800 (0).  (not work any more)
        int nsflag; ///> NS/BH formation (3 normally Delayed)
        int wdflag; ///> WD formation (1 normally)
        int piflag; ///> (Pulsational) Pair Instability (1 normally)
    } flags_;

//    extern struct{
//        int psflag; ///> PPSN condition
//        int kmech;  ///> kick mechanism
//        int ecflag; ///> ECS switcher
//    } flags2_;

    extern struct{
        double pts1;  ///> time step of MS (0.05)
        double pts2;  ///> time step of GB, CHeB, AGB, HeGB (0.01)
        double pts3;  ///> time step of HG, HeMS (0.02)
    } points_;

    extern struct{
        int idum2;
        int iy;
        int ir[32];
    } rand3_;

    //! function for initial metallicity parameters
    void zcnsts_(double* z, double* zpars);

    //!function for collison matrix
    void instar_();

    //! function for evolving one star
    void evolv1_(int* kw, double* mass, double* mt, double* r, double* lum, double* mc, double* rc, double* menv, double* renv, double* ospin,
                 double* epoch, double* tm, double* tphys, double* tphysf, double* dtp, double* z, double* zpars, double* vkick);

    //! function for evolving one binary
    void evolv2_(int* kw, double* mass, double* mt, double* r, double* lum, double* mc, double* rc, double* menv, double* renv, double* ospin,
                 double* epoch, double* tm, double* tphys, double* tphysf, double* dtp, double* z, double* zpars, 
                 double* period, double* ecc, double* bse_event, double* vkick);

    void star_(int* kw, double* mass, double* mt, double* tm, double* tn, double* tscls, double* lums, double* GB, double* zpars);

    void deltat_(int* kw, double* age, double* tm, double* tn, double* tscls, double* dt, double* dtr);

    void mix_(double* m0, double* mt, double* age, int* kw, double* zpars, int* krol);

    void printconst_();
}

//! MOBSE star parameter for saving
struct StarParameter{
    long long int kw;       ///> stellar type
    double m0;    ///> Initial stellar mass in solar units
    double mt;    ///> Current mass in solar units (used for R)
    double r;     ///> Stellar radius in solar units
    double mc;    ///> core mass in solar units 
    double rc;    ///> core radius in solar units (output)
    double ospin;  ///> spin of star
    double epoch;  ///> starting time of one evolution phase, age = tphys - epoch
    double tphys;  ///> physical evolve time in Myr
    double lum;    ///> Landmark luminosities 

    //! initial zero age main sequence
    /*!
      @param[in] _mass: initial mass
      @param[in] _kw: initial type (default: 1: MS)
      @param[in] _ospin: initial spin (default: 0.0)
      @param[in] _epoch: initial age for the given type (default: 0.0)
     */
    void initial(double _mass, int _kw=1, double _ospin=0.0, double _epoch=0.0) {
        kw = _kw;
        m0 = _mass;
        mt = _mass;
        r  = 0.0;
        mc = 0.0;
        rc = 0.0;
        ospin = _ospin;
        epoch = _epoch;
        tphys = _epoch;
    }

    //! write class data with ASCII format
    /*! @param[in] _fout: file IO for write
     */
    void writeAscii(FILE* fp) const{
        fprintf(fp, "%lld %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e ",
                this->kw, this->m0, this->mt, this->r, this->mc, this->rc, this->ospin, this->epoch, this->tphys, this->lum);
    }

    //! read class data with ASCII format
    /*! @param[in] _fin: file IO for read
     */
    void readAscii(FILE* fp) {
        int rcount=fscanf(fp, "%lld %lf %lf %lf %lf %lf %lf %lf %lf %lf ",
                          &this->kw, &this->m0, &this->mt, &this->r, &this->mc, &this->rc, &this->ospin, &this->epoch, &this->tphys, & this->lum);
        if(rcount<10) {
            std::cerr<<"Error: Data reading fails! requiring data number is 10, only obtain "<<rcount<<".\n";
            abort();
        }
    }

    //! for print in one line
    void print(std::ostream & fout) const{
        fout<<" type= "<<kw
            <<" m0[M*]= "<<m0
            <<" m[M*]= "<<mt
            <<" rad[R*]= "<<r
            <<" mc[M*]= "<<mc
            <<" rc[M*]= "<<rc
            <<" spin= "<<ospin
            <<" epoch= "<<epoch
            <<" t[myr]= "<<tphys
            <<" lum[L*]= "<<lum;
    }

    //! print titles of class members using column style
    /*! print titles of class members in one line for column style
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
     */
    static void printColumnTitle(std::ostream & _fout, const int _width=20) {
        _fout<<std::setw(_width)<<"s_type"
             <<std::setw(_width)<<"s_mass0[M*]"
             <<std::setw(_width)<<"s_mass[M*]"
             <<std::setw(_width)<<"s_rad[R*]"
             <<std::setw(_width)<<"s_mcore[M*]"
             <<std::setw(_width)<<"s_rcore[R*]"
             <<std::setw(_width)<<"s_spin"
             <<std::setw(_width)<<"s_epoch[Myr]"
             <<std::setw(_width)<<"s_time[Myr]"
             <<std::setw(_width)<<"s_lum[L*]";
    }    

    //! print data of class members using column style
    /*! print data of class members in one line for column style. Notice no newline is printed at the end
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
     */
    void printColumn(std::ostream & _fout, const int _width=20) const{
        _fout<<std::setw(_width)<<kw
             <<std::setw(_width)<<m0
             <<std::setw(_width)<<mt
             <<std::setw(_width)<<r
             <<std::setw(_width)<<mc
             <<std::setw(_width)<<rc
             <<std::setw(_width)<<ospin
             <<std::setw(_width)<<epoch
             <<std::setw(_width)<<tphys
             <<std::setw(_width)<<lum;
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
        _fout<<std::setw(_offset)<<" "<<counter<<". s_type: MOBSE stellar type\n";
        counter++;
        _fout<<std::setw(_offset)<<" "<<counter<<". s_mass0: initial mass at each evolution stage [Msun]\n";
        counter++;
        _fout<<std::setw(_offset)<<" "<<counter<<". s_mass: current mass [Msun]\n";
        counter++;
        _fout<<std::setw(_offset)<<" "<<counter<<". s_rad: stellar radius [Rsun]\n";        
        counter++;
        _fout<<std::setw(_offset)<<" "<<counter<<". s_mcore: stellar core mass [Msun]\n";
        counter++;
        _fout<<std::setw(_offset)<<" "<<counter<<". s_rcore: stellar core radius [Rsun]\n";        
        counter++;
        _fout<<std::setw(_offset)<<" "<<counter<<". s_spin: stellar rotation\n";        
        counter++;
        _fout<<std::setw(_offset)<<" "<<counter<<". s_epoch: time offset at each evolution stage [Myr]\n";
        counter++;
        _fout<<std::setw(_offset)<<" "<<counter<<". s_time: physical time [Myr]\n";
        counter++;
        _fout<<std::setw(_offset)<<" "<<counter<<". s_lum: luminosity [Lsun]\n";
        return counter;
    }
};

//! MOBSE star parameter for output
struct StarParameterOut{
    long long int kw0;       ///> original type before evolution
    double dtmiss; ///> required evolution time - actually evolved time
    double menv;  ///> mass of convective envelope 
    double renv;  ///> radius of convective envelope
    double tm;   ///> Main sequence lifetime
    double vkick[4]; ///> kick velocity for NS/BH formation
    double dm;   ///> mass loss

    StarParameterOut(): kw0(0), dtmiss(0.0), menv(0.0), renv(0.0), tm(0.0), vkick{0.0}, dm(0.0) {}

    //! print titles of class members using column style
    /*! print titles of class members in one line for column style
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
     */
    static void printColumnTitle(std::ostream & _fout, const int _width=20) {
        _fout<<std::setw(_width)<<"type0"
             <<std::setw(_width)<<"m_CE[M*]"
             <<std::setw(_width)<<"r_CE[R*]"
             <<std::setw(_width)<<"t_MS[Myr]"
             <<std::setw(_width)<<"vk.x[km/s]"
             <<std::setw(_width)<<"vk.y[km/s]"
             <<std::setw(_width)<<"vk.z[km/s]"
             <<std::setw(_width)<<"vk[km/s]"
             <<std::setw(_width)<<"m_loss[M*]";
    }
    
    //! print data of class members using column style
    /*! print data of class members in one line for column style. Notice no newline is printed at the end
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
     */
    void printColumn(std::ostream & _fout, const int _width=20) const{
        _fout<<std::setw(_width)<<kw0
             <<std::setw(_width)<<menv
             <<std::setw(_width)<<renv
             <<std::setw(_width)<<tm
             <<std::setw(_width)<<vkick[0]
             <<std::setw(_width)<<vkick[1]
             <<std::setw(_width)<<vkick[2]
             <<std::setw(_width)<<vkick[3]
             <<std::setw(_width)<<dm;
    }

    //! for print in one line
    void print(std::ostream & fout) const{
        fout<<" type0= "<<kw0
            <<" menv[M*]= "<<menv
            <<" renv[R*]= "<<renv
            <<" t_MS[Myr]= "<<tm
            <<" vk.x[km/s]= "<<vkick[0]
            <<" vk.y[km/s]= "<<vkick[1]
            <<" vk.z[km/s]= "<<vkick[2]
            <<" vk[km/s]= "<<vkick[3]
            <<" m_loss[M*]= "<<dm;
    }

};

//! MOBSE event recorder class
class BinaryEvent{
public:
    double record[10][9];

    void recordInitial(const StarParameter& _p1, const StarParameter& _p2, const double _semi, const double _ecc) {
        const int init_index = getEventIndexInit();
        record[0][init_index] = std::min(_p1.tphys, _p2.tphys);
        record[1][init_index] = _p1.mt;
        record[2][init_index] = _p2.mt;
        record[3][init_index] = _p1.kw;
        record[4][init_index] = _p2.kw;
        record[5][init_index] = _semi;
        record[6][init_index] = _ecc;
        record[7][init_index] = _p1.r;
        record[8][init_index] = _p2.r;
        record[9][init_index] = 0.0;
    }
    
    int getEventNMax() const {
        return 8;
    }
    
    int getEventIndexInit() const {
        return 8;
    }

    int getType(const int index) const {
        return int(record[9][index]);
    }

    //! for print in one line
    void print(std::ostream & fout, const int index) const{
        fout<<" t[Myr]= "<<record[0][index]
            <<" m1[M*]= "<<record[1][index]
            <<" m2[M*]= "<<record[2][index]
            <<" type1= "<<int(record[3][index])
            <<" type2= "<<int(record[4][index])
            <<" semi[R*]= "<<record[5][index]
            <<" ecc= "<<record[6][index]
            <<" rad1[R*]= "<<record[7][index]
            <<" rad2[R*]= "<<record[8][index];
    }

    //! print titles of class members using column style
    /*! print titles of class members in one line for column style
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
     */
    static void printColumnTitle(std::ostream & _fout, const int _width=20) {
        _fout<<std::setw(_width)<<"t[Myr]"
             <<std::setw(_width)<<"m1[M*]"
             <<std::setw(_width)<<"m2[M*]"
             <<std::setw(_width)<<"type1"
             <<std::setw(_width)<<"type2"
             <<std::setw(_width)<<"semi[R*]"
             <<std::setw(_width)<<"ecc"
             <<std::setw(_width)<<"rad1[R*]"
             <<std::setw(_width)<<"rad2[R*]"
             <<std::setw(_width)<<"btype";
    }

    //! print data of class members using column style
    /*! print data of class members in one line for column style. Notice no newline is printed at the end
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
      @param[in] _index: event index to print
     */
    void printColumn(std::ostream & _fout, const int _index, const int _width=20) const{
        for (int i=0; i<3; i++) _fout<<std::setw(_width)<<record[i][_index];
        for (int i=3; i<5; i++) _fout<<std::setw(_width)<<int(record[i][_index]);
        for (int i=5; i<9; i++) _fout<<std::setw(_width)<<record[i][_index];
        _fout<<std::setw(_width)<<int(record[9][_index]);
    }
};

//! IO parameters for MOBSE manager
class IOParamsMOBSE{
public:
    IOParamsContainer input_par_store;
    IOParams<double> neta;
    IOParams<double> bwind;
    IOParams<double> hewind;
    IOParams<double> alpha;
    IOParams<double> lambda;
    IOParams<double> beta;
    IOParams<double> xi;
    IOParams<double> bhwacc;
    IOParams<double> epsnov;
    IOParams<double> eddfac;
    IOParams<double> gamma;
    //IOParams<double> mxns;
    IOParams<double> sigma1;
    IOParams<double> sigma2;
    IOParams<long long int> ceflag;
    IOParams<long long int> tflag;
    //IOParams<long long int> ifflag;
    IOParams<long long int> wdflag;
    IOParams<long long int> bhflag;
    IOParams<long long int> nsflag;
    IOParams<long long int> piflag;
    //IOParams<int> kmech;
    //IOParams<int> ecflag;
    IOParams<double> pts1;
    IOParams<double> pts2;
    IOParams<double> pts3;
    IOParams<long long int> idum;
    IOParams<double> tscale;
    IOParams<double> rscale;
    IOParams<double> mscale;
    IOParams<double> vscale;
    IOParams<double> z;

    bool print_flag;

    IOParamsMOBSE(): input_par_store(),
                   neta  (input_par_store, 0.5,     "mobse-neta",   "Reimers mass-loss coefficent [neta*4x10^-13]"),
                   bwind (input_par_store, 0.0,     "mobse-wind",   "Binary enhanced mass loss parameter; inactive for single"),
                   hewind(input_par_store, 1.0,     "mobse-hewind", "Helium star mass loss factor"),
                   //mxns  (input_par_store, 1.0, "Helium star mass loss factor"),
                   alpha (input_par_store, 3.0,     "mobse-alpha",  "Common-envelope efficiency parameter"),
                   lambda(input_par_store, 0.1,     "mobse-lambda", "Binding energy factor for common envelope evolution"),
                   beta  (input_par_store, 0.125,   "mobse-beta",   "wind velocity factor: proportional to vwind**2"),
                   xi    (input_par_store, 1.0,     "mobse-xi",       "wind accretion efficiency factor"),
                   bhwacc(input_par_store, 1.5,     "mobse-bwacc",  "Bondi-Hoyle wind accretion factor"),
                   epsnov(input_par_store, 0.001,   "mobse-epsnov", "The fraction of accreted matter retained in nova eruption"),
                   eddfac(input_par_store, 1.0,     "mobse-eddfac", "Eddington limit factor for mass transfer"),
                   gamma (input_par_store, -1.0,    "mobse-gamma",  "Angular momentum factor for mass lost during Roche"),
                   sigma1 (input_par_store, 265.0,  "mobse-sigma1", "Dispersion in the Maxwellian for the CCSN kick speed [km/s]"),
                   sigma2 (input_par_store, 265.0,  "mobse-sigma2", "Dispersion in the Maxwellian for the ECSN kick speed [km/s]"),
                   ceflag(input_par_store, 0,       "mobse-cflag",  "if =3, activates de Kool common-envelope model"),
                   tflag (input_par_store, 1,       "mobse-tflag",  "if >0, activates tidal circularisation"),
                   //ifflag(input_par_store, 2,   "if > 0 uses WD IFMR of HPE, 1995, MNRAS, 272, 800"),
                   wdflag(input_par_store, 1,       "mobse-wdflag", "if >0, uses WD IFMR of HPE, 1995, MNRAS, 272, 800"),
                   bhflag(input_par_store, 3,       "mobse-bhflag", "BH kick option: 0: no kick; 1: same as NS; 2: scaled by fallback; 3: Giacobbo&Mapelli (2020)"),
                   nsflag(input_par_store, 3,       "mobse-nsflag", "NS/BH foramtion options: 0: original SSE; 2: Belczynski (2008); 2: Fryer (2012) rapid SN; 2: Fryer (2012) delayed SN; 4: Belczynski (2008); 5: no SN explosion"),
                   piflag(input_par_store, 1,       "mobse-piflag", "PPSN condition (Spera et al. 2015)"),
                   //psflag(input_par_store, 1,  "PPSN condition (Belczynski 2016): 0: no PPSN; 1: strong; (Leung 2019): 2: moderate; 3: weak"),
                   //kmech (input_par_store, 1,  "Kick mechanism: 1: standard momentum-conserving; 2: convection-asymmetry-driven; 3: collapse-asymmerty-driven; 4: neutrino driven"),
                   //ecflag(input_par_store, 1,  "if >0, ECS is switched on"),
                   pts1  (input_par_store, 0.05,    "mobse-pts1",   "time step of MS"),
                   pts2  (input_par_store, 0.01,    "mobse-pts2",   "time step of GB, CHeB, AGB, HeGB"),
                   pts3  (input_par_store, 0.02,    "mobse-pts3",   "time step of HG, HeMS"),
                   idum  (input_par_store, 1234,    "mobse-idum",   "random number seed used by the kick routine"),
                   tscale(input_par_store, 1.0,     "mobse-tscale", "Time scale factor from input data unit (IN) to Myr (time[Myr]=time[IN]*tscale)"),
                   rscale(input_par_store, 1.0,     "mobse-rscale", "Radius scale factor from input data unit (IN) to Rsun (r[Rsun]=r[IN]*rscale)"),
                   mscale(input_par_store, 1.0,     "mobse-msclae", "Mass scale factor from input data unit (IN) to Msun (m[Msun]=m[IN]*mscale)"),
                   vscale(input_par_store, 1.0,     "mobse-vsclae",  "Velocity scale factor from input data unit(IN) to km/s (v[km/s]=v[IN]*vscale)"),
                   z     (input_par_store, 0.001,   "mobse-metallicity",    "Metallicity"),
                   print_flag(false) {}

    //! reading parameters from GNU option API
    /*!
      @param[in] argc: number of options
      @param[in] argv: string of options
      @param[in] opt_used_pre: already used option number from previous reading, use to correctly count the remaining argument number
      \return -1 if help is used; else the used number of argv
     */
    int read(int argc, char *argv[], const int opt_used_pre=0) {
        static int sse_flag=-1;
        const struct option long_options[] = {
            {neta.key,   required_argument, &sse_flag, 0},  
            {bwind.key,  required_argument, &sse_flag, 1},  
            {hewind.key, required_argument, &sse_flag, 2},  
          //{mxns",   required_argument, &sse_flag, 3}, 
            {alpha.key,  required_argument, &sse_flag, 22},
            {lambda.key, required_argument, &sse_flag, 23},
            {beta.key,   required_argument, &sse_flag, 24},
            {xi.key,     required_argument, &sse_flag, 25},
            {bhwacc.key, required_argument, &sse_flag, 26},
            {epsnov.key, required_argument, &sse_flag, 27},
            {eddfac.key, required_argument, &sse_flag, 28},
            {gamma.key,  required_argument, &sse_flag, 29},
            {sigma1.key, required_argument, &sse_flag, 4},
            {sigma2.key, required_argument, &sse_flag, 5},
            {ceflag.key, required_argument, &sse_flag, 6},
            {tflag.key,  required_argument, &sse_flag, 7},
          //{ifflag", required_argument, &sse_flag, 7},
            {wdflag.key, required_argument, &sse_flag, 8},
            {bhflag.key, required_argument, &sse_flag, 9}, 
            {nsflag.key, required_argument, &sse_flag, 10}, 
            {piflag.key, required_argument, &sse_flag, 11},
          //{psflag", required_argument, &sse_flag, 11},
          //{kmech",  required_argument, &sse_flag, 12},
          //{ecflag", required_argument, &sse_flag, 13},
            {pts1.key,   required_argument, &sse_flag, 14},
            {pts2.key,   required_argument, &sse_flag, 15},       
            {pts3.key,   required_argument, &sse_flag, 16},
            {idum.key,   required_argument, &sse_flag, 17}, 
            {tscale.key, required_argument, &sse_flag, 18},
            {rscale.key, required_argument, &sse_flag, 19},
            {mscale.key, required_argument, &sse_flag, 20},
            {vscale.key, required_argument, &sse_flag, 21},
            {z.key,      required_argument, 0, 'z'},
            {"help",     no_argument,       0, 'h'},
            {0,0,0,0}
        };

        int opt_used=opt_used_pre;
        int copt;
        int option_index;
        std::string fname_par;
        optind = 0;
        while ((copt = getopt_long(argc, argv, "-z:p:h", long_options, &option_index)) != -1) 
            switch (copt) {
            case 0:
                switch (sse_flag) {
                case 0:
                    neta.value = atof(optarg);
                    if(print_flag) neta.print(std::cout);
                    opt_used+=2;
                    break;
                case 1:
                    bwind.value = atof(optarg);
                    if(print_flag) bwind.print(std::cout);
                    opt_used+=2;
                    break;
                case 2:
                    hewind.value = atof(optarg);
                    if(print_flag) hewind.print(std::cout);
                    opt_used+=2;
                    break;
                //case 3:
                //    mxns.value = atof(optarg);
                //    if(print_flag) mxns.print(std::cout);
                //    break;
                case 4:
                    sigma1.value = atof(optarg);
                    if(print_flag) sigma1.print(std::cout);
                    opt_used+=2;
                    break;
                case 5:
                    sigma2.value = atof(optarg);
                    if(print_flag) sigma2.print(std::cout);
                    opt_used+=2;
                    break;
                case 6:
                    ceflag.value = atof(optarg);
                    if(print_flag) ceflag.print(std::cout);
                    opt_used+=2;
                    break;
                case 7:
                    tflag.value = atof(optarg);
                    if(print_flag) tflag.print(std::cout);
                    opt_used+=2;
                    break;
                case 8:
                    wdflag.value = atof(optarg);
                    if(print_flag) wdflag.print(std::cout);
                    opt_used+=2;
                    break;
                case 9:
                    bhflag.value = atof(optarg);
                    if(print_flag) bhflag.print(std::cout);
                    opt_used+=2;
                    break;
                case 10:
                    nsflag.value = atof(optarg);
                    if(print_flag) nsflag.print(std::cout);
                    opt_used+=2;
                    break;
                case 11:
                    piflag.value = atof(optarg);
                    if(print_flag) piflag.print(std::cout);
                    opt_used+=2;
                    break;
                case 14:
                    pts1.value = atof(optarg);
                    if(print_flag) pts1.print(std::cout);
                    opt_used+=2;
                    break;
                case 15:
                    pts2.value = atof(optarg);
                    if(print_flag) pts2.print(std::cout);
                    opt_used+=2;
                    break;
                case 16:
                    pts3.value = atof(optarg);
                    if(print_flag) pts3.print(std::cout);
                    opt_used+=2;
                    break;
                case 17:
                    idum.value = atof(optarg);
                    if(print_flag) idum.print(std::cout);
                    opt_used+=2;
                    break;
                case 18:
                    tscale.value = atof(optarg);
                    if(print_flag) tscale.print(std::cout);
                    opt_used+=2;
                    break;
                case 19:
                    rscale.value = atof(optarg);
                    if(print_flag) rscale.print(std::cout);
                    opt_used+=2;
                    break;
                case 20:
                    mscale.value = atof(optarg);
                    if(print_flag) mscale.print(std::cout);
                    opt_used+=2;
                    break;
                case 21:
                    vscale.value = atof(optarg);
                    if(print_flag) vscale.print(std::cout);
                    opt_used+=2;
                    break;
                case 22:
                    alpha.value = atof(optarg);
                    if(print_flag) alpha.print(std::cout);
                    opt_used+=2;
                    break;
                case 23:
                    lambda.value = atof(optarg);
                    if(print_flag) lambda.print(std::cout);
                    opt_used+=2;
                    break;
                case 24:
                    beta.value = atof(optarg);
                    if(print_flag) beta.print(std::cout);
                    opt_used+=2;
                    break;
                case 25:
                    xi.value = atof(optarg);
                    if(print_flag) xi.print(std::cout);
                    opt_used+=2;
                    break;
                case 26:
                    bhwacc.value = atof(optarg);
                    if(print_flag) bhwacc.print(std::cout);
                    opt_used+=2;
                    break;
                case 27:
                    epsnov.value = atof(optarg);
                    if(print_flag) epsnov.print(std::cout);
                    opt_used+=2;
                    break;
                case 28:
                    eddfac.value = atof(optarg);
                    if(print_flag) eddfac.print(std::cout);
                    opt_used+=2;
                    break;
                case 29:
                    gamma.value = atof(optarg);
                    if(print_flag) gamma.print(std::cout);
                    opt_used+=2;
                    break;
                default:
                    break;
                }
                break;
            case 'z':
                z.value = atof(optarg);
                if(print_flag) z.print(std::cout);
                opt_used+=2;
                break;
            case 'p':
                fname_par = optarg;
                if(print_flag) {
                    std::string fbse_par = fname_par+".mobse"; 
                    FILE* fpar_in;
                    if( (fpar_in = fopen(fbse_par.c_str(),"r")) == NULL) {
                        fprintf(stderr,"Error: Cannot open file %s.\n", fbse_par.c_str());
                        abort();
                    }
                    input_par_store.readAscii(fpar_in);
                    fclose(fpar_in);
                }
                opt_used+=2;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
                input_par_store.mpi_broadcast();
                PS::Comm::barrier();
#endif
                break;
            case 'h':
                if(print_flag){
                    std::cout<<"MOBSE options:"<<std::endl;
                    input_par_store.printHelp(std::cout, 2, 10, 23);
                }
                return -1;
            case '?':
                opt_used +=2;
                break;
            default:
                break;
            }

        if(print_flag) std::cout<<"----- Finish reading input options of MOBSE -----\n";

        return opt_used;
    }    
};

//! MOBSE interface manager
class MOBSEManager{
public:
    double z, zpars[20]; ///> metallicity parameters
    double tscale; ///> time scaling factor from NB to Myr (t[Myr]=t[NB]*tscale)
    double rscale; ///> radius scaling factor from NB to Rsun
    double mscale; ///> mass scaling factor from NB to Msun
    double vscale; ///> velocity scaling factor from NB to km/s
    const double year_to_day; ///> year to day 
    const char* single_type[16]; ///> name of single type from single star
    const char* binary_type[14]; ///> name of binary type return from binary evolv2, notice if it is -1, it indicate the end of record

    MOBSEManager(): z(0.0), zpars{0}, tscale(0.0), rscale(0.0), mscale(0.0), vscale(0.0), year_to_day(3.6525e8),
                  single_type{"LMS", "MS", "HG", "GB", "CHeB", "FAGB", "SAGB", "HeMS", "HeHG", "HeGB", "HeWD", "COWD", "ONWD", "NS", "BH", "SN"},
                  binary_type{"Unset",               //0
                              "Initial",             //1
                              "Type_change",         //2
                              "Start_Roche",         //3
                              "End_Roche",           //4
                              "Contact",             //5
                              "Start_Symbiotic",     //6
                              "End_Symbiotic",       //7
                              "Common_envelope",     //8
                              "Giant",               //9
                              "Coalescence",         //10
                              "Blue_straggler",      //11
                              "No_remain",           //12
                              "Disrupt"              //13
                              } {}
    

    bool checkParams() {
        assert(z>0.0);
        assert(tscale>0.0);
        assert(rscale>0.0);
        assert(mscale>0.0);
        assert(vscale>0.0);
        return true;
    }

    //! dump MOBSE rand constant to file
    void dumpRandConstant(const char* _fname) {
        FILE* fin;
        if( (fin = fopen(_fname,"w")) == NULL) {
            fprintf(stderr,"Error: Cannot open file %s.\n", _fname);
            abort();
        }
        fprintf(fin, "%d %d %d ", value3_.idum, rand3_.idum2, rand3_.iy);
        for (int i=0; i<32; i++) fprintf(fin, "%d ", rand3_.ir[i]);
        fprintf(fin, "\n");
        fclose(fin);
    }

    //! read MOBSE rand constant from file
    void readRandConstant(const char* _fname) {
        FILE* fin;
        if( (fin = fopen(_fname,"r")) == NULL) {
            printf("Not found.\n");
        }
        else {
            int rcount = fscanf(fin, "%d %d %d ", &value3_.idum, &rand3_.idum2, &rand3_.iy);
            for (int i=0; i<32; i++) rcount += fscanf(fin, "%d ", &rand3_.ir[i]);
            if(rcount<35) {
                std::cerr<<"Error: Data reading fails! requiring data number is 35, only obtain "<<rcount<<".\n";
                abort();
            }
            fclose(fin);
        }
    }
        
    //! print terminal Logo
    static void printLogo(std::ostream & fout) {
        fout<<"---------------------------------------\n"
            <<"             ╔╦╗╔═╗╔╗ ╔═╗╔═╗\n"
            <<"             ║║║║ ║╠╩╗╚═╗║╣ \n"
            <<"             ╩ ╩╚═╝╚═╝╚═╝╚═╝\n"
            <<"---------------------------------------"<<std::endl;
            fout<<" Online document: https://mobse-webpage.netlify.app/\n"
            <<std::endl;
    }

    //! print reference to cite
    static void printReference(std::ostream & fout, const int offset=4) {
        for (int i=0; i<offset; i++) fout<<" ";
        fout<<"MOBSE: Giacobbo N., Mapelli M. & Spera M., 2018, MNRAS, 474, 2959\n";
        fout<<"\t \t (Online document: https://mobse-webpage.netlify.app/)\n"
        <<std::endl;
    }
    
    bool isMassTransfer(const int _binary_type) {
        return (_binary_type>=3&&_binary_type<=9);
    }

    //! notice kick priority is higher than others
    bool isKick(const int _binary_type) {
        return (_binary_type==13);
    }

    bool isMerger(const int _binary_type) {
        return (_binary_type>=10&&_binary_type<=12);
    }

    bool isDisrupt(const int _binary_type) {
        return (_binary_type==14);
    }

    //! initial MOBSE global parameters
    void initial(const IOParamsMOBSE& _input, const bool _print_flag=false) {
        // common block
        value1_.neta  = _input.neta.value;
        value1_.bwind = _input.bwind.value;
        value1_.hewind= _input.hewind.value;

        value2_.alpha  = _input.alpha.value;
        value2_.lambda = _input.lambda.value;
        
        value4_.sigma1  = _input.sigma1.value;
        value4_.sigma2  = _input.sigma2.value;
        value4_.mxns  = 1.8;
        if (_input.nsflag.value>0) value4_.mxns = 3.;
        value4_.bhflag = _input.bhflag.value;

        value5_.beta = _input.beta.value;
        value5_.xi   = _input.xi.value;
        value5_.bhwacc = _input.bhwacc.value;
        value5_.epsnov = _input.epsnov.value;
        value5_.eddfac = _input.eddfac.value;
        value5_.gamma  = _input.gamma.value;

        flags_.ceflag = _input.ceflag.value;
        flags_.tflag  = _input.tflag.value;
        //flags_.ifflag = _input.ifflag.value;
        flags_.wdflag = _input.wdflag.value;
        flags_.nsflag = _input.nsflag.value;
        flags_.piflag = _input.piflag.value;

        //flags2_.kmech  = _input.kmech.value;
        //flags2_.ecflag = _input.ecflag.value;

        points_.pts1 = _input.pts1.value;
        points_.pts2 = _input.pts2.value;
        points_.pts3 = _input.pts3.value;

        tscale = _input.tscale.value;
        rscale = _input.rscale.value;
        mscale = _input.mscale.value;
        vscale = _input.vscale.value;

        // Set parameters which depend on the metallicity 
        z = _input.z.value;
        zcnsts_(&z, zpars);
        value3_.idum = (_input.idum.value>0)? -_input.idum.value: _input.idum.value;

        // collision matrix
        instar_();

        if (_print_flag) {
            printconst_();
            std::cout<<"z: "<<z<<" zpars: ";
            for (int i=0;i<20;i++) std::cout<<zpars[i]<<" ";
            std::cout<<std::endl;
        }

    }

    //! get current mass in NB unit
    double getMass(StarParameter& _star) {
        return _star.mt/mscale;
    }
    
    //! get mass loss in NB unit
    double getMassLoss(StarParameterOut& _out) {
        return _out.dm/mscale;
    }

    //! get merger radius in NB unit
    double getMergerRadius(StarParameter& _star) {
        // use stellar radius as merger radius
        return _star.r/rscale;
    }

    //! get evolved Time in NB unit
    double getTime(StarParameter& _star) {
        return _star.tphys/tscale;
    }

    //! get the difference of required finishing time and actually evolved time in NB unit
    double getDTMiss(StarParameterOut& _out) {
        return _out.dtmiss/tscale;
    }

    //! print type change 
    void printTypeChange(std::ostream& _fout, StarParameter& _star, StarParameterOut& _out, const int _width=4) {
        assert(_out.kw0>=0&&_out.kw0<16);
        assert(_star.kw>=0&&_star.kw<16);
        _fout<<" "<<std::setw(_width)<<single_type[_out.kw0]<<" -> "<<std::setw(_width)<<single_type[_star.kw];
    }

    //! print MOBSE event
    void printBinaryEvent(std::ostream& _fout, const BinaryEvent& _bin_event) {
        int nmax = _bin_event.getEventNMax();
        for (int k=0; k<nmax; k++) {
            int type = _bin_event.getType(k);
            if(type>0) {
                _fout<<std::setw(16)<<binary_type[type];
                _bin_event.print(_fout, k);
                _fout<<std::endl;
            }
            else if(type<0) break;
        }
    }

    //! print MOBSE event one
    void printBinaryEventOne(std::ostream& _fout, const BinaryEvent& _bin_event, const int k) {
        int type = _bin_event.getType(k);
        assert(type>=0&&type<14);
        _fout<<std::setw(16)<<binary_type[type];
        if (k==0) _bin_event.print(_fout, _bin_event.getEventIndexInit());
        else _bin_event.print(_fout, k-1);
        _bin_event.print(_fout, k);
    }

    //! print BSE event one
    void printBinaryEventColumnOne(std::ostream& _fout, const BinaryEvent& _bin_event, const int k, const int _width=20) {
        int type = _bin_event.getType(k);
        assert(type>=0&&type<14);
        _fout<<std::setw(16)<<binary_type[type]
             <<std::setw(_width)<<type;
        if (k==0) _bin_event.printColumn(_fout, _bin_event.getEventIndexInit(), _width);
        else _bin_event.printColumn(_fout, k-1, _width);
        _bin_event.printColumn(_fout, k, _width);
    }

    //! get velocity change in NB unit
    /*!
      @param[in] _dv: 3-D array to record velocity change
      \return value of velocity change
     */
    double getVelocityChange(double* dv, StarParameterOut& _out) {
        for (int k=0; k<3; k++) dv[k] = _out.vkick[k]/vscale;
        return _out.vkick[3]/vscale;
    }

    //! call evolv1 for single star
    /*!
      @param[in,out] _star: star parameter
      @param[out] _out: output parameter from evolv1
      @param[in] _dt_nb: physical time step in NB unit to evolve
      \return event flag: -1: error, 0: normal, 1: type change, 2: velocity kick
     */
    int evolveStar(StarParameter& _star, StarParameterOut& _out, const double _dt_nb) {
        double tphysf = _dt_nb*tscale + _star.tphys;
        double dtp=tphysf*100.0+1000.0;
        _out.dm = _star.mt;
        _out.kw0 = _star.kw;
        int kw = _star.kw;
        evolv1_(&kw, &_star.m0, &_star.mt, &_star.r, 
                &_star.lum, &_star.mc, &_star.rc, &_out.menv, &_out.renv, 
                &_star.ospin, &_star.epoch, 
                &_out.tm, &_star.tphys, &tphysf, &dtp, &z, zpars, _out.vkick);
                _star.kw = kw;
        _out.dm = _star.mt - _out.dm;
        _out.dtmiss = tphysf - _star.tphys;

        if (_star.kw<0) {
            _star.kw = -_star.kw;
            return -1; // error 
        }
        else if (_out.vkick[3]>0) return 2; // kick
        else if (_star.kw!=_out.kw0) return 1; // kw change
        else return 0;
    }

    //! call evolv2 for a binary
    /*!
      @param[in,out] _star1: star parameter of first
      @param[in,out] _star2: star parameter of second
      @param[out] _out1: output parameter of first from evolv2
      @param[out] _out2: output parameter of second from evolv2
      @param[out] _bse_event: MOBSE event record (bpp array)
      @param[in] _semi: semi-major axis, only used to record initial semi
      @param[in,out] _period: period of binary in NB unit, used for MOBSE
      @param[in,out] _ecc: eccentricity of binary, used for MOBSE
      @param[in] _dt_nb: physical time step in Myr to evolve
      \return error flag: -1: error, 0: normal
     */
    int evolveBinary(StarParameter& _star1, StarParameter& _star2, StarParameterOut& _out1, StarParameterOut& _out2, 
                     double& _semi, double& _period, double& _ecc, BinaryEvent& _bse_event, const double _dt_nb) {
        double tphys = std::max(_star1.tphys, _star2.tphys);
        double tphysf = _dt_nb*tscale + tphys;
        double dtp=tphysf*100.0+1000.0;
        double period_days = _period*tscale*year_to_day;
        double semi_rsun = _semi*rscale;
        // in case two component have different tphys, evolve to the same time first
        int event_flag = 0 ;
        _out1.dm = _star1.mt;
        _out2.dm = _star2.mt;

        // backup initial state
        _bse_event.recordInitial(_star1, _star2, semi_rsun, _ecc);

        if (_star1.tphys<tphys) event_flag = evolveStar(_star1, _out1, tphys);
        if (event_flag<0) return event_flag;
        if (_star2.tphys<tphys) event_flag = evolveStar(_star2, _out2, tphys);
        if (event_flag<0) return event_flag;
        
        int kw[2];
        double m0[2],mt[2],r[2],lum[2],mc[2],rc[2],menv[2],renv[2],ospin[2],epoch[2],tm[2],vkick[8];
        for (int k =0; k<4; k++) {
            vkick[k]  = _out1.vkick[k];
            vkick[k+4]= _out2.vkick[k];
        }

        kw[0] = _star1.kw;
        m0[0] = _star1.m0;
        mt[0] = _star1.mt;
        r[0]  = _star1.r;
        mc[0] = _star1.mc;
        rc[0] = _star1.rc;
        ospin[0] = _star1.ospin;
        epoch[0] = _star1.epoch;
 
        kw[1] = _star2.kw;
        m0[1] = _star2.m0;
        mt[1] = _star2.mt;
        r[1]  = _star2.r;
        mc[1] = _star2.mc;
        rc[1] = _star2.rc;
        ospin[1] = _star2.ospin;
        epoch[1] = _star2.epoch;
        
        evolv2_(kw, m0, mt, r, lum, mc, rc, menv, renv, ospin, epoch, tm, &tphys, &tphysf, &dtp, &z, zpars, &period_days, &_ecc, _bse_event.record[0], vkick);
        _period = period_days/year_to_day/tscale;

        _star1.kw = kw[0];
        _star1.m0 = m0[0];
        _star1.mt = mt[0];
        _star1.r  = r[0];
        _star1.mc = mc[0];
        _star1.rc = rc[0];
        _star1.ospin  = ospin[0];
        _star1.epoch  = epoch[0];
        _star1.tphys  = tphys;
        _star1.lum    = lum[0];

        _star2.kw = kw[1];
        _star2.m0 = m0[1];
        _star2.mt = mt[1];
        _star2.r  = r[1];
        _star2.mc = mc[1];
        _star2.rc = rc[1];
        _star2.ospin  = ospin[1];
        _star2.epoch  = epoch[1];
        _star2.tphys  = tphys;
        _star2.lum    = lum[1];

        _out1.menv = menv[0];
        _out1.renv = renv[0];
        _out1.tm = tm[0];
        _out1.dm = _star1.mt - _out1.dm;
        _out1.dtmiss = tphysf - _star1.tphys;

        _out2.menv = menv[1];
        _out2.renv = renv[1];
        _out2.tm = tm[1];
        _out2.dm = _star2.mt - _out2.dm;
        _out2.dtmiss = tphysf - _star2.tphys;

        for (int k=0; k<4; k++) _out1.vkick[k]=vkick[k];
        for (int k=0; k<4; k++) _out2.vkick[k]=vkick[k+4];

        if (kw[0]<0||kw[1]<0||(_star1.mt<0&&_star1.kw==15)||(_star2.mt<0&&_star2.kw==15)) {
            kw[0] = abs(kw[0]);
            kw[1] = abs(kw[1]);
            return -1; // error case
        }
        //else if (vkick[3]>0||vkick[7]>0) return 3; // kick
        //else if (isDisrupt(_binary_type)) return 4; // disrupt without kick
        //else if (isMerger(_binary_type)) return 5; // Merger
        //else if (isMassTransfer(_binary_type)) return 2; // orbit change
        //else if (_binary_type>0) return 1; // type change
        else return 0;
    }

    //! merge two star using MOBSE mix function, star 2 will becomes zero mass
    void merge(StarParameter& _star1, StarParameter& _star2) {
        double m0[2],mt[2],age[2];
        int kw[2],krol[2];

        kw[0] = _star1.kw;
        m0[0] = _star1.m0;
        mt[0] = _star1.mt;
        age[0]= _star1.tphys-_star1.epoch;

        kw[1] = _star2.kw;
        m0[1] = _star2.m0;
        mt[1] = _star2.mt;
        age[1]= _star2.tphys-_star2.epoch;

        //! This is necessary for the collision matrix
        krol[0] = _star1.kw;
        krol[1] = _star2.kw;

        mix_(m0,mt,age,kw,zpars,krol);

        _star1.kw = kw[0];
        _star1.m0 = m0[0];
        _star1.mt = mt[0];
        _star1.epoch = _star1.tphys - age[0];

        _star2.kw = kw[1];
        _star2.m0 = m0[1];
        _star2.mt = mt[1];
        _star2.epoch = _star2.tphys - age[1];
    }

    //! get next time step to check in Myr
    double getTimeStep(StarParameter& _star) {
        if (_star.kw==15) return 1.0e30/tscale; // give very large value to avoid evolve

        double tm, tn, tscls[20], lums[10], gb[10], dtm, dtr;
        
        // obtain star parameters
        int kw = _star.kw;
        star_(&kw, &_star.m0, &_star.mt, &tm, &tn, tscls, lums, gb, zpars);

        // get next step
        double age = _star.tphys-_star.epoch;
        deltat_(&kw, &age, &tm, &tn, tscls, &dtm, &dtr);
        _star.kw = kw;

        //assert(dtr>0.0);
        //assert(dtm>0.0);

        return std::min(dtr, dtm)/tscale;
    }
    
};
