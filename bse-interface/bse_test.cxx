#include <iostream>
#include <fstream>
#include <getopt.h>
#include <cmath>
#include <vector>
#include <string>
#include "bse_interface.h"
#include "../src/io.hpp"
#include "../parallel-random/rand_interface.hpp"
#define ASSERT assert
#include "../src/gw_kick.hpp"

const int WRITE_WIDTH=23;
const int WRITE_PRECISION=14;

struct BinaryBase{
    double m1,m2,semi,period,ecc,period0,ecc0;
    int kw1,kw2;
    double tphys;
    StarParameter star[2];
    StarParameterOut out[2];
    BinaryEvent bse_event;

    void readAscii(FILE* fp) {
        int rcount=fscanf(fp, "%lf %lf %d %d %lf %lf %lf",
                          &m1, &m2, &kw1, &kw2, &period, &ecc, &tphys);
        if(rcount<4) {
            std::cerr<<"Error: Data reading fails! requiring data number is 4, only obtain "<<rcount<<".\n";
            abort();
        }
    }
};

int main(int argc, char** argv){

    int arg_label;
    int width=13;
    int n=5000;
    double m_min=0.08, m_max=150.0;
    double time=100.0;
    double time0=0.0;
    double dtmin=1.0;
    std::vector<StarParameter> star;
    std::vector<double> mass0;
    std::vector<BinaryBase> bin;
    std::vector<BinaryBase> hyb;
    std::string fprint_name;
    std::string fbin_name;
    std::string fsin_name;
    std::string fhyb_name;
    bool read_mass_flag = false;
    bool always_output_flag = false;
    std::string bse_prefix = BSEManager::getBSEOutputFilenameSuffix();
    std::string sse_prefix = BSEManager::getSSEOutputFilenameSuffix();

    auto printHelp= [&]() {
        std::cout<<"The tool to evolve single stars or binaries using "<<BSEManager::getBSEName()<<std::endl;
#ifdef MOBSE
        BSEManager::printLogo(std::cout);
#endif
        BSEManager::printReference(std::cout);
        std::cout<<"Usage: petar"<<bse_prefix<<" [options] [initial mass of stars, can be multiple values]\n"
                 <<"       * If no initial mass or no single/binary/hyperbolic table (-s, -b or -m) is provided,\n"
                 <<"         N single stars (-n) with equal mass interal in Log scale will be evolved.\n"
                 <<"       * When a single/binary/hyperbolic table is provided (-s, -b or -m), times and types of stars can be set individually.\n"
                 <<"       * When only one star or binary is provided, the single and binary types change are printed;\n"
                 <<"         if -o is used, data of every step is printed.\n"
                 <<"       * The default unit set is: Msun, Myr. If the input data have different units (referred as [IN]),\n"
                 <<"         please modify the scaling fators.\n"
                 <<"Options:\n"
                 <<"    -n [I]: number of stars when evolve an IMF ("<<n<<")\n"
                 <<"    -i [D]: start time ("<<time0<<") [IN]\n"
                 <<"    -t [D]: finish time ("<<time<<") [IN]\n"
                 <<"        --mmin [D]: mimimum mass ("<<m_min<<") [M*]\n"
                 <<"        --mmax [D]: maximum mass ("<<m_max<<") [M*]\n"
                 <<"    -d [D]: minimum time step to call SSE/BSE evolution functions ("<<dtmin<<")[IN]\n"
                 <<"    -s [S]: a file of single table: First line: number of single (unit:IN); After: mass, type, time per line\n"
                 <<"    -b [S]: a file of binary table: First line: number of binary (unit:IN); After: m1, m2, type1, type2, period, ecc, time per line\n"
                 <<"    -m [S]: a file of hyperbolic orbit table for checking merger: First line: number of orbit (unit:IN); After: m1, m2, type1, type2, semi, ecc, time\n"
                 <<"    -w [I]: print column width ("<<width<<")\n"
                 <<"    -o    : print evolution data every time when stellar evolution function is called (maximum time step by -d)\n"
                 <<"            * If this option is not used, only output data when single or binary types change.\n"
                 <<"    -f [S]: the prefix of data file that save the stellar evolution data, if not given, this file is not generated\n"
                 <<"            * The files [prefix].["<<sse_prefix<<"/"<<bse_prefix<<"].type_change store single or binary types change; if -o is used, save data every step\n"
                 <<"            * The files [prefix].["<<sse_prefix<<"/"<<bse_prefix<<"].sn_kick store supernovae kick events\n"
                 <<"            * The files [prefix].["<<bse_prefix<<"].gw_kick store gravitational wave kick events\n"
                 <<"    -h    : help\n";
    };

    IOParamsBSE bse_io;
    IOParamsRand rand_io;
    opterr = 0;

    // reset optind
    static int long_flag=-1;
    static struct option long_options[] = {
        {"mmin", required_argument, &long_flag, 0},
        {"mmax", required_argument, &long_flag, 1},
        {0,0,0,0}
    };

    int opt_used = 0;
    optind=0;
    int option_index;
    bool help_flag=false;
    while ((arg_label = getopt_long(argc, argv, "-n:i:m:t:w:d:s:b:of:h", long_options, &option_index)) != -1)
        switch (arg_label) {
        case 0:
            switch (long_flag) {
            case 0:
                m_min = atof(optarg);
                std::cout<<"min mass: "<<m_min<<std::endl;
                opt_used+=2;
                break;
            case 1:
                m_max = atof(optarg);
                std::cout<<"max mass: "<<m_max<<std::endl;
                opt_used+=2;
                break;
            default:
                break;
            }
            break;
        case 'n':
            n = atof(optarg);
            std::cout<<"N: "<<n<<std::endl;
            opt_used+=2;
            break;
        case 'i':
            time0 = atof(optarg);
            std::cout<<"start time: "<<time0<<std::endl;
            opt_used+=2;
            break;
        case 't':
            time = atof(optarg);
            std::cout<<"finish time: "<<time<<std::endl;
            opt_used+=2;
            break;
        case 'w':
            width = atoi(optarg);
            std::cout<<"print width: "<<width<<std::endl;
            opt_used+=2;
            break;
        case 'd':
            dtmin = atof(optarg);
            std::cout<<"minimum time step "<<dtmin<<std::endl;
            opt_used+=2;
            break;
        case 's':
            fsin_name = optarg;
            read_mass_flag = true;
            std::cout<<"Single data file "<<fsin_name<<std::endl;
            opt_used+=2;
            break;
        case 'b':
            fbin_name = optarg;
            read_mass_flag = true;
            std::cout<<"Binary data file "<<fbin_name<<std::endl;
            opt_used+=2;
            break;
        case 'm':
            fhyb_name= optarg;
            read_mass_flag = true;
            std::cout<<"hyperbolic data file "<<fhyb_name<<std::endl;
            opt_used+=2;
            break;
        case 'o':
            always_output_flag = true;
            std::cout<<"Always output data every call of bse function"<<std::endl;
            opt_used++;
            break;
        case 'f':
            fprint_name = optarg;
            std::cout<<"Output data file "<<fprint_name<<std::endl;
            opt_used+=2;
            break;
        case 'h':
            printHelp();
            help_flag=true;
            break;
        case '?':
            opt_used +=2;
            break;
        default:
            break;
        }        

    bse_io.print_flag = true;
    bse_io.read(argc,argv);
    rand_io.print_flag = true;
    rand_io.read(argc,argv);

    if (help_flag) return 0;

    BSEManager bse_manager;
    RandomManager rand_manager;
    GWKick gw_kick;
    bse_manager.initial(bse_io, true);
    rand_manager.initialAll(rand_io);
    rand_manager.printRandSeeds(std::cout);
    gw_kick.vscale = bse_manager.vscale;
    gw_kick.speed_of_light = bse_manager.getSpeedOfLight()*gw_kick.vscale;
    assert(bse_manager.checkParams());
    assert(gw_kick.checkParams());

    // argc is 1 no input is given
    opt_used++;
    if (opt_used<argc) {
        // read initial mass list
        while (opt_used<argc) {
            mass0.push_back(atof(argv[opt_used++]));
        }
        n = mass0.size();
        star.resize(n);
        for (int i=0; i<n; i++) {
            star[i].initial(mass0[i]*bse_manager.mscale);
        }
        read_mass_flag = true;
    }

    // if no mass is read, use a mass range
    if (!read_mass_flag) {
        double dm_factor = exp((log(m_max) - log(m_min))/n);

        mass0.resize(n);
        star.resize(n);

        mass0[0] = m_min;
        star[0].initial(m_min);
        for (int i=1; i<n; i++) {
            mass0[i] = mass0[i-1]*dm_factor;
            star[i].initial(mass0[i]);
        }
    }

    if (fsin_name!="") {
        FILE* fsin;
        if( (fsin = fopen(fsin_name.c_str(),"r")) == NULL) {
            fprintf(stderr,"Error: Cannot open file %s.\n", fsin_name.c_str());
            abort();
        }
        int n0 = mass0.size();
        int nadd;
        int rcount=fscanf(fsin, "%d", &nadd);
        if(rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            abort();
        }
        n = n0+nadd;
        mass0.resize(n);
        star.resize(n);
        for (int k=n0; k<n; k++) {
            int type;
            double tphys;
            rcount=fscanf(fsin, "%lf %d %lf", &mass0[k], &type, &tphys);
            if(rcount<3) {
                std::cerr<<"Error: Data reading fails! requiring data number is 3, only obtain "<<rcount<<".\n";
                abort();
            }
            star[k].initial(mass0[k]*bse_manager.mscale, type, 0.0, tphys*bse_manager.tscale);
        }
    }

    if (fbin_name!="") {
        FILE* fbin;
        if( (fbin = fopen(fbin_name.c_str(),"r")) == NULL) {
            fprintf(stderr,"Error: Cannot open file %s.\n", fbin_name.c_str());
            abort();
        }
        int nb;
        int rcount=fscanf(fbin, "%d", &nb);
        if(rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            abort();
        }
        for (int k=0; k<nb; k++) {
            BinaryBase bink;
            bink.readAscii(fbin);
            bink.star[0].initial(bink.m1*bse_manager.mscale);
            bink.star[1].initial(bink.m2*bse_manager.mscale);
            bink.tphys *= bse_manager.tscale;
            bink.star[0].tphys = bink.tphys;
            bink.star[1].tphys = bink.tphys;
            bink.star[0].kw = bink.kw1;
            bink.star[1].kw = bink.kw2;
            bin.push_back(bink);
            // initial
        }
    }

    if (fhyb_name!="") {
        FILE* fhyb;
        if( (fhyb = fopen(fhyb_name.c_str(),"r")) == NULL) {
            fprintf(stderr,"Error: Cannot open file %s.\n", fhyb_name.c_str());
            abort();
        }
        int nb;
        int rcount=fscanf(fhyb, "%d", &nb);
        if(rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            abort();
        }
        double pc_to_rsun = 44334448.006896;
        for (int k=0; k<nb; k++) {
            BinaryBase hybk;
            hybk.readAscii(fhyb);
            hybk.star[0].initial(hybk.m1*bse_manager.mscale);
            hybk.star[1].initial(hybk.m2*bse_manager.mscale);
            hybk.tphys *= bse_manager.tscale;
            hybk.star[0].tphys = hybk.tphys;
            hybk.star[1].tphys = hybk.tphys;
            hybk.star[0].kw = hybk.kw1;
            hybk.star[1].kw = hybk.kw2;
            hybk.semi = hybk.period*pc_to_rsun;
            hybk.period = 0.0;
            hybk.period0 = 0.0;
            hybk.ecc0 = hybk.ecc;
            hyb.push_back(hybk);
            // initial
        }
    }

    auto printBinaryTitle=[&](std::ostream & _fout) {
        _fout<<std::setw(width)<<"mass0_1[M*]"
             <<std::setw(width)<<"mass0_2[M*]"
             <<std::setw(width)<<"P0[days]"
             <<std::setw(width)<<"ecc0"
             <<std::setw(width)<<"P[days]"
             <<std::setw(width)<<"ecc";
        StarParameter::printColumnTitle(_fout, width);
        StarParameterOut::printColumnTitle(_fout, width);
        StarParameter::printColumnTitle(_fout, width);
        StarParameterOut::printColumnTitle(_fout, width);
        _fout<<std::endl;
    };

    auto printBinary=[&](std::ostream & _fout, BinaryBase& _bin){
        _fout<<std::setw(width)<<_bin.m1*bse_manager.mscale;
        _fout<<std::setw(width)<<_bin.m2*bse_manager.mscale;
        _fout<<std::setw(width)<<_bin.period0*bse_manager.tscale*3.6524e8;
        _fout<<std::setw(width)<<_bin.ecc0;
        _fout<<std::setw(width)<<_bin.period*bse_manager.tscale*3.6524e8;
        _fout<<std::setw(width)<<_bin.ecc;
        for (int k=0; k<2; k++) {
            _bin.star[k].printColumn(_fout, width);
            _bin.out[k].printColumn(_fout, width);
        }
        _fout<<std::endl;
    };

    auto printSingleColumnTitle=[&](std::ostream & _fout) {
        _fout<<std::setw(width)<<"Mass_init[Msun]";
        StarParameter::printColumnTitle(_fout, width);
        StarParameterOut::printColumnTitle(_fout, width);
        _fout<<std::endl;
    };

    auto printSingleColumn=[&](std::ostream & _fout, double& _mass0, StarParameter& _star, StarParameterOut& _out) {
        _fout<<std::setw(width)<<_mass0*bse_manager.mscale;
        _star.printColumn(_fout, width);
        _out.printColumn(_fout, width);
        _fout<<std::endl;
    };

    auto printSingle=[&](std::ostream & _fout, StarParameter& _star, StarParameterOut& _out) {
        _star.print(_fout);
        _out.print(_fout);
        _fout<<std::endl;
    };

    bool output_flag=false;
    if (fprint_name!="") output_flag = true;

    // first check whether binary exist
    if (bin.size()>0) {

        // open output file
        std::ofstream fout_bse_type, fout_bse_sn, fout_gw_kick;
        if (output_flag) {
            std::string bse_suffix=BSEManager::getBSEOutputFilenameSuffix();
            fout_bse_type.open((fprint_name+bse_suffix+std::string(".type_change")).c_str(), std::ofstream::out);
            fout_bse_sn.open((fprint_name+bse_suffix+std::string(".sn_kick")).c_str(), std::ofstream::out);
            fout_gw_kick.open((fprint_name+bse_suffix+std::string(".gw_kick")).c_str(), std::ofstream::out);
            fout_bse_type<<std::setprecision(WRITE_PRECISION);
            fout_bse_sn<<std::setprecision(WRITE_PRECISION);
            fout_gw_kick<<std::setprecision(WRITE_PRECISION);
        }

        int nbin = bin.size();
#pragma omp parallel for schedule(dynamic)
        for (int i=0; i<nbin; i++) {
            // output file

            bool kick_print_flag[2]={false,false};
            // evolve
            double tend = time*bse_manager.tscale;
            bin[i].period0 = bin[i].period;
            bin[i].ecc0 = bin[i].ecc;
            int bin_type_last=0;
            bool initial_flag = true;
            while (bse_manager.getTime(bin[i].star[0])<tend) {
                // time step
                //double dt1 = bse_manager.getTimeStep(bin[i].star[0]);
                //double dt2 = bse_manager.getTimeStep(bin[i].star[1]);
                //double dt = std::min(dt1,dt2);
                double dt = bse_manager.getTimeStepBinary(bin[i].star[0],bin[i].star[1],bin[i].semi,bin[i].ecc,bin_type_last);
                dt = std::max(dt,dtmin);
                dt = std::min(tend-bse_manager.getTime(bin[i].star[0]), dt);
                if (initial_flag) {
                    dt = 0;
                    initial_flag = false;
                }
                double mtot = bin[i].star[0].mt + bin[i].star[1].mt;
                double period_myr = bin[i].period*bse_manager.tscale;
                double G= 0.00449830997959438;
                const double PI = 4.0*atan(1.0);
                double pc_to_rsun = 44334448.006896;
                bin[i].semi = std::pow(period_myr*period_myr*G*mtot/(4*PI*PI),1.0/3.0)*pc_to_rsun;

                StarParameter p1_star_bk = bin[i].star[0];
                StarParameter p2_star_bk = bin[i].star[1];
                
                // evolve function
                int event_flag=bse_manager.evolveBinary(bin[i].star[0],bin[i].star[1],bin[i].out[0],bin[i].out[1],bin[i].semi,bin[i].period,bin[i].ecc,bin[i].bse_event, bin_type_last, dt);
                    

                // error
                if (event_flag<0) {
                    std::cerr<<"Error! ";
                    std::cerr<<" BID="<<i+1<<" ";
                    std::cerr<<" semi[R*]: "
                             <<bin[i].semi
                             <<" ecc: "<<bin[i].ecc
                             <<" period[days]: "<<bin[i].period;
                    std::cerr<<" Init: Star1: ";
                    p1_star_bk.print(std::cerr);
                    std::cerr<<" Star2: ";
                    p2_star_bk.print(std::cerr);
                    std::cerr<<" final: Star1: ";
                    bin[i].star[0].print(std::cerr);
                    bin[i].out[0].print(std::cerr);
                    std::cerr<<" Star2: ";
                    bin[i].star[1].print(std::cerr);
                    bin[i].out[1].print(std::cerr);
                    std::cerr<<std::endl;
                    abort();
                }

                int nmax = bin[i].bse_event.getEventNMax();
                int bin_type_init = bin[i].bse_event.getType(bin[i].bse_event.getEventIndexInit());
                bin_type_last = bin_type_init;
                int merger_event_index = -1; // record binary event index for merger, if no merger, is -1
                for (int k=0; k<nmax; k++) {
                    int binary_type = bin[i].bse_event.getType(k);
                    if (binary_type>0) {
                        bool first_event = (k==0);
                        if ( (first_event&&bin_type_init!=binary_type) || !first_event || always_output_flag) {
                            //if (!(bin_type_last==11&&binary_type==3)) {// avoid repeating printing Start Roche and BSS
                            if (output_flag) {
#pragma omp critical
                                {
                                    bse_manager.printBinaryEventColumnOne(fout_bse_type, bin[i].bse_event, k, WRITE_WIDTH, false);
                                    fout_bse_type<<std::setw(WRITE_WIDTH)<<2*i+1
                                                 <<std::setw(WRITE_WIDTH)<<2*i+2
                                                 <<std::setw(WRITE_WIDTH)<<0
                                                 <<std::setw(WRITE_WIDTH)<<0
                                                 <<std::endl;
                                }
                            }
                            // print data
                            if (nbin==1) {
                                //std::cout<<" BID="<<i+1<<" index="<<k<<" "<<" dt="<<dt;
                                bse_manager.printBinaryEventOne(std::cout, bin[i].bse_event, k);
                                std::cout<<" dt="<<dt<<std::endl;
                            }
                        }
                        bin_type_last = binary_type;
                        if (bse_manager.isMerger(binary_type)) {
                            if (merger_event_index==-1) merger_event_index = i; // avoid save index twice
                        }
                    }
                    else if (binary_type<0) {
                        if (always_output_flag) {
                            if (k==0) bin[i].bse_event.print(std::cout,bin[i].bse_event.getEventIndexInit());
                            else bin[i].bse_event.print(std::cout,k-1);
                            std::cout<<std::endl;
                        }
                        break;
                    }
                }

                // check SN event
                for (int k=0; k<2; k++) {
                    double dv[4];
                    dv[3] = bse_manager.getVelocityChange(dv,bin[i].out[k]);
                    if (dv[3]>0&&!kick_print_flag[k]) {
#pragma omp critical 
                        {
                            fout_bse_sn<<std::setw(WRITE_WIDTH)<<2*i+1
                                       <<std::setw(WRITE_WIDTH)<<2*i+2
                                       <<std::setw(WRITE_WIDTH)<<k+1
                                       <<std::setw(WRITE_WIDTH)<<dv[3]*bse_manager.vscale;
                            bin[i].star[k].printColumn(fout_bse_sn, WRITE_WIDTH);
                            fout_bse_sn<<std::endl;
                        }
                        if (nbin==1) {
                            std::cout<<"SN kick, BID="<<i+1<<" Member="<<k+1<<" vkick[IN]="<<dv[3]<<" ";
                            bin[i].star[k].print(std::cout);
                            std::cout<<std::endl;
                        }
                        kick_print_flag[k]=true;
                    }
                }
                
                // check GW merger
                if (merger_event_index>=0) {
                    Float m1 = bse_manager.getMass(bin[i].star[0]);
                    Float m2 = bse_manager.getMass(bin[i].star[1]);
                    ASSERT(m1==0.0 || m2==0.0);
                    ASSERT(!(m1==0.0 && m2==0.0));
    
                    int type1, type2;    
                    Float q, m1_pre, m2_pre, mf, semi_pre, ecc_pre; // mass ratio
                    auto& bin_event = bin[i].bse_event;
                        
                    if (merger_event_index==0) {
                        type1 = p1_star_bk.kw;
                        type2 = p2_star_bk.kw;
                        m1_pre = bse_manager.getMass(p1_star_bk);
                        m2_pre = bse_manager.getMass(p2_star_bk);
                        q = m1_pre/m2_pre;
                        if (q>1) q = 1/q;
                        semi_pre = bin[i].semi;
                        ecc_pre = bin[i].ecc;
                    }
                    else {
                        type1 = bin_event.getType1(merger_event_index-1);
                        type2 = bin_event.getType2(merger_event_index-1);
                        m1_pre = bin_event.getMass1(merger_event_index-1);
                        m2_pre = bin_event.getMass2(merger_event_index-1);
                        ASSERT(m1_pre>0 && m2_pre>0);
                        q = bin_event.getMassRatio(merger_event_index-1);
                        semi_pre = bin_event.getSemi(merger_event_index-1);
                        ecc_pre = bin_event.getEcc(merger_event_index-1);    
                    }
                    Float tmerge = bin_event.getTime(merger_event_index);
                    // GW with NS or BH binaries
                    if (type1>=13 && type1<=14 && type2>=13 && type2<=14) {
                        
                        std::array<Float, 3> chi1 = gw_kick.uniformPointsInsideSphere(0.8);
                        std::array<Float, 3> chi2 = gw_kick.uniformPointsInsideSphere(0.8);
                        
                        // save kick to output                    
                        int kick_index;
                        if (m1==0.0) {
                            mf = m2;
                            kick_index = 1;    
                        }
                        else {
                            mf = m1;
                            kick_index = 0;
                        }
                        Float L[3] = {0.0, 0.0, 1.0};
                        Float dr[3] = {1.0, 0.0, 0.0};
                        Float vkick[3];
                        gw_kick.calcKickVel(vkick, chi1.data(), chi2.data(), L, dr, q);
#pragma omp critical      
                        {
                            fout_gw_kick<<std::setw(WRITE_WIDTH)<<2*i+1
                                        <<std::setw(WRITE_WIDTH)<<2*i+2
                                        <<std::setw(WRITE_WIDTH)<<kick_index+1
                                        <<std::setw(WRITE_WIDTH)<<vkick[0]*gw_kick.vscale
                                        <<std::setw(WRITE_WIDTH)<<vkick[1]*gw_kick.vscale
                                        <<std::setw(WRITE_WIDTH)<<vkick[2]*gw_kick.vscale
                                        <<std::setw(WRITE_WIDTH)<<tmerge
                                        <<std::setw(WRITE_WIDTH)<<m1_pre
                                        <<std::setw(WRITE_WIDTH)<<m2_pre
                                        <<std::setw(WRITE_WIDTH)<<mf
                                        <<std::setw(WRITE_WIDTH)<<semi_pre
                                        <<std::setw(WRITE_WIDTH)<<ecc_pre
                                        <<std::setw(WRITE_WIDTH)<<chi1[0]
                                        <<std::setw(WRITE_WIDTH)<<chi1[1]
                                        <<std::setw(WRITE_WIDTH)<<chi1[2]
                                        <<std::setw(WRITE_WIDTH)<<chi2[0]
                                        <<std::setw(WRITE_WIDTH)<<chi2[1]
                                        <<std::setw(WRITE_WIDTH)<<chi2[2]
                                        <<std::setw(WRITE_WIDTH)<<L[0]
                                        <<std::setw(WRITE_WIDTH)<<L[1]
                                        <<std::setw(WRITE_WIDTH)<<L[2]
                                        <<std::setw(WRITE_WIDTH)<<dr[0]
                                        <<std::setw(WRITE_WIDTH)<<dr[1]
                                        <<std::setw(WRITE_WIDTH)<<dr[2]
                                        <<std::endl;
                        }
                        if (nbin==1) {
                            std::cout<<"GW kick, time[Myr] = "<<tmerge<<" m1[M*] ="<<m1_pre<<" m2[M*] ="<<m2_pre
                                    <<" chi1 ="<<chi1[0]<<" "<<chi1[1]<<" "<<chi1[2]
                                    <<" chi2 ="<<chi2[0]<<" "<<chi2[1]<<" "<<chi2[2]
                                    <<" vkick[km/s]="<<vkick[0]<<" "<<vkick[1]<<" "<<vkick[2]
                                    <<std::endl;
                        }   
                    }   
                }
                else {
                    assert(bse_manager.getMass(bin[i].star[0])>0 && bse_manager.getMass(bin[i].star[1])>0);     
                }

                double dt_miss = bse_manager.getDTMiss(bin[i].out[0]);

                if (dt_miss!=0.0&&bin[i].star[0].kw>=15&&bin[i].star[1].kw>=15) break;
                
                if (bin[i].star[0].kw>=15) {
                    mass0.push_back(bin[i].star[1].m0/bse_manager.mscale);
                    star.push_back(bin[i].star[1]);
                    break;
                }
                if (bin[i].star[1].kw>=15) {
                    mass0.push_back(bin[i].star[0].m0/bse_manager.mscale);
                    star.push_back(bin[i].star[0]);
                    break;
                }
                 
                if (bse_manager.isDisrupt(bin_type_last)) {
                    mass0.push_back(bin[i].star[0].m0/bse_manager.mscale);
                    star.push_back(bin[i].star[0]);
                    mass0.push_back(bin[i].star[1].m0/bse_manager.mscale);
                    star.push_back(bin[i].star[1]);
                    break;
                }
            }
        }

        printBinaryTitle(std::cout);
        for (int i=0; i<nbin; i++) printBinary(std::cout, bin[i]);

        if (output_flag) {
            fout_bse_type.close();
            fout_bse_sn.close();
            fout_gw_kick.close();
        }
    }

    if (star.size()>0) {

        // output file open
        std::ofstream fout_sse_type, fout_sse_sn;
        if (output_flag) {
            std::string sse_suffix=BSEManager::getSSEOutputFilenameSuffix();
            fout_sse_type.open((fprint_name+sse_suffix+std::string(".type_change")).c_str(), std::ofstream::out);
            fout_sse_sn.open((fprint_name+sse_suffix+std::string(".sn_kick")).c_str(), std::ofstream::out);
            fout_sse_type<<std::setprecision(WRITE_PRECISION);
            fout_sse_sn<<std::setprecision(WRITE_PRECISION);
        }

        StarParameterOut output[star.size()];

#pragma omp parallel for schedule(dynamic)
        for (size_t i=0; i<star.size(); i++) {
            //int error_flag = bse_manager.evolveStar(star[i],output[i],time);
            bool kick_print_flag=false;
            double tend = time*bse_manager.tscale;
            bool initial_flag = true;
            while (bse_manager.getTime(star[i])<tend) {
                double dt = std::max(bse_manager.getTimeStepStar(star[i]),dtmin);
                dt = std::min(tend-bse_manager.getTime(star[i]), dt);
                StarParameter star_bk = star[i];
                if (initial_flag) {
                    dt = 0;
                    initial_flag = false;
                }
                int event_flag=bse_manager.evolveStar(star[i],output[i],dt);

                // error 
                if (event_flag<0) {
                    std::cerr<<"SSE Error: ID= "<<i+1<<" mass0[IN]="<<mass0[i]<<" ";
                    star[i].print(std::cerr);
                    output[i].print(std::cerr);
                    std::cerr<<std::endl;
                    abort();
                }

                double dv[4];
                dv[3] = bse_manager.getVelocityChange(dv, output[i]);

                if ((event_flag>0 && event_flag<2) || (event_flag==2 && !kick_print_flag) || always_output_flag) {
                    if (output_flag) {
#pragma omp critical 
                        {
                            fout_sse_type<<std::setw(WRITE_WIDTH)<<i+1;
                            star_bk.printColumn(fout_sse_type, WRITE_WIDTH);
                            star[i].printColumn(fout_sse_type, WRITE_WIDTH);
                            fout_sse_type<<std::endl;
                        }
                    }
                    if (star.size()==1) printSingle(std::cout,  star[i], output[i]);
                }
                if(event_flag==2 && !kick_print_flag) {
                    if (dv[3]>0&&!kick_print_flag) {
#pragma omp critical 
                        {
                            fout_sse_sn<<std::setw(WRITE_WIDTH)<<i+1
                                       <<std::setw(WRITE_WIDTH)<<dv[3]*bse_manager.vscale;
                            star[i].printColumn(fout_sse_sn, WRITE_WIDTH);
                            fout_sse_sn<<std::endl;
                        }

                        if (star.size()==1) {
                            std::cout<<"SN kick, vkick[IN]="<<dv[3]<<" ";
                            star[i].print(std::cout);
                            std::cout<<std::endl;
                        }
                        kick_print_flag=true;
                    }
                }
                double dt_miss = bse_manager.getDTMiss(output[i]);
                if (dt_miss!=0.0&&star[i].kw>=15) break;
            }
        }
        
        printSingleColumnTitle(std::cout);
        for (size_t i=0; i<star.size(); i++) {
            printSingleColumn(std::cout, mass0[i], star[i], output[i]);
        }

        if (output_flag) {
            fout_sse_type.close();
            fout_sse_sn.close();
        }
    }

    if (hyb.size()>0) {
        int nhyb = hyb.size();
        printBinaryTitle(std::cout);
#pragma omp parallel for schedule(dynamic)
        for (int i=0; i<nhyb; i++) {
            bse_manager.evolveStar(hyb[i].star[0],hyb[i].out[0],0);
            bse_manager.evolveStar(hyb[i].star[1],hyb[i].out[1],0);

            bse_manager.merge(hyb[i].star[0],hyb[i].star[1],hyb[i].out[0],hyb[i].out[1],hyb[i].semi,hyb[i].ecc);

            printBinary(std::cout,hyb[i]);
        }
    }

    return 0;
}
