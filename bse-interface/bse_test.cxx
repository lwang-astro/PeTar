#include <iostream>
#include <fstream>
#include <getopt.h>
#include <cmath>
#include <vector>
#include <string>
#include "bse_interface.h"
#include "../src/io.hpp"

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
    std::string fprint_name;
    std::string fbin_name;
    std::string fsin_name;
    bool read_mass_flag = false;

    auto printHelp= [&]() {
        std::cout<<"The tool to evolve single stars or binaries using SSE/BSE\n"
                 <<"Usage: petar.bse [options] [initial mass of stars, can be multiple values]\n"
                 <<"       If no initial mass or no single/binary table (-s or -b) is provided, N single stars (-n) with equal mass interal in Log scale will be evolved\n"
                 <<"       When single table is provided, times and types can be set individually\n"
                 <<"       The default unit set is: Msun, Myr. If input data have different units [IN], please modify the scaling fators\n"
                 <<"Options:\n"
                 <<"    -n [I]: number of stars when evolve an IMF ("<<n<<")\n"
                 <<"    -i [D]: start time ("<<time0<<") [IN]\n"
                 <<"    -t [D]: finish time ("<<time<<") [IN]\n"
                 <<"        --mmin [D]: mimimum mass ("<<m_min<<") [M*]\n"
                 <<"        --mmax [D]: maximum mass ("<<m_max<<") [M*]\n"
                 <<"    -d [D]: minimum time step ("<<dtmin<<")[IN]\n"
                 <<"    -s [S]: a file of single table: First line: number of single (unit:IN); After: mass, type, time per line\n"
                 <<"    -b [S]: a file of binary table: First line: number of binary (unit:IN); After: m1, m2, type1, type2, period, ecc, time per line\n"
                 <<"    -w [I]: print column width ("<<width<<")\n"
                 <<"    -o [S]: a file to output data every step\n"
                 <<"    -h    : help\n";
    };

    IOParamsBSE bse_io;
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
    while ((arg_label = getopt_long(argc, argv, "-n:i:t:w:d:s:b:o:h", long_options, &option_index)) != -1)
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
        case 'o':
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

    if (help_flag) return 0;

    BSEManager bse_manager;
    bse_manager.initial(bse_io,true);
    assert(bse_manager.checkParams());

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

    auto printSingleTitle=[&](std::ostream & _fout) {
        _fout<<std::setw(width)<<"Mass_init[Msun]";
        StarParameter::printColumnTitle(_fout, width);
        StarParameterOut::printColumnTitle(_fout, width);
        _fout<<std::endl;
    };

    auto printSingle=[&](std::ostream & _fout, double& _mass0, StarParameter& _star, StarParameterOut& _out) {
        _fout<<std::setw(width)<<_mass0*bse_manager.mscale;
        _star.printColumn(_fout, width);
        _out.printColumn(_fout, width);
        _fout<<std::endl;
    };

    // first check whether binary exist
    if (bin.size()>0) {
        int nbin = bin.size();
#pragma omp parallel for schedule(dynamic)
        for (int i=0; i<nbin; i++) {
            // output file
            std::ofstream fout;
            bool print_flag=false;
            if (fprint_name!="") {
                std::string fi = fprint_name+std::string(".b.")+std::to_string(i);
                fout.open(fi.c_str(),std::ofstream::out);
                printBinaryTitle(fout);
                print_flag = true;
            }

            bool kick_print_flag[2]={false,false};
            // evolve
            double tend = time*bse_manager.tscale;
            bin[i].period0 = bin[i].period;
            bin[i].ecc0 = bin[i].ecc;
            int bin_type_init=0;
            int bin_type_last=1;
            while (bse_manager.getTime(bin[i].star[0])<tend) {
                // time step
                //double dt1 = bse_manager.getTimeStep(bin[i].star[0]);
                //double dt2 = bse_manager.getTimeStep(bin[i].star[1]);
                //double dt = std::min(dt1,dt2);
                double dt = bse_manager.getTimeStepBinary(bin[i].star[0],bin[i].star[1],bin[i].semi,bin[i].ecc,bin_type_last);
                dt = std::max(dt,dtmin);
                dt = std::min(tend-bse_manager.getTime(bin[i].star[0]), dt);
                double mtot = bin[i].star[0].mt + bin[i].star[1].mt;
                double period_myr = bin[i].period*bse_manager.tscale;
                double G= 0.00449830997959438;
                const double PI = 4.0*atan(1.0);
                double pc_to_rsun = 44334448.006896;
                bin[i].semi = std::pow(period_myr*period_myr*G*mtot/(4*PI*PI),1.0/3.0)*pc_to_rsun;
                
                // evolve function
                int error_flag=bse_manager.evolveBinary(bin[i].star[0],bin[i].star[1],bin[i].out[0],bin[i].out[1],bin[i].semi,bin[i].period,bin[i].ecc,bin[i].bse_event, bin_type_init, dt);
                int nmax = bin[i].bse_event.getEventNMax();
                for (int k=0; k<nmax; k++) {
                    int binary_type = bin[i].bse_event.getType(k);
                    if (binary_type>0) {
                        bin_type_last = binary_type;
#pragma omp critical
                        {
                            std::cout<<" ID="<<i<<" index="<<k<<" "<<" dt="<<dt;
                            bse_manager.printBinaryEventOne(std::cout, bin[i].bse_event, k);
                            std::cout<<std::endl;
                        }
                    }
                    else if (binary_type<0) break;
                }
                bin_type_init = bin_type_last;

                for (int k=0; k<2; k++) {
                    double dv[4];
                    dv[3] = bse_manager.getVelocityChange(dv,bin[i].out[k]);
                    if (dv[3]>0&&!kick_print_flag[k]) {
#pragma omp critical 
                        {
                            std::cout<<"SN kick, i="<<i<<" vkick[IN]="<<dv[3]<<" ";
                            bin[i].star[k].print(std::cout);
                            std::cout<<std::endl;
                        }
                        kick_print_flag[k]=true;
                    }
                }
                if (error_flag<0) {
#pragma omp critical 
                    {
                        std::cerr<<"Error: i="<<i<<" mass0[IN]="<<bin[i].m1<<" "<<bin[i].m2<<" period[IN]="<<bin[i].period<<" ecc[IN]="<<bin[i].ecc
                                 <<std::endl;
                        std::cerr<<"Star 1:";
                        bin[i].star[0].print(std::cerr);
                        std::cerr<<"\nStar 2:";
                        bin[i].star[1].print(std::cerr);
                        std::cerr<<std::endl;
                        std::cerr<<std::endl;
                    }
                }
                double dt_miss = bse_manager.getDTMiss(bin[i].out[0]);

                if (print_flag) printBinary(fout,bin[i]);
                if (dt_miss!=0.0&&bin[i].star[0].kw>=15&&bin[i].star[1].kw>=15) break;
            }
            fout.close();
        }

        printBinaryTitle(std::cout);
        for (int i=0; i<nbin; i++) printBinary(std::cout, bin[i]);
    }

    if (star.size()>0) {
        StarParameterOut output[star.size()];

#pragma omp parallel for schedule(dynamic)
        for (size_t i=0; i<star.size(); i++) {
            // output file
            std::ofstream fout;
            bool print_flag=false;
            if (fprint_name!="") {
                std::string fi = fprint_name+std::string(".s.")+std::to_string(i);
                fout.open(fi.c_str(),std::ofstream::out);
                printSingleTitle(fout);
                print_flag = true;
            }

            //int error_flag = bse_manager.evolveStar(star[i],output[i],time);
            bool kick_print_flag=false;
            double tend = time*bse_manager.tscale;
            while (bse_manager.getTime(star[i])<tend) {
                double dt = std::max(bse_manager.getTimeStepStar(star[i]),dtmin);
                dt = std::min(tend-bse_manager.getTime(star[i]), dt);
                int error_flag=bse_manager.evolveStar(star[i],output[i],dt);
                double dv[4];
                dv[3] = bse_manager.getVelocityChange(dv, output[i]);
                if (dv[3]>0&&!kick_print_flag) {
#pragma omp critical 
                    {
                        std::cout<<"SN kick, i="<<i<<" vkick[IN]="<<dv[3]<<" ";
                        star[i].print(std::cout);
                        std::cout<<std::endl;
                    }
                    kick_print_flag=true;
                }
                if (error_flag<0) {
#pragma omp critical 
                    {
                        std::cerr<<"Error: i="<<i<<" mass0[IN]="<<mass0[i]<<" ";
                        star[i].print(std::cerr);
                        std::cerr<<std::endl;
                    }
                }
                double dt_miss = bse_manager.getDTMiss(output[i]);
                if (print_flag) printSingle(fout, mass0[i], star[i], output[i]);
                if (dt_miss!=0.0&&star[i].kw>=15) break;
            }
            fout.close();
        }

        printSingleTitle(std::cout);
        for (size_t i=0; i<star.size(); i++) {
            printSingle(std::cout, mass0[i], star[i], output[i]);
        }
    }

    return 0;
}
