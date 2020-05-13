#include <iostream>
#include <getopt.h>
#include <cmath>
#include <vector>
#include "bse_interface.h"
#include "../src/io.hpp"

int main(int argc, char** argv){

    int arg_label;
    int width=20;
    int n=5000;
    double m_min=0.08, m_max=150.0;
    double time=100.0;
    double dtmin=1.0;
    std::vector<double> mass0;

    auto printHelp= [&]() {
        std::cout<<"BSE_test [options] [initial mass of stars, can be multiple values, if not set, evolve an IMF with equal mass interal in Log scale]\n"
                 <<"    -s [D]: mimimum mass ("<<m_min<<")\n"
                 <<"    -e [D]: maximum mass ("<<m_max<<")\n"
                 <<"    -n [I]: number of stars when evolve an IMF ("<<n<<")\n"
                 <<"    -t [D]: evolve time ("<<time<<")[Myr]\n"
                 <<"    -d [D]: minimum time step ("<<dtmin<<")[Myr]\n"
                 <<"    -w [I]: print column width ("<<width<<")\n"
                 <<"    -h    : help\n";
    };

    IOParamsBSE bse_io;
    opterr = 0;
    bse_io.print_flag = true;
    int opt_used = bse_io.read(argc,argv);

    // reset optind
    optind=1;
    static struct option long_options[] = {{0,0,0,0}};

    int option_index;
    while ((arg_label = getopt_long(argc, argv, "s:e:n:t:d:w:h", long_options, &option_index)) != -1)
        switch (arg_label) {
        case 's':
            m_min = atof(optarg);
            std::cout<<"min mass: "<<m_min<<std::endl;
            break;
        case 'e':
            m_max = atof(optarg);
            std::cout<<"max mass: "<<m_max<<std::endl;
            break;
        case 'n':
            n = atof(optarg);
            std::cout<<"N: "<<n<<std::endl;
            break;
        case 't':
            time = atof(optarg);
            std::cout<<"finish time[Myr]: "<<time<<std::endl;
            break;
        case 'w':
            width = atoi(optarg);
            std::cout<<"print width: "<<width<<std::endl;
            break;
        case 'd':
            dtmin = atof(optarg);
            std::cout<<"minimum time step "<<dtmin<<std::endl;
            break;
        case 'h':
            printHelp();
            return 0;
        case '?':
            opt_used--;
            break;
        default:
            break;
        }        

    // argc and optind are 1 when no input is given
    opt_used += optind;
    // read initial mass list
    bool read_mass_flag = false;
    if (opt_used<argc) {
        while (opt_used<argc) 
            mass0.push_back(atof(argv[opt_used++]));
        read_mass_flag = true;
    }

    BSEManager bse_manager;

    bse_manager.initial(bse_io);

    assert(bse_manager.checkParams());

    if (read_mass_flag) n = mass0.size();
    StarParameter star[n];
    StarParameterOut output[n];

    // if no mass is read, use IMF
    if (!read_mass_flag) {
        double dm_factor = exp((log(m_max) - log(m_min))/n);

        mass0.push_back(m_min);
        for (int i=1; i<n; i++) {
            mass0.push_back(mass0.back()*dm_factor);
        }
    }

    // initial parameter
    for (int i=0; i<n; i++) {
        star[i].initial(mass0[i]);
    }

#pragma omp parallel for schedule(dynamic)
    for (int i=0; i<n; i++) {
        bse_manager.evolveStar(star[i],output[i],time);
        //while (star[i].tphys/bse_manager.tscale<time) {
        //    double dt = std::max(bse_manager.getTimeStep(star[i]),dtmin);
        //    dt = std::min(time-star[i].tphys/bse_manager.tscale, dt);
        //    double dt_miss=bse_manager.evolveStar(star[i],output[i],dt);
        //    if (dt_miss!=0.0&&star[i].kw>=15) break;
        //}
    }

    std::cout<<std::setw(width)<<"Mass_init[Msun]";
    StarParameter::printColumnTitle(std::cout, width);
    StarParameterOut::printColumnTitle(std::cout, width);
    std::cout<<std::endl;

    for (int i=0; i<n; i++) {
        std::cout<<std::setw(width)<<mass0[i];
        star[i].printColumn(std::cout, width);
        output[i].printColumn(std::cout, width);
        std::cout<<std::endl;
    }

    return 0;
}
