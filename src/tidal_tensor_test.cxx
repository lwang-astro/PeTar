#include <iostream>
#include <cstdio>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <cassert>
#include <unistd.h>
#include <particle_simulator.hpp>
#include "tidal_tensor.hpp"
#include "soft_ptcl.hpp"
#include "soft_force.hpp"
#include "static_variables.hpp"
#include <getopt.h>
#ifdef GALPY
#include "galpy_interface.h"
#endif

#ifdef GALPY
void addExtAcc(FPSoft& ptcl, double* pos_offset, GalpyManager& galpy_manager) {
        double pos[3] = {ptcl.pos[0] + pos_offset[0],
                         ptcl.pos[1] + pos_offset[1],
                         ptcl.pos[2] + pos_offset[2]};
        double acc[3]={0.0,0.0,0.0};
        PS::F64 pot=0.0;
        galpy_manager.calcAccPot(acc, pot, 0, pos, &ptcl.pos[0]);
        ptcl.acc[0] += acc[0];
        ptcl.acc[1] += acc[1];
        ptcl.acc[2] += acc[2];
        ptcl.pot_tot += pot;
        ptcl.pot_soft += pot;
        ptcl.pot_ext += pot;
}
#endif

void printAcc(FPSoft* ptcl, FPSoft& pcm, int n, TidalTensor& tt) {
    std::cout<<std::setw(3)<<"I"
             <<std::setw(PRINT_WIDTH*3+1)<<"pos"
             <<std::setw(PRINT_WIDTH*3+1)<<"acc(tt)"
             <<std::setw(PRINT_WIDTH*3+1)<<"acc(measured)"
             <<std::setw(PRINT_WIDTH*3+1)<<"diff"
             <<std::endl;
    for (int i=0; i<n; i++) {
        PS::F64vec acc;
        acc.x = acc.y = acc.z = 0.0;
        PS::F64vec pos = ptcl[i].pos - pcm.pos;
        tt.eval(&acc.x, pos);
        std::cout<<std::setw(3)<<i
                 <<std::setw(PRINT_WIDTH)<<pos.x
                 <<std::setw(PRINT_WIDTH)<<pos.y
                 <<std::setw(PRINT_WIDTH)<<pos.z
                 <<" "
                 <<std::setw(PRINT_WIDTH)<<acc.x
                 <<std::setw(PRINT_WIDTH)<<acc.y
                 <<std::setw(PRINT_WIDTH)<<acc.z
                 <<" "
                 <<std::setw(PRINT_WIDTH)<<ptcl[i].acc.x
                 <<std::setw(PRINT_WIDTH)<<ptcl[i].acc.y
                 <<std::setw(PRINT_WIDTH)<<ptcl[i].acc.z
                 <<" "
                 <<std::setw(PRINT_WIDTH)<<acc.x-ptcl[i].acc.x
                 <<std::setw(PRINT_WIDTH)<<acc.y-ptcl[i].acc.y
                 <<std::setw(PRINT_WIDTH)<<acc.z-ptcl[i].acc.z
                 <<std::endl;
    }
}

int main(int argc, char **argv){

#ifdef GALPY
    IOParamsGalpy galpy_io;
    bool unit_astro_flag = false;
#endif
    int arg_label;
    opterr = 0;
    int opt_used = 0;
    optind = 0;
    bool help_flag = false;
    std::cout<<std::setprecision(11);

    double gravitational_constant = 1.0;

    static struct option long_options[] = {
        {0,0,0,0}
    };
    int option_index;
    
    while ((arg_label = getopt_long(argc, argv, "-uh", long_options, &option_index)) != -1)
        switch (arg_label) {
        case 'u':
#ifdef GALPY
            unit_astro_flag = true;
#endif
            std::cout<<"Use the astronomical unit set (Myr, pc, Msun)\n";
            opt_used ++;
            gravitational_constant = 0.00449830997959438; // pc^3/(Msun*Myr^2)
            break;
        case 'h':
            std::cout<<"petar.tt.test [options] [data filename]\n"
                     <<"The data files content:\n"
                     <<"1st line: N_t, N_p, xt, yt, zt, rscale\n"
                     <<"    N_t: number of particles for measuring tidal force\n"
                     <<"    N_p: number of particles to provide soft force\n"
                     <<"    xt, yt, zt: center position of the tidal tensor box\n"
                     <<"    rscale: distance scale for creating TT box\n"
#ifdef GALPY
                     <<"2nd line: position and velocity offsets of the center referring to the galactic center\n"
#endif
                     <<"Following lines: particle data (first N_t particles for measurement; then N_p particles for forces)\n"
                     <<"Each particle line contains: \n";
            ParticleBase::printTitleWithMeaning(std::cout,0,13);
            std::cout<<"Options:\n"
                     <<"    -h    : help\n";
            help_flag=true;
            break;
        case '?':
            opt_used +=2;
            break;
        default:
            break;
        }


#ifdef GALPY
    galpy_io.print_flag = true;
    galpy_io.read(argc,argv);
#endif
    
    if (help_flag) return 0;

#ifdef GALPY
    GalpyManager galpy_manager;
    if (unit_astro_flag) galpy_io.setStdUnit();
    galpy_manager.initial(galpy_io,true);
#endif

    opt_used ++;
    std::string filename;
    if (opt_used<argc) {
        filename=argv[argc-1];
        std::cout<<"Reading data file name: "<<filename<<std::endl;
    }
    else {
        std::cerr<<"Error: no input data file\n";
        abort();
    }

    // open data file
    FILE* fin;
    if ( (fin = fopen(filename.c_str(),"r")) == NULL) {
        fprintf(stderr,"Error: Cannot open input file %s.\n",filename.c_str());
        abort();
    }

    int n_check, n_pert; // number of particles
    double r_scale;
    FPSoft ptcl_binary_cm;

    int rcount=fscanf(fin, "%d %d %lf %lf %lf %lf", &n_check, &n_pert, &ptcl_binary_cm.pos.x, &ptcl_binary_cm.pos.y, &ptcl_binary_cm.pos.z, &r_scale);
    if (rcount<6) {
        std::cerr<<"Error: reading first line data fail, required 6, given"<<rcount<<std::endl;
        abort();
    }

#ifdef GALPY
    double pos_off[3], vel_off[3];
    fscanf(fin, "%lf %lf %lf %lf %lf %lf", &pos_off[0], &pos_off[1], &pos_off[2], &vel_off[0], &vel_off[1], &vel_off[2]);
#endif

    int n_tot = n_check + n_pert;
    ParticleBase ptcl[n_tot];
    for (int i=0; i<n_tot; i++) {
        ptcl[i].readAscii(fin);
    }
    //ptcl_binary_cm.mass = ptcl[0].mass + ptcl[1].mass;
    //ptcl_binary_cm.pos = (ptcl[0].mass*ptcl[0].pos + ptcl[1].mass*ptcl[1].pos)/ptcl_binary_cm.mass;
    //ptcl_binary_cm.vel = (ptcl[0].mass*ptcl[0].vel + ptcl[1].mass*ptcl[1].vel)/ptcl_binary_cm.mass;

    // create box
    int n_tt = TidalTensor::getParticleN();
    FPSoft ptcl_tt[n_tt];
    FPSoft ptcl_check[n_check];
    TidalTensor::createTidalTensorMeasureParticles(ptcl_tt, ptcl_binary_cm, r_scale);
    
    int Nepi=n_tt + n_check + 1;
    int Nepj=n_pert;

    // print particle
    std::cout<<"Tidal tensor box and center:\n";
    for (int i=0; i<n_tt+1; i++) std::cout<<"I"<<i<<" "<<ptcl_tt[i].pos<<std::endl;
    std::cout<<"Measuring points:\n";
    for (int i=0; i<n_check; i++) {
        ptcl_check[i].pos = ptcl[i].pos;
        ptcl_check[i].mass = ptcl[i].mass;

        std::cout<<"I"<<i<<" "<<ptcl_check[i].pos<<std::endl;
    }

    if (Nepj>0) {
        EPISoft epi[Nepi];
        for (int i=0; i<n_tt; i++) {
            epi[i].id = i+1;
            epi[i].pos = ptcl_tt[i].pos;
            epi[i].r_search = r_scale*2;
            ptcl_tt[i].group_data.artificial.setParticleTypeToSingle();
        }
        epi[n_tt].id = 9;
        epi[n_tt].pos = ptcl_binary_cm.pos;
        epi[n_tt].r_search  = r_scale*2;

        for (int i=0; i<n_check; i++) {
            int k = i+n_tt+1;
            epi[k].id = k+1;
            epi[k].pos = ptcl[i].pos;
            epi[k].r_search = r_scale*2;
        }

        EPJSoft epj[Nepj];
        for (int i=0; i<Nepj; i++) {
            epj[i].id = n_check + n_tt + 1 + i;
            epj[i].mass=ptcl[i+n_check].mass;
            epj[i].pos= ptcl[i+n_check].pos;
        }

        // calculate force
        ForceSoft force_sp[Nepi]; //8: box; last cm
        for (int i=0; i<Nepi; i++) force_sp[i].acc = PS::F64vec(0,0,0);
        CalcForceEpSpMonoNoSimd f_ep_sp;
        ForceSoft::grav_const = gravitational_constant;
        EPISoft::eps = 0.0;
        
        std::cout<<"EPJ mass,pos: (perturber) \n";
        for (int i=0; i<Nepj; i++) std::cout<<"I"<<i<<" "<<epj[i].pos<<std::endl;
    
        f_ep_sp(epi, Nepi, epj, Nepj, force_sp);
        for (int i=0; i<n_tt; i++) ptcl_tt[i].copyFromForce(force_sp[i]);
        ptcl_binary_cm.copyFromForce(force_sp[n_tt]);
        for (int i=0; i<n_check; i++) ptcl_check[i].copyFromForce(force_sp[i+n_tt+1]);
    }

#ifdef GALPY
    addExtAcc(ptcl_binary_cm, pos_off, galpy_manager);
    for (int i=0; i<n_tt; i++) addExtAcc(ptcl_tt[i], pos_off, galpy_manager);
    for (int i=0; i<n_check; i++) {
        addExtAcc(ptcl_check[i], pos_off, galpy_manager);
        ptcl_check[i].acc -= ptcl_binary_cm.acc;
    }
#endif

    // print acc at box
    std::cout<<"Original acc at tensor box:\n";
    for (int i=0; i<n_tt; i++) 
        std::cout<<"I"<<i<<" "<<ptcl_tt[i].acc<<std::endl;

    TidalTensor::subtractCMForce(ptcl_tt, ptcl_binary_cm);
    TidalTensor tt;
    tt.fit(ptcl_tt, ptcl_binary_cm, r_scale);
    
    // print matrix
    std::cout<<"Tidal tensor:\n";
    tt.print(std::cout,PRINT_WIDTH);

    // evolve force at box and compare with original force;
    std::cout<<"Check acc at box sample points\n";
    printAcc(ptcl_tt, ptcl_binary_cm, n_tt, tt);

    std::cout<<"Check acc at center\n";
    ptcl_binary_cm.acc = PS::F64vec(0.0);
    printAcc(&ptcl_binary_cm, ptcl_binary_cm, 1, tt);
    
    // evolve force at the measuring points
    std::cout<<"Check acc at measuring points\n";
    printAcc(ptcl_check, ptcl_binary_cm, n_check, tt);

    return 0;
}

