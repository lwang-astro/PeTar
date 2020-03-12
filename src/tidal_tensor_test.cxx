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

#define PRINT_WIDTH 14

int main(int argc, char **argv){

    int n_opt=0;
    int arg_label;
    double r_scale;

    while ((arg_label = getopt(argc, argv, "r:h")) != -1)
        switch (arg_label) {
        case 'r':
            r_scale = atof(optarg);
            break;
        case 'h':
            std::cout<<"tidal_tensor_test.out [options] [particle data]\n"
                     <<"particle data file contain three particles: first two is binary, third is distant particle\n"
                     <<"options:\n"
                     <<"      -r [float]: distance scale for creating TT box (1.0)\n"
                     <<"      -h:         help\n";
            return 0;
        default:
            break;
        }

    std::string filename;
    if (argc-n_opt>1) filename=argv[argc-1];
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

    ParticleBase ptcl[3];
    FPSoft ptcl_binary_cm;
    for (int i=0; i<3; i++) {
        ptcl[i].readAscii(fin);
    }
    ptcl_binary_cm.mass = ptcl[0].mass + ptcl[1].mass;
    ptcl_binary_cm.pos = (ptcl[0].mass*ptcl[0].pos + ptcl[0].mass*ptcl[0].pos)/ptcl_binary_cm.mass;
    ptcl_binary_cm.vel = (ptcl[0].mass*ptcl[0].vel + ptcl[0].mass*ptcl[0].vel)/ptcl_binary_cm.mass;
    
    // create box
    FPSoft ptcl_tt[8];
    TidalTensor::createTidalTensorMeasureParticles(ptcl_tt, ptcl_binary_cm, r_scale);
    
    int Nepi=9, Nepj=1;
    EPISoft epi[Nepi];
    for (int i=0; i<8; i++) {
        epi[i].id = i+1;
        epi[i].pos = ptcl_tt[i].pos;
        epi[i].r_search = r_scale*2;
        ptcl_tt[i].group_data.artificial.setParticleTypeToSingle();
    }
    epi[8].id = 9;
    epi[8].pos = ptcl_binary_cm.pos;
    epi[8].r_search  = r_scale*2;

    EPJSoft epj[Nepj];
    epj[0].id = 10;
    epj[0].mass=ptcl[2].mass;
    epj[0].pos= ptcl[2].pos;

    // calculate force
    ForceSoft force_sp[Nepi]; //8: box; last cm
    for (int i=0; i<Nepi; i++) force_sp[i].acc = PS::F64vec(0,0,0);
    CalcForceEpSpMonoNoSimd f_ep_sp;
    // print particle
    std::cout<<"EPI pos: (8+cm) \n";
    for (int i=0; i<Nepi; i++) std::cout<<"I"<<i<<" "<<epi[i].pos<<std::endl;
    std::cout<<"EPJ mass,pos: (perturber) \n";
    std::cout<<epj[0].mass<<" "<<epj[0].pos<<std::endl;
    
    f_ep_sp(epi, Nepi, epj, Nepj, force_sp);
    for (int i=0; i<8; i++) ptcl_tt[i].copyFromForce(force_sp[i]);
    ptcl_binary_cm.copyFromForce(force_sp[8]);

    TidalTensor::subtractCMForce(ptcl_tt, ptcl_binary_cm);
    TidalTensor tt;
    tt.fit(ptcl_tt, ptcl_binary_cm, r_scale);
    
    // print matrix
    tt.print(std::cout,PRINT_WIDTH);
    
    // evolve force at box and compare with original force;
    std::cout<<"Check acc at box sample points\n";
    std::cout<<std::setw(3)<<"I"
             <<std::setw(PRINT_WIDTH*3+1)<<"pos"
             <<std::setw(PRINT_WIDTH*3+1)<<"acc(tt)"
             <<std::setw(PRINT_WIDTH*3+1)<<"acc(measured)"
             <<std::setw(PRINT_WIDTH*3+1)<<"diff"
             <<std::endl;
    for (int i=0; i<8; i++) {
        PS::F64vec acc;
        acc.x = acc.y = acc.z = 0.0;
        PS::F64vec pos = ptcl_tt[i].pos - ptcl_binary_cm.pos;
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
                 <<std::setw(PRINT_WIDTH)<<ptcl_tt[i].acc.x
                 <<std::setw(PRINT_WIDTH)<<ptcl_tt[i].acc.y
                 <<std::setw(PRINT_WIDTH)<<ptcl_tt[i].acc.z
                 <<" "
                 <<std::setw(PRINT_WIDTH)<<acc.x-ptcl_tt[i].acc.x
                 <<std::setw(PRINT_WIDTH)<<acc.y-ptcl_tt[i].acc.y
                 <<std::setw(PRINT_WIDTH)<<acc.z-ptcl_tt[i].acc.z
                 <<std::endl;
    }

    return 0;
}

