#include <iostream>
#include <cstdio>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <cassert>
#include <random>

#include "pn.hpp"
#include "pn_BH.h"
#include "astro_units.hpp"
#include <particle_simulator.hpp>
#include "particle_base.hpp"
#include "Common/Float.h"
#include "Common/binary_tree.h"

int main(int argc, char **argv){
    
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.1,1);    
        
    ParticleBase p1,p2;

    COMM::Binary bin;
    
    bin.semi = 1e-9;
    bin.ecc = 0.1;
    bin.incline = 0.0;
    bin.rot_horizon = 0.0;
    bin.rot_self = 0.0;
    bin.t_peri = 0.0;
    bin.period = 0.0;
    bin.ecca = 3.14;
    bin.m1 = 10.0;
    bin.m2 = 20.0;
    bin.r = 0.0;
    bin.stab = 0;

    const Float G = G_ASTRO;
    bin.calcParticles(p1, p2, G);

    p1.printColumnTitle(std::cout);
    p2.printColumnTitle(std::cout);
    std::cout<<std::endl;
    p1.printColumn(std::cout);
    p2.printColumn(std::cout);

    std::cout<<std::endl;

    PostNewtonian pn;
    pn.speed_of_light = 2.99792458e5*KMS_TO_PCMYR;
    pn.gravitational_constant = G;

    bool used_pn_order[6] = {true};
    used_pn_order[5] = false;

    Float a1[3], a2[3], ad1[3], ad2[3], s1[3], s2[3];
    pn.calcAccJerkPN(a1, a2, ad1, ad2, NULL, NULL, p1, p2);

    Float a1r[3], a2r[3], ad1r[3], ad2r[3], s1r[3], s2r[3];
    calc_force_pn_BH(p1.mass, p1.pos, p1.vel, s1r, p2.mass, p2.pos, p2.vel, s2r, pn.speed_of_light, used_pn_order, 0, a1r, ad1r, a2r, ad2r);

    int width = 10;
    std::cout<<std::setw(width)<<"a1.x"
             <<std::setw(width)<<"a1.y"
             <<std::setw(width)<<"a1.z";
    std::cout<<std::setw(width)<<"a2.x"
             <<std::setw(width)<<"a2.y"
             <<std::setw(width)<<"a2.z";
    std::cout<<std::endl;

    std::cout<<std::setw(width)<<a1[0]
             <<std::setw(width)<<a1[1]
             <<std::setw(width)<<a1[2];
    std::cout<<std::setw(width)<<a2[0]
             <<std::setw(width)<<a2[1]
             <<std::setw(width)<<a2[2];
    std::cout<<std::endl;

    std::cout<<std::setw(width)<<a1r[0]
             <<std::setw(width)<<a1r[1]
             <<std::setw(width)<<a1r[2];
    std::cout<<std::setw(width)<<a2r[0]
             <<std::setw(width)<<a2r[1]
             <<std::setw(width)<<a2r[2];
    std::cout<<std::endl;

    std::cout<<std::setw(width)<<"ad1.x"
             <<std::setw(width)<<"ad1.y"
             <<std::setw(width)<<"ad1.z";
    std::cout<<std::setw(width)<<"ad2.x"
             <<std::setw(width)<<"ad2.y"
             <<std::setw(width)<<"ad2.z";
    std::cout<<std::endl;

    std::cout<<std::setw(width)<<ad1[0]
             <<std::setw(width)<<ad1[1]
             <<std::setw(width)<<ad1[2];
    std::cout<<std::setw(width)<<ad2[0]
             <<std::setw(width)<<ad2[1]
             <<std::setw(width)<<ad2[2];
    std::cout<<std::endl;
    
    std::cout<<std::setw(width)<<ad1r[0]
             <<std::setw(width)<<ad1r[1]
             <<std::setw(width)<<ad1r[2];
    std::cout<<std::setw(width)<<ad2r[0]
             <<std::setw(width)<<ad2r[1]
             <<std::setw(width)<<ad2r[2];
    std::cout<<std::endl;

    std::cout<<std::setw(width)<<"s1.x"
             <<std::setw(width)<<"s1.y"
             <<std::setw(width)<<"s1.z";
    std::cout<<std::setw(width)<<"s2.x"
             <<std::setw(width)<<"s2.y"
             <<std::setw(width)<<"s2.z";
    std::cout<<std::endl;

    std::cout<<std::setw(width)<<s1[0]
             <<std::setw(width)<<s1[1]
             <<std::setw(width)<<s1[2];
    std::cout<<std::setw(width)<<s2[0]
             <<std::setw(width)<<s2[1]
             <<std::setw(width)<<s2[2];
    std::cout<<std::endl;

    std::cout<<std::setw(width)<<s1r[0]
             <<std::setw(width)<<s1r[1]
             <<std::setw(width)<<s1r[2];
    std::cout<<std::setw(width)<<s2r[0]
             <<std::setw(width)<<s2r[1]
             <<std::setw(width)<<s2r[2];
    std::cout<<std::endl;

    
    return 0;
}
