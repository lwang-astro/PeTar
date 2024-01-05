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

    Float a1[3], a2[3], ad1[3], ad2[3], s1[3], s2[3];
    pn.calcAccJerkPN(a1, a2, ad1, ad2, NULL, NULL, p1, p2);

}
