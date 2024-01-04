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

    const Float G = G_ASTRO;

}
