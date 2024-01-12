#include <iostream>
#include <cstdio>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <cassert>
#include <random>
#include <cassert>

#define ASSERT(expr) assert(expr)

#include "pn.hpp"
#include "pn_BH.h"
#include "astro_units.hpp"
#include <particle_simulator.hpp>
#include "particle_base.hpp"
#include "Common/Float.h"
#include "Common/binary_tree.h"

class ParticleSpin: public ParticleBase{
public:
    Float spin[3];
};

int main(int argc, char **argv){
    
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.1,1);    
        
    ParticleSpin p1,p2;

    COMM::Binary bin;
    
    bin.semi = 1e-9;
    bin.ecc = 0.1;
    bin.incline = 2.0;
    bin.rot_horizon = 1.0;
    bin.rot_self = 1.0;
    bin.t_peri = 0.0;
    bin.period = 0.0;
    bin.ecca = 3.14;
    bin.m1 = 10.0;
    bin.m2 = 20.0;
    bin.r = 0.0;
    bin.stab = 0;

    const Float G = G_ASTRO;
    bin.calcParticles(p1, p2, G);
    p1.spin[0] = 0.5;
    p1.spin[1] = 0.3;
    p1.spin[2] = 0.2;
    p2.spin[0] = 0.3;
    p2.spin[1] = 0.5;
    p2.spin[2] = -0.2;

    p1.spin[0] = 0.0;
    p1.spin[1] = 0.0;
    p1.spin[2] = 0.0;
    p2.spin[0] = 0.0;
    p2.spin[1] = 0.0;
    p2.spin[2] = 0.0;

    p1.printColumnTitle(std::cout);
    p2.printColumnTitle(std::cout);
    std::cout<<std::endl;
    p1.printColumn(std::cout);
    p2.printColumn(std::cout);

    std::cout<<std::endl;

    PostNewtonian pn;
    pn.speed_of_light = 2.99792458e5*KMS_TO_PCMYR;
    pn.gravitational_constant = G;
    pn.precession_criterion = 1e-9;

    Float a1[6][3], a2[6][3], ad1[6][3], ad2[6][3], s1[3]={0.0}, s2[3]={0.0};
    bool used_pn_orders[6];

    bin.calcOrbit(p1,p2, G);
    pn.setUsedPNOrders(used_pn_orders, bin.r, bin.m1+bin.m2);
    used_pn_orders[5] = true;
    std::cout<<"Used orders: "
             <<"  PN1 "<<used_pn_orders[0]
             <<"  PN2 "<<used_pn_orders[1]
             <<"  PN2.5 "<<used_pn_orders[2]
             <<"  PN3 "<<used_pn_orders[3]
             <<"  PN3.5 "<<used_pn_orders[4]
             <<"  Spin "<<used_pn_orders[5]<<std::endl;
    pn.calcAccJerkPN(a1, a2, ad1, ad2, s1, s2, p1, p2, used_pn_orders);

    int usedOrNot[6] = {1,1,1,1,1,1};
    Float a1r[6][3], a2r[6][3], ad1r[6][3], ad2r[6][3], s1r[3], s2r[3];
    for (int i=0; i<3; i++) {
        s1r[i] = p1.spin[i];
        s2r[i] = p2.spin[i];
    }

    // unit scale, let G = 1
    Float vscale = std::sqrt(G);
    Float tscale = 1/std::sqrt(G);
    Float fscale = vscale/tscale;
    Float jscale = vscale/tscale/tscale;
    //Float Gfactor[6] = {1, G, G*G, G*G*std::sqrt(G), G*G*G, G*G*G*std::sqrt(G)};
    //Float jfactor[6] = {1, G*std::sqrt(G), G*G, G*G*std::sqrt(G), G*G*G, G*G*G*std::sqrt(G)};

    Float vs1[3] = {p1.vel[0]/vscale, p1.vel[1]/vscale, p1.vel[2]/vscale};
    Float vs2[3] = {p2.vel[0]/vscale, p2.vel[1]/vscale, p2.vel[2]/vscale};
    
    calc_force_pn_BH(p1.mass, &p1.pos[0], vs1, s1r,
                     p2.mass, &p2.pos[0], vs2, s2r,
                     pn.speed_of_light/vscale, usedOrNot, 1, a1r, ad1r, a2r, ad2r);

    int width = 15;
    std::string pn_name[6]={"N","pn1", "pn2", "pn2.5", "pn3", "pn3.5"};

    for (int i=0; i<5; i++) {
        std::cout<<std::setw(width)<<"";
        std::cout<<std::setw(width)<<"a1_"+pn_name[i]+".x"
                 <<std::setw(width)<<"a1_"+pn_name[i]+".y"
                 <<std::setw(width)<<"a1_"+pn_name[i]+".z";
        std::cout<<std::setw(width)<<"a2_"+pn_name[i]+".x"
                 <<std::setw(width)<<"a2_"+pn_name[i]+".y"
                 <<std::setw(width)<<"a2_"+pn_name[i]+".z";
        std::cout<<std::endl;

        std::cout<<std::setw(width)<<"Res";
        std::cout<<std::setw(width)<<a1[i][0]
                 <<std::setw(width)<<a1[i][1]
                 <<std::setw(width)<<a1[i][2];
        std::cout<<std::setw(width)<<a2[i][0]
                 <<std::setw(width)<<a2[i][1]
                 <<std::setw(width)<<a2[i][2];
        std::cout<<std::endl;

        for (int k=0; k<3; k++) {
            a1r[i][k] *= fscale;
            a2r[i][k] *= fscale;
        }

        std::cout<<std::setw(width)<<"Ref";
        std::cout<<std::setw(width)<<a1r[i][0]
                 <<std::setw(width)<<a1r[i][1]
                 <<std::setw(width)<<a1r[i][2];
        std::cout<<std::setw(width)<<a2r[i][0]
                 <<std::setw(width)<<a2r[i][1]
                 <<std::setw(width)<<a2r[i][2];
        std::cout<<std::endl;

        Float da[6];
        da[0] = (a1r[i][0]-a1[i][0])/(a1r[i][0]+1e-64); 
        da[1] = (a1r[i][1]-a1[i][1])/(a1r[i][1]+1e-64); 
        da[2] = (a1r[i][2]-a1[i][2])/(a1r[i][2]+1e-64); 
        da[3] = (a2r[i][0]-a2[i][0])/(a2r[i][0]+1e-64); 
        da[4] = (a2r[i][1]-a2[i][1])/(a2r[i][1]+1e-64); 
        da[5] = (a2r[i][2]-a2[i][2])/(a2r[i][2]+1e-64); 

        std::cout<<std::setw(width)<<"diff(should=0)";
        std::cout<<std::setw(width)<<da[0] 
                 <<std::setw(width)<<da[1] 
                 <<std::setw(width)<<da[2];
        std::cout<<std::setw(width)<<da[3] 
                 <<std::setw(width)<<da[4] 
                 <<std::setw(width)<<da[5];
        std::cout<<std::endl;

        Float diff_max = 0;
        for (int k=0; k<6; k++) {
            diff_max = std::max(std::fabs(da[k]), diff_max);
        }
        if (diff_max>ROUND_OFF_ERROR_LIMIT*10) {
            std::cerr<<"Test failed! difference > round off error\n";
            abort();
        }

    }

    for (int i=0; i<5; i++) {
        std::cout<<std::setw(width)<<"";
        std::cout<<std::setw(width)<<"ad1_"+pn_name[i]+".x"
                 <<std::setw(width)<<"ad1_"+pn_name[i]+".y"
                 <<std::setw(width)<<"ad1_"+pn_name[i]+".z";
        std::cout<<std::setw(width)<<"ad2_"+pn_name[i]+".x"
                 <<std::setw(width)<<"ad2_"+pn_name[i]+".y"
                 <<std::setw(width)<<"ad2_"+pn_name[i]+".z";
        std::cout<<std::endl;

        std::cout<<std::setw(width)<<"Res";
        std::cout<<std::setw(width)<<ad1[i][0]
                 <<std::setw(width)<<ad1[i][1]
                 <<std::setw(width)<<ad1[i][2];
        std::cout<<std::setw(width)<<ad2[i][0]
                 <<std::setw(width)<<ad2[i][1]
                 <<std::setw(width)<<ad2[i][2];
        std::cout<<std::endl;

        for (int k=0; k<3; k++) {
            ad1r[i][k] *= jscale;
            ad2r[i][k] *= jscale;
        }

        std::cout<<std::setw(width)<<"Ref";
        std::cout<<std::setw(width)<<ad1r[i][0]
                 <<std::setw(width)<<ad1r[i][1]
                 <<std::setw(width)<<ad1r[i][2];
        std::cout<<std::setw(width)<<ad2r[i][0]
                 <<std::setw(width)<<ad2r[i][1]
                 <<std::setw(width)<<ad2r[i][2];
        std::cout<<std::endl;

        Float da[6];
        da[0] = (ad1r[i][0]-ad1[i][0])/(ad1r[i][0]+1e-64); 
        da[1] = (ad1r[i][1]-ad1[i][1])/(ad1r[i][1]+1e-64); 
        da[2] = (ad1r[i][2]-ad1[i][2])/(ad1r[i][2]+1e-64); 
        da[3] = (ad2r[i][0]-ad2[i][0])/(ad2r[i][0]+1e-64); 
        da[4] = (ad2r[i][1]-ad2[i][1])/(ad2r[i][1]+1e-64); 
        da[5] = (ad2r[i][2]-ad2[i][2])/(ad2r[i][2]+1e-64); 

        std::cout<<std::setw(width)<<"diff(should=0)";
        std::cout<<std::setw(width)<<da[0] 
                 <<std::setw(width)<<da[1] 
                 <<std::setw(width)<<da[2];
        std::cout<<std::setw(width)<<da[3] 
                 <<std::setw(width)<<da[4] 
                 <<std::setw(width)<<da[5];
        std::cout<<std::endl;

        Float diff_max = 0;
        for (int k=0; k<6; k++) {
            diff_max = std::max(std::fabs(da[k]), diff_max);
        }
        if (diff_max>ROUND_OFF_ERROR_LIMIT*1e2) {
            std::cerr<<"Test failed! difference > round off error\n";
            abort();
        }
    
    }

    std::cout<<std::setw(width)<<"";
    std::cout<<std::setw(width)<<"s1.x"
             <<std::setw(width)<<"s1.y"
             <<std::setw(width)<<"s1.z";
    std::cout<<std::setw(width)<<"s2.x"
             <<std::setw(width)<<"s2.y"
             <<std::setw(width)<<"s2.z";
    std::cout<<std::endl;

    std::cout<<std::setw(width)<<"res";
    std::cout<<std::setw(width)<<s1[0]
             <<std::setw(width)<<s1[1]
             <<std::setw(width)<<s1[2];
    std::cout<<std::setw(width)<<s2[0]
             <<std::setw(width)<<s2[1]
             <<std::setw(width)<<s2[2];
    std::cout<<std::endl;

    Float dspin[6];
    dspin[0] = (s1r[0]-p1.spin[0]-s1[0])/(s1[0]+1e-64); 
    dspin[1] = (s1r[1]-p1.spin[1]-s1[1])/(s1[1]+1e-64); 
    dspin[2] = (s1r[2]-p1.spin[2]-s1[2])/(s1[2]+1e-64); 
    dspin[3] = (s2r[0]-p2.spin[0]-s2[0])/(s2[0]+1e-64); 
    dspin[4] = (s2r[1]-p2.spin[1]-s2[1])/(s2[1]+1e-64); 
    dspin[5] = (s2r[2]-p2.spin[2]-s2[2])/(s2[2]+1e-64); 
        

    std::cout<<std::setw(width)<<"diff";
    std::cout<<std::setw(width)<<dspin[0]
             <<std::setw(width)<<dspin[1]
             <<std::setw(width)<<dspin[2];
    std::cout<<std::setw(width)<<dspin[0]
             <<std::setw(width)<<dspin[1]
             <<std::setw(width)<<dspin[2];
    std::cout<<std::endl;

    Float diff_max = 0;
    for (int k=0; k<6; k++) {
        diff_max = std::max(std::fabs(dspin[k]), diff_max);
    }
    if (diff_max>ROUND_OFF_ERROR_LIMIT*10) {
        std::cerr<<"Test failed! difference "<<diff_max<<" > round off error "<<ROUND_OFF_ERROR_LIMIT*10<<std::endl;
        abort();
    }
    
    return 0;
}
