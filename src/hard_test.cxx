#include <iostream>
#include <cstdio>
#include <iomanip>
#include <string>
#include <unordered_map>
#include <particle_simulator.hpp>
#define HARD_DEBUG_PRINT_FEQ 1

#include <getopt.h>
#include "io.hpp"
#include "hard_assert.hpp"
#include "cluster_list.hpp"
#include "hard.hpp"
#include "soft_ptcl.hpp"
#include "static_variables.hpp"

#ifndef NAN_CHECK
#define NAN_CHECK(val) assert((val) == (val));
#endif

class IOParamsHardTest{
public:
    // IO parameters
    IOParamsContainer input_par_store;
    IOParams<PS::F64> time;
    IOParams<PS::F64> rin;
    IOParams<PS::F64> rout;
    IOParams<PS::F64> rsearch;
    IOParams<PS::F64> rbin;
    IOParams<PS::F64> dt_limit;
    IOParams<PS::F64> eta;
    IOParams<PS::F64> eps;
    IOParams<PS::S64> n;
#ifdef ORBIT_SAMPLING
    IOParams<PS::S64> n_split;
#endif
    IOParams<std::string> fname_inp;

    bool print_flag;

    IOParamsHardTest(): input_par_store(),
                        time (input_par_store, 0.0, "t", "finishing time"),
                        rin  (input_par_store, 0.0, "r-in", "changeover inner boundary"),
                        rout (input_par_store, 0.0, "r-out","changeover outer boundary"),
                        rsearch(input_par_store, 0.0, "r-search","search radius"),
                        rbin (input_par_store, 0.0, "r-bin","group radius"),
                        dt_limit(input_par_store, 0.0, "s", "maximum hermite time step"),
                        eta  (input_par_store, 0.0, "eta", "eta for hermite time step"),
                        eps  (input_par_store, 0.0, "eps", "softening length"),
                        n    (input_par_store, 0, "n", "number of particles"),
#ifdef ORBIT_SAMPLING
                        n_split (input_par_store, 4,  "n-split", "Number of binary sample points for tree perturbation force"),
#endif
                        fname_inp(input_par_store, "__NONE__", "snap-filename", "Input data file", NULL, false), 
                        print_flag(false) {}
    
    int read(int argc, char *argv[], const int opt_used_pre=0) {
        static int flag=-1;
        static struct option long_options[] = {
            {rin.key,     required_argument, &flag, 1},
            {rout.key,    required_argument, &flag, 2},
            {rsearch.key, required_argument, &flag, 3},
            {rbin.key,    required_argument, &flag, 4}, 
            {eta.key,     required_argument, &flag, 6}, 
            {eps.key,     required_argument, &flag, 7}, 
#ifdef ORBIT_SAMPLING
            {n_split.key, required_argument, &flag, 8}, 
#endif
            {"help",      no_argument, 0, 'h'},
            {0,0,0,0}
        };
        
        int opt_used = opt_used_pre;
        int copt;
        int option_index;
        optind = 0; // reset getopt
        while ((copt = getopt_long(argc, argv, "-t:s:n:h", long_options, &option_index)) != -1) 
            switch (copt) {
            case 0:
                switch (flag) {
                case 1:
                    rin.value = atof(optarg);
                    if(print_flag) rin.print(std::cout);
                    opt_used += 2;
                    break;
                case 2:
                    rout.value = atof(optarg);
                    if(print_flag) rout.print(std::cout);
                    opt_used += 2;
                    break;
                case 3:
                    rsearch.value = atof(optarg);
                    if(print_flag) rsearch.print(std::cout);
                    opt_used += 2;
                    break;
                case 4:
                    rbin.value = atof(optarg);
                    if(print_flag) rbin.print(std::cout);
                    opt_used += 2;
                    break;
                case 6:
                    eta.value = atof(optarg);
                    if(print_flag) eta.print(std::cout);
                    opt_used += 2;
                    break;
                case 7:
                    eps.value = atof(optarg);
                    if(print_flag) eps.print(std::cout);
                    opt_used += 2;
                    break;
#ifdef ORBIT_SAMPLING
                case 8:
                    n_split.value = atoi(optarg);
                    if(print_flag) n_split.print(std::cout);
                    opt_used += 2;
                    break;
#endif
                default:
                    break;
                }
                break;
            case 't':
                time.value = atof(optarg);
                if(print_flag) time.print(std::cout);
                opt_used += 2;
                break;
            case 's':
                dt_limit.value = atof(optarg);
                if(print_flag) dt_limit.print(std::cout);
                opt_used += 2;
                break;
            case 'n':
                n.value = atoi(optarg);
                if(print_flag) n.print(std::cout);
                opt_used += 2;
                break;
            case 'h':
                if(print_flag) {
                    std::cout<<"The tool to test hard integrator\n"
                             <<"Usage: petar.hard.test [options] [data filename]\n"
                             <<"       Data file content per line:\n";
                    ParticleBase::printTitleWithMeaning(std::cout,0,13);
                    std::cout<<"Options:\n";
                    input_par_store.printHelp(std::cout, 2, 10, 23);
                }
                return -1;
            case '?':
                opt_used +=2;
                break;
            default:
                break;
            }
        opt_used ++;
        //std::cout<<"Opt used:"<<opt_used<<std::endl;
        if (opt_used<argc) {
            fname_inp.value =argv[argc-1];
            if(print_flag) std::cout<<"Reading data file name: "<<fname_inp.value<<std::endl;
        }

        if(print_flag) std::cout<<"----- Finish reading input options -----\n";

        return opt_used-1;
    }

    //! check paramters
    bool checkParams() {
        assert(time.value>=0.0);
        assert(rin.value>0.0);
        assert(rout.value>rin.value);
        assert(rsearch.value>rout.value);
        assert(rbin.value>0.0&&rbin.value<rin.value);
        assert(eta.value>0.0);
        assert(eps.value>=0.0);
        assert(dt_limit.value>0.0);
        assert(n.value>0);
#ifdef ORBIT_SAMPLING
        assert(n_split.value>=0);
#endif
        return true;
    }    
};


int main(int argc, char** argv)
{
    IOParamsHardTest io;
    io.print_flag = true;
    int opt_used = io.read(argc,argv);
    if (opt_used==-1) return 0;
    io.checkParams();
    
    // open data file
    FILE* fin;
    if ( (fin = fopen(io.fname_inp.value.c_str(),"r")) == NULL) {
        fprintf(stderr,"Error: Cannot open input file %s.\n",io.fname_inp.value.c_str());
        abort();
    }

    PS::F64 &rsearch = io.rsearch.value;
    PS::F64 &rin = io.rin.value;
    PS::F64 &rout = io.rout.value;
    PS::F64 &rbin = io.rbin.value;
    PS::F64 &eta = io.eta.value;
    PS::F64 &eps = io.eps.value;
    PS::F64 &dt_limit = io.dt_limit.value;
    PS::F64 &time = io.time.value;
    PS::S64 &N = io.n.value;
#ifdef ORBIT_SAMPLING
    PS::S64 &n_split = io.n_split.value;
#endif
    
    ParticleBase pin;
    PS::ReallocatableArray<PS::S32> p_list;
    PS::ReallocatableArray<PS::S32> n_cluster;
    PS::ReallocatableArray<PS::S32> np;
    p_list.resizeNoInitialize(N);
    n_cluster.resizeNoInitialize(1);
    n_cluster[0] = N;

    PS::ParticleSystem<FPSoft> sys;

    PS::F64 m_average =0;
    for (int i=0; i<N; i++) {
        pin.readAscii(fin);
        ChangeOver co;
        GroupDataDeliver gdata;
        sys.addOneParticle(FPSoft(Ptcl(pin, rsearch, i+1, gdata, co), 0, i));
        p_list[i]=i;
        m_average = pin.mass;
    }
    m_average = pin.mass/N;
    Ptcl::r_search_min = rsearch;
    Ptcl::search_factor = 3;
    Ptcl::r_group_crit_ratio = rbin/rin;
    Ptcl::mean_mass_inv = 1.0/m_average;

    for (int i=0; i<N; i++) {
        sys[i].changeover.setR(sys[i].mass*Ptcl::mean_mass_inv, rin, rout);
        sys[i].calcRSearch(dt_limit);
    }

    PS::F64 time_sys = 0.0;
    PS::F64 dt_min_hard = 1.0;
    for (PS::S32 i=0;i<40;i++) dt_min_hard *= 0.5;

    // system hard paramters
    HardManager hard_manager;
    hard_manager.setDtRange(dt_limit, 40);
    hard_manager.setEpsSq(eps);
    hard_manager.setGravitationalConstant(1.0);
    hard_manager.r_in_base = rin;
    hard_manager.r_out_base = rout;
#ifdef HARD_CHECK_ENERGY
    hard_manager.energy_error_max = 1e-4;
#else
    hard_manager.energy_error_max = NUMERIC_FLOAT_MAX;
#endif
    hard_manager.n_step_per_orbit = 8;
    hard_manager.ap_manager.r_tidal_tensor = rbin;
    hard_manager.ap_manager.id_offset = N;
#ifdef ORBIT_SAMPLING
    hard_manager.ap_manager.orbit_manager.setParticleSplitN(input_parameters.n_split.value);
#endif
    hard_manager.h4_manager.step.eta_4th = eta;
    hard_manager.h4_manager.step.eta_2nd = 0.01*eta;
    hard_manager.h4_manager.step.calcAcc0OffsetSq(m_average, rout);
    hard_manager.ar_manager.energy_error_relative_max = 1e-8;
    hard_manager.ar_manager.step_count_max = 1e6;
    // set symplectic order
    hard_manager.ar_manager.step.initialSymplecticCofficients(-6);
    hard_manager.ar_manager.slowdown_pert_ratio_ref = 1e-4;
    hard_manager.ar_manager.slowdown_timescale_max = dt_limit;
#ifdef SLOWDOWN_MASSRATIO
    hard_manager.ar_manager.slowdown_mass_ref = m_average;
#endif
    hard_manager.ar_manager.interrupt_detection_option = 0;
  
    // check consistence of paramters
    hard_manager.checkParams();

    fclose(fin);

    SystemHard sys_hard;
    sys_hard.manager = &hard_manager;
    sys_hard.allocateHardIntegrator(1);
    sys_hard.setPtclForIsolatedMultiClusterOMP(sys, p_list, n_cluster);

    PS::S32 n_sys = sys.getNumberOfParticleLocal();
    sys_hard.findGroupsAndCreateArtificialParticlesOMP<PS::ParticleSystem<FPSoft>, FPSoft>(sys, dt_limit);
  
    PS::ReallocatableArray<PS::S32> mass_modify_list;

    // correct change over
    auto& hard_ptcl = sys_hard.getPtcl();
    int n =hard_ptcl.size();
    for (int i=0; i<n; i++) hard_ptcl[i].changeover.updateWithRScale();

    // recover mass
    for(int i=0; i<n; i++) {
        auto& pi_artificial = hard_ptcl[i].group_data.artificial;
        if(pi_artificial.isMember()){
            const PS::S64 cm_adr = hard_ptcl[i].getParticleCMAddress(); // notice status is negative 
#ifdef HARD_DEBUG
            assert(pi_artificial.getMassBackup()>0); 
            assert(cm_adr>0);
#endif
            hard_ptcl[i].mass = pi_artificial.getMassBackup();
        }
    }

    std::cout<<std::setprecision(WRITE_PRECISION);

    sys_hard.setTimeOrigin(time_sys);
    while(time_sys < time){
        fprintf(stderr,"Time = %e\n", time_sys+dt_limit);
        sys_hard.driveForMultiClusterOMP(dt_limit, &sys[0]);
        sys_hard.writeBackPtclForMultiCluster(sys, mass_modify_list);
        time_sys += dt_limit;
        sys.setNumberOfParticleLocal(n_sys);
        sys_hard.setTimeOrigin(time_sys);
        sys_hard.ARC_substep_sum = 0;
        sys_hard.findGroupsAndCreateArtificialParticlesOMP<PS::ParticleSystem<FPSoft>, FPSoft>(sys, dt_limit);
        for (int i=0; i<n; i++) hard_ptcl[i].changeover.updateWithRScale();
        // recover mass
        for(int i=0; i<n; i++) {
            auto& pi_artificial = hard_ptcl[i].group_data.artificial;
            if(pi_artificial.isMember()){
                const PS::S64 cm_adr= hard_ptcl[i].getParticleCMAddress(); // notice status is negative 
#ifdef HARD_DEBUG
                assert(cm_adr>0);
                assert(pi_artificial.getMassBackup()>0); 
#endif
                hard_ptcl[i].mass = pi_artificial.getMassBackup();
            }
        }
    } 

    return 0;
}

