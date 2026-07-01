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
#include "kickdriftstep.hpp"
#include "profile.hpp"

#ifndef NAN_CHECK
#define NAN_CHECK(val) assert((val) == (val));
#endif

class IOParamsHardTest{
public:
    // IO parameters
    IOParamsContainer input_par_store;
    IOParams<PS::F64> time;
    IOParams<PS::F64> r_in_over_out;
    IOParams<PS::F64> r_out;
    IOParams<PS::F64> r_search_min;
    IOParams<PS::F64> dt_soft;
    IOParams<PS::S64> n_bin;
    IOParams<PS::S64> write_style;
    IOParams<std::string> fname_inp;
    IOParams<std::string> fname_snap;

    bool print_flag;

    IOParamsHardTest(): input_par_store(),
                        time (input_par_store, 1.0, "t", "finishing time"),
                        r_in_over_out (input_par_store, 0.1, "r-ratio", "ratio of changeover inner over outer radius"),
                        r_out         (input_par_store, 0.1, "r",   "changeover outer boundary (r_out)"),
                        r_search_min  (input_par_store, 0, "r-search-min","search radius, if 0, use 1.5*r_out"),
                        dt_soft       (input_par_store, 0.125, "s", "soft time step (max hermite step), will be regularized to 2^-n"),
                        n_bin         (input_par_store, 0, "b", "number of binaries"),
                        write_style   (input_par_store, 2,   "w", "Data file writing style; 0: no output; 1: write snapshots separately; 2. write snapshots in one status file"),
                        fname_inp     (input_par_store, "__NONE__", "snap-filename", "Input data file", NULL, false), 
                        fname_snap    (input_par_store, "data", "f", "Prefix of filenames for output data: [prefix].**"),
                        print_flag(false) {}
    
    int read(int argc, char *argv[], const int opt_used_pre=0) {
        static int flag=-1;
        static struct option long_options[] = {
            {r_in_over_out.key, required_argument, &flag, 1},
            {r_search_min.key,  required_argument, &flag, 2},
            {"help",      no_argument, 0, 'h'},
            {0,0,0,0}
        };
        
        int opt_used = opt_used_pre;
        int copt;
        int option_index;
        optind = 0; // reset getopt
        while ((copt = getopt_long(argc, argv, "-t:b:r:s:w:f:h", long_options, &option_index)) != -1) 
            switch (copt) {
            case 0:
                switch (flag) {
                case 1:
                    r_in_over_out.value = atof(optarg);
                    if(print_flag) r_in_over_out.print(std::cout);
                    opt_used += 2;
                    break;
                case 2:
                    r_search_min.value = atof(optarg);
                    if(print_flag) r_search_min.print(std::cout);
                    opt_used += 2;
                    break;
                default:
                    break;
                }
                break;
            case 't':
                time.value = atof(optarg);
                if(print_flag) time.print(std::cout);
                opt_used += 2;
                break;
            case 'r':
                r_out.value = atof(optarg);
                if(print_flag) r_out.print(std::cout);
                opt_used += 2;
                break;
            case 's':
                dt_soft.value = atof(optarg);
                if(print_flag) dt_soft.print(std::cout);
                opt_used += 2;
                break;
            case 'b':
                n_bin.value = atoi(optarg);
                if(print_flag) n_bin.print(std::cout);
                opt_used += 2;
                break;
            case 'w':
                write_style.value = atoi(optarg);
                if(print_flag) write_style.print(std::cout);
                opt_used += 2;
                break;
            case 'f':
                fname_snap.value = optarg;
                if(print_flag) fname_snap.print(std::cout);
                opt_used += 2;
                break;
            case 'h':
                if(print_flag) {
                    std::cout<<"The tool to test hard integrator\n"
                             <<"Usage: petar.hard.test [options] [data filename]\n"
                             <<"       filename: initial or restart data (snapshot) filename\n"
                             <<"       Data file content:\n"
                             <<"            First line: \n"
                             <<"             1. File_ID: 0 for initialization, else for restarting\n"
                             <<"             2. N_particle \n"
                             <<"             3. Time\n"
                             <<"            Following lines:\n";
                    FPSoft::printTitleWithMeaning(std::cout,0,13);
                    std::cout<<"Options:\n";
                    input_par_store.printHelp(std::cout,true);
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
        assert(r_in_over_out.value>0.0 && r_in_over_out.value <1.0);
        assert(r_out.value>0.0);
        assert(r_search_min.value>=r_out.value);
        assert(dt_soft.value>0.0);
        return true;
    }    
};

class HardTestProfile{
public:
    Tprofile total;
    Tprofile findGroups;
    Tprofile integration;
    Tprofile output;
    const PS::S32 n_profile;

    HardTestProfile(): total         (Tprofile("Total      ")),
                       findGroups    (Tprofile("FindGroups ")),
                       integration   (Tprofile("Integration")),
                       output        (Tprofile("Output     ")),
                       n_profile(4) {} 

	void print(std::ostream & fout, const PS::F64 time_sys, const PS::S64 n_loop=1){
        fout<<"Time: "<<time_sys<<std::endl;
        
        for(PS::S32 i=0; i<n_profile; i++) {
            Tprofile* iptr = (Tprofile*)this+i;
            iptr->print(fout, n_loop);
        }
    }

    void dump(std::ostream & fout, const PS::S64 n_loop=1, const PS::S32 width=PROFILE_PRINT_WIDTH) const {
        for(PS::S32 i=0; i<n_profile; i++) {
            Tprofile* iptr = (Tprofile*)this+i;
            iptr->dump(fout, n_loop, width);
        }
    }

    void dumpName(std::ostream & fout, const PS::S32 width=PROFILE_PRINT_WIDTH) const {
        for(PS::S32 i=0; i<n_profile; i++) {
            Tprofile* iptr = (Tprofile*)this+i;
            iptr->dumpName(fout, width);
        }
    }

    void clear(){
        for(PS::S32 i=0; i<n_profile; i++) {
            Tprofile* iptr = (Tprofile*)this+i;
            iptr->reset();
        }
    }

} profile;

class HardCounts{
public:
    NumCounter sdar_substep_sum;
    NumCounter sdar_tsyn_step_sum;
    NumCounter H4_step_sum;
    NumCounter H4_force_sum;
    NumCounter n_neighbor_zero;
    NumCounter sdar_n_groups;
    NumCounter sdar_n_groups_iso;
    NumCounter sdar_n_groups_new;
    NumCounter sdar_n_groups_end;
    NumCounter sdar_n_groups_arti_change;
    const PS::S32 n_counter;

    HardCounts():
                 sdar_substep_sum  (NumCounter("AR_step_sum")),
                 sdar_tsyn_step_sum(NumCounter("AR_tsyn_sum")),
                 H4_step_sum       (NumCounter("H4_step_sum")),
                 H4_force_sum      (NumCounter("H4_int_sum ")),
                 n_neighbor_zero   (NumCounter("H4_no_NB   ")),
                 sdar_n_groups     (NumCounter("Stab_Ngroup")),
                 sdar_n_groups_iso (NumCounter("Iso_Ngroup ")),
                 sdar_n_groups_new (NumCounter("Form_Ngroup")),
                 sdar_n_groups_end (NumCounter("End_Ngroup ")),
                 sdar_n_groups_arti_change(NumCounter("Modf_Ngroup")),
                 n_counter(10) 
                 {}

    void dump(std::ostream & fout, const PS::S64 n_loop=1, const PS::S32 print_part=0, const PS::S32 width=PROFILE_PRINT_WIDTH) const{
        for(PS::S32 i=0; i<n_counter; ++i) {
            NumCounter* iptr = (NumCounter*)this+i;
            iptr->dump(fout, n_loop, width);
        }
    }

    void dumpName(std::ostream & fout, const PS::S32 print_part=0, const PS::S32 width=PROFILE_PRINT_WIDTH) const{
        for(PS::S32 i=0; i<n_counter; ++i) {
            NumCounter* iptr = (NumCounter*)this+i;
            iptr->dumpName(fout, width);
        }
    }
    
    void clear() {
        for(PS::S32 i=0; i<n_counter; ++i) {
            NumCounter* iptr = (NumCounter*)this+i;
            *iptr = 0;
        }
    }                     


} n_count;

int main(int argc, char** argv)
{
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL    
    omp_set_nested(1);
#endif

    // IO parameters    
    std::vector<IOParamsContainer*> all_pars;
        
    IOParamsHardTest main_parameters;
    IOParamsHard hard_parameters;
    all_pars.push_back(&main_parameters.input_par_store);
    all_pars.push_back(&hard_parameters.input_par_store);

    std::vector<std::string> known_options;
    known_options.push_back("help");
    known_options.push_back("h");
    FindUndefinedOptions(all_pars, argc, argv, &known_options);
    
    main_parameters.print_flag = true;
    hard_parameters.print_flag = true;
    opterr = 0; // do not print getopt error messages
    int opt_used = main_parameters.read(argc,argv);
    hard_parameters.read(argc,argv,false,true);

    if (opt_used==-1) return 0;

    // get parameters    
    PS::F64 &r_search_min = main_parameters.r_search_min.value;
    PS::F64 &r_in_over_out = main_parameters.r_in_over_out.value;
    PS::F64 &r_out = main_parameters.r_out.value;
    PS::F64 r_in = r_out*r_in_over_out;
    PS::F64 &dt_soft = main_parameters.dt_soft.value;
    PS::F64 &time = main_parameters.time.value;

    if (r_search_min==0) r_search_min = 1.5*r_out;
    dt_soft = regularTimeStep(dt_soft);

    main_parameters.checkParams();
    
    // Particle system    
    FileHeader file_header;
    PS::ParticleSystem<FPSoft> sys;
    
    // reading data    
    sys.readParticleAscii(main_parameters.fname_inp.value.c_str(), file_header);
    PS::S32 N = sys.getNumberOfParticleLocal();
    hard_parameters.id_offset.value = N+1;

    // initial particle parameters    
    PS::F64 m_average =0;
    for (int i=0; i<N; i++) m_average += sys[i].mass;
    m_average /= N;

    Ptcl::mean_mass_inv = 1.0/m_average;
    Ptcl::r_search_min = r_search_min;
    Ptcl::search_factor = 3;

    // initial changeover radii    
    PS::S32 n_bin = main_parameters.n_bin.value;
    for (int i=0; i<2*n_bin; i+=2) {
        sys[i].changeover.setR(sys[i].mass*Ptcl::mean_mass_inv, r_in, r_out);
        sys[i+1].changeover.setR(sys[i+1].mass*Ptcl::mean_mass_inv, r_in, r_out);
        PS::F64vec vave = (sys[i].mass*sys[i].vel + sys[i+1].mass*sys[i+1].vel)/(sys[i].mass+sys[i+1].mass);
        PS::F64 v = sqrt(vave[0]*vave[0]+vave[1]*vave[1]+vave[2]*vave[2]);
        sys[i].r_search = std::max(r_search_min, v*dt_soft*Ptcl::search_factor + sys[i].changeover.getRout());
        sys[i+1].r_search = std::max(r_search_min, v*dt_soft*Ptcl::search_factor + sys[i+1].changeover.getRout());
    }
    for (int i=2*n_bin; i<N; i++) {
        sys[i].changeover.setR(sys[i].mass*Ptcl::mean_mass_inv, r_in, r_out);
        sys[i].calcRSearch(dt_soft);
    }

    // status
    PS::F64 time_sys = 0.0;
    Status stat;
    stat.time = time_sys;
    stat.n_real_loc = N;
    stat.n_real_glb = N;
    
    // set r_group
    // dt_max/2^20 as period
    //if (hard_parameters.r_group.value==-1) {
    //    PS::F64 r_group = COMM::Binary::periodToSemi(dt_soft/1024, m_average, hard_parameters.gravitational_constant.value);
    //    hard_parameters.r_group.value = r_group;
    //    hard_parameters.r_search_group.value = r_group*1.5;
    //}

    // system hard paramters
    HardManager hard_manager;
    hard_manager.initial(hard_parameters, stat, true, false, m_average, 0.0, r_out, r_in_over_out, dt_soft);

    // print parameters
    std::cout<<"r_out="<<r_out<<std::endl
             <<"r_in="<<r_in<<std::endl
             <<"r_search_min="<<r_search_min<<std::endl  
             <<"dt_soft="<<dt_soft<<std::endl
             <<"r_group="<<hard_parameters.r_group.value<<std::endl
             <<"r_search_group="<<hard_parameters.r_search_group.value<<std::endl
             <<"hermite_dt_max="<<hard_manager.h4_manager.step.getDtMax()<<std::endl
             <<"hermite_dt_min="<<hard_manager.h4_manager.step.getDtMin()<<std::endl;

    // system hard    
    SystemHard sys_hard;
    sys_hard.manager = &hard_manager;
    sys_hard.setTimeOrigin(time_sys);

#ifdef HARD_DEBUG_PRINT
    sys_hard.output_filename_prefix = main_parameters.fname_snap.value;
#endif

    // initial system hard particles, only one cluster
    PS::ReallocatableArray<PS::S32> p_list;
    p_list.resizeNoInitialize(N);
    for (int i=0; i<N; i++) p_list[i]=i;

    PS::ReallocatableArray<PS::S32> n_cluster;
    n_cluster.resizeNoInitialize(1);
    n_cluster[0] = N;

    sys_hard.setPtclForIsolatedMultiClusterOMP(sys, p_list, n_cluster);

    auto& hard_ptcl = sys_hard.getPtcl();
    int n_hard = hard_ptcl.size();

    // output style    
    //std::cout<<std::setprecision(WRITE_PRECISION);
    PS::S32 write_style = main_parameters.write_style.value;

    std::ofstream fstatus;
    if (write_style==2) {
        std::string fname = main_parameters.fname_snap.value + ".status";    
        fstatus.open(fname.c_str(),std::ofstream::out);
    }

    // integration loop    
    PS::ReallocatableArray<PS::S32> mass_modify_list;
    while(time_sys < time){

        // find groups
        fprintf(stderr,"Time = %e\n", time_sys+dt_soft);
        profile.total.start();

        // update artificial particles
        profile.findGroups.start();        
        sys_hard.findGroupsAndCreateArtificialParticlesOMP<PS::ParticleSystem<FPSoft>, FPSoft>(sys, dt_soft);

        // correct change over
        for (int i=0; i<n_hard; i++) hard_ptcl[i].changeover.updateWithRScale();
        
        // recover mass
        for(int i=0; i<n_hard; i++) {
            auto& pi_artificial = hard_ptcl[i].group_data.artificial;
            if(pi_artificial.isMember()){
#ifdef HARD_DEBUG
                const PS::S64 cm_adr= hard_ptcl[i].getParticleCMAddress(); // notice status is negative 
                assert(cm_adr>0);
                assert(pi_artificial.getMassBackup()>0); 
#endif
                hard_ptcl[i].mass = pi_artificial.getMassBackup();
            }
        }
        profile.findGroups.barrier();
        profile.findGroups.end();

        // integration
        profile.integration.start();
        sys_hard.driveForMultiClusterOMP(dt_soft, &sys[0]);
        profile.integration.barrier();
        profile.integration.end();

#ifdef PROFILE
        n_count.sdar_substep_sum  += sys_hard.sdar_substep_sum;
        n_count.sdar_tsyn_step_sum+= sys_hard.sdar_tsyn_step_sum;
        n_count.H4_step_sum       += sys_hard.H4_step_sum;
        n_count.H4_force_sum      += sys_hard.H4_force_sum;
        n_count.sdar_n_groups     += sys_hard.sdar_n_groups;
        n_count.sdar_n_groups_iso += sys_hard.sdar_n_groups_iso;
        
        sys_hard.sdar_substep_sum =0;
        sys_hard.sdar_tsyn_step_sum=0;
        sys_hard.H4_step_sum =0;
        sys_hard.H4_force_sum=0;
        sys_hard.sdar_n_groups=0;
        sys_hard.sdar_n_groups_iso=0;
#endif
#ifdef HARD_COUNT_NO_NEIGHBOR
        n_count.n_neighbor_zero  += sys_hard.n_neighbor_zero;
        sys_hard.n_neighbor_zero =0;
#endif
        n_count.sdar_n_groups_new += sys_hard.sdar_n_groups_new;
        n_count.sdar_n_groups_end += sys_hard.sdar_n_groups_end;
        n_count.sdar_n_groups_arti_change += sys_hard.sdar_n_groups_arti_change;
        sys_hard.sdar_n_groups_new =0;
        sys_hard.sdar_n_groups_end =0;
        sys_hard.sdar_n_groups_arti_change =0;

        profile.output.start();        
        sys_hard.writeBackPtclForMultiCluster(sys, mass_modify_list);

        // reset number, remove artificial particles        
        sys.setNumberOfParticleLocal(N);

        // update time        
        time_sys += dt_soft;
        sys_hard.setTimeOrigin(time_sys);        
        stat.time = time_sys;

        // output
        if (write_style==1) {
            file_header.time = stat.time;
            file_header.nfile++;
            std::string fname = main_parameters.fname_snap.value + "." + std::to_string(file_header.nfile);    
            sys.writeParticleAscii(fname.c_str(), file_header);
        } 
        else if (write_style==2) {
            // write snapshot with one line
            stat.printColumnAscii(fstatus, WRITE_WIDTH);
            for (int i=0; i<stat.n_real_loc; i++) sys[i].printColumnAscii(fstatus, WRITE_WIDTH);
            fstatus<<std::endl;
        }    
        profile.output.barrier();
        profile.output.end();

        profile.total.barrier();
        profile.total.end();

        profile.dumpName(std::cout);
        std::cout<<std::endl;
        profile.dump(std::cout, 1);
        std::cout<<std::endl;

        profile.clear();        

        n_count.dumpName(std::cout);
        std::cout<<std::endl;
        n_count.dump(std::cout,1);
        std::cout<<std::endl;

        n_count.clear();
    } 

    return 0;
}
