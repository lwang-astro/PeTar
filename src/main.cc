#ifdef P3T_64BIT
#define CALC_EP_64bit
#define CALC_SP_64bit
#define RSQRT_NR_EPJ_X4
#define RSQRT_NR_SPJ_X4
#elif P3T_MIXBIT
#define CALC_EP_64bit
#define RSQRT_NR_EPJ_X4
#else
#define RSQRT_NR_EPJ_X2
#endif //P3T_64BIT

#ifdef FORCE_CHECK
#define FORCE_DIRECT
#endif

#if defined(INTRINSIC_K) || defined(INTRINSIC_X86)
#define INTRINSIC
#endif

#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<sstream>
//#include<unistd.h>
#include<getopt.h>
#ifdef USE_C03
#include<map>
#else //USE_C03
#include<unordered_map>
#endif //USE_C03
#include<particle_simulator.hpp>
#include"soft.hpp"
#include"hard.hpp"
#include"kepler.hpp"
#include"io.hpp"
#include"init.hpp"
#include"integrate.hpp"
#include"domain.hpp"
#include"AR.h" /// include AR.h (L.Wang)
//#include"cluster.hpp"
#include"cluster_list.hpp"
#ifdef PROFILE
#include"profile.hpp"
#endif

#ifdef USE_QUAD
typedef PS::TreeForForceLong<ForceSoft, EPISoft, EPJSoft>::QuadrupoleWithSymmetrySearch TreeForce; 
#else
typedef PS::TreeForForceLong<ForceSoft, EPISoft, EPJSoft>::MonopoleWithSymmetrySearch TreeForce;
//typedef PS::TreeForForceLong<ForceSoft, EPISoft, EPJSoft>::MonopoleWithScatterSearch Tree;
#endif
typedef PS::ParticleSystem<FPSoft> SystemSoft;

// For neighbor searching
typedef PS::TreeForForceShort<ForceSoft, EPISoft, EPJSoft>::Symmetry TreeNB;

int main(int argc, char *argv[]){
    std::cout<<std::setprecision(PRINT_PRECISION);
    std::cerr<<std::setprecision(PRINT_PRECISION);
    PS::Initialize(argc, argv);

    PS::S32 my_rank = PS::Comm::getRank();
    const PS::S32 n_proc = PS::Comm::getNumberOfProc();

#ifdef PROFILE
    if(my_rank==0) {
        std::cout<<"----- Parallelization information -----\n";
        std::cout<<"MPI processors: "<<n_proc<<std::endl;
        std::cout<<"OMP threads:    "<<PS::Comm::getNumberOfThread()<<std::endl;
    }
    
    SysProfile profile;
    SysCounts  n_count;
    SysCounts  n_count_sum;
    PsProfile  ps_profile;

    std::ofstream fprofile;
    std::string rank_str;
    std::stringstream atmp;
    atmp<<my_rank;
    atmp>>rank_str;
    PS::S64 dn_loop = 0;

#endif
    
    // initial parameters
    IOParams<PS::F64> ratio_r_cut  (0.1,  "r_in / r_out");
    IOParams<PS::F64> theta        (0.3,  "Openning angle theta");
    IOParams<PS::S32> n_leaf_limit (20,   "Tree leaf number limit", "optimized value shoudl be slightly >=11+N_bin_sample (20)");
#ifdef USE__AVX512
    IOParams<PS::S32> n_group_limit(1024, "Tree group number limit", "optimized for x86-AVX512 (2048)");    
#else
    IOParams<PS::S32> n_group_limit(512,  "Tree group number limit", "optimized for x86-AVX2 (512)");
#endif
    IOParams<PS::S32> n_smp_ave    (100,  "Average target number of sample particles per process");
    IOParams<PS::S32> n_split      (8,    "Number of binary sample points for tree perturbation force");
    IOParams<PS::S64> n_bin        (0,    "Number of primordial binaries (assume binaries ID=1,2*n_bin)");
    IOParams<PS::F64> time_end     (10.0, "Finishing time");
    IOParams<PS::F64> eta          (0.1,  "Hermite time step coefficient eta");
    IOParams<PS::S64> n_glb        (16384,"Total number of particles (this will suppress reading snapshot data and use Plummer model generator without binary)");
    IOParams<PS::F64> dt_soft      (0.0,  "Tree timestep","0.1*r_out/sigma_1D");
    IOParams<PS::F64> dt_snp       (0.0625,"Output time interval of particle dataset");
    IOParams<PS::F64> search_factor(1.0,  "Neighbor searching coefficient");
    IOParams<PS::F64> dt_limit_hard_factor(4.0,  "Limit of tree time step/hard time step");
    IOParams<PS::S32> dt_min_hermite_index(40,   "Power index n for the smallest time step (0.5^n) allowed in Hermite integrator");
    IOParams<PS::S32> dt_min_arc_index    (64,   "Power index n for the smallest time step (0.5^n) allowed in ARC integrator");
    IOParams<PS::F64> dt_err_pert  (1e-6, "Time synchronization maximum (relative) error for perturbed ARC integrator");
    IOParams<PS::F64> dt_err_soft  (1e-3, "Time synchronization maximum (relative) error for no-perturber (only soft perturbation) ARC integrator");
    IOParams<PS::F64> e_err_arc    (1e-10,"Maximum energy error allown for ARC integrator");
    IOParams<PS::F64> eps          (0.0,  "Softerning eps");
    IOParams<PS::F64> r_out        (0.0,  "Transit function outer boundary radius", "<m>/sigma_1D^2*/ratio_r_cut");
    IOParams<PS::F64> r_bin        (0.0,  "Maximum binary radius criterion", "0.8*r_in");
    IOParams<PS::F64> sd_factor    (1e-6, "Slowdown perturbation criterion");
    IOParams<PS::S32> data_format  (1,    "Data read(r)/write(w) format BINARY(B)/ASCII(A): Writing off: r-B(5), r-A(4); Writing on: r-B/w-A (3); r-A/w-B (2); rw-A (1); rw-B (0)");
    IOParams<std::string> fname_snp("data","Prefix filename of dataset: [prefix].[File ID]");

    // PS::F64 g_min = 1e-6;

    // reading parameters
    bool reading_flag=true;
    bool app_flag=false; // appending data flag

    static struct option long_options[] = {
        {"n-split", required_argument, 0, 0},           //0 
        {"search-factor", required_argument, 0, 0},     //1 
        {"dt-max-factor", required_argument, 0, 0},     //2
        {"dt-min-hermite", required_argument, 0, 0},    //3
        {"dt-min-arc", required_argument, 0, 0},        //4
        {"dt-err-pert", required_argument, 0, 0},       //5
        {"dt-err-soft", required_argument, 0, 0},       //6
        {"energy-err-arc", required_argument, 0, 0},    //7
        {"soft-eps", required_argument, 0, 0},          //8
        {"slowdown-factor", required_argument, 0, 0},   //9
        {"help",no_argument, 0, 'h'},                   //10
        {0,0,0,0}
    };

    int copt;
    int option_index;
    if(my_rank == 0) std::cout<<"----- input option -----\n";
    while ((copt = getopt_long(argc, argv, "i:at:D:o:r:R:b:B:N:n:l:s:T:E:f:h", long_options, &option_index)) != -1) 
        switch (copt) {
        case 0:
            switch(option_index){
            case 0:
                n_split.value = atoi(optarg);
                if(my_rank == 0) n_split.print(std::cout);
                assert(n_split.value>=8);
                break;
            case 1:
                search_factor.value = atof(optarg);
                if(my_rank == 0) search_factor.print(std::cout);
                assert(search_factor.value>0.0);
                break;
            case 2:
                dt_limit_hard_factor.value = atof(optarg);
                if(my_rank == 0) dt_limit_hard_factor.print(std::cout);
                assert(dt_limit_hard_factor.value > 0.0);
                break;
            case 3:
                dt_min_hermite_index.value = atoi(optarg);
                if(my_rank == 0) dt_min_hermite_index.print(std::cout);
                assert(dt_min_hermite_index.value > 0);
                break;
            case 4:
                dt_min_arc_index.value = atoi(optarg);
                if(my_rank == 0) dt_min_arc_index.print(std::cout);
                assert(dt_min_arc_index.value > 0);
                break;
            case 5:
                dt_err_pert.value = atof(optarg);
                if(my_rank == 0) dt_err_pert.print(std::cout);
                assert(dt_err_pert.value > 0.0);
                break;
            case 6:
                dt_err_soft.value = atof(optarg);
                if(my_rank == 0) dt_err_soft.print(std::cout);
                assert(dt_err_soft.value > 0.0);
                break;
            case 7:
                e_err_arc.value = atof(optarg);
                if(my_rank == 0) e_err_arc.print(std::cout);
                assert(e_err_arc.value > 0.0);
                break;
            case 8:
                eps.value = atof(optarg);
                if(my_rank == 0) eps.print(std::cout);
                assert(eps.value>=0.0);
                break;
            case 9:
                sd_factor.value = atof(optarg);
                if(my_rank == 0) sd_factor.print(std::cout);
                assert(sd_factor.value>0.0);
                break;
            default:
                if(my_rank == 0) std::cerr<<"Unknown option. check '-h' for help.\n";
                abort();
            }
            break;
        case 'i':
            data_format.value = atoi(optarg);
            if(my_rank == 0) data_format.print(std::cout);
            assert(data_format.value>=0||data_format.value<=3);
            break;
        case 'a':
            app_flag=true;
            break;
        case 'b':
            r_bin.value = atof(optarg);
            if(my_rank == 0) r_bin.print(std::cout);
            assert(r_bin.value>0.0);
            break;
        case 'B':
            n_bin.value = atoi(optarg);
            if(my_rank == 0) n_bin.print(std::cout);
            assert(n_bin.value>=0);
            break;
        case 'T':
            theta.value = atof(optarg);
            if(my_rank == 0) theta.print(std::cout);
            assert(theta.value>=0.0);
            break;
        case 't':
            time_end.value = atof(optarg);
            if(my_rank == 0) time_end.print(std::cout);
            assert(time_end.value>=0.0);
            break;
        case 'E':
            eta.value = atof(optarg);
            if(my_rank == 0) eta.print(std::cout);
            assert(eta.value>0.0);
            break;
        case 'n':
            n_group_limit.value = atoi(optarg);
            if(my_rank == 0) n_group_limit.print(std::cout);
            assert(n_group_limit.value>0);
            break;
        case 'N':
            reading_flag=false;
            n_glb.value = atol(optarg);
            if(my_rank == 0) n_glb.print(std::cout);
            assert(n_glb.value>0);
            break;
        case 's':
            n_smp_ave.value = atoi(optarg);
            if(my_rank == 0) n_smp_ave.print(std::cout);
            assert(n_smp_ave.value>0.0);
            break;
        case 'D':
            dt_soft.value = atof(optarg);
            if(my_rank == 0) dt_soft.print(std::cout);
            assert(dt_soft.value>0.0);
            break;
        case 'o':
            dt_snp.value = atof(optarg);
            if(my_rank == 0) dt_snp.print(std::cout);
            assert(dt_snp.value>0.0);
            break;
        case 'l':
            n_leaf_limit.value = atoi(optarg);
            if(my_rank == 0) n_leaf_limit.print(std::cout);
            assert(n_leaf_limit.value>0);
            break;
        case 'r':
            ratio_r_cut.value = atof(optarg);
            if(my_rank == 0) ratio_r_cut.print(std::cout);
            assert(ratio_r_cut.value>0.0);
            assert(ratio_r_cut.value<1.0);
            break;
        case 'R':
            r_out.value = atof(optarg);
            if(my_rank == 0) r_out.print(std::cout);
            assert(r_out.value>0.0);
            break;
        case 'f':
            fname_snp.value = optarg;
            if(my_rank == 0) fname_snp.print(std::cout);
            break;
        case 'h':
            if(my_rank == 0){
                std::cout<<"Usage: nbody.out [option] [filename]"<<std::endl;
                std::cout<<"       Option defaulted values are shown after ':'\n"<<std::endl;
                std::cout<<"  -i: [I] "<<data_format<<std::endl;
                std::cout<<"          File content:\n"
                         <<"            First line: \n"
                         <<"             1. File_ID: 0 for initialization, else for restarting\n"
                         <<"             2. N_particle \n"
                         <<"             3. N_offset: for naming special particle ID, should be > N_particle\n"
                         <<"             4. Time\n"
                         <<"             5. Tree time step (set 0 for input file)\n"
                         <<"             6. n_split (set 0 for input file)\n"
                         <<"            Following lines:\n"
                         <<"             1. mass\n"
                         <<"             2. position[3]\n"
                         <<"             5. velocity[3]\n"
                         <<"             8. r_search (0.0)\n"
                         <<"             9. mass_backup (0.0)\n"
                         <<"             10. ID (>0,unique)\n"
                         <<"             11. status (0)\n"
                         <<"             12. Acceleration[3] (0.0)\n"
                         <<"             15. Potential (0.0)\n"
                         <<"             16. N_neighbor (0)\n"
                         <<"             (*) show initialization values which should be used together with FILE_ID = 0"<<std::endl;
                std::cout<<"  -a:     data output style (except snapshot) becomes appending, defaulted: replace"<<std::endl;
                std::cout<<"  -t: [F] "<<time_end<<std::endl;
                std::cout<<"  -D: [F] "<<dt_soft<<std::endl;
                std::cout<<"  -o: [F] "<<dt_snp<<std::endl;
                std::cout<<"        --dt-max-factor:   [F] "<<dt_limit_hard_factor<<std::endl;
                std::cout<<"        --dt-min-hermite:  [I] "<<dt_min_hermite_index<<std::endl;
                std::cout<<"        --dt-min-arc:      [I] "<<dt_min_arc_index<<std::endl;
                std::cout<<"        --dt-err-pert:     [F] "<<dt_err_pert<<std::endl;
                std::cout<<"        --dt-err-soft:     [F] "<<dt_err_soft<<std::endl;
                std::cout<<"  -r: [F] "<<ratio_r_cut<<std::endl;
                std::cout<<"  -R: [F] "<<r_out<<std::endl;
                std::cout<<"  -b: [F] "<<r_bin<<std::endl;
                std::cout<<"  -B: [I] "<<n_bin<<std::endl;
                std::cout<<"  -N: [I] "<<n_glb<<std::endl;
                std::cout<<"  -n: [I] "<<n_group_limit<<std::endl;
                std::cout<<"  -l: [I] "<<n_leaf_limit<<std::endl;
                std::cout<<"  -s: [I] "<<n_smp_ave<<std::endl;
                std::cout<<"        --n-split:         [I] "<<n_split<<std::endl;
                std::cout<<"  -T: [F] "<<theta<<std::endl;
                std::cout<<"  -E: [F] "<<eta<<std::endl;
                std::cout<<"        --search-factor:   [F] "<<search_factor<<std::endl;
                std::cout<<"        --energy-err-arc:  [F] "<<e_err_arc<<std::endl;
                std::cout<<"        --slowdown-factor: [F] "<<sd_factor<<std::endl;
                std::cout<<"        --soft-eps:        [F] "<<eps<<std::endl;
                std::cout<<"  -f: [S] "<<fname_snp<<std::endl;
                std::cout<<"  -h(--help):               print help"<<std::endl;
                std::cout<<"*** PS: r_in : transit function inner boundary radius\n"
                         <<"        r_out: transit function outer boundary radius\n"
                         <<"        sigma: half-mass radius velocity dispersion\n"
                         <<"        n_bin: number of primordial binaries\n"
                         <<"        <m>  : averaged mass"<<std::endl;
            }
            PS::Finalize();
            return 0;
        }
    
    PS::F64 time_sys = 0.0;
    PS::S32 n_loc;

#ifdef PROFILE
    std::string fproname=fname_snp.value+".prof.rank."+rank_str;
    if(app_flag) fprofile.open(fproname.c_str(),std::ofstream::out|std::ofstream::app);
    else         fprofile.open(fproname.c_str(),std::ofstream::out);
#endif
    
    SystemSoft system_soft;
    system_soft.initialize();
    system_soft.setAverageTargetNumberOfSampleParticlePerProcess(n_smp_ave.value);
    
    FileHeader file_header;
    if (reading_flag) {
      char* sinput=argv[argc-1];
      if(data_format.value==1||data_format.value==2||data_format.value==4)
          system_soft.readParticleAscii(sinput, file_header);
      else
          system_soft.readParticleBinary(sinput, file_header);
      time_sys = file_header.time;
      PS::Comm::broadcast(&time_sys, 1, 0);
      PS::Comm::broadcast(&file_header, 1, 0);
      n_glb.value = system_soft.getNumberOfParticleGlobal();
      n_loc = system_soft.getNumberOfParticleLocal();
      //      for(PS::S32 i=0; i<n_loc; i++) system_soft[i].id = i;
      if(my_rank==0)
          std::cout<<"Reading file "<<sinput<<std::endl
                   <<"N_tot = "<<n_glb.value<<"\nN_loc = "<<n_loc<<std::endl;
    }
    else {
        SetParticlePlummer(system_soft, n_glb.value, n_loc);
        file_header.nfile = 0;
        time_sys = file_header.time = 0.0;
        file_header.n_body = n_glb.value;
        file_header.id_offset = n_glb.value;
    }

    bool restart_flag = file_header.nfile; // nfile = 0 is assumed as initial data file

    PS::Comm::barrier();

    // status information output
    std::ofstream fstatus;
    if(my_rank==0) {
        if(app_flag) fstatus.open((fname_snp.value+".status").c_str(),std::ofstream::out|std::ofstream::app);
        else         fstatus.open((fname_snp.value+".status").c_str(),std::ofstream::out);
        fstatus<<std::setprecision(WRITE_PRECISION);
    }

    Status stat;

    // tree time step and n_split
    if (restart_flag) {
        if(dt_soft.value!=file_header.dt_soft&&dt_soft.value>0) 
            std::cerr<<"Warning: tree time step cannot be changed for restarting, the value from the data file ("<<file_header.dt_soft<<") will be used\n";
        if(n_split.value!=file_header.n_split)
            std::cerr<<"Warning: n_split cannot be changed for restarting, the value from the data file ("<<file_header.n_split<<") will be used\n";
        dt_soft.value = file_header.dt_soft;
        n_split.value = file_header.n_split;
    }
    else {
        file_header.dt_soft = dt_soft.value;
        file_header.n_split = n_split.value;        
        if(my_rank==0&&app_flag==false) {
            stat.dumpName(fstatus,WRITE_WIDTH);
            fstatus<<std::endl;
        }
    }
    
    PS::F64 r_in, m_average, v_disp, r_search_min;
    GetR(system_soft, r_in, r_out.value, r_bin.value, r_search_min, m_average, dt_soft.value, v_disp, search_factor.value, ratio_r_cut.value, n_bin.value, restart_flag);

//    EPISoft::r_out = r_out;
    EPISoft::r_in  = r_in;
    EPISoft::eps   = eps.value;
    EPISoft::r_out = EPJSoft::r_out = FPSoft::r_out = r_out.value;
    Ptcl::search_factor = search_factor.value;
    Ptcl::r_search_min = r_search_min;
    Ptcl::mean_mass_inv = 1.0/m_average;
//    const PS::F64 r_oi_inv = 1.0/(r_out - r_in);
//    EPJSoft::r_search_min = r_out*search_factor;
//    EPJSoft::m_average = m_average;
  
    // ID safety check
    if (!restart_flag) {
        for (PS::S32 i=0; i<n_loc; i++) {
            PS::S64 id = system_soft[i].id;
            if(id<=0) {
                std::cerr<<"Error: for initial data, the id should always larger than zero. current index i = "<<i<<", id = "<<system_soft[i].id<<"!"<<std::endl;
                abort();
            }

            // for binary, research depend on v_disp
            //if(id<=2*n_bin.value) system_soft[i].r_search = std::max(r_search_min*std::sqrt(system_soft[i].mass*Ptcl::mean_mass_inv),v_disp*dt_soft.value*search_factor.value);
            if(id<=2*n_bin.value) system_soft[i].r_search = std::max(r_search_min,v_disp*dt_soft.value*search_factor.value);
            else system_soft[i].calcRSearch(dt_soft.value);
        }
    }
    
    if(my_rank == 0) {
        std::cout<<"----- Parameter list: -----\n";
        std::cout<<" m_average    = "<<m_average      <<std::endl
                 <<" r_in         = "<<r_in           <<std::endl
                 <<" r_out        = "<<r_out.value    <<std::endl
                 <<" r_bin        = "<<r_bin.value    <<std::endl
                 <<" r_search_min = "<<r_search_min   <<std::endl
                 <<" vel_disp     = "<<v_disp         <<std::endl
                 <<" dt_soft      = "<<dt_soft.value  <<std::endl;
    }

// Finish dataset initialization


    const PS::F32 coef_ema = 0.2;
    PS::DomainInfo dinfo;
    PS::F32 domain_decompose_weight=1.0;
    dinfo.initialize(coef_ema);
    dinfo.decomposeDomainAll(system_soft);

    PS::F64ort * pos_domain = new PS::F64ort[n_proc];
    for(PS::S32 i=0; i<n_proc; i++) pos_domain[i] = dinfo.getPosDomain(i);

    system_soft.exchangeParticle(dinfo); 

    n_loc = system_soft.getNumberOfParticleLocal();

#pragma omp parallel for
    for(PS::S32 i=0; i<n_loc; i++){
        system_soft[i].rank_org = my_rank;
        system_soft[i].adr = i;
    }

    TreeNB tree_nb;
    tree_nb.initialize(n_glb.value, theta.value, n_leaf_limit.value, n_group_limit.value);
    TreeForce tree_soft;
    tree_soft.initialize(n_glb.value, theta.value, n_leaf_limit.value, n_group_limit.value);

    SystemHard system_hard_one_cluster;
    PS::F64 dt_limit_hard = dt_soft.value/dt_limit_hard_factor.value;
    PS::F64 dt_min_hermite = 1.0;
    PS::F64 dt_min_arc = 1.0;
    for (PS::S32 i=0;i<dt_min_hermite_index.value;i++) dt_min_hermite *= 0.5;
    for (PS::S32 i=0;i<dt_min_arc_index.value;i++) dt_min_arc *= 0.5;
    system_hard_one_cluster.setParam(r_bin.value, r_out.value, r_in, eps.value, dt_limit_hard, dt_min_hermite, eta.value, time_sys, sd_factor.value, file_header.id_offset, n_split.value);
    // system_hard_one_cluster.setARCParam();
    SystemHard system_hard_isolated;
    system_hard_isolated.setParam(r_bin.value, r_out.value, r_in, eps.value,  dt_limit_hard, dt_min_hermite, eta.value, time_sys, sd_factor.value, file_header.id_offset, n_split.value);
    system_hard_isolated.setARCParam(e_err_arc.value, dt_err_pert.value, dt_err_soft.value, dt_min_arc);
    SystemHard system_hard_connected;
    system_hard_connected.setParam(r_bin.value, r_out.value, r_in, eps.value, dt_limit_hard, dt_min_hermite, eta.value, time_sys, sd_factor.value, file_header.id_offset, n_split.value);
    system_hard_connected.setARCParam(e_err_arc.value, dt_err_pert.value, dt_err_soft.value, dt_min_arc);

    PS::ReallocatableArray<PS::S32> remove_list;
                                       
    SearchCluster search_cluster;
    search_cluster.initialize();


    if(!restart_flag) {
        file_header.n_body = system_soft.getNumberOfParticleGlobal();
        file_header.dt_soft= 0.0;
        std::string fname = fname_snp.value+"."+std::to_string(file_header.nfile);
        file_header.dt_soft= dt_soft.value;
    }

    //if(my_rank==0) {
    //    std::cout<<"----- Initial status -----\n";
    //    stat.eng_init.print(std::cout);
    //    stat.dump(fstatus,WRITE_WIDTH);
    //    fstatus<<std::endl;
    //}
#ifdef MAIN_DEBUG
    FILE* fout;
    if ( (fout = fopen("nbody.dat","w")) == NULL) {
        fprintf(stderr,"Error: Cannot open file nbody.dat\n");
        abort();
    }
#endif


    bool first_step_flag = true;
    bool output_flag = false;
    PS::S64 n_loop = 0;
    PS::F64 dt_kick = dt_soft.value;
/// Main loop
    while(time_sys < time_end.value){

#ifdef PROFILE
        profile.tot.start();

        profile.tree_nb.start();
#endif
        // >1. Tree for neighbor searching ----------------------------------------
#ifndef USE_SIMD
        tree_nb.calcForceAllAndWriteBack(CalcForceEpEpWithLinearCutoffNoSIMD(), system_soft, dinfo);
#else
        tree_nb.calcForceAllAndWriteBack(CalcForceEpEpWithLinearCutoffSimd(), system_soft, dinfo);
#endif

#ifdef PROFILE
        profile.tree_nb.end();
        profile.search_cluster.start();
#endif
        // >2.1 search clusters ----------------------------------------
        search_cluster.searchNeighborOMP<SystemSoft, Tree, EPJSoft>
            (system_soft, tree_nb, r_out.value, r_in, pos_domain, EPISoft::eps*EPISoft::eps);

        search_cluster.searchClusterLocal();
        search_cluster.setIdClusterLocal();

        // >2.2 Send/receive connect cluster
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
        search_cluster.connectNodes(pos_domain,tree_soft);
        search_cluster.setIdClusterGlobalIteration();
        search_cluster.sendAndRecvCluster(system_soft);
#endif

        // record real particle n_loc/glb
        n_loc = system_soft.getNumberOfParticleLocal();
        n_glb = system_soft.getNumberOfParticleGlobal();

        // >2.3 Find ARC groups and create artificial particles
        search_cluster.findGroupsAndCreateArtificalParticles(system_soft);

        // update n_glb, n_loc for all
        PS::S64 n_loc_all = system_soft.getNumberOfParticleLocal();
        PS::S64 n_glb_all = system_soft.getNumberOfParticleGlobal();
            
        // >2.4 set adr/rank for artificial particles in GPS
#pragma omp parallel for
        for(PS::S32 i=n_loc; i<n_loc_new; i++){
            system_soft[i].rank_org = my_rank;
            system_soft[i].adr = i;
            assert(system_soft[i].id<0&&system_soft[i].status<0);
        }

        // >3 Tree for force ----------------------------------------
#ifdef PROFILE
        profile.search_cluster.end();

        profile.soft_tot.start();
        profile.tree_soft.start();

        tree_soft.clearNumberOfInteraction();
        tree_soft.clearTimeProfile();
#endif
        
#ifndef USE_SIMD
        tree_soft.calcForceAllAndWriteBack(CalcForceEpEpWithLinearCutoffNoSIMD(),
#ifdef USE_QUAD
                                           CalcForceEpSpQuadNoSimd(),
#else
                                           CalcForceEpSpNoSIMD(),
#endif
                                           system_soft,
                                           dinfo);
#else
        tree_soft.calcForceAllAndWriteBack(CalcForceEpEpWithLinearCutoffSimd(),
#ifdef USE_QUAD
                                           CalcForceEpSpQuadSimd(),
#else
                                           CalcForceEpSpMonoSimd(),
#endif
                                           system_soft,
                                           dinfo);
#endif

#ifdef PROFILE
        n_count.ep_ep_interact     += tree_soft.getNumberOfInteractionEPEPLocal();
        n_count_sum.ep_ep_interact += tree_soft.getNumberOfInteractionEPEPGlobal();
        n_count.ep_sp_interact     += tree_soft.getNumberOfInteractionEPSPLocal();
        n_count_sum.ep_sp_interact += tree_soft.getNumberOfInteractionEPSPGlobal(); 

        ps_profile += tree_soft.getTimeProfile();
        domain_decompose_weight = ps_profile.calc_force;
        profile.tree_soft.end();
        profile.soft_tot.end();
#endif

        // for first step
        if(first_step_flag) {
            first_step_flag = false;

            // update status
            stat.time = time_sys;
            stat.N = n_glb;
            stat.N_all = n_glb_all;

            // calculate initial energy
            stat.eng_init.clear();
            stat.eng_init.calc(&system_soft[0], system_soft.getNumberOfParticleLocal(), 0.0);
            stat.eng_init.getSumMultiNodes();
            stat.eng_now = stat.eng_init;
             
            // print status
            if(my_rank==0) {
                std::cout<<std::endl;
                stat.print(std::cout);
                stat.dump(fstatus, WRITE_WIDTH);
                fstatus<<std::endl;
            }
            
            output_flag = false;

            // set kick step
            dt_kick = dt_soft.value*0.5;
        }
        else {
            // output step, reduce step by half 
            if( fmod(time_sys, dt_snp.value) == 0.0) {
                output_flag = true;
                dt_kick = dt_soft.value*0.5;
            } 
            else {
                output_flag = false;
                // double kick step
                dt_kick = dt_soft.value;
            }
        }

#ifdef PROFILE
        profile.soft_tot.start();
        profile.kick.start();
#endif
        // >4. kick  ----------------------------------------
        Kick(system_soft, tree_soft, dt_kick);
#ifdef PROFILE
        profile.kick.end();
        profile.soft_tot.end();
#endif

        // output information
        if( fmod(time_sys, dt_snp.value) == 0.0){

#ifdef MAIN_DEBUG
            write_p(fout, time_sys, dt_soft.value*0.5, system_soft, stat.eng_now, stat.eng_diff);
//        //output
//        PS::S32 ntot = system_soft.getNumberOfParticleLocal();
//        fout<<std::setprecision(17)<<time_sys<<" ";
//        for (PS::S32 i=0;i<ntot;i++){
//          fout<<system_soft[i].mass<<" ";
//          for (PS::S32 k=0;k<3;k++) fout<<system_soft[i].pos[k]<<" ";
//          for (PS::S32 k=0;k<3;k++) fout<<system_soft[i].vel[k]<<" ";
//          fout<<system_soft[i].pot_tot<<" ";
//          fout<<0.5*system_soft[i].mass*system_soft[i].vel*system_soft[i].vel<<" ";
//        }
//        fout<<std::endl;
#endif

            // update status
            stat.time = time_sys;
            stat.N = n_glb;
            
            stat.eng_now.clear();
            stat.eng_now.calc(&system_soft[0], system_soft.getNumberOfParticleLocal(), dt_soft.value*0.5);
            stat.eng_now.getSumMultiNodes();
        
            stat.eng_diff = stat.eng_now - stat.eng_init;

            // print status
            if(my_rank==0) {
                std::cout<<std::endl;
                stat.print(std::cout);
                stat.dump(fstatus, WRITE_WIDTH);
                fstatus<<std::endl;
            }
            
#ifdef PROFILE
            //const int NProc=PS::Comm::getNumberOfProc();

            if(my_rank==0) {
                std::cout<<std::setprecision(5);
                std::cout<<"Tree step number: "<<dn_loop<<std::endl;

                std::cout<<"**** Wallclock time per step (local):\n";
                //std::cout<<std::setw(PRINT_WIDTH)<<"Rank";
                profile.dumpName(std::cout,PRINT_WIDTH);
                std::cout<<std::endl;

                //std::cout<<std::setw(PRINT_WIDTH)<<my_rank;
                profile.dump(std::cout,PRINT_WIDTH,dn_loop);
                std::cout<<std::endl;

                std::cout<<"**** FDPS time profile (local):\n";
                ps_profile.dumpName(std::cout,PRINT_WIDTH);
                std::cout<<std::endl;
                ps_profile.dump(std::cout,PRINT_WIDTH,dn_loop);
                std::cout<<std::endl;

                std::cout<<"**** Number per step (global):\n";
                n_count_sum.dumpName(std::cout,PRINT_WIDTH);
                std::cout<<std::endl;
                n_count_sum.dump(std::cout,PRINT_WIDTH,dn_loop);
                std::cout<<std::endl;
                
                std::cout<<"**** Number of members in clusters (local):\n";
                n_count.printHist(std::cout,PRINT_WIDTH,dn_loop);
            }
                    
//#ifdef ARC_ERROR
//            for (int i=1;i<20;i++)  system_hard_isolated.N_count[i] = PS::Comm::getSum(system_hard_isolated.N_count[i]);
//            if(my_rank==0) {
//                std::cerr<<"NHist: ";
//                for (PS::S32 i=0;i<20;i++) {
//                    std::cerr<<system_hard_isolated.N_count[i]<<" ";
//                    system_hard_isolated.N_count[i] = 0;
//                }
//                std::cerr<<std::endl;
//            }
//#endif
            

            fprofile<<std::setprecision(WRITE_PRECISION);
            fprofile<<std::setw(WRITE_WIDTH)<<my_rank;
            fprofile<<std::setw(WRITE_WIDTH)<<time_sys
                    <<std::setw(WRITE_WIDTH)<<dn_loop
                    <<std::setw(WRITE_WIDTH)<<n_loc;
            profile.dump(fprofile, WRITE_WIDTH, dn_loop);
            ps_profile.dump(fprofile, WRITE_WIDTH, dn_loop);
            n_count.dump(fprofile, WRITE_WIDTH, dn_loop);
            fprofile<<std::endl;

#endif

            // data output
            file_header.n_body = n_glb;
            file_header.time = time_sys;
            file_header.nfile++;
            std::string fname = fname_snp.value+"."+std::to_string(file_header.nfile);
            if (data_format.value==1||data_format.value==3)
                system_soft.writeParticleAscii(fname.c_str(), file_header);
            else if(data_format.value==0||data_format.value==2)
                system_soft.writeParticleBinary(fname.c_str(), file_header);


#ifdef PROFILE            
            profile.clear();
            ps_profile.clear();
            n_count.clear();
            n_count_sum.clear();
            dn_loop=0;

            // second half kick
            profile.soft_tot.start();
            profile.kick.start();
#endif
            Kick(system_soft, tree_soft, dt_soft.value);

#ifdef PROFILE
            profile.kick.end();
            profile.soft_tot.end();
#endif
        }

#endif
        profile.search_cluster.start();
#endif

        // Send receive cluster particles after kick
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
        search_cluster.connectNodes(pos_domain,tree_soft);
        search_cluster.setIdClusterGlobalIteration();
        search_cluster.sendAndRecvCluster(system_soft);
#endif

#ifdef PROFILE
        profile.search_cluster.end();

        // >5. Hard integration --------------------------------------
        profile.hard_tot.start();
#endif

        ////// set time
        system_hard_one_cluster.setTimeOrigin(time_sys);
        system_hard_isolated.setTimeOrigin(time_sys);
        system_hard_connected.setTimeOrigin(time_sys);
        ////// set time
        

#ifdef PROFILE
        profile.hard_single.start();
#endif
        ////// integrater one cluster
        system_hard_one_cluster.initializeForOneCluster(search_cluster.getAdrSysOneCluster().size());
        system_hard_one_cluster.setPtclForOneCluster(system_soft, search_cluster.getAdrSysOneCluster());
        system_hard_one_cluster.driveForOneCluster(dt_soft.value);
        system_hard_one_cluster.writeBackPtclForOneClusterOMP(system_soft, search_cluster.getAdrSysOneCluster());
        ////// integrater one cluster
#ifdef PROFILE
        profile.hard_single.end();
#endif
        ////////////////

        /////////////
#ifdef PROFILE
        profile.hard_isolated.start();
#endif
        // integrate multi cluster A
        system_hard_isolated.setPtclForIsolatedMultiCluster(system_soft, search_cluster.adr_sys_multi_cluster_isolated_, search_cluster.n_ptcl_in_multi_cluster_isolated_);
        system_hard_isolated.driveForMultiClusterOMP<SystemSoft, FPSoft>(dt_soft.value,system_soft,first_step_flag);
        system_hard_isolated.writeBackPtclForMultiCluster(system_soft, search_cluster.adr_sys_multi_cluster_isolated_,remove_list);
        // integrate multi cluster A
#ifdef PROFILE
        profile.hard_isolated.end();
#endif
        /////////////

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
        /////////////
#ifdef PROFILE
        PS::Comm::barrier();
        profile.hard_connected.start();
#endif
        // integrate multi cluster B
        system_hard_connected.setPtclForConnectedCluster(system_soft, search_cluster.mediator_sorted_id_cluster_, search_cluster.ptcl_recv_);
        system_hard_connected.driveForMultiClusterOMP<SystemSoft, FPSoft>(dt_soft.value,system_soft,first_step_flag);
        search_cluster.writeAndSendBackPtcl(system_soft, system_hard_connected.getPtcl(), remove_list);
        // integrate multi cluster B
#ifdef PROFILE
        profile.hard_connected.end();
#endif        

#endif
        
//#ifdef MAIN_DEBUG
//        n_loc = system_soft.getNumberOfParticleLocal();
//        for(PS::S32 i=0; i<n_loc; i++){
//            if(system_soft[i].id<0&&system_soft[i].status<0) {
//                std::cerr<<"Error! Ghost detected in system_soft after remove, i="<<i<<std::endl;
//                abort();
//            }
//        }
//#endif


        /////////////
        // Remove ghost particles
        system_soft.removeParticle(remove_list.getPointer(), remove_list.size());
        // reset particle number
        system_soft.setNumberOfParticleLocal(n_loc-remove_list.size());
        remove_list.resizeNoInitialize(0);

#ifdef PROFILE
        profile.hard_tot.end();
        /////////////


        // > 6. Domain decomposition
#ifdef PROFILE
        profile.domain_ex_ptcl.start();
#endif
        // Domain decomposition, parrticle exchange and force calculation
        if(n_loop % 16 == 0) {
            dinfo.decomposeDomainAll(system_soft,domain_decompose_weight);
            //std::cout<<"rank: "<<my_rank<<" weight: "<<domain_decompose_weight<<std::endl;
        }
        system_soft.exchangeParticle(dinfo);

        n_loc = system_soft.getNumberOfParticleLocal();
#pragma omp parallel for
        for(PS::S32 i=0; i<n_loc; i++){
            system_soft[i].rank_org = my_rank;
            system_soft[i].adr = i;
        }
        
#ifdef PROFILE
        profile.domain_ex_ptcl.end();
        profile.tot.end();
#endif

        time_sys += dt_soft.value;

#ifdef PROFILE
        // profile analysis
        PS::S32 n_hard_single     = system_hard_one_cluster.getPtcl().size();
        PS::S32 n_hard_isolated   = system_hard_isolated.getPtcl().size();
        PS::S32 n_hard_connected  = system_hard_connected.getPtcl().size();

        n_count.hard_single      += n_hard_single;
        n_count.hard_isolated    += n_hard_isolated;
        n_count.hard_connected   += n_hard_connected;

        n_count_sum.hard_single      += PS::Comm::getSum(n_hard_single);
        n_count_sum.hard_isolated    += PS::Comm::getSum(n_hard_isolated);
        n_count_sum.hard_connected   += PS::Comm::getSum(n_hard_connected);
                                           
        PS::S64 ARC_substep_sum   = system_hard_isolated.ARC_substep_sum;
        PS::S64 ARC_n_groups      = system_hard_isolated.ARC_n_groups;
        n_count.ARC_substep_sum  += ARC_substep_sum;
        n_count.ARC_n_groups     += ARC_n_groups;

        n_count_sum.ARC_substep_sum  += PS::Comm::getSum(ARC_substep_sum);
        n_count_sum.ARC_n_groups     += PS::Comm::getSum(ARC_n_groups);

        system_hard_isolated.ARC_substep_sum = 0;
        system_hard_isolated.ARC_n_groups = 0;
                                           
        n_count.cluster_count(1, n_hard_single);

        const PS::S32  n_isolated_cluster = system_hard_isolated.getNCluster();
        n_count.cluster_isolated += n_isolated_cluster;
        n_count_sum.cluster_isolated += PS::Comm::getSum(n_isolated_cluster);

        const PS::S32* isolated_cluster_n_list = system_hard_isolated.getClusterNList();
        for (PS::S32 i=0; i<n_isolated_cluster; i++) n_count.cluster_count(isolated_cluster_n_list[i]);

        const PS::S32  n_connected_cluster = system_hard_connected.getNCluster();
        n_count.cluster_connected += n_connected_cluster;
        n_count_sum.cluster_connected += PS::Comm::getSum(n_connected_cluster);
        const PS::S32* connected_cluster_n_list = system_hard_connected.getClusterNList();
        for (PS::S32 i=0; i<n_connected_cluster; i++) n_count.cluster_count(connected_cluster_n_list[i]);

        dn_loop++;
#endif
        n_loop++;
    }

//#ifdef ARC_ERROR
//    std::cout<<"Hist[1]"<<system_hard_isolated.N_count[1]/(PS::F64)n_loop<<std::endl;
//    for (int i=1;i<20;i++)  system_hard_isolated.N_count[i] = PS::Comm::getSum(system_hard_isolated.N_count[i]);
//    if(my_rank==0) {
//        std::cout<<"NHist: ";
//        for (PS::S32 i=0;i<20;i++) std::cout<<system_hard_isolated.N_count[i]/(PS::F64)n_loop<<" ";
//        std::cout<<std::endl;
//    }
//#endif
#ifdef MAIN_DEBUG
    fclose(fout);
#endif

    PS::Finalize();
    return 0;
}

