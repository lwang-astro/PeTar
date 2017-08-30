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
#include<unistd.h>
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
#include"profile.hpp"
#include"domain.hpp"
#include"AR.h" /// include AR.h (L.Wang)
//#include"cluster.hpp"
#include"cluster_list.hpp"

#ifdef USE_QUAD
typedef PS::TreeForForceLong<ForceSoft, EPISoft, EPJSoft>::QuadrupoleWithSymmetrySearch Tree; 
#else
typedef PS::TreeForForceLong<ForceSoft, EPISoft, EPJSoft>::MonopoleWithSymmetrySearch Tree;
//typedef PS::TreeForForceLong<ForceSoft, EPISoft, EPJSoft>::MonopoleWithScatterSearch Tree;
#endif
typedef PS::ParticleSystem<FPSoft> SystemSoft;

#ifdef MAIN_DEBUG
// flag: 1: c.m; 2: individual; 
template<class Teng, class Tsys>
void write_p(FILE* fout, const PS::F64 time, const PS::F64 dt_soft, const Tsys& p, const Teng &et, const Teng &ediff) {
    fprintf(fout,"%20.14e ",time);
    fprintf(fout,"%20.14e %20.14e %20.14e %20.14e %20.14e %20.14e %20.14e %20.14e %20.14e %20.14e %20.14e %20.14e ",
            ediff.tot/et.tot, et.kin, et.pot, et.tot,
            ediff.Lt/et.Lt, ediff.L[0]/et.Lt, ediff.L[1]/et.Lt, ediff.L[2]/et.Lt,
               et.Lt,    et.L[0],    et.L[1],    et.L[2]);
    for (int i=0; i<p.getNumberOfParticleLocal(); i++) {
        if(p[i].status>0||p[i].id<0) continue;
        PS::F64 mi = p[i].mass;
        PS::F64vec vi = p[i].vel;
        if(p[i].status!=0) {
            mi = p[i].mass_bk;
            vi += p[i].acc * dt_soft;
        }
        fprintf(fout,"%20.14e %20.14e %20.14e %20.14e %20.14e %20.14e %20.14e ", 
                mi, p[i].pos[0], p[i].pos[1], p[i].pos[2], 
                vi[0], vi[1], vi[2]);
    }
    fprintf(fout,"\n");
}
#endif

int main(int argc, char *argv[]){
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);
    PS::Initialize(argc, argv);

    PS::S32 my_rank = PS::Comm::getRank();

    Wtime profile;

	PS::S64 n_ptcl_hard_one_cluster = 0;
	PS::S64 n_ptcl_hard_isolated_cluster = 0;
	PS::S64 n_ptcl_hard_nonisolated_cluster = 0;

    // initial parameters
    IOParams<PS::F64> ratio_r_cut  (0.1,  "r_in / r_out");
    IOParams<PS::F64> theta        (0.3,  "openning angle theta");
    IOParams<PS::S32> n_leaf_limit (8,    "tree leaf number limit");
    IOParams<PS::S32> n_group_limit(64,   "tree group number limit");
    IOParams<PS::S32> n_smp_ave    (100,  "average target number of sample particles per process");
    IOParams<PS::S32> n_split      (8,    "number of binary sample points for tree perturbation force");
    IOParams<PS::S64> n_bin        (0,    "number of primordial binaries (assume binaries ID=1,2*n_bin)");
    IOParams<PS::F64> time_end     (10.0, "finishing time");
    IOParams<PS::F64> eta          (0.1,  "Hermite time step coefficient eta");
    IOParams<PS::S64> n_glb        (16384,"Total number of particles");
    IOParams<PS::F64> dt_soft      (0.0,  "Tree timestep","0.1*r_in/sigma");
    IOParams<PS::F64> dt_snp       (0.0625,"Output time interval of particle dataset");
    IOParams<PS::F64> search_factor(1.0,  "neighbor searching coefficient");
    IOParams<PS::F64> dt_limit_hard_factor(4.0, "limit of tree time step/hard time step");
    IOParams<PS::F64> eps          (0.0,  "softerning eps");
    IOParams<PS::F64> r_out        (0.0,  "transit function outer boundary radius", "3.0*<m>/sigma^2");
    IOParams<PS::F64> r_bin        (0.0,  "maximum binary radius criterion", "0.1*r_in");
    IOParams<PS::F64> sd_factor    (1e-8, "Slowdown perturbation criterion");
    IOParams<PS::S32> data_format  (1,    "Data read(r)/write(w) format BINARY(B)/ASCII(A)","r-B/w-A (3); r-A/w-B (2); rw-A (1); rw-B (0)");
    IOParams<std::string> fname_snp("data","Prefix filename of dataset: [prefix].[File ID]");

    // PS::F64 g_min = 1e-6;

    // reading parameters
    int c;
    bool reading_flag=false;

    while((c=getopt(argc,argv,"i:b:B:T:t:e:E:n:N:s:S:d:D:o:l:r:R:X:p:f:h")) != -1){
        switch(c){
        case 'i':
            reading_flag=true;
            data_format.value = atoi(optarg);
            if(my_rank == 0) data_format.print(std::cerr);
            assert(data_format.value>=0||data_format.value<=3);
            break;
        case 'b':
            r_bin.value = atof(optarg);
            if(my_rank == 0) r_bin.print(std::cerr);
            assert(r_bin.value>0.0);
            break;
        case 'B':
            n_bin.value = atoi(optarg);
            if(my_rank == 0) n_bin.print(std::cerr);
            assert(n_bin.value>=0);
            break;
        case 'T':
            theta.value = atof(optarg);
            if(my_rank == 0) theta.print(std::cerr);
            assert(theta.value>=0.0);
            break;
        case 't':
            time_end.value = atof(optarg);
            if(my_rank == 0) time_end.print(std::cerr);
            assert(time_end.value>=0.0);
            break;
        case 'e':
            eps.value = atof(optarg);
            if(my_rank == 0) eps.print(std::cerr);
            assert(eps.value>=0.0);
            break;
        case 'E':
            eta.value = atof(optarg);
            if(my_rank == 0) eta.print(std::cerr);
            assert(eta.value>0.0);
            break;
        case 'n':
            n_group_limit.value = atoi(optarg);
            if(my_rank == 0) n_group_limit.print(std::cerr);
            assert(n_group_limit.value>0);
            break;
        case 'N':
            n_glb.value = atol(optarg);
            if(my_rank == 0) n_glb.print(std::cerr);
            assert(n_glb.value>0);
            break;
        case 's':
            n_smp_ave.value = atoi(optarg);
            if(my_rank == 0) n_smp_ave.print(std::cerr);
            assert(n_smp_ave.value>0.0);
            break;
        case 'S':
            search_factor.value = atof(optarg);
            if(my_rank == 0) search_factor.print(std::cerr);
            assert(search_factor.value>0.0);
            break;
        case 'd':
            sd_factor.value = atof(optarg);
            if(my_rank == 0) sd_factor.print(std::cerr);
            assert(sd_factor.value>0.0);
            break;
        case 'D':
            dt_soft.value = atof(optarg);
            if(my_rank == 0) dt_soft.print(std::cerr);
            assert(dt_soft.value>0.0);
            break;
        case 'o':
            dt_snp.value = atof(optarg);
            if(my_rank == 0) dt_snp.print(std::cerr);
            assert(dt_snp.value>0.0);
            break;
        case 'l':
            n_leaf_limit.value = atoi(optarg);
            if(my_rank == 0) n_leaf_limit.print(std::cerr);
            assert(n_leaf_limit.value>0);
            break;
        case 'r':
            ratio_r_cut.value = atof(optarg);
            if(my_rank == 0) ratio_r_cut.print(std::cerr);
            assert(ratio_r_cut.value>0.0);
            assert(ratio_r_cut.value<1.0);
            break;
        case 'R':
            r_out.value = atof(optarg);
            if(my_rank == 0) r_out.print(std::cerr);
            assert(r_out.value>0.0);
            break;
        case 'X':
            dt_limit_hard_factor.value = atof(optarg);
            if(my_rank == 0) dt_limit_hard_factor.print(std::cerr);
            assert(dt_limit_hard_factor.value > 0.0);
            break;
        case 'p':
            n_split.value = atoi(optarg);
            if(my_rank == 0) n_split.print(std::cerr);
            assert(n_split.value>=8);
            break;
        case 'f':
            fname_snp.value = optarg;
            if(my_rank == 0) fname_snp.print(std::cerr);
            break;
        case 'h':
            std::cerr<<"Usage: nbody.out [option] [filename]"<<std::endl;
            std::cerr<<"       Option defaulted values are shown after ':'\n"<<std::endl;
            std::cerr<<"  -i: [I] enable reading data file: disabled with Plummer model"<<std::endl;
            std::cerr<<"          "<<data_format<<std::endl;
            std::cerr<<"          File content:\n"
                     <<"            First line: \n"
                     <<"             1. File_ID: 0 for initialization, else for restarting\n"
                     <<"             2. N_particle \n"
                     <<"             3. N_offset: for naming special particle ID, should be > N_particle\n"
                     <<"             4. Time\n"
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
            std::cerr<<"  -b: [F] "<<r_bin<<std::endl;
            std::cerr<<"  -B: [I] "<<n_bin<<std::endl;
            std::cerr<<"  -T: [F] "<<theta<<std::endl;
            std::cerr<<"  -t: [F] "<<time_end<<std::endl;
            std::cerr<<"  -e: [F] "<<eps<<std::endl;
            std::cerr<<"  -E: [F] "<<eta<<std::endl;
            std::cerr<<"  -n: [I] "<<n_group_limit<<std::endl;
            std::cerr<<"  -N: [I] "<<n_glb<<std::endl;
            std::cerr<<"  -s: [I] "<<n_smp_ave<<std::endl;
            std::cerr<<"  -S: [F] "<<search_factor<<std::endl;
            std::cerr<<"  -d: [F] "<<sd_factor<<std::endl;
            std::cerr<<"  -D: [F] "<<dt_soft<<std::endl;
            std::cerr<<"  -o: [F] "<<dt_snp<<std::endl;
            std::cerr<<"  -l: [I] "<<n_leaf_limit<<std::endl;
            std::cerr<<"  -r: [F] "<<ratio_r_cut<<std::endl;
            std::cerr<<"  -R: [F] "<<r_out<<std::endl;
            std::cerr<<"  -X: [F] "<<dt_limit_hard_factor<<std::endl;
            std::cerr<<"  -p: [I] "<<n_split<<std::endl;
            std::cerr<<"  -f: [S] "<<fname_snp<<std::endl;
            std::cerr<<"*** PS: r_in : transit function inner boundary radius\n"
                     <<"        r_out: transit function outer boundary radius\n"
                     <<"        sigma: half-mass radius velocity dispersion\n"
                     <<"        n_bin: number of primordial binaries\n"
                     <<"        <m>  : averaged mass"<<std::endl;
            PS::Finalize();
            return 0;
        }
    }

    PS::F64 time_sys = 0.0;
    PS::S32 n_loc;

    SystemSoft system_soft;
    system_soft.initialize();
    system_soft.setAverageTargetNumberOfSampleParticlePerProcess(n_smp_ave.value);
    
    FileHeader file_header;
    if (reading_flag) {
      char* sinput=argv[argc-1];
      if(data_format.value==1||data_format.value==2)
          system_soft.readParticleAscii(sinput, file_header);
      else
          system_soft.readParticleBinary(sinput, file_header);
      time_sys = file_header.time;
      PS::Comm::broadcast(&time_sys, 1, 0);
      n_glb.value = system_soft.getNumberOfParticleGlobal();
      n_loc = system_soft.getNumberOfParticleLocal();
      //      for(PS::S32 i=0; i<n_loc; i++) system_soft[i].id = i;
      std::cerr<<"Reading file "<<sinput<<std::endl
               <<"N_tot = "<<n_glb.value<<"\nN_loc = "<<n_loc<<std::endl;
    }
    else SetParticlePlummer(system_soft, n_glb.value, n_loc, time_sys);

    bool restart_flag = file_header.nfile; // nfile = 0 is assumed as initial data file

    PS::Comm::barrier();

    // tree time step and n_split
    if (restart_flag) {
        if(dt_soft.value!=file_header.dt_soft) 
            std::cerr<<"Warning: tree time step cannot be changed for restarting, the value from the data file ("<<file_header.dt_soft<<") will be used\n";
        if(n_split.value!=file_header.n_split)
            std::cerr<<"Warning: n_split cannot be changed for restarting, the value from the data file ("<<file_header.n_split<<") will be used\n";
        dt_soft.value = file_header.dt_soft;
        n_split.value = file_header.n_split;
    }
    else {
        file_header.dt_soft = dt_soft.value;
        file_header.n_split = n_split.value;        
    }
    
    PS::F64 r_in, m_average, v_disp, r_search_min;
    GetR(system_soft, r_in, r_out.value, r_bin.value, r_search_min, m_average, dt_soft.value, v_disp, search_factor.value, ratio_r_cut.value, n_bin.value);
//    EPISoft::r_out = r_out;
//    EPISoft::r_in  = r_in;
    EPISoft::eps   = eps.value;
    EPISoft::r_out = EPJSoft::r_out = FPSoft::r_out = r_out.value;
    Ptcl::search_factor = search_factor.value;
    Ptcl::r_search_min = r_search_min;
//    EPJSoft::r_search_min = r_out*search_factor;
//    EPJSoft::m_average = m_average;
  
    // ID safety check
    if (!restart_flag) {
        for (PS::S32 i=0; i<n_loc; i++) {
            if(system_soft[i].id<=0) {
                std::cerr<<"Error: for initial data, the id should always larger than zero. current index i = "<<i<<", id = "<<system_soft[i].id<<"!"<<std::endl;
                abort();
            }
            
            system_soft[i].calcRSearch(dt_soft.value);
        }
    }
    
    if(my_rank == 0) {
        std::cerr<<"Parameter list:\n";
        std::cerr<<" m_average    = "<<m_average      <<std::endl
                 <<" r_in         = "<<r_in           <<std::endl
                 <<" r_out        = "<<r_out.value    <<std::endl
                 <<" r_search_min = "<<r_search_min   <<std::endl
                 <<" vel_disp     = "<<v_disp         <<std::endl
                 <<" dt_soft      = "<<dt_soft.value  <<std::endl;
    }

    // set r_search
    //if(n_bin>n_loc) {
    //    if (n_loc%2) std::cerr<<"Warning! Binary number is larger than local particle numbers, but n_loc is odd, which means the last binary is splitted to two nodes"<<std::endl;
    //    int ndiv = 2*n_bin/n_loc;
    //    if(my_rank<ndiv) n_bin = n_loc/2;
    //    else if(my_rank==ndiv) n_bin = 2*n_bin % n_loc;
    //    else n_bin = 0;
    //}
    //else {
    //    if(my_rank>0) n_bin = 0;
    //}

    //if(n_bin>0) SetBinaryRout(system_soft, n_bin, g_min, r_in, r_out, m_average);
    //SetSingleRout(system_soft, n_loc, 2*n_bin, r_out);

    const PS::F32 coef_ema = 0.2;
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);
    dinfo.decomposeDomainAll(system_soft);
    const PS::S32 n_proc = PS::Comm::getNumberOfProc();
    PS::F64ort * pos_domain = new PS::F64ort[n_proc];
    for(PS::S32 i=0; i<n_proc; i++) pos_domain[i] = dinfo.getPosDomain(i);

    system_soft.exchangeParticle(dinfo); 

    n_loc = system_soft.getNumberOfParticleLocal();

#pragma omp parallel for
    for(PS::S32 i=0; i<n_loc; i++){
        system_soft[i].rank_org = my_rank;
        system_soft[i].adr = i;
    }

    Tree tree_soft;
    tree_soft.initialize(n_glb.value, theta.value, n_leaf_limit.value, n_group_limit.value);
    tree_soft.calcForceAllAndWriteBack(CalcForceEpEpWithLinearCutoffNoSIMD(),
#ifdef USE_QUAD
                                       CalcForceEpSpQuadNoSimd(),
#else
                                       CalcForceEpSpNoSIMD(),
#endif
                                       system_soft,
                                       dinfo);

    SystemHard system_hard_one_cluster;
    PS::F64 dt_limit_hard = dt_soft.value/dt_limit_hard_factor.value;
    system_hard_one_cluster.setParam(r_bin.value, r_out.value, r_in, eps.value, dt_limit_hard, eta.value, time_sys, sd_factor.value, file_header.id_offset, n_split.value);
    // system_hard_one_cluster.setARCParam();
    SystemHard system_hard_isolated;
    system_hard_isolated.setParam(r_bin.value, r_out.value, r_in, eps.value,  dt_limit_hard, eta.value, time_sys, sd_factor.value, file_header.id_offset, n_split.value);
    system_hard_isolated.setARCParam();
    SystemHard system_hard_conected;
    system_hard_conected.setParam(r_bin.value, r_out.value, r_in, eps.value, dt_limit_hard, eta.value, time_sys, sd_factor.value, file_header.id_offset, n_split.value);
    system_hard_conected.setARCParam();

    SearchCluster search_cluster;
    search_cluster.initialize();
    search_cluster.searchNeighborAndCalcHardForceOMP<SystemSoft, Tree, EPJSoft>
      (system_soft, tree_soft, r_out.value, r_in, pos_domain, EPISoft::eps*EPISoft::eps);

    if(!restart_flag) {
        search_cluster.searchClusterLocal();
        search_cluster.setIdClusterLocal();
        search_cluster.conectNodes(pos_domain,tree_soft);
        search_cluster.setIdClusterGlobalIteration();
        search_cluster.sendAndRecvCluster(system_soft);

        system_hard_isolated.setPtclForIsolatedMultiCluster(system_soft, search_cluster.adr_sys_multi_cluster_isolated_, search_cluster.n_ptcl_in_multi_cluster_isolated_);
        system_hard_isolated.initialMultiClusterOMP<SystemSoft,FPSoft>(system_soft, dt_soft.value);
        system_hard_isolated.writeBackPtclForMultiCluster(system_soft, search_cluster.adr_sys_multi_cluster_isolated_);
    
        system_hard_conected.setPtclForConectedCluster(system_soft, search_cluster.mediator_sorted_id_cluster_, search_cluster.ptcl_recv_);
        system_hard_conected.initialMultiClusterOMP<SystemSoft,FPSoft>(system_soft, dt_soft.value);
        search_cluster.writeAndSendBackPtcl(system_soft, system_hard_conected.getPtcl());

#pragma omp parallel for
        for(PS::S32 i=0; i<n_loc; i++){
            system_soft[i].rank_org = my_rank;
            system_soft[i].adr = i;
        }
        tree_soft.calcForceAllAndWriteBack(CalcForceEpEpWithLinearCutoffNoSIMD(),
#ifdef USE_QUAD
                                           CalcForceEpSpQuadNoSimd(),
#else
                                           CalcForceEpSpNoSIMD(),
#endif
                                           system_soft,
                                           dinfo);

        search_cluster.searchNeighborAndCalcHardForceOMP<SystemSoft, Tree, EPJSoft>
            (system_soft, tree_soft, r_out.value, r_in, pos_domain, EPISoft::eps*EPISoft::eps);

        file_header.n_body = system_soft.getNumberOfParticleGlobal();
        file_header.dt_soft= 0.0;
        std::string fname = fname_snp.value+"."+std::to_string(file_header.nfile);
        if (data_format.value==1||data_format.value==3)
            system_soft.writeParticleAscii(fname.c_str(), file_header);
        else
            system_soft.writeParticleBinary(fname.c_str(), file_header);
        file_header.dt_soft= dt_soft.value;
    }    

    EnergyAndMomemtum eng_init, eng_now, eng_diff;

    eng_init.clear();
    eng_init.calc(&system_soft[0], system_soft.getNumberOfParticleLocal(), 0.0);
    eng_init.getSumMultiNodes();

    if(my_rank==0) eng_init.dump(std::cerr);
    eng_now = eng_init;


#ifdef MAIN_DEBUG
    FILE* fout;
    if ( (fout = fopen("nbody.dat","w")) == NULL) {
        fprintf(stderr,"Error: Cannot open file nbody.dat\n");
        abort();
    }
    write_p(fout, time_sys, 0.0, system_soft, eng_now, eng_diff);
#endif
    std::ofstream fprofile;
    if(my_rank==0) fprofile.open("profile.out");

    PS::S64 n_loop = 0;
    PS::S64 dn_loop = 0;
    bool first_step_flag = true;
    while(time_sys < time_end.value){
      
        profile.tot.start();
        ////////////////
        ////// 1st kick
        Kick(system_soft, tree_soft, dt_soft.value*0.5);
        ////// 1st kick
        ////////////////
        
        ////////////////
        profile.hard_tot.start();
        ////// set time
        system_hard_one_cluster.setTimeOrigin(time_sys);
        system_hard_isolated.setTimeOrigin(time_sys);
        system_hard_conected.setTimeOrigin(time_sys);
        ////// set time
        ////////////////
        
        ////////////////
        ////// search cluster
        search_cluster.searchClusterLocal();
        search_cluster.setIdClusterLocal();
        search_cluster.conectNodes(pos_domain,tree_soft);
        search_cluster.setIdClusterGlobalIteration();
        search_cluster.sendAndRecvCluster(system_soft);
        ////// search cluster
        ////////////////

        ////////////////
        profile.hard_single.start();
        ////// integrater one cluster
        system_hard_one_cluster.initializeForOneCluster(search_cluster.getAdrSysOneCluster().size());
        system_hard_one_cluster.setPtclForOneCluster(system_soft, search_cluster.getAdrSysOneCluster());
        system_hard_one_cluster.driveForOneCluster(dt_soft.value);
        system_hard_one_cluster.writeBackPtclForOneCluster(system_soft, search_cluster.getAdrSysOneCluster());
        ////// integrater one cluster
        profile.hard_single.end();
        ////////////////

        /////////////
        profile.hard_isolated.start();
        // integrate multi cluster A
        system_hard_isolated.setPtclForIsolatedMultiCluster(system_soft, search_cluster.adr_sys_multi_cluster_isolated_, search_cluster.n_ptcl_in_multi_cluster_isolated_);
        system_hard_isolated.driveForMultiClusterOMP<SystemSoft, FPSoft>(dt_soft.value,system_soft,first_step_flag);
        system_hard_isolated.writeBackPtclForMultiCluster(system_soft, search_cluster.adr_sys_multi_cluster_isolated_);
        // integrate multi cluster A
        profile.hard_isolated.end();
        /////////////

        /////////////
        profile.hard_connected.start();
        // integrate multi cluster B
        system_hard_conected.setPtclForConectedCluster(system_soft, search_cluster.mediator_sorted_id_cluster_, search_cluster.ptcl_recv_);
        system_hard_conected.driveForMultiClusterOMP<SystemSoft, FPSoft>(dt_soft.value,system_soft,first_step_flag);
        search_cluster.writeAndSendBackPtcl(system_soft, system_hard_conected.getPtcl());
        // integrate multi cluster B
        profile.hard_connected.end();

        first_step_flag = false;
        profile.hard_tot.end();
        /////////////

        /////////////
        profile.soft_tot.start();
        // Domain decomposition, parrticle exchange and force calculation

        if(n_loop % 16 == 0) dinfo.decomposeDomainAll(system_soft);
        system_soft.exchangeParticle(dinfo);
        n_loc = system_soft.getNumberOfParticleLocal();

#pragma omp parallel for
        for(PS::S32 i=0; i<n_loc; i++){
            system_soft[i].rank_org = my_rank;
            system_soft[i].adr = i;
        }
        tree_soft.calcForceAllAndWriteBack(CalcForceEpEpWithLinearCutoffNoSIMD(),
#ifdef USE_QUAD
                                           CalcForceEpSpQuadNoSimd(),
#else
                                           CalcForceEpSpNoSIMD(),
#endif
                                           system_soft,
                                           dinfo);
        profile.search_cluster.start();
        search_cluster.searchNeighborAndCalcHardForceOMP<SystemSoft, Tree, EPJSoft>
          (system_soft, tree_soft, r_out.value, r_in, pos_domain, EPISoft::eps*EPISoft::eps);
        profile.search_cluster.end();
        
        profile.soft_tot.end();

        // Domain decomposition, parrticle exchange and force calculation
        /////////////

        ////////////////
        ////// 2nd kick
        Kick(system_soft, tree_soft, dt_soft.value*0.5);
        time_sys += dt_soft.value;
        ////// 2nd kick
        ////////////////

        eng_now.clear();
        eng_now.calc(&system_soft[0], system_soft.getNumberOfParticleLocal(), dt_soft.value*0.5);
        eng_now.getSumMultiNodes();
        
        profile.tot.end();
        /////////////

        eng_diff = eng_now - eng_init;
        PS::S64 n_glb = system_soft.getNumberOfParticleGlobal();
        PS::S32 n_one_cluster         = system_hard_one_cluster.getPtcl().size();
        PS::S32 n_isolated_cluster    = system_hard_isolated.getPtcl().size();
        PS::S32 n_nonisolated_cluster = system_hard_conected.getPtcl().size();
        n_ptcl_hard_one_cluster         += PS::Comm::getSum(n_one_cluster);
        n_ptcl_hard_isolated_cluster    += PS::Comm::getSum(n_isolated_cluster);
        n_ptcl_hard_nonisolated_cluster += PS::Comm::getSum(n_nonisolated_cluster);
        
        dn_loop++;

//#ifdef ARC_ERROR
//        system_hard_isolated.N_count[0] += PS::Comm::getSum(n_one_cluster);
//#endif
        
        if( fmod(time_sys, dt_snp.value) == 0.0){
#ifdef MAIN_DEBUG
            write_p(fout, time_sys, dt_soft.value*0.5, system_soft, eng_now, eng_diff);
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

            //eng_diff.dump(std::cerr);
            if(my_rank==0) {
                std::cerr<<"n_loop= "<<n_loop<<std::endl;
                std::cerr<<"n_glb= "<<n_glb<<std::endl;
                std::cerr<<"Time= "<<time_sys<<" Enow-Einit="<<eng_diff.tot<<" (Enow-Einit)/Einit= "<<eng_diff.tot/eng_init.tot
//#ifdef ARC_ERROR
//                         <<" ARC_error_relative="<<system_hard_isolated.ARC_error_relative+system_hard_conected.ARC_error_relative<<" ARC_error="<<system_hard_isolated.ARC_error+system_hard_conected.ARC_error
//#endif
                         <<std::endl;
                eng_now.dump(std::cerr);
                profile.print(std::cerr,time_sys,dn_loop);
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
            
            if(my_rank==0) {
                fprofile<<std::setprecision(PRINT_PRECISION);
                fprofile<<std::setw(PRINT_WIDTH)<<time_sys
                        <<std::setw(PRINT_WIDTH)<<n_loop
                        <<std::setw(PRINT_WIDTH)<<n_glb;
                profile.dump(fprofile, PRINT_WIDTH, dn_loop);
                fprofile<<std::setw(PRINT_WIDTH)<<(PS::F64)n_ptcl_hard_one_cluster/dn_loop
                        <<std::setw(PRINT_WIDTH)<<(PS::F64)n_ptcl_hard_isolated_cluster/dn_loop
                        <<std::setw(PRINT_WIDTH)<<(PS::F64)n_ptcl_hard_nonisolated_cluster/dn_loop
                        <<std::setw(PRINT_WIDTH)<<eng_diff.tot/eng_init.tot
//#ifdef ARC_ERROR
//                        <<std::setw(PRINT_WIDTH)<<system_hard_isolated.ARC_error+system_hard_conected.ARC_error
//#endif              
                        <<std::endl;
            }
            
            file_header.n_body = system_soft.getNumberOfParticleGlobal();
            file_header.time = time_sys;
            file_header.nfile++;
            std::string fname = fname_snp.value+"."+std::to_string(file_header.nfile);
            if (data_format.value==1||data_format.value==3)
                system_soft.writeParticleAscii(fname.c_str(), file_header);
            else
                system_soft.writeParticleBinary(fname.c_str(), file_header);

//            if (n_bin>0) {
//              char bout[99] = "bin.";
//              sprintf(&bout[4],"%lld",file_header.nfile);
//              FILE* bfout=fopen(bout,"w");
//              for (PS::S64 k=0; k<n_bin; k++) {
//                PS::F64 ax,ecc,inc,OMG,omg,tperi,anomaly;
//                anomaly=PosVel2OrbParam(ax,ecc,inc,OMG,omg,tperi,
//                                        system_soft[2*k].pos, system_soft[2*k+1].pos,
//                                        system_soft[2*k].vel, system_soft[2*k+1].vel,
//                                        system_soft[2*k].mass,system_soft[2*k+1].mass);
//                FPSoft pcm;
//                calc_center_of_mass(pcm, &(system_soft[2*k]), 2);
//                
//                fprintf(bfout,"%lld %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e\n",
//                        k, ax, ecc, inc, OMG, omg, tperi, anomaly, system_soft[2*k].mass, system_soft[2*k+1].mass,
//                        pcm.pos[0],pcm.pos[1],pcm.pos[2],pcm.vel[0],pcm.vel[1],pcm.vel[2]);
//              }
//              fclose(bfout);
//            }
            
            profile.clear();
            n_ptcl_hard_one_cluster = n_ptcl_hard_isolated_cluster = n_ptcl_hard_nonisolated_cluster = 0;

            dn_loop=0;
        }


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

