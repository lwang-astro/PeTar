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

#ifdef ENERGY_CHECK
//#define ENERGY_DIRECT
#endif

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
#include"profile.hpp"
#include"domain.hpp"
#include"AR.h" /// include AR.h (L.Wang)
//#include"cluster.hpp"
#include"cluster_list.hpp"

template<class Tpsys, class Ttree>
void Kick(Tpsys & system,
          const Ttree & tree,
          const PS::F64 dt){
    const PS::S32 n = system.getNumberOfParticleLocal();
#pragma omp parallel for
    for(int i=0; i<n; i++){
	system[i].vel  += system[i].acc * dt;
    }
}

template<class Tpsys, class Ttree>
void Drift(Tpsys & system,
           const Ttree & tree,
           const PS::F64 dt){
    const PS::S32 n = system.getNumberOfParticleLocal();
#pragma omp parallel for
    for(int i=0; i<n; i++){
        //if(tree.getForce(i).n_ngb <= 0){
	if(system[i].n_ngb <= 0){
            system[i].pos  += system[i].vel * dt;
        }
    }
}

#ifdef USE_QUAD
typedef PS::TreeForForceLong<ForceSoft, EPISoft, EPJSoft>::QuadrupoleWithScatterSearch Tree; 
#else
typedef PS::TreeForForceLong<ForceSoft, EPISoft, EPJSoft>::MonopoleWithScatterSearch Tree;
#endif
typedef PS::ParticleSystem<FPSoft> SystemSoft;


int main(int argc, char *argv[]){
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);
    PS::Initialize(argc, argv);

	PS::F64 wtime_tot = 0.0;
	PS::F64 wtime_tot_offset = 0.0;
	PS::F64 wtime_hard_tot = 0.0;
	PS::F64 wtime_hard_tot_offset = 0.0;
	PS::F64 wtime_hard_1st_block = 0.0;
	PS::F64 wtime_hard_1st_block_offset = 0.0;
	PS::F64 wtime_hard_2nd_block = 0.0;
	PS::F64 wtime_hard_2nd_block_offset = 0.0;
	PS::F64 wtime_hard_3rd_block = 0.0;
	PS::F64 wtime_hard_3rd_block_offset = 0.0;
	PS::F64 wtime_soft_tot = 0.0;
	PS::F64 wtime_soft_tot_offset = 0.0;
    //	PS::F64 wtime_soft_force = 0.0;
	PS::F64 wtime_soft_search_neighbor_offset = 0.0;
	PS::F64 wtime_soft_search_neighbor = 0.0;

	PS::S64 n_ptcl_hard_one_cluster = 0;
	PS::S64 n_ptcl_hard_isolated_cluster = 0;
	PS::S64 n_ptcl_hard_nonisolated_cluster = 0;

//#ifdef CALC_HARD_ENERGY
//    PS::F64 dEerr_1body_loc = 0.0;
//    PS::F64 dEerr_1body_glb = 0.0;
//    PS::F64 dEerr_2body_loc = 0.0;
//    PS::F64 dEerr_2body_glb = 0.0;
//    PS::F64 dEerr_mbody_loc = 0.0;
//    PS::F64 dEerr_mbody_glb = 0.0;
//#endif
    PS::F64 ratio_r_cut = 0.5;
    PS::F64 time_sys = 0.0;
    PS::F64 theta = 0.4;
    PS::S32 n_leaf_limit = 8;
    PS::S32 n_group_limit = 64;
    PS::S32 n_smp_ave = 100;
    PS::F64 time_end = 1000.0;
    PS::F64 dt_soft = 1.0/256.0;
    //    PS::F64 eta = 0.1;
    //    PS::F64 eta_s = eta * 0.1;
    //    char dir_name[1024];
    PS::S64 n_glb = 16384;
    PS::S64 n_bin = 0;
    PS::F64 dt_snp = 1.0 / 16.0;
    PS::F64 search_factor = 0.1; 
    PS::F64 dt_limit_hard_factor = 4.0;
    PS::F64 eps = 1e-8;
    PS::F64 r_out = 0.0;
    PS::F64 g_min = 1e-6;
    int c;
    bool reading_flag=false;

    while((c=getopt(argc,argv,"id:t:T:e:n:N:b:s:S:g:l:r:R:X:h")) != -1){
        switch(c){
        case 'i':
            reading_flag=true;
            break;
//        case 'o':
//            sprintf(dir_name, optarg);
//            std::cerr<<"dir_name="<<dir_name<<std::endl;
//            break;
//        case 'd':
//            dt_soft = 1.0 / atof(optarg);
//            std::cerr<<"tree time step="<<dt_soft<<std::endl;
//            break;
        case 't':
            theta = atof(optarg);
            std::cerr<<"tree openning angle theta="<<theta<<std::endl;
            break;
        case 'T':
            time_end = atof(optarg);
            std::cerr<<"finishing time="<<time_end<<std::endl;
            break;
        case 'e':
            eps = atof(optarg);
            std::cerr<<"softening="<<eps<<std::endl;
            break;
        case 'n':
            n_group_limit = atoi(optarg);
            std::cerr<<"n_group_limit="<<n_group_limit<<std::endl;
            break;
        case 'N':
            n_glb = atol(optarg);
            std::cerr<<"Total number of particles="<<n_glb<<std::endl;
            break;
        case 'b':
            n_bin = atol(optarg);
            std::cerr<<"Binary number="<<n_bin<<std::endl;
            break;
        case 's':
            n_smp_ave = atoi(optarg);
            std::cerr<<"n_smp_ave="<<n_smp_ave<<std::endl;
            break;
        case 'S':
            search_factor = atof(optarg);
            std::cerr<<"neighbor searching factor="<<search_factor<<std::endl;
            break;
        case 'g':
            g_min = atof(optarg);
            std::cerr<<"Perturber parameter gamma minimum="<<g_min<<std::endl;
            break;
        case 'l':
            n_leaf_limit = atoi(optarg);
            std::cerr<<"n_leaf_limit="<<n_leaf_limit<<std::endl;
            break;
        case 'r':
            ratio_r_cut = atof(optarg);
            std::cerr<<"r_in/r_out="<<ratio_r_cut<<std::endl;
            break;
        case 'R':
            r_out = atof(optarg);
            std::cerr<<"r_out="<<r_out<<std::endl;
            break;
        case 'X':
            dt_limit_hard_factor = atof(optarg);
            std::cerr<<"soft (tree) time step/hard time step="<<dt_limit_hard_factor<<std::endl;
            assert(dt_limit_hard_factor > 0.0);
            break;
        case 'h':
            std::cerr<<"  -i:     enable reading data file (default: disabled with Plummer model)"<<std::endl;
            //              std::cerr<<"o: dir name of output"<<std::endl;
            //            std::cerr<<"  -d: [F] tree time step (default: "<<dt_soft<<")"<<std::endl;
            std::cerr<<"  -t: [F] openning angle theta (default: "<<theta<<")"<<std::endl;
            std::cerr<<"  -T: [F] finishing time (default: "<<time_end<<")"<<std::endl;
            std::cerr<<"  -e: [F] softening parameter (default: "<<eps<<")"<<std::endl;
            std::cerr<<"  -n: [I] n_group_limit (default: "<<n_group_limit<<")"<<std::endl;
            std::cerr<<"  -N: [I] total number of particles if no -i used (default: "<<n_glb<<")"<<std::endl;
            std::cerr<<"  -b: [I] binary number (default: "<<n_bin<<")"<<std::endl;
            std::cerr<<"  -s: [I] n_smp_ave (default: "<<n_smp_ave<<")"<<std::endl;
            std::cerr<<"  -S: [F] neighbor searching factor (default: "<<search_factor<<")"<<std::endl;
            std::cerr<<"  -g: [F} perturber parameter gamma minimum (default: "<<g_min<<")"<<std::endl;
            std::cerr<<"  -l: [I] n_leaf_limit (default: "<<n_leaf_limit<<")"<<std::endl;
            std::cerr<<"  -r: [F] r_in / r_out (default: "<<ratio_r_cut<<")"<<std::endl;
            std::cerr<<"  -R: [F] r_out (default: autodetermine)"<<std::endl;
            std::cerr<<"  -X: [F] soft (tree) time step/hard time step (default: "<<dt_limit_hard_factor<<")"<<std::endl;
            PS::Finalize();
            return 0;
        }
    }

    SystemSoft system_soft;
    system_soft.initialize();
    system_soft.setAverageTargetNumberOfSampleParticlePerProcess(n_smp_ave);
    PS::S32 n_loc;
    
    FileHeader file_header;
    if (reading_flag) {
      char* sinput=argv[argc-1];
      system_soft.readParticleAscii(sinput, file_header);
      time_sys = file_header.time;
      PS::Comm::broadcast(&time_sys, 1, 0);
      n_glb = system_soft.getNumberOfParticleGlobal();
      n_loc = system_soft.getNumberOfParticleLocal();
      //      for(PS::S32 i=0; i<n_loc; i++) system_soft[i].id = i;
    }
    else SetParticlePlummer(system_soft, n_glb, n_loc, time_sys);

    PS::Comm::barrier();

    PS::F64 r_in, m_average, v_disp;
    GetR(r_in, r_out, m_average, dt_soft, v_disp, system_soft, ratio_r_cut);
    EPISoft::r_out = r_out;
    EPISoft::r_in  = r_in;
    EPISoft::eps   = eps;
    EPJSoft::r_search_min = r_out*search_factor;
    EPJSoft::m_average = m_average;
    
    std::cerr<<"m_average= "<<EPJSoft::m_average<<" r_out="<<EPISoft::r_out<<" r_in="<<EPISoft::r_in<<" dt_soft="<<dt_soft<<std::endl;

    PS::S32 my_rank = PS::Comm::getRank();

    // set r_search
    if (my_rank==0) {
      if(n_bin>0) SetBinaryRsearch(system_soft, n_bin, g_min, r_in);
      SetSingleRsearch(system_soft, n_glb, n_bin, r_in);
    }

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
    tree_soft.initialize(n_glb, theta, n_leaf_limit, n_group_limit);
    tree_soft.calcForceAllAndWriteBack(CalcForceEpEpWithLinearCutoffNoSIMD(),
                                       CalcForceEpSpNoSIMD(),
                                       system_soft,
                                       dinfo);

    SystemHard system_hard_one_cluster;
    system_hard_one_cluster.setParam(EPISoft::r_out, EPISoft::r_in, eps, dt_soft/dt_limit_hard_factor, time_sys, g_min);
    // system_hard_one_cluster.setARCParam();
    SystemHard system_hard_isolated;
    system_hard_isolated.setParam(EPISoft::r_out, EPISoft::r_in, eps,  dt_soft/dt_limit_hard_factor, time_sys, g_min);
    system_hard_isolated.setARCParam();
    SystemHard system_hard_conected;
    system_hard_conected.setParam(EPISoft::r_out, EPISoft::r_in, eps, dt_soft/dt_limit_hard_factor, time_sys, g_min);
    system_hard_conected.setARCParam();

    SearchCluster search_cluster;
    search_cluster.initialize();
    search_cluster.searchNeighborAndCalcHardForceOMP<SystemSoft, Tree, EPJSoft>
      (system_soft, tree_soft, EPISoft::r_out, EPISoft::r_in, pos_domain, EPISoft::eps*EPISoft::eps);

    Energy eng_init, eng_now;

    eng_init.calc(system_soft, true);

    eng_init.dump(std::cerr);
    eng_now = eng_init;

#ifdef DEBUG_OUTPUT
    std::ofstream fout;
	fout.open("pdata");
#endif
    std::ofstream fprofile;
    fprofile.open("profile.out");

    PS::S32 n_loop = 0;
    PS::S32 dn_loop = 0;
    while(time_sys < time_end){
      
        wtime_tot_offset = PS::GetWtime();
        ////////////////
        ////// 1st kick
        Kick(system_soft, tree_soft, dt_soft*0.5);
        ////// 1st kick
        ////////////////
        
        ////////////////
        wtime_hard_tot_offset = PS::GetWtime();
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
        search_cluster.conectNodes(pos_domain);
        search_cluster.setIdClusterGlobalIteration();
        search_cluster.sendAndRecvCluster(system_soft);
        ////// search cluster
        ////////////////

        ////////////////
        wtime_hard_1st_block_offset = PS::GetWtime();
        ////// integrater one cluster
        system_hard_one_cluster.initializeForOneCluster(search_cluster.getAdrSysOneCluster().size());
        system_hard_one_cluster.setPtclForOneCluster(system_soft, search_cluster.getAdrSysOneCluster());
        system_hard_one_cluster.driveForOneCluster(dt_soft);
        system_hard_one_cluster.writeBackPtclForOneCluster(system_soft, search_cluster.getAdrSysOneCluster());
        ////// integrater one cluster
        wtime_hard_1st_block += PS::GetWtime() - wtime_hard_1st_block_offset;
        ////////////////


        /////////////
        wtime_hard_2nd_block_offset = PS::GetWtime();
        // integrate multi cluster A
        system_hard_isolated.setPtclForIsolatedMultiCluster(system_soft,
                                                            search_cluster.adr_sys_multi_cluster_isolated_,
                                                            search_cluster.n_ptcl_in_multi_cluster_isolated_);
        system_hard_isolated.driveForMultiClusterOMP(dt_soft);
        system_hard_isolated.writeBackPtclForMultiCluster(system_soft, search_cluster.adr_sys_multi_cluster_isolated_);
        // integrate multi cluster A
        wtime_hard_2nd_block += PS::GetWtime() - wtime_hard_2nd_block_offset;
        /////////////

        /////////////
        wtime_hard_3rd_block_offset = PS::GetWtime();
        // integrate multi cluster B
        system_hard_conected.setPtclForConectedCluster(system_soft, search_cluster.mediator_sorted_id_cluster_, search_cluster.ptcl_recv_);
        system_hard_conected.driveForMultiClusterOMP(dt_soft);
        search_cluster.writeAndSendBackPtcl(system_soft, system_hard_conected.getPtcl());
        // integrate multi cluster B
        wtime_hard_3rd_block += PS::GetWtime() - wtime_hard_3rd_block_offset;
        wtime_hard_tot += PS::GetWtime() - wtime_hard_tot_offset;
        /////////////

        /////////////
        wtime_soft_tot_offset = PS::GetWtime();
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
                                           CalcForceEpSpNoSIMD(),
                                           system_soft,
                                           dinfo);
        wtime_soft_search_neighbor_offset = PS::GetWtime();
        
        search_cluster.searchNeighborAndCalcHardForceOMP<SystemSoft, Tree, EPJSoft>
          (system_soft, tree_soft, EPISoft::r_out, EPISoft::r_in, pos_domain, EPISoft::eps*EPISoft::eps);
        
        wtime_soft_search_neighbor += PS::GetWtime() - wtime_soft_search_neighbor_offset;
        
        wtime_soft_tot += PS::GetWtime() - wtime_soft_tot_offset;

        // Domain decomposition, parrticle exchange and force calculation
        /////////////

        ////////////////
        ////// 2nd kick
        Kick(system_soft, tree_soft, dt_soft*0.5);
        time_sys += dt_soft;
        ////// 2nd kick
        ////////////////

        eng_now.calc(system_soft, true);
        
        wtime_tot += PS::GetWtime() - wtime_tot_offset;
        /////////////

        Energy eng_diff = eng_now.calcDiff(eng_init);
        PS::S64 n_glb = system_soft.getNumberOfParticleGlobal();
        PS::S32 n_one_cluster         = system_hard_one_cluster.getPtcl().size();
        PS::S32 n_isolated_cluster    = system_hard_isolated.getPtcl().size();
        PS::S32 n_nonisolated_cluster = system_hard_conected.getPtcl().size();
        n_ptcl_hard_one_cluster         += PS::Comm::getSum(n_one_cluster);
        n_ptcl_hard_isolated_cluster    += PS::Comm::getSum(n_isolated_cluster);
        n_ptcl_hard_nonisolated_cluster += PS::Comm::getSum(n_nonisolated_cluster);
        
        dn_loop++;
        
#ifdef ARC_ERROR
        system_hard_isolated.N_count[0] += PS::Comm::getSum(n_one_cluster);
#endif
        
#ifdef DEBUG_OUTPUT
        //output
        PS::S32 ntot = system_soft.getNumberOfParticleLocal();
        fout<<std::setprecision(17)<<time_sys<<" ";
        for (PS::S32 i=0;i<ntot;i++){
          fout<<system_soft[i].mass<<" ";
          for (PS::S32 k=0;k<3;k++) fout<<system_soft[i].pos[k]<<" ";
          for (PS::S32 k=0;k<3;k++) fout<<system_soft[i].vel[k]<<" ";
          fout<<system_soft[i].pot_tot<<" ";
          fout<<0.5*system_soft[i].mass*system_soft[i].vel*system_soft[i].vel<<" ";
        }
        fout<<std::endl;
#endif
        
        if( fmod(time_sys, dt_snp) == 0.0 ){
            //eng_diff.dump(std::cerr);
            std::cerr<<"n_loop= "<<n_loop<<std::endl;
            std::cerr<<"n_glb= "<<n_glb<<std::endl;
            std::cerr<<"Time= "<<time_sys<<" Enow-Einit="<<eng_diff.tot<<" (Enow-Einit)/Einit= "<<eng_diff.tot/eng_init.tot
#ifdef ARC_ERROR
                     <<" ARC_error_relative="<<system_hard_isolated.ARC_error_relative+system_hard_conected.ARC_error_relative<<" ARC_error="<<system_hard_isolated.ARC_error+system_hard_conected.ARC_error
#endif
                     <<std::endl;
            eng_now.dump(std::cerr);
            std::cerr<<"wtime_tot/dn_loop= "<<wtime_tot/dn_loop<<" dn_loop= "<<dn_loop<<std::endl;
            std::cerr<<"wtime_hard_tot/dn_loop= "<<wtime_hard_tot/dn_loop<<std::endl
                     <<"wtime_hard_1st_block/dn_loop= "<<wtime_hard_1st_block/dn_loop<<std::endl
                     <<"wtime_hard_2nd_block/dn_loop= "<<wtime_hard_2nd_block/dn_loop<<std::endl
                     <<"wtime_hard_3rd_block/dn_loop= "<<wtime_hard_3rd_block/dn_loop<<std::endl;
            std::cerr<<"wtime_soft_tot/dn_loop= "<<wtime_soft_tot/dn_loop<<std::endl
                     <<"wtime_soft_search_neighbor/dn_loop= "<<wtime_soft_search_neighbor/dn_loop<<std::endl;
            std::cerr<<"n_ptcl_hard_one_cluster/dn_loop= "         <<(PS::F64)n_ptcl_hard_one_cluster/dn_loop<<std::endl
                     <<"n_ptcl_hard_isolated_cluster/dn_loop= "   <<(PS::F64)n_ptcl_hard_isolated_cluster/dn_loop<<std::endl
                     <<"n_ptcl_hard_nonisolated_cluster/dn_loop= "<<(PS::F64)n_ptcl_hard_nonisolated_cluster/dn_loop<<std::endl;

#ifdef ARC_ERROR
            std::cerr<<"NHist: ";
            for (PS::S32 i=0;i<20;i++) std::cerr<<system_hard_isolated.N_count[i]<<" ";
            std::cerr<<std::endl;
#endif
            
            fprofile<<time_sys<<" "
                    <<n_loop<<" "
                    <<n_glb<<" "
                    <<wtime_tot/dn_loop<<" "
                    <<wtime_hard_tot/dn_loop<<" "
                    <<wtime_hard_1st_block/dn_loop<<" "
                    <<wtime_hard_2nd_block/dn_loop<<" "
                    <<wtime_hard_3rd_block/dn_loop<<" "
                    <<wtime_soft_tot/dn_loop<<" "
                    <<wtime_soft_search_neighbor/dn_loop<<" "
                    <<(PS::F64)n_ptcl_hard_one_cluster/dn_loop<<" "
                    <<(PS::F64)n_ptcl_hard_isolated_cluster/dn_loop<<" "
                    <<(PS::F64)n_ptcl_hard_nonisolated_cluster/dn_loop<<" "
                    <<eng_diff.tot/eng_init.tot<<" "
#ifdef ARC_ERROR
                    <<system_hard_isolated.ARC_error+system_hard_conected.ARC_error<<" "
#endif              
                    <<std::endl;
            
            file_header.time = time_sys;
            file_header.nfile++;
            char sout[99] = "data.";
            sprintf(&sout[5],"%lld",file_header.nfile);
            WriteFile<SystemSoft, FileHeader, FPSoft>(system_soft, sout, file_header);
//            system_soft.writeParticleAscii(sout, file_header);

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
            
            wtime_tot = 0.0;
            wtime_hard_tot = wtime_hard_1st_block = wtime_hard_2nd_block = wtime_hard_3rd_block = 0.0;
            wtime_soft_tot = wtime_soft_search_neighbor = 0.0;
            n_ptcl_hard_one_cluster = n_ptcl_hard_isolated_cluster = n_ptcl_hard_nonisolated_cluster = 0;

            dn_loop=0;
        }

        n_loop++;
    }

#ifdef ARC_ERROR    
    std::cout<<"NHist: ";
    for (PS::S32 i=0;i<20;i++) std::cout<<system_hard_isolated.N_count[i]/(PS::F64)n_loop<<" ";
    std::cout<<std::endl;
#endif

    PS::Finalize();
    return 0;
}

