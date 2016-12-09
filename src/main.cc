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
#include"class.hpp"
#include"hard.hpp"
#include"kepler.hpp"
#include"io.hpp"
#include"profile.hpp"
#include"domain.hpp"
#include"AR.h" /// include AR.h (L.Wang)
//#include"cluster.hpp"
#include"cluster_list.hpp"

/*
template<class Tsys, class Ttree>
void CalcForce(PS::DomainInfo & dinfo, Tsys & system, Ttree & tree){
    const PS::F64vec root_cen(0.0);
    const PS::F64 root_len = GetRootFullLenght(system, root_cen);
    tree.setParticleLocalTree(system);
    tree.setRootCell(root_len, root_cen);
    tree.mortonSortLocalTreeOnly();
    tree.linkCellLocalTreeOnly();
    tree.calcMomentLocalTreeOnly();
    tree.exchangeLocalEssentialTree(dinfo);
    tree.setLocalEssentialTreeToGlobalTree();
    tree.mortonSortGlobalTreeOnly();
    tree.linkCellGlobalTreeOnly();
    tree.calcMomentGlobalTreeOnly();
    tree.makeIPGroup();
#ifdef USE_QUAD
    tree.calcForce(CalcForceEpEpWithLinearCutoff(), CalcForceEPSPQuad(), true);
#else
    tree.calcForce(CalcForceEpEpWithLinearCutoff(), CalcForceEPSPMono(), true);
#endif
    const PS::S32 n_loc = system.getNumberOfParticleLocal();
#pragma omp parallel for
    for(PS::S32 i=0; i<n_loc; i++){
        system[i].copyFromForce(tree.getForce(i));
    }
}
*/

PS::F64 GetQuantizedValue(const PS::F64 & val){
    static const PS::F64 inv_log2 = 1.0 / log(2.0);
    PS::F64 log2val = log(val) * inv_log2;
    log2val = (log2val > 0.0) ? log2val : log2val - 1.0;
    PS::S32 power = (PS::S32)(log2val);
    return pow(2.0, (PS::F64)(power));

}

bool GetFlagSnpWithEnergy(const char sinput[]){
    bool flag = true;
    if(PS::Comm::getRank() == 0){
	std::ifstream fin;
	fin.open(sinput);
	std::cout<<"sinput:"<<sinput<<std::endl;
	std::string line;
	getline(fin, line);
	std::stringstream ss(line);
	std::string str;
	PS::S32 ncnt = 0;
	while(ss>>str){
	    std::cout<<"str: "<<str<<std::endl;
	    ncnt++;
	}
	std::cout<<"line="<<line<<std::endl;
	std::cout<<"ncnt="<<ncnt<<std::endl;
	if(ncnt == 2){ flag = false; }
	else if(ncnt == 12){ flag = true; }
	else{
	    std::cerr<<"input file is wrong format"<<std::endl;
	    PS::Abort();
	}
    }
    PS::Comm::broadcast(&flag, 1, 0);
    return flag;
}

template<class Tpsys>
PS::F64 GetRootFullLenght(const Tpsys & psys, const PS::F64vec & cen){
    PS::S64 nloc = psys.getNumberOfParticleLocal();
    PS::F64 len_loc_max = 0.0;
    for(PS::S32 i=0; i<nloc; i++){
	PS::F64vec dr = psys[i].pos - cen;
	for(PS::S32 k=0; k<3; k++){
	    if(len_loc_max < dr[k]) len_loc_max = dr[k];
	}
    }
    return 2.1*fabs(PS::Comm::getMaxValue(len_loc_max));
}

class FileHeader{
public:
    PS::S64 n_body;
    PS::F64 time;
    FileHeader(){
        n_body = 0;
        time = 0.0;
    }
    FileHeader(const PS::S64 n, const PS::F64 t){
        n_body = n;
        time = t;
    }
    PS::S32 readAscii(FILE * fp){
        fscanf(fp, "%lld\t%lf\n", &n_body, &time);
	std::cout<<"n_body="<<n_body<<" time="<<time<<std::endl;
        return n_body;
    }
    void writeAscii(FILE* fp) const{
        fprintf(fp, "%lld\t%lf\n", n_body, time);
    }
};

class FileHeaderWithEnergy{
public:
    PS::S64 n_body;
    PS::F64 time;
    Energy eng_init;
    Energy eng_now;
    FileHeaderWithEnergy(){
        n_body = 0;
        time = 0.0;
        eng_init.clear();
        eng_now.clear();
    }
    FileHeaderWithEnergy(const PS::S64 n, const PS::F64 t, const Energy & e_i, const Energy & e_n){
        n_body = n;
        time = t;
	eng_init = e_i;
	eng_now = e_n;
    }
    PS::S32 readAscii(FILE * fp){
        fscanf(fp, "%lld%lf %lf%lf%lf %lf%lf%lf\n", 
	       &n_body, &time, 
	       &eng_init.kin, &eng_init.pot, &eng_init.tot,
	       &eng_now.kin,  &eng_now.pot,  &eng_now.tot);
	std::cout<<"n_body="<<n_body<<" time="<<time<<std::endl;
        return n_body;
    }
    void writeAscii(FILE* fp) const{
        fprintf(fp, "%lld\t%lf\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n", 
		n_body, time, 
		eng_init.kin, eng_init.pot, eng_init.tot,
		eng_now.kin,  eng_now.pot,  eng_now.tot);
    }

};

template<class Tpsys>
PS::F64 GetMassMax(const Tpsys & system, const PS::S64 n){
    PS::F64 m_max_loc = -1.0;
    for(PS::S64 i=0; i<n; i++){
        if(m_max_loc < system[i].mass) m_max_loc = system[i].mass;
    }
    return PS::Comm::getMaxValue(m_max_loc);
}

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


template<class Tpsys>
PS::F64 GetRSearch(const Tpsys & system_soft,
		   const PS::F64 r_out,
		   const PS::F64 search_factor,
		   const PS::F64 dt){
    const PS::S32 n_loc = system_soft.getNumberOfParticleLocal();
    PS::F64vec vel_cm_loc = 0.0;
    PS::F64 mass_cm_loc = 0.0;
    for(PS::S32 i=0; i<n_loc; i++){
	mass_cm_loc += system_soft[i].mass;
	vel_cm_loc += system_soft[i].mass * system_soft[i].vel;
    }
    PS::F64    mass_cm_glb = PS::Comm::getSum(mass_cm_loc);
    PS::F64vec vel_cm_glb  = PS::Comm::getSum(vel_cm_loc);
    vel_cm_glb /= mass_cm_glb;
    PS::F64 vel_sq_loc = 0.0;
    for(PS::S32 i=0; i<n_loc; i++){
	PS::F64vec dv = system_soft[i].vel - vel_cm_glb;
	vel_sq_loc += dv * dv;
    }
    const PS::S32    n_glb      = PS::Comm::getSum(n_loc);
    const PS::F64 vel_sq_glb = PS::Comm::getSum(vel_sq_loc);
    const PS::F64    vel_disp   = sqrt(vel_sq_glb) / n_glb;
    return r_out + vel_disp*dt*search_factor;
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
	PS::F64 wtime_soft_force = 0.0;
	PS::F64 wtime_soft_search_neighbor_offset = 0.0;
	PS::F64 wtime_soft_search_neighbor = 0.0;

	PS::S64 n_ptcl_hard_one_cluster = 0;
	PS::S64 n_ptcl_hard_isolated_cluster = 0;
	PS::S64 n_ptcl_hard_nonisolated_cluster = 0;

#ifdef CALC_HARD_ENERGY
    PS::F64 dEerr_1body_loc = 0.0;
    PS::F64 dEerr_1body_glb = 0.0;
    PS::F64 dEerr_2body_loc = 0.0;
    PS::F64 dEerr_2body_glb = 0.0;
    PS::F64 dEerr_mbody_loc = 0.0;
    PS::F64 dEerr_mbody_glb = 0.0;
#endif
    PS::F64 ratio_r_cut = 0.1;
    PS::F64 time_sys = 0.0;
    PS::F64 theta = 0.4;
    PS::S32 n_leaf_limit = 8;
    PS::S32 n_group_limit = 64;
    PS::S32 n_smp_ave = 100;
    PS::F64 time_end = 1000.0;
    PS::F64 dt_soft = 1.0/256.0;
    PS::F64 eta = 0.1;
    PS::F64 eta_s = eta * 0.1;
    char dir_name[1024];
    PS::S64 n_glb = 16384;
    PS::F64 dt_snp = 1.0 / 16.0;
    PS::F64 search_factor = 3.0; 
    PS::F64 dt_limit_hard_factor = 4.0;
    PS::S64 snp_id = 0;
    PS::F64 r_out = 1.0 / 256.0;
    PS::F64 eps = 1.0 / 64.0;
    int c;
#ifdef READ_FILE
    char sinput[2048];
    while((c=getopt(argc,argv,"i:I:o:d:D:E:t:T:n:s:S:l:X:h")) != -1){
        switch(c){
        case 'i':
            sprintf(sinput, optarg);
            std::cerr<<"sinput="<<sinput<<std::endl;
            break;
        case 'I':
            snp_id = atoi(optarg);
            std::cerr<<"snp_id="<<snp_id<<std::endl;
            break;
        case 'o':
            sprintf(dir_name, optarg);
            std::cerr<<"dir_name="<<dir_name<<std::endl;
            break;
        case 'd':
            dt_soft = 1.0 / atof(optarg);
            std::cerr<<"dt_soft="<<dt_soft<<std::endl;
            break;
        case 'D':
            dt_snp = atof(optarg);
            std::cerr<<"dt_snp="<<dt_snp<<std::endl;
            break;
        case 'E':
            eta = atof(optarg);
            std::cerr<<"eta="<<eta<<std::endl;
            break;
        case 't':
            theta = atof(optarg);
            std::cerr<<"theta="<<theta<<std::endl;
            break;
        case 'T':
            time_end = atof(optarg);
            std::cerr<<"time_end="<<time_end<<std::endl;
            break;
        case 'n':
            n_group_limit = atoi(optarg);
            std::cerr<<"n_group_limit="<<n_group_limit<<std::endl;
            break;
        case 's':
            n_smp_ave = atoi(optarg);
            std::cerr<<"n_smp_ave="<<n_smp_ave<<std::endl;
            break;
        case 'S':
            search_factor = atoi(optarg);
            std::cerr<<"search_factor="<<search_factor<<std::endl;
            break;
        case 'l':
            n_leaf_limit = atoi(optarg);
            std::cerr<<"n_leaf_limit="<<n_leaf_limit<<std::endl;
            break;
        case 'X':
	    dt_limit_hard_factor = atof(optarg);
            std::cerr<<"dt_limit_hard_factor="<<dt_limit_hard_factor<<std::endl;
	    assert(dt_limit_hard_factor > 0.0);
            break;
        case 'h':
            std::cerr<<"i: input_file"<<std::endl;
            std::cerr<<"I: snp_id"<<std::endl;
            std::cerr<<"o: dir name of output"<<std::endl;
            std::cerr<<"d: inv_dt (dafult 16 ~ 0.01yr  )"<<std::endl;
            std::cerr<<"D: dt_snp (dafult 100 ~ 16yr  )"<<std::endl;
            std::cerr<<"E: eta (dafult 0.1)"<<std::endl;
            std::cerr<<"t: theta (dafult: 0.5)"<<std::endl;
            std::cerr<<"T: time_end (dafult: 6000 ~ 1000yr)"<<std::endl;
            std::cerr<<"n: n_group_limit (dafult: 64.0)"<<std::endl;
            std::cerr<<"s: n_smp_ave (dafult: 100)"<<std::endl;
            std::cerr<<"S: search_factor (dafult: 3.0)"<<std::endl;
            std::cerr<<"l: n_leaf_limit (dafult: 8)"<<std::endl;
            std::cerr<<"X: dt_limit_hard_factor(dafult: 4.0 -> dt_limit_hard = dt_soft/4.0)"<<std::endl;
            PS::Finalize();
            return 0;
        }
    }
#else // not READ_FILE
    PS::S32 seed = 0;
    while((c=getopt(argc,argv,"e:o:d:D:E:t:T:n:N:s:S:l:r:r:x:X:h")) != -1){
        switch(c){
        case 'e':
            eps = atof(optarg);
            break;
        case 'o':
            sprintf(dir_name,optarg);
            break;
        case 'd':
            dt_soft = 1.0 / atof(optarg);
            break;
        case 'D':
            dt_snp = atof(optarg);
            break;
        case 'E':
            eta = atof(optarg);
            break;
        case 't':
            theta = atof(optarg);
            break;
        case 'T':
            time_end = atof(optarg);
            break;
        case 'n':
            n_group_limit = atoi(optarg);
            break;
        case 'N':
            n_glb = atol(optarg);
            break;
        case 's':
            n_smp_ave = atoi(optarg);
            break;
        case 'S':
            search_factor = atoi(optarg);
            break;
        case 'l':
            n_leaf_limit = atoi(optarg);
            break;
        case 'r':
            r_out = atof(optarg);
        case 'x':
            seed = atoi(optarg);
            break;
        case 'X':
            dt_limit_hard_factor = atof(optarg);
	    assert(dt_limit_hard_factor > 0.0);
            break;
        case 'h':
	    std::cerr<<"e: softening (dafule 1/64)"<<std::endl;
            std::cerr<<"o: dir name of output"<<std::endl;
            std::cerr<<"d: inv_dt (dafult 16 ~ 0.01yr  )"<<std::endl;
            std::cerr<<"D: dt_snp (dafult 100 ~ 16yr  )"<<std::endl;
            std::cerr<<"E: eta (dafult 0.1)"<<std::endl;
            std::cerr<<"t: theta (dafult: 0.5)"<<std::endl;
            std::cerr<<"T: time_end (dafult: 8000 ~ 1yr)"<<std::endl;
            std::cerr<<"n: n_group_limit (dafult: 64.0)"<<std::endl;
            std::cerr<<"N: n_glb (dafult: 16384)"<<std::endl;
            std::cerr<<"s: n_smp_ave (dafult: 100)"<<std::endl;
            std::cerr<<"S: search_factor (dafult: 3.0)"<<std::endl;
            std::cerr<<"l: n_leaf_limit (dafult: 8)"<<std::endl;
            std::cerr<<"r: r_out (dafult: 1/256)"<<std::endl;
            std::cerr<<"x: seed (dafult: 0)"<<std::endl;
            std::cerr<<"X: dt_limit_hard_factor(dafult: 4.0 -> dt_limit_hard = dt_soft/4.0)"<<std::endl;
            PS::Finalize();
            return 0;
        }
    }
#endif//READ_FILE

    SystemSoft system_soft;
    system_soft.initialize();
    system_soft.setAverageTargetNumberOfSampleParticlePerProcess(n_smp_ave);
    PS::S32 n_loc;

    SetParticlePlummer(system_soft, n_glb, n_loc, time_sys);

    PS::Comm::barrier();

    PS::F64 r_in = r_out * ratio_r_cut;
    PS::F64 r_search = GetRSearch(system_soft, r_out, search_factor, dt_soft);
    std::cerr<<"r_search= "<<r_search<<std::endl;
    EPJSoft::r_search = r_search;
    EPISoft::r_out    = r_out;
    EPISoft::r_in     = r_in;
    EPISoft::eps      = eps;

    const PS::F32 coef_ema = 0.2;
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);
    dinfo.decomposeDomainAll(system_soft);
    const PS::S32 n_proc = PS::Comm::getNumberOfProc();
    PS::F64ort * pos_domain = new PS::F64ort[n_proc];
    for(PS::S32 i=0; i<n_proc; i++) pos_domain[i] = dinfo.getPosDomain(i);

    system_soft.exchangeParticle(dinfo); 

    n_loc = system_soft.getNumberOfParticleLocal();

    PS::S32 my_rank = PS::Comm::getRank();
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
    system_hard_one_cluster.setParam(r_out, r_in, dt_soft/dt_limit_hard_factor, time_sys);
    system_hard_one_cluster.setARCParam();
    SystemHard system_hard_isolated;
    system_hard_isolated.setParam(r_out, r_in, dt_soft/dt_limit_hard_factor, time_sys);
    system_hard_isolated.setARCParam();
    SystemHard system_hard_conected;
    system_hard_conected.setParam(r_out, r_in, dt_soft/dt_limit_hard_factor, time_sys);
    system_hard_conected.setARCParam();

    SearchCluster search_cluster;
    search_cluster.initialize();
    search_cluster.searchNeighborAndCalcHardForceOMP<SystemSoft, Tree, EPJSoft>
	(system_soft, tree_soft, r_out, r_in, pos_domain, EPISoft::eps*EPISoft::eps);

    Energy eng_init, eng_now;
#ifndef READ_FILE
    eng_init.calc(system_soft, true);
#endif
    eng_init.dump(std::cerr);
    eng_now = eng_init;

    PS::S32 n_loop = 0;
    PS::S32 dn_loop = 0;
    while(time_sys < time_end){
        wtime_tot_offset = PS::GetWtime();
        ////////////////
        ////// 1st kick
        Kick(system_soft, tree_soft, dt_soft*0.5);

        wtime_hard_tot_offset = PS::GetWtime();
        
        system_hard_one_cluster.setTimeOrigin(time_sys);
        system_hard_isolated.setTimeOrigin(time_sys);
        system_hard_conected.setTimeOrigin(time_sys);
        ////// 1st kick
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
        ////// integrater one cluster
        wtime_hard_1st_block_offset = PS::GetWtime();

        system_hard_one_cluster.initializeForOneCluster(search_cluster.getAdrSysOneCluster().size());

        system_hard_one_cluster.setPtclForOneCluster(system_soft, search_cluster.getAdrSysOneCluster());
        system_hard_one_cluster.driveForOneCluster(dt_soft);
        system_hard_one_cluster.writeBackPtclForOneCluster(system_soft, search_cluster.getAdrSysOneCluster());

        wtime_hard_1st_block += PS::GetWtime() - wtime_hard_1st_block_offset;

        ////// integrater one cluster
        ////////////////


        //std::cerr<<"check a"<<std::endl;

        /////////////
        // integrate multi cluster A
        wtime_hard_2nd_block_offset = PS::GetWtime();

        system_hard_isolated.setPtclForIsolatedMultiCluster(system_soft,
                                                            search_cluster.adr_sys_multi_cluster_isolated_,
                                                            search_cluster.n_ptcl_in_multi_cluster_isolated_);
        system_hard_isolated.driveForMultiCluster(dt_soft);
        system_hard_isolated.writeBackPtclForMultiCluster(system_soft, search_cluster.adr_sys_multi_cluster_isolated_);

        wtime_hard_2nd_block += PS::GetWtime() - wtime_hard_2nd_block_offset;
        // integrate multi cluster A
        /////////////

        //std::cerr<<"check b"<<std::endl;

        /////////////
        // integrate multi cluster B
        wtime_hard_3rd_block_offset = PS::GetWtime();

        system_hard_conected.setPtclForConectedCluster(system_soft, search_cluster.mediator_sorted_id_cluster_, search_cluster.ptcl_recv_);
        system_hard_conected.driveForMultiClusterOMP(dt_soft);
        search_cluster.writeAndSendBackPtcl(system_soft, system_hard_conected.getPtcl());

        wtime_hard_3rd_block += PS::GetWtime() - wtime_hard_3rd_block_offset;
        // integrate multi cluster B
        /////////////

        //std::cerr<<"check c"<<std::endl;

        wtime_hard_tot += PS::GetWtime() - wtime_hard_tot_offset;
        /////////////
        // Domain decomposition, parrticle exchange and force calculation

        wtime_soft_tot_offset = PS::GetWtime();

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
            (system_soft, tree_soft, r_out, r_in, pos_domain, EPISoft::eps*EPISoft::eps);
        
        wtime_soft_search_neighbor += PS::GetWtime() - wtime_soft_search_neighbor_offset;
        
        wtime_soft_tot += PS::GetWtime() - wtime_soft_tot_offset;

        // Domain decomposition, parrticle exchange and force calculation
        /////////////

        //std::cerr<<"check d"<<std::endl;

        ////////////////
        ////// 2nd kick
        Kick(system_soft, tree_soft, dt_soft*0.5);
        time_sys += dt_soft;
        ////// 2nd kick
        ////////////////

        eng_now.calc(system_soft, true);
        
        wtime_tot += PS::GetWtime() - wtime_tot_offset;

        Energy eng_diff = eng_now.calcDiff(eng_init);
        PS::S64 n_glb = system_soft.getNumberOfParticleGlobal();
        PS::S32 n_one_cluster         = system_hard_one_cluster.getPtcl().size();
        PS::S32 n_isolated_cluster    = system_hard_isolated.getPtcl().size();
        PS::S32 n_nonisolated_cluster = system_hard_conected.getPtcl().size();
        n_ptcl_hard_one_cluster         += PS::Comm::getSum(n_one_cluster);
        n_ptcl_hard_isolated_cluster    += PS::Comm::getSum(n_isolated_cluster);
        n_ptcl_hard_nonisolated_cluster += PS::Comm::getSum(n_nonisolated_cluster);
        
        dn_loop++;

        if( fmod(time_sys, dt_snp) == 0.0 ){
            //eng_diff.dump(std::cerr);
            std::cerr<<"n_loop= "<<n_loop<<std::endl;
            std::cerr<<"n_glb= "<<n_glb<<std::endl;
            std::cerr<<"time_sys= "<<time_sys<<" (Enow-Einit)/Einit= "<<eng_diff.tot/eng_init.tot<<std::endl;
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
            wtime_tot = 0.0;
            wtime_hard_tot = wtime_hard_1st_block = wtime_hard_2nd_block = wtime_hard_3rd_block = 0.0;
            wtime_soft_tot = wtime_soft_search_neighbor = 0.0;
            n_ptcl_hard_one_cluster = n_ptcl_hard_isolated_cluster = n_ptcl_hard_nonisolated_cluster = 0;

            dn_loop=0;
        }

        n_loop++;
    }


    PS::Finalize();
    return 0;
}

