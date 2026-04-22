#include <iostream>
#include <iomanip>
#include <string>
#include<getopt.h>

#include <particle_simulator.hpp>
#define HARD_DEBUG_PRINT_FEQ 1

#include "io.hpp"
#include "hard_assert.hpp"
#include "cluster_list.hpp"
#include "hard.hpp"
#include "soft_ptcl.hpp"
#include "static_variables.hpp"
#include "status.hpp"

#if defined(BSE_BASE) || defined(DISK_STAR_MERGER)
#include "../parallel-random/rand_io.hpp"
#include "../parallel-random/rand_interface.hpp"
#endif

class IOParamsHardDebug {
public:
    IOParamsContainer input_par_store;
    IOParams<PS::S64> mode;
    IOParams<PS::S64> n_crit_ptcl;
    IOParams<PS::F64> tstart;
    IOParams<PS::F64> tend;
    IOParams<PS::S64> istart;
    IOParams<PS::S64> iend;
    IOParams<PS::S64> n_crit_group;
    IOParams<PS::S64> n_crit_arti;
#ifdef SOFT_PERT
    IOParams<PS::S64> soft_perturbation;
#endif
    IOParams<std::string> fname_par;
    IOParams<std::string> fname_dump;
    bool print_flag;

    IOParamsHardDebug(): input_par_store(),
                         mode(input_par_store, 0, "m", "running mode; 0: evolve system to time_end; 1: stability check"),
                         n_crit_ptcl(input_par_store, 0, "n", "if >0, only do integration when particle number matches the given value"),
                         tstart(input_par_store, -1.0, "tstart", "if >0 only do integration when physical time >= tstart"),
                         tend(input_par_store, -1.0, "tend", "if >0 only do integration when physical time < tend"),
                         istart(input_par_store, -1, "istart", "if >0 only do integration when dump index >= istart (counting from 1)"),
                         iend(input_par_store, -1, "iend", "if >0 only do integration when dump index < iend (counting from 1)"),
                         n_crit_group(input_par_store, 0, "n-crit-group", "if >0 only do integration when group number matches the given value"),
                         n_crit_arti(input_par_store, 0, "n-crit-arti", "if >0 only do integration when artificial particle number matches the given value"),
#ifdef SOFT_PERT
                         soft_perturbation(input_par_store, 1, "soft-perturbation", "switch of soft perturbation (tidal tensor); 1: enable, 0: suppress"),
#endif
                         fname_par(input_par_store, "input.par", "p", "Parameter file prefix; related files are <prefix>.hard/.bse/.dsm/.exthard/.galpy/.rand", NULL, false),
                         fname_dump(input_par_store, "hard_dump", "dump-filename", "Hard dump data filename", NULL, false),
                         print_flag(false)
    {}

    int read(int argc, char* argv[], const int opt_used_pre=0) {
        static int debug_flag=-1;
        static struct option long_options[] = {
            {tstart.key, required_argument, &debug_flag, 0},
            {tend.key, required_argument, &debug_flag, 1},
            {istart.key, required_argument, &debug_flag, 2},
            {iend.key, required_argument, &debug_flag, 3},
            {n_crit_group.key, required_argument, &debug_flag, 4},
            {n_crit_arti.key, required_argument, &debug_flag, 5},
#ifdef SOFT_PERT
            {soft_perturbation.key, required_argument, &debug_flag, 6},
#endif
            {"help", no_argument, 0, 'h'},
            {0,0,0,0}
        };

        int opt_used = opt_used_pre;
        int copt;
        int option_index;
        optind = 0;
        while ((copt = getopt_long(argc, argv,
                                   "-m:n:p:h",
                                   long_options, &option_index)) != -1) {
            switch (copt) {
            case 0:
                switch (debug_flag) {
                case 0:
                    tstart.value = atof(optarg);
                    if (print_flag) tstart.print(std::cout);
                    opt_used += 2;
                    break;
                case 1:
                    tend.value = atof(optarg);
                    if (print_flag) tend.print(std::cout);
                    opt_used += 2;
                    break;
                case 2:
                    istart.value = atoi(optarg);
                    if (print_flag) istart.print(std::cout);
                    opt_used += 2;
                    break;
                case 3:
                    iend.value = atoi(optarg);
                    if (print_flag) iend.print(std::cout);
                    opt_used += 2;
                    break;
                case 4:
                    n_crit_group.value = atoi(optarg);
                    if (print_flag) n_crit_group.print(std::cout);
                    opt_used += 2;
                    break;
                case 5:
                    n_crit_arti.value = atoi(optarg);
                    if (print_flag) n_crit_arti.print(std::cout);
                    opt_used += 2;
                    break;
#ifdef SOFT_PERT
                case 6:
                    soft_perturbation.value = atoi(optarg);
                    if (print_flag) soft_perturbation.print(std::cout);
                    opt_used += 2;
                    break;
#endif
                default:
                    break;
                }
                break;
            case 'm':
                mode.value = atoi(optarg);
                if (print_flag) mode.print(std::cout);
                opt_used += 2;
                break;
            case 'n':
                n_crit_ptcl.value = atoi(optarg);
                if (print_flag) n_crit_ptcl.print(std::cout);
                opt_used += 2;
                break;
            case 'p':
                fname_par.value = optarg;
                if (print_flag) fname_par.print(std::cout);
                opt_used += 2;
                break;
            case 'h':
                if (print_flag) {
                    std::cout<<"A tool to integrate a dumped cluster of neighbor particles using particle-particle method (Hermite/SDAR)\n"
                         <<"Usage: petar.hard.debug [options] [hard parameter filename (defaulted: input.par or input.par.hard)] [dumped data filename (defaulted: hard_dump)]\n"
                         <<"   dumped data file: Hard dump file from petar simulation\n"
                         <<"---- Main Options---- :\n"
                         <<"    -h (--help):  help\n";
                        input_par_store.printHelp(std::cout, true, false);
                }
                return -1;
            case '?':
                opt_used += 2;
                break;
            default:
                break;
            }
        }

        opt_used++;
        if (opt_used < argc) {
            fname_dump.value = argv[argc-1];
            if (print_flag) std::cout<<"Reading dump data filename: "<<fname_dump.value<<std::endl;
        }
        return opt_used-1;
    }

    bool checkParams() {
        assert(mode.value==0 || mode.value==1);
        assert(n_crit_ptcl.value >= 0 && n_crit_group.value >= 0 && n_crit_arti.value >= 0);
        assert(istart.value==-1 || istart.value>=1);
        assert(iend.value==-1 || iend.value>=1);
    #ifdef SOFT_PERT
        assert(soft_perturbation.value==0 || soft_perturbation.value==1);
    #endif
        assert(!(istart.value>0 && iend.value>0 && iend.value<istart.value));
        assert(!(tstart.value>0 && tend.value>0 && tend.value<=tstart.value));
        return true;
    }

};

int main(int argc, char **argv){
    IOParamsHardDebug debug_io;
    IOParamsHard hard_io;
#ifdef STELLAR_EVOLUTION
#ifdef BSE_BASE
    IOParamsBSE bse_io;
    std::string bse_name = BSEManager::getBSEName();
    std::string fsse_suffix = BSEManager::getSSEOutputFilenameSuffix();
    std::string fbse_suffix = BSEManager::getBSEOutputFilenameSuffix();
#elif DISK_STAR_MERGER
    IOParamsDiskStarMerger dsm_io;
#endif
#if defined(BSE_BASE) || defined(DISK_STAR_MERGER)
    IOParamsRand rand_io;
    uint64_t seed_override = 0;
    bool has_rand_seed_option = false;
#endif
#endif
#ifdef EXTERNAL_HARD
    IOParamsExternalHard ext_hard_io;
#ifdef GALPY
    IOParamsGalpy galpy_io;
#endif
#endif

    opterr = 0;

    auto has_suffix = [](const std::string& _str, const std::string& _suffix) {
        return _str.size() >= _suffix.size() && _str.compare(_str.size()-_suffix.size(), _suffix.size(), _suffix) == 0;
    };

    // Inject default -p early so explicit CLI options that appear later keep precedence.
    std::vector<std::string> adjusted_args;
    adjusted_args.reserve(argc + 2);
    bool has_p_option = false;
    adjusted_args.push_back(argv[0]);
    for (int i=1; i<argc; i++) adjusted_args.push_back(argv[i]);
    for (int i=1; i<argc; i++) {
        if (adjusted_args[i] == "-p" && i+1<argc) {
            has_p_option = true;
            if (has_suffix(adjusted_args[i+1], ".hard")) {
                adjusted_args[i+1] = adjusted_args[i+1].substr(0, adjusted_args[i+1].size()-5);
            }
        }
#if defined(BSE_BASE) || defined(DISK_STAR_MERGER)
        if (adjusted_args[i] == "--rand-seed" || adjusted_args[i] == "--rand-seedfile") has_rand_seed_option = true;
#endif
    }
    if (!has_p_option) {
        adjusted_args.insert(adjusted_args.begin()+1, "-p");
        adjusted_args.insert(adjusted_args.begin()+2, debug_io.fname_par.value);
    }
    std::vector<char*> adjusted_argv;
    adjusted_argv.reserve(adjusted_args.size());
    for (auto& item: adjusted_args) adjusted_argv.push_back(item.data());

    debug_io.print_flag=true;    
    const bool help_flag = (debug_io.read(adjusted_args.size(), adjusted_argv.data()) == -1);
    if (!help_flag) debug_io.checkParams();

    hard_io.print_flag = true;
    hard_io.read(adjusted_args.size(), adjusted_argv.data(), false);
#ifdef BSE_BASE
    bse_io.print_flag = true;
    bse_io.read(adjusted_args.size(), adjusted_argv.data(), false);
#endif
#ifdef DISK_STAR_MERGER
    dsm_io.print_flag = true;
    dsm_io.read(adjusted_args.size(), adjusted_argv.data(), false);
#endif
#if defined(BSE_BASE) || defined(DISK_STAR_MERGER)
    rand_io.print_flag = true;
    rand_io.read(adjusted_args.size(), adjusted_argv.data(), false);
    if (has_rand_seed_option && rand_io.seed.value > 0) seed_override = static_cast<uint64_t>(rand_io.seed.value);
#endif
#ifdef EXTERNAL_HARD
    ext_hard_io.print_flag = true;
    ext_hard_io.read(adjusted_args.size(), adjusted_argv.data(), false);
#ifdef GALPY
    galpy_io.print_flag = true;
    galpy_io.read(adjusted_args.size(), adjusted_argv.data(), false);
#endif
#endif

    if (help_flag) return 0;
    
    const std::string filename = debug_io.fname_dump.value;
    const std::string fhardpar = debug_io.fname_par.value;
    const std::string fhardpar_ascii = has_suffix(fhardpar, ".hard") ? fhardpar.substr(0, fhardpar.size()-5) : fhardpar;
    const std::string fhardpar_report = fhardpar_ascii + ".hard";

    std::cerr<<"Reading dump file:"<<filename<<std::endl;
    std::cerr<<"Hard manager parameter file:"<<fhardpar_report<<std::endl;
    std::cerr<<"Hard manager read mode:ascii-only"<<std::endl;

    std::cout<<std::setprecision(WRITE_PRECISION);

    std::vector<IOParamsContainer*> all_pars;
    all_pars.push_back(&debug_io.input_par_store);
    all_pars.push_back(&hard_io.input_par_store);
#ifdef BSE_BASE
    all_pars.push_back(&bse_io.input_par_store);
#endif
#ifdef DISK_STAR_MERGER
    all_pars.push_back(&dsm_io.input_par_store);
#endif
#if defined(BSE_BASE) || defined(DISK_STAR_MERGER)
    all_pars.push_back(&rand_io.input_par_store);
#endif
#ifdef EXTERNAL_HARD
    ext_hard_io.appendInputParamStores(all_pars);
#ifdef GALPY
    all_pars.push_back(&galpy_io.input_par_store);
#endif
#endif
    std::vector<std::string> known_options;
    known_options.push_back("help");
    known_options.push_back("h");
    FindUndefinedOptions(all_pars, static_cast<int>(adjusted_argv.size()), adjusted_argv.data(), &known_options);

    HardManager hard_manager;
    Status stat;
    hard_manager.status = &stat;

#ifdef STELLAR_EVOLUTION
#ifdef BSE_BASE
    std::cerr<<bse_name<<" parameters source prefix:"<<fhardpar_ascii<<std::endl;
#else // BSE_BASE
#endif 
#ifdef DISK_STAR_MERGER
    std::cerr<<"DSM parameters source prefix:"<<fhardpar_ascii<<std::endl;
#endif
#endif //STELLAR_EVOLUTION

#ifdef EXTERNAL_HARD
    std::cerr<<"External hard parameters source prefix:"<<fhardpar_ascii<<std::endl;

#ifdef GALPY
    std::cerr<<"Galpy parameters source prefix:"<<fhardpar_ascii<<std::endl;

#endif
#endif        

    // Reconstruct hard manager from IO parameters instead of binary state.
    // With the latest initial() behavior, acc_offset_sq is persisted after first auto-computation.
    if (hard_io.acc_offset_sq.value < 0.0) {
        std::cerr<<"Error: hard_io.acc_offset_sq < 0 in ASCII mode. Please provide a hard parameter file with persisted hermite-acc-offset-sq.\n";
        abort();
    }

    if (hard_io.dt_max_hermite.value <=0.0) {
        std::cerr<<"Error: hard_io.dt_max_hermite <= 0 in ASCII mode. Please provide a hard parameter file with proper hermite-dt-max.\n";
        abort();
    }
    stat.time = 0.0;
    // bse_manager/disk_star_merger_manager are initialized inside HardManager::initial.
#ifdef BSE_BASE
    hard_manager.initial(hard_io, bse_io, stat, 1, false);
#elif DISK_STAR_MERGER
    hard_manager.initial(hard_io, dsm_io, stat, 1, false);
#else
    hard_manager.initial(hard_io, stat, 1, false);
#endif

#ifdef STELLAR_EVOLUTION
#ifdef BSE_BASE
    if (hard_manager.ar_manager.interaction.stellar_evolution_write_flag) {
        hard_manager.ar_manager.interaction.fout_sse.open((filename+fsse_suffix).c_str(), std::ofstream::out);
        hard_manager.ar_manager.interaction.fout_bse.open((filename+fbse_suffix).c_str(), std::ofstream::out);
        hard_manager.ar_manager.interaction.fout_sse<<std::setprecision(WRITE_PRECISION);
        hard_manager.ar_manager.interaction.fout_bse<<std::setprecision(WRITE_PRECISION);
    }
#else
    if (hard_manager.ar_manager.interaction.interrupt_detection_option>0) {
        hard_manager.ar_manager.interaction.fout_interrupt.open((filename+".interrupt").c_str(), std::ofstream::out);
        hard_manager.ar_manager.interaction.fout_interrupt<<std::setprecision(WRITE_PRECISION);
    }
#endif
#endif


#ifdef ADJUST_GROUP_PRINT
    if (hard_manager.h4_manager.group_info_output.isWriteEnabled()) {
        const bool binary_flag = (hard_io.adjust_group_write_option.value==2);
        hard_manager.h4_manager.group_info_output.setup(filename+".group", false, binary_flag, WRITE_PRECISION);
    }
#endif

    hard_manager.checkParams();
    hard_manager.print(std::cerr);

    std::FILE* fp = std::fopen(filename.c_str(),"rb");
    if (fp==NULL) {
        std::cerr<<"Error: filename "<<filename<<" cannot be open!\n";
        abort();
    }

    HardDump hard_dump;
    int ncount = 0;
    
    while (true) {
        int c = fgetc(fp);
        if (c == EOF) break;
        ungetc(c, fp);
        hard_dump.readOneClusterBinary(fp);
#ifdef EXTERNAL_HARD
        hard_manager.center.readBinary(fp);
#endif

#if defined(BSE_BASE) || defined(DISK_STAR_MERGER)
        if (seed_override!=0) hard_dump.rand_manager.initialFromSeed(seed_override, 0);
#endif

        ncount++;

        // skip if particle/group/artificial particle number not match n_crit_**
        if (debug_io.n_crit_ptcl.value>0 && hard_dump.n_ptcl != debug_io.n_crit_ptcl.value) continue;
        if (debug_io.n_crit_group.value>0 && hard_dump.n_group != debug_io.n_crit_group.value) continue;
        if (debug_io.n_crit_arti.value>0 && hard_dump.n_arti != debug_io.n_crit_arti.value) continue;
        // skip if time is out of range
        if (debug_io.tstart.value>0 && hard_dump.time_offset < debug_io.tstart.value) continue;
        if (debug_io.tend.value>0 && hard_dump.time_offset >= debug_io.tend.value) continue;
        // skip if dump index is out of range
        if (ncount < debug_io.istart.value) continue;
        if (debug_io.iend.value>0 && ncount>debug_io.iend.value) continue;

        std::cerr<<"Dump "<<ncount<<"\nTime: "<<hard_dump.time_offset<<std::endl;
#if defined(BSE_BASE) || defined(DISK_STAR_MERGER)
        hard_dump.rand_manager.printRandSeeds(std::cerr);
#endif

        stat.time = hard_dump.time_offset;
        stat.pcm.mass = hard_dump.gcm_mass;
        stat.pcm.pos = hard_dump.gcm_pos;
        stat.pcm.vel = hard_dump.gcm_vel;

        std::cerr<<"Global CM: mass="<<stat.pcm.mass<<" pos="<<stat.pcm.pos<<" vel="<<stat.pcm.vel<<std::endl;  

#ifdef SOFT_PERT
        if (debug_io.soft_perturbation.value==0) {
            std::cerr<<"Suppress soft perturbation\n";
            if (hard_dump.n_group>0) {
                // if no artificial particles, continue
                if (hard_dump.ptcl_arti_bk.getPointer()!=NULL) {
                    // set all tidal tensor force to zero
                    for (int i=0; i<hard_dump.n_group; i++) {
                        int offset= i*hard_manager.ap_manager.getArtificialParticleN();
                        auto* pi = &(hard_dump.ptcl_arti_bk[offset]);
                        //auto* pcm = ap_manager.getCMParticles(pi);
                        auto* ptt = hard_manager.ap_manager.getTidalTensorParticles(pi);
                        for (int j=0; j<hard_manager.ap_manager.getTidalTensorParticleN(); j++) {
                            ptt[j].acc = PS::F64vec(0.0);
                        }
                    }
                }
            }
        }
#endif

        // running mode
        if (debug_io.mode.value==0) {
            //SystemHard sys;
            //sys.manager = &hard_manager;
            //sys.allocateHardIntegrator();

            // change ARC parameters
            //sys.driveForMultiClusterImpl(hard_dump.ptcl_bk.getPointer(), hard_dump.n_ptcl, hard_dump.ptcl_arti_bk.getPointer(), hard_dump.n_group, hard_dump.time_end, 0);

#ifdef EXTERNAL_HARD
#ifdef GALPY
            GalpyManager galpy_manager;
            std::string galpy_conf_filename = filename+".galpy";
            galpy_manager.initial(galpy_io, stat.time, galpy_conf_filename, true, true, false);
            hard_manager.h4_manager.interaction.ext_force.initial(ext_hard_io, galpy_manager, stat, hard_manager.center, hard_manager.center_id, true);
#else
            hard_manager.h4_manager.interaction.ext_force.initial(ext_hard_io, stat.time, hard_manager.center, hard_manager.center_id, true);
#endif
#endif
#ifdef STELLAR_EVOLUTION
            hard_manager.ar_manager.interaction.time_interrupt_max = hard_dump.time_end; // set max interrupt time to integration end time, to avoid unexpected interruption after integration end time due to inaccurate time offset in hard dump
#endif
            HardIntegrator hard_int;
#ifdef HARD_DEBUG_PRINT
            hard_int.output_filename_prefix = filename;
#endif
            auto* ptcl_artificial_ptr =  hard_dump.ptcl_arti_bk.getPointer();
            if (hard_dump.n_arti == 0) ptcl_artificial_ptr = NULL; // if no artificial particle, avoid reading artificial data from last hard_dump
            hard_int.initial(hard_dump.ptcl_bk.getPointer(), hard_dump.n_ptcl, ptcl_artificial_ptr, hard_dump.n_group, hard_dump.n_member_in_group.getPointer(), &hard_manager, hard_dump.time_offset, hard_dump.time_end);

            hard_int.integrateToTime(hard_dump.time_end);
            hard_int.driftClusterAndArtificialCMAndWriteBack(hard_dump.time_end, ptcl_artificial_ptr, hard_dump.n_group);

        }
        // test stability
        else if (debug_io.mode.value==1) {
            typedef H4::ParticleH4<PtclHard> PtclH4;

            SearchGroupCandidate<PtclH4> group_candidate;
            auto* ptcl = hard_dump.ptcl_bk.getPointer();
            PS::S32 n_ptcl = hard_dump.n_ptcl;

            group_candidate.searchAndMerge(ptcl, n_ptcl);

            PS::ReallocatableArray<PtclH4> ptcl_new;
            PS::S32 n_group_in_cluster;

            SystemHard sys;
            sys.manager = &hard_manager;

            PS::ReallocatableArray<COMM::BinaryTree<PtclH4,COMM::Binary>> binary_table;
            PS::ReallocatableArray<SystemHard::GroupIndexInfo> n_member_in_group;
            PS::ReallocatableArray<PS::S32> i_cluster_changeover_update;
            // generate artificial particles, stability test is included
            sys.findGroupsAndCreateArtificialParticlesOneCluster(0, ptcl, n_ptcl, ptcl_new, binary_table, n_group_in_cluster, n_member_in_group, i_cluster_changeover_update, group_candidate, hard_dump.time_end);
        }
    } 

    fclose(fp);

#ifdef STELLAR_EVOLUTION
    auto& interaction = hard_manager.ar_manager.interaction;
#ifdef BSE_BASE
    if (interaction.fout_sse.is_open()) interaction.fout_sse.close();
    if (interaction.fout_bse.is_open()) interaction.fout_bse.close();
#else
    if (interaction.fout_interrupt.is_open()) interaction.fout_interrupt.close();
#endif
#endif
    return 0;
}
