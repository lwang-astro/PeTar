#include <iostream>
#include <iomanip>
#include <cstdio>
#include <string>
#include <cstring>
#include <vector>
#include <getopt.h>
#include <algorithm>
#include <cmath>
#include <cinttypes>
#include <ctime>
#include <cassert>

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

// -----------------------------------------------------------------------
//  IOParamsDump2Test: CLI parameters for the conversion tool
// -----------------------------------------------------------------------
class IOParamsDump2Test {
public:
    IOParamsContainer input_par_store;
    IOParams<PS::S64> n_crit_ptcl;
    IOParams<PS::F64> tstart;
    IOParams<PS::F64> tend;
    IOParams<PS::S64> istart;
    IOParams<PS::S64> iend;
    IOParams<PS::S64> n_crit_group;
    IOParams<PS::S64> n_crit_arti;
    IOParams<std::string> fname_par;
    IOParams<std::string> fname_dump;
    bool print_flag;

    IOParamsDump2Test(): input_par_store(),
                         n_crit_ptcl(input_par_store, 0, "n", "if >0, only convert dumps where particle number matches"),
                         tstart(input_par_store, -1.0, "tstart", "if >0, only convert dumps where physical time >= tstart"),
                         tend(input_par_store, -1.0, "tend", "if >0, only convert dumps where physical time < tend"),
                         istart(input_par_store, -1, "istart", "if >0, only convert dumps with index >= istart (1-based)"),
                         iend(input_par_store, -1, "iend", "if >0, only convert dumps with index < iend (1-based)"),
                         n_crit_group(input_par_store, 0, "n-crit-group", "if >0, only convert dumps where group number matches"),
                         n_crit_arti(input_par_store, 0, "n-crit-arti", "if >0, only convert dumps where artificial particle number matches"),
                         fname_par(input_par_store, "data.par", "p", "Parameter file prefix; reads <prefix>.hard"),
                         fname_dump(input_par_store, "hard_dump", "dump-filename", "Hard dump data filename"),
                         print_flag(false)
    {}

    int read(int argc, char* argv[], const int opt_used_pre=0) {
        static int flag=-1;
        static struct option long_options[] = {
            {tstart.key, required_argument, &flag, 0},
            {tend.key,   required_argument, &flag, 1},
            {istart.key, required_argument, &flag, 2},
            {iend.key,   required_argument, &flag, 3},
            {n_crit_group.key, required_argument, &flag, 4},
            {n_crit_arti.key,  required_argument, &flag, 5},
            {"help", no_argument, 0, 'h'},
            {0,0,0,0}
        };

        int opt_used = opt_used_pre;
        int copt;
        int option_index;
        optind = 0;
        while ((copt = getopt_long(argc, argv,
                                   "-n:p:h",
                                   long_options, &option_index)) != -1) {
            switch (copt) {
            case 0:
                switch (flag) {
                case 0: tstart.value = atof(optarg); if(print_flag) tstart.print(std::cout); opt_used += 2; break;
                case 1: tend.value   = atof(optarg); if(print_flag) tend.print(std::cout);   opt_used += 2; break;
                case 2: istart.value = atoi(optarg); if(print_flag) istart.print(std::cout); opt_used += 2; break;
                case 3: iend.value   = atoi(optarg); if(print_flag) iend.print(std::cout);   opt_used += 2; break;
                case 4: n_crit_group.value = atoi(optarg); if(print_flag) n_crit_group.print(std::cout); opt_used += 2; break;
                case 5: n_crit_arti.value  = atoi(optarg); if(print_flag) n_crit_arti.print(std::cout);  opt_used += 2; break;
                default: break;
                }
                break;
            case 'n':
                n_crit_ptcl.value = atoi(optarg);
                if(print_flag) n_crit_ptcl.print(std::cout);
                opt_used += 2;
                break;
            case 'p':
                fname_par.value = optarg;
                if(print_flag) std::cout<<"Parameter prefix: "<<fname_par.value<<std::endl;
                opt_used += 2;
                break;
            case 'h':
                if(print_flag) {
                    std::cout<<"Convert PeTar hard dump files to petar.hard.test-compatible snapshots\n"
                             <<"Usage: petar.dump2test [options] [dump filename]\n"
                             <<"   Options:\n"
                             <<"    -p <prefix>   Parameter file prefix (reads <prefix>.hard)\n"
                             <<"    -n <N>        Only convert dumps with exactly N particles\n"
                             <<"    --tstart <T>  Only convert dumps with time >= T\n"
                             <<"    --tend <T>    Only convert dumps with time < T\n"
                             <<"    --istart <I>  Only convert dumps with index >= I (1-based)\n"
                             <<"    --iend <I>    Only convert dumps with index < I (1-based)\n"
                             <<"    --n-crit-group <N> Only convert dumps with exactly N groups\n"
                             <<"    --n-crit-arti <N>  Only convert dumps with exactly N artificial particles\n"
                             <<"    -h (--help)   Print this help\n";
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
            if(print_flag) std::cout<<"Dump filename: "<<fname_dump.value<<std::endl;
        }
        return opt_used-1;
    }

    bool checkParams() {
        assert(n_crit_ptcl.value >= 0 && n_crit_group.value >= 0 && n_crit_arti.value >= 0);
        assert(istart.value==-1 || istart.value>=1);
        assert(iend.value==-1 || iend.value>=1);
        assert(!(istart.value>0 && iend.value>0 && iend.value<istart.value));
        assert(!(tstart.value>0 && tend.value>0 && tend.value<=tstart.value));
        return true;
    }
};

// -----------------------------------------------------------------------
//  Helper: check filename suffix
// -----------------------------------------------------------------------
static bool has_suffix(const std::string& str, const std::string& suffix) {
    return str.size() >= suffix.size() &&
           str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

// -----------------------------------------------------------------------
//  Infer parameter prefix from dump filename
//  e.g. "data.hard_dump_t12345_M0_O0_c1" -> prefix "data.par"
// -----------------------------------------------------------------------
static bool infer_par_from_dump(const std::string& dump_name,
                                std::string& par_prefix) {
    std::size_t dot_pos = dump_name.find('.');
    if (dot_pos != std::string::npos && dot_pos > 0) {
        par_prefix = dump_name.substr(0, dot_pos) + ".par";
        return true;
    }
    return false;
}

// -----------------------------------------------------------------------
//  Write one cluster's particles as an ASCII snapshot for petar.hard.test
//
//  Format (matching FPSoft::readAscii without GROUP_DATA_WRITE_ARTIFICIAL,
//  without COLLECT_SP_ACC, without EXTERNAL_POT_IN_PTCL):
//
//    FileHeader line:   nfile n_body time
//    Per particle:      mass pos.x pos.y pos.z vel.x vel.y vel.z binary_state
//                       r_search id group_data1 group_data2
//                       changeover.r_in changeover.r_out
//                       acc.x acc.y acc.z pot_tot pot_soft n_ngb
// -----------------------------------------------------------------------
static void writeSnapshotAscii(const HardDump& dump,
                               const std::string& snap_filename,
                               int cluster_index) {

    FILE* fp = fopen(snap_filename.c_str(), "w");
    if (!fp) {
        std::cerr << "Error: Cannot open output file " << snap_filename << "\n";
        abort();
    }

    // -- FileHeader: nfile=0, n_body, time --
    fprintf(fp, "0 %d %26.17e\n", dump.n_ptcl, dump.time_offset);

    // -- Write each particle --
    for (int i = 0; i < dump.n_ptcl; i++) {
        const auto& src = dump.ptcl_bk[i];

        // 1) Determine the real mass
        PS::F64 real_mass;
        if (src.group_data.artificial.isSingle()) {
            real_mass = src.mass;
        } else if (src.group_data.artificial.isMember()) {
            real_mass = src.group_data.artificial.getMassBackup();
        } else {
            // CM or other type — unlikely in dump for real particles
            real_mass = src.mass;
        }

        // 2) Determine changeover — apply pending r_scale_next
        PS::F64 r_in_val  = src.changeover.getRin();
        PS::F64 r_out_val = src.changeover.getRout();
        if (src.changeover.r_scale_next != 1.0) {
            r_in_val  *= src.changeover.r_scale_next;
            r_out_val *= src.changeover.r_scale_next;
        }

        // 3) Write ParticleBase fields
        fprintf(fp, "%26.17e %26.17e %26.17e %26.17e "
                    "%26.17e %26.17e %26.17e 0 ",
                real_mass,
                src.pos.x, src.pos.y, src.pos.z,
                src.vel.x, src.vel.y, src.vel.z);

        // 4) Write Ptcl fields (r_search, id, group_data=0,0)
        fprintf(fp, "%26.17e %" PRId64 " 0 0 ",
                src.r_search, src.id);

        // 5) Write changeover
        fprintf(fp, "%26.17e %26.17e ",
                r_in_val, r_out_val);

        // 6) Write FPSoft fields (acc=0, pot=0, pot_soft=0, n_ngb=0)
        fprintf(fp, "0 0 0 0 0 0\n");
    }

    fclose(fp);
    std::cerr << "  Wrote snapshot: " << snap_filename
              << "  (" << dump.n_ptcl << " particles, "
              << dump.n_group << " groups)" << std::endl;
}

// -----------------------------------------------------------------------
//  Compute base r_out and r_in/r_out ratio from particle changeovers
//
//  In the HardManager::initial flow:
//    changeover.setR(mass * mean_mass_inv, r_in_base, r_out_base)
//  which does:
//    m_fac3 = max(pow(mass * mean_mass_inv, 1/3), 1.0)
//    r_in  = m_fac3 * r_in_base
//    r_out = m_fac3 * r_out_base
//
//  For particles with mass <= mean_mass, m_fac3 = 1, so
//    r_in_min = r_in_base,  r_out_min = r_out_base
//
//  However, the dump particles may already have mass=0 for members
//  (backup mass was used to compute the changeover originally).
//  We use the RESTORED mass for this calculation.
// -----------------------------------------------------------------------
struct BaseRadii {
    PS::F64 r_out_base;
    PS::F64 r_in_over_out_base;
};

static BaseRadii computeBaseRadii(const HardDump& dump) {
    PS::F64 r_in_min  = PS::LARGE_FLOAT;
    PS::F64 r_out_min = PS::LARGE_FLOAT;

    for (int i = 0; i < dump.n_ptcl; i++) {
        const auto& p = dump.ptcl_bk[i];
        PS::F64 rin  = p.changeover.getRin();
        PS::F64 rout = p.changeover.getRout();
        if (rout > 0 && rout < r_out_min) r_out_min = rout;
        if (rin  > 0 && rin  < r_in_min)  r_in_min  = rin;
    }

    // Fallback if no valid radii found (shouldn't happen with valid data)
    if (r_out_min >= PS::LARGE_FLOAT * 0.5) r_out_min = 1.0;
    if (r_in_min  >= PS::LARGE_FLOAT * 0.5) r_in_min  = r_out_min * 0.1;

    BaseRadii br;
    br.r_out_base          = r_out_min;
    br.r_in_over_out_base  = r_in_min / r_out_min;
    return br;
}

// -----------------------------------------------------------------------
//  Print the recommended petar.hard.test command
// -----------------------------------------------------------------------
static void printCommand(const std::string& dump_basename,
                         const std::string& par_prefix,
                         const std::string& snap_filename,
                         int cluster_index,
                         const HardDump& dump,
                         const BaseRadii& br,
                         const IOParamsHard& hard_io,
                         bool r_search_min_valid) {

    std::cout << "\n"
              << "# ============================================================\n"
              << "# Cluster " << cluster_index << " converted from: " << dump_basename << "\n"
              << "#   n_ptcl=" << dump.n_ptcl
              << "  n_group=" << dump.n_group
              << "  n_arti=" << dump.n_arti << "\n"
              << "#   time_offset=" << dump.time_offset
              << "  time_end=" << dump.time_end << "\n"
              << "#   par prefix: " << par_prefix << "\n"
              << "# ============================================================\n"
              << "#\n"
              << "# Suggested petar.hard.test command:\n"
              << "#\n";

    // output filename prefix = snap filename without .snap suffix
    std::string out_prefix = snap_filename;
    if (out_prefix.size() >= 5 && out_prefix.substr(out_prefix.size()-5) == ".snap")
        out_prefix.resize(out_prefix.size() - 5);

    std::cout << "petar.hard.test -p " << par_prefix
              << " -f " << out_prefix
              << " -t " << std::setprecision(17) << dump.time_end
              << " -r " << br.r_out_base
              << " --r-ratio " << br.r_in_over_out_base
              << " -s " << hard_io.dt_max_hermite.value;

    if (r_search_min_valid)
        std::cout << " --r-search-min " << Ptcl::r_search_min;

    std::cout << " " << snap_filename << "\n"
              << "#\n"
              << std::endl;
}

// =======================================================================
//  Main
// =======================================================================
int main(int argc, char** argv) {

    // ---- Parse CLI arguments ----
    IOParamsDump2Test io;
    IOParamsHard hard_io;

    // Auto-inject -p unless explicitly provided (same pattern as hard_debug.cxx)
    opterr = 0;
    std::vector<std::string> adjusted_args;
    adjusted_args.reserve(argc + 2);
    bool has_p_option = false;
    bool has_help_option = false;
    adjusted_args.push_back(argv[0]);
    for (int i = 1; i < argc; i++) adjusted_args.push_back(argv[i]);
    for (int i = 1; i < argc; i++) {
        if (adjusted_args[i] == "-p" && i+1 < argc) {
            has_p_option = true;
            if (has_suffix(adjusted_args[i+1], ".hard")) {
                adjusted_args[i+1] = adjusted_args[i+1].substr(0, adjusted_args[i+1].size() - 5);
            }
        }
        if (adjusted_args[i] == "-h" || adjusted_args[i] == "--help")
            has_help_option = true;
    }
    if (!has_help_option && !has_p_option) {
        std::string par_prefix = io.fname_par.value;
        if (adjusted_args.size() > 1) {
            const std::string& dump_candidate = adjusted_args.back();
            if (!dump_candidate.empty() && dump_candidate[0] != '-') {
                std::string inferred;
                if (infer_par_from_dump(dump_candidate, inferred))
                    par_prefix = inferred;
            }
        }
        adjusted_args.insert(adjusted_args.begin()+1, "-p");
        adjusted_args.insert(adjusted_args.begin()+2, par_prefix);
    }

    std::vector<std::vector<char>> adjusted_argbuf;
    adjusted_argbuf.reserve(adjusted_args.size());
    std::vector<char*> adjusted_argv;
    adjusted_argv.reserve(adjusted_args.size());
    for (const auto& item : adjusted_args) {
        adjusted_argbuf.emplace_back(item.begin(), item.end());
        adjusted_argbuf.back().push_back('\0');
        adjusted_argv.push_back(adjusted_argbuf.back().data());
    }

    io.print_flag = true;
    const bool help_flag = (io.read((int)adjusted_argv.size(), adjusted_argv.data()) == -1);
    if (!help_flag) io.checkParams();

    hard_io.print_flag = true;
    hard_io.read((int)adjusted_argv.size(), adjusted_argv.data(), false);

    if (help_flag) return 0;

    const std::string filename       = io.fname_dump.value;
    const std::string fpval = io.fname_par.value;
    const std::string fhardpar_ascii = has_suffix(fpval, ".hard")
        ? fpval.substr(0, fpval.size() - 5)
        : fpval;
    const std::string fhardpar_report = fhardpar_ascii + ".hard";

    std::cerr << "Dump file:       " << filename << std::endl;
    std::cerr << "Hard par file:   " << fhardpar_report << std::endl;

    // ---- Read .hard parameter file ----
    // IOParamsHard::read with -p flag in argv already loaded the
    // .hard file's contents into hard_io.input_par_store.
    // Now verify that key parameters are present / usable.

    if (hard_io.dt_max_hermite.value <= 0.0) {
        std::cerr << "Error: hermite-dt-max <= 0. Cannot infer dt_soft for hard.test.\n"
                  << "Please ensure the .hard file contains a valid hermite-dt-max.\n";
        abort();
    }

    // ---- Open dump file ----
    std::FILE* fp = std::fopen(filename.c_str(), "rb");
    if (!fp) {
        std::cerr << "Error: Cannot open dump file " << filename << "\n";
        abort();
    }

    // ---- Iterate over clusters ----
    int ncount = 0;
    int n_converted = 0;

    while (true) {
        int c = fgetc(fp);
        if (c == EOF) break;
        ungetc(c, fp);

        HardDump hard_dump;
        hard_dump.readOneClusterBinary(fp);

#ifdef EXTERNAL_HARD
        // Skip EXTERNAL_HARD center data that follows each cluster in dump
        // (just advance the file pointer — we don't need the data)
        FPSoft dummy_center;
        dummy_center.readBinary(fp);
#endif

        ncount++;

        // ---- Filter criteria ----
        if (io.n_crit_ptcl.value > 0 && hard_dump.n_ptcl != io.n_crit_ptcl.value) continue;
        if (io.n_crit_group.value > 0 && hard_dump.n_group != io.n_crit_group.value) continue;
        if (io.n_crit_arti.value  > 0 && hard_dump.n_arti  != io.n_crit_arti.value) continue;
        if (io.tstart.value > 0 && hard_dump.time_offset < io.tstart.value) continue;
        if (io.tend.value   > 0 && hard_dump.time_offset >= io.tend.value) continue;
        if (ncount < io.istart.value) continue;
        if (io.iend.value > 0 && ncount >= io.iend.value) continue;

        // ---- Restore static particle parameters from dump ----
        // These were set by HardDump::readOneClusterBinary
        // Ptcl::search_factor, Ptcl::r_search_min, Ptcl::mean_mass_inv,
        // PtclHard::r_group_over_in, PtclHard::r_search_group_over_in

        // ---- Compute base r_out, r_in/r_out ratio from particle changeovers ----
        BaseRadii br = computeBaseRadii(hard_dump);

        // ---- Generate snapshot filename ----
        std::string snap_filename = filename + ".cluster" + std::to_string(ncount) + ".snap";

        // ---- Write ASCII snapshot ----
        writeSnapshotAscii(hard_dump, snap_filename, ncount);

        // ---- Print recommended command ----
        bool r_search_min_valid = (Ptcl::r_search_min > 0.0);
        printCommand(filename, fhardpar_ascii, snap_filename, ncount, hard_dump, br,
                     hard_io, r_search_min_valid);

        n_converted++;
    }

    fclose(fp);

    std::cerr << "\n"
              << "=== Conversion complete ===\n"
              << "  Total clusters in dump: " << ncount << "\n"
              << "  Converted:              " << n_converted << "\n"
              << "  Run the printed petar.hard.test command(s) above.\n"
              << std::endl;

    return 0;
}
