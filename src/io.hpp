#pragma once
#include "Common/io.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdio>
#include <cstring>

#define PRINT_WIDTH 15
#define PRINT_PRECISION 7

// stringify a macro's token text as-is
#define STRINGIFY(x) #x
// expand macros in x first, then stringify the result (e.g. XSTRINGIFY(G_ASTRO))
#define XSTRINGIFY(x) STRINGIFY(x)

// Bring COMM IO types into global namespace (PeTar code uses them without prefix)
using COMM::IOParams;
using COMM::IOParamsContainer;
using COMM::IOParamsPrintHelp;
using COMM::FindUndefinedOptions;

//! configure-time interrupt mode (exact token recovered from compile macros;
//! 'off' collapses to 'none' since no stellar-evolution columns exist then)
inline const char* getInterruptModeName() {
#ifdef BSEEMP
    return "bseEmp";
#elif defined(MOBSE)
    return "mobse";
#elif defined(SEVN)
    return "sevn";
#elif defined(BSE_BASE)
    return "bse";
#elif defined(DISK_STAR_MERGER)
    return "dsm";
#elif defined(STELLAR_EVOLUTION)
    return "merger";
#else
    return "none";
#endif
}

//! configure-time external potential mode
inline const char* getExternalModeName() {
#ifdef GALPY
    return "galpy";
#elif defined(AGAMA)
    return "agama";
#else
    return "none";
#endif
}

inline void writeOutputMeta(const std::string& _prefix, const int _data_format, const std::string& _version, const int _write_style = 1) {
    const char* interrupt = getInterruptModeName();
    const char* external = getExternalModeName();
    std::string family = std::string("petar");
    if (std::string(interrupt) != "none") family += std::string(".") + interrupt;
    if (std::string(external) != "none") family += std::string(".") + external;
    FILE* fmeta = fopen((_prefix + ".meta").c_str(), "w");
    if (fmeta == NULL) {
        std::cerr<<"Error: cannot open "<<_prefix<<".meta for writing.\n";
        return;
    }
    fprintf(fmeta, "# output-meta/1: output schema tag; see doc/output-schema-tagging-plan.md\n");
    fprintf(fmeta, "format = output-meta/1\n");
    fprintf(fmeta, "producer = petar\n");
    fprintf(fmeta, "petar_version = %s\n", _version.c_str());
    fprintf(fmeta, "binary_family = %s\n", family.c_str());
    fprintf(fmeta, "snapshot_format = %s\n", (_data_format == 1 || _data_format == 3) ? "ascii" : "binary");
    fprintf(fmeta, "snapshot_container = %s\n", (_write_style == 4) ? "pack" : "file");
    fprintf(fmeta, "layout.interrupt_mode = %s\n", interrupt);
    fprintf(fmeta, "layout.external_mode = %s\n", external);
    fprintf(fmeta, "layout.use_mpfrc = %d\n",
#ifdef PETAR_USE_MPFRC
            1
#else
            0
#endif
            );
    fprintf(fmeta, "layout.collect_sp_acc = %d\n",
#ifdef COLLECT_SP_ACC
            1
#else
            0
#endif
            );
#ifdef BSE_BASE
    // BSE-family star.spin width: 3 columns since 2024-12 (ospin[3] in bse_interface.h)
    fprintf(fmeta, "layout.spin_3d = true\n");
#endif
    fprintf(fmeta, "layout.float_type = %s\n",
#ifdef PETAR_USE_MPFRC
            "float128"
#else
            "float64"
#endif
            );
    fclose(fmeta);
}


class FileHeader{
public:
    long long int nfile;  // file id
    long long int n_body;
    double time;
#ifdef RECORD_CM_IN_HEADER
    PS::F64vec pos_offset;
    PS::F64vec vel_offset;
#endif
    FileHeader(){
        n_body = 0;
        time = 0.0;
#ifdef RECORD_CM_IN_HEADER
        pos_offset = PS::F64vec(0.0);
        vel_offset = PS::F64vec(0.0);
#endif
    }
#ifdef RECORD_CM_IN_HEADER
    FileHeader(const long long int ni, const long long int n, const double t, const PS::F64vec& pos, const PS::F64vec& vel){
        nfile = ni;
        n_body = n;
        time = t;
        pos_offset = pos;
        vel_offset = vel;
    }
    int readAscii(FILE * fp){
        int rcount=fscanf(fp, "%lld %lld %lf %lf %lf %lf %lf %lf %lf\n", &nfile, &n_body, &time, &pos_offset.x,&pos_offset.y, &pos_offset.z, &vel_offset.x, &vel_offset.y, &vel_offset.z);
        if (rcount<9) {
          std::cerr<<"Error: cannot read header, please check your data file header!\n";
          abort();
        }
        //std::cout<<"Number of particles ="<<n_body<<";  Time="<<time<<std::endl;
        return n_body;
    }

    void writeAscii(FILE* fp) const{
        fprintf(fp, "%lld %lld %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e\n", nfile, n_body, time, pos_offset.x, pos_offset.y, pos_offset.z, vel_offset.x, vel_offset.y, vel_offset.z);
    }
#else
    FileHeader(const long long int ni, const long long int n, const double t) {
        nfile = ni;
        n_body = n;
        time = t;
    }

    int readAscii(FILE * fp){
        int rcount=fscanf(fp, "%lld %lld %lf\n", &nfile, &n_body, &time);
        if (rcount<3) {
          std::cerr<<"Error: cannot read header, please check your data file header!\n";
          abort();
        }
        //std::cout<<"Number of particles ="<<n_body<<";  Time="<<time<<std::endl;
        return n_body;
    }

    void writeAscii(FILE* fp) const{
        fprintf(fp, "%lld %lld %26.17e\n", nfile, n_body, time);
    }
#endif
    int readBinary(FILE* fp){
        size_t rcount=fread(this, sizeof(FileHeader), 1, fp);
        if(rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is "<<1<<" bytes, only obtain "<<rcount<<" bytes.\n";
            abort();
        }
        //std::cout<<"Number of particles ="<<n_body<<";  Time="<<time<<std::endl;
        return n_body;
    }

    void writeBinary(FILE* fp) const{
        fwrite(this, sizeof(FileHeader), 1, fp);
    }

    int readBinaryStream(std::istream& _fin){
        _fin.read(reinterpret_cast<char*>(this), sizeof(FileHeader));
        if(!_fin) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1 bytes.\n";
            abort();
        }
        return n_body;
    }

    void printColumnBinary(std::ostream& _fout) const{
        _fout.write(reinterpret_cast<const char*>(this), sizeof(FileHeader));
    }

};