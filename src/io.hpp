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

// Bring COMM IO types into global namespace (PeTar code uses them without prefix)
using COMM::IOParams;
using COMM::IOParamsContainer;
using COMM::IOParamsPrintHelp;
using COMM::FindUndefinedOptions;


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