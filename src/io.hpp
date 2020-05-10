#pragma once
#include <vector>

#define PRINT_WIDTH 18
#define PRINT_PRECISION 14

//! IO Params container
class IOParamsContainer{
    std::vector<double*> d_f64;
    std::vector<long long int*> d_s64;
    std::vector<int*> d_s32;
    std::vector<std::string*> d_str;
    
public:
    
    void store(double* _item) {
        d_f64.push_back(_item);
    }

    void store(long long int* _item) {
        d_s64.push_back(_item);
    }

    void store(int* _item) {
        d_s32.push_back(_item);
    }

    void store(std::string* _item) {
        d_str.push_back(_item);
    }

    void writeAscii(FILE *_fout) {
        for(size_t i=0; i<d_f64.size(); i++) fprintf(_fout, "%26.15e ", *d_f64[i]);
        for(size_t i=0; i<d_s64.size(); i++) fprintf(_fout, "%lld ",    *d_s64[i]);
        for(size_t i=0; i<d_s32.size(); i++) fprintf(_fout, "%d ",      *d_s32[i]);
        for(size_t i=0; i<d_str.size(); i++) fprintf(_fout, "%s ",d_str[i]->c_str());
        fprintf(_fout,"\n");
    }
    
    void readAscii(FILE *_fin) {
        size_t rcount=0;
        for(size_t i=0; i<d_f64.size(); i++) {
            rcount=fscanf(_fin, "%lf ", d_f64[i]);
            if (rcount<1) {
                std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                abort();
            }
        }
        for(size_t i=0; i<d_s64.size(); i++) {
            rcount=fscanf(_fin, "%lld ", d_s64[i]);
            if (rcount<1) {
                std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                abort();
            }
        }
        for(size_t i=0; i<d_s32.size(); i++) {
            rcount=fscanf(_fin, "%d ", d_s32[i]);
            if (rcount<1) {
                std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                abort();
            }
        }
        for(size_t i=0; i<d_str.size(); i++) {
            char dtmp[1024];
            rcount=fscanf(_fin, "%s ", dtmp);
            if (rcount<1) {
                std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                abort();
            }
            *d_str[i] = dtmp;
        }
    }

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
    void mpi_broadcast() {
        for(size_t i=0; i<d_f64.size(); i++) PS::Comm::broadcast(d_f64[i], 1, 0);
        for(size_t i=0; i<d_s64.size(); i++) PS::Comm::broadcast(d_s64[i], 1, 0);
        for(size_t i=0; i<d_s32.size(); i++) PS::Comm::broadcast(d_s32[i], 1, 0);
        for(size_t i=0; i<d_str.size(); i++) {
            size_t str_size=d_str[i]->size();
            PS::Comm::broadcast(&str_size, 1, 0);
            char stmp[str_size+1];
            std::strcpy(stmp, d_str[i]->c_str());
            PS::Comm::broadcast(stmp, str_size, 0);
            *d_str[i]=std::string(stmp);
        }
    }
#endif
};

// IO Params
template <class Type>
struct IOParams{
    Type value;
    const char* name;
    const char* defaulted;

    IOParams(IOParamsContainer& _ioc, const Type& _value, const char* _name, const char* _defaulted=NULL): value(_value), name(_name), defaulted(_defaulted)  {
        _ioc.store(&value);
    }

    void print(std::ostream& os) const{
        os<<name<<":   "<<value<<std::endl;
    }
};

template <class Type>
std::ostream& operator <<(std::ostream& os, const IOParams<Type>& par) {
    if (par.defaulted!=NULL) os<<par.name<<": "<<par.defaulted;
    else os<<par.name<<": "<<par.value;
    return os;
}


class FileHeader{
public:
    long long int nfile;  // file id
    long long int n_body;
    double time;
    FileHeader(){
        n_body = 0;
        time = 0.0;
    }
    FileHeader(const long long int ni, const long long int n, const double t){
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
        std::cout<<"Number of particles ="<<n_body<<";  Time="<<time<<std::endl;
        return n_body;
    }

    int readBinary(FILE* fp){
        size_t rcount=fread(this, sizeof(FileHeader), 1, fp);
        if(rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is "<<1<<" bytes, only obtain "<<rcount<<" bytes.\n";
            abort();
        }
        std::cout<<"Number of particles ="<<n_body<<";  Time="<<time<<std::endl;
        return n_body;
    }

    void writeAscii(FILE* fp) const{
        fprintf(fp, "%lld %lld %26.17e\n", nfile, n_body, time);
    }

    void writeBinary(FILE* fp) const{
        fwrite(this, sizeof(FileHeader), 1, fp);
    }
};
