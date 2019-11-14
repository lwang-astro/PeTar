#pragma once

#define PRINT_WIDTH 18
#define PRINT_PRECISION 14

//! IO Params container
class IOParamsContainer{
    PS::ReallocatableArray<PS::F64*> d_f64;
    PS::ReallocatableArray<PS::S64*> d_s64;
    PS::ReallocatableArray<PS::S32*> d_s32;
    PS::ReallocatableArray<std::string*> d_str;
    
public:
    
    void store(PS::F64* _item) {
        d_f64.push_back(_item);
    }

    void store(PS::S64* _item) {
        d_s64.push_back(_item);
    }

    void store(PS::S32* _item) {
        d_s32.push_back(_item);
    }

    void store(std::string* _item) {
        d_str.push_back(_item);
    }

    void writeAscii(FILE *_fout) {
        for(PS::S32 i=0; i<d_f64.size(); i++) fprintf(_fout, "%26.15e ", *d_f64[i]);
        for(PS::S32 i=0; i<d_s64.size(); i++) fprintf(_fout, "%lld ",    *d_s64[i]);
        for(PS::S32 i=0; i<d_s32.size(); i++) fprintf(_fout, "%d ",      *d_s32[i]);
        for(PS::S32 i=0; i<d_str.size(); i++) fprintf(_fout, "%s ",d_str[i]->c_str());
        fprintf(_fout,"\n");
    }
    
    void readAscii(FILE *_fin) {
        size_t rcount=0;
        for(PS::S32 i=0; i<d_f64.size(); i++) {
            rcount=fscanf(_fin, "%lf ", d_f64[i]);
            if (rcount<1) {
                std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                abort();
            }
        }
        for(PS::S32 i=0; i<d_s64.size(); i++) {
            rcount=fscanf(_fin, "%lld ", d_s64[i]);
            if (rcount<1) {
                std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                abort();
            }
        }
        for(PS::S32 i=0; i<d_s32.size(); i++) {
            rcount=fscanf(_fin, "%d ", d_s32[i]);
            if (rcount<1) {
                std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                abort();
            }
        }
        for(PS::S32 i=0; i<d_str.size(); i++) {
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
        for(PS::S32 i=0; i<d_f64.size(); i++) PS::Comm::broadcast(d_f64[i], 1, 0);
        for(PS::S32 i=0; i<d_s64.size(); i++) PS::Comm::broadcast(d_s64[i], 1, 0);
        for(PS::S32 i=0; i<d_s32.size(); i++) PS::Comm::broadcast(d_s32[i], 1, 0);
        for(PS::S32 i=0; i<d_str.size(); i++) {
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
    PS::S64 nfile;  // file id
    PS::S64 n_body;
    //PS::S64 id_offset; // file id offset for add new artificial particles, should be larger than the maximum file id
    PS::F64 time;
    //PS::F64 dt_soft;   // tree time step should be recorded for restarting (soft kick of cm)
    //PS::S64 n_split;   // n_split is also needed for restarting (soft kick of cm)
    FileHeader(){
        n_body = 0;
        time = 0.0;
    }
    FileHeader(const PS::S64 ni, const PS::S64 n, const PS::F64 t){
        nfile = ni;
        n_body = n;
        time = t;
    }
    PS::S32 readAscii(FILE * fp){
        PS::S32 rcount=fscanf(fp, "%lld %lld %lf\n", &nfile, &n_body, &time);
        if (rcount<3) {
          std::cerr<<"Error: cannot read header, please check your data file header!\n";
          abort();
        }
        std::cout<<"Number of particles ="<<n_body<<";  Time="<<time<<std::endl;
        return n_body;
    }

    PS::S32 readBinary(FILE* fp){
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
