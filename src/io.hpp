#pragma once
#include <iostream>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>

#define PRINT_WIDTH 15
#define PRINT_PRECISION 7

// print format parameters
struct IOParamsPrintHelp{
    int offset_short_key;
    int offset_long_key;
    int width_key;

    IOParamsPrintHelp(const int _offset_short_key, const int _offset_long_key, const int _width_key): 
        offset_short_key(_offset_short_key), offset_long_key(_offset_long_key), width_key(_width_key) {}

    static void printTypeShortNameDescription(std::ostream& os) {
        os<<"I: 64bit integer; F: 64bit floating; S: string\n";
    }

    static char getValueTypeShortName(const long long int& value) {
        return 'I';
    }

    static char getValueTypeShortName(const double& value) {
        return 'F';
    }

    static char getValueTypeShortName(const std::string& value) {
        return 'S';
    }
};


// IO Params
template <class Type>
struct IOParams{
    Type value;
    const char* key;
    const char* description;
    const char* defaulted;
    bool print_help_flag;

    template <class TContainer>
    IOParams(TContainer& _ioc, const Type& _value, const char* _key, const char* _description, const char* _defaulted=NULL, const bool _print_help_flag=true): value(_value), key(_key), description(_description), defaulted(_defaulted), print_help_flag(_print_help_flag)  {
        _ioc.store(_key, this);
    }

    void print(std::ostream& os) const{
        os<<description<<":   "<<value<<std::endl;
    }
    
    void printHelp(std::ostream& os, const IOParamsPrintHelp& _align, const bool print_short_flag, const bool always_print=false) const {
        if (print_help_flag || always_print) {
            //key
            bool print_flag = false;
            int multiline_width = 0;

            if (strlen(key)==1 && print_short_flag) {
                os<<std::setw(_align.offset_short_key)<<"-"<<key<<"  ";
                //  <<std::setw(_align.width_key+_align.offset_long_key-_align.offset_short_key-1)<<" "; 
                multiline_width = _align.offset_short_key+5;
                print_flag = true;
            }
            else if (strlen(key)>1 && !print_short_flag) {
                os<<std::setw(_align.offset_long_key)<<"--"
                  <<std::left<<std::setw(_align.width_key)<<key<<std::right;   
                multiline_width = _align.offset_long_key+_align.width_key+5;
                print_flag = true;
            }

            if (print_flag) {
                // type
                os<<"["<<IOParamsPrintHelp::getValueTypeShortName(value)<<"] ";

                // description            
                std::vector<std::string> tokens;
                std::stringstream ss(description);
                std::string token;

                while (std::getline(ss, token, ';')) {
                    tokens.push_back(token);
                }

                if (tokens.size()>1) {
                    if (defaulted!=NULL) os<<tokens[0]<<": "<<defaulted<<std::endl;
                    else os<<tokens[0]<<": "<<value<<std::endl;
                    for (size_t i=1; i<tokens.size(); i++) 
                        os<<std::setw(multiline_width)<<" "<<tokens[i]<<std::endl;
                }
                else os<<*this<<std::endl;
            }
        }
    }
};

template <class Type>
std::ostream& operator <<(std::ostream& os, const IOParams<Type>& par) {
    if (par.defaulted!=NULL) os<<par.description<<": "<<par.defaulted;    
    else os<<par.description<<": "<<par.value;
    return os;
}

//! IO Params container
class IOParamsContainer{
    struct char_cmp {
        bool operator () (const char *a,const char *b) const {
            return strcmp(a,b)<0;
        }
    };
    std::map<const char*, IOParams<double>*, char_cmp> d_f64;
    std::map<const char*, IOParams<long long int>*, char_cmp> d_i64;
    std::map<const char*, IOParams<std::string>*, char_cmp> d_str;
    
public:

    void store(const char* _name, IOParams<double>* _item) {
        d_f64[_name] = _item;
    }

    void store(const char* _name, IOParams<long long int>* _item) {
        d_i64[_name] = _item;
    }

    void store(const char* _name, IOParams<std::string>* _item) {
        d_str[_name] = _item;
    }

    void writeAscii(FILE *_fout) {
        for(auto iter = d_f64.begin(); iter!=d_f64.end(); iter++) fprintf(_fout, "%c %s %26.15e\n", IOParamsPrintHelp::getValueTypeShortName(iter->second->value), iter->first, iter->second->value);
        for(auto iter = d_i64.begin(); iter!=d_i64.end(); iter++) fprintf(_fout, "%c %s %lld\n",    IOParamsPrintHelp::getValueTypeShortName(iter->second->value), iter->first, iter->second->value);
        for(auto iter = d_str.begin(); iter!=d_str.end(); iter++) fprintf(_fout, "%c %s %s\n",      IOParamsPrintHelp::getValueTypeShortName(iter->second->value), iter->first, iter->second->value.c_str());
    }
    
    void readAscii(FILE *_fin) {
        
        char key_name[1024];
        while (!feof(_fin)) {
            size_t rcount=0;
            char type_id;
            rcount=fscanf(_fin, "%c ", &type_id);
            if (rcount<1) {
                std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                abort();
            }

            switch (type_id) {
            case 'F':
            {
                double dtmp;
                rcount=fscanf(_fin, "%s %lf\n", key_name, &dtmp);
                if (rcount<2) {
                    std::cerr<<"Error: Data reading fails! requiring data number is 2, only obtain "<<rcount<<".\n";
                    abort();
                }
                auto search = d_f64.find(key_name);
                if (search == d_f64.end())
                    std::cerr<<"Warning: parameter name key "<<key_name<<" is not found!\n";
                else 
                    search->second->value = dtmp;
                break;
            }
            case 'I':
            {
                long long int dtmp;
                rcount=fscanf(_fin, "%s %lld\n", key_name, &dtmp);
                if (rcount<2) {
                    std::cerr<<"Error: Data reading fails! requiring data number is 2, only obtain "<<rcount<<".\n";
                    abort();
                }
                auto search = d_i64.find(key_name);
                if (search == d_i64.end()) 
                    std::cerr<<"Warning: parameter name key "<<key_name<<" is not found!\n";
                else 
                    search->second->value = dtmp;
                break;
            }
            case 'S':
            {
                char dtmp[1024];
                rcount=fscanf(_fin, "%s %s\n", key_name, dtmp);
                if (rcount<2) {
                    std::cerr<<"Error: Data reading fails! requiring data number is 2, only obtain "<<rcount<<".\n";
                    abort();
                }
                auto search = d_str.find(key_name);
                if (search == d_str.end()) 
                    std::cerr<<"Warning: parameter name key "<<key_name<<" is not found!\n";
                else 
                    search->second->value = dtmp;
                break;
            }
            default:
                std::cerr<<"Warning: parameter type not found, given "<<type_id<<", should be one of F, L, I, S\n";
                break;
            }
        }
    }

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
    void mpi_broadcast() {
        for(auto iter=d_f64.begin(); iter!=d_f64.end(); iter++) PS::Comm::broadcast(&(iter->second->value), 1, 0);
        for(auto iter=d_i64.begin(); iter!=d_i64.end(); iter++) PS::Comm::broadcast(&(iter->second->value), 1, 0);
        for(auto iter=d_str.begin(); iter!=d_str.end(); iter++) {
            size_t str_size=iter->second->value.size();
            PS::Comm::broadcast(&str_size, 1, 0);
            char stmp[str_size+1];
            std::strcpy(stmp, iter->second->value.c_str());
            PS::Comm::broadcast(stmp, str_size, 0);
            iter->second->value = std::string(stmp);
        }
    }
#endif

    void print(std::ostream& os) const{
        for(auto iter=d_f64.begin(); iter!=d_f64.end(); iter++) os<<iter->first<<": "<<iter->second->value<<std::endl;
        for(auto iter=d_i64.begin(); iter!=d_i64.end(); iter++) os<<iter->first<<": "<<iter->second->value<<std::endl;
        for(auto iter=d_str.begin(); iter!=d_str.end(); iter++) os<<iter->first<<": "<<iter->second->value<<std::endl;
    }

    void printHelp(std::ostream& os, const bool print_format_info_flag=true, const bool print_all_flag=false, const int _offset_short_key=3, const int _offset_long_key=4, const int _width_key=23) const{
        IOParamsPrintHelp print_help(_offset_short_key, _offset_long_key, _width_key);
        if (print_format_info_flag) {
            os<<"** Default values are shown after ':'\n"
              <<"** The char in [] indicates argument type: ";
            print_help.printTypeShortNameDescription(os);
        }
        // short key
        for(auto iter=d_f64.begin(); iter!=d_f64.end(); iter++) iter->second->printHelp(os, print_help, true, print_all_flag);
        for(auto iter=d_i64.begin(); iter!=d_i64.end(); iter++) iter->second->printHelp(os, print_help, true, print_all_flag);
        for(auto iter=d_str.begin(); iter!=d_str.end(); iter++) iter->second->printHelp(os, print_help, true, print_all_flag);
        // long key
        for(auto iter=d_f64.begin(); iter!=d_f64.end(); iter++) iter->second->printHelp(os, print_help, false, print_all_flag);
        for(auto iter=d_i64.begin(); iter!=d_i64.end(); iter++) iter->second->printHelp(os, print_help, false, print_all_flag);
        for(auto iter=d_str.begin(); iter!=d_str.end(); iter++) iter->second->printHelp(os, print_help, false, print_all_flag);
    }
};


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
};
