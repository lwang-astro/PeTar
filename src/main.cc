#include "petar.hpp"
#ifdef GPERF_PROFILE
#include <gperftools/profiler.h>
#endif
#ifdef PETAR_USE_MPFRC
#include <mpreal.h>
#endif

int main(int argc, char *argv[]){

#ifdef NAN_CHECK_DEBUG
    assert(std::isnan(NAN));
#endif

    PeTar::initialFDPS(argc,argv);

#ifdef PETAR_USE_MPFRC
    mpfr::mpreal::set_default_prec(128); // Set precision to 128 bits
#endif

    PeTar petar;
    
    PS::S32 iread = petar.readParameters(argc,argv);
    if (iread<0) {
        PeTar::finalizeFDPS();
        return 0;
    }

    auto& inp = petar.input_parameters;

    if (inp.fname_inp.value=="__Plummer") petar.generatePlummer();
    //else if (inp.fname_inp.value!="__KeplerDisk") petar.generateKeplerDisk();
    else petar.readDataFromFile();

    petar.initialParameters();

    petar.initialStep();

#ifdef GPERF_PROFILE
    std::string rank_str;
    std::stringstream atmp;
    atmp<<petar.my_rank;
    atmp>>rank_str;
    std::string fproname=petar.input_parameters.fname_snp.value+".gperf.out.r"+rank_str;
    ProfilerStart(fproname.c_str());
#endif

#if 0
    PS::F64 dt_break = inp.dt_snp.value;
    PS::F64 dt_end = inp.time_end.value;
    PS::S32 n_loop = dt_end/dt_break;
    PS::F64 time_break = 0.0;
    
    for (int i=0; i<n_loop; i++) {   
        time_break += dt_break;
        int n_interupt = 1;
        while(n_interupt>0) n_interupt = petar.integrateToTime();
    }
#else
    int n_interupt = 1;
    while(n_interupt>0) n_interupt = petar.integrateToTime();
    
#endif

#ifdef GPERF_PROFILE
    ProfilerStop();
#endif

    PeTar::finalizeFDPS();

    return 0;

}
