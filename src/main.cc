#include "petar.hpp"

int main(int argc, char *argv[]){

    PeTar petar;

    petar.initialFDPS(argc,argv);
    
    PS::S32 iread = petar.readParameters(argc,argv);
    if (iread<0) return 0;

    auto& inp = petar.input_parameters;

    if (inp.n_glb.value==0) petar.readDataFromFile();
    else petar.generatePlummer();

    petar.initialParameters();

    petar.initialStep();

#if 0
    PS::F64 dt_break = inp.dt_snp.value;
    PS::F64 dt_end = inp.time_end.value;
    PS::S32 n_loop = dt_end/dt_break;
    PS::F64 time_break = 0.0;
    
    for (int i=0; i<n_loop; i++) {   
        time_break += dt_break;
        petar.evolveToTime(time_break);
    }
#else
    petar.evolveToTime();
#endif

    return 0;

}
