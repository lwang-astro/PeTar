#include "petar.hpp"

int main(int argc, char *argv[]){

    PeTar petar;

    petar.initialFDPS(argc,argv);
    
    PS::S32 iread = petar.readParameters(argc,argv);
    if (iread<0) return 0;

    auto& inp = petar.input_parameters;

    PS::S32 reading_style = inp.reading_style.value;
    if (reading_style==1) petar.readDataFromFile();
    else if (reading_style==2) petar.generatePlummer();

    petar.initialParameters();

    petar.initialStep();

#if 1
    PS::F64 dt_break = petar.input_parameters.dt_snp.value;
    PS::F64 dt_end = petar.input_parameters.time_end.value;
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
