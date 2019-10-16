#include "petar.hpp"

int main(int argc, char *argv[]){

    PeTar petar;

    PS::S32 iread = petar.readParameters(argc,argv);
    if (iread<0) return 0;

    petar.evolveToTime(petar.time_end);

    return 0;

}
