#pragma once
#include <iostream>
#include <particle_simulator.hpp>
#include "soft_ptcl.hpp"

#ifdef GPU_PROFILE
#include "profile.hpp"
extern struct GPUProfile{
public:
    Tprofile copy;
    Tprofile send;
    Tprofile recv;
    Tprofile calc;
    const PS::S32 n_profile;

    GPUProfile(): 
        copy (Tprofile("copy       ")),
        send (Tprofile("send       ")),
        recv (Tprofile("receive    ")),
        calc (Tprofile("calc_force ")),
        n_profile(4) {}

	void print(std::ostream & fout, const PS::F64 time_sys, const PS::S64 n_loop=1){
        fout<<"Time: "<<time_sys<<std::endl;
        
        for(PS::S32 i=0; i<n_profile; i++) {
            Tprofile* iptr = (Tprofile*)this+i;
            iptr->print(fout, n_loop);
        }
    }

    void dump(std::ostream & fout, const PS::S32 width=20, const PS::S64 n_loop=1) const {
        for(PS::S32 i=0; i<n_profile; i++) {
            Tprofile* iptr = (Tprofile*)this+i;
            iptr->dump(fout, width, n_loop);
        }
    }

    void dumpName(std::ostream & fout, const PS::S32 width=20) {
        for(PS::S32 i=0; i<n_profile; i++) {
            Tprofile* iptr = (Tprofile*)this+i;
            iptr->dumpName(fout, width);
        }
    }
    
    void clear(){
        for(PS::S32 i=0; i<n_profile; i++) {
            Tprofile* iptr = (Tprofile*)this+i;
            iptr->reset();
        }
    }

} gpu_profile;
#endif

#ifdef USE_QUAD
#define SPJSoft PS::SPJQuadrupoleInAndOut
#else
#define SPJSoft PS::SPJMonopoleInAndOut
#endif

PS::S32 DispatchKernelWithSP(const PS::S32 tag,
                             const PS::S32 n_walk,
                             const EPISoft ** epi,
                             const PS::S32 *  n_epi,
                             const EPJSoft ** epj,
                             const PS::S32 *  n_epj,
                             const SPJSoft ** spj,
                             const PS::S32  * n_spj);

PS::S32 RetrieveKernel(const PS::S32 tag,
                       const PS::S32 n_walk,
                       const PS::S32 * ni,
                       ForceSoft      ** force);
