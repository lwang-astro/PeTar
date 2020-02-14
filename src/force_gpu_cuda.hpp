#pragma once
#include <iostream>
#include <particle_simulator.hpp>
#include "soft_ptcl.hpp"
#include "soft_force.hpp"

PS::S32 DispatchKernelWithSP(const PS::S32 tag,
                             const PS::S32 n_walk,
                             const EPISoft ** epi,
                             const PS::S32 *  n_epi,
                             const EPJSoft ** epj,
                             const PS::S32 *  n_epj,
                             const PS::SPJMonopoleInAndOut ** spj,
                             const PS::S32  * n_spj);

PS::S32 RetrieveKernel(const PS::S32 tag,
                       const PS::S32 n_walk,
                       const PS::S32 * ni,
                       ForceSoft      ** force);
