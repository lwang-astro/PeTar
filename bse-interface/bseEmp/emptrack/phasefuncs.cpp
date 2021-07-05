#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>

typedef float     F32;
typedef double    F64;
typedef int       S32;
typedef long long S64;

#include "empfuncs.hpp"

#if defined(__cplusplus)
extern "C" {
#endif

    void followAGBPhase(F64 * aj,
			F64 * mass,
			F64 * mt,
			F64 * lum,
			F64 * r,
			F64 * rg,
			F64 * mcbagb,
			F64 * mc,
			F64 * mcx,
			F64 * mcmax);

    void calcTimestepAGBPhase(S32 * kw,
			      F64 * aj,
			      F64 * mass,
			      F64 * tn,
			      F64 * pts3,
			      F64 * dt,
			      F64 * dtr);

#if defined(__cplusplus)
}
#endif

void followAGBPhase(F64 * aj,
		    F64 * mass,
		    F64 * mt,
		    F64 * lum,
		    F64 * r,
		    F64 * rg,
		    F64 * mcbagb,
		    F64 * mc,
		    F64 * mcx,
		    F64 * mcmax) {
    static const F64 TinyValue = 1e-14;

    *mcbagb = getHeCoreMassBAGBTime(mass);
    *mc     = *mcbagb;

    *mcmax = getCOCoreMassEndTime(mass);
    if(*aj < getEndTime(mass)) {
	*mcx = *mcmax * 0.99;
    } else {
	*mcx = *mcmax;
    }

    F64 tBAGB = getHeITime(mass) + getTimeIntervalOfCHeBPhase(mass);
    F64 tBlue = getEndTimeOfBluePhase(mass);
    F64 tFin  = getEndTime(mass);
    if(askBlueOrRed(aj, mass) || askAllBlueOrNot(mt)) {
	F64 dt    = (tBlue - tBAGB != 0.) ? (*aj - tBAGB) : 0.;
	F64 tau   = dt / (tBlue - tBAGB + TinyValue);
	F64 tau3  = tau * tau * tau;
	F64 lBAGB = getLuminosityBAGBTime(mass);
	F64 rBAGB = getRadiusEndTimeOfBlueCHeBPhase(mt, mt, &lBAGB);
	F64 lBlue = getLuminosityEndTimeOfBluePhase(mass);
	F64 ltemp = pow(10., lBlue);
	F64 rBlue = (getEndTimeOfBluePhase(mt) == getEndTime(mt)) ?
	    getRadiusEndTime(mt) : getRadiusRedPhase(mt, &ltemp);
	F64 lpower = lBAGB + tau3 * (lBlue - lBAGB);
	F64 rpower = rBAGB + tau3 * (rBlue - rBAGB);
	*lum = pow(10., lpower);
	*r   = pow(10., rpower);
	*rg  = getRadiusRedPhase(mt, lum);
    } else {
	F64 tau   = (*aj - tBlue) / (tFin - tBlue);
	F64 lBlue = getLuminosityEndTimeOfBluePhase(mass);
	F64 lFin  = getLuminosityEndTime(mass);
	F64 lpower = lBlue + tau  * (lFin - lBlue);
	*lum = pow(10., lpower);
	F64 rpower = getRadiusRedPhase(mt, lum);
	*r   = pow(10., rpower);
	*rg  = *r;
    }

    return;
}

void calcTimestepAGBPhase(S32 * kw,
			  F64 * aj,
			  F64 * mass,
			  F64 * tn,
			  F64 * pts3,
			  F64 * dt,
			  F64 * dtr) {
    if((*kw) == 5) {
	F64 tBAGB = getHeITime(mass) + getTimeIntervalOfCHeBPhase(mass);
	F64 tFin  = getEndTime(mass);
	*dt  = (*pts3) * (tFin - tBAGB);
//	*dtr = std::min((*tn), tFin) - (*aj);
	*dtr = tFin - (*aj);
    }
    return;
}

