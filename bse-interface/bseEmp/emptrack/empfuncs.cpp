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

namespace ExtremeMetalPoors {
    F64 LowerLimitOfMassInScopeOfApplication = 8.;
    F64 UpperLimitOfMassInScopeOfApplication = 1e5;

    const  S64 NumberOfCoefficientForTime       = 4;
    //const  S64 NumberOfCoefficientForLuminosity = 6;
    const  S64 NumberOfCoefficientForLuminosity = 7;
    const  S64 NumberOfCoefficientForRadius     = 6;
    const  S64 NumberOfCoefficientForCoreMass   = 3;
    const  S64 NumberOfCoefficientForMSLuminosity = 3;
    const  S64 NumberOfCoefficientForMSRadius     = 4;
    const  S64 NumberOfCoefficientForAGBRadius    = 2;
    static bool InUse = false;
    static F64 zeta   = -100;

    static F64 UpperLimitOfMassHG          = 0.;
    static F64 UpperLimitOfMassBlueLoop    = 0.;
    static F64 LowerLimitOfMassAllBlue     = 0.;
    static F64 UpperLimitOfMassAllBlue     = 0.;
    static F64 UpperLimitOfMassAllBlueCHeB = 0.;
    static F64 atime[NumberOfCoefficientForTime];
    static F64 btime[NumberOfCoefficientForTime];
    static F64 ctim2[NumberOfCoefficientForTime];
    static F64 dtime[NumberOfCoefficientForTime];
    static F64 alums[NumberOfCoefficientForLuminosity];
    static F64 blums[NumberOfCoefficientForLuminosity];
    static F64 clums[NumberOfCoefficientForLuminosity];
    static F64 dlums[NumberOfCoefficientForLuminosity];
    static F64 arads[NumberOfCoefficientForRadius];
    static F64 brads[NumberOfCoefficientForRadius];
    static F64 crads[NumberOfCoefficientForRadius];
    static F64 drads[NumberOfCoefficientForRadius];
    static F64 acore[NumberOfCoefficientForCoreMass];
    static F64 bcore[NumberOfCoefficientForCoreMass];
    static F64 ccore[NumberOfCoefficientForCoreMass];
    static F64 dcore[NumberOfCoefficientForCoreMass];
    static F64 msalums[NumberOfCoefficientForMSLuminosity];
    static F64 msblums[NumberOfCoefficientForMSLuminosity];
    static F64 msclums[NumberOfCoefficientForMSLuminosity];
    static F64 msdlums[NumberOfCoefficientForMSLuminosity];
    static F64 msarads[NumberOfCoefficientForMSRadius];
    static F64 msbrads[NumberOfCoefficientForMSRadius];
    static F64 mscrads[NumberOfCoefficientForMSRadius];
    static F64 msdrads[NumberOfCoefficientForMSRadius];
    static F64 agbarads[NumberOfCoefficientForAGBRadius];
    static F64 agbbrads[NumberOfCoefficientForAGBRadius];

    F64 CriticalMassMassive = 160.;
    const  S64 NumberOfCoefficientForTimeMassive       = 4;
    const  S64 NumberOfCoefficientForLuminosityMassive = 6;
    const  S64 NumberOfCoefficientForRadiusMassive     = 5;
    const  S64 NumberOfCoefficientForCoreMassMassive   = 3;
    const  S64 NumberOfCoefficientForMSRadiusMassive   = 4;
    static F64 atimeMassive[NumberOfCoefficientForTimeMassive];
    static F64 btimeMassive[NumberOfCoefficientForTimeMassive];
    static F64 ctimeMassive[NumberOfCoefficientForTimeMassive];
    static F64 dtimeMassive[NumberOfCoefficientForTimeMassive];
    static F64 alumsMassive[NumberOfCoefficientForLuminosityMassive];
    static F64 blumsMassive[NumberOfCoefficientForLuminosityMassive];
    static F64 clumsMassive[NumberOfCoefficientForLuminosityMassive];
    static F64 dlumsMassive[NumberOfCoefficientForLuminosityMassive];
    static F64 aradsMassive[NumberOfCoefficientForRadiusMassive];
    static F64 bradsMassive[NumberOfCoefficientForRadiusMassive];
    static F64 cradsMassive[NumberOfCoefficientForRadiusMassive];
    static F64 dradsMassive[NumberOfCoefficientForRadiusMassive];
    static F64 acoreMassive[NumberOfCoefficientForCoreMassMassive];
    static F64 bcoreMassive[NumberOfCoefficientForCoreMassMassive];
    static F64 ccoreMassive[NumberOfCoefficientForCoreMassMassive];
    static F64 dcoreMassive[NumberOfCoefficientForCoreMassMassive];
    static F64 msaradsMassive[NumberOfCoefficientForMSRadiusMassive];
    static F64 msbradsMassive[NumberOfCoefficientForMSRadiusMassive];
    static F64 mscradsMassive[NumberOfCoefficientForMSRadiusMassive];
    static F64 msdradsMassive[NumberOfCoefficientForMSRadiusMassive];

    bool UseK2strOrNot = false;
    const  S64 NumberOfCoefficientForK2str = 2;
    static F64 ak2str[NumberOfCoefficientForK2str];
    static F64 bk2str[NumberOfCoefficientForK2str];

    inline F64 calcTauOne(F64 tau) {
	return std::min(1., tau);
    }

    inline F64 calcTauTwo(F64 tau) {
	static const F64 eps = 0.01;
	F64 value = (tau - (1. - eps)) / eps;
	return std::max(0., std::min(1., value));
    }

    inline F64 getTimeAndTimeInterval(F64 mass,
				      S64 index) {
	F64 minv1 = 1. / mass;
	F64 minv2 = minv1 * minv1;
	F64 minv3 = minv1 * minv2;
	if(mass <= CriticalMassMassive) {
	    return atime[index] + btime[index] * minv1 + ctim2[index] * minv2 + dtime[index] * minv3;
	} else {
	    return atimeMassive[index] + btimeMassive[index] * minv1;
	}
    }
    
    inline F64 getLuminosity(F64 mass,
			     S64 index) {
	F64 lmass1 = log10(mass);
	F64 lmass2 = lmass1 * lmass1;
	F64 lmass3 = lmass1 * lmass2;
	if(mass <= CriticalMassMassive) {
	    return alums[index] + blums[index] * lmass1 + clums[index] * lmass2 + dlums[index] * lmass3;
	} else {
	    index = (index != 6) ? index : 5;
	    return alumsMassive[index] + blumsMassive[index] * lmass1 + clumsMassive[index] * lmass2 + dlumsMassive[index] * lmass3;
	}
    }

    inline F64 getRadius(F64 mass,
			 S64 index) {
	F64 lmass1 = log10(mass);
	F64 lmass2 = lmass1 * lmass1;
	F64 lmass3 = lmass1 * lmass2;
	if(mass <= CriticalMassMassive) {
	    return arads[index] + brads[index] * lmass1 + crads[index] * lmass2 + drads[index] * lmass3;
	} else {
	    index = (index != 5) ? index : 4;
	    return aradsMassive[index] + bradsMassive[index] * lmass1 + cradsMassive[index] * lmass2 + dradsMassive[index] * lmass3;
	}
    }

    inline F64 getCoreMass(F64 mass,
			   S64 index) {
	F64 lmass1 = log10(mass);
	F64 lmass2 = lmass1 * lmass1;
	F64 lmass3 = lmass1 * lmass2;
#if 0
	F64 power = acore[index] + bcore[index] * lmass1 + ccore[index] * lmass2 + dcore[index] * lmass3;
#else
	F64 power = 0.;
	if(mass <= CriticalMassMassive) {
	    power = acore[index] + bcore[index] * lmass1 + ccore[index] * lmass2 + dcore[index] * lmass3;
	} else {
	    power = acoreMassive[index] + bcoreMassive[index] * lmass1 + ccoreMassive[index] * lmass2 + dcoreMassive[index] * lmass3;
	}
#endif
	return pow(10., power);
    }

    inline F64 getMSLuminosity(F64 mass,
			       S64 index) {
	F64 value = 0.;
	if(index != NumberOfCoefficientForMSLuminosity - 1) {
	    F64 lmass1 = log10(mass);
	    F64 lmassi = 1. / lmass1;
	    F64 lmass2 = lmass1 * lmass1;
	    value = msalums[index] * lmassi + msblums[index]
		+ msclums[index] * lmass1 + msdlums[index] * lmass2;
	} else {
	    F64 lmass1 = log10(mass);
	    F64 lmass2 = lmass1 * lmass1;
	    F64 lmass3 = lmass1 * lmass2;
	    value = msalums[index] + msblums[index] * lmass1
		+ msclums[index] * lmass2 + msdlums[index] * lmass3;
	}
	return value;
    }

    inline F64 getMSRadius(F64 mass,
			   S64 index) {
	F64 value = 0.;
	if(index != NumberOfCoefficientForMSRadius - 1) {
	    F64 lmass1 = log10(mass);
	    F64 lmassi = 1. / lmass1;
	    F64 lmass2 = lmass1 * lmass1;
	    value = msarads[index] * lmassi + msbrads[index]
		+ mscrads[index] * lmass1 + msdrads[index] * lmass2;
	    /*
	    if(mass <= CriticalMassMassive) {
		value = msarads[index] * lmassi + msbrads[index]
		    + mscrads[index] * lmass1 + msdrads[index] * lmass2;
	    } else {
		value = msaradsMassive[index] * lmassi + msbradsMassive[index]
		    + mscradsMassive[index] * lmass1 + msdradsMassive[index] * lmass2;
	    }
	    */
	} else {
	    F64 lmass1 = log10(mass);
	    F64 lmass2 = lmass1 * lmass1;
	    F64 lmass3 = lmass1 * lmass2;
	    value = msarads[index] + msbrads[index] * lmass1
		+ mscrads[index] * lmass2 + msdrads[index] * lmass3;
	    /*
	    if(mass <= CriticalMassMassive) {
		value = msarads[index] + msbrads[index] * lmass1
		    + mscrads[index] * lmass2 + msdrads[index] * lmass3;
	    } else {
		value = msaradsMassive[index] + msbradsMassive[index] * lmass1
		    + mscradsMassive[index] * lmass2 + msdradsMassive[index] * lmass3;
	    }
	    */
	}
	return value;
    }

    inline F64 getQuantityEndTimeOfBluePhase(F64 mass,
					     F64 qBlue,
					     F64 qBAGB,
					     F64 qFin) {
	F64 tHeI  = getTimeAndTimeInterval(mass, 0);
	F64 tHef  = getTimeAndTimeInterval(mass, 1);
	F64 tBlue = getTimeAndTimeInterval(mass, 2);
#if 0 // 20/01/30
	if(tBlue < tHeI) {
	    fprintf(stderr, "Error: tBlue is smaller than tHeI.\n");
	    exit(0);
	} else if(LowerLimitOfMassAllBlue <= mass && mass < UpperLimitOfMassAllBlue)  {
	    qBlue = qFin;
	} else if((mass < UpperLimitOfMassAllBlueCHeB && tBlue <= tHeI + tHef)
		|| (mass >= UpperLimitOfMassAllBlueCHeB && tBlue >= tHeI + tHef)) {
	    qBlue = qBAGB;
	}
#else
	if(LowerLimitOfMassAllBlue <= mass && mass < UpperLimitOfMassAllBlue)  {
	    qBlue = qFin;
	} else if((mass < UpperLimitOfMassAllBlueCHeB && tBlue <= tHeI + tHef)
		|| (mass >= UpperLimitOfMassAllBlueCHeB && tBlue >= tHeI + tHef)) {
	    qBlue = qBAGB;
	}
#endif
	return qBlue;
    }

    inline F64 getAGBRadius(F64 mt,
			    S64 index) {
	F64 lmt = log10(mt);
	return agbarads[index] + agbbrads[index] * lmt;
    }

    inline F64 interpolateCoefficient(F64 zeta,
				      F64 val1,
				      F64 val2) {
	return ((val1 - val2) * (zeta - (-2.)) + val2);
    }

    void setMetallicityInterpolation(F64 zeta) {
	F64 atimetmp1[NumberOfCoefficientForTime];
	F64 btimetmp1[NumberOfCoefficientForTime];
	F64 ctimetmp1[NumberOfCoefficientForTime];
	F64 dtimetmp1[NumberOfCoefficientForTime];
	F64 alumstmp1[NumberOfCoefficientForLuminosity];
	F64 blumstmp1[NumberOfCoefficientForLuminosity];
	F64 clumstmp1[NumberOfCoefficientForLuminosity];
	F64 dlumstmp1[NumberOfCoefficientForLuminosity];
	F64 aradstmp1[NumberOfCoefficientForRadius];
	F64 bradstmp1[NumberOfCoefficientForRadius];
	F64 cradstmp1[NumberOfCoefficientForRadius];
	F64 dradstmp1[NumberOfCoefficientForRadius];
	F64 acoretmp1[NumberOfCoefficientForCoreMass];
	F64 bcoretmp1[NumberOfCoefficientForCoreMass];
	F64 ccoretmp1[NumberOfCoefficientForCoreMass];
	F64 dcoretmp1[NumberOfCoefficientForCoreMass];
	F64 msalumstmp1[NumberOfCoefficientForMSLuminosity];
	F64 msblumstmp1[NumberOfCoefficientForMSLuminosity];
	F64 msclumstmp1[NumberOfCoefficientForMSLuminosity];
	F64 msdlumstmp1[NumberOfCoefficientForMSLuminosity];
	F64 msaradstmp1[NumberOfCoefficientForMSRadius];
	F64 msbradstmp1[NumberOfCoefficientForMSRadius];
	F64 mscradstmp1[NumberOfCoefficientForMSRadius];
	F64 msdradstmp1[NumberOfCoefficientForMSRadius];
	F64 agbaradstmp1[NumberOfCoefficientForAGBRadius];
	F64 agbbradstmp1[NumberOfCoefficientForAGBRadius];
	F64 atimeMassivetmp1[NumberOfCoefficientForTimeMassive];
	F64 btimeMassivetmp1[NumberOfCoefficientForTimeMassive];
	F64 ctimeMassivetmp1[NumberOfCoefficientForTimeMassive];
	F64 dtimeMassivetmp1[NumberOfCoefficientForTimeMassive];
	F64 alumsMassivetmp1[NumberOfCoefficientForLuminosityMassive];
	F64 blumsMassivetmp1[NumberOfCoefficientForLuminosityMassive];
	F64 clumsMassivetmp1[NumberOfCoefficientForLuminosityMassive];
	F64 dlumsMassivetmp1[NumberOfCoefficientForLuminosityMassive];
	F64 aradsMassivetmp1[NumberOfCoefficientForRadiusMassive];
	F64 bradsMassivetmp1[NumberOfCoefficientForRadiusMassive];
	F64 cradsMassivetmp1[NumberOfCoefficientForRadiusMassive];
	F64 dradsMassivetmp1[NumberOfCoefficientForRadiusMassive];
	F64 acoreMassivetmp1[NumberOfCoefficientForCoreMassMassive];
	F64 bcoreMassivetmp1[NumberOfCoefficientForCoreMassMassive];
	F64 ccoreMassivetmp1[NumberOfCoefficientForCoreMassMassive];
	F64 dcoreMassivetmp1[NumberOfCoefficientForCoreMassMassive];
	F64 msaradsMassivetmp1[NumberOfCoefficientForMSRadiusMassive];
	F64 msbradsMassivetmp1[NumberOfCoefficientForMSRadiusMassive];
	F64 mscradsMassivetmp1[NumberOfCoefficientForMSRadiusMassive];
	F64 msdradsMassivetmp1[NumberOfCoefficientForMSRadiusMassive];

	F64 atimetmp2[NumberOfCoefficientForTime];
	F64 btimetmp2[NumberOfCoefficientForTime];
	F64 ctimetmp2[NumberOfCoefficientForTime];
	F64 dtimetmp2[NumberOfCoefficientForTime];
	F64 alumstmp2[NumberOfCoefficientForLuminosity];
	F64 blumstmp2[NumberOfCoefficientForLuminosity];
	F64 clumstmp2[NumberOfCoefficientForLuminosity];
	F64 dlumstmp2[NumberOfCoefficientForLuminosity];
	F64 aradstmp2[NumberOfCoefficientForRadius];
	F64 bradstmp2[NumberOfCoefficientForRadius];
	F64 cradstmp2[NumberOfCoefficientForRadius];
	F64 dradstmp2[NumberOfCoefficientForRadius];
	F64 acoretmp2[NumberOfCoefficientForCoreMass];
	F64 bcoretmp2[NumberOfCoefficientForCoreMass];
	F64 ccoretmp2[NumberOfCoefficientForCoreMass];
	F64 dcoretmp2[NumberOfCoefficientForCoreMass];
	F64 msalumstmp2[NumberOfCoefficientForMSLuminosity];
	F64 msblumstmp2[NumberOfCoefficientForMSLuminosity];
	F64 msclumstmp2[NumberOfCoefficientForMSLuminosity];
	F64 msdlumstmp2[NumberOfCoefficientForMSLuminosity];
	F64 msaradstmp2[NumberOfCoefficientForMSRadius];
	F64 msbradstmp2[NumberOfCoefficientForMSRadius];
	F64 mscradstmp2[NumberOfCoefficientForMSRadius];
	F64 msdradstmp2[NumberOfCoefficientForMSRadius];
	F64 agbaradstmp2[NumberOfCoefficientForAGBRadius];
	F64 agbbradstmp2[NumberOfCoefficientForAGBRadius];
	F64 atimeMassivetmp2[NumberOfCoefficientForTimeMassive];
	F64 btimeMassivetmp2[NumberOfCoefficientForTimeMassive];
	F64 ctimeMassivetmp2[NumberOfCoefficientForTimeMassive];
	F64 dtimeMassivetmp2[NumberOfCoefficientForTimeMassive];
	F64 alumsMassivetmp2[NumberOfCoefficientForLuminosityMassive];
	F64 blumsMassivetmp2[NumberOfCoefficientForLuminosityMassive];
	F64 clumsMassivetmp2[NumberOfCoefficientForLuminosityMassive];
	F64 dlumsMassivetmp2[NumberOfCoefficientForLuminosityMassive];
	F64 aradsMassivetmp2[NumberOfCoefficientForRadiusMassive];
	F64 bradsMassivetmp2[NumberOfCoefficientForRadiusMassive];
	F64 cradsMassivetmp2[NumberOfCoefficientForRadiusMassive];
	F64 dradsMassivetmp2[NumberOfCoefficientForRadiusMassive];
	F64 acoreMassivetmp2[NumberOfCoefficientForCoreMassMassive];
	F64 bcoreMassivetmp2[NumberOfCoefficientForCoreMassMassive];
	F64 ccoreMassivetmp2[NumberOfCoefficientForCoreMassMassive];
	F64 dcoreMassivetmp2[NumberOfCoefficientForCoreMassMassive];
	F64 msaradsMassivetmp2[NumberOfCoefficientForMSRadiusMassive];
	F64 msbradsMassivetmp2[NumberOfCoefficientForMSRadiusMassive];
	F64 mscradsMassivetmp2[NumberOfCoefficientForMSRadiusMassive];
	F64 msdradsMassivetmp2[NumberOfCoefficientForMSRadiusMassive];
	
#ifdef GENEVAMODEL
	char idir[1024] = "ffgeneva";
	fprintf(stderr, " Use Geneva interpolation\n");
#else
	char idir[1024] = "ffbonn";
	fprintf(stderr, " Use Bonn interpolation\n");
#endif
	char ifile[1024];
	sprintf(ifile, "./%s/metal1e-1.dat", idir);
	FILE * fp;
	fp = fopen(ifile, "r");
	assert(fp);
	fscanf(fp, "%lf", &UpperLimitOfMassHG);
	fscanf(fp, "%lf", &UpperLimitOfMassBlueLoop);
	fscanf(fp, "%lf", &LowerLimitOfMassAllBlue);
	fscanf(fp, "%lf", &UpperLimitOfMassAllBlue);
	fscanf(fp, "%lf", &UpperLimitOfMassAllBlueCHeB);
	for(S64 i = 0; i < NumberOfCoefficientForTime; i++) {
	    fscanf(fp, "%lf%lf%lf%lf", &atimetmp1[i], &btimetmp1[i], &ctimetmp1[i], &dtimetmp1[i]);
	}
	for(S64 i = 0; i < NumberOfCoefficientForLuminosity; i++) {
	    fscanf(fp, "%lf%lf%lf%lf", &alumstmp1[i], &blumstmp1[i], &clumstmp1[i], &dlumstmp1[i]);
	}
	for(S64 i = 0; i < NumberOfCoefficientForRadius; i++) {
	    fscanf(fp, "%lf%lf%lf%lf", &aradstmp1[i], &bradstmp1[i], &cradstmp1[i], &dradstmp1[i]);
	}
	for(S64 i = 0; i < NumberOfCoefficientForCoreMass; i++) {
	    fscanf(fp, "%lf%lf%lf%lf", &acoretmp1[i], &bcoretmp1[i], &ccoretmp1[i], &dcoretmp1[i]);
	}
	for(S64 i = 0; i < NumberOfCoefficientForMSLuminosity; i++) {
	    fscanf(fp, "%lf%lf%lf%lf", &msalumstmp1[i], &msblumstmp1[i], &msclumstmp1[i], &msdlumstmp1[i]);
	}
	for(S64 i = 0; i < NumberOfCoefficientForMSRadius; i++) {
	    fscanf(fp, "%lf%lf%lf%lf", &msaradstmp1[i], &msbradstmp1[i], &mscradstmp1[i], &msdradstmp1[i]);
	}
	for(S64 i = 0; i < NumberOfCoefficientForAGBRadius; i++) {
	    fscanf(fp, "%lf%lf", &agbaradstmp1[i], &agbbradstmp1[i]);
	}
	for(S64 i = 0; i < NumberOfCoefficientForTimeMassive; i++) {
	    fscanf(fp, "%lf%lf%lf%lf", &atimeMassivetmp1[i], &btimeMassivetmp1[i], &ctimeMassivetmp1[i], &dtimeMassivetmp1[i]);
	}
	for(S64 i = 0; i < NumberOfCoefficientForLuminosityMassive; i++) {
	    fscanf(fp, "%lf%lf%lf%lf", &alumsMassivetmp1[i], &blumsMassivetmp1[i], &clumsMassivetmp1[i], &dlumsMassivetmp1[i]);
	}	
	for(S64 i = 0; i < NumberOfCoefficientForRadiusMassive; i++) {
	    fscanf(fp, "%lf%lf%lf%lf", &aradsMassivetmp1[i], &bradsMassivetmp1[i], &cradsMassivetmp1[i], &dradsMassivetmp1[i]);
	}	
	for(S64 i = 0; i < NumberOfCoefficientForCoreMassMassive; i++) {
	    fscanf(fp, "%lf%lf%lf%lf", &acoreMassivetmp1[i], &bcoreMassivetmp1[i], &ccoreMassivetmp1[i], &dcoreMassivetmp1[i]);
	}
	fclose(fp);

	sprintf(ifile, "./%s/metal1e-2.dat", idir);
	fp = fopen(ifile, "r");
	assert(fp);
	fscanf(fp, "%lf", &UpperLimitOfMassHG);
	fscanf(fp, "%lf", &UpperLimitOfMassBlueLoop);
	fscanf(fp, "%lf", &LowerLimitOfMassAllBlue);
	fscanf(fp, "%lf", &UpperLimitOfMassAllBlue);
	fscanf(fp, "%lf", &UpperLimitOfMassAllBlueCHeB);
	for(S64 i = 0; i < NumberOfCoefficientForTime; i++) {
	    fscanf(fp, "%lf%lf%lf%lf", &atimetmp2[i], &btimetmp2[i], &ctimetmp2[i], &dtimetmp2[i]);
	}
	for(S64 i = 0; i < NumberOfCoefficientForLuminosity; i++) {
	    fscanf(fp, "%lf%lf%lf%lf", &alumstmp2[i], &blumstmp2[i], &clumstmp2[i], &dlumstmp2[i]);
	}
	for(S64 i = 0; i < NumberOfCoefficientForRadius; i++) {
	    fscanf(fp, "%lf%lf%lf%lf", &aradstmp2[i], &bradstmp2[i], &cradstmp2[i], &dradstmp2[i]);
	}
	for(S64 i = 0; i < NumberOfCoefficientForCoreMass; i++) {
	    fscanf(fp, "%lf%lf%lf%lf", &acoretmp2[i], &bcoretmp2[i], &ccoretmp2[i], &dcoretmp2[i]);
	}
	for(S64 i = 0; i < NumberOfCoefficientForMSLuminosity; i++) {
	    fscanf(fp, "%lf%lf%lf%lf", &msalumstmp2[i], &msblumstmp2[i], &msclumstmp2[i], &msdlumstmp2[i]);
	}
	for(S64 i = 0; i < NumberOfCoefficientForMSRadius; i++) {
	    fscanf(fp, "%lf%lf%lf%lf", &msaradstmp2[i], &msbradstmp2[i], &mscradstmp2[i], &msdradstmp2[i]);
	}
	for(S64 i = 0; i < NumberOfCoefficientForAGBRadius; i++) {
	    fscanf(fp, "%lf%lf", &agbaradstmp2[i], &agbbradstmp2[i]);
	}
	for(S64 i = 0; i < NumberOfCoefficientForTimeMassive; i++) {
	    fscanf(fp, "%lf%lf%lf%lf", &atimeMassivetmp2[i], &btimeMassivetmp2[i], &ctimeMassivetmp2[i], &dtimeMassivetmp2[i]);
	}
	for(S64 i = 0; i < NumberOfCoefficientForLuminosityMassive; i++) {
	    fscanf(fp, "%lf%lf%lf%lf", &alumsMassivetmp2[i], &blumsMassivetmp2[i], &clumsMassivetmp2[i], &dlumsMassivetmp2[i]);
	}	
	for(S64 i = 0; i < NumberOfCoefficientForRadiusMassive; i++) {
	    fscanf(fp, "%lf%lf%lf%lf", &aradsMassivetmp2[i], &bradsMassivetmp2[i], &cradsMassivetmp2[i], &dradsMassivetmp2[i]);
	}	
	for(S64 i = 0; i < NumberOfCoefficientForCoreMassMassive; i++) {
	    fscanf(fp, "%lf%lf%lf%lf", &acoreMassivetmp2[i], &bcoreMassivetmp2[i], &ccoreMassivetmp2[i], &dcoreMassivetmp2[i]);
	}
	fclose(fp);

	for(S64 i = 0; i < NumberOfCoefficientForTime; i++) {
	    atime[i] = interpolateCoefficient(zeta, atimetmp1[i], atimetmp2[i]);
	    btime[i] = interpolateCoefficient(zeta, btimetmp1[i], btimetmp2[i]);
	    ctim2[i] = interpolateCoefficient(zeta, ctimetmp1[i], ctimetmp2[i]);
	    dtime[i] = interpolateCoefficient(zeta, dtimetmp1[i], dtimetmp2[i]);
	}
	for(S64 i = 0; i < NumberOfCoefficientForLuminosity; i++) {
	    alums[i] = interpolateCoefficient(zeta, alumstmp1[i], alumstmp2[i]);
	    blums[i] = interpolateCoefficient(zeta, blumstmp1[i], blumstmp2[i]);
	    clums[i] = interpolateCoefficient(zeta, clumstmp1[i], clumstmp2[i]);
	    dlums[i] = interpolateCoefficient(zeta, dlumstmp1[i], dlumstmp2[i]);
	}
	for(S64 i = 0; i < NumberOfCoefficientForRadius; i++) {
	    arads[i] = interpolateCoefficient(zeta, aradstmp1[i], aradstmp2[i]);
	    brads[i] = interpolateCoefficient(zeta, bradstmp1[i], bradstmp2[i]);
	    crads[i] = interpolateCoefficient(zeta, cradstmp1[i], cradstmp2[i]);
	    drads[i] = interpolateCoefficient(zeta, dradstmp1[i], dradstmp2[i]);
	}
	for(S64 i = 0; i < NumberOfCoefficientForCoreMass; i++) {
	    acore[i] = interpolateCoefficient(zeta, acoretmp1[i], acoretmp2[i]);
	    bcore[i] = interpolateCoefficient(zeta, bcoretmp1[i], bcoretmp2[i]);
	    ccore[i] = interpolateCoefficient(zeta, ccoretmp1[i], ccoretmp2[i]);
	    dcore[i] = interpolateCoefficient(zeta, dcoretmp1[i], dcoretmp2[i]);
	}
	for(S64 i = 0; i < NumberOfCoefficientForMSLuminosity; i++) {
	    msalums[i] = interpolateCoefficient(zeta, msalumstmp1[i], msalumstmp2[i]);
	    msblums[i] = interpolateCoefficient(zeta, msblumstmp1[i], msblumstmp2[i]);
	    msclums[i] = interpolateCoefficient(zeta, msclumstmp1[i], msclumstmp2[i]);
	    msdlums[i] = interpolateCoefficient(zeta, msdlumstmp1[i], msdlumstmp2[i]);
	}
	for(S64 i = 0; i < NumberOfCoefficientForMSRadius; i++) {
	    msarads[i] = interpolateCoefficient(zeta, msaradstmp1[i], msaradstmp2[i]);
	    msbrads[i] = interpolateCoefficient(zeta, msbradstmp1[i], msbradstmp2[i]);
	    mscrads[i] = interpolateCoefficient(zeta, mscradstmp1[i], mscradstmp2[i]);
	    msdrads[i] = interpolateCoefficient(zeta, msdradstmp1[i], msdradstmp2[i]);
	}
	for(S64 i = 0; i < NumberOfCoefficientForAGBRadius; i++) {
	    agbarads[i] = interpolateCoefficient(zeta, agbaradstmp1[i], agbaradstmp2[i]);
	    agbbrads[i] = interpolateCoefficient(zeta, agbbradstmp1[i], agbbradstmp2[i]);
	}
	for(S64 i = 0; i < NumberOfCoefficientForTimeMassive; i++) {
	    atimeMassive[i] = interpolateCoefficient(zeta, atimeMassivetmp1[i], atimeMassivetmp2[i]);
	    btimeMassive[i] = interpolateCoefficient(zeta, btimeMassivetmp1[i], btimeMassivetmp2[i]);
	    ctimeMassive[i] = interpolateCoefficient(zeta, ctimeMassivetmp1[i], ctimeMassivetmp2[i]);
	    dtimeMassive[i] = interpolateCoefficient(zeta, dtimeMassivetmp1[i], dtimeMassivetmp2[i]);
	}
	for(S64 i = 0; i < NumberOfCoefficientForLuminosityMassive; i++) {
	    alumsMassive[i] = interpolateCoefficient(zeta, alumsMassivetmp1[i], alumsMassivetmp2[i]);
	    blumsMassive[i] = interpolateCoefficient(zeta, blumsMassivetmp1[i], blumsMassivetmp2[i]);
	    clumsMassive[i] = interpolateCoefficient(zeta, clumsMassivetmp1[i], clumsMassivetmp2[i]);
	    dlumsMassive[i] = interpolateCoefficient(zeta, dlumsMassivetmp1[i], dlumsMassivetmp2[i]);
	}	
	for(S64 i = 0; i < NumberOfCoefficientForRadiusMassive; i++) {
	    aradsMassive[i] = interpolateCoefficient(zeta, aradsMassivetmp1[i], aradsMassivetmp2[i]);
	    bradsMassive[i] = interpolateCoefficient(zeta, bradsMassivetmp1[i], bradsMassivetmp2[i]);
	    cradsMassive[i] = interpolateCoefficient(zeta, cradsMassivetmp1[i], cradsMassivetmp2[i]);
	    dradsMassive[i] = interpolateCoefficient(zeta, dradsMassivetmp1[i], dradsMassivetmp2[i]);
	}	
	for(S64 i = 0; i < NumberOfCoefficientForCoreMassMassive; i++) {
	    acoreMassive[i] = interpolateCoefficient(zeta, acoreMassivetmp1[i], acoreMassivetmp2[i]);
	    bcoreMassive[i] = interpolateCoefficient(zeta, bcoreMassivetmp1[i], bcoreMassivetmp2[i]);
	    ccoreMassive[i] = interpolateCoefficient(zeta, ccoreMassivetmp1[i], ccoreMassivetmp2[i]);
	    dcoreMassive[i] = interpolateCoefficient(zeta, dcoreMassivetmp1[i], dcoreMassivetmp2[i]);
	}
#ifdef K2_YOSHIDA
	char afile[1024];
	sprintf(afile, "./%s/additional1e-1.dat", idir);
	fp = fopen(afile, "r");
	if(fp != NULL) {
	    F64 ak2strtmp1[NumberOfCoefficientForK2str], bk2strtmp1[NumberOfCoefficientForK2str];
	    F64 ak2strtmp2[NumberOfCoefficientForK2str], bk2strtmp2[NumberOfCoefficientForK2str];
	    fscanf(fp, "%lf%lf", &ak2strtmp1[0], &bk2strtmp1[0]);
	    fscanf(fp, "%lf%lf", &ak2strtmp1[1], &bk2strtmp1[1]);
	    fclose(fp);
	    sprintf(afile, "./%s/additional1e-2.dat", idir);
	    fp = fopen(afile, "r");
	    assert(fp);
	    fscanf(fp, "%lf%lf", &ak2strtmp2[0], &bk2strtmp2[0]);
	    fscanf(fp, "%lf%lf", &ak2strtmp2[1], &bk2strtmp2[1]);
	    fclose(fp);
	    for(S64 i = 0; i < NumberOfCoefficientForK2str; i++) {
		ak2str[i] = interpolateCoefficient(zeta, ak2strtmp1[i], ak2strtmp2[i]);
		bk2str[i] = interpolateCoefficient(zeta, bk2strtmp1[i], bk2strtmp2[i]);
	    }
	    fprintf(stderr, " Use additional file\n");
	    UseK2strOrNot = true;
	}
#endif
    }

    void checkTimeSequence(F64 mass) {
	F64 tHeI  = getHeITime(&mass);
	F64 tBAGB = tHeI + getTimeIntervalOfCHeBPhase(&mass);
	F64 tBlue = getEndTimeOfBluePhase(&mass);
	F64 tFin  = getEndTime(&mass);
	if(mass < LowerLimitOfMassAllBlue) {
	    if(!(tBAGB <= tBlue && tBlue <= tFin)) {
		fprintf(stderr, "Error: m: %e tBAGB: %e tBlue: %e tFin: %e\n", mass, tBAGB, tBlue, tFin);
		assert(NULL);
	    }
	} else if(mass < UpperLimitOfMassAllBlue) {
	    if(!(tBAGB <= tBlue && tBlue == tFin)) {
		fprintf(stderr, "Error: m: %e tBAGB: %e tBlue: %e tFin: %e\n", mass, tBAGB, tBlue, tFin);
		assert(NULL);
	    }
	} else if(mass < UpperLimitOfMassAllBlueCHeB) {
	    if(!(tBAGB <= tBlue && tBlue <= tFin)) {
		fprintf(stderr, "Error: m: %e tBAGB: %e tBlue: %e tFin: %e\n", mass, tBAGB, tBlue, tFin);
		assert(NULL);
	    }
	} else {
	    if(!(tBlue <= tBAGB && tBAGB <= tFin)) {
		fprintf(stderr, "Error: m: %e tBlue: %e tBAGB: %e tFin: %e\n", mass, tBlue, tBAGB, tFin);
		assert(NULL);
	    }
	}
    }

    void checkLuminosity(F64 mass) {
	F64 tHeI  = getHeITime(&mass);
	F64 tBAGB = tHeI + getTimeIntervalOfCHeBPhase(&mass);
	F64 tBlue = getEndTimeOfBluePhase(&mass);
	F64 tFin  = getEndTime(&mass);

	F64 lTMS  = getLuminosityTMSTime(&mass);
	F64 lHeI  = getLuminosityHeITime(&mass);
	F64 lBAGB = getLuminosityBAGBTime(&mass);
	F64 lBlue = getLuminosityEndTimeOfBluePhase(&mass);
	F64 lFin  = getLuminosityEndTime(&mass);

	if(mass >= UpperLimitOfMassHG) {
	    assert(lTMS == lHeI);
	}
	if(tBAGB == tBlue) {
	    assert(lBAGB == lBlue);
	}
	if(tBAGB == tFin) {
	    assert(lBAGB == lFin);
	}
	/*
	if(tBlue == tFin) {
	    assert(lBlue == lFin);
	}
	*/

    }

    void checkRadius(F64 mass) {
	F64 tHeI  = getHeITime(&mass);
	F64 tBAGB = tHeI + getTimeIntervalOfCHeBPhase(&mass);
	F64 tFin  = getEndTime(&mass);

	F64 rTMS  = getRadiusTMSTime(&mass);
	F64 rHeI  = getRadiusHeITime(&mass);
	F64 rMin  = getMinimumRadiusCHeBPhase(&mass);
	F64 rBAGB = getRadiusEndTimeOfBlueCHeBPhase(&mass, &mass, &mass);
	F64 rFin  = getRadiusEndTime(&mass);

	if(mass >= UpperLimitOfMassHG) {
	    assert(rTMS == rHeI);
	}
	if(mass >= UpperLimitOfMassBlueLoop) {
	    assert(rMin == rHeI);
	}

	if(tBAGB == tFin) {
	    assert(rBAGB == rFin);
	}

    }

    void checkMSHGPhase(char * mode) {
	bool flag = (strcmp(mode, "MSHGLuminosity") == 0) ? true : false;
	const S64 nsample = 14;
	static F64 m[nsample] = {8., 10., 13., 16., 20., 25., 32.,
				 40., 50., 65., 80., 100., 125., 160.};
	for(S64 i = 0; i < nsample; i++) {
	    char ifile[1024];
	    if(flag) {
		sprintf(ifile, "png/m%03lld_lums.dat", (S64)m[i]);
	    } else {
		sprintf(ifile, "png/m%03lld_rads.dat", (S64)m[i]);
	    }
	    F64 alpha = flag ? getAlphaOfLuminosityMSPhase(&m[i]) : getAlphaOfRadiusMSPhase(&m[i]);
	    F64 beta  = flag ? getBetaOfLuminosityMSPhase(&m[i])  : getBetaOfRadiusMSPhase(&m[i]);
	    F64 gamma = flag ? 0.                                 : getGammaOfRadiusMSPhase(&m[i]);
	    F64 delta = flag ? getDifferentialLuminosityOfHook(&m[i]) : getDifferentialRadiusOfHook(&m[i]);
	    F64 qzams = flag ? getLuminosityZAMSTime(&m[i]) : getRadiusZAMSTime(&m[i]);
	    F64 qtms  = flag ? getLuminosityTMSTime(&m[i])  : getRadiusTMSTime(&m[i]);
	    F64 qHeI  = flag ? getLuminosityHeITime(&m[i])  : getRadiusHeITime(&m[i]);
	    F64 powsm = flag ?  2. :  3.;
	    F64 powmd = flag ? 20. : 10.;
	    F64 powlg = 40.;
	    F64 tHeI  = getHeITime(&m[i]);
	    F64 tTMS  = tHeI * getRatioOfTMSTimeToHeITime(&m[i]);
	    F64 dtim  = tHeI / 1024.;
	    FILE * fp = fopen(ifile, "w");
	    assert(fp);
	    for(F64 tim = 0.; tim <= tHeI; tim += dtim) {
		if(tim < tTMS) {
		    F64 tau   = tim / tTMS;
		    F64 tausm = pow(tau, powsm);
		    F64 taumd = pow(tau, powmd);
		    F64 taulg = pow(tau, powlg);
		    F64 tauone = calcTauOne(tau);
		    F64 tautwo = calcTauTwo(tau);
		    F64 tauonesm = pow(tauone, powsm);
		    F64 tautwosm = pow(tautwo, powsm);
		    F64 qtime = qzams + alpha * tau + beta * taumd + gamma * taulg
			+ (qtms - qzams - alpha - beta - gamma) * tausm
			- delta * (tauonesm - tautwosm);
		    fprintf(fp, "%+e %+e\n", tim, qtime);
		} else {
		    F64 tau     = (tim - tTMS) / (tHeI - tTMS);
		    F64 qtime   = qtms + tau * (qHeI - qtms);
		    fprintf(fp, "%+e %+e\n", tim, qtime);
		}
	    }
	    fclose(fp);
	}
    }

    void checkAGBRadius() {
	const S64 nsample = 14;
	static F64 m[nsample] = {8., 10., 13., 16., 20., 25., 32.,
				 40., 50., 65., 80., 100., 125., 160.};
	for(S64 i = 0; i < nsample; i++) {
	    if(LowerLimitOfMassAllBlue <= m[i] && m[i] < UpperLimitOfMassAllBlue) {
		continue;
	    }
	    F64 lBlue = getLuminosityEndTimeOfBluePhase(&m[i]);
	    F64 lFin  = getLuminosityEndTime(&m[i]);
	    F64 dllog = (lFin - lBlue) / 256.;
	    char ifile[1024];
	    sprintf(ifile, "png/m%03lld_ragb.dat", (S64)m[i]);
	    FILE * fp = fopen(ifile, "w");
	    assert(fp);
	    for(F64 logl = lBlue; logl < lFin; logl += dllog) {
		F64 lum = pow(10., logl);
		F64 rad = getRadiusRedPhase(&m[i], &lum);
		fprintf(fp, "%+e %+e\n", logl, rad);
	    }
	    fclose(fp);
	}
    }

};

namespace EMP = ExtremeMetalPoors;

bool askInUseOrNot() {
    using namespace EMP;
    return InUse;
}

bool askInScopeOfApplication(F64 * mass) {
    using namespace EMP;
    bool inScope = false;
    if(*mass >= LowerLimitOfMassInScopeOfApplication) {
	inScope = true;
    }
    return inScope;
}

void setLowerLimitOfMass(F64 * mass) {
    using namespace EMP;
    LowerLimitOfMassInScopeOfApplication = *mass;
    return;
}

F64 getLowerLimitOfMass() {
    return EMP::LowerLimitOfMassInScopeOfApplication;
}

F64 getUpperLimitOfMass() {
    return EMP::UpperLimitOfMassInScopeOfApplication;
}

F64 getMetallicity() {
    using namespace EMP;
#ifdef MARIGOMODEL
    return 0.;
#else
    return 0.02 * pow(10., zeta);
#endif
}

F64 getWindMetallicity() {
    using namespace EMP;
    static bool first = true;
    const  F64  czeta = -3.;
/*
#if 1
    if(first && zeta < czeta) {
	fprintf(stderr, "Caution: Wind metallicity is 10^%.f Zsun\n", czeta);
	first = false;
    }
    F64 tzeta = (zeta < czeta) ? czeta : zeta;
#else
    if(first) {
	fprintf(stderr, "Caution: Wind metallicity is 10^%.f Zsun\n", czeta);
	first = false;
    }
    F64 tzeta = czeta;
#endif
    F64 zwind =  0.02 * pow(10., tzeta);
*/
    if(first) {
	fprintf(stderr, "Caution: Wind metallicity is 10^%.f Zsun\n", zeta);
	first = false;
    }
    F64 zwind =  0.02 * pow(10., zeta);
    return zwind;
}

bool askAllBlueOrNot(F64 * mt) {
    using namespace EMP;
    if(LowerLimitOfMassAllBlue <= *mt && *mt < UpperLimitOfMassAllBlue) {
	return true;
    } else{
	return false;
    }
}

bool askBlueOrRed(F64 * aj,
		  F64 * mass) {
    F64 tBlue = getEndTimeOfBluePhase(mass);
    F64 tFin  = getEndTime(mass);
    if((*aj) < tBlue || tBlue == tFin) {
	return true;
    } else {
	return false;
    }
}

bool askBlueOrRed2(F64 * lumpersun,
		   F64 * radpersun) {
    const F64 sigm = 5.670e-5;
    const F64 lsun = 3.828e33;
    const F64 rsun = 6.963e10;
    F64 lum  = (*lumpersun) * lsun;
    F64 rad  = (*radpersun) * rsun;
    F64 rad2 = rad * rad;
    F64 logteff = 0.25*log10(lum / (4. * M_PI * sigm * rad2));
    if(logteff >= 3.65) {
	return true;
    } else {
	return false;
    }
}

bool askRadiativeOrNot(S32 * kw,
		       F64 * aj,
		       F64 * mass) {
    /*
    if(*kw == 4 || (*kw == 5 && *mass <= 50.)) { 
	return true;
    } else {
	return false;
    }
    */
    if((*kw == 4 || *kw == 5) && askBlueOrRed(aj, mass)) {
	return true;
    } else {
	return false;
    }
}

bool askRadiativeOrNot2(S32 * kw,
			F64 * lumpersun,
			F64 * radpersun) {
    if((*kw == 4 || *kw == 5) && askBlueOrRed2(lumpersun, radpersun)) {
	return true;
    } else {
	return false;
    }
}

F64 getCriticalMassRatio(S32 * kw,
			 F64 * aj,
			 F64 * mass,
			 F64 * massc) {
    F64 qc = 0.;

#ifdef KINUGAWACE
    if(*kw == 1) {
	qc = 2.;
    } else if(*kw == 2) {
	qc = 4.;
    } else if(*kw == 3 || *kw == 5 || *kw == 6) {
	qc = 0.362 + 1. / (3. * (1. - *massc / *mass));
    } else if(*kw == 4) {
	qc = 4.;
    } else if(*kw == 7) {
	qc = 1.7;
    } else if(*kw == 8 || *kw == 9) {
	qc = 3.5;
    }
#else
    if(*kw == 1) {
	qc = 2.;
    } else if(*kw == 2) {
	qc = 4.;
    } else if(*kw == 3 || *kw == 6
	      || ((*kw == 4 || *kw == 5) && (!askBlueOrRed(aj, mass)))) {
	qc = 0.362 + 1. / (3. * (1. - *massc / *mass));
    } else if((*kw == 4 || *kw == 5) && askBlueOrRed(aj, mass)) {
	qc = 4.;
    } else if(*kw == 7) {
	qc = 1.7;
    } else if(*kw == 8 || *kw == 9) {
	qc = 3.5;
    }
#endif

    return qc;
}

F64 getCriticalMassRatio2(S32 * kw,
			  F64 * lumpersun,
			  F64 * radpersun,
			  F64 * mass,
			  F64 * massc) {
    static bool FirstCall = true;
    static const bool OriginalFormula = true;
    F64 qc = 0.;
    if(OriginalFormula) {
	if(FirstCall) {
	    fprintf(stderr, "**** CriticalMassRatio is Original! ****\n");
	    FirstCall = false;
	}
	if(*kw == 1) {
	    qc = 3.;
	} else if(*kw == 2) {
	    qc = 4.;
	} else if(*kw == 3 || *kw == 6
		  || ((*kw == 4 || *kw == 5) && (!askBlueOrRed2(lumpersun, radpersun)))) {
	    qc = 0.362 + 1. / (3. * (1. - *massc / *mass));
	} else if((*kw == 4 || *kw == 5) && askBlueOrRed2(lumpersun, radpersun)) {
	    qc = 3.;
	} else if(*kw == 7) {
	    qc = 3.;
	} else if(*kw == 8 || *kw == 9) {
	    qc = 0.784;
	}
    } else {
	if(*kw == 1) {
	    qc = 2.;
	} else if(*kw == 2) {
	    qc = 4.;
	} else if(*kw == 3 || *kw == 6
		  || ((*kw == 4 || *kw == 5) && (!askBlueOrRed2(lumpersun, radpersun)))) {
	    qc = 0.362 + 1. / (3. * (1. - *massc / *mass));
	} else if((*kw == 4 || *kw == 5) && askBlueOrRed2(lumpersun, radpersun)) {
	    qc = 4.;
	} else if(*kw == 7) {
	    qc = 1.7;
	} else if(*kw == 8 || *kw == 9) {
	    qc = 3.5;
	}
    }

    return qc;
}

bool askCommonEnvelopeOrNot(S32 * kw,
			    F64 * aj,
			    F64 * mass,
			    F64 * q,
			    F64 * qc,
			    F64 * radx,
			    F64 * radc) {
    bool cce = false;

#ifdef KINUGAWACE
    if((((*kw == 3) || (*kw == 5 && *mass > 50.) || (*kw == 8) || (*kw == 9))
	&& ((*q > *qc) || (*radx <= *radc)))
       || ((*kw == 2 || *kw == 4) && (*q > *qc))
       || (*kw == 5 && *mass <= 50. && *q > 4.)) {
        cce = true;
    }
#else
    if((((*kw == 2) || ((*kw == 4 || *kw == 5) && askBlueOrRed(aj, mass))) && (*q > *qc))
       || (((*kw == 3) || (*kw == 8) || (*kw == 9)
            || ((*kw == 4 || *kw == 5) && (!askBlueOrRed(aj, mass))))
           && ((*q > *qc) || (*radx <= *radc)))) {
        cce = true;
    } else if((*kw == 4 && askBlueOrRed(aj, mass) && *q > *qc)
              || (*kw == 5 && askBlueOrRed(aj, mass) && *q > *qc)
              || (*kw == 4 && (!askBlueOrRed(aj, mass)) && (*q > *qc || *radx <= *radc))
              || (*kw == 5 && (!askBlueOrRed(aj, mass)) && (*q > *qc || *radx <= *radc))) {
        assert(NULL);
    }
#endif
    return cce;
}

bool askCommonEnvelopeOrNot2(S32 * kw,
			     F64 * lumpersun,
			     F64 * radpersun,
			     F64 * mass,
			     F64 * q,
			     F64 * qc,
			     F64 * radx,
			     F64 * radc) {
    bool cce = false;

    if((((*kw == 2) || ((*kw == 4 || *kw == 5) && askBlueOrRed2(lumpersun, radpersun))) && (*q > *qc))
       || (((*kw == 3) || (*kw == 8) || (*kw == 9)
            || ((*kw == 4 || *kw == 5) && (!askBlueOrRed2(lumpersun, radpersun))))
           && ((*q > *qc) || (*radx <= *radc)))) {
        cce = true;
    } else if((*kw == 4 && askBlueOrRed2(lumpersun, radpersun) && *q > *qc)
              || (*kw == 5 && askBlueOrRed2(lumpersun, radpersun) && *q > *qc)
              || (*kw == 4 && (!askBlueOrRed2(lumpersun, radpersun)) && (*q > *qc || *radx <= *radc))
              || (*kw == 5 && (!askBlueOrRed2(lumpersun, radpersun)) && (*q > *qc || *radx <= *radc))) {
        assert(NULL);
    }
    return cce;
}

void setMetallicity(F64 * _zeta) {
    using namespace EMP;
    assert(sizeof(F32) == 4);
    assert(sizeof(F64) == 8);
    assert(sizeof(S32) == 4);
    assert(sizeof(S64) == 8);
    if(InUse) {
	return;
    }
    InUse = true;
    zeta  = *_zeta;
#ifdef MARIGOMODEL
    char ifile[1024];
    sprintf(ifile, "./ffmarigo/metal1e0.dat");
    fprintf(stderr, " Use Marigo model\n");    
#else
#ifdef GENEVAMODEL
    char idir[1024] = "ffgeneva";
    fprintf(stderr, " Use Geneva model log(Z/Zsun) = %.1f\n", zeta);
    CriticalMassMassive = 1e9;
    UpperLimitOfMassInScopeOfApplication = 1e3;
#else
    char idir[1024] = "ffbonn";
    fprintf(stderr, " Use Bonn model log(Z/Zsun) = %.1f\n", zeta);
#endif
    char ifile[1024], afile[1024];
    if(zeta == -8.) {
	sprintf(ifile, "./%s/metal1e-8.dat", idir);
	sprintf(afile, "./%s/additional1e-8.dat", idir);
    } else if(zeta == -6.) {
	sprintf(ifile, "./%s/metal1e-6.dat", idir);
	sprintf(afile, "./%s/additional1e-6.dat", idir);
	CriticalMassMassive = 1e9;
	UpperLimitOfMassInScopeOfApplication = 1e3;
    } else if(zeta == -5.) {
	sprintf(ifile, "./%s/metal1e-5.dat", idir);
	sprintf(afile, "./%s/additional1e-5.dat", idir);
	CriticalMassMassive = 1e9;
	UpperLimitOfMassInScopeOfApplication = 1e3;
    } else if(zeta == -4.) {
	sprintf(ifile, "./%s/metal1e-4.dat", idir);
	sprintf(afile, "./%s/additional1e-4.dat", idir);
	CriticalMassMassive = 1e9;
	UpperLimitOfMassInScopeOfApplication = 1e3;
    } else if(zeta == -2.) {
	sprintf(ifile, "./%s/metal1e-2.dat", idir);
	sprintf(afile, "./%s/additional1e-2.dat", idir);
    } else if(zeta == -1.) {
	sprintf(ifile, "./%s/metal1e-1.dat", idir);
	sprintf(afile, "./%s/additional1e-1.dat", idir);
	LowerLimitOfMassInScopeOfApplication = 16.;
    } else {
#if 0
	fprintf(stderr, "Error: Not yet implemented Z/Zsun = %+e\n", zeta);
	exit(0);
#else
	if(-2 < zeta && zeta < -1) {
	    EMP::setMetallicityInterpolation(zeta);
	    LowerLimitOfMassInScopeOfApplication = 16.;
#ifdef GENEVAMODEL
	    CriticalMassMassive = 1e9;
	    UpperLimitOfMassInScopeOfApplication = 1e3;
#endif
	    return;
	} else {
	    fprintf(stderr, "Error: Not yet implemented Z/Zsun = %+e\n", zeta);
	    exit(0);
	}
#endif
    }
#endif
    FILE * fp = fopen(ifile, "r");
    assert(fp);
    fscanf(fp, "%lf", &UpperLimitOfMassHG);
    fscanf(fp, "%lf", &UpperLimitOfMassBlueLoop);
    fscanf(fp, "%lf", &LowerLimitOfMassAllBlue);
    fscanf(fp, "%lf", &UpperLimitOfMassAllBlue);
    fscanf(fp, "%lf", &UpperLimitOfMassAllBlueCHeB);
    for(S64 i = 0; i < NumberOfCoefficientForTime; i++) {
	fscanf(fp, "%lf%lf%lf%lf", &atime[i], &btime[i], &ctim2[i], &dtime[i]);
    }
    for(S64 i = 0; i < NumberOfCoefficientForLuminosity; i++) {
	fscanf(fp, "%lf%lf%lf%lf", &alums[i], &blums[i], &clums[i], &dlums[i]);
    }
    for(S64 i = 0; i < NumberOfCoefficientForRadius; i++) {
	fscanf(fp, "%lf%lf%lf%lf", &arads[i], &brads[i], &crads[i], &drads[i]);
    }
    for(S64 i = 0; i < NumberOfCoefficientForCoreMass; i++) {
	fscanf(fp, "%lf%lf%lf%lf", &acore[i], &bcore[i], &ccore[i], &dcore[i]);
    }
    for(S64 i = 0; i < NumberOfCoefficientForMSLuminosity; i++) {
	fscanf(fp, "%lf%lf%lf%lf", &msalums[i], &msblums[i], &msclums[i], &msdlums[i]);
    }
    for(S64 i = 0; i < NumberOfCoefficientForMSRadius; i++) {
	fscanf(fp, "%lf%lf%lf%lf", &msarads[i], &msbrads[i], &mscrads[i], &msdrads[i]);
    }
    for(S64 i = 0; i < NumberOfCoefficientForAGBRadius; i++) {
	fscanf(fp, "%lf%lf", &agbarads[i], &agbbrads[i]);
    }
#ifndef GENEVAMODEL
    if(zeta == -8. || zeta == -2. || zeta == -1.) {
	for(S64 i = 0; i < NumberOfCoefficientForTimeMassive; i++) {
	    fscanf(fp, "%lf%lf%lf%lf", &atimeMassive[i], &btimeMassive[i], &ctimeMassive[i], &dtimeMassive[i]);
	}
	for(S64 i = 0; i < NumberOfCoefficientForLuminosityMassive; i++) {
	    fscanf(fp, "%lf%lf%lf%lf", &alumsMassive[i], &blumsMassive[i], &clumsMassive[i], &dlumsMassive[i]);
	}	
	for(S64 i = 0; i < NumberOfCoefficientForRadiusMassive; i++) {
	    fscanf(fp, "%lf%lf%lf%lf", &aradsMassive[i], &bradsMassive[i], &cradsMassive[i], &dradsMassive[i]);
	}	
	for(S64 i = 0; i < NumberOfCoefficientForCoreMassMassive; i++) {
	    fscanf(fp, "%lf%lf%lf%lf", &acoreMassive[i], &bcoreMassive[i], &ccoreMassive[i], &dcoreMassive[i]);
	}
    }
#endif
    fclose(fp);
#ifdef K2_YOSHIDA
    fp = fopen(afile, "r");
    if(fp != NULL) {
	fscanf(fp, "%lf%lf", &ak2str[0], &bk2str[0]);
	fscanf(fp, "%lf%lf", &ak2str[1], &bk2str[1]);
	fclose(fp);
	fprintf(stderr, " Use additional file\n");
	UseK2strOrNot = true;
    }
#endif
}

F64 getRatioOfTMSTimeToHeITime(F64 * mass) {
    using namespace EMP;
    F64 frac = 0.;
    if(*mass < UpperLimitOfMassHG) {
	frac = 0.99;
    } else {
	frac = 1.;
    }
    return frac;
}

F64 getHeITime(F64 * mass) {
    using namespace EMP;
    F64 tHeI = getTimeAndTimeInterval(*mass, 0);
    F64 tFin = getTimeAndTimeInterval(*mass, 3);
    if(tHeI > tFin) {
	fprintf(stderr, "Error: tHeI is larger than tFin.\n");
	fprintf(stderr, "Error: mass= %e tHeI= %e tFin= %e\n", *mass, tHeI, tFin);
	exit(0);
    }
    return tHeI;
}

F64 getTimeIntervalOfCHeBPhase(F64 * mass) {
    using namespace EMP;
    F64 tHeI = getTimeAndTimeInterval(*mass, 0);
    F64 tHef = getTimeAndTimeInterval(*mass, 1);
    F64 tFin = getTimeAndTimeInterval(*mass, 3);
    if(tHeI + tHef >= tFin) {
#if not defined MARIGOMODEL
	/*
	fprintf(stderr, "Error: tBAGB is larger than tFin.\n");
	fprintf(stderr, "Error: Mass  is %+e.\n", *mass);
	fprintf(stderr, "Error: tBAGB is %+e.\n", tHeI + tHef);
	fprintf(stderr, "Error: tFin  is %+e.\n", tFin);
	exit(0);
	*/
	tHef = (tFin - tHeI) * 0.99;
#else
	tHef = (tFin - tHeI) * 0.99;
#endif
    }
#ifdef MARIGOMODEL
    {
	F64 logtHe = 0.;
	if((*mass) <= 12.0) {
	    logtHe = 57.5925e0 - 15.1429e0 * (*mass) + 1.45846e0 * (*mass) * (*mass)
		- 0.0463627e0 * (*mass) * (*mass) * (*mass);
	} else {
	    logtHe = 6.043e0 - 0.0331597e0 * (*mass) + 0.000516053e0 * (*mass) * (*mass)
		- 2.988e-6 * (*mass) * (*mass) * (*mass);
	}
	tHef = pow(10., logtHe) * 1e-6;
    }
#endif
    return tHef;
}

F64 getEndTimeOfBluePhase(F64 * mass) {
    using namespace EMP;
#if 0 // 200205
    F64 tBAGB = getTimeAndTimeInterval(*mass, 0) + getTimeAndTimeInterval(*mass, 1);
#else
    F64 ttmp1 = getTimeAndTimeInterval(*mass, 1) / getTimeAndTimeInterval(*mass, 0);
    F64 ttmp2 = ttmp1 * getTimeAndTimeInterval(*mass, 0);
    F64 tBAGB = getTimeAndTimeInterval(*mass, 0) + ttmp2;
#endif
    F64 tBlue = getTimeAndTimeInterval(*mass, 2);
    F64 tFin  = getTimeAndTimeInterval(*mass, 3);
    tBlue = getQuantityEndTimeOfBluePhase(*mass, tBlue, tBAGB, tFin);
    if(tBlue >= tFin
       && (*mass < LowerLimitOfMassAllBlue || UpperLimitOfMassAllBlue <= *mass)) {
	//assert(tBAGB < tFin);
	tBlue = tBAGB + (tFin - tBAGB) * 0.99;
    }
    return tBlue;
}

F64 getEndTime(F64 * mass) {
    using namespace EMP;
    return getTimeAndTimeInterval(*mass, 3);
}

F64 getLuminosityZAMSTime(F64 * mass) {
    using namespace EMP;
    return getLuminosity(*mass, 0);
}

F64 getLuminosityTMSTime(F64 * mass) {
    using namespace EMP;
    F64 lTMS = getLuminosity(*mass, 1);
    F64 lHeI = getLuminosity(*mass, 2);
    if(*mass >= UpperLimitOfMassHG) {
	lTMS = lHeI;
    }
    return lTMS;
}

F64 getLuminosityHeITime(F64 * mass) {
    using namespace EMP;
    return getLuminosity(*mass, 2);
}

F64 getLuminosityBAGBTime(F64 * mass) {
    using namespace EMP;
    return getLuminosity(*mass, 3);
}

F64 getLuminosityEndTimeOfBluePhase(F64 * mass) {
    using namespace EMP;

    F64 lBAGB = getLuminosity(*mass, 3);
    F64 lBlue = getLuminosity(*mass, 4);
    F64 lFin  = getLuminosity(*mass, 5);
    lBlue = getQuantityEndTimeOfBluePhase(*mass, lBlue, lBAGB, lFin);
    if(LowerLimitOfMassAllBlue <= *mass && *mass < UpperLimitOfMassAllBlue)  {
	lBlue = getLuminosity(*mass, 4);
    }
    return lBlue;
}

F64 getLuminosityEndTime(F64 * mass) {
    using namespace EMP;
#if 0
    return getLuminosity(*mass, 5);
#else
    F64 lFin;
    if(*mass < LowerLimitOfMassAllBlue) {
	lFin = getLuminosity(*mass, 5);
    } else {
	lFin = getLuminosity(*mass, 6);
    }
    return lFin;
#endif
}

F64 getRadiusZAMSTime(F64 * mass) {
    using namespace EMP;
    return getRadius(*mass, 0);
}

F64 getRadiusTMSTime(F64 * mass) {
    using namespace EMP;
    F64 rTMS = getRadius(*mass, 1);
    F64 rHeI = getRadius(*mass, 2);
    if(*mass >= UpperLimitOfMassHG) {
	rTMS = rHeI;
    }
    // Patch for no red MS 20/09/30
    /*
    if(*mass > CriticalMassMassive) {
	F64 lTMS    = pow(10., getLuminosityTMSTime(mass));
	F64 rTMSred = getRadiusRedPhase(mass, &lTMS);
	//F64 rTMStmp = log10(0.9 * pow(10., rTMSred));
	F64 rTMStmp = log10(0.1 * pow(10., rTMSred));
	rTMS = (rTMS <= rTMStmp) ? rTMS : rTMStmp;
    }
    */
    return rTMS;
}

F64 getRadiusHeITime(F64 * mass) {
    using namespace EMP;
    return getRadius(*mass, 2);
}

F64 getMinimumRadiusCHeBPhase(F64 * mass) {
    using namespace EMP;
    F64 rMin = getRadius(*mass, 3);
    F64 rHeI = getRadius(*mass, 2);
    //if(*mass >= UpperLimitOfMassBlueLoop) {
    if(*mass >= UpperLimitOfMassBlueLoop || rMin > rHeI) {
	rMin = getRadius(*mass, 2);
    }
    assert(rMin <= rHeI);
    return rMin;
}

F64 getRadiusRedPhase(F64 * mt,
		      F64 * lum) {
    using namespace EMP;
    F64 a   = getAGBRadius(*mt, 0);
    F64 b   = getAGBRadius(*mt, 1);
    return a + b * log10(*lum);
}

F64 getRadiusEndTimeOfBlueCHeBPhase(F64 * mass,
				    F64 * mt,
				    F64 * lum) {
    using namespace EMP;
    F64 rBlueCHeB = getRadius(*mass, 4);
    if(*mass >= UpperLimitOfMassAllBlueCHeB) {
#ifdef TESTMODE
	rBlueCHeB = 0.;
#else
	rBlueCHeB = getRadiusRedPhase(mt, lum);
#endif
    }
#ifdef MARIGOMODEL
    F64 m1 = *mass;
    F64 m2 = m1 * m1;
    F64 m3 = m2 * m1;
    F64 m4 = m3 * m1;
    if(m1 < 70.) {
	rBlueCHeB = -0.0723005e0 + 0.0814329e0  * m1 - 0.00252995e0 * m2
	    + 5.88465e-5 * m3 - 4.28501e-7 * m4;
    } else {
	rBlueCHeB =  2.93694e0   + 0.00271667e0 * m1;
    }    
#endif
    return rBlueCHeB;
}

F64 getRadiusEndTime(F64 * mass) {
    using namespace EMP;
    F64 rFin = getRadius(*mass, 5);
    if(*mass < LowerLimitOfMassAllBlue || UpperLimitOfMassAllBlue <= *mass) {
	rFin = 0.;
    }
    return rFin;
}

#if not defined MARIGOMODEL
F64 getHeCoreMassHeITime(F64 * mass) {
    using namespace EMP;
    return getCoreMass(*mass, 0);
}

F64 getHeCoreMassBAGBTime(F64 * mass) {
    using namespace EMP;
    return getCoreMass(*mass, 1);
}

F64 getCOCoreMassEndTime(F64 * mass) {
    using namespace EMP;
    return getCoreMass(*mass, 2);
}
#else
F64 getHeCoreMassHeITime(F64 * mass) {
    F64 mc = 0.;
    if(*mass < 15.) {
	mc = -0.47466e0 + 0.184148e0 * pow(*mass, 1.13274e0);
    } else {
	mc = -2.3546e0  + 0.271582e0 * pow(*mass, 1.12392e0);
    }
    return mc;
}

F64 getHeCoreMassBAGBTime(F64 * mass) {
    F64 m1 = *mass;
    F64 m2 = m1 * m1;
    F64 m3 = m2 * m1;
    F64 m4 = m3 * m1;
    F64 mc = -0.00584516e0 + 0.131569e0 * m1 + 0.00993475e0 * m2 - 0.000112405e0* m3
	+ 4.60669e-7 * m4;
    return mc;
}

F64 getCOCoreMassEndTime(F64 * mass) {
    //F64 m1 = *mass;
    F64 m1 = *mass <= 100. ? *mass : 100.;
    F64 m2 = m1 * m1;
    F64 m3 = m2 * m1;
    F64 m4 = m3 * m1;
    F64 m5 = m4 * m1;
    F64 mc = 0.618397e0 - 0.057395e0 * m1 + 0.0173053e0 * m2 - 0.000312008e0 * m3
	+ 2.99858e-6 * m4 - 1.12942e-8*m5;
    return mc;
}
#endif

F64 getAlphaOfLuminosityMSPhase(F64 * mass) {
    using namespace EMP;
    return getMSLuminosity(*mass, 0);
}

F64 getBetaOfLuminosityMSPhase(F64 * mass) {
    using namespace EMP;
    return getMSLuminosity(*mass, 1);
}

F64 getDifferentialLuminosityOfHook(F64 * mass) {
    using namespace EMP;
    return getMSLuminosity(*mass, 2);
}

F64 getAlphaOfRadiusMSPhase(F64 * mass) {
    using namespace EMP;
    return getMSRadius(*mass, 0);
}

F64 getBetaOfRadiusMSPhase(F64 * mass) {
    using namespace EMP;
    return getMSRadius(*mass, 1);
}

F64 getGammaOfRadiusMSPhase(F64 * mass) {
    using namespace EMP;
    return getMSRadius(*mass, 2);
}

F64 getDifferentialRadiusOfHook(F64 * mass) {
    using namespace EMP;
    return getMSRadius(*mass, 3);
}

F64 getLuminosityMSPhaseMassive(F64 * mass,
				F64 * tau) {
    using namespace EMP;
    F64 lzams  = pow(10., getLuminosityZAMSTime(mass));
    F64 ltms   = pow(10., getLuminosityTMSTime(mass));
    F64 lum    = lzams - (lzams - ltms) * *tau;
    return log10(lum);
}

F64 getLuminosityMSPhaseIntermediate(F64 * mass,
				     F64 * tau,
				     F64 * tau1,
				     F64 * tau2,
				     F64 * taumin,
				     F64 * eta) {
    using namespace EMP;
    F64 mcrit1 = CriticalMassMassive;
    F64 mcrit2 = mcrit1 * 2;
    F64 logl2  = getLuminosityMSPhaseMassive(&mcrit2, tau);
    F64 dell   = getDifferentialLuminosityOfHook(&mcrit1);
    F64 dtau   = pow(*tau1, 2.) - pow(*tau2, 2.);
    F64 alpha  = getAlphaOfLuminosityMSPhase(&mcrit1);
    F64 beta   = getBetaOfLuminosityMSPhase(&mcrit1);
    F64 lzams  = pow(10., getLuminosityZAMSTime(&mcrit1));
    F64 ltms   = pow(10., getLuminosityTMSTime(&mcrit1));
    F64 lx     = log10(ltms / lzams);
    F64 xx     = 0.;
    if(*tau > *taumin) {
	xx = alpha * *tau + beta * pow(*tau, *eta)
	    + (lx - alpha - beta) * pow(*tau, 2.) - dell * dtau;
    } else {
	xx = alpha * *tau + (lx - alpha) * pow(*tau, 2.) - dell * dtau;
    }
    F64 logl1  = log10(lzams * pow(10., xx));
    F64 loglx  = (logl2 - logl1) * (*mass - mcrit1) / (mcrit2 - mcrit1) + logl1;
    /////// keep continuity from type 1 to type 4 (21/09/26)
    F64 logly  = getLuminosityMSPhaseMassive(mass, tau);
    loglx = (1. - (*tau)) * loglx + (*tau) * logly;
    ///////
    return loglx;
}

F64 getRadiusMSPhaseMassive(F64 * mass,
			    F64 * tau) {
    using namespace EMP;
    F64 rzams  = getRadiusZAMSTime(mass);
    F64 rtms   = getRadiusTMSTime(mass);
    F64 rad    = rzams + (rtms - rzams) * (*tau) * (*tau) * (*tau) * (*tau) * (*tau);
    return rad;
}

F64 getRadiusMSPhaseIntermediate(F64 * mass,
				 F64 * tau,
				 F64 * tau1,
				 F64 * tau2,
				 F64 * taumin) {
    using namespace EMP;
    F64 mcrit1 = CriticalMassMassive;
    F64 mcrit2 = mcrit1 * 2;
    F64 logr2  = getRadiusMSPhaseMassive(&mcrit2, tau);
    F64 delr   = getDifferentialRadiusOfHook(&mcrit1);
    F64 dtau   = pow(*tau1, 3.) - pow(*tau2, 3.);
    F64 alpha  = getAlphaOfRadiusMSPhase(&mcrit1);
    F64 beta   = getBetaOfRadiusMSPhase(&mcrit1);
    F64 gamma  = getGammaOfRadiusMSPhase(&mcrit1);
    F64 rzams  = pow(10., getRadiusZAMSTime(&mcrit1));
    F64 rtms   = pow(10., getRadiusTMSTime(&mcrit1));
    F64 rx     = log10(rtms / rzams);
    F64 xx     = 0.;
    if(*tau > *taumin) {
	xx = alpha * *tau + beta * pow(*tau, 10.) + gamma * pow(*tau, 40.)
	    + (rx - alpha - beta - gamma) * pow(*tau, 3.) - delr * dtau;
    } else {
	xx = alpha * *tau + (rx - alpha) * pow(*tau, 3.) - delr * dtau;
    }
    F64 logr1  = log10(rzams * pow(10., xx));
    F64 logrx  = (logr2 - logr1) * (*mass - mcrit1) / (mcrit2 - mcrit1) + logr1;
    /////// keep continuity from type 1 to type 4 (21/09/26)
    F64 logry  = getRadiusMSPhaseMassive(mass, tau);
    logrx = (1. - (*tau)) * logrx + (*tau) * logry;
    ///////
    return logrx;
}

F64 getK2strOfBluePhase(F64 * mt,
			F64 * mc,
			F64 * rzams,
			F64 * rtms,
			F64 * rad,
			F64 * k2e) {
    using namespace EMP;
    if(UseK2strOrNot && *mt >= LowerLimitOfMassInScopeOfApplication) {
	F64 lmass  = log10(*mt);
	F64 power1 = ak2str[0] + bk2str[0] * lmass;
	if(*mc == 0.) {
	    *k2e = 0.1 * pow((*rad / *rzams), power1);
	} else {
	    F64 k2etms = 0.1 * pow((*rtms / *rzams), power1);
	    F64 power2 = ak2str[1] + bk2str[1] * lmass;
	    *k2e = k2etms * pow((*rad / *rtms), power2);
	}
    }
    return *k2e;
}

F64 getConvectiveCoreRadiusOfBluePhase(F64 * mt) {
    using namespace EMP;
    F64 rconv = pow(10., 0.05 * (zeta - 1.)) * pow((*mt / 10.), 0.8);
    return rconv;
}

#ifdef TESTMODE
int main(int argc,
	 char ** argv) {

    if (argc != 3) {
	fprintf(stderr, "%s <Z/Zsun> <mode[=time/lums/rads/core/MSHGLuminosity/MSHGRadius/AGBRadius]>\n", argv[0]);
	exit(0);
    }

    F64 zeta = atof(argv[1]);

    setMetallicity(&zeta);

    if(strcmp(argv[2], "MSHGLuminosity") == 0 || strcmp(argv[2], "MSHGRadius") == 0) {
	EMP::checkMSHGPhase(argv[2]);
	exit(0);
    }

    if(strcmp(argv[2], "AGBRadius") == 0) {
	EMP::checkAGBRadius();
	exit(0);
    }

    for(F64 m = 8.; m < 200.; m *= 1.03) {
#if not defined MARIGOMODEL
	EMP::checkTimeSequence(m);
	EMP::checkLuminosity(m);
	EMP::checkRadius(m);
#endif
	if(strcmp(argv[2], "time") == 0) {
	    printf("%+e %+e %+e %+e %+e\n", m, getHeITime(&m), getTimeIntervalOfCHeBPhase(&m),
		   getEndTimeOfBluePhase(&m), getEndTime(&m));
	} else if(strcmp(argv[2], "lums") == 0) {
	    printf("%+e %+e %+e %+e %+e %+e %+e\n", m, getLuminosityZAMSTime(&m), getLuminosityTMSTime(&m),
		   getLuminosityHeITime(&m), getLuminosityBAGBTime(&m),
		   getLuminosityEndTimeOfBluePhase(&m), getLuminosityEndTime(&m));
	} else if(strcmp(argv[2], "rads") == 0) {
	    printf("%+e %+e %+e %+e %+e %+e %+e\n", m, getRadiusZAMSTime(&m), getRadiusTMSTime(&m),
		   getRadiusHeITime(&m), getMinimumRadiusCHeBPhase(&m),
		   getRadiusEndTimeOfBlueCHeBPhase(&m, &m, &m), getRadiusEndTime(&m));
	} else if(strcmp(argv[2], "core") == 0) {
	    printf("%+e %+e %+e %+e\n", m, getHeCoreMassHeITime(&m),
		   getHeCoreMassBAGBTime(&m), getCOCoreMassEndTime(&m));
	} else {
	    fprintf(stderr, "Error: Mode %s is not found.\n", argv[2]);
	    exit(0);
	}
    }

    return 0;
}
#endif
