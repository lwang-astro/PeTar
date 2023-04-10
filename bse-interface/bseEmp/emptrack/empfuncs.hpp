#if defined(__cplusplus)
extern "C" {
#endif
    void setMetallicity(F64 * _zeta);

    bool askAllBlueOrNot(F64 * mt);
    bool askInUseOrNot();
    bool askInScopeOfApplication(F64 * mass);
    void setLowerLimitOfMass(F64 * mass);
    F64 getLowerLimitOfMass();
    F64 getUpperLimitOfMass();
    F64 getMetallicity();
    F64 getWindMetallicity();
    F64 getCriticalMassMassive();
    bool askBlueOrRed(F64 * aj,
		      F64 * mass);
    bool askBlueOrRed2(F64 * lumpersun,
		       F64 * radpersun);
    bool askBlueOrRed3(F64 * lumpersun,
		       F64 * radpersun,
		       F64 * mt);
    bool askRadiativeOrNot(S32 * kw,
			   F64 * aj,
			   F64 * mass);
    bool askRadiativeOrNot2(S32 * kw,
			    F64 * lumpersun,
			    F64 * radpersun);
    bool askRadiativeOrNot3(S32 * kw,
			    F64 * lumpersun,
			    F64 * radpersun,
			    F64 * mt);
    F64 getCriticalMassRatio(S32 * kw,
			     F64 * aj,
			     F64 * mass,
			     F64 * massc);
    F64 getCriticalMassRatio2(S32 * kw,
			      F64 * lumpersun,
			      F64 * radpersun,
			      F64 * mass,
			      F64 * massc,
			      S32 * preventCe);
    F64 getCriticalMassRatio3(S32 * kw,
			      F64 * lumpersun,
			      F64 * radpersun,
			      F64 * mass,
			      F64 * massc,
			      S32 * preventCe);
    bool askCommonEnvelopeOrNot(S32 * kw,
				F64 * aj,
				F64 * mass,
				F64 * q,
				F64 * qc,
				F64 * radx,
				F64 * radc);
    bool askCommonEnvelopeOrNot2(S32 * kw,
				 F64 * lumpersun,
				 F64 * radpersun,
				 F64 * mass,
				 F64 * q,
				 F64 * qc,
				 F64 * radx,
				 F64 * radc);
    bool askCommonEnvelopeOrNot3(S32 * kw,
				 F64 * lumpersun,
				 F64 * radpersun,
				 F64 * mass,
				 F64 * q,
				 F64 * qc,
				 F64 * radx,
				 F64 * radc);

    F64 getRatioOfTMSTimeToHeITime(F64 * mass);
    F64 getHeITime(F64 * mass);
    F64 getTimeIntervalOfCHeBPhase(F64 * mass);
    F64 getEndTimeOfBluePhase(F64 * mass);
    F64 getEndTime(F64 * mass);
    
    F64 getLuminosityZAMSTime(F64 * mass);
    F64 getLuminosityTMSTime(F64 * mass);
    F64 getDifferentialLuminosityOfHook(F64 * mass);
    F64 getAlphaOfLuminosityMSPhase(F64 * mass);
    F64 getBetaOfLuminosityMSPhase(F64 * mass);
    F64 getLuminosityHeITime(F64 * mass);
    F64 getLuminosityBAGBTime(F64 * mass);
    F64 getLuminosityEndTimeOfBluePhase(F64 * mass);
    F64 getLuminosityEndTime(F64 * mass);
    
    F64 getRadiusZAMSTime(F64 * mass);
    F64 getRadiusTMSTime(F64 * mass);
    F64 getDifferentialRadiusOfHook(F64 * mass);
    F64 getAlphaOfRadiusMSPhase(F64 * mass);
    F64 getBetaOfRadiusMSPhase(F64 * mass);
    F64 getGammaOfRadiusMSPhase(F64 * mass);
    F64 getRadiusHeITime(F64 * mass);
    F64 getMinimumRadiusCHeBPhase(F64 * mass);
    F64 getRadiusRedPhase(F64 * mt,
			  F64 * lum);
    F64 getRadiusEndTimeOfBlueCHeBPhase(F64 * mass,
					F64 * mt,
					F64 * lum);
    F64 getRadiusEndTime(F64 * mass);
    
    F64 getHeCoreMassHeITime(F64 * mass);
    F64 getHeCoreMassBAGBTime(F64 * mass);
    F64 getCOCoreMassEndTime(F64 * mass);

    F64 getTMSMassFromHeCoreMassHeITime(F64 * mc);
    F64 getTMSMassFromHeCoreMassBAGBTime(F64 * mc);

    F64 getLuminosityMSPhaseMassive(F64 * mass,
				    F64 * tau);
    F64 getLuminosityMSPhaseIntermediate(F64 * mass,
					 F64 * tau,
					 F64 * tau1,
					 F64 * tau2,
					 F64 * taumin,
					 F64 * eta);
    F64 getRadiusMSPhaseMassive(F64 * mass,
				F64 * tau);
    F64 getRadiusMSPhaseIntermediate(F64 * mass,
				     F64 * tau,
				     F64 * tau1,
				     F64 * tau2,
				     F64 * taumin);
    F64 getK2strOfBluePhase(F64 * mt,
			    F64 * mc,
			    F64 * rzams,
			    F64 * rtms,
			    F64 * rad,
			    F64 * k2e);

    F64 getConvectiveCoreRadiusOfBluePhase(F64 * mt);

    void setMetallicityInZUsingInBSE(F64 * zbse);
    F64 getMetallicityInZUsingInBSE();

#if defined(__cplusplus)
}
#endif
