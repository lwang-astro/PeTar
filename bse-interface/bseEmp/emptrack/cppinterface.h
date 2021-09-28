      interface 
      subroutine setMetallicity(zeta) bind(c, name='setMetallicity')
      import
      implicit none
      real(c_double) zeta
      end subroutine

      function askInUseOrNot() bind(c, name='askInUseOrNot')
      import
      implicit none
      logical(c_bool) askInUseOrNot
      end function
      
      function askInScopeOfApplication(mass) 
     &     bind(c, name='askInScopeOfApplication')
      import
      implicit none
      real(c_double) mass
      logical(c_bool) askInScopeOfApplication
      end function

      subroutine setLowerLimitOfMass(mass)
     &     bind(c, name='setLowerLimitOfMass')
      import
      implicit none
	  real(c_double) mass
      end subroutine
      
      function getLowerLimitOfMass() 
     &     bind(c, name='getLowerLimitOfMass')
      import
      implicit none
      real(c_double) getLowerLimitOfMass
      end function

      function getUpperLimitOfMass() 
     &     bind(c, name='getUpperLimitOfMass')
      import
      implicit none
      real(c_double) getUpperLimitOfMass
      end function

      function getMetallicity() bind(c, name='getMetallicity')
      import
      implicit none
      real(c_double) getMetallicity
      end function

      function getWindMetallicity() bind(c, name='getWindMetallicity')
      import
      implicit none
      real(c_double) getWindMetallicity
      end function

      function askAllBlueOrNot(mt) bind(c, name='askAllBlueOrNot')
      import
      implicit none
      real(c_double) mt
      logical(c_bool) askAllBlueOrNot
      end function
      
      function askBlueOrRed(aj,mass) bind(c, name='askBlueOrRed')
      import
      implicit none
      real(c_double) aj
      real(c_double) mass
      logical(c_bool) askBlueOrRed
      end function
      
      function askBlueOrRed2(lumpersun,radpersun)
     &     bind(c, name='askBlueOrRed2')
      import
      implicit none
      real(c_double) lumpersun
      real(c_double) radpersun
      logical(c_bool) askBlueOrRed2
      end function
      
      function askRadiativeOrNot(kw,aj,mass)
     & bind(c, name='askRadiativeOrNot')
      import
      implicit none
      integer(c_int) kw
      real(c_double) aj
      real(c_double) mass
      logical(c_bool) askRadiativeOrNot
      end function
      
      function askRadiativeOrNot2(kw,lumpersun,radpersun)
     & bind(c, name='askRadiativeOrNot2')
      import
      implicit none
      integer(c_int) kw
      real(c_double) lumpersun
      real(c_double) radpersun
      logical(c_bool) askRadiativeOrNot2
      end function

      function getCriticalMassRatio(kw,aj,mass,massc) 
     & bind(c, name='getCriticalMassRatio')
      import
      implicit none
      integer(c_int) kw
      real(c_double) aj
      real(c_double) mass
      real(c_double) massc
      real(c_double) getCriticalMassRatio
      end function

      function getCriticalMassRatio2(kw,lumpersun,radpersun,mass,massc) 
     & bind(c, name='getCriticalMassRatio2')
      import
      implicit none
      integer(c_int) kw
      real(c_double) lumpersun
      real(c_double) radpersun
      real(c_double) mass
      real(c_double) massc
      real(c_double) getCriticalMassRatio2
      end function

      function askCommonEnvelopeOrNot(kw,aj,mass,q,qc,radx,radc) 
     &     bind(c, name='askCommonEnvelopeOrNot')
      import
      implicit none
      integer(c_int) kw
      real(c_double) aj
      real(c_double) mass
      real(c_double) q
      real(c_double) qc
      real(c_double) radx
      real(c_double) radc
      logical(c_bool) askCommonEnvelopeOrNot
      end function

      function askCommonEnvelopeOrNot2(kw,lumpersun,radpersun,
     &     mass,q,qc,radx,radc) 
     &     bind(c, name='askCommonEnvelopeOrNot2')
      import
      implicit none
      integer(c_int) kw
      real(c_double) lumpersun
      real(c_double) radpersun
      real(c_double) mass
      real(c_double) q
      real(c_double) qc
      real(c_double) radx
      real(c_double) radc
      logical(c_bool) askCommonEnvelopeOrNot2
      end function
      
      function getRatioOfTMSTimeToHeITime(mass)
     &     bind(c, name='getRatioOfTMSTimeToHeITime')
      import
      implicit none
      real(c_double) mass
      real(c_double) getRatioOfTMSTimeToHeITime
      end function

      function getHeITime(mass) bind(c, name='getHeITime')
      import
      implicit none
      real(c_double) mass
      real(c_double) getHeITime
      end function

      function getTimeIntervalOfCHeBPhase(mass)
     &     bind(c, name='getTimeIntervalOfCHeBPhase')
      import
      implicit none
      real(c_double) mass
      real(c_double) getTimeIntervalOfCHeBPhase
      end function

      function getEndTimeOfBluePhase(mass) 
     &     bind(c, name='getEndTimeOfBluePhase')
      import
      implicit none
      real(c_double) mass
      real(c_double) getEndTimeOfBluePhase
      end function

      function getEndTime(mass) 
     &     bind(c, name='getEndTime')
      import
      implicit none
      real(c_double) mass
      real(c_double) getEndTime
      end function

      function getLuminosityZAMSTime(mass)
     &     bind(c, name='getLuminosityZAMSTime')
      import
      implicit none
      real(c_double) mass
      real(c_double) getLuminosityZAMSTime
      end function

      function getLuminosityTMSTime(mass)
     &     bind(c, name='getLuminosityTMSTime')
      import
      implicit none
      real(c_double) mass
      real(c_double) getLuminosityTMSTime
      end function

      function getDifferentialLuminosityOfHook(mass)
     &     bind(c, name='getDifferentialLuminosityOfHook')
      import
      implicit none
      real(c_double) mass
      real(c_double) getDifferentialLuminosityOfHook
      end function

      function getAlphaOfLuminosityMSPhase(mass)
     &     bind(c, name='getAlphaOfLuminosityMSPhase')
      import
      implicit none
      real(c_double) mass
      real(c_double) getAlphaOfLuminosityMSPhase
      end function

      function getBetaOfLuminosityMSPhase(mass)
     &     bind(c, name='getBetaOfLuminosityMSPhase')
      import
      implicit none
      real(c_double) mass
      real(c_double) getBetaOfLuminosityMSPhase
      end function

      function getLuminosityHeITime(mass) 
     &     bind(c, name='getLuminosityHeITime')
      import
      implicit none
      real(c_double) mass
      real(c_double) getLuminosityHeITime
      end function

      function getLuminosityBAGBTime(mass) 
     &     bind(c, name='getLuminosityBAGBTime')
      import
      implicit none
      real(c_double) mass
      real(c_double) getLuminosityBAGBTime
      end function

      function getRadiusZAMSTime(mass) bind(c, name='getRadiusZAMSTime')
      import
      implicit none
      real(c_double) mass
      real(c_double) getRadiusZAMSTime
      end function

      function getRadiusTMSTime(mass) bind(c, name='getRadiusTMSTime')
      import
      implicit none
      real(c_double) mass
      real(c_double) getRadiusTMSTime
      end function

      function getDifferentialRadiusOfHook(mass)
     &     bind(c, name='getDifferentialRadiusOfHook')
      import
      implicit none
      real(c_double) mass
      real(c_double) getDifferentialRadiusOfHook
      end function

      function getAlphaOfRadiusMSPhase(mass)
     &     bind(c, name='getAlphaOfRadiusMSPhase')
      import
      implicit none
      real(c_double) mass
      real(c_double) getAlphaOfRadiusMSPhase
      end function

      function getBetaOfRadiusMSPhase(mass)
     &     bind(c, name='getBetaOfRadiusMSPhase')
      import
      implicit none
      real(c_double) mass
      real(c_double) getBetaOfRadiusMSPhase
      end function

      function getGammaOfRadiusMSPhase(mass)
     &     bind(c, name='getGammaOfRadiusMSPhase')
      import
      implicit none
      real(c_double) mass
      real(c_double) getGammaOfRadiusMSPhase
      end function

      function getRadiusHeITime(mass)
     &     bind(c, name='getRadiusHeITime')
      import
      implicit none
      real(c_double) mass
      real(c_double) getRadiusHeITime
      end function

      function getMinimumRadiusCHeBPhase(mass)
     &     bind(c, name='getMinimumRadiusCHeBPhase')
      import
      implicit none
      real(c_double) mass
      real(c_double) getMinimumRadiusCHeBPhase
      end function

      function getRadiusRedPhase(mass, lum) 
     &     bind(c, name='getRadiusRedPhase')
      import
      implicit none
      real(c_double) mass
      real(c_double) lum
      real(c_double) getRadiusRedPhase
      end function

      function getRadiusEndTimeOfBlueCHeBPhase(mass, mt, lum) 
     &     bind(c, name='getRadiusEndTimeOfBlueCHeBPhase')
      import
      implicit none
      real(c_double) mass
      real(c_double) mt
      real(c_double) lum
      real(c_double) getRadiusEndTimeOfBlueCHeBPhase
      end function

      function getHeCoreMassHeITime(mass)
     &     bind(c, name='getHeCoreMassHeITime')
      import
      implicit none
      real(c_double) mass
      real(c_double) getHeCoreMassHeITime
      end function

      function getHeCoreMassBAGBTime(mass)
     &     bind(c, name='getHeCoreMassBAGBTime')
      import
      implicit none
      real(c_double) mass
      real(c_double) getHeCoreMassBAGBTime
      end function

      function getCOCoreMassEndTime(mass)
     &     bind(c, name='getCOCoreMassEndTime')
      import
      implicit none
      real(c_double) mass
      real(c_double) getCOCoreMassEndTime
      end function

      subroutine followAGBPhase(aj, mass, mt, lum, 
     &     r, rg, mcbagb, mc, mcx, mcmax)
     &     bind(c, name='followAGBPhase')
      import
      implicit none
      real(c_double) aj
      real(c_double) mass
      real(c_double) mt
      real(c_double) lum
      real(c_double) r
      real(c_double) rg
      real(c_double) mcbagb
      real(c_double) mc
      real(c_double) mcx
      real(c_double) mcmax
      end subroutine

      subroutine calcTimestepAGBPhase(kw, aj, mass, tn, pts3,
     &     dtm, dtr)
     &     bind(c, name='calcTimestepAGBPhase')
      import
      implicit none
      integer(c_int) kw
      real(c_double) aj
      real(c_double) mass
      real(c_double) tn
      real(c_double) pts3
      real(c_double) dtm
      real(c_double) dtr
      end subroutine

      function askInUseOfConvectiveCore()
     & bind(c, name='askInUseOfConvectiveCore')
      import
      implicit none
      logical(c_bool) askInUseOfConvectiveCore
      end function

      subroutine initializeConvectiveCore(zeta)
     & bind(c, name='initializeConvectiveCore')
      import
      implicit none
      real(c_double) zeta
      end subroutine

      function getRadiusOfConvectiveCore(kw,aj,mass)
     & bind(c, name='getRadiusOfConvectiveCore')
      import
      implicit none
      integer(c_int) kw
      real(c_double) aj
      real(c_double) mass
      real(c_double) getRadiusOfConvectiveCore
      end function

      function getLuminosityMSPhaseMassive(mass,tau)
     & bind(c, name='getLuminosityMSPhaseMassive')
      import
      implicit none
      real(c_double) mass
      real(c_double) tau
      real(c_double) getLuminosityMSPhaseMassive
      end function

      function getLuminosityMSPhaseIntermediate(mass,tau,
     &     tau1,tau2,taumin,eta)
     & bind(c, name='getLuminosityMSPhaseIntermediate')
      import
      implicit none
      real(c_double) mass
      real(c_double) tau
      real(c_double) tau1
      real(c_double) tau2
      real(c_double) taumin
      real(c_double) eta
      real(c_double) getLuminosityMSPhaseIntermediate
      end function

      function getRadiusMSPhaseMassive(mass,tau)
     & bind(c, name='getRadiusMSPhaseMassive')
      import
      implicit none
      real(c_double) mass
      real(c_double) tau
      real(c_double) getRadiusMSPhaseMassive
      end function

      function getRadiusMSPhaseIntermediate(mass,tau,
     &     tau1,tau2,taumin)
     & bind(c, name='getRadiusMSPhaseIntermediate')
      import
      implicit none
      real(c_double) mass
      real(c_double) tau
      real(c_double) tau1
      real(c_double) tau2
      real(c_double) taumin
      real(c_double) getRadiusMSPhaseIntermediate
      end function

      function getK2strOfBluePhase(mt,mc,rzams,rtms,rad,k2e)
     & bind(c, name='getK2strOfBluePhase')
      import
      implicit none
      real(c_double) mt
      real(c_double) mc
      real(c_double) rzams
      real(c_double) rtms
      real(c_double) rad
      real(c_double) k2e
      real(c_double) getK2strOfBluePhase
      end function

      function getConvectiveCoreRadiusOfBluePhase(mt)
     & bind(c, name='getConvectiveCoreRadiusOfBluePhase')
      import
      implicit none
      real(c_double) mt
      real(c_double) getConvectiveCoreRadiusOfBluePhase
      end function

      end interface
      
