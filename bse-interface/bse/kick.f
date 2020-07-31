***
      SUBROUTINE kick(kw,m1,m1n,m2,ecc,sep,jorb,v,fbfac,fbtot,mco,ecs)
      implicit none
*
      integer kw,k
      INTEGER idum
      COMMON /VALUE3/ idum
      INTEGER idum2,iy,ir(32)
      COMMON /RAND3/ idum2,iy,ir
      integer bhflag,mxns
      real*8 m1,m2,m1n,ecc,sep,jorb,ecc2
      real*8 pi,twopi,gmrkm,yearsc,rsunkm
      parameter(yearsc=3.1557d+07,rsunkm=6.96d+05)
      real*8 mm,em,dif,der,del,r
      real*8 u1,u2,vk,v(4),s,theta,phi
      real*8 sphi,cphi,stheta,ctheta,salpha,calpha
      real*8 vr,vr2,vk2,vn2,hn2
      real*8 mu,cmu,vs(3),v1,v2,mx1,mx2
      real*8 FBFAC,FBTOT,MCO
      integer ECS
      REAL*8 CONVF,MNUEFF,MNSTYP
      INTEGER psflag,kmech,ecflag
      real*8 sigma
      real*8 DISP0,ECSIG,WDSIG1,WDSIG2,WDKMAX
      COMMON /VALUE4/ sigma,mxns,bhflag
      real ran3,xx
      external ran3
      COMMON /FLAGS2/ psflag,kmech,ecflag
*
*       Choose the kick settings.
*       Adopt the Maxwellian of 'sigma' (input file) 
*       for regular NSs/BHs (no fallback), with EC kicks from a Maxwellian with a lower 
*       peak (ECSIG > 0, ECSIG = 0.0 => EC kicks are zero) and BH/NS kicks scaled by fallback.
*       WD kicks are lowered (if WDSIG1,WDSIG2,WDKMAX > 0, otherwise WD kicks are zero).
*
      ECSIG = 3.D0
      WDSIG1 = 2.D0
      WDSIG2 = 2.D0
      WDKMAX = 6.D0
*******
*     Chris L. Fryer. Oct 2018.
*
* CONVF: convective boost factor larger CO core masses 
*        in the case of convection-asymmerty-driven
*        kick mechanism (typical range: 2.0-10.0)
*
* MNUEFF: in case of neutrino-driven kick mechanism, the 
*         effective remnant mass beyond which the neutrino emission does not
*         enhance significantly as the remnant (baryonic) mass is increased
*         (typical range: 5.0-10.0 Msun)
*
* MNSTYP: typical mass of a neutron star with the input dispersion velocity 'DISP'
*
* KMECH: kick mechanism. 1 = standard momentum-conserving,
*                        2 = convection-asymmetry-driven,
*                        3 = collapse-asymmerty-driven,
*                        4 = neutrino driven
*
* It is assumed that one of these four mechanisms is the primary driver
* of SN natal kicks that we observe
*
      CONVF = 5.0D0
      MNUEFF = 7.0D0
      MNSTYP = 1.4D0
*******
*
*     ECSIG = 0.D0
*     WDSIG1 = 2.D0
*     WDSIG2 = 2.D0
*     WDKMAX = 6.D0
*
      do k = 1,3
         vs(k) = 0.d0
      enddo
*     if(kw.eq.14.and.bhflag.eq.0) goto 95
*
      pi = ACOS(-1.d0)
      twopi = 2.d0*pi
* Conversion factor to ensure velocities are in km/s using mass and
* radius in solar units.
      gmrkm = 1.906125d+05
*
* Find the initial separation by randomly choosing a mean anomaly.
      if(sep.gt.0.d0.and.ecc.ge.0.d0)then
         xx = RAN3(idum)
         mm = xx*twopi
         em = mm
 2       dif = em - ecc*SIN(em) - mm
         if(ABS(dif/mm).le.1.0d-04) goto 3
         der = 1.d0 - ecc*COS(em)
         del = dif/der
         em = em - del
         goto 2
 3       continue
         r = sep*(1.d0 - ecc*COS(em))
*
* Find the initial relative velocity vector.
         salpha = SQRT((sep*sep*(1.d0-ecc*ecc))/(r*(2.d0*sep-r)))
         calpha = (-1.d0*ecc*SIN(em))/SQRT(1.d0-ecc*ecc*COS(em)*COS(em))
         vr2 = gmrkm*(m1+m2)*(2.d0/r - 1.d0/sep)
         vr = SQRT(vr2)
      else
         vr = 0.d0
         vr2 = 0.d0
         salpha = 0.d0
         calpha = 0.d0
      endif
*
* Generate Kick Velocity using Maxwellian Distribution (Phinney 1992).
      DISP0 = MAX(sigma,0.D0)
      IF(KW.EQ.10.OR.KW.EQ.11) DISP0 = MAX(WDSIG1,0.D0)
      IF(KW.EQ.12) DISP0 = MAX(WDSIG2,0.D0)
      IF(ECS.EQ.1)THEN
         IF(ECSIG.LT.-0.01)THEN
            DISP0 = DISP0*ABS(ECSIG)
         ELSE
            DISP0 = MAX(ECSIG,0.D0)
         ENDIF
*         write(6,*) "ECS-NS formation MASS KS DISP: ",
*     &   m1n, KW, DISP0
      ENDIF
* Use Henon's method for pairwise components (Douglas Heggie 22/5/97).
      do 20 k = 1,2
         u1 = RAN3(idum)
         u2 = RAN3(idum)
* Generate two velocities from polar coordinates S & THETA.
         s = DISP0*SQRT(-2.d0*LOG(1.d0 - u1))
         theta = twopi*u2
         v(2*k-1) = s*COS(theta)
         v(2*k) = s*SIN(theta)
 20   continue
      vk2 = v(1)**2 + v(2)**2 + v(3)**2
      vk = SQRT(vk2)
      v(4) = vk
* Impose the maximum WD kick velocity. 
      IF(KW.GE.10.AND.KW.LE.12.AND.vk.GT.WDKMAX)THEN
         vk = WDKMAX
         vk2 = vk*vk 
      ENDIF
* Kick modification for BH/heavy-NS according to 'bhflag'  
      if((kw.eq.14.and.bhflag.eq.0).or.kw.lt.0)then
         vk2 = 0.d0
         vk = 0.d0
         if(kw.lt.0) kw = 13
      endif
      IF (BHFLAG.GT.1) THEN
         IF (KW.EQ.13.OR.KW.EQ.14) THEN
********* Exclude ECS-NS from the treatment *********
          IF(ECS.EQ.0)THEN
********* Nutrino-driven kick *****************
            IF (KMECH.EQ.4) THEN
                vk = vk*(MIN(m1n,MNUEFF)/m1n)
            ELSE
********* Standard momentum-conserving kick *********
                vk = vk*(1.0D0 - FBFAC)
********* Covection-asymmetry-driven kick ********
            IF ((KMECH.EQ.2).AND.(MCO.LE.3.5D0))
     &          vk = vk*(MNSTYP/m1n)
            IF ((KMECH.EQ.2).AND.(MCO.GT.3.5D0))
     &          vk = vk*(MNSTYP/m1n)*CONVF
********* Collapse-asymmetry-driven kick ********
            IF ((KMECH.EQ.3).AND.(MCO.LE.3.0D0))
     &          vk = vk*(MNSTYP/m1n)
            IF ((KMECH.EQ.3).AND.(MCO.GT.3.0D0))
     &          vk = vk*(MNSTYP/m1n)*0.1D0
            ENDIF
          ENDIF
            vk2 = vk*vk
*            WRITE(6,*) "Heavy-NS/BH formation", 
*     &" (mechanism/fallback control) MASS0 MASS KS",
*     &" FBFAC FBTOT MCO VKICK KMECH: ",
*     &           m1, m1n, KW, FBFAC, FBTOT, MCO, vk, KMECH
          
         ENDIF     
      ENDIF

*     rescale base on new vk
      do k = 1,3
         v(k) = v(k)*vk/v(4)
      end do
      v(4) = vk

      sphi = -1.d0 + 2.d0*u1
      phi = ASIN(sphi)
      cphi = COS(phi)
      stheta = SIN(theta)
      ctheta = COS(theta)
*     WRITE(66,*)' KICK VK PHI THETA ',vk,phi,theta
      if(sep.le.0.d0.or.ecc.lt.0.d0) goto 90
*
* Determine the magnitude of the new relative velocity.
      vn2 = vk2+vr2-2.d0*vk*vr*(ctheta*cphi*salpha-stheta*cphi*calpha)
* Calculate the new semi-major axis.
      sep = 2.d0/r - vn2/(gmrkm*(m1n+m2))
      sep = 1.d0/sep
*     if(sep.le.0.d0)then
*        ecc = 1.1d0
*        goto 90
*     endif
* Determine the magnitude of the cross product of the separation vector
* and the new relative velocity.
      v1 = vk2*sphi*sphi
      v2 = (vk*ctheta*cphi-vr*salpha)**2
      hn2 = r*r*(v1 + v2)
* Calculate the new eccentricity.
      ecc2 = 1.d0 - hn2/(gmrkm*sep*(m1n+m2))
      ecc2 = MAX(ecc2,0.d0)
      ecc = SQRT(ecc2)
* Calculate the new orbital angular momentum taking care to convert
* hn to units of Rsun^2/yr.
      jorb = (m1n*m2/(m1n+m2))*SQRT(hn2)*(yearsc/rsunkm)
* Determine the angle between the new and old orbital angular
* momentum vectors.
      cmu = (vr*salpha-vk*ctheta*cphi)/SQRT(v1 + v2)
      mu = ACOS(cmu)
* Calculate the components of the velocity of the new centre-of-mass.
 90   continue
      if(ecc.le.1.0)then
* Calculate the components of the velocity of the new centre-of-mass.
         mx1 = vk*m1n/(m1n+m2)
         mx2 = vr*(m1-m1n)*m2/((m1n+m2)*(m1+m2))
         vs(1) = mx1*ctheta*cphi + mx2*salpha
         vs(2) = mx1*stheta*cphi + mx2*calpha
         vs(3) = mx1*sphi
      else
* Calculate the relative hyperbolic velocity at infinity (simple method).
         sep = r/(ecc-1.d0)
*        cmu = SQRT(ecc-1.d0)
*        mu = ATAN(cmu)
         mu = ACOS(1.d0/ecc)
         vr2 = gmrkm*(m1n+m2)/sep
         vr = SQRT(vr2)
         vs(1) = vr*SIN(mu)
         vs(2) = vr*COS(mu)
         vs(3) = 0.d0
         ecc = MIN(ecc,99.99d0)
      endif
*
*
      RETURN
      END
***
