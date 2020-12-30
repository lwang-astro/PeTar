***
      SUBROUTINE MERGE(kstar,mass0,mass,rad,massc,radc,menv,renv,ospin,
     &     age,semi,ecc,vkick,zpars)
*
*
*     Hyperbolic merger
*

      INTEGER kstar(2),I1,I2
      REAL*8  mass0(2),mass(2),rad(2),massc(2),radc(2),age(2)
      REAL*8  semi,peri
      REAL*8  TSCLS(20),LUMS(10),GB(10),zpars(20),TM,TN
      REAL*8  LUM,Q,RL1,RL2
      REAL*8  MENV(2),RENV(2),K2,ospin(2),jspin(2),jorb,vkick(8)
      real*8 FBFAC,FBTOT,MCO,k3,EPS,TOL
      integer ECS
      PARAMETER (EPS = 1.0D-06, TOL = 1.0D-04, k3=0.21d0)
      REAL*8 RL
      EXTERNAL RL
      LOGICAL coel

      I1 = 1
      I2 = 2
*       Determine indices for primary & secondary star (donor & accretor).
      Q = mass(I1)/mass(I2)
      PERI = SEMI*(1-ecc)
      RL1 = RL(Q)*ABS(PERI)
*       Evaluate Roche radius for the second star.
      Q = 1.0/Q
      RL2 = RL(Q)*ABS(PERI)
*
*       Compare scaled Roche radii when choosing the primary.
      IF (RAD(I1)/RL1.LT.RAD(I2)/RL2) THEN
          I1 = 2
          I2 = 1
          RL1 = RL2
      END IF

      IF(RAD(I1).GE.RL1)THEN

*     Convert Roche radius, radius and initial & current mass to SU.
         
         DO k=1,2 
            CALL star(kstar(k),mass0(k),mass(k),tm,tn,tscls,
     &           lums,GB,zpars)
            CALL hrdiag(mass0(k),age(k),mass(k),tm,tn,tscls,lums,
     &           GB,zpars,rad(k),lum,kstar(k),massc(k),radc(k),
     &           menv(k),renv(k),k2,fbfac,fbtot,mco,ecs)
            jspin(k) = ospin(k)*(k2*rad(k)*rad(k)*(mass(k)-massc(k))
     &           +k3*radc(k)*radc(k)*massc(k))
         END DO

         CALL comenv(mass0(i1),mass(i1),massc(i1),age(i1),jspin(i1),
     &        kstar(i1),mass0(i2),mass(i2),massc(i2),age(i2),
     &        jspin(i2),kstar(i2),zpars,ecc,semi,jorb,
     &        vkick(4*(i1-1)+1),vkick(4*(i2-1)+1),coel)
         
      ELSE
         CALL MIX(mass0,mass,age,kstar,zpars)
      END IF

      RETURN

      END
