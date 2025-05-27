***
      SUBROUTINE TRFLOW(kstar,mass0,mass,rad,massc,radc,age,
     &     dtr,semi,ecc,zpars)
      implicit none
*
*
*       Time until Roche overflow.
*       --------------------------
*
*
      INTEGER kstar(2),I1,I2,IT,KW
      REAL*8  mass0(2),mass(2),rad(2),massc(2),radc(2),age(2),dtr
      REAL*8  semi,aursun,ecc
      REAL*8  TSCLS(20),LUMS(10),GB(10),zpars(20),TM,TN
      REAL*8  M0,M1,LUM,MC,RC,Q,RL1,RL2
      REAL*8  MENV,RENV,K2
      REAL*8  AJ,RHIGH,THIGH,T,RR,RX,TLOW,RLOW,DEL,DER,EPS,TOL
      real*8 FBFAC,FBTOT,MCO
      integer ECS
      PARAMETER (EPS=1.0D-06, TOL=1.0D-04, aursun=2.1493370240907944d2)
      REAL*8 RL,RTMSF,RGBF,RGBDF,RAGBF,RAGBDF,RZHEF
      EXTERNAL RL,RTMSF,RGBF,RGBDF,RAGBF,RAGBDF,RZHEF
      REAL*8 RSJ,CM
* Tanikawa's prescription for BH spin
      real*8 jspin(2),aspin(2)
*
*
      I1 = 1
      I2 = 2
*       Determine indices for primary & secondary star (donor & accretor).
      Q = mass(I1)/mass(I2)
      RL1 = RL(Q)*SEMI
*       Evaluate Roche radius for the second star.
      Q = 1.0/Q
      RL2 = RL(Q)*SEMI
*
*       Compare scaled Roche radii when choosing the primary.
      IF (RAD(I1)/RL1.LT.RAD(I2)/RL2) THEN
          I1 = 2
          I2 = 1
          RL1 = RL2
      END IF
*     IF(TEV(J1).LE.0.0)THEN
*        WRITE(6,99)J1,NAME(J1),NAME(I),KSTAR(J1),TEV(J1)+TOFF,TTOT
*99      FORMAT(' TRFLOW WARNING! ',3I6,I4,2F10.3)
*     ENDIF
*
*       Exit with large interval if semi-major axis is negative.
      IF(SEMI.LE.0.D0)THEN
          DTR = 1.0D+10
          GOTO 210
      ENDIF
*
*       Convert Roche radius, radius and initial & current mass to SU.
      M0 = MASS0(I1)
      M1 = MASS(I1)
      RSJ = RAD(I1)
      MC = MASSC(I1)
      RC = RADC(I1)
      KW = KSTAR(I1)
      CALL star(KW,M0,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
      CALL hrdiag(M0,AGE(I1),M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &            RSJ,LUM,KW,MC,RC,MENV,RENV,K2,fbfac,fbtot,mco,ecs,
     &     jspin(i1),aspin(i1))
*
* If the star already fills its Roche lobe or will fill after peri-center, return one period
*
      IF(RSJ.GE.RL1.OR.SEMI*(1-ECC).LT.(RAD(I1)+RAD(I2)))THEN
         DTR = semi/aursun*sqrt(semi/(aursun*(MASS(I1)+MASS(I2))))
         DTR = DTR/1.0d6
*         DTR = 0.D0
*         RAD(I1) = RS1
         GOTO 210
      ENDIF
*
*       Exit with large interval if primary has terminated its evolution.
      IF(KSTAR(I1).GE.10)THEN
          DTR = 1.0D+10
          GOTO 210
      ENDIF
*
* If the star is on the MS then see if it will fill its Roche Lobe
* during this stage.
*
      T = 2.0D+10
      IF(KW.LE.1.OR.KW.EQ.7)THEN
         AJ = 0.99D0*TM
         CALL hrdiag(M0,AJ,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &        RHIGH,LUM,KW,MC,RC,MENV,RENV,K2,fbfac,fbtot,mco,ecs,
     &        jspin(i1),aspin(i1))
         IF(RHIGH.GT.RL1.AND.AJ.GT.AGE(I1))THEN
* Find the overflow time using binary chopping.
            TLOW = AGE(I1)
            RLOW = RSJ
            THIGH = AJ
            IT = 0
   10       T = 0.5D0*(THIGH + TLOW)
            CALL hrdiag(M0,T,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &           RR,LUM,KW,MC,RC,MENV,RENV,K2,fbfac,fbtot,mco,ecs,
     &           jspin(i1),aspin(i1))
            IF(RR.LT.RL1)THEN
               TLOW = T
               RLOW = RR
            ELSE
* Stop when 0.1% accuracy is achieved.
               IT = IT + 1
*               IF(IT.EQ.25)THEN
*                  IF(ABS(RR-RL1).GT.0.1)THEN
*                     WRITE(38,*)' TRFLOW KW1: ',rr,rls,t,tm
*                  ENDIF
*               ENDIF
               IF(ABS(RR - RL1).LE.0.001*RR.OR.IT.GT.25) GOTO 200
               THIGH = T
               RHIGH = RR
            ENDIF
            GOTO 10
         ENDIF
      ENDIF
*
* See if the star will fill its Roche Lobe on the HG if not already
* past this stage.
*
      IF(KW.LE.2)THEN
         AJ = MIN(TN,TSCLS(1))*(1.D0-EPS)
         CALL hrdiag(M0,AJ,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &        RHIGH,LUM,KW,MC,RC,MENV,RENV,K2,fbfac,fbtot,mco,ecs,
     &        jspin(i1),aspin(i1))
         IF(RHIGH.GT.RL1.AND.AJ.GT.AGE(I1))THEN
* Solve for the overflow time.
            RR = RTMSF(M0)
            T = LOG(RL1/RR)/LOG(RHIGH/RR)
            T = TM + T*(TSCLS(1) - TM)
            GOTO 200
         ENDIF
         IF(TN.LE.TSCLS(1)) GOTO 200
      ENDIF
*
* If the star is pre-helium ignition see if it will fill its Roche
* lobe before helium ignition.
*
      IF(KW.LE.3)THEN
         AJ = MIN(TN,TSCLS(2))*(1.D0-EPS)
         CALL hrdiag(M0,AJ,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &        RHIGH,LUM,KW,MC,RC,MENV,RENV,K2,fbfac,fbtot,mco,ecs,
     &        jspin(i1),aspin(i1))
         IF(RHIGH.GE.RL1)THEN
* Solve for the luminosity corresponding to the Roche radius and
* then for the overflow time.
            IT = 0
            LUM = LUMS(3)
 30         IT = IT + 1
            DEL = RGBF(M1,LUM) - RL1
*            IF(IT.EQ.25.AND.ABS(DEL).GT.0.1)THEN
*               if(rank.eq.0)then
*               WRITE(38,*)' TRFLOW KW3: ',rgbf(m1,lum),rls,lum,lums(3)
*               end if
*            ENDIF
            IF(ABS(DEL/RL1).LE.TOL.OR.IT.EQ.25) GOTO 40
            DER = RGBDF(M1,LUM)
            LUM = LUM - DEL/DER
            GOTO 30
 40         CONTINUE
            IF(LUM.LE.LUMS(6))THEN
               T = TSCLS(4) - (1.D0/((GB(5)-1.D0)*GB(1)*GB(4)))*
     &                      ((GB(4)/LUM)**((GB(5)-1.D0)/GB(5)))
            ELSE
               T = TSCLS(5) - (1.D0/((GB(6)-1.D0)*GB(1)*GB(3)))*
     &                      ((GB(3)/LUM)**((GB(6)-1.D0)/GB(6)))
            ENDIF
            GOTO 200
         ELSE
* If a low mass star has not yet filled its Roche Lobe then it will
* not do it until after the Helium Flash. Thus we don't let it go
* any further until it has actually become a type 4.
            IF(M0.LT.ZPARS(2).OR.TN.LE.TSCLS(2)) GOTO 200
         ENDIF
      ENDIF
*
* Check for overflow during the CHeB stage.
*
      IF(KW.EQ.4)THEN
*
         IF(TN.LT.(TSCLS(2)+TSCLS(3)))THEN
            T = MAX(TSCLS(2),AGE(I1))
            AJ = T + 0.5D0*(TN - T)
            CALL hrdiag(M0,AJ,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &           RHIGH,LUM,KW,MC,RC,MENV,RENV,K2,fbfac,fbtot,mco,ecs,
     &           jspin(i1),aspin(i1))
            IF(RHIGH.LT.RL1)THEN
* If the evolution is due to end during CHeB then the pertubation
* functions for small envelope mass will take effect at some point
* causing the radius to decrease. We assume that this point is after
* AJ and quit with T as an underestimate of the overflow time.
               T = AJ
               GOTO 200
            ENDIF
         ELSE
            AJ = (TSCLS(2)+TSCLS(3))*(1.D0-EPS)
            CALL hrdiag(M0,AJ,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &           RHIGH,LUM,KW,MC,RC,MENV,RENV,K2,fbfac,fbtot,mco,ecs,
     &           jspin(i1),aspin(i1))
            IF(RHIGH.LT.RL1) GOTO 50
         ENDIF
* Find the overflow time using binary chopping.
         TLOW = AGE(I1)
         RLOW = RSJ
         THIGH = AJ
         IT = 0
 60      T = 0.5D0*(THIGH + TLOW)
         CALL hrdiag(M0,T,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &        RR,LUM,KW,MC,RC,MENV,RENV,K2,fbfac,fbtot,mco,ecs,
     &        jspin(i1),aspin(i1))
         IF(RR.LT.RL1)THEN
            TLOW = T
            RLOW = RR
         ELSE
* Stop when 0.1% accuracy is achieved.
            IT = IT + 1
*            IF(IT.EQ.25)THEN
*               IF(rank.eq.0.and.ABS(RR-RL1).GT.0.1)THEN
*                  WRITE(38,*)' TRFLOW KW4: ',rr,rls,tscls(2),t,tscls(3)
*               ENDIF
*            ENDIF
            IF(ABS(RR - RL1).LE.0.001*RR.OR.IT.GT.25) GOTO 200
            THIGH = T
            RHIGH = RR
         ENDIF
         GOTO 60
 50      CONTINUE
      ENDIF
*
* The star now has only until the end of the AGB phase to fill its lobe.
*
      IF(KW.LE.6)THEN
* Solve for the luminosity corresponding to the Roche radius and
* then for the overflow time.
         IT = 0
         IF(KW.EQ.6)THEN
            LUM = LUMS(8)
         ELSE
            LUM = LUMS(7)
         ENDIF
         RR = RAGBF(M1,LUM,ZPARS(2))
         IF(RR.GT.RL1)THEN
* In this case the solution should already have been found. If it hasn't
* then most likely the envelope is small and the pertubation functions
* are taking effect.
            T = TN
            GOTO 200
         ENDIF
 70      IT = IT + 1
         RR = RAGBF(M1,LUM,ZPARS(2))
         DEL = RR - RL1
*         IF(IT.EQ.25.AND.ABS(DEL).GT.0.1)THEN
*            if(rank.eq.0) WRITE(38,*)' TRFLOW KW6: ',rr,rls,lum,lums(7)
*         ENDIF
         IF(ABS(DEL/RL1).LE.TOL.OR.IT.EQ.25) GOTO 80
         DER = RAGBDF(M1,LUM,ZPARS(2))
         LUM = LUM - DEL/DER
         GOTO 70
 80      CONTINUE
         IF(LUM.LE.LUMS(8))THEN
            IF(LUM.LE.LUMS(6))THEN
               T = TSCLS(7) - (1.D0/((GB(5)-1.D0)*GB(8)*GB(4)))*
     &                      ((GB(4)/LUM)**((GB(5)-1.D0)/GB(5)))
            ELSE
               T = TSCLS(8) - (1.D0/((GB(6)-1.D0)*GB(8)*GB(3)))*
     &                      ((GB(3)/LUM)**((GB(6)-1.D0)/GB(6)))
            ENDIF
         ELSE
            IF(LUM.LE.LUMS(6))THEN
               T = TSCLS(10) - (1.D0/((GB(5)-1.D0)*GB(2)*GB(4)))*
     &                       ((GB(4)/LUM)**((GB(5)-1.D0)/GB(5)))
            ELSE
               T = TSCLS(11) - (1.D0/((GB(6)-1.D0)*GB(2)*GB(3)))*
     &                       ((GB(3)/LUM)**((GB(6)-1.D0)/GB(6)))
            ENDIF
         ENDIF
         GOTO 200
      ENDIF
*
      IF(KW.LE.9)THEN
*
* The star has until the end of the Helium GB to fill its lobe.
*
         LUM = (RL1/0.08D0)**(4.0/3.0)
         CM = 2.0D-03*M1**2.5/(2.D0 + M1**5)
         IF(CM*LUM.GE.100.D0)THEN
            DTR = 1.0D+10
            GOTO 210
         ENDIF
         RX = RZHEF(M1)
         RR = RX*(LUM/LUMS(2))**0.2 + 
     &            0.02D0*(EXP(CM*LUM) - EXP(CM*LUMS(2)))
         IF(RR.LT.RL1)THEN
* Find the overflow luminosity using binary chopping.
            TLOW = LUM
            RLOW = RR
 89         LUM = 2.0*LUM
            RR = RX*(LUM/LUMS(2))**0.2 + 
     &               0.02D0*(EXP(CM*LUM) - EXP(CM*LUMS(2)))
            IF(RR.LE.RL1) GOTO 89
            THIGH = LUM
            RHIGH = RR
            IT = 0
 90         LUM = 0.5D0*(THIGH + TLOW)
            RR = RX*(LUM/LUMS(2))**0.2 + 
     &               0.02D0*(EXP(CM*LUM) - EXP(CM*LUMS(2)))
            IF(RR.LT.RL1)THEN
               TLOW = LUM
               RLOW = RR
            ELSE
* Stop when 0.1% accuracy is achieved.
               IT = IT + 1
               IF(IT.EQ.25)THEN
*                  IF(ABS(RR-RL1).GT.0.1)THEN
*                     if(rank.eq.0)WRITE(38,*)' TRFLOW KW9: ',rr,rls,t,tm
*                  ENDIF
               ENDIF
               IF(ABS(RR - RL1).LE.0.001*RR.OR.IT.GT.25) GOTO 95
               THIGH = LUM
               RHIGH = RR
            ENDIF
            GOTO 90
 95         CONTINUE
         ENDIF
         IF(LUM.LE.LUMS(6))THEN
            T = TSCLS(4) - (1.D0/((GB(5)-1.D0)*GB(8)*GB(4)))*
     &                   ((GB(4)/LUM)**((GB(5)-1.D0)/GB(5)))
         ELSE
            T = TSCLS(5) - (1.D0/((GB(6)-1.D0)*GB(8)*GB(3)))*
     &                   ((GB(3)/LUM)**((GB(6)-1.D0)/GB(6)))
         ENDIF
      ENDIF
*
 200  CONTINUE
      IF(T.GE.TN)THEN
         DTR = 1.0D+10
      ELSE
         DTR = T - AGE(I1)
*         TT = (T+EPOCH(I1))/TSTAR
*         DTR = TT - TIME
*         IF(DTR.LT.0.D0.AND.ITEV)THEN
*            TEV(J1) = MIN(TEV(J1),TIME)
*            TEV(J2) = TEV(J1)
*         ENDIF
         DTR = MAX(DTR,0.D0)
*         IF(TEV(J1).LT.TEV0(J1))THEN
*            if(rank.eq.0) WRITE (6,101)J1,KSTAR(J1),TEV(J1)+TOFF,TTOT
* 101        FORMAT(' TRFLOW WARNING! TEV<TEV0 ',I6,I4,2F10.3)
*         ENDIF
      ENDIF
*
 210  RETURN
      END
***
