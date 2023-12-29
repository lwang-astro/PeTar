***
      SUBROUTINE MIX(M0,M,AJ,KS,ZPARS)
*
*     Author : J. R. Hurley
*     Date :   7th July 1998
*
*       Evolution parameters for mixed star.
*       ------------------------------------
*
      implicit none
*
      INTEGER KS(2),I1,I2,K1,K2,KW,ICASE
      INTEGER KTYPE(0:14,0:14)
      COMMON /TYPES/ KTYPE
      REAL*8 M0(2),M(2),AJ(2),ZPARS(20)
      REAL*8 TSCLS(20),LUMS(10),GB(10),TMS1,TMS2,TMS3,TN
      REAL*8 M01,M02,M03,M1,M2,M3,AGE1,AGE2,AGE3,MC3,MCH
      PARAMETER(MCH=1.44D0)

****
      REAL*8 ftzacc
      PARAMETER (ftzacc=0.5D0)
****
      REAL*8 sigma, mxns
      INTEGER bhflag
      COMMON /VALUE4/ sigma, mxns, bhflag
* Tanikawa's BH model (prevent HGCE)
      logical mixgiant
      logical nohgce
      integer krol(2)
      common /bseemp/ nohgce,krol
*
*
* Tanikawa's BH model (prevent HGCE)
      icase = ktype(krol(1),krol(2))
*
*       Define global indices with body #I1 being most evolved.
      IF(KS(1).GE.KS(2))THEN
          I1 = 1
          I2 = 2
      ELSE
          I1 = 2
          I2 = 1
      END IF
*
*       Specify case index for collision treatment.
      K1 = KS(I1)
      K2 = KS(I2)
* Tanikawa's BH model (prevent HGCE)
!      ICASE = KTYPE(K1,K2)
*
*     if(icase.gt.100) WRITE(66,*)' MIX ERROR ICASE>100 ',icase,k1,k2
*
*       Determine evolution time scales for first star.
      M01 = M0(I1)
      M1 = M(I1)
      AGE1 = AJ(I1)
      CALL star(K1,M01,M1,TMS1,TN,TSCLS,LUMS,GB,ZPARS)
*
*       Obtain time scales for second star.
      M02 = M0(I2)
      M2 = M(I2)
      AGE2 = AJ(I2)
      CALL star(K2,M02,M2,TMS2,TN,TSCLS,LUMS,GB,ZPARS)
*
*       Check for planetary systems - defined as HeWDs and low-mass WDs!
      IF(K1.EQ.10.AND.M1.LT.0.05)THEN
         ICASE = K2
         IF(K2.LE.1)THEN
            ICASE = 1
            AGE1 = 0.D0
         ENDIF
      ELSEIF(K1.GE.11.AND.M1.LT.0.5.AND.ICASE.EQ.6)THEN
         ICASE = 9
      ENDIF
      IF(K2.EQ.10.AND.M2.LT.0.05)THEN
         ICASE = K1
         IF(K1.LE.1)THEN
            ICASE = 1
            AGE2 = 0.D0
         ENDIF
      ENDIF
*
*       Specify total mass.
      M3 = M1 + M2
      M03 = M01 + M02
      KW = ICASE
      AGE3 = 0.d0
*
*       Restrict merged stars to masses less than 100 Msun. 
*      IF(M3.GE.100.D0)THEN
*         M3 = 99.D0
*         M03 = MIN(M03,M3)
*      ENDIF
*
*       Evaluate apparent age and other parameters.
*
      IF(ICASE.EQ.1)THEN
*       Specify new age based on complete mixing.
         IF(K1.EQ.7) KW = 7
         CALL star(KW,M03,M3,TMS3,TN,TSCLS,LUMS,GB,ZPARS)
         AGE3 = 0.1d0*TMS3*(AGE1*M1/TMS1 + AGE2*M2/TMS2)/M3
* Tanikawa's BH model (prevent HGCE)
      else if(nohgce .and. icase.eq.2)then
         call star(kw,m03,m3,tms3,tn,tscls,lums,gb,zpars)
         age3 = 0.1d0*tms3*(age1*m1/tms1+age2*m2/tms2)/m3
*
      ELSEIF(ICASE.EQ.3.OR.ICASE.EQ.6.OR.ICASE.EQ.9)THEN
         MC3 = M1
         CALL gntage(MC3,M3,KW,ZPARS,M03,AGE3)
* Tanikawa's BH model (prevent HGCE)
!      ELSEIF(ICASE.EQ.4)THEN
      else if(mixgiant(nohgce,icase))then
*
         MC3 = M1
         AGE3 = AGE1/TMS1
         CALL gntage(MC3,M3,KW,ZPARS,M03,AGE3)
      ELSEIF(ICASE.EQ.7)THEN
         CALL star(KW,M03,M3,TMS3,TN,TSCLS,LUMS,GB,ZPARS)
         AGE3 = TMS3*(AGE2*M2/TMS2)/M3
      ELSEIF(ICASE.LE.12)THEN
*       Ensure that a new WD has the initial mass set correctly.
         M03 = M3
         IF(ICASE.LT.12.AND.M3.GE.MCH)THEN
            M3 = 0.D0
            KW = 15
         ENDIF
      ELSEIF(ICASE.EQ.13.OR.ICASE.EQ.14)THEN
*       Set unstable Thorne-Zytkow object with fast mass loss of envelope 
*       unless the less evolved star is a WD, NS or BH. 
         IF(K2.LT.10)THEN
            M03 = M1
            M3 = M1
**** Condition to add fraction of mass to BH in case of
**** main-seuqence BH collision         
            IF((K1.EQ.14).AND.(K2.LE.1)) M3 = M1 + (ftzacc*M2)
         ENDIF
         IF(ICASE.EQ.13.AND.M3.GT.MXNS) KW = 14
      ELSEIF(ICASE.EQ.15)THEN
         M3 = 0.D0
      ELSEIF(ICASE.GT.100)THEN
*       Common envelope case which should only be used after COMENV.
         KW = K1
         AGE3 = AGE1
         M3 = M1
         M03 = M01
      ELSE
*       This should not be reached.
        KW = 1
        M03 = M3
      ENDIF
*
* Put the result in *1.
*
      KS(1) = KW
      KS(2) = 15
      M(1) = M3
      M(2) = 0.D0
      M0(1) = M03
      AJ(1) = AGE3
*
      RETURN
      END
***
      logical function mixgiant(nohgce,icase)
      implicit none
      logical nohgce
      integer icase
      
      mixgiant = .false.
      if(nohgce)then
         if(icase.eq.4 .or. icase.eq.5)then
            mixgiant = .true.
         endif
      else
         if(icase.eq.4)then
            mixgiant = .true.
         endif
      endif

      return
      end
