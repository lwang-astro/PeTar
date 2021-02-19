***
      SUBROUTINE trdot(kw,m0,m1,rm,mc,rc,age,dt,dtr,zpars)
      implicit none
*
      INTEGER kw,kw0,it
      REAL*8 age,tm,tn,tscls(20),LUMS(10),GB(10),zpars(20)
      REAL*8 dt,dtr
      REAL*8 M0,M1,RM,LUM,MC,MC1,RC,RM0,AGE0,M10
      REAL*8 menv,renv,k2
      REAL*8 dr,dtdr
      REAL*8 pts1,pts2,pts3,eps,alpha2,tol
      real*8 FBFAC,FBTOT,MCO
      integer ECS
      COMMON /POINTS/ pts1,pts2,pts3
      PARAMETER(eps=1.0d-06,alpha2=0.09d0,tol=1.0d-10)
*
      kw0 = kw
      IF(M1.LE.0.0) M1 = rm
      M10 = M1
      CALL star(KW,M0,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
      CALL hrdiag(M0,AGE,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &            RM,LUM,KW,MC,RC,MENV,RENV,K2,fbfac,fbtot,mco,ecs)

*       Quit if there is a change of type at the current TEV.
      if((kw0.le.6.and.kw.gt.6).or.
     &     (kw0.le.9.and.kw.gt.9))then
         m1 = m10
         dt = 0.d0
         dtr= 0.d0
         return
      endif
***

*     Base new time scale for changes in radius & mass on stellar type.
*
      if(kw.le.1)then
         dt = pts1*tm
         dtr = tm - age
      elseif(kw.eq.2)then
         dt = pts1*(tscls(1) - tm)
         dtr = tscls(1) - age
      elseif(kw.eq.3)then
         if(age.lt.tscls(6))then
            dt = pts2*(tscls(4) - age)
         else
            dt = pts2*(tscls(5) - age)
         endif
         dtr = MIN(tscls(2),tn) - age
      elseif(kw.eq.4)then
         dt = pts2*tscls(3)
         dtr = MIN(tn,tscls(2) + tscls(3)) - age
      elseif(kw.eq.5)then
         if(age.lt.tscls(9))then
            dt = pts3*(tscls(7) - age)
         else
            dt = pts3*(tscls(8) - age)
         endif
         dtr = MIN(tn,tscls(13)) - age
      elseif(kw.eq.6)then
         if(age.lt.tscls(12))then
            dt = pts3*(tscls(10) - age)
         else
            dt = pts3*(tscls(11) - age)
         endif
         dt = MIN(dt,0.005d0)
         dtr = tn - age
      elseif(kw.eq.7)then
         dt = pts1*tm
         dtr = tm - age
      elseif(kw.eq.8.or.kw.eq.9)then
         if(age.lt.tscls(6))then
            dt = pts2*(tscls(4) - age)
         else
            dt = pts2*(tscls(5) - age)
         endif
         dtr = tn - age
      else
*        dt = MAX(0.1d0,age*10.d0)
         dt = MAX(0.1d0,dt*10.d0)
         dt = MIN(dt,5.0d+02)
         dtr = dt
      endif
*

* Record radius.
*
      rm0 = rm
      if(kw.ge.10) goto 30
      age0 = age
      kw0 = kw
      mc1 = mc
*
* Check for type change.
*
      it = 0
      if((dtr-dt).le.tol)then
*
* Check final radius for too large a jump.
*
         age = MAX(age,age*(1.d0-eps) + dtr)
         CALL hrdiag(M0,AGE,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &               RM,LUM,KW,MC1,RC,MENV,RENV,K2,fbfac,fbtot,mco,ecs)
         dr = rm - rm0
         if(ABS(dr).gt.0.1*rm0)then
            dt = dtr - age0*eps
            dtdr = dt/ABS(dr)
            dt = alpha2*MAX(rm,rm0)*dtdr
            goto 20
         else
            dt = dtr
            goto 30
         endif
      endif
*
* Limit to a 10% increase assuming no further mass loss
* and thus that the pertubation functions due to small envelope mass
* will not change the radius.
*
   20 age = age0 + dt
      mc1 = mc
      CALL hrdiag(M0,AGE,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &            RM,LUM,KW,MC1,RC,MENV,RENV,K2,fbfac,fbtot,mco,ecs)
      dr = rm - rm0
      it = it + 1
      if(it.eq.20.and.kw.eq.4) goto 30
* Bug fix, do not iterate if the star is compact object, thanks to Ataru Takikawa
      if(kw.ge.13) goto 30
      IF(IT.GT.30)THEN
         WRITE (6,22) IT, KW, M0, DR, RM0
   22    FORMAT (' BSE DANGER!    deltat: Iter KW M0 DR RM0 ',
     &        2I4,1P,3E10.2)
         goto 30
      ENDIF
      if(ABS(dr).gt.0.1*rm0)then
         dtdr = dt/ABS(dr)
         dt = alpha2*MAX(rm0,rm)*dtdr
         if(it.ge.20) dt = 0.5d0*dt
         goto 20
      endif
*
 30   continue
*
      RETURN
      END
***
