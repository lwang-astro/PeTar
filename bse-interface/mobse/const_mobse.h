*
* const_mobse.h
*
      INTEGER idum
      COMMON /VALUE3/ idum
      INTEGER idum2,iy,ir(32)
      COMMON /RAND3/ idum2,iy
      INTEGER ktype(0:14,0:14)
      COMMON /TYPES/ ktype
      INTEGER ceflag,tflag,ifflag,nsflag,wdflag,piflag
      COMMON /FLAGS/ ceflag,tflag,ifflag,nsflag,wdflag,piflag
*
      INTEGER bhflag
      REAL*8 sigma1,sigma2,mxns
      REAL*8 neta,bwind,hewind,alpha1,lambda
      REAL*8 beta,xi,acc2,epsnov,eddfac,gamma
      COMMON /VALUE1/ neta,bwind,hewind
      COMMON /VALUE2/ alpha1,lambda
      COMMON /VALUE4/ sigma1,sigma2,mxns,bhflag
      COMMON /VALUE5/ beta,xi,acc2,epsnov,eddfac,gamma
      REAL*8 pts1,pts2,pts3
      COMMON /POINTS/ pts1,pts2,pts3
*      REAL*8 dmmax,drmax
*      COMMON /TSTEPC/ dmmax,drmax
*      REAL scm(50000,14),spp(20,3)
*      COMMON /SINGLE/ scm,spp
*      REAL bcm(50000,34),bpp(200,33)
*      COMMON /BINARY/ bcm,bpp
*