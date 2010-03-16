C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
C
      SUBROUTINE PR_KIN
C
C     This program calculates all kinematic variabls for Real Prim. Exp.
C                 ETA Primakoff
C
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
      COMMON/GCKINE/IKINE,PKINE(10),ITRA,ISTAK,IVERT,IPART,ITRTYP
     +      ,NAPART(5),AMASS,CHARGE,TLIFE,VERT(3),PVERT(4),IPAOLD
C
      COMMON/KINEM1/EINI,TPI0,FIPI0G,EPI0SP,TRECSP,TKINRM
      COMMON/KINEM2/PPREL(3)
      COMMON/KINEM3/EPI0LF,PPI0LF(3),EG1LF,PG1LF(3),EG2LF,PG2LF(3)
C
      COMMON/TMCAR1/FTPNOR(1:150,0:10000),tpi0ms(10000)
C
      COMMON/TMCAR4/egammn,egammx,denest
      COMMON/BRMST1/egamak
C
      DIMENSION RNDM(2)
C
      DOUBLE PRECISION egammx,egammn,denest
C
      DOUBLE PRECISION C1,A1,A2,A3
      DOUBLE PRECISION BMASS,PI0MAS,TPI0DP,TRECDP,TKINRC,ENERRC
      DOUBLE PRECISION PPI0DP,EPI0DP
      DOUBLE PRECISION EGBEDP,DISCRI,U,ARECOI,TCONS1
C
      DOUBLE PRECISION PI,TWOPI,PIBY2,DEGRAD,RADDEG,EMASS
C
      PARAMETER (PI=3.14159265358979324)
      PARAMETER (TWOPI=6.28318530717958648)
      PARAMETER (PIBY2=1.57079632679489662)
      PARAMETER (DEGRAD=0.0174532925199432958)
      PARAMETER (RADDEG=57.2957795130823209)
      PARAMETER (EMASS=0.0005109990615)
C
CC for pio exp.      PARAMETER (PI0MAS=0.1349764d0)
C
      PARAMETER (PI0MAS=0.54745d0)
C
      U = 0.931502D0
C
C He4      ARECOI = 4.00260D0
C
      ARECOI = 1.00726568D0                    ! for  Proton
C
      BMASS = ARECOI*U
C 
      EINI  = egamak
C
      EGBEDP = DBLE(egamak)
C
C                                         Neew Monte Carlo for teta pi0
C                         selects the prod. angle vs. cross section
      CALL ANGLSP(tpisel)
C
      TPI0DP = DBLE(tpisel)                  ! in degrees
C
      C1=(2.D0*EGBEDP*BMASS+PI0MAS*PI0MAS)/2.D0
      TCONS1=(EGBEDP+BMASS)**2
      A1=EGBEDP*EGBEDP*(DCOS(DEGRAD*TPI0DP))**2-TCONS1
      A2=2.D0*EGBEDP*C1*DCOS(DEGRAD*TPI0DP)
      A3=C1*C1-PI0MAS*PI0MAS*(EGBEDP+BMASS)**2
C
      DISCRI = (A2*A2-4.D0*A1*A3)
      IF(DISCRI.LE.0.D0)THEN
      WRITE(6,*)'Discriminator under SQRT is LT than 0.'
C
      ELSE
         PPI0DP=(-A2-DSQRT(DISCRI))/(2.D0*A1)
C
      EPI0DP=DSQRT(PPI0DP*PPI0DP+PI0MAS*PI0MAS)
      TKINRC=EGBEDP-EPI0DP
      ENERRC=BMASS+TKINRC
      PMOMRC=DSQRT(ENERRC*ENERRC-BMASS*BMASS)
      TRECDP=RADDEG*(DASIN((PPI0DP*DSIN(DEGRAD*TPI0DP))/PMOMRC))
C
         TPI0 = tpisel               ! in degrees
C
C                              sampling the fipi0 in 0-360 degrees
C
      CALL RANLUX(RNDM,2)
C
      FIPI0R = 2.*3.14159265*RNDM(2)        ! in radians
      FIPI0G = 57.2957795*FIPI0R
C
         EPI0LF=SNGL(EPI0DP)
         EPI0SP=EPI0LF
         PPI0 = SNGL(PPI0DP)
         TRECSP=SNGL(TRECDP)
         TKINRM=1000.*SNGL(TKINRC)          ! T kin of Rec. Nuc. in MeV
C
         PPI0LF(1)=PPI0*SIN(DEGRAD*TPI0)*COS(FIPI0R)
         PPI0LF(2)=PPI0*SIN(DEGRAD*TPI0)*SIN(FIPI0R)
         PPI0LF(3)=PPI0*COS(DEGRAD*TPI0)
C
         ENDIF
C
C
      CALL HF1(1,EINI,1.)
      CALL HF1(2,TPI0,1.)
      CALL HF1(3,FIPI0G,1.)
      CALL HF1(4,EPI0LF,1.)
      CALL HF1(5,TRECSP,1.)
      CALL HF1(6,TKINRM,1.)
C
C
         CALL ETA_CM
C
C
      RETURN
      END
C
C
C......................
C
C
C........................................................................
C
      SUBROUTINE ANGLSP(tpisel)
C
C     This subr. samples the theta ETA  for Real Prim. Exp.
C........................................................................
C
      COMMON/TMCAR1/FTPNOR(1:150,0:10000),tpi0ms(10000)
      COMMON/TMCAR4/egammn,egammx,denest
      COMMON/BRMST1/egamak
C
      DOUBLE PRECISION egammx,egammn,denest
C
C
C    definition of initial g beam energy bin for theta pi0 sampling, new
C
      igebin = 1 + IDINT((DBLE(egamak)-egammn)/denest)
C
      if(igebin.GT.200)write(6,*)'? S. wrong in ANGLSP,igebin=',igebin 
C
C         Sampling of Theta pi0 angle according cross section*domega  ******
C
      CALL RANLUX(URN,1)
C
      k=0
  1   k=k+1
C
      IF(k.GT.10000) go to 2
C
      IF(FTPNOR(igebin,k).LT.URN) go to 1
C
      tpisel = tpi0ms(k)            ! selected teta pi0 angle in degrees
C
      RETURN
C
 2    write(6,*)'Subr. ANGLSP**, The URN bigger than Max of FTPNOR ??'
C
      RETURN
      END
C








