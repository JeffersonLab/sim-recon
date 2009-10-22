C          sum of all amplituds on proton, 121208
C            in eta lifetime experiment
C
CC      FUNCTION etaall(thpi0)
C
      SUBROUTINE CRSEC6(egama,thpi0,sumall)
C
C      DOUBLE PRECISION EMMU,PMASS,AVO
C
      DOUBLE PRECISION thpi0
C
      DOUBLE PRECISION C1,A1,A2,A3
      DOUBLE PRECISION AMASS,PI0MAS,TPI0DP,TRECDP,TKINRC,ENERRC
      DOUBLE PRECISION t,Delta2
      DOUBLE PRECISION PPI0DP,EPI0DP,PI0BET
      DOUBLE PRECISION EGBEDP,DISCRI,U,ARECOI,TCONS1
C
      DOUBLE PRECISION DCONS3,DCONS1,DCONS2,DCONS4,DCONS5
      DOUBLE PRECISION DCON11,DCON12,DCON14,DCON25,DCON26
      DOUBLE PRECISION BCOULM,BNCOHE
      DOUBLE PRECISION TCOULM,TNCOHE
C
      DOUBLE PRECISION qGeV,qMeV,qfm,ALPHA,Z,Gpi02g
      DOUBLE PRECISION q,F
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
      PARAMETER (PI0MAS=0.54745)
C
      U = 0.931502D0
C
CPb      Z=82.00D0                            ! Pb nucleus
CPb      ARECOI = 207.976627D0                ! for Pb^208
C
CC12      Z= 6                                    ! for C(6-12)
CC12      ARECOI = 12.0000                        ! 
C
C
CHe4      Z= 2
CHe4      ARECOI = 4.002603D0                  ! for  He^4
C
      Z= 1.0D0
      ARECOI = 1.00726568D0                    ! for  Proton
C
C
      AMASS = ARECOI*U
C
      ALPHA = 1.D0/137.0359895D0           !Fine structure constant.
C
CC      EGBEDP = 11.000D0                        ! Initial g Energy in GeV
C
      EGBEDP = DBLE(egama)                ! Tagged g Energy in GeV
C
      TPI0DP = DBLE(thpi0)                  ! in degrees
C
      C1=(2.D0*EGBEDP*AMASS+PI0MAS*PI0MAS)/2.D0
      TCONS1=(EGBEDP+AMASS)**2
      A1=EGBEDP*EGBEDP*(DCOS(DEGRAD*TPI0DP))**2-TCONS1
      A2=2.D0*EGBEDP*C1*DCOS(DEGRAD*TPI0DP)
      A3=C1*C1-PI0MAS*PI0MAS*(EGBEDP+AMASS)**2
C
      DISCRI = (A2*A2-4.D0*A1*A3)
      IF(DISCRI.LE.0.D0)THEN
C
      WRITE(6,*)'Discriminator under SQRT is LT than 0.'
C
      ELSE
         PPI0DP=(-A2-DSQRT(DISCRI))/(2.D0*A1)
C
         EPI0DP=DSQRT(PPI0DP*PPI0DP+PI0MAS*PI0MAS)
         TKINRC=EGBEDP-EPI0DP
         ENERRC=AMASS+TKINRC
         PMOMRC=DSQRT(ENERRC*ENERRC-AMASS*AMASS)
         TRECDP=RADDEG*(DASIN((PPI0DP*DSIN(DEGRAD*TPI0DP))/PMOMRC))
C
         DCS21 = DCOS(DEGRAD*TPI0DP)
         t = PI0MAS*PI0MAS-2.D0*(EGBEDP*EPI0DP-EGBEDP*PPI0DP*DCS21)
C
         TPI0 = thpi0                ! in degrees
C
C
CC      FIPI0R = 2.*3.14159265*RNDM(2)        ! in radians
CC      FIPI0G = 57.2957795*FIPI0R
C
         EPI0LF=SNGL(EPI0DP)
         EPI0SP=EPI0LF
         PPI0 = SNGL(PPI0DP)
         TRECSP=SNGL(TRECDP)
         TKINRM=1000.*SNGL(TKINRC)          ! T kin of Rec. Nuc. in MeV
C
         Delta2 = -2.D0*AMASS*TKINRC        ! the same as t in Mandelstamm
         ENDIF
C
C
       qGeV = DSQRT(2.D0*AMASS*TKINRC)
       qMeV =1000.D0*DSQRT(2.D0*AMASS*TKINRC)
       qfm  =qMeV/197.327053D0                            ! in 1/fm
C
       q = qfm
C
C for He4      CALL FFACT(q,F)            ! new for test
C                         for proton
      CALL FFACT(qGeV,F)
C
C
CCC       Gpi02g= 1.0D-9                  ! is set to 1 ev
CCC       BCOULM= 7.92D0
C
       Gpi02g= 1.0D-6                  ! is set Kev to 1 Gev
CC       BCOULM= 0.343
       BCOULM= 0.500
C
C
C                            for eta
CCC       BNCOHE= 0.25
C
CC       BNCOHE= 1.02        ! changed on 05/19/00
C       BNCOHE= 0.5
C                      from Laget, should be ~3.0
       BNCOHE= 2.0
C
CE=4.4       BFIINT= 1.18       ! rad's which is 67.61 degrees
C
CCC       BFIINT= 1.01          ! rad's which is 57.61 degrees, E=5.6 GeV
C
       BFIINT= -1.00             ! for eta
C
C
       PI0BET= PPI0DP/EPI0DP
C
       DCONS3= 0.389379D0               ! GeV-2 to mb
C
C  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C                                       Squere of Primakoff amplitude
C
       DCONS1=8.D0*ALPHA*Z*Z*Gpi02g*PI0BET**3*EGBEDP**4
       DCONS2=((DSIN(DEGRAD*TPI0DP))**2*F*F)/(PI0MAS**3*Delta2**2)
C
       dsdopr=SNGL(DCONS3*BCOULM*DCONS1*DCONS2)       ! in millibarn/sr
C
C >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C                                       ! Nuclear Coherent
C
CC       DCONS4=100.D0*EGBEDP**2                 ! in microbarn (K in GeV)
       DCONS4=16.D0*EGBEDP**2 
       DCONS5=(DSIN(DEGRAD*TPI0DP))**2*F*F
C
       dsdonc=SNGL(BNCOHE*ARECOI**2*DCONS4*DCONS5)/1000. !in milibarn/sr
C
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C                                ! Interference between pr and nc
C                                       Squere of Primakoff amplitude
C
       DCON11=DSQRT(8.D0*ALPHA*Z*Z*Gpi02g)
       DCON12=EGBEDP**2*(DSIN(DEGRAD*TPI0DP))*F/ABS(Delta2)
       DCON14=DSQRT(PI0BET**3/PI0MAS**3)
C
       TCOULM=DSQRT(DCONS3)*DCON11*DCON12*DCON14   ! in sqrt(millibarn)
C
C
C                                     ! squere of Nucl Coherent
C
       DCON25=DSQRT(16.D0*EGBEDP**2/1000.D0) !LinSQRT(microbarn KinGeV)
       DCON26=(DSIN(DEGRAD*TPI0DP))*F
C
       TNCOHE=ARECOI*DCON25*DCON26              !in sqrt(milibarn)
C
C
       con101= SNGL(DSQRT(BCOULM))
       con102= SNGL(TCOULM)
       con103= SNGL(DSQRT(BNCOHE))
       con104= SNGL(TNCOHE)
C
       sintrf=2.*con101*con102*con103*con104*COS(BFIINT)
C
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C                                    Nuclear Incoherent
C
CE=4.4       BNINCO = 4.10
       BNINCO = 6.3                     ! E=5.6 GeV
C                      for eta
       BNINCO = 2.26
C
CC       dsdoic=BNINCO*SNGL(DSQRT(DSQRT(ARECOI**3)))/1000.  !milib/sr
C
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
C                                    Sum of all parts
C
CCCCC       etaall=dsdopr+sintrf+dsdonc+dsdoic
CCC       etaall=dsdopr+sintrf+dsdonc
C
       sumall=dsdopr+sintrf+dsdonc
C
      return
      end
C
C
C..........
C
C
C
      SUBROUTINE FFACT(q,F)
C
      DOUBLE PRECISION q,F
C
C   Form Factor for Proton
C
C He4      F = (1.D0-(a*a*q*q)**6)*dexp(-b*b*q*q)
C
      F = 1.D0/(1.D0+q*q/0.7D0)**2
C
C
      return
      end
C








