C             this subroutine samples the eta->g+g in center of mass
C
      SUBROUTINE ETA_CM
C
C
C
      COMMON/GCKINE/IKINE,PKINE(10),ITRA,ISTAK,IVERT,IPART,ITRTYP
     +      ,NAPART(5),AMASS,CHARGE,TLIFE,VERT(3),PVERT(4),IPAOLD
C
      COMMON/KINEM1/EINI,TPI0,FIPI0G,EPI0SP,TRECSP,TKINRM
      COMMON/KINEM2/PPREL(3)
      COMMON/KINEM3/EPI0LF,PPI0LF(3),EG1LF,PG1LF(3),EG2LF,PG2LF(3)
C
      COMMON/KING12/TG1CMG,FG1CMG
      COMMON/KING23/PSIG12
C
Cn      COMMON/G1SAMP/CMG1MM(2)
C
      DIMENSION PG1CM(3),PG2CM(3)
      DIMENSION RNDM(2)
C
      PARAMETER (PI=3.1415926)
      PARAMETER (TWOPI=6.283185)
      PARAMETER (PIBY2=1.57079632)
      PARAMETER (DEGRAD=0.01745329)
      PARAMETER (RADDEG=57.29577951)
      PARAMETER (EMASS=0.00051099)
C
      PARAMETER (EMPI0=0.54745)            ! For eta experiment
C
C
C
C        Eta -> g + g  in Eta Cener Macc
C
      EG1CM = EMPI0/2.
      EG2CM = EMPI0/2.
C
      PG1CF = EG1CM
      PG2CF = EG2CM
C
C      Random  sampling of eta -> g + g   decay in CM
C
C
C   General formula for translation (0,1) to (a,b) is
C       c=RNDM*(b-a)+a    where RNDM from (0,1) and c will be in (a,b)
C
C      TG1CM in (0,180 deg) and FG1CM in (0,180 deg),<= they are identical
C      and because d(solid angle)= d(cos(Teta))*d(Fi)
C      and   COS(0)= +1 and COS(180)= -1.
C
      TG1MIN= 0.00000
      TG1MAX= 180.00000
C
      CSMIN = COS(DEGRAD*TG1MIN)
      CSMAX = COS(DEGRAD*TG1MAX)
C
      CALL RANLUX(RNDM,2)
C
      RANDCS= RNDM(1)*(CSMAX-CSMIN)+CSMIN
      TG1CM = ACOS(RANDCS)
C
      TG1CMG= RADDEG*TG1CM  
C
      FG1CM = PI*RNDM(2)
C
      FG1CMG = RADDEG*FG1CM
C                              for first gama
C
      PG1CM(1) = PG1CF*SIN(TG1CM)*COS(FG1CM)
      PG1CM(2) = PG1CF*SIN(TG1CM)*SIN(FG1CM)
      PG1CM(3) = PG1CF*COS(TG1CM)
C
C                              for second gamma 
C
      PG2CM(1) = - PG1CM(1)
      PG2CM(2) = - PG1CM(2)
      PG2CM(3) = - PG1CM(3)
C
C                                  Lorentz Transformation to the Lab Frame
C
C                   ETA  CM Parameters
C
      GAMMA = EPI0LF/EMPI0
      BETAX = PPI0LF(1)/EPI0LF
      BETAY = PPI0LF(2)/EPI0LF
      BETAZ = PPI0LF(3)/EPI0LF
C
C                                          For First Gamma
C
      SCALAR= BETAX*PG1CM(1)+BETAY*PG1CM(2)+BETAZ*PG1CM(3)
C
      EG1LF    = GAMMA*(EG1CM + SCALAR)
C
      TERM     = GAMMA*(EG1CM+GAMMA*SCALAR/(GAMMA+1.))
C
      PG1LF(1) =  PG1CM(1) + BETAX*TERM
      PG1LF(2) =  PG1CM(2) + BETAY*TERM
      PG1LF(3) =  PG1CM(3) + BETAZ*TERM
C
C                                          For Second Gamma
C
      SCALAR= BETAX*PG2CM(1)+BETAY*PG2CM(2)+BETAZ*PG2CM(3)
C
      EG2LF    = GAMMA*(EG2CM + SCALAR)
C
      TERM     = GAMMA*(EG2CM+GAMMA*SCALAR/(GAMMA+1.))
C
      PG2LF(1) =  PG2CM(1) + BETAX*TERM
      PG2LF(2) =  PG2CM(2) + BETAY*TERM
      PG2LF(3) =  PG2CM(3) + BETAZ*TERM
C
C      Opening Angle betwean g1 and g2 in Lab fram
C
      PG1PG2 = PG1LF(1)*PG2LF(1)+PG1LF(2)*PG2LF(2)+PG1LF(3)*PG2LF(3)
      CSPSI  = PG1PG2/(EG1LF*EG2LF)
      PSIG12 = RADDEG*ACOS(CSPSI)        ! in degrees
C
      CALL HF1(7,PSIG12,1.)
C
      RETURN
      END
C





























