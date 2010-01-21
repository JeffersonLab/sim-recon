C
C**************************************************************************
C                                                                         *
C-- Date: 10/18/09                                                        *
C-- Author : A. Gasparian                                                 *
C                                                                         *
C      INTERACTEV VERSION OF                                              *
C......GEANT SIMULATION PROG. FOR HALL D Real Promakoff Experiment      *
C.                                               GEANT VERSION   3.21     *
C.                                               PROG. VERSION 1.01     *
C              Eta Primakoff generator for Hall D project                 *
C                                                                         *
C.*************************************************************************
C.....
C.
C
C****************************************************************************
*-- Author :A. Gasparian                                                    *
*-- DATE   :08/29/97                                                         *
*                                                                           *
*                                                                           *
*                                                                           *
C.               To open FFREAD and HBOOK files                             *
*                                                                           *
*****************************************************************************
C
      SUBROUTINE UFILES
C
      OPEN(UNIT=4,FILE='eta_p_gen.dat',STATUS='OLD')
C
C       the file below is for Pawl to generate events in Matt's format
CDL      open(unit=11, file='Prim_eta_Prt.dat', status='new')
C
      RETURN
      END
C.
C.....
C.
C
C$LIST ON
C.*************************************************************************
*-- Author :A. Gasparian                                                  *
*-- DATE   :9/04/97                                                       *
*                                                                         *
C.    *                                                                   *
C.    *                                                                   *
C.    *        To initialise GEANT  prog. and read data cards           *
C.    *                                                                   *
*                                                                         *
C.*************************************************************************
C
      SUBROUTINE UGINIT
C.
C.
      COMMON/GCFLAG/IDEBUG,IDEMIN,IDEMAX,ITEST,IDRUN,IDEVT,IEORUN
     +        ,IEOTRI,IEVENT,ISWIT(10),IFINIT(20),NEVENT,NRNDM(2)
      COMMON/GCUNIT/LIN,LOUT,NUNITS,LUNITS(5)
      INTEGER LIN,LOUT,NUNITS,LUNITS
      COMMON/GCMAIL/CHMAIL
      CHARACTER*132 CHMAIL
C
      COMMON/GCLIST/NHSTA,NGET ,NSAVE,NSETS,NPRIN,NGEOM,NVIEW,NPLOT
     +       ,NSTAT,LHSTA(20),LGET (20),LSAVE(20),LSETS(20),LPRIN(20)
     +             ,LGEOM(20),LVIEW(20),LPLOT(20),LSTAT(20)
C
      COMMON/GCKINE/IKINE,PKINE(10),ITRA,ISTAK,IVERT,IPART,ITRTYP
     +      ,NAPART(5),AMASS,CHARGE,TLIFE,VERT(3),PVERT(4),IPAOLD
C
      PARAMETER (KWBANK=99000,KWWORK=8200)
      COMMON/GCBANK/NZEBRA,GVERSN,ZVERSN,IXSTOR,IXDIV,IXCONS,FENDQ(16)
     +             ,LMAIN,LR1,WS(KWBANK)
      COMMON/GCLINK/JDIGI ,JDRAW ,JHEAD ,JHITS ,JKINE ,JMATE ,JPART
     +      ,JROTM ,JRUNG ,JSET  ,JSTAK ,JGSTAT,JTMED ,JTRACK,JVERTX
     +      ,JVOLUM,JXYZ  ,JGPAR ,JGPAR2,JSKLT
C
      COMMON/PHBMEN/EGFRAC(2)
      COMMON/EQUIPH/QEQUPH
      COMMON/EQHORS/THOURS
      COMMON/TMCAR4/egammn,egammx,denest
      DOUBLE PRECISION egammx,egammn,denest
C
      COMMON/TMCAR2/CRSSUM,CRSOME
      COMMON/TMCAR3/YIELD
      COMMON/TMCAR9/totpht
C
      INTEGER t
C
      DIMENSION PAR(8)
C
C             Open user files
C
      CALL UFILES
C
C             Initialise GEANT
C
      CALL GINIT
C
      CALL FFKEY('TSAM',THOURS,1,'REAL')
C                                          for 45 days on LH2
Cn      THOURS = 1080.000
C                                Tagged beam equivalent Photon Number in sec
      CALL FFKEY('QEQU',QEQUPH,1,'REAL')
C                                        Tagged gamma accepted fractions
      CALL FFKEY('EGFR',EGFRAC,2,'REAL')
C
      CALL FFSET('LINP',4)
      CALL GFFGO
      CLOSE(4)
C
C       rundom namber generater initialization
C
      t = TIME()
      CALL RLUXGO(4,t,0,0)
C
C
C
C
C........      Initalization of data structure to be used for the eta polar angle sampling
C              vs. cross sections. 
C
      write(6,*)'?????????????????????????????? CALL YILINT ????'
      CALL YILINT
C
      write(6,*)'???????????????????????????? FINISH CALL YILINT ????'
C
C                                        from 120708
C120708      for the 30 cm He4 target
C      N(He4) = 5.644x10^23 /cm^2
C      N(hour)= 3600*1.*E+7x5.644*E+23{1.E-6*YIELD in mb*Nphotons}*1.E-27
C             = 20.3184
C
C121408      EVENTH= 20.3184*THOURS*YIELD
C                            Br. Ratio eta -- g+g    39.25+/- 0.31 %
C121408      EVETGG = EVENTH*0.3925
C
C120708      for the 30 cm LH2 target,  work below
C      N(H) = 1.28x10^24 /cm^2
C      N(hour)= 3600*1.*E+7x12.8*E+23{1.E-6*YIELD in mb*Nphotons}*1.E-27
C             = 46.080
C
      EVENTH= 46.080*THOURS*YIELD
C                            Br. Ratio eta -- g+g    39.25+/- 0.31 %
      EVETGG = EVENTH*0.3925
C
      write(6,*)'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(6,*)'>>>>  Integral (Cr.Sect.*dteta) = ',CRSSUM,' mb/sr'
      write(6,*)'>>>>  Integral (Cr.Sect.*domega) = ',CRSOME/1.E+6,' mb'
      write(6,*)'>>>>  For Time= ',THOURS,' hours'
C      write(6,*)'>>>>  For 5.644 x 10^23 he4/cm^2 target '
      write(6,*)'>>>>  For 12.8 x 10^23 proton/cm^2 target, (30 cm) '
      Eg1=SNGL(egammn)
      Eg2=SNGL(egammx)
      write(6,*)'>>>>  For Eg int from ',Eg1,' to ',Eg2, ' GeV'
C      write(6,*)'>>>>  and for ',QEQUPH,' microA  el. beam '
      write(6,*)'>>>>  and for ',QEQUPH,' eq. ph. beam '
      write(6,*)'>>>>  total photons = ', totpht,' x10^7'
      write(6,*)'>>>>  You must sample TOTAL ETA EVENT = ',EVENTH
      write(6,*)'>>>>  OR must sample TOTAL eta - g+g EVENT = ',EVETGG
      write(6,*)'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
C
C
      PRINT*,'------------ Calling GZINIT'
      CALL GZINIT
C      PRINT*,'------------ Calling GDINIT'
C      CALL GDINIT
      PRINT*,'------------ Calling GPART'
      CALL GPART
      PRINT*,'------------'
C
C             Prints version number
C
      WRITE(LOUT,1000)
C
C
C              Geometry and materials description
C
      CALL GMATE           ! geant default materials  #1-16

C
C ****************   Defines  tracking media parameters        **********
C.
      NMED1 = 1                ! AIR, NOT OPTIC               # 1 
C
C...       Tracking medium parameters for detectors
C.
C                           !! negative value of parameters meens that
C                              even in case IGAUTO=0, automatic calculation
C                              still takes place
      IFIELD =  0
      FIELDM =  0.
      TMAXFD =  10.
      STEMAX =  10.
      DEEMAX =  0.1             ! FROM GEXAM1
      EPSIL  = 0.001
      STMIN  = 0.01             ! IN GEXAM8 = -0.01
C.
      ISVOL = 1                   ! VOLUME IS SENSITIVE
C.
      CALL GSTMED( 1,'DEFAULT MEDIUM AIR$ ' ,15 ,  0 , IFIELD,
     *             FIELDM,TMAXFD,STEMAX,DEEMAX, EPSIL, STMIN, 0 , 0 )
C......................................................................
C
C.
C...  PARAMETERS OF MOTHER VOLUME 'HALD'
C... Hall D is described as a 'BOX ' with half sizes =1320.x2610.x2700. cm^3
C... hall D referense frame
C... OZ ALONG BEAM, OX HORIZONTAL(NON DISPERSION), OY- VERTICAL 
C
      PAR(1) =  500.000
      PAR(2) =  500.000
      PAR(3) =  2000.000
C                 Temporarily Hall D is an air #1, not mag.
C
      CALL GSVOLU( 'HALD','BOX ', 1,PAR,3,IVOL)               !  HALL D
C
C                      Close geometry banks
C
      CALL GGCLOS
C.
C........................................................
C
      CALL GLOOK('MATE',LPRIN,NPRIN,IM)
      CALL GLOOK('TMED',LPRIN,NPRIN,IT)
      CALL GLOOK('VOLU',LPRIN,NPRIN,IV)
      IF(IM.NE.0)CALL GPRINT('MATE',0)
      IF(IT.NE.0)CALL GPRINT('TMED',0)
      IF(IV.NE.0)CALL GPRINT('VOLU',0)
C
C              Energy loss and cross-sections initialisations
C
      CALL GPHYSI
C
C
C             Define  histograms
C
CDL      CALL UHINIT
C
 1000 FORMAT(/,'R. ETA Pr. VERS 0.01, :10/18/09, by A. Gasparian',/)
C
      RETURN
      END
C.
C....
C.
C
C
C.*************************************************************************
*     Author :A. Gasparian                                                *
*     date : 10/18/09                                                     *
*                                                                         *
*.           Sampling of Real Primakoff Events                      *
*                                                                         *
C.*************************************************************************
C.
      SUBROUTINE YILINT
C
C
      COMMON/GCKINE/IKINE,PKINE(10),ITRA,ISTAK,IVERT,IPART,ITRTYP
     +      ,NAPART(5),AMASS,CHARGE,TLIFE,VERT(3),PVERT(4),IPAOLD
C
      COMMON/PHBMEN/EGFRAC(2)
      COMMON/EQUIPH/QEQUPH
C
      COMMON/TMCAR1/FTPNOR(1:150,0:10000),tpi0ms(10000)
C
      COMMON/TMCAR2/CRSSUM,CRSOME
      COMMON/TMCAR3/YIELD
      COMMON/TMCAR4/egammn,egammx,denest
      COMMON/TMCAR9/totpht
C
C
      DOUBLE PRECISION egammx,egammn,denest
C
      DOUBLE PRECISION PI,TWOPI,PIBY2,DEGRAD
C
      PARAMETER (PI=3.141592653)
      PARAMETER (TWOPI=6.2831853071)
      PARAMETER (PIBY2=1.57079632679)
      PARAMETER (DEGRAD=0.01745329251)
C 
      nenedv = 150
C
      egammx = DBLE(EGFRAC(2))*DBLE(PKINE(1))
      egammn = DBLE(EGFRAC(1))*DBLE(PKINE(1))
C
      denest = (egammx - egammn)/DBLE(nenedv)
C
      YIELD = 0.000000
C
      do i = 1, nenedv
C
      engamk = SNGL(egammn)+i*SNGL(denest)-SNGL(denest/2.D0)
C
      CALL TOCRS2(engamk)
C
C
C from subr TOCRS2 we now have CRSOME = integr. crsec*domega in 0-? degree
C for the carrent g energy => engamk
C Now we have to calculate Number of Photons in the [(engamk-denest),engamk]
C   see page #36 book #2
C  QEQUPH is the number of equivalent photons in second in units 1.E+7
C
C
      phnumb = QEQUPH*(DLOG(DBLE(engamk))-DLOG(DBLE(engamk-denest)))
C
      totpht= totpht+phnumb
C
      YIELD = YIELD+phnumb*CRSOME
C
C                    new
      ienbin = i
      CALL FSUMMC(ienbin,engamk)
C
C      write(6,*)'?? in YILINT  i = ',i,' engamk=',engamk
C
      enddo
C
C
      RETURN
      END
C
C
C.....
C
C
C.*************************************************************************
*     Author :A. Gasparian                                                *
*     date : 10/20/97                                                     *
*                                                                         *
*.           For sampling of Real Primakoff Events                      *
*    now for each event we prepareing unique FTPNOR(i) for Thetapi0 sampl.*
*                                                                         *
C.*************************************************************************
C.
      SUBROUTINE FSUMMC(ienbin,engamk)
C
C
      COMMON/GCKINE/IKINE,PKINE(10),ITRA,ISTAK,IVERT,IPART,ITRTYP
     +      ,NAPART(5),AMASS,CHARGE,TLIFE,VERT(3),PVERT(4),IPAOLD
C
C
      COMMON/TMCAR1/FTPNOR(1:150,0:10000),tpi0ms(10000)
C
      DOUBLE PRECISION FTPROB(0:10000)
      DOUBLE PRECISION tpi0mx,tpi0mn,dpi0st,tpi0k
      DOUBLE PRECISION UPCOS,DNCOS,DEOMEG
C
      DOUBLE PRECISION PI,TWOPI,PIBY2,DEGRAD
C
      PARAMETER (PI=3.141592653)
      PARAMETER (TWOPI=6.2831853071)
      PARAMETER (PIBY2=1.57079632679)
      PARAMETER (DEGRAD=0.01745329251)
C 
C
      ntetdv = 10000
      tpi0mx = DBLE(PKINE(9))                  ! degree
      tpi0mn = 0.000000D0                      ! degree
C
      dpi0st = (tpi0mx - tpi0mn)/DBLE(ntetdv)
C
      FTPROB(0) = 0.000000D0
C
      do i = 1, ntetdv
C
         tpi0k = tpi0mn + DBLE(i)*dpi0st-dpi0st/2.D0
C
         tpi0ms(i) = SNGL(tpi0k)
C
C           here, energy of photon comes from YILINT, for each bin
C
         CALL CRSEC6(engamk,tpi0k,sumall)
C
C
         UPCOS = DCOS(DEGRAD*tpi0k)
         DNCOS = DCOS(DEGRAD*(tpi0k-dpi0st))
C
         DEOMEG=2.D0*PI*(DNCOS-UPCOS)
C
         FTPROB(i) = FTPROB(i-1) + DBLE(sumall)*DEOMEG   ! in mbarn
C
      enddo
C
      do i = 1, ntetdv
C 
         FTPNOR(ienbin,i) = SNGL(FTPROB(i)/FTPROB(ntetdv))
C
      enddo
C
C
C
      RETURN
      END
C
C
C.....
C
CC$LIST ON
C
C.*************************************************************************
*     Author :A. Gasparian                                                *
*     date : 9/19/97                                                      *
*                                                                         *
*.           Total sampling of Real Primakoff Events                      *
*                                                                         *
C.*************************************************************************
C.
      SUBROUTINE TOCRS2(engamk)
C
C
      COMMON/GCKINE/IKINE,PKINE(10),ITRA,ISTAK,IVERT,IPART,ITRTYP
     +      ,NAPART(5),AMASS,CHARGE,TLIFE,VERT(3),PVERT(4),IPAOLD
C
C
      COMMON/TMCAR1/FTPNOR(1:150,0:10000),tpi0ms(10000)
C
      COMMON/TMCAR2/CRSSUM,CRSOME
C
      DOUBLE PRECISION tpi0mx,tpi0mn,dpi0st,tpi0k
      DOUBLE PRECISION UPCOS,DNCOS,DEOMEG
C
      DOUBLE PRECISION PI,TWOPI,PIBY2,DEGRAD
C
      PARAMETER (PI=3.141592653)
      PARAMETER (TWOPI=6.2831853071)
      PARAMETER (PIBY2=1.57079632679)
      PARAMETER (DEGRAD=0.01745329251)
C
C
      ntetdv = 10000
      tpi0mx = DBLE(PKINE(9))         ! degree
      tpi0mn = 0.000000D0             ! degree
C
      dpi0st = (tpi0mx - tpi0mn)/DBLE(ntetdv)
C
      CRSSUM = 0.000000
      CRSOME = 0.000000
C
      do i = 1, ntetdv
C
      tpi0k = tpi0mn + DBLE(i)*dpi0st-dpi0st/2.D0
C
C         new total cross section subroutine, with variable Egamma
C
      CALL CRSEC6(engamk,tpi0k,sumall)
C
C
      CRSSUM = CRSSUM + sumall*SNGL(dpi0st)               ! in milib
C
C   here I add a term 1.D6 to the domega for accuracy, it will be taken
C   in account in the call-ing subroutine
C
      UPCOS = 1.D6*DCOS(DEGRAD*(tpi0k+dpi0st/2.D0))
      DNCOS = 1.D6*DCOS(DEGRAD*(tpi0k-dpi0st/2.D0))
C
      DEOMEG=2.D0*PI*(DNCOS-UPCOS)
C
      CRSOME=CRSOME+sumall*SNGL(DEOMEG)                   ! in milibarn
C
      enddo
C
C
      RETURN
      END
C
C
C.....
C
C
C$LIST OFF
C
C.*************************************************************************
*-- Author :A. Gasparian                                                  *
*-- DATE   :9/02/97                                                       *
*                                                                         *
C     *                                                                   *
C     ***         To book  Histograms for EVENT development               *
C     *                                                                   *
*                                                                         *
C.*************************************************************************
C
      SUBROUTINE UHINIT
C
      COMMON/GCLINK/JDIGI ,JDRAW ,JHEAD ,JHITS ,JKINE ,JMATE ,JPART
     +      ,JROTM ,JRUNG ,JSET  ,JSTAK ,JGSTAT,JTMED ,JTRACK,JVERTX
     +      ,JVOLUM,JXYZ  ,JGPAR ,JGPAR2,JSKLT
C
      COMMON/GCFLAG/IDEBUG,IDEMIN,IDEMAX,ITEST,IDRUN,IDEVT,IEORUN
     +        ,IEOTRI,IEVENT,ISWIT(10),IFINIT(20),NEVENT,NRNDM(2)
      COMMON/GCFLAX/BATCH, NOLOG
      LOGICAL BATCH, NOLOG
C
      COMMON/THARGT/THRPOS(3)
C
      COMMON/GCKINE/IKINE,PKINE(10),ITRA,ISTAK,IVERT,IPART,ITRTYP
     +      ,NAPART(5),AMASS,CHARGE,TLIFE,VERT(3),PVERT(4),IPAOLD
C
      CHARACTER*6 COORDN(8)
C
      CHARACTER*15 histfl
C
      CHARACTER*2 subnam(100)
      DATA subnam /'01','02','03','04','05','06','07','08','09','10',
     +'11','12','13','14','15','16','17','18','19','20',
     +'21','22','23','24','25','26','27','28','29','30',
     +'31','32','33','34','35','36','37','38','39','40',
     +'41','42','43','44','45','46','47','48','49','50',
     +'51','52','53','54','55','56','57','58','59','60',
     +'61','62','63','64','65','66','67','68','69','70',
     +'71','72','73','74','75','76','77','78','79','80',
     +'81','82','83','84','85','86','87','88','89','90',
     +'91','92','93','94','95','96','97','98','99','XX'/
C
C below are only for the reference
C
      DATA COORDN/'X_g1  ', 'Y_g1  ','X_g2  ', 'Y_g2  ',
     +'invmft', 'thetft', 'thetex', 'invmas'/
C
Ct     +'Tkirec','PL1G1X','PL1G1Y','PL1G2X','PL1G2Y','G1_reg','G2_reg',
Ct     +'Eg1_LF','Eg2_LF','X_g1  ','Y_g1  ','Z_g1  ','X_g2  ','Y_g2  ',
Ct     +'Z_g2  ','t_reco','efm_pi','Texppi','Fexppi',
Ct     +'md#_g1','md#_g2','Psi_al','Psi_ex','Eg1_ex','Eg2_ex'/
C
      histfl='Gen_prt_'//subnam(IDRUN)//'.hist'
C
      CALL HROPEN(1,'HBOOK',histfl,'N',1024,ISTAT)
C
      IF(ISTAT.NE.0) GO TO 99
C
      CALL HBOOKN(101,'COOR. IN PL.',8,'//HBOOK',9950,COORDN)   ! Ntupel #101
C

      CALL HBOOK1(1,'Incident g-beam Energy BINS(GeV)$'
     *,100,0.0,12.0, 0.0)
      CALL HBOOK1(2,'Eta Polar Angle  (deg.)$'
     *,100,-0.1,3.9, 0.0)
C
      CALL HBOOK1(21,'Eta Solid Angle  (deg.)$'
     *,100,-0.1,3.9, 0.0)
C
      CALL HBOOK1(3,'Eta Azimutal Angle  (deg.)$'
     *,100,-10.0,370.0, 0.0)
      CALL HBOOK1(4,'Eta Energy  (GeV)$'
     *,100,6.0,12.0, 0.0)
      CALL HBOOK1(5,'Recoil Nucleus Polar Angle  (deg.)$'
     *,100,0.0,100.0, 0.0)
      CALL HBOOK1(6,'Recoil Nucleus Kinetic Energy (MeV)$'
     *,100,-0.1,19.9, 0.0)
      CALL HBOOK1(7,'Two gamma Opening Angle (Degr.)$'
     *,100,0.0,60.0, 0.0)
C
      CALL HBOOK1(11,'X-position of the vertex on Target$'
     *,400,-2.0,2.0, 0.0)
      CALL HBOOK1(12,'Y-position of the vertex on Target$'
     *,400,-2.0,2.0, 0.0)
      ZCENTER =THRPOS(3)+20.00
      CALL HBOOK1(13,'Z-position of the vertex on Target$'
     *,100,-(ZCENTER-50.0),(ZCENTER+50.0), 0.0)
C
      CALL HBOOK1(31,'Random number gener               $'
     *,100,-3.0           ,3.0           , 0.0)
C
      CALL HBOOK1(41,'Scattered-g1, PgX           (GeV)$'
     *,100,0.0,12.0, 0.0)
      CALL HBOOK1(42,'Scattered-g1, PgY           (GeV)$'
     *,100,0.0,12.0, 0.0)
      CALL HBOOK1(43,'Scattered-g1, PgZ           (GeV)$'
     *,100,0.0,12.0, 0.0)
C
      CALL HBOOK1(44,'Scattered-g2, PgX           (GeV)$'
     *,100,0.0,12.0, 0.0)
      CALL HBOOK1(45,'Scattered-g2, PgY           (GeV)$'
     *,100,0.0,12.0, 0.0)
      CALL HBOOK1(46,'Scattered-g2, PgZ           (GeV)$'
     *,100,0.0,12.0, 0.0)
C
        CALL HBIGBI(0,4)
C
 99   CONTINUE
C
      RETURN
      END
C.
C..........
C.
C
C.************************************************************************
*__ Author: A. Gasparian                                                 *
*-- DATE : 9/09/94                                                       *
*                                                                        *
C.  *       Routine to generate one hadronic interaction                 *
C.  *                                                                    *
C.  *    ==>Called by :GTHADR,GTNEUT                                     *
C.  *                                                                    *
*                                                                        *
C.************************************************************************
C
      SUBROUTINE GUHADR
C.
C.
      COMMON/GCPHYS/IPAIR,SPAIR,SLPAIR,ZINTPA,STEPPA
     +             ,ICOMP,SCOMP,SLCOMP,ZINTCO,STEPCO
     +             ,IPHOT,SPHOT,SLPHOT,ZINTPH,STEPPH
     +             ,IPFIS,SPFIS,SLPFIS,ZINTPF,STEPPF
     +             ,IDRAY,SDRAY,SLDRAY,ZINTDR,STEPDR
     +             ,IANNI,SANNI,SLANNI,ZINTAN,STEPAN
     +             ,IBREM,SBREM,SLBREM,ZINTBR,STEPBR
     +             ,IHADR,SHADR,SLHADR,ZINTHA,STEPHA
     +             ,IMUNU,SMUNU,SLMUNU,ZINTMU,STEPMU
     +             ,IDCAY,SDCAY,SLIFE ,SUMLIF,DPHYS1
     +             ,ILOSS,SLOSS,SOLOSS,STLOSS,DPHYS2
     +             ,IMULS,SMULS,SOMULS,STMULS,DPHYS3
     +             ,IRAYL,SRAYL,SLRAYL,ZINTRA,STEPRA
C.
C.
C.
C          GHEISHA only if IHADR<3   (default)
C          GHEISHA and HADRIN/NUCRIN if IHADR=3
C
      IF (IHADR.NE.4) THEN
         CALL GHEISH
      ELSE
         CALL FLUFIN
      ENDIF
      END
C
C.
C.....
C.
C
C.***************************************************************************
*__ Autor: A. Gasparian                                                     *
*-- DATE : 7/05/94                                                          *
*                                                                           *
C.    *                                                                     *
C.    *   Routine to compute Hadron  inter-n probabilities                  *
C.    *                                                                     *
C.    * ==>Called by : GTHADR,GTNEUT                                        *
C.    *                                                                     *
*                                                                           *
C.***************************************************************************
C
      SUBROUTINE GUPHAD
C.
C.
      COMMON/GCPHYS/IPAIR,SPAIR,SLPAIR,ZINTPA,STEPPA
     +             ,ICOMP,SCOMP,SLCOMP,ZINTCO,STEPCO
     +             ,IPHOT,SPHOT,SLPHOT,ZINTPH,STEPPH
     +             ,IPFIS,SPFIS,SLPFIS,ZINTPF,STEPPF
     +             ,IDRAY,SDRAY,SLDRAY,ZINTDR,STEPDR
     +             ,IANNI,SANNI,SLANNI,ZINTAN,STEPAN
     +             ,IBREM,SBREM,SLBREM,ZINTBR,STEPBR
     +             ,IHADR,SHADR,SLHADR,ZINTHA,STEPHA
     +             ,IMUNU,SMUNU,SLMUNU,ZINTMU,STEPMU
     +             ,IDCAY,SDCAY,SLIFE ,SUMLIF,DPHYS1
     +             ,ILOSS,SLOSS,SOLOSS,STLOSS,DPHYS2
     +             ,IMULS,SMULS,SOMULS,STMULS,DPHYS3
     +             ,IRAYL,SRAYL,SLRAYL,ZINTRA,STEPRA
C.
      IF (IHADR.NE.4) THEN
         CALL GPGHEI
      ELSE
         CALL FLDIST
      ENDIF
      END
C
C
C.....
C.
C
C.*************************************************************************
*____ Author :A. Gasparian                                                *
*____ date : 8/23/94                                                      *
*                                                                         *
*.           Generates Kinematics for primary track                       *
*                                                                         *
C.*************************************************************************
C.
      SUBROUTINE GUKINE
C.
      COMMON/GCLINK/JDIGI ,JDRAW ,JHEAD ,JHITS ,JKINE ,JMATE ,JPART
     +      ,JROTM ,JRUNG ,JSET  ,JSTAK ,JGSTAT,JTMED ,JTRACK,JVERTX
     +      ,JVOLUM,JXYZ  ,JGPAR ,JGPAR2,JSKLT
C
      COMMON/GCFLAG/IDEBUG,IDEMIN,IDEMAX,ITEST,IDRUN,IDEVT,IEORUN
     +        ,IEOTRI,IEVENT,ISWIT(10),IFINIT(20),NEVENT,NRNDM(2)
      COMMON/GCFLAX/BATCH, NOLOG
      LOGICAL BATCH, NOLOG
C
      COMMON/GCKINE/IKINE,PKINE(10),ITRA,ISTAK,IVERT,IPART,ITRTYP
     +      ,NAPART(5),AMASS,CHARGE,TLIFE,VERT(3),PVERT(4),IPAOLD
C
      COMMON/KINEM1/EINI,TPI0,FIPI0G,EPI0SP,TRECSP,TKINRM
      COMMON/KINEM2/PPREL(3)
      COMMON/KINEM3/EPI0LF,PPI0LF(3),EG1LF,PG1LF(3),EG2LF,PG2LF(3)
C
      COMMON/PHBMEN/EGFRAC(2)
      COMMON/EQUIPH/QEQUPH
C
      COMMON/BRMST1/egamak
C
      DOUBLE PRECISION PI,TWOPI,PIBY2,DEGRAD,RADDEG,CLIGHT,BIG,EMASS
      DOUBLE PRECISION EMMU,PMASS,AVO
C
      PARAMETER (PI=3.141592653)
      PARAMETER (TWOPI=6.2831853071)
      PARAMETER (PIBY2=1.57079632679)
      PARAMETER (DEGRAD=0.01745329251)
      PARAMETER (RADDEG=57.2957795130)
      PARAMETER (CLIGHT=29979245800.)
      PARAMETER (BIG=10000000000.)
      PARAMETER (EMASS=0.0005109990615)
      PARAMETER (EMMU=0.105658387)
      PARAMETER (PMASS=0.9382723128)
      PARAMETER (AVO=0.60221367)
C
C
      DIMENSION VERTEX(3),PLAB(3)
      SAVE VERTEX,PLAB
      DATA VERTEX/3*0./
      DATA PLAB  /3*0./
C
C
C           Kinematics for Real eta Primakoff
C
      IK=IKINE
C
C 
C  PKINE(1) => is the initial Electron beam energy in GeV 
C  PKINE(9) => is the Eta meson maximum Polar angle in degees
C  First one here is IKINE => the type of the particle ( be careful )
C 
C                sampling of photon energy egamak, for each event
C
      CALL BRENER
C
C                With the initial photon energy the eta polar angle will be sampled
C                and all kinematical variables including the eta -> g+g decay channel
C                will be calculated.
C
      CALL PR_KIN
C
C
C      Dave,  here:  EPI0LF,PPI0LF(3) are the eta energy and three momenta 
C                    EG1LF,PG1LF(3)  energy and three momenta of the 1th gamma
C                    EG2LF,PG2LF(3)  energy and three momenta of the 2th gamma
C
C
C                 below are the X, Y and Z widths of the beam (are +/-)
C
      xwidth = 0.0
      ywidth = 0.0
      zwidth = 0.0
C
Cn      CALL RANLUX(RNDM,2)
C
      VERTEX(1) = 0.000
      VERTEX(2) = 0.000
      VERTEX(3) = 50.000
C
Cn      CALL RANLUX(RNDM,2)
C
Cn      CALL HF1(11,VERTEX(1),1.)
Cn      CALL HF1(12,VERTEX(2),1.)
Cn      CALL HF1(13,VERTEX(3),1.)
C                                        First Gamma  from eta -> g + g
C
      CALL GSVERT(VERTEX,0,0,0,0,NVTX)
C
      CALL GSKINE(PG1LF, 1 ,NVTX,0,0,NT)
C
C                                       Second Gamma from  pi0 -> g + g
C
      CALL GSVERT(VERTEX,NT,NT,0,0,NVTX)
C
      CALL GSKINE(PG2LF, 1 ,NVTX,0,0,NT)
C
      CALL HF1(41,PG1LF(1),1.)
      CALL HF1(42,PG1LF(2),1.)
      CALL HF1(43,PG1LF(3),1.)
      CALL HF1(44,PG2LF(1),1.)
      CALL HF1(45,PG2LF(2),1.)
      CALL HF1(46,PG2LF(3),1.)
C    ........................................................
C                                  writing events for Pawel like Matt, 10/20/2009
C
C                    IDRUN is the run mumber from *.dat file
      Id_run = IDRUN
      loop = IEVENT
      Npartl = 2
C
CDL        write(11,1)Id_run, loop, Npartl
CDL 1      format(I5,1x,I8,1x,I8)
C
C                      the 1th particle: first gama from eta --> g+g

CDL        write(11,6)
CDL 6      format('1', ' 1', 1x, '0.000')
CDL        write(11,7)PG1LF(1),PG1LF(2),PG1LF(3), EG1LF
CDL 7      format(' 0',1x, F10.6,1x,F10.6,1x,F10.6, 1x, F10.6)
C
C                     the 2th particle: second decay g from eta-->g+g
C
CDL        write(11,8)
CDL 8      format('2', ' 1', 1x, '0.000')
CDL        write(11,9)PG2LF(1),PG2LF(2),PG2LF(3),EG2LF
CDL 9      format(' 0',1x, F10.6,1x,F10.6,1x,F10.6, 1x, F10.6)

C     ................................................... end of Matt
C ***          Kinematics debug (controlled by ISWIT(1) )
C
      IF(IDEBUG.EQ.1) THEN
        IF(ISWIT(1).EQ.1) THEN
          CALL GPRINT('VERT',0)
          CALL GPRINT('KINE',0)
        ENDIF
      ENDIF
C.
C.
      RETURN
      END
C
C.
C.....
C.
C
C.*************************************************************************
*____ Author :A. Gasparian                                                *
*____ date : 10/21/97                                                     *
*                                                                         *
*.        Samples gamma with energy according bremsstrahlung beam         *
*         assuming Bethe Gaitler form of Energy                           *
*                                                                         *
C.*************************************************************************
C.
      SUBROUTINE BRENER
C
C
      COMMON/GCKINE/IKINE,PKINE(10),ITRA,ISTAK,IVERT,IPART,ITRTYP
     +      ,NAPART(5),AMASS,CHARGE,TLIFE,VERT(3),PVERT(4),IPAOLD
C
      COMMON/PHBMEN/EGFRAC(2)
      COMMON/EQUIPH/QEQUPH
C
      COMMON/TMCAR4/egammn,egammx,denest
      COMMON/BRMST1/egamak
C
      DOUBLE PRECISION egammx,egammn,denest
C
C                         PKINE(1) => is the initial Electron beam energy in GeV 
C
      CALL RANLUX(RTL,1)
C
C               see page # 38 in the book #2
C
      DRL = DBLE(RTL)
      egamak=SNGL(DEXP(DRL*(DLOG(egammx)-DLOG(egammn))+DLOG(egammn)))
C
C
C
      RETURN
      END
C
C.
C.....
C.
C
C.************************************************************************
*____ Author :A. Gasparian                                               *
*____ date : 8/23/94                                                     *
*                                                                        *
C.    *                                                                  *
C.    *       Routine to control tracking of one event                   *
C.    *                                                                  *
C.    *       Called by GRUN                                             *
C.    *                                                                  *
*.************************************************************************
C
      SUBROUTINE GUTREV
C.
      COMMON/GCFLAG/IDEBUG,IDEMIN,IDEMAX,ITEST,IDRUN,IDEVT,IEORUN
     +        ,IEOTRI,IEVENT,ISWIT(10),IFINIT(20),NEVENT,NRNDM(2)
      COMMON/GCKINE/IKINE,PKINE(10),ITRA,ISTAK,IVERT,IPART,ITRTYP
     +      ,NAPART(5),AMASS,CHARGE,TLIFE,VERT(3),PVERT(4),IPAOLD
      COMMON/ANGL12/SNTEL,CSTEL,SNTHD,CSTHD
      COMMON/ANGL13/SNTET,CSTET,SPMSET
C
C.
C.
      CALL GTREVE
C.
      RETURN
      END
C.
C
C......
C.
C.************************************************************************
*____ Author :A. Gasparian                                               *
*____ date : 7/12/94                                                     *
*                                                                        *
C.    *                                                                  *
C.    *       This routine called at the end of each tracking step       *
C.    *       INWVOL is different from 0 when the track has reached      *
C.    *       ISTOP is different from 0 if the track has stopped         *
C.    *                                                                  *
C.                                                                       *
C.************************************************************************
C
      SUBROUTINE GUSTEP
C
      COMMON/GCLINK/JDIGI ,JDRAW ,JHEAD ,JHITS ,JKINE ,JMATE ,JPART
     +      ,JROTM ,JRUNG ,JSET  ,JSTAK ,JGSTAT,JTMED ,JTRACK,JVERTX
     +      ,JVOLUM,JXYZ  ,JGPAR ,JGPAR2,JSKLT
      COMMON/GCFLAG/IDEBUG,IDEMIN,IDEMAX,ITEST,IDRUN,IDEVT,IEORUN
     +        ,IEOTRI,IEVENT,ISWIT(10),IFINIT(20),NEVENT,NRNDM(2)
      COMMON/GCFLAX/BATCH, NOLOG
      COMMON/GCVOLU/NLEVEL,NAMES(15),NUMBER(15),
     +LVOLUM(15),LINDEX(15),INFROM,NLEVMX,NLDEV(15),LINMX(15),
     +GTRAN(3,15),GRMAT(10,15),GONLY(15),GLX(3)
      LOGICAL BATCH, NOLOG
C
      COMMON/GCTMED/NUMED,NATMED(5),ISVOL,IFIELD,FIELDM,TMAXFD,STEMAX
     +      ,DEEMAX,EPSIL,STMIN,CFIELD,PREC,IUPD,ISTPAR,NUMOLD
      COMMON/GCTLIT/THRIND,PMIN,DP,DNDL,JMIN,ITCKOV,IMCKOV,NPCqOV

C
      COMMON/GCKINE/IKINE,PKINE(10),ITRA,ISTAK,IVERT,IPART,ITRTYP
     +      ,NAPART(5),AMASS,CHARGE,TLIFE,VERT(3),PVERT(4),IPAOLD
C
C
      COMMON/GCSETS/IHSET,IHDET,ISET,IDET,IDTYPE,NVNAME,NUMBV(20)
C
      INTEGER MXGKIN
      PARAMETER (MXGKIN=100)
      COMMON/GCKING/KCASE,NGKINE,GKIN(5,MXGKIN),
     +                           TOFD(MXGKIN),IFLGK(MXGKIN)
      INTEGER MXPHOT
      PARAMETER (MXPHOT=800)
      COMMON/GCKIN2/NGPHOT,XPHOT(11,MXPHOT)
      COMMON/GCKIN3/GPOS(3,MXGKIN)
      INTEGER       KCASE,NGKINE ,IFLGK
      INTEGER       NGPHOT
      REAL          GKIN,TOFD,GPOS
C
      PARAMETER (MAXMEC=30)
      COMMON/GCTRAK/VECT(7),GETOT,GEKIN,VOUT(7),NMEC,LMEC(MAXMEC)
     + ,NAMEC(MAXMEC),NSTEP ,MAXNST,DESTEP,DESTEL,SAFETY,SLENG
     + ,STEP  ,SNEXT ,SFIELD,TOFG  ,GEKRAT,UPWGHT,IGNEXT,INWVOL
     + ,ISTOP ,IGAUTO,IEKBIN, ILOSL, IMULL,INGOTO,NLDOWN,NLEVIN
     + ,NLVSAV,ISTORY
C
      COMMON/MAGCUT/CUTMAG(10),MGNSTP
C
      COMMON/GCNUM/NMATE ,NVOLUM,NROTM,NTMED,NTMULT,NTRACK,NPART
     +            ,NSTMAX,NVERTX,NHEAD,NBIT
C
      CHARACTER*4 NAMES
C.
      CALL GDCXYZ                   ! SHOWING THE TRECKS INTERACTIVLY
C.
Cn      CALL ACC
C
C.......,,,,,,,
C.
C.             Debug event
C
CT      IF(IDEBUG.EQ.1) THEN
C
CT 50   CALL GDEBUG
C
CT         IF(ISWIT(1).EQ.1) THEN
CT         CALL GSXYZ
C
CT         ENDIF
CT      ENDIF
C
C
      RETURN
      END
C
C......
C.
C
C.************************************************************************
*____ Autor :A. Gasparian                                                *
*____ data : 7/21/94                                                     *
*                                                                        *
*                                                                        *
C.************************************************************************
C
      SUBROUTINE GUOUT
C
      COMMON/GCFLAG/IDEBUG,IDEMIN,IDEMAX,ITEST,IDRUN,IDEVT,IEORUN
     +        ,IEOTRI,IEVENT,ISWIT(10),IFINIT(20),NEVENT,NRNDM(2)
      COMMON/GCKINE/IKINE,PKINE(10),ITRA,ISTAK,IVERT,IPART,ITRTYP
     +      ,NAPART(5),AMASS,CHARGE,TLIFE,VERT(3),PVERT(4),IPAOLD
      COMMON/GCFLAX/BATCH, NOLOG
      LOGICAL BATCH, NOLOG
C
C
      RETURN
      END
C.
C
C.........
C
C
C
C.**************************************************************************
*____ Author :A. Gasparian                                                 *
*____ date : 8/23/94                                                       *
*                                                                          *
C.    *                                                                    *
C.    *      Termination routine to print histograms and statistics        *
C.    *                                                                    *
*                                                                          *
C.**************************************************************************
C
      SUBROUTINE UGLAST
C
C
      COMMON/GCFLAG/IDEBUG,IDEMIN,IDEMAX,ITEST,IDRUN,IDEVT,IEORUN
     +        ,IEOTRI,IEVENT,ISWIT(10),IFINIT(20),NEVENT,NRNDM(2)
C
      COMMON/ACCCAL/NHOLMC,NOUTMC
      COMMON/ACCCA2/N1gpwo
C
      COMMON/GCFLAX/BATCH, NOLOG
      LOGICAL BATCH, NOLOG
C
      CALL GLAST
C
      CALL HROUT(0,ICYCLE,' ')
      CALL HREND('HBOOK')
C
C         for the file opened for the event generator
         close(unit=11)
C
C
      RETURN
      END
C.
C......
C.
