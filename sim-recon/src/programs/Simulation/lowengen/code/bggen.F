C
C---  Simulates low energy (<3 GeV) photoproduction by a coherent Bremsstrahlung beam
C     Reaction: gamma+p
C     Includes: a) calculation of the coherent+incoherent photon energy spectra
C               b) parametrization for the cross sections:
C                  - full (formula fit to data)
C                  - p pi0, n pi+ - using SAID
C                  - p 2pi, n 2pi, p eta, p 3pi, n 3pi (formula fit to data)
C               c)  simulation if unbiased (equal weight) events in a given beam energy range: 
C                  - beam energy simulated (beam spectrum times the total cross section)
C                  - the process is chosen randomly accordingly to the their probabilities 
C
C---  Input: file "run.input" containing: 
C                      the number of events, the beam energy range and the distance to the collimator
C            file "particle.dat" contains a table for particle masses (GEANT numbering)
C
C
      PROGRAM BGGEN
C
      IMPLICIT NONE
C
      INCLUDE 'bg_partc.inc'
      INCLUDE 'bg_evec.inc'
C
      INTEGER mxpawc
      PARAMETER (mxpawc=1000000)
      COMMON/PAWC/ HMEM(mxpawc)
      REAL HMEM
C
      COMMON/QUEST/ IQUEST(100)
      INTEGER       IQUEST
C
      REAL      enelec,enpeak,beamlim(2),zcollim,am,wg,qq
C
      CHARACTER cline*132
C
      INTEGER lun,lout,lunnt,i,j,nevent,lenc,iost,iev,ip,nb,ipro,iwrnt
     +       ,kd(4),lundat,idd
      INTEGER id,lrec,icycle,idnt
      INTEGER    mxpro
      PARAMETER (mxpro=10)
      INTEGER nevpro(mxpro)
      INTEGER kftr(2,MXTRA)    ! for HDDS
      REAL    pptr(4,MXTRA)    ! for HDDS
      CHARACTER cnampro(mxpro)*16
      DATA cnampro/'p pi0           '
     +            ,'n pi+           '
     +            ,'p pi+ pi-       '
     +            ,'p rho0          '
     +            ,'Delta++ pi-     '
     +            ,'p pi0 pi0       '
     +            ,'n pi+ pi0       '
     +            ,'p eta           '
     +            ,'p pi+ pi- pi0   '
     +            ,'n pi+ pi+ pi-   '
     +            /
C
C
C     ------------------------------------------------------------------
C
      lun=9
      lout=6
      lunnt=2
      lundat=21
      id=10000
C
      CALL HLIMIT(mxpawc)
C
C---   Read the input file
C
      OPEN(lun,FILE='run.input',STATUS='OLD',IOSTAT=iost
     +       ,FORM='FORMATTED')
      IF(iost.NE.0) THEN
         WRITE(lout,*) ' *** ERROR: Missing file run.input'
         GO TO 999
      ENDIF
      nevent=0
 5    READ(lun,FMT='(A)',IOSTAT=iost) cline
      IF(iost.EQ.0) THEN
         IF(cline(1:1).NE.'*'.AND.cline(1:1).NE.'C') THEN
            READ(cline,*) nevent,iwrnt,enelec,enpeak,beamlim,zcollim
         ENDIF
         GO TO 5
      ENDIF
      IF(nevent.EQ.0) THEN
         WRITE(lout,*) ' *** ERROR: Reading file run.input: '
     +            ,'event number=',nevent
         GO TO 999
      ENDIF
      CLOSE(lun)
C
      IF(beamlim(1).LT.0.01.OR.beamlim(2).GT.12..OR.
     +   beamlim(1).GE.beamlim(2)) THEN
         WRITE(lout,*) ' *** ERROR: out of range - beammom=',beamlim
         GO TO 999
      ENDIF
C
C---   Read the particle masses (GEANT numbering)
C
      DO ip=1,MXPART
         AM_PART(ip)=0.
         WG_PART(ip)=0.
         DO i=1,4
            KD_PART(i,ip)=0
         ENDDO
      ENDDO
C
      OPEN(lun,FILE='particle.dat',STATUS='OLD',IOSTAT=iost
     +       ,FORM='FORMATTED')
      IF(iost.NE.0) THEN
         WRITE(lout,*) ' *** ERROR: Missing file particle.dat'
         GO TO 999
      ENDIF
 10   READ(lun,FMT='(A)',IOSTAT=iost) cline
      IF(iost.EQ.0) THEN
C
         IF(cline(1:1).NE.'*'.AND.cline(1:1).NE.'C') THEN
            READ(cline,*) ip,am,wg,kd
C            write(6,*) ip,am,wg,kd
            IF(ip.GT.0.AND.ip.LE.MXPART) THEN
               AM_PART(ip)=am
               WG_PART(ip)=wg
               DO i=1,4
                  KD_PART(i,ip)=kd(i)
               ENDDO
            ELSE
               WRITE(lout,*) ' --- ERROR: Reading file particle.dat ',
     +              'GEANT index is out of range ',ip 
            ENDIF
         ENDIF
C
         GO TO 10
C
      ELSE IF(iost.GT.0) THEN
         WRITE(lout,*) ' *** ERROR: Reading file particle.dat'
         GO TO 999
      ENDIF
      CLOSE(lun)
C
C---   Initialize the beam spectrum and the partial processes
C
      CALL BGGEN_INI(id,enelec,enpeak,beamlim(1),zcollim)  
C
C---   Initialize the NTUPLE
C
      IF(iwrnt.NE.0) THEN
         lrec=2048
         IQUEST(10)=128000
         CALL HROPEN(lunnt,'bgkin','bggen.nt','N',lrec,iost)
         IF(iost.NE.0) THEN
            WRITE(6,*)'*** ERROR opening NTUPLE, iost=',iost
            GO TO 999
         ENDIF
         idnt=1
         CALL HBNT(idnt,'Kinematics',' ')
         CALL HBNAME(idnt,'run',IEVPROC,'iproc')
         CALL HBNAME(idnt,'run',NTRA,'np[0,20]')
         CALL HBNAME(idnt,'run',ITYPTR(1) ,'ityp(np)')
         CALL HBNAME(idnt,'run',AMTRA(1) ,'am(np)')
         CALL HBNAME(idnt,'run',PTRA(1,1),'ptra(3,np)')
         CALL HBNAME(idnt,'run',NDECTR(1) ,'ndec(np)')
C
      ENDIF
C
C---   Open output file
C
      OPEN(lundat,FILE='bggen.dat',STATUS='UNKNOWN'
     +       ,FORM='UNFORMATTED')
      CALL OPEN_HDDM_OUTPUT('bggen.hddm')
C
      DO i=1,mxpro
         nevpro(i)=0
      ENDDO
C
      DO iev=1,nevent
C
         CALL BGGEN_EVE(id)
C
         IF(NTRA.GT.2) THEN
            ipro=IEVPROC
            IF(ipro.GT.0.AND.ipro.LE.mxpro) THEN
               nevpro(ipro)=nevpro(ipro)+1
            ENDIF
         ENDIF
C
         WRITE(lundat) iev,IEVPROC,NTRA
     +       ,(ITYPTR(i),AMTRA(i),NDECTR(i),(PTRA(j,i),j=1,3),i=1,NTRA)
         IF(iwrnt.NE.0) THEN
            CALL HFNT(idnt)
         ENDIF
C
         idd=0
         DO i=1,NTRA
            IF(NDECTR(i).EQ.0) THEN
               idd=idd+1
               kftr(1,idd)=ITYPTR(i)
               kftr(2,idd)=NDECTR(i)
               qq=0.
               DO j=1,3
                  pptr(j,idd)=PTRA(j,i)
                  qq=qq+pptr(j,idd)**2
               ENDDO
               pptr(4,idd)=SQRT(qq+AMTRA(i)**2)
            ENDIF
         ENDDO
         CALL WRITE_HDDM_EVENT(iev,IEVPROC,idd,kftr(1,1),pptr(1,1))
C
      ENDDO
C
      WRITE(6,*) ' Events Simulated: ',nevent
      IF(nevent.GT.0) THEN
         WRITE(6,2000) 
 2000    FORMAT('  process ',16X,'   events     fraction')
         WRITE(6,2010) (i,cnampro(i),nevpro(i)
     +                ,nevpro(i)/REAL(nevent)*100.,i=1,mxpro) 
 2010    FORMAT(3X,I4,2X,A16,3X,I7,5X,F5.1,' %')
      ENDIF
C
      CLOSE(lundat)
C
      IF(iwrnt.NE.0) THEN
        i=0
        CALL HCDIR('//bgkin',' ')
        CALL HROUT(idnt,i,' ')
        CALL HREND('bgkin')
        CLOSE(UNIT=lunnt)
      ENDIF
C
C-- Write the histograms
C
      lrec=1024
      CALL HROPEN(lun,'HISOUT','bggen.his','N',lrec,iost)
      CALL HROUT(0,icycle,'    ')
      CALL HREND('HISOUT')
      CLOSE(UNIT=lun)
      CALL CLOSE_HDDM_OUTPUT
C
 999  CONTINUE
C
      END

