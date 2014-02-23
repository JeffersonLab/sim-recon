C
C---  Runs PYTHIA
C     Reaction: gamma+p
C---  Input: file "run.input" containing: 
C                      the number of events and the beam energy 
C            file "pythia-geant.dat" contains a table for PYTHIA<->GEANT particle code conversion
C            file "pythia.dat" - redefinition of PYTHIA parameters (from HERMES, adapted to GLUEX)           
C
C
      PROGRAM RUN_PYTH
C
      IMPLICIT NONE
C
      INTEGER PYK,PYCOMP
      DOUBLE PRECISION PYP
      EXTERNAL PYK,PYP,PYCOMP
C
      REAL      beammom
      DOUBLE PRECISION dbeam
      CHARACTER chbeam*6, text*80,cpar*100
C
      INTEGER    mxpgeant,mxpkc
      PARAMETER (mxpgeant=100,mxpkc=1000)
      INTEGER iplund(mxpgeant)   ! PYTHIA particle codes (KF)
      INTEGER kcgean(-mxpkc:mxpkc)  ! GEANT code for the PYTHIA internal code KC (with sign)
C
      INTEGER lun,lout,i,j,nevent,lenc,kf,kc,ks,iost,iev,nlnd,nlusf,ilnd
C
      INTEGER mxlnd
      PARAMETER (mxlnd=200)
      REAL     plund(4,mxlnd)
      INTEGER  klund(6,mxlnd)
C
C     ------------------------------------------------------------------
C
      lun=9
      lout=6
C
C---   Read the input file
C
      OPEN(lun,FILE='run.input',STATUS='OLD',IOSTAT=iost
     +       ,FORM='FORMATTED')
      IF(iost.NE.0) THEN
         WRITE(lout,*) ' *** ERROR: Missing file run.input'
         GO TO 999
      ENDIF
      READ(lun,*,IOSTAT=iost) nevent,beammom 
      IF(iost.NE.0) THEN
         WRITE(lout,*) ' *** ERROR: Reading file run.input'
         GO TO 999
      ENDIF
      CLOSE(lun)
C
      IF(beammom.LT.3.OR.beammom.GT.1.E5) THEN
         WRITE(lout,*) ' *** ERROR: out of range - beammom=',beammom
         GO TO 999
      ENDIF
C
C---   Read the GEANT<->PYTHIA particle table
C      
C
      DO i=1,mxpgeant
         iplund(i)=0
      ENDDO
      DO i=-mxpkc,mxpkc
         kcgean(i)=0
      ENDDO
C
      OPEN(lun,FILE='pythia-geant.map',STATUS='OLD',IOSTAT=iost
     +       ,FORM='FORMATTED')
      IF(iost.NE.0) THEN
         WRITE(lout,*) ' *** ERROR: Missing file pythia-geant.map'
         GO TO 999
      ENDIF
 10   READ(lun,'(A)',IOSTAT=iost) cpar
      IF(iost.EQ.0) THEN
C
         lenc=0
         DO i=1,LEN_TRIM(cpar)
            IF(cpar(i:i).EQ.'!') GO TO 20
            lenc=i
         ENDDO
 20      CONTINUE
         IF(lenc.GE.3) THEN
            READ(cpar(1:lenc),*) j,kf
         ENDIF
         IF(j.GT.0.AND.j.LE.mxpgeant) THEN
            IF(kf.NE.0) THEN
               kc=PYCOMP(kf)
               IF(kc.GT.0.AND.kc.LE.mxpkc) THEN
                  IF(kf.LT.0) kc=-kc
                  kcgean(kc)=j
                  kc=ABS(kc)
                  iplund(j) =kf
C                  write(lout,FMT='(10I8)') j,kf,kc
C
C---             Forbid the decays for particles with GEANT code
C
                  WRITE(cpar,1000) kc,0
 1000             FORMAT('MDCY(',I4,',1)=',I2)
                  CALL PYGIVE(cpar)
               ENDIF
C
            ENDIF
         ENDIF
C
         GO TO 10
C
      ELSE IF(iost.GT.0) THEN
         WRITE(lout,*) ' *** ERROR: Reading file pythia-geant.map'
         GO TO 999
      ENDIF
      CLOSE(lun)
 
C---   Read the pythia settings for JLab energies
C      
      OPEN(lun,FILE='pythia.dat',STATUS='OLD',IOSTAT=iost
     +       ,FORM='FORMATTED')
      IF(iost.NE.0) THEN
         WRITE(lout,*) ' *** ERROR: Missing file pythia.dat'
         GO TO 999
      ENDIF
 30   READ(lun,'(A)',IOSTAT=iost) cpar
      IF(iost.EQ.0) THEN
         CALL PYGIVE(cpar)
         GO TO 30
      ELSE IF(iost.GT.0) THEN
         WRITE(lout,*) ' *** ERROR: Reading file pythia.dat'
         GO TO 999
      ENDIF
      CLOSE(lun)
C
C---   Open output file
C
      CALL OPEN_HDDM_OUTPUT('pytout.hddm')
      OPEN(lun,FILE='pytout.dat',STATUS='UNKNOWN'
     +       ,FORM='UNFORMATTED')
C
C---   Initialize PYTHIA
C
      dbeam=DBLE(beammom)
C
      CALL PYINIT('FIXT','gamma','p+',dbeam)
C
C     Uncomment the following to enable Ko decays in pythia
C      CALL pythia6_setdecay(.true., .true.)
C
      DO iev=1,nevent
C
         CALL PYEVNT
         CALL PYEDIT(15)
         nlnd=PYK(0,1)
C
         DO ilnd=1,MIN(nlnd,mxlnd)
            DO i=1,5
               klund(i,ilnd)=PYK(ilnd,i)
            ENDDO
            klund(6,ilnd)=0
            kf=klund(2,ilnd)
            kc=PYCOMP(kf)
            IF(kf.LT.0) kc=-kc
            IF(ABS(kc).LE.mxpkc) THEN
               klund(6,ilnd)=kcgean(kc)
            ENDIF
            DO i=1,3
               plund(i,ilnd)=REAL(PYP(ilnd,i))
            ENDDO
            plund(4,ilnd)=REAL(PYP(ilnd,5))
C
         ENDDO
C
         IF(nlnd.GT.0) THEN
			   CALL WRITE_HDDM_EVENT(iev, beammom, nlnd, klund, plund)
            WRITE(lun) iev,beammom,nlnd
     +              ,((klund(i,ilnd),i=1,6) 
     +               ,(plund(i,ilnd),i=1,4)
     +               ,ilnd=1,nlnd) 
         ENDIF
C
      ENDDO
C
      CALL CLOSE_HDDM_OUTPUT
      CLOSE(lun)
      CALL PYSTAT(1)
C--
C
 999  CONTINUE
C
      END

