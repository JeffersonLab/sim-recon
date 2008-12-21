      SUBROUTINE BG_NTUP_INI(IERR)
C
C---      Initialize the ntuple
C
      IMPLICIT NONE
      INTEGER IERR
C
      COMMON/QUEST/ IQUEST(100)
      INTEGER       IQUEST
C
      INCLUDE 'bg_ctrl.inc'
      INCLUDE 'bg_evec.inc'
C
      INTEGER lrec,idnt,iost
      CHARACTER cpar*3
C
      IERR=1
      lrec=2048
      IQUEST(10)=128000
      CALL HROPEN(LUNWR(3),'bgkin','bggen.nt','N',lrec,iost)
      IF(iost.NE.0) THEN
         WRITE(6,*)'*** ERROR opening NTUPLE, iost=',iost
         GO TO 999
      ENDIF
      idnt=1
      CALL HBNT(idnt,'Kinematics',' ')
      CALL HBNAME(idnt,'run',IEVPROC,'iproc')
      CALL HBNAME(idnt,'run',ITPIN(1,1) ,'itypin(2,2)')
      CALL HBNAME(idnt,'run',AMIN(1) ,'amin(2)')
      CALL HBNAME(idnt,'run',PIN(1,1),'pin(3,2)')
      WRITE(cpar,FMT='(I3)') MXTRA
      CALL HBNAME(idnt,'run',NTRA,'np[0,'//cpar//']')
      CALL HBNAME(idnt,'run',ITPTRA(1,1) ,'ityp(6,np)')
      CALL HBNAME(idnt,'run',AMTRA(1) ,'am(np)')
      CALL HBNAME(idnt,'run',PTRA(1,1),'ptra(3,np)')
C
      IERR=0
 999  RETURN
      END
