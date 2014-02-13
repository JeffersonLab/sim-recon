C
C---  Example: read the file produced with BGG
C
      PROGRAM READ_BGG
C
      IMPLICIT NONE
C
C
      REAL      beammom
C
      INTEGER lun,lout,i,j,nevent,lenc,iost,iev,ntra
C
      INTEGER mxtra
      PARAMETER (mxtra=10)
      REAL     am(mxtra),pp(3,mxtra)
      INTEGER  iproc,ityp(mxtra),idec(mxtra)
C
C     ------------------------------------------------------------------
C
      lun=9
      lout=6
C
C
C---   Open file
C
      OPEN(lun,FILE='bggen.dat',STATUS='OLD',IOSTAT=iost
     +       ,FORM='UNFORMATTED')
      IF(iost.NE.0) THEN
         WRITE(lout,*) ' *** ERROR: Missing file bggen.dat'
         GO TO 999
      ENDIF
C
      nevent=0
C
 10   READ(lun,IOSTAT=iost) iev,iproc,ntra
     +    ,(ityp(i),am(i),idec(i),(pp(j,i),j=1,3),i=1,MIN(mxtra,ntra)) 
      IF(iost.EQ.0) THEN
         nevent=nevent+1
         IF(nevent.LT.50) THEN
            WRITE(lout,2000) nevent,iev,iproc
 2000       FORMAT(' event ',2I10,4X,' process:',I2/
     +            ,'   #      GEANT#   mass       decay           '
     +            ,'Momentum')
            DO i=1,ntra
               WRITE(lout,2010) i,ityp(i),am(i),idec(i),(pp(j,i),j=1,3)
 2010          FORMAT(I4,5X,I4,6X,F6.3,5X,I3,5X,3F10.4)
            ENDDO
         ENDIF
C
         GO TO 10
      ENDIF
C
      WRITE(lout,*) ' Read events =',nevent
C
      CLOSE(lun)
C--
C
 999  CONTINUE
C
      END

