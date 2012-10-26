      SUBROUTINE XSEC_PRI
C
C---     Print a cross section table in nb
C
      IMPLICIT NONE
C
      INTEGER iproc,np,ip,iost,np1
      REAL    e1,e2,ef1,ef2,eb,xs,xsec(3000),epn(3000)
      REAL XSEC_JPSI
C
 10   WRITE(6,1000)
 1000 FORMAT('Enter the input numbers: IPROC E1 E2 NP ',/,
     +  '  where IPROC=1 for J/psi',/,
     +  '        E1,E2 - min,max of the energy range, say 8. 12.',/
     +  '        NP - number of points including E1,E2')
      READ(5,*,IOSTAT=iost) iproc,e1,e2,np
      IF(iost.NE.0) THEN
         WRITE(6,*) ' *** ERROR: Reading input string'
         GO TO 10
      ENDIF
      IF(iproc.LT.1.OR.iproc.GT.1) THEN
         WRITE(6,*) 'Wrong number ',iproc
         GO TO 10
      ENDIF
      IF(e1.GE.e2.OR.e1.LT.0..OR.e2.GT.1000.) THEN
         WRITE(6,*) 'Wrong values e1,e2 ',e1,e2
         GO TO 10
      ENDIF
      IF(np.LT.2.OR.np.GT.3000) THEN
         WRITE(6,*) 'Wrong values np ',np
         GO TO 10
      ENDIF
C
      np1=0
      DO ip=1,np
         eb=e1+(e2-e1)/(np-1)*(ip-1)
         IF(iproc.EQ.1) THEN
            xs=XSEC_JPSI(eb)
         ENDIF
         IF(xs.GT.1.E-20) THEN
            np1=np1+1
            IF(np1.EQ.1) ef1=eb
            ef2=eb
            xsec(np1)=xs
            epn(np1)=eb
         ENDIF
      ENDDO
C
      WRITE(6,FMT='(4X,2F10.4,2X,I4)') ef1,ef2,np1
      WRITE(6,FMT='(4X,10E11.4)') (xsec(ip),ip=1,np1)
C      WRITE(6,FMT='(4X,10E11.4)') (epn(ip),ip=1,np1)
C
 999  CONTINUE
      END
C
