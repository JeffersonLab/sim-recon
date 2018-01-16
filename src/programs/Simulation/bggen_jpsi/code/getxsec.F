      REAL FUNCTION GETXSEC(E)
C
C---    Return the cross section (in barns) from the input table
C---    Interpolated 
C
      IMPLICIT NONE
      REAL E,DIVDIF
      INTEGER ipol
C
      INCLUDE 'bg_reac.inc'
C
      GETXSEC=0.
C      write(6,*) E,ELREAC
      IF(E.LT.ELREAC(1).OR.E.GT.ELREAC(2)) GO TO 999
C
      ipol=3
      GETXSEC=DIVDIF(XSREAC(1),ESREAC(1),NPXREAC,E,ipol)*1.E-9 ! nb -->b
C      i=(E-ELREAC(1))/((ELREAC(2)-ELREAC(1))/(NPXREAC-1))+0.1
C      i=i+1
C      write(6,*) E,i,ESREAC(i),XSREAC(i)
C
 999  RETURN
      END
