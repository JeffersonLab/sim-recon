C
      DOUBLE PRECISION FUNCTION PYR(IX)
      IMPLICIT NONE
      INTEGER IX
      REAL a
C
      IX=0
      CALL RANLUX(a,1)
      PYR=DBLE(a)
      RETURN
      END

