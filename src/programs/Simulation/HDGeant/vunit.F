*
* $Id$
*
* $Log$
* Revision 1.2  2001/11/12 23:03:02  brash
* Found a small bug in my original modification to this routine, which
* exists only here for the alpha BTW.  The bug was that I was checking for small
* values of rotation matrix elements, to avoid fpe's, and I should have been
* checking for small ABSOLUTE values, since of course the matrix elements can
* be negative.
*
* Revision 1.1  2001/10/04 06:19:36  brash
* This routine is part of the cern library, but suffered arithmetic traps on the alpha.  I have modified to trap these situations. EJB
*
* Revision 1.1.1.1  1996/02/15 17:50:18  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE VUNIT (A,B,N)
C
C CERN PROGLIB# F121    VUNIT           .VERSION KERNFOR  1.0   710701
C ORIG. 01/07/71
C
      DIMENSION A(*),B(*)
C
      IF (N.LE.0)  RETURN
      XX= 0.
      DO 9 I=1,N
C EJB - added the following line to avoid arithmetic exceptions on the
C alpha.  Seems to work.
      if(abs(a(i)).lt.1.0E-10)a(i)=0.0
C
    9 XX= XX + A(I)*A(I)
      IF (XX.LE.0.)  RETURN
      YY=1./ SQRT (XX)
      DO 19 I= 1,N
   19 B(I)=A(I) * YY
      RETURN
      END
