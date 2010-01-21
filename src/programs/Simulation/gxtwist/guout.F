*
* $Id$
*
*
* Revision 1.1.1.1  1995/10/24 10:21:52  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.35  by  S.Giani
*-- Author :
      SUBROUTINE GUOUT
*
C.    ******************************************************************
C.    *                                                                *
C.    *       User routine called at the end of each event.            *
C.    *                                                                *
C.    ******************************************************************
C.
C.
#include "geant321/gcomis.inc"
C.    ------------------------------------------------------------------
C.

      integer iskip
      integer nseen
      call gelh_outp(iskip)
      nseen = loadOutput();

#ifdef DISABLE_OUTPUT
      if (0.gt.1) then
#elif WRITE_ONLY_IF_SOMETHING_WAS_SEEN
      if (iseen.gt.0) then
#else
      if (iseen.ge.0) then
#endif
         call flushOutput();
      endif

      END
