*
* $Id$
*
* $Log$
* Revision 1.2  2001/07/24 05:37:16  jonesrt
* First working prototype of hits package -rtj
*
* Revision 1.1  2001/07/10 18:05:46  jonesrt
* imported several of the gu*.F user subroutines for Hall D customization -rtj
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
*
      call loadOutput();
      call flushOutput();
      END
