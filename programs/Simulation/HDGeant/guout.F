*
* $Id$
*
* $Log$
* Revision 1.3  2001/10/29 17:39:23  jonesrt
* - added mc truth info to output event for internal track/photon generators
* - added special code for background studies, selected by the conditional
*   #define BACKGROUND_STUDIES   (in gustep.F)
* - added conditional code to disable normal event output for bg studies, using
*   #define DISABLE_OUTPUT    (in guout.F)
* Both of the above defines are disabled in the distribution code by default.
* -rtj-
*
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
* #define DISABLE_OUTPUT 1

#ifndef DISABLE_OUTPUT
      call loadOutput();
      call flushOutput();
#endif

      END
