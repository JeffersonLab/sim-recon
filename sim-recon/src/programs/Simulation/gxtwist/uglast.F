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
      SUBROUTINE UGLAST
*
************************************************************************
*                                                                      *
*            Termination routine to print histograms and statistics    *
*                                                                      *
************************************************************************
#include "geant321/gcomis.inc"
*
*     -----------------------------------------------------------------
*
      call gelh_last()
      CALL GLAST
*
*             Close HIGZ
*
      CALL IGEND
*
      call HROUT(0,icycle,' ')
      call HREND('RZfile')
      call closeOutput()
      END
