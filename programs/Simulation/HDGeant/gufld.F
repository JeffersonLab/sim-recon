*
* $Id$
*
* $Log$
* Revision 1.2  2001/07/24 05:37:16  jonesrt
* First working prototype of hits package -rtj
*
* Revision 1.1  2001/07/10 18:05:45  jonesrt
* imported several of the gu*.F user subroutines for Hall D customization -rtj
*
* Revision 1.1.1.1  1995/10/24 10:21:46  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.25  by  S.Giani
*-- Author :
      SUBROUTINE GUFLD(X,F)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       User routine to compute the magnetic field F             *
C.    *       at space point X                                         *
C.    *                                                                *
C.    *    ==>Called by : GHELIX,GRKUTA                                *
C.    *                                                                *
C.    ******************************************************************
C.
#include <geant321/gctmed.inc>
C.
      DIMENSION X(3),F(3)
C.
C.    ------------------------------------------------------------------
C.
      F(1)=0.
      F(2)=0.
      if (X(3).lt.0) then
        F(3)=FIELDM * (100/(100-X(3)))**3
      elseif (X(3).lt.400) then
        F(3)=FIELDM
      else
        F(3)=FIELDM * (100/(X(3)-300))**3
      endif
C
      END
