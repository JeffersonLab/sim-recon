*
* $Id$
*
* $Log$
* Revision 1.3  2001/07/24 05:37:16  jonesrt
* First working prototype of hits package -rtj
*
* Revision 1.2  2001/07/15 07:31:37  jonesrt
* HDGeant now supportskinematic input from Monte Carlo generators
* via the routines in hddmInput.c  -rtj
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
      SUBROUTINE GUKINE
*
************************************************************************
*                                                                      *
*             Generates Kinematics for primary tracks                  *
*                                                                      *
************************************************************************
*
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcscan.inc"
#include "geant321/gcomis.inc"
*
      DIMENSION VERTEX(3),PLAB(3)
      DIMENSION RNDM(2)
*
*     -----------------------------------------------------------------
*
*              Try input from MonteCarlo generator first
*
      itry = nextInput()
      if (itry .eq. 0) then
         itry = loadInput()
      elseif (itry .eq. 9) then
 
         VERTEX(1)=VSCAN(1)
         VERTEX(2)=VSCAN(2)
         VERTEX(3)=VSCAN(3)
         IF(IKINE.GT.100)THEN
            IK=IKINE-100
            THETA=PKINE(2)*DEGRAD
            PHI=PKINE(3)*DEGRAD
         ELSE
            IK=IKINE
            CALL GRNDM(RNDM,2)
            THETA=PI*RNDM(1)
            PHI=TWOPI*RNDM(2)
         ENDIF

         PLAB(1) = PKINE(1)*SIN(THETA)*COS(PHI)
         PLAB(2) = PKINE(1)*SIN(THETA)*SIN(PHI)
         PLAB(3) = PKINE(1)*COS(THETA)

         CALL GSVERT(VERTEX,0,0,0,0,NVERT)
         CALL GSKINE(PLAB,IK,NVERT,0,0,NT)

      else
         ieorun = 1
         ieotri = 1
         return
      endif
*
*              Kinematic debug (controled by ISWIT(1))
*
      IF(IDEBUG.EQ.1.AND.ISWIT(1).EQ.1) THEN
        CALL GPRINT('VERT',0)
        CALL GPRINT('KINE',0)
      ENDIF
*
  99  END
