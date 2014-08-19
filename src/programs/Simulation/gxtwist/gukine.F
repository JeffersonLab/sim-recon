*
* $Id$
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
#include "geant321/gcunit.inc"
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcscan.inc"
#include "geant321/gcomis.inc"
#include "cobrems.inc"
#include "halo.inc"
*
      DIMENSION VERTEX(4),PLAB(5)
      DIMENSION RNDM(20)

#if 0
      real tgen
      real unif01(100)
      integer i,j
      character*20 pname
      integer nubuf
      real ubuf(99)
#endif

      integer idhalo
      parameter (idhalo=9876)
      real haloxy
      external haloxy
      logical hexist
      external hexist
*
*     -----------------------------------------------------------------
*
      UPWGHT = 1
      ISTORY = 0

      ev = IDEVT
      do i=1,10
        ev = ev/10.
        if (ev.lt.10) goto 2
      enddo
    2 if (int(ev).eq.ev) then
        write(LOUT,*) IDEVT," events simulated"
      endif
*
*              Try input from MonteCarlo generator first
*
      itry = nextInput()
      if (itry .eq. 0) then
         itry = loadInput()
      elseif (itry .ne. 9) then
         ieorun = 1
         ieotri = 1
         return
*
*              Try coherent bremsstrahlung beam generation next
*
      elseif (E.gt.0) then
         call beamgen(0.)
         call storeInput(IDRUN,IDEVT,1);
*
*              If all else fails, do automatic single-track generation
*
      else
         TOFG=0
         VERTEX(1)=VSCAN(1)
         VERTEX(2)=VSCAN(2)
         VERTEX(3)=VSCAN(3)
         IF (IKINE.GT.100) THEN
            IK=IKINE-100
            CALL GRNDM(RNDM,3)
            PABS=PKINE(1)+PKINE(4)*(RNDM(1)-0.5)
            THETA=(PKINE(2)+PKINE(5)*(RNDM(2)-0.5))*DEGRAD
            PHI=(PKINE(3)+PKINE(6)*(RNDM(3)-0.5))*DEGRAD
         ELSE
            IK=IKINE
            CALL GRNDM(RNDM,2)
            PABS=PKINE(1)
            THETA=PI*RNDM(1)
            PHI=TWOPI*RNDM(2)
         ENDIF
         PLAB(1) = PABS*SIN(THETA)*COS(PHI)
         PLAB(2) = PABS*SIN(THETA)*SIN(PHI)
         PLAB(3) = PABS*COS(THETA)
*
*   If the incident track is on the z axis then simulate the actual
*   electron beam profile, including a central gaussian and a halo
*   modeled according to CASA technical note JLAB-TN-06-048.
*
         if (vertex(1).eq.0 .and. vertex(2).eq.0 .and.
     +       plab(1).eq.0 .and. plab(2).eq.0) then
            if (.not.hexist(idhalo)) then
               call hbook2(idhalo-1,'halo work histogram',
     +                     150,-1.5,1.5,150,-1.5,1.5,0.)
               call hbook2(idhalo,'beam y vs x',
     +                     150,-1.5,1.5,150,-1.5,1.5,0.)
               call hbook2(idhalo+1,'beam px vs x',
     +                     150,-1.5,1.5,150,-3.0,3.0,0.)
               call hbook2(idhalo+2,'beam py vs y',
     +                     150,-1.5,1.5,150,-3.0,3.0,0.)
               call hbook2(idhalo+3,'beam py versus px',
     +                     150,-3.0,3.0,150,-3.0,3.0,0.)
               do ix=1,150
                  x = (ix-75.5)/50.
                  do iy=1,150
                     y = (iy-75.5)/50.
                     call hfill(idhalo-1,x,y,haloxy(x*1e-2,y*1e-2,1))
                  enddo
               enddo
            endif
            call grndm(rndm,1)
            if (rndm(1).lt.fhalo/0.52) then
               call hrndm2(idhalo-1,vertex(1),vertex(2))
               call grndm(rndm,2)
               phig = rndm(1)*TWOPI
               rhog = sqrt(-2*log(rndm(2)))
               plab(1) = plab(3)*(0.2*rhog*cos(phig)-vertex(1))/8000.
               plab(2) = plab(3)*(0.1*rhog*sin(phig)-vertex(2))/4000.
            else
               call grndm(rndm,4)
               phig = rndm(1)*TWOPI
               rhog = sqrt(-2*log(rndm(2)))
               thetaX = (emitx/spot)*rhog*cos(phig)
               thetaY = (emity/spot)*rhog*sin(phig)
               phig = rndm(3)*TWOPI
               rhog = sqrt(-2*log(rndm(4)))
               vertex(1) = (spot*rhog*cos(phig)-thetaX*D)*1e2
               vertex(2) = (spot*rhog*sin(phig)-thetaY*D)*1e2
               plab(1) = plab(3)*thetaX
               plab(2) = plab(3)*thetaY
            endif
            call hfill(idhalo,vertex(1),vertex(2),1.)
            call hfill(idhalo+1,vertex(1),plab(1)*1e3,1.)
            call hfill(idhalo+2,vertex(2),plab(2)*1e3,1.)
            call hfill(idhalo+3,plab(1)*1e3,plab(2)*1e3,1.)
         endif

         CALL GSVERT(VERTEX,0,0,0.0,0,NVERT)
         CALL GSKINE(PLAB,IK,NVERT,0,0,NT)

         call storeInput(IDRUN,IDEVT,NT);

      endif
*
*              Kinematic debug (controled by ISWIT(1))
*
      IF(IDEBUG.EQ.1.AND.ISWIT(1).EQ.1) THEN
        CALL GPRINT('VERT',0)
        CALL GPRINT('KINE',0)
      ENDIF
*
      END
