*
* $Id$
*
* $Log$
* Revision 1.7  2001/10/30 11:52:36  jonesrt
* - fixed bug in gukine.F in coherent beam simulation
*   where variable spot was in meters but treated as if it were cm    -rtj-
*
* Revision 1.6  2001/10/29 17:39:23  jonesrt
* - added mc truth info to output event for internal track/photon generators
* - added special code for background studies, selected by the conditional
*   #define BACKGROUND_STUDIES   (in gustep.F)
* - added conditional code to disable normal event output for bg studies, using
*   #define DISABLE_OUTPUT    (in guout.F)
* Both of the above defines are disabled in the distribution code by default.
* -rtj-
*
* Revision 1.5  2001/08/02 03:08:05  jonesrt
* Now the BEAM data card is supported, with correct generation of
* coherent bremsstrahlung radiation. -rtj
*
* Revision 1.4  2001/07/27 21:04:09  jonesrt
* With this release, HDGeant version 1.0 is now in beta.  -rtj
*
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
#include "geant321/gcunit.inc"
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcscan.inc"
#include "geant321/gcomis.inc"
*
      DIMENSION VERTEX(3),PLAB(3)
      DIMENSION RNDM(20)
*
      include 'cobrems.inc'
      integer coSwitch
      real freqMaximum
      common /coherentGen/coSwitch,freqMaximum
      data coSwitch,freqMaximum/0,0/
      save /coherentGen/
      real xMinimum,Theta02,CollimPos
      parameter (xMinimum=1e-2)
      parameter (Theta02=1.8)
      parameter (CollimPos=-1500.0)
*
*     -----------------------------------------------------------------
*
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
         call GRNDM(RNDM,5)
         phim = RNDM(1)*TWOPI
         rhom = mospread*sqrt(-2*log(RNDM(2)))
         thxMosaic = rhom*cos(phim)
         thyMosaic = rhom*sin(phim)
         phib = RNDM(3)*TWOPI
         varBeam = (emit/spot)**2 + sigma2MS(t*RNDM(4))
         rhob = sqrt(-2*varBeam*log(RNDM(5)))
         thxBeam = rhob*cos(phib)
         thyBeam = rhob*sin(phib)
         cos45 = 1/sqrt(2d0)
         rotate(1,1) = 0
         rotate(1,2) = cos45			!point (1,0,0) along beam
         rotate(1,3) = -cos45			!point (0,1,1) vertically
         rotate(2,1) = 0
         rotate(2,2) = cos45
         rotate(2,3) = cos45
         rotate(3,1) = 1
         rotate(3,2) = 0
         rotate(3,3) = 0
         call rotmat(rotate,thxBeam-thx-thxMosaic,0d0,0d0)
         call rotmat(rotate,0d0,thyBeam-thy-thyMosaic,0d0)
         do i=1,10000
            call GRNDM(RNDM,4)
            phi = RNDM(1)*TWOPI
            x = xMinimum**RNDM(2)
            coSwitch = coSwitch+1
            if (coSwitch/2*2.eq.coSwitch) then  !try coherent generation
               freq = dNcdxdp(x,phi)
               f = freq*RNDM(3)
               do ip=1,q2points
                  if (f.le.q2weight(ip)) then
                     theta2 = q2theta2(ip)
                     goto 5
                  endif
               enddo
    5          continue
               freq = freq*TWOPI*Theta02
            else                                  !try incoherent generation
               theta2 = Theta02
     +                  *RNDM(3)/(1-RNDM(3)+1e-30)
               freq = dNidxdt2(x,theta2)
               freq = freq*(Theta02+theta2)**2
            endif
            freq = x*freq
            if (freq.gt.freqMaximum) then
               freqMaximum = freq*1.2
            elseif (freq.ge.freqMaximum*RNDM(4)) then
               goto 50
            endif
         enddo
         write(5,*) 'Gukine: photon beam generator failed, giving up'

   50    continue
         theta = sqrt(theta2)*(me/E)
         thetaX = thxBeam+theta*cos(phi)
         thetaY = thyBeam+theta*sin(phi)
         PLAB(3) = E*x
         PLAB(1) = PLAB(3)*thetaX
         PLAB(2) = PLAB(3)*thetaY
         call GRNDM(RNDM,2)
         phic = RNDM(1)*TWOPI
         rhoc = spot*sqrt(-2*log(RNDM(2)))
         VERTEX(1) = (rhoc*cos(phic) + D*thetaX)*100
         VERTEX(2) = (rhoc*sin(phic) + D*thetaY)*100
         VERTEX(3) = CollimPos
         call GSVERT(VERTEX,0,0,0,0,NVERT)
         call GSKINE(PLAB,1,NVERT,0,0,NT)
         call storeInput(IDRUN,IDEVT,NT);
*
*              If all else fails, do automatic single-track generation
*
      else

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
