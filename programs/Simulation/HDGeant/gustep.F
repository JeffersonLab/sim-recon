*
* $Id$
*
* $Log$
* Revision 1.8  2001/12/19 02:34:55  jonesrt
* Fixed the names of sensitive volumes in the save-hits section of gustep.F,
* also added support for the r="float" parameter of <cdcPoint> in hitCDC.c.
* -rtj-
*
* Revision 1.7  2001/10/29 17:39:23  jonesrt
* - added mc truth info to output event for internal track/photon generators
* - added special code for background studies, selected by the conditional
*   #define BACKGROUND_STUDIES   (in gustep.F)
* - added conditional code to disable normal event output for bg studies, using
*   #define DISABLE_OUTPUT    (in guout.F)
* Both of the above defines are disabled in the distribution code by default.
* -rtj-
*
* Revision 1.6  2001/07/27 21:04:09  jonesrt
* With this release, HDGeant version 1.0 is now in beta.  -rtj
*
* Revision 1.5  2001/07/24 05:37:16  jonesrt
* First working prototype of hits package -rtj
*
* Revision 1.4  2001/07/19 23:25:49  jonesrt
* numerous new files as I develop the prototype hits libraries -rtj
*
* Revision 1.3  2001/07/17 22:38:40  jonesrt
* Adding hits registry in gustep -rtj
*
* Revision 1.2  2001/07/15 07:31:37  jonesrt
* HDGeant now supports kinematic input from Monte Carlo generators
* via the routines in hddmInput.c  -rtj
*
* Revision 1.1  2001/07/10 18:05:47  jonesrt
* imported several of the gu*.F user subroutines for Hall D customization -rtj
*
* Revision 1.1.1.1  1995/10/24 10:21:52  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.35  by  S.Giani
*-- Author :
      SUBROUTINE GUSTEP
*
************************************************************************
*                                                                      *
*             User routine called at the end of each tracking step     *
*             MEC   is the mechanism origin of the step                *
*             INWVOL is different from 0 when the track has reached    *
*                    a volume boundary                                 *
*             ISTOP is different from 0 if the track has stopped       *
*                                                                      *
************************************************************************
*
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gcomis.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gctmed.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcflag.inc"
*
*     -----------------------------------------------------------------
*

* #define BACKGROUND_STUDIES 1

#ifdef BACKGROUND_STUDIES
      integer type
      real xv(4),Etot
      common /bgNtuple/type,xv,Etot
      character*80 bgntdef
      parameter (bgntdef='type:I,xv(4):R,Etot:R')
      integer bgnt
      parameter (bgnt=10)
      logical HEXIST
      external HEXIST
#endif

#define TOP_CERENKOV_EFFICIENCY 0.25

      character*4 cnames(15)
      equivalence (NAMES(1),cnames(1))
      real rnd(100)

      real xin(4),xout(4),pin(5),pout(5),dEsum
      save xin,xout,pin,pout,dEsum

      CALL GDEBUG

* Place any secondaries generated during this step onto the stack

      do i=1,NGKINE
        itypa  = GKIN(5,i)
        if (itypa.ne.4) call GSKING(i)
      enddo

* For explicit Cerenkov generation, apply an inefficiency factor

      if (NGPHOT.gt.100) then
        call GRANOR(rnd(1),rnd(2))
        sigma = sqrt(NGPHOT * (TOP_CERENKOV_EFFICIENCY
     +                      * (1 - TOP_CERENKOV_EFFICIENCY)))
        NGPHOT = TOP_CERENKOV_EFFICIENCY*NGPHOT + rnd(1)*sigma + 0.5
        call GSKPHO(0)
      elseif (NGPHOT.gt.0) then
        call GRNDM(rnd,NGPHOT)
        do i=1,NGPHOT
          if (rnd(i).le.TOP_CERENKOV_EFFICIENCY) then
            call GSKPHO(i)
          endif
        enddo
      endif

* Stop wimpy charged particles that are taking forever to range out

      if ((NSTEP.ge.9999).and.(CHARGE.ne.0)) then
        DESTEP = GEKIN
        ISTOP = 1
      endif

#ifdef BACKGROUND_STUDIES
      if (.not.HEXIST(bgnt)) then
        call HBNT(bgnt,'background particles','')
        call HBNAME(bgnt,'tracks',type,bgntdef)
      endif
      if (INWVOL.eq.1) then
        z = VECT(3)
        r = sqrt(VECT(1)**2 + VECT(2)**2)
        if (cnames(NLEVEL).eq.'COL1') then
          ISTOP = 1
          return
        elseif ((cnames(NLEVEL).eq.'UWIT' .and. abs(z-0.02).lt.0.1)
     +    .or.  (cnames(NLEVEL).eq.'VRTX' .and. abs(r-4.95).lt.0.1)
     +    .or.  (cnames(NLEVEL).eq.'CDCI' .and. abs(r-15.0).lt.0.1)
     +    .or.  (cnames(NLEVEL).eq.'DC12' .and. abs(r-37.0).lt.0.1)
     +    .or.  (cnames(NLEVEL).eq.'CDCO' .and. abs(r-59.0).lt.0.1)
     +    .or.  (cnames(NLEVEL).eq.'BCAL' .and. abs(r-65.0).lt.0.1)
     +    .or.  (cnames(NLEVEL).eq.'FDC ' .and. abs(z-224.).lt.1.0)
     +    .or.  (cnames(NLEVEL).eq.'CERE' .and. abs(z-410.).lt.1.0)
     +    .or.  (cnames(NLEVEL).eq.'FCAL' .and. abs(z-575.).lt.1.0))
     +  then
          xv(1) = VECT(1)
          xv(2) = VECT(2)
          xv(3) = VECT(3)
          xv(4) = r
          type = IPART
          Etot = GETOT
          call HFNT(bgnt)
        endif
      endif
#endif

* If not a sensitive volume then exit here

      if (ISVOL.eq.0) return

* Inside sensitive medium: accumulate info about track segment

      if (INWVOL.eq.1) then
        xin(1) = VECT(1)
        xin(2) = VECT(2)
        xin(3) = VECT(3)
        xin(4) = TOFG
        pin(1) = VECT(4)
        pin(2) = VECT(5)
        pin(3) = VECT(6)
        pin(4) = GETOT
        pin(5) = VECT(7)
        dEsum = 0
        return
      elseif ((ISTOP.ne.0).or.(INWVOL.eq.2)) then
        xout(1) = VECT(1)
        xout(2) = VECT(2)
        xout(3) = VECT(3)
        xout(4) = TOFG
        pout(1) = VECT(4)
        pout(2) = VECT(5)
        pout(3) = VECT(6)
        pout(4) = GETOT
        pout(5) = VECT(7)
        dEsum = dEsum + DESTEP
      else
        dEsum = dEsum + DESTEP
        return
      endif

* At end of track segment in sensitive medium: register hit

      if (CHARGE.ne.0) then
        if (cnames(NLEVEL).eq.'VRTX') then		! start counter
          call hitStartCntr(xin,xout,pin,pout,dEsum,ITRA)
        elseif ((cnames(NLEVEL).eq.'STRA').or.		! CDC straight straw
     +          (cnames(NLEVEL).eq.'STLA')) then	! CDC stereo straw
          call hitCentralDC(xin,xout,pin,pout,dEsum,ITRA)
        elseif ((cnames(NLEVEL).eq.'FDCA').or.		! FDC anode drift cell
     +          (cnames(NLEVEL).eq.'FDCC')) then	! FDC cathode strip
          call hitForwardDC(xin,xout,pin,pout,dEsum,ITRA)
        elseif (cnames(NLEVEL).eq.'CGAS') then		! Cerenkov counter
          call hitCerenkov(xin,xout,pin,pout,dEsum,ITRA)
        elseif (cnames(NLEVEL).eq.'BCAM') then		! barrel calorimeter
          call hitBarrelEMcal(xin,xout,pin,pout,dEsum,ITRA)
        elseif ((cnames(NLEVEL).eq.'FTOC').or.		! forward TOF counter
     +          (cnames(NLEVEL).eq.'FTOS')) then
          call hitForwardTOF(xin,xout,pin,pout,dEsum,ITRA)
        elseif (cnames(NLEVEL).eq.'LGBL') then		! forward calorimeter
          call hitForwardEMcal(xin,xout,pin,pout,dEsum,ITRA)
        endif
      endif

      END
