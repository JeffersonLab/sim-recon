*
* $Id$
*
* $Log$
* Revision 1.3  2001/07/17 22:38:40  jonesrt
* Adding hits registry in gustep -rtj
*
* Revision 1.2  2001/07/15 07:31:37  jonesrt
* HDGeant now supportskinematic input from Monte Carlo generators
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
#include "geant321/gcking.inc"
#include "geant321/gcomis.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gctrak.inc"
*
*     -----------------------------------------------------------------
*

#define TOP_CERENKOV_EFFICIENCY 0.25

      character*4 cnames(15)
      equivalence (NAMES(1),cnames(1))
      real rnd(100)

      real xin(4),xout(4),pin(4),pout(4),dEsum
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

* If not a sensitive volume then exit here

      if (ISVOL.eq.0) return

* Inside sensitive medium: accumulate info about track segment

      if (INWVOL.eq.1) then
        xin(1) = VECT(1)
        xin(2) = VECT(2)
        xin(3) = VECT(3)
        xin(4) = TOFG
        pin(1) = VECT(4)*VECT(7)
        pin(2) = VECT(5)*VECT(7)
        pin(3) = VECT(6)*VECT(7)
        pin(4) = GETOT
        dEsum = 0
        return
      elseif (INWVOL.eq.2) then
        xout(1) = VECT(1)
        xout(2) = VECT(2)
        xout(3) = VECT(3)
        xout(4) = TOFG
        pout(1) = VECT(4)*VECT(7)
        pout(2) = VECT(5)*VECT(7)
        pout(3) = VECT(6)*VECT(7)
        pout(4) = GETOT
        dEsum = dEsum + DESTEP
      else
        dEsum = dEsum + DESTEP
        return
      endif

* At end of track segment in sensitive medium: register hit

      if (cnames(NLEVEL).eq.'BCAM') then	! barrel calorimeter
      elseif (cnames(NLEVEL).eq.'CSTI') then	! CDC inner cathode strips
      elseif (cnames(NLEVEL).eq.'CSTO') then	! CDC outer cathode strips
      elseif (cnames(NLEVEL).eq.'STRA') then	! CDC drift straw
      elseif (cnames(NLEVEL).eq.'STRA') then	! CDC drift straw
      elseif (cnames(NLEVEL).eq.'CGAS') then	! Cerenkov counter
      elseif (cnames(NLEVEL).eq.'CGAS') then	! Cerenkov counter
      elseif (cnames(NLEVEL).eq.'FDCA') then	! FDC anode drift cell
      elseif (cnames(NLEVEL).eq.'FDCC') then	! FDC cathode strip
      elseif (cnames(NLEVEL).eq.'LGBL') then	! forward calorimeter
      elseif (cnames(NLEVEL).eq.'FTOF') then	! forward TOF counter
      elseif (cnames(NLEVEL).eq.'VRTX') then	! start counter
      endif

      END
