      subroutine savehits
*
************************************************************************
*                                                                      *
* savehits: dispatches to hits registry functions for each detector    *
*           subsystem in the simulation.                               *
*                                                                      *
************************************************************************
*
#include "geant321/gckine.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gctmed.inc"
#include "geant321/gctrak.inc"

      character*4 cnames(15)
      equivalence (NAMES(1),cnames(1))

      real xin(4),xout(4),pin(5),pout(5),dEsum
      save xin,xout,pin,pout,dEsum

      if (ISTOP.ne.0) then		! particle stops
        continue
      elseif (INWVOL.eq.2) then		! particle exits current volume
        continue
      elseif (INWVOL.eq.1) then		! particle enters new volume
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
      else
        dEsum = dEsum + DESTEP
        return
      endif

* At end of track segment in sensitive medium: register hit

      dEsum = dEsum + DESTEP
      xout(1) = VECT(1)
      xout(2) = VECT(2)
      xout(3) = VECT(3)
      xout(4) = TOFG
      pout(1) = VECT(4)
      pout(2) = VECT(5)
      pout(3) = VECT(6)
      pout(4) = GETOT
      pout(5) = VECT(7)
      if ((cnames(NLEVEL).eq.'STRC').or.		! vertex cylinder
     +    (cnames(NLEVEL).eq.'STRP')) then		! vertex endcap
        if (dEsum.gt.0) then
          call hitStartCntr(xin,xout,pin,pout,dEsum,ITRA,ISTAK,ISTORY)
        endif
      elseif ((cnames(NLEVEL).eq.'STRA').or.		! CDC straight straw
     +        (cnames(NLEVEL).eq.'STLA')) then	! CDC stereo straw
        if (dEsum.gt.0) then
          call hitCentralDC(xin,xout,pin,pout,dEsum,ITRA,ISTAK,ISTORY)
        endif
      elseif (cnames(NLEVEL).eq.'FDCA') then		! FDC anode drift cell
        if (dEsum.gt.0) then
          call hitForwardDC(xin,xout,pin,pout,dEsum,ITRA,ISTAK,ISTORY)
        endif
      elseif (cnames(NLEVEL).eq.'CSEC') then		! Cerenkov counter
        if (dEsum.gt.0) then
          call hitCerenkov(xin,xout,pin,pout,dEsum,ITRA,ISTAK,ISTORY)
        endif
      elseif ((cnames(NLEVEL).eq.'FTOC').or.		! forward TOF counter
     +        (cnames(NLEVEL).eq.'FTOS')) then
        if (dEsum.gt.0) then
          call hitForwardTOF(xin,xout,pin,pout,dEsum,ITRA,ISTAK,ISTORY)
        endif
      elseif ((cnames(NLEVEL)(1:3).eq.'BCM').or.	! BCal segment
     +        (cnames(NLEVEL)(1:3).eq.'BM1')) then
        if ((dEsum.gt.0).or.(ISTORY.ne.1)) then
          call hitBarrelEMcal(xin,xout,pin,pout,dEsum,ITRA,ISTAK,ISTORY)
          ISTORY = 1  ! this particle has entered the BCal (inherited trait)
        endif
      elseif (cnames(NLEVEL).eq.'LGBL') then		! forward calorimeter
        if ((dEsum.ne.0).or.(ISTORY.ne.2)) then
          call hitForwardEMcal
     +                       (xin,xout,pin,pout,dEsum,ITRA,ISTAK,ISTORY)
          ISTORY = 2  ! this particle has entered the FCal (inherited trait)
        endif
      elseif ((cnames(NLEVEL).eq.'UPVP').or.		! UPV channel
     +        (cnames(NLEVEL).eq.'UPVC')) then
        if ((dEsum.ne.0).or.(ISTORY.ne.3)) then
          call hitUpstreamEMveto
     +                       (xin,xout,pin,pout,dEsum,ITRA,ISTAK,ISTORY)
          ISTORY = 3  ! this particle has entered the UPV (inherited trait)
        endif
      endif
      END
