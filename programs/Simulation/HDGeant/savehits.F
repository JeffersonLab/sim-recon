      subroutine savehits
*
************************************************************************
*                                                                      *
* savehits: dispatches to hits registry functions for each detector    *
*           subsystem in the simulation.                               *
*                                                                      *
************************************************************************
*
* changes: Wed Jun 20 13:19:56 EDT 2007 B. Zihlmann 
*          add ipart to the function calls hitxxxxxxx() 
*

#include "geant321/gckine.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gctmed.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcphys.inc"

#define MAKE_HITS_WHEN_ILOSS_0    1

      character*4 cnames(15)
      equivalence (NAMES(1),cnames(1))

      real xin(4),xout(4),pin(5),pout(5),dEsum
      save xin,xout,pin,pout,dEsum

      integer nlevel_in,inwvol_in
      nlevel_in = NLEVEL
      inwvol_in = INWVOL

* At this point we hide any inner volumes that should not interrupt the
* formation of a single hit by pretending that we are still in the mother;
* the true geometry information is restored before we exit.

      if (cnames(NLEVEL).eq.'FDCA') then ! do not split hits around wire plane
        if (INWVOL.eq.2.and.INFROM.eq.0) then
          INWVOL = 0 ! continue, do not finalize hit
        elseif (INWVOL.eq.1.and.INFROM.eq.1) then
          INWVOL = 0 ! append to previous hit segment
        endif
      elseif (cnames(NLEVEL).eq.'FDCW') then ! discard energy lost in wires
        goto 9
      endif

* Accumulate hit information and save it when particle exits the volume

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
        goto 9
      else
        dEsum = dEsum + DESTEP
        goto 9
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
#ifdef MAKE_HITS_WHEN_ILOSS_0
        if ((ILOSS.eq.0).and.(CHARGE.ne.0)) then
          dEsum = 1e-3   ! 1 MeV in plastic
        endif
#endif
        if (dEsum.gt.0) then
          call hitStartCntr(xin,xout,pin,pout,dEsum,ITRA,ISTAK,ISTORY,
     >          IPART)
        endif
      elseif ((cnames(NLEVEL).eq.'STRA').or.		! CDC straight straw
     +        (cnames(NLEVEL).eq.'STLA')) then	! CDC stereo straw
#ifdef MAKE_HITS_WHEN_ILOSS_0
        if ((ILOSS.eq.0).and.(CHARGE.ne.0)) then
          dEsum = 1e-5   ! 10 KeV in gas
        endif
#endif
        if (dEsum.gt.0) then
          call hitCentralDC(xin,xout,pin,pout,dEsum,ITRA,ISTAK,ISTORY,
     >          IPART)
        endif
      elseif (cnames(NLEVEL).eq.'FDCA') then		! FDC anode drift cell
#ifdef MAKE_HITS_WHEN_ILOSS_0
        if ((ILOSS.eq.0).and.(CHARGE.ne.0)) then
          dEsum = 1e-5   ! 10 KeV in gas
        endif
#endif
        if (dEsum.gt.0) then
          call hitForwardDC(xin,xout,pin,pout,dEsum,ITRA,ISTAK,ISTORY,
     >          IPART)
        endif
      elseif (cnames(NLEVEL).eq.'CERW') then		! Cerenkov truth
        if (dEsum.gt.0) then
          call hitCerenkov(xin,xout,pin,pout,dEsum,ITRA,ISTAK,ISTORY,
     >          IPART)
        endif
      elseif (cnames(NLEVEL).eq.'CPMT') then		! Cerenkov counter
        if ((dEsum.gt.0).and.(IPART.eq.50)) then
          call hitCerenkov(xin,xout,pin,pout,-dEsum,ITRA,ISTAK,ISTORY,
     >          IPART)
        endif
      elseif ((cnames(NLEVEL).eq.'FTOC').or.		! forward TOF counter
     +        (cnames(NLEVEL).eq.'FTOH')) then
#ifdef MAKE_HITS_WHEN_ILOSS_0
        if ((ILOSS.eq.0).and.(CHARGE.ne.0)) then
          dEsum = 1e-2   ! 10 MeV in plastic
        endif
#endif
        if (dEsum.gt.0) then
          call hitForwardTOF(xin,xout,pin,pout,dEsum,ITRA,ISTAK,ISTORY,
     >          IPART)
        endif
      elseif ((cnames(NLEVEL)(1:3).eq.'BCM').or.	! BCal segment
     +        (cnames(NLEVEL)(1:3).eq.'BM1')) then
#ifdef MAKE_HITS_WHEN_ILOSS_0
        if ((ILOSS.eq.0).and.(CHARGE.ne.0)) then
          dEsum = 1e-2  ! 10 MeV in the calorimeter module
        endif
#endif
        if ((dEsum.gt.0).or.(ISTORY.ne.1)) then
          call hitBarrelEMcal(xin,xout,pin,pout,dEsum,ITRA,ISTAK,ISTORY,
     >          IPART)
          ISTORY = 1  ! this particle has entered the BCal (inherited trait)
        endif
      elseif (cnames(NLEVEL).eq.'LGBL') then		! forward calorimeter
#ifdef MAKE_HITS_WHEN_ILOSS_0
        if ((ILOSS.eq.0).and.(CHARGE.ne.0)) then
          dEsum = 1e-2   ! 10 MeV in the calorimeter block
        endif
#endif
        if ((dEsum.ne.0).or.(ISTORY.ne.2)) then
          call hitForwardEMcal
     +                       (xin,xout,pin,pout,dEsum,ITRA,ISTAK,ISTORY,
     >          IPART)
          ISTORY = 2  ! this particle has entered the FCal (inherited trait)
        endif
      elseif ((cnames(NLEVEL).eq.'UPVP').or.		! UPV channel
     +        (cnames(NLEVEL).eq.'UPVC')) then
#ifdef MAKE_HITS_WHEN_ILOSS_0
        if ((ILOSS.eq.0).and.(CHARGE.ne.0)) then
          dEsum = 1e-2   ! 10 MeV in the calorimeter paddle
        endif
#endif
        if ((dEsum.ne.0).or.(ISTORY.ne.3)) then
          call hitUpstreamEMveto
     +                       (xin,xout,pin,pout,dEsum,ITRA,ISTAK,ISTORY,
     >          IPART)
          ISTORY = 3  ! this particle has entered the UPV (inherited trait)
        endif
      endif

* Restore the true geometry flags and return

    9 NLEVEL = nlevel_in
      INWVOL = inwvol_in
      END
