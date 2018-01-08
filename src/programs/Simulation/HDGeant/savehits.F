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
*           Oct 11 2012, yqiang, add hits in RICH volume: RAWN RDCD

#include "geant321/gckine.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gctmed.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcphys.inc"

#define MAKE_HITS_WHEN_ILOSS_0    1

      character*4 cnames(15)
      equivalence (NAMES(1),cnames(1))

      integer nlevellast,nameslast,numberlast
      common /gcvolulast/nlevellast,nameslast(15),numberlast(15)
      character*4 cnameslast(15)
      equivalence (nameslast(1),cnameslast(1))
      save /gcvolulast/
      data nlevellast/0/

      logical hitopen
      save hitopen
      data hitopen/.false./

      real xin(4),xout(4),pin(5),pout(5),dEsum
      save xin,xout,pin,pout,dEsum

      integer isame
      integer level,nlevel_in

C  Reason 1 for being here: entry to a sensitive volume
      if (INWVOL.eq.1) then
        if (hitopen) then
          if (NLEVEL.eq.nlevellast) then
            do level=NLEVEL,1,-1
              if (NAMES(level).ne.nameslast(level).or.
     +            NUMBER(level).ne.numberlast(level)) then
                write(6,*) 'Warning in savehits: unsaved hit',
     +                     ' information found from sensitive volume ',
     +                     cnameslast(level),', hit discarded.'
                hitopen = .false.
                goto 3
              endif
            enddo
            return
          endif
    3     continue
        endif

*	Initialize a new hit when a particle enters a sensitive volume,
*	unless it is just coming out from an inner daughter volume,
*	in which case it should just continue to add to the existing hit.

        if (istop.eq.0) then
          do level=1,NLEVEL
            nameslast(level) = NAMES(level)
            numberlast(level) = NUMBER(level)
          enddo
          nlevellast = NLEVEL
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
          hitopen = .true.
        endif
        return

C  Reason 2 for being here: stepping inside a sensitive volume
      elseif (INWVOL.eq.0) then
        dEsum = dEsum + DESTEP
        if (ISTOP.eq.0) then
          return
        endif
        xout(1) = VECT(1)
        xout(2) = VECT(2)
        xout(3) = VECT(3)
        xout(4) = TOFG
        pout(1) = VECT(4)
        pout(2) = VECT(5)
        pout(3) = VECT(6)
        pout(4) = GETOT
        pout(5) = VECT(7)

C  Reason 3 for being here: exiting from *ANY* volume
      elseif (ISVOL.eq.1) then
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
      endif

* 	At this point, we know that the particle is either leaving the
*	current volume or stopping inside it.

C  If no hit is being recorded then nothing to do
      if (.not.hitopen) then
        return
      endif

*  If about to enter a daughter volume then nothing to do
      if (INWVOL.eq.2.and.INGOTO.ne.0) then
        return
      endif

*	Otherwise, we may be exiting a sensitive region or one of its contents;
*	if a content then we have no way of knowing whether we are at the same
*	time leaving the outer sensitive volume or merely passing through one
*	of its interior interfaces.  At this point there is no way to find out
*	except to do the rather expensive call to GINVOL for the current point.
*	To do that we will need to save and restore the geometry state info.

      if (INWVOL.eq.2.and.NLEVEL.gt.nlevellast) then
        call GSCVOL
        NLEVEL = nlevellast
        call GINVOL(VECT,isame)
        call GFCVOL
        if (isame.ne.0) then
          return
        endif
      endif

*	At this point we back out to the level of the sensitive volume
*	to save the hits information, then make sure that the true
*	geometry information is restored before we exit.

      do level=nlevellast,1,-1
        if (level.gt.NLEVEL.or.
     +      NAMES(level).ne.nameslast(level).or.
     +      NUMBER(level).ne.numberlast(level)) then
          write(6,*) 'Warning in savehits: unsaved hit',
     +               ' information found from sensitive volume ',
     +               cnameslast(level),', hit discarded.'
          hitopen = .false.
          return
        endif
      enddo
      nlevel_in = NLEVEL
      NLEVEL = nlevellast

c      write(6,*)cnames(NLEVEL)
*	At end of track segment in sensitive medium: register hit
      if (cnames(NLEVEL)(1:3).eq.'PTS') then ! TPOL silicon detector sectors
#ifdef MAKE_HITS_WHEN_ILOSS_0
        if ((ILOSS.eq.0).and.(CHARGE.ne.0)) then
          dEsum = 0.6e-3   ! 600 keV in silicon
        endif
#endif
        if (dEsum.gt.0) then
          call hitTPOL(xin,xout,pin,pout,dEsum,ITRA,ISTAK,ISTORY,
     >          IPART)
        endif
      elseif (cnames(NLEVEL)(1:3).eq.'PSF') then ! Pair spectrometer fine counters
#ifdef MAKE_HITS_WHEN_ILOSS_0
        if ((ILOSS.eq.0).and.(CHARGE.ne.0)) then
          dEsum = 1e-3   ! 1 MeV in plastic
        endif
#endif
        if (dEsum.gt.0) then
          call hitPS(xin,xout,pin,pout,dEsum,ITRA,ISTAK,ISTORY,
     >          IPART)
        endif
      elseif (cnames(NLEVEL)(1:3).eq.'PSC') then ! PS coarse counters
#ifdef MAKE_HITS_WHEN_ILOSS_0
        if ((ILOSS.eq.0).and.(CHARGE.ne.0)) then
          dEsum = 1e-3   ! 1 MeV in plastic
        endif
#endif
        if (dEsum.gt.0) then
          call hitPSC(xin,xout,pin,pout,dEsum,ITRA,ISTAK,ISTORY,
     >          IPART)
        endif
      elseif (cnames(NLEVEL).eq.'STRC') then  ! start counter paddle
#ifdef MAKE_HITS_WHEN_ILOSS_0
        if ((ILOSS.eq.0).and.(CHARGE.ne.0)) then
          dEsum = 1e-3   ! 1 MeV in plastic
        endif
#endif
        if (dEsum.gt.0) then
          call hitStartCntr(xin,xout,pin,pout,dEsum,ITRA,ISTAK,ISTORY,
     >          IPART)
        endif
      elseif ((cnames(NLEVEL).eq.'STRA').or.  ! CDC straight straw
     +        (cnames(NLEVEL).eq.'STLA').or.  ! CDC stereo straw
     +        (cnames(NLEVEL).eq.'STLB')) then ! CDC close-packed stereo straw
#ifdef MAKE_HITS_WHEN_ILOSS_0
        if ((ILOSS.eq.0).and.(CHARGE.ne.0)) then
          dEsum = 1e-5   ! 10 KeV in gas
        endif
#endif
        if (dEsum.gt.0) then
          call hitCentralDC(xin,xout,pin,pout,dEsum,ITRA,ISTAK,ISTORY,
     >          IPART)
        endif
      elseif (cnames(NLEVEL)(1:3).eq.'FDA') then ! FDC anode drift cell
#ifdef MAKE_HITS_WHEN_ILOSS_0
        if ((ILOSS.eq.0).and.(CHARGE.ne.0)) then
          dEsum = 1e-5   ! 10 KeV in gas
        endif
#endif
        if (dEsum.gt.0) then
          call hitForwardDC(xin,xout,pin,pout,dEsum,ITRA,ISTAK,ISTORY,
     >          IPART)
        endif
      elseif (cnames(NLEVEL).eq.'CERW') then  ! Cerenkov truth
        if (dEsum.gt.0) then
          call hitCerenkov(xin,xout,pin,pout,dEsum,ITRA,ISTAK,ISTORY,
     >          IPART)
        endif
      elseif (cnames(NLEVEL).eq.'CPPC') then  ! Cerenkov counter
        if ((dEsum.gt.0).and.(IPART.eq.50)) then
          call hitCerenkov(xin,xout,pin,pout,-dEsum,ITRA,ISTAK,ISTORY,
     >          IPART)
        endif
      elseif (cnames(NLEVEL).eq.'RAWN') then ! DIRC window: truth
        if (dEsum.gt.0) then
          call hitDirc(xin,xout,pin,pout,dEsum,ITRA,ISTAK,ISTORY,
     >          IPART)
        endif
      elseif (cnames(NLEVEL).eq.'RDCD') then ! DIRC Cathode: counter
        if ((dEsum.gt.0).and.(IPART.eq.50)) then
          call hitDirc(xin,xout,pin,pout,-dEsum,ITRA,ISTAK,ISTORY,
     >          IPART)
        endif
      elseif ((cnames(NLEVEL).eq.'FTOC').or.  ! forward TOF counter
     +        (cnames(NLEVEL).eq.'FTOX').or.  ! forward TOF counter
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
      elseif (cnames(NLEVEL)(1:2).eq.'BM') then  ! BCal segment
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
      elseif (cnames(NLEVEL).eq.'GCAL') then  ! GCal segment
#ifdef MAKE_HITS_WHEN_ILOSS_0
        if ((ILOSS.eq.0).and.(CHARGE.ne.0)) then
          dEsum = 1e-2  ! 10 MeV in the calorimeter module
        endif
#endif
        if ((dEsum.gt.0).or.(ISTORY.ne.1)) then
          call hitGapEMcal(xin,xout,pin,pout,dEsum,ITRA,ISTAK,ISTORY,
     >          IPART)
          ISTORY = 1  ! this particle has entered the GCal (inherited trait)
        endif
      elseif (cnames(NLEVEL).eq.'LGBL'.or.cnames(NLEVEL).eq.'LTBK') then  ! forward calorimeter
#ifdef MAKE_HITS_WHEN_ILOSS_0
        if ((ILOSS.eq.0).and.(CHARGE.ne.0)) then
          dEsum = 1e-2   ! 10 MeV in the calorimeter block
        endif
#endif
        if ((dEsum.ne.0).or.(ISTORY.ne.2)) then
          call hitForwardEMcal(xin,xout,pin,pout,dEsum,ITRA,ISTAK,
     >          ISTORY,IPART,0)
          ISTORY = 2  ! this particle has entered the FCal (inherited trait)
        endif
      elseif (cnames(NLEVEL).eq.'LGLG') then  ! forward calorimeter light guides
        if (dEsum.gt.0) then
          call hitForwardEMcal(xin,xout,pin,pout,dEsum,ITRA,ISTAK,
     >          ISTORY,IPART,1)
        endif
      elseif (cnames(NLEVEL).eq.'LTBL') then  ! Compton calorimeter
#ifdef MAKE_HITS_WHEN_ILOSS_0
        if ((ILOSS.eq.0).and.(CHARGE.ne.0)) then
          dEsum = 1e-2   ! 10 MeV in the calorimeter block
        endif
#endif
        if ((dEsum.ne.0).or.(ISTORY.ne.2)) then
          call hitComptonEMcal
     +                       (xin,xout,pin,pout,dEsum,ITRA,ISTAK,ISTORY,
     >          IPART)
          ISTORY = 2  ! this particle has entered the CCal (inherited trait)
        endif
      elseif ((cnames(NLEVEL).eq.'UPVP').or.  ! UPV channel
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

* Mark this hit as closed, so that it does not get registered twice

      hitopen = .false.

* If this is the end of tracking for this particle, reset the saved state

      if (ISTOP.ne.0) then
        nlevellast = 0
      endif

* Restore the true geometry state and return

      NLEVEL = nlevel_in
      END
