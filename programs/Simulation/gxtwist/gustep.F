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
#include "geant321/gcnum.inc"
#include "hdtrackparams.inc"
*
*     -----------------------------------------------------------------
*
      character*4 cnames(15)
      equivalence (NAMES(1),cnames(1))
      character*4 chfrom
      save chfrom
      data chfrom/'NULL'/
      real rnd(100)

      real vertx(3),tofgx,ubuf(99),xint(3,999)
      integer ntbeamx,nttargx,nubuf

      integer ptype
      real xin(4),xout(4),pin(5),pout(5),dEsum,x0(4),p0(5),p1(5),ppol
      common /nt1/ ptype,xin,xout,pin,pout,dEsum,x0,p0,p1,ppol
      character*80 ntform
      parameter (ntform='kind:i,xin(4):r,xout(4):r,'
     +                //'pin(5):r,pout(5):r,dEsum:r,'
     +                //'x0(4):r,p0(5):r,p1(5):r,ppol:r')
      integer ntid1,ntid2,ntid3
      parameter(ntid1=1,ntid2=2,ntid3=3)

      logical hexist
      external hexist   

      if (.not.hexist(ntid1)) then
          call hbnt(ntid1,'microscope hits',' ')
          call hbname(ntid1,'hits',ptype,ntform)
          call hbnt(ntid2,'fixed array hits',' ')
          call hbname(ntid2,'hits',ptype,ntform)
          call hbnt(ntid3,'endpoint array hits',' ')
          call hbname(ntid3,'hits',ptype,ntform)
      endif
 
      CALL GDEBUG

*#define GENERATE_BUT_DO_NOT_TRACK 1
#if GENERATE_BUT_DO_NOT_TRACK
      istop = 1
      return
#endif

* Implement an importance-sampling cascade scheme
#if TUNL_BACKSTREAMING_CASCADE_FACTOR
      if (INWVOL.eq.1.and.NSTEP.gt.0) then
          if (ISTORY.lt.2.and.
     +        NLEVEL.ge.2.and.cnames(2).eq.'AREA'.and.
     +        .not.(NLEVEL.ge.3.and.cnames(3).eq.'TUNL')) then
              if (chfrom.eq.'TUNL') then
                  isave = ISTORY
                  ISTORY = 2
                  call replicate(TUNL_BACKSTREAMING_CASCADE_FACTOR)
*                 ISTORY = isave
              endif
          elseif (ISTORY.eq.0.and.cnames(NLEVEL).eq.'TUNL') then
              if (chfrom.eq.'EDHS') then
                  ISTORY = 1
                  call replicate(EDHS_BACKSTREAMING_CASCADE_FACTOR)
*                 ISTORY = 0
              endif
          endif
      elseif (INWVOL.eq.2) then
          if ((NLEVEL.ge.2.and.cnames(2).eq.'TUNL').or.
     +        (NLEVEL.ge.3.and.cnames(3).eq.'TUNL')) then
              if ((NLEVEL.ge.3.and.cnames(3).eq.'EDHS').or.
     +            (NLEVEL.ge.4.and.cnames(4).eq.'EDHS')) then
                  chfrom = 'EDHS'
              else
                  chfrom = 'TUNL'
              endif
          else
              chfrom = cnames(NLEVEL)
          endif
      else
          chfrom = cnames(NLEVEL)
      endif
#endif

* Place any secondaries generated during this step onto the stack

      if (nosecondaries.eq.0) then
        do i=1,NGKINE
          itypa  = GKIN(5,i)
          if (itypa.ne.4) call GSKING(i)
        enddo
      endif

* Stop wimpy charged particles that are taking forever to range out

      if ((NSTEP.ge.9999).and.(CHARGE.ne.0)) then
        DESTEP = GEKIN
        ISTOP = 1
      endif

* Stop at exit from the tagger area

c     if (NLEVEL.eq.1) then
c       STOP = 1
c     endif

* If not a sensitive volume then exit here

      if (NTMULT.eq.1.and.NSTEP.eq.0) then
        x0(1) = VECT(1)
        x0(2) = VECT(2)
        x0(3) = VECT(3)
        x0(4) = TOFG
        p0(1) = VECT(4)
        p0(2) = VECT(5)
        p0(3) = VECT(6)
        p0(4) = GETOT
        p0(5) = VECT(7)
        call gfvert(1,vertx,ntbeamx,nttargx,tofgx,ubuf,nubuf)
        ppol = ubuf(1)
        p1(1) = ubuf(2)
        p1(2) = ubuf(3)
        p1(3) = ubuf(4)
        p1(4) = ubuf(5)
        p1(5) = ubuf(6)
      elseif (ISVOL.eq.0) then
        return
      endif

* Inside sensitive medium: accumulate info about track segment

      if (ISTOP.ne.0) then		! particle stops
        continue
      elseif (INWVOL.eq.2) then		! particle exits current volume
        continue
      elseif (INWVOL.eq.1) then		! particle enters new volume
        ptype = ipart
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
      if (CNAMES(NLEVEL).eq.'MSFI') then
        call hfnt(ntid1)
      elseif (CNAMES(NLEVEL)(1:2).eq.'FX') then
        call hfnt(ntid2)
      elseif (CNAMES(NLEVEL)(1:4).eq.'ENDP') then
        call hfnt(ntid3)
      elseif (CNAMES(NLEVEL).eq.'XTAL') then
        x0(1) = VECT(1)
        x0(2) = VECT(2)
        x0(3) = VECT(3)
        x0(4) = TOFG
        p0(1) = VECT(4)
        p0(2) = VECT(5)
        p0(3) = VECT(6)
        p0(4) = GETOT
        p0(5) = VECT(7)
        call gfvert(1,vertx,ntbeamx,nttargx,tofgx,ubuf,nubuf)
        ppol = ubuf(1)
        p1(1) = ubuf(2)/ubuf(6)
        p1(2) = ubuf(3)/ubuf(6)
        p1(3) = ubuf(4)/ubuf(6)
        p1(4) = ubuf(5)
        p1(5) = ubuf(6)
      endif
      END

      subroutine replicate(count)
      integer count
      real psave(3),xsave(3),wsave
      integer i
#undef CERNLIB_GEANT321_GCKINE_INC
#undef CERNLIB_GEANT321_GCTRAK_INC
#include "geant321/gckine.inc"
#include "geant321/gctrak.inc"
      if (count.le.1) return
      do i=1,3
        xsave(i) = VERT(i)
        psave(i) = PVERT(i)
        VERT(i) = VECT(i)
        PVERT(i) = VECT(i+3)*VECT(7)
      enddo
      wsave = UPWGHT
      UPWGHT = count-1
      call GSSTAK(0)
      UPWGHT = wsave
      do i=1,3
        VERT(i) = xsave(i)
        PVERT(i) = psave(i)
      enddo
      end
