      subroutine gamate(alo,aup,yes)
c
c
c check if material contains material with
c A between alo and aup
c      subroutine gamat(alo,aup,yes)
      implicit none
      logical yes
      integer *4 i,nmix
      real *4 alo,aup
#include "geant321/gcmate.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gcbank.inc"
c
      yes=.true.
      if(q(jma+6).lt.aup+0.5.and.q(jma).gt.alo-0.5) return
      nmix=abs(q(jma+11))
      if(nmix.ge.2) then
         do 1000 i=1,nmix
            if(q(jmixt+i).lt.aup+0.5.and.q(jmixt+i).gt.alo-0.5) return
 1000    continue
      endif
      yes=.false.
      return
      end
