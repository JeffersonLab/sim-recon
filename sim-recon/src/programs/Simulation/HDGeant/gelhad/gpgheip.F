
      subroutine gpgheip(idvec,ppvec,kinergy)
c
c *** compute distance to next hadronic interaction point ***
c *** this routine is an interface to gheisha8 ***
c *** nve 06-apr-1988 cern geneva ***
c
c called by : guphad (user routine)
c origin : f.carminati
c
#include "geant321/gcflag.inc"
#include "geant321/gcbank.inc"
#include "geant321/gckine.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcmate.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcphys.inc"
#include "geant321/gcjloc.inc"
c --- gheisha commons ---
#include "gelhad/ghcdes/prntfl.inc"
      real *4 ppvec,kinergy
      integer *4 idvec
c
c --- initialize relevant gheisha variables at first pass ---
      if (ifinit(4) .eq. 0) call gheini
c
      if (z .lt. 1.0) go to 1000
c
      kk=abs(q(jma+11))
      if (kk .gt. 1) go to 10
c
      sig=ghesig(ppvec,kinergy,a,a,z,1.0,1,dens,0.0,idvec)
      go to 20
c
 10   continue
      qcor=0.0
      if(jtm.gt.0) then
         lnve=lq(jtm)
         if (lnve .gt. 0) qcor=q(lnve+26)
      endif
      sig=ghesig(ppvec,kinergy,a,q(jmixt+1),q(jmixt+kk+1),
     $           q(jmixt+2*kk+1),kk,dens,qcor,idvec)
c
 20   continue
      if (sig .le. 0.0) go to 1000
      shadr=zintha/sig
      if (nprt(9)) print 2000,kk,sig,shadr
 2000 format(' *gpghei* kk,sig,shadr = ',i3,1x,2(g12.5,1x))
      go to 9999
c
c --- ensure no interaction in current medium ---
c
 1000 continue
      shadr=big
      if (nprt(9)) print 2001,kk,sig,shadr
 2001 format(' *gpghei* === no interaction in current medium ==='/
     $ ' kk,sig,shadr = ',i3,1x,2(g12.5,1x))
c
 9999 continue
      end
