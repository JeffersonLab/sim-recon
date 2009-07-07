c Mon Oct  3 16:35:06 PDT 1994
c routine to compute photo-production xsection. A. Snyder
c Tue Nov 28 10:10:58 PST 1995
c modified to use improved xsection routine |sigmag.F|
c Wed Jan 24 16:34:07 PST 1996
c modified to handle 0 xsections. A. Snyder
c
      subroutine gpsig(sigma)
      implicit none
c
c /geant321/ commons
#include "geant321/gcbank.inc"
#include "geant321/gcmate.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gctrak.inc"
C
c /gelsig/ variables
      real *4 sigma
      real *4 temp        !working space
      integer *4 nmix     !number of components in current material
      real *4 amix        !A of mixture component
      real *4 zmix        !Z of mixture component
       real *4 wmix        !weight of mixture component
      integer *4 i        !dummy index
c
      nmix=abs(q(jma+11))
      if(nmix.eq.1) then !not a mixture ?
c
       call sigmag(a,z,getot,sigma)
       if(sigma.eq.0.0) return
       sigma=sigma*1000.0 !mb->mub
       sigma=sigma*6.022e-07*dens/a !convert to cm**-1
c
      else               !a mixture
c
       sigma=0.0
       do 1000 i=1,nmix !loop over components of mixture
        amix=q(jmixt+i)
        zmix=q(jmixt+1*nmix+i)
        call sigmag(amix,zmix,getot,temp)
        temp=temp*1000.0 !mb->mub
        wmix=q(jmixt+2*nmix+i)
        sigma=sigma+temp*wmix
1000   continue
       if(sigma.eq.0.0) return
       sigma=sigma*6.022e-07*dens/a !convert to cm**-1
c
      endif !nmix.eq.1
c
      return
      end
