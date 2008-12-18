c Tue Oct  4 10:35:00 PDT 1994
c subtract /da/ and /dz/ from /a/ and /a/ for material /material/
c A. Snyder
c
      subroutine grmate(da,dz)
c
c /geant321/ commons
c
      implicit none
#include "geant321/gcmate.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gcbank.inc"
C
C
c
c /grmate/ variables
c
      real *4 da                !amount to reduce /a/
      real *4 dz                !amount to reduce /z/
      integer *4 nmix           !number of mixtures
      integer *4 i              !loop index
c
      q(jma+6)=q(jma+6)-da*sign(1.0,q(jma+6))
      q(jma+7)=q(jma+7)-dz*sign(1.0,q(jma+7))
      a=abs(q(jma+6))
      z=abs(q(jma+7))
c
      nmix=abs(q(jma+11))
      if(nmix.le.1) return !pure ?
      do 1000 i=1,nmix
       q(jmixt+i)=q(jmixt+i)-da*sign(1.0,q(jmixt+i))
       q(jmixt+nmix+i)=q(jmixt+nmix+i)-dz*sign(1.0,q(jmixt+nmix+i))
1000  continue
c
      return
      end





