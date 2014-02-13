      subroutine gmmate(mass)
c
c /geant321/ commons
c
      implicit none
#include "geant321/gcmate.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gcbank.inc"
c
      real *4 mass,anumber
      anumber=abs(q(jma+6))
      mass=anumber*0.939
      return
      end
