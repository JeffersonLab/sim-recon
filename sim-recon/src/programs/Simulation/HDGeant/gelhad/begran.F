* Random Number Generators for BEGET, JETSET, etc...
* should be based on GEANT

* Doug Wright, Feb 1994
* Fred Kral, September 23, 1994 - Beget v021 now handles all but begran.F.
*                               - This routine added to bbforce.F.
*                               - This fragment moved from bbrandom.F.
 
******************************
* BEGET
      real function begran()
****************************** 
      implicit none
#include "geant321/gcunit.inc"
c      integer*4 i
      real*4    r
      logical*4 first
      save first
      data first /.true./
 
      if (first) then
        first = .false.
        write (lout,1000)
      end if
 1000 format (/'   >> this version of begran (BEGET random numbers) ',
     + ' uses GEANT grndm'/)

      call grndm(r,1)
      begran = r

      return
      end
