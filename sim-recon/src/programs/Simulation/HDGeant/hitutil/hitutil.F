      subroutine transformCoord(xin,cin,xout,cout)
      real xin(3), xout(3)
      character*(*) cin, cout
#include <geant321/gcvolu.inc>
      character*4 cnames(15)
      equivalence (cnames(1),NAMES(1))
      integer level,saveLevel
      integer levelIn,levelOut
      real xglobal(3)
c
      if (cin.eq.'global') then
        levelIn = 1
      elseif (cin.eq.'local') then
        levelIn = NLEVEL
      else
        do level=1,NLEVEL-1
          if (cin.eq.cnames(level)) goto 10
        enddo
   10   levelIn = level
      endif 
      if (cout.eq.'global') then
        levelOut = 1
      elseif (cout.eq.'local') then
        levelOut = NLEVEL
      else
        do level=1,NLEVEL-1
          if (cout.eq.cnames(level)) goto 20
        enddo
   20   levelOut = level
      endif 
      if (levelIn.eq.levelOut) then
        xout(1) = xin(1)
        xout(2) = xin(2)
        xout(3) = xin(3)
      elseif (levelIn.eq.1) then
        saveLevel = NLEVEL
        NLEVEL = levelOut
        call gmtod(xin,xout,1)
        NLEVEL = saveLevel
      elseif (levelOut.eq.1) then
        saveLevel = NLEVEL
        NLEVEL = levelIn
        call gdtom(xin,xout,1)
        NLEVEL = saveLevel
      else
        saveLevel = NLEVEL
        NLEVEL = levelIn
        call gdtom(xin,xglobal,1)
        NLEVEL = levelOut
        call gmtod(xglobal,xout,1)
        NLEVEL = saveLevel
      endif
      end
