      real function sigmagamma(e)
      implicit none
      real *4 e
      integer *4 mode/0/ !mode=0=>homemade based on pdb
c
      real *4 elo,ehi,sigmalo,sigmahi,lelo,lehi,le
      integer *4 i
c
      real *4 sigmatable(8),etable(8)
      data sigmatable/120.0,522.0,522.0,209.0,266.0,152.0,143.0,122.0/
      data etable/0.2,0.3,0.35,0.475,0.7,1.5,2.0,10./
c
      logical init/.false./
c
      if(.not.init) then
       write(6,*) 'mode=',mode,'in /sigmagamma/'
       init=.true.
      endif !.not.init
c
      if(mode.eq.0) then !homebrew?
c
      sigmagamma=sigmatable(1)
      if(e.lt.etable(1)) return !off bottom ?
      sigmagamma=sigmatable(8)
      if(e.gt.etable(8)) return !off top
c
c find bin
      do 1000 i=1,7
       if(e.ge.etable(i).and.e.lt.etable(i+1)) then
        elo=etable(i)
        ehi=etable(i+1)
        sigmalo=sigmatable(i)
        sigmahi=sigmatable(i+1)
        go to 1099
       endif !e.ge.etable(i).and.e.lt.etable(i+1)
1000  continue
1099  continue
c
c interpolate
      lelo=alog(elo)
      lehi=alog(ehi)
      le=alog(e)
      sigmagamma=sigmalo+(le-lelo)*(sigmahi-sigmalo)/(lehi-lelo)
c
      endif !mode.eq.0
c
      return
      end
