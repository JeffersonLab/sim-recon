      Subroutine Gelh_vrfy ! verify parameters from database/ffread cards
*
*..History:
*    Kral      1/18/96 - Use this when the db common block is used (default).
*                      - Clear event counts.
*    Kral      1/19/96 - Add interaction count.
*
      implicit none
#include "gnbase/gelhad_db.inc"
#include "gelhad/gelhadused.inc"
      logical print/.true./
      save print
      logical init/.false./
      save init

*..Clear output
      nevtot_gelh = 0
      nevhad_gelh = 0
      nevout_gelh = 0
      ninthad_gelh = 0

*..Print verification
      if (print) then
        if (.not. init) then
          init = .true.
          if (jphadr_gelhad .ne. 0) then
            write (6, *)
            write (6, *) 
     $        'GELH_VRFY:  GELHAD hadronic interactions are on'
            write (6, *) '  jphadr  = ', jphadr_gelhad
            write (6, *) '  ecut    = ', ecut_gelhad
            write (6, *) '  scale   = ', scale_gelhad
            write (6, *) '  mode    = ', mode_gelhad
            write (6, *) '  ethresh = ', ethresh_gelhad
            write (6, *)
          endif
        endif
      endif

      return
      end
