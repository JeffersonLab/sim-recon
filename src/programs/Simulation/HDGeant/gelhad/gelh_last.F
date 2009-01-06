      Subroutine Gelh_last ! statistics printed
*
*..History:
*    Kral      1/18/96 - Print event counts.
*    Kral      1/19/96 - Add interaction count.
*
      implicit none
#include "gnbase/gelhad_db.inc"
#include "gelhad/gelhadused.inc"
      logical print/.true./
      save print
      logical init/.false./
      save init

*..Print summary
      if (print) then
        if (.not. init) then
          init = .true.
          if (jphadr_gelhad .ne. 0) then
            write (6, *)
            write (6, *) 
     $        'GELH_LAST:  GELHAD event summary'
            write (6, *) '  nevtot  = ', nevtot_gelh
            write (6, *) '  nevhad  = ', nevhad_gelh,
     $        REAL(nevhad_gelh)/REAL(MAX(nevtot_gelh, 1))
            write (6, *) '  nevout  = ', nevout_gelh
            write (6, *) '  ninthad = ', ninthad_gelh
            write (6, *)
          endif
        endif
      endif

      return
      end
