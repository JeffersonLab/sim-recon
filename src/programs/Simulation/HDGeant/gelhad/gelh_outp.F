#include "gnbase/Flags.h"
      Subroutine Gelh_outp(Iskip)
******-****************************************************************-*******
*
*..Description:
*    Called after each event. For analysis and output flagging.
*    If ISKIP is set to 1 by this routine, event is NOT output (1=output veto).
*
*    Called from gnbbg/gucall.F.
*
*..Author:       Frederic Kral
*..History:
*    Kral     12/24/95 - GELHAD User Output Routine selects gelhad events only.
*    Kral      1/18/96 - Use jphadr switch in new gelhad.db scheme.
*    Kral      1/19/96 - Print number of interactions per event.
*                      - Option jphadr = 3 only saves if ==1 interaction/event.
*                      - Print skipped events.
*    Kral      1/19/96 - Add interaction count.
*    RT Jones  7/10/02 - Suppress special BaBar code using conditional BABAR
*
******-****************************************************************-*******
      IMPLICIT NONE
*..Include:
#ifdef _GELH_
# ifdef BABAR
#include "gnbase/run_db.inc"
# endif
#include "geant321/gcflag.inc"
C   from gelhad V00-00-03 and onwards
#include "gelhad/gelhadused.inc"
C   from gelhad V00-00-05 and onwards
#include "gnbase/gelhad_db.inc"
#endif 
*..Input:
*..Output:
      Integer*4 Iskip           ! Output veto
*..Local:
#ifdef _GELH_
#ifdef BABAR
      Logical*1 First/.true./
      Save First
#endif 
#endif

*=====-================================================================-=======

*..Clear to not skip
      Iskip = 0

#ifdef _GELH_
# ifdef BABAR
      If (First) then
        First = .false.
        If (jphadr_gelhad .ge. 2) then
          Write (6, *)
          CALL UT_ERR_FORM('W', 1, 'GELH_OUTP',
     $      '(''Modified to only save GELHAD beam crossings'')',
     $      0, 0, 0, 0, 0, 0)
          Write (6, *) '*****************************************'
          Write (6, *) 'GELH_OUTP warns: bad mode for normal runs'
          Write (6, *) 'If normal run, set gelhad.db / jphadr < 2'
          If (jphadr_gelhad .ge. 3) then
          Write (6, *) 'This run only saves when ==1 interaction.'
          Endif                 ! (jphadr_gelhad >= 3)
          Write (6, *) '*****************************************'
          Write (6, *)
        Elseif (jphadr_gelhad .ge. 1) then
          CALL UT_ERR_FORM('I', 2, 'GELH_OUTP',
     $      '(''GELHAD on - all beam crossings saved as usual'')',
     $      0, 0, 0, 0, 0, 0)
        Endif                   ! (jphadr_gelhad >= 2)
      Endif                     ! (First)

*..   Skipping logic
      If (jphadr_gelhad .ge. 3) then
*..     GELHAD on with jphadr =3, skip `events' without ==1 GELHAD interaction
        If (gelhadused .and. ngelhperev .eq. 1) then
          if (ev_debug.gt.1) print *,' GELH_OUTP:  Event ',ievent,IDEVT,
     +      '  saved, has ==1 interaction'
        Else
          Iskip = 1             ! SKIP THIS EVENT if not gelhad
        Endif
      Elseif (jphadr_gelhad .ge. 2) then
*..     GELHAD on with jphadr =2, skip `events' without GELHAD interaction
        If (gelhadused) then
          if (ev_debug.gt.1) print *,' GELH_OUTP:  Event ',ievent,IDEVT,
     +      '  saved, has n interactions = ', ngelhperev
        Else
          Iskip = 1             ! SKIP THIS EVENT if not gelhad
        Endif
      Endif                     ! (jphadr_gelhad)

*..Print skips as well
      if (iskip .eq. 1 .and. ev_debug .gt. 1)
     +  print *,' GELH_OUTP:  Event ',ievent,IDEVT,
     +  '  skipped, n = ', ngelhperev,
     +  '  gelhadused = ', gelhadused

# endif

*..Count output
      nevtot_gelh = nevtot_gelh + 1
      if (gelhadused)   nevhad_gelh = nevhad_gelh + 1
      if (iskip .eq. 0) nevout_gelh = nevout_gelh + 1
      if (gelhadused)   ninthad_gelh = ninthad_gelh + ngelhperev
#endif 

c 999  Return
      Return
      End
