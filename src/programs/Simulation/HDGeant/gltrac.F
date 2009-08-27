*
* $Id: gltrac.F,v 1.1 2006/04/15 04:38:38 jonesrt Exp $
*
* $Log: gltrac.F,v $
* Revision 1.1  2006/04/15 04:38:38  jonesrt
* gltrac.F, gsstak.F
*    - replacements for geant321 library functions that enable stacking of
*      secondaries with repeat counts and saving of ISTORY on the stack,
*      both for the purposes of enabling cascaded simulations. [rtj]
* gustep.F
*    - new code to support electron beam dump simulations with a two-level
*      cascade to enhance the statistics of dump-related backgrounds. [rtj]
* taggerCoords2.xls
*    - updates to the dimensions of the building and electron beam dump. [rtj]
* hdds/Spectrometer.xml, hdds/TaggerArea.xml
*    - new geometry description including the electron beam dump attached
*      by a corridor to the tagger building. [rtj]
* hdds/Makefile hdds/ElectronDump.xml
*    - added new document to describe the electron beam dump geometry [rtj]
* gxtwist, gxtwist++, hdgeant [deleted]
*    - binary files removed from repository [rtj]
*
* Revision 1.1.1.1  1995/10/24 10:21:41  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/04 13/12/94  15.36.22  by  S.Giani
*-- Author :
      SUBROUTINE GLTRAC
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    SUBR. GLTRAC                                                *
C.    *                                                                *
C.    *   Extracts next track from stack JSTAK and prepares commons    *
C.    *    /GCTRAK/, /GCKINE/ and /GCVOLU/                             *
C.    *                                                                *
C.    *   Called by : GTREVE                                           *
C.    *   Authors   : R.Brun, F.Bruyant                                *
C.    *                                                                *
C.    ******************************************************************
C.
#define USE_UPWGHT_AS_REPEAT_COUNT 1
*
#include "geant321/gcbank.inc"
#include "geant321/gckine.inc"
#include "geant321/gcnum.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcphys.inc"
#include "geant321/gcstak.inc"
#include "geant321/gctmed.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcvolu.inc"
      DIMENSION RNDM(5)
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION P2,GETOTD,GEKIND
      DOUBLE PRECISION PXD,PYD,PZD,ONE,HNORM,DAMASS,PP
#endif
      PARAMETER (ONE=1)
C.
C.    ------------------------------------------------------------------
*** Restore USE_UPWGHT_AS_REPEAT_COUNT November 18 2006
***  I added a reset ISTORY=3 to the top of gukine() and turned the
***  USE_UPWGHT_AS_REPEAT_COUNT back on.  Notice to other users: if you
***  turn this off, you will get hundreds of truth tags per event.
*** richard.t.jones@uconn.edu
***
*** Disabled use of UPWGHT November 6,2006
*** I have commented out the #define of USE_UPWGHT_AS_REPEAT_COUNT
*** above because there appears to be a problem with it turned on.
*** The symptom is that after the first dozen events or so, every
*** single call to hitUpstreamEMveto(...) passes ISTORY=3 thereby
*** disabling the UPVTruthShowers in the output. Disabling the
*** UPWGHT feature results in consistent generation of UPVTruthShowers.
*** I will leave it to Richard to track down the true source of
*** the problem so that use of UPWGHT can be re-enabled.
*** David Lawrence

*** Modification introduced March 26, 2006
***  There is a  "user word" UPWGHT that is associated with each particle
***  on the temporary stack.  In the standard usage this word is a priority
***  that is used to select the order in which particles are tracked, in
***  conjunction with the SORD control card.  In this modification I change
***  the meaning of UPWGHT to represent a repeat count for the stacked
***  particle.  That is, each time a particle is retrieved from the stack
***  its value of UPWGHT on the stack is decremented and it is removed from
***  the stack only when its UPWGHT reaches zero.  This behaviour is useful
***  in implementing an importance sampling scheme.  Note that the default
***  value of UPWGHT is 1 so this modification has no effect unless user
***  code in gustep() or elsewhere overwrites its value.  The SORD card
***  will have no effect if USE_UPWGHT_AS_REPEAT_COUNT is in effect.
*** richard.t.jones@uconn.edu

*
* *** Extract next track from stack JSTAK
*
#ifndef USE_UPWGHT_AS_REPEAT_COUNT
      IF(ISTORD.EQ.1) THEN
*
* *** User ordering of tracks if requested
         CALL GSTORD
      ENDIF
#endif
      ISTAK = IQ(JSTAK+1)
      IQ(JSTAK+1) = ISTAK -1
      JST = JSTAK +NWSTAK*IQ(JSTAK+1) +3
#ifdef USE_UPWGHT_AS_REPEAT_COUNT
      IF (Q(JST+12).GT.1) THEN
         IQ(JSTAK+1) = ISTAK
      ENDIF
#endif
      ITRA   = IQ(JST+1)
      IF (ITRA.LT.0) THEN
         ITRA = -ITRA
      ELSE
*
*        This is a new track. We set to zero the stack number and
*        update the vertex number
         ISTAK = 0
         JK=LQ(JKINE-ITRA)
         IVERT=Q(JK+6)
      ENDIF
      IPART  = IQ(JST+2)
      DO 60 I = 1,3
         VERT(I) = Q(JST+3+I)
        PVERT(I) = Q(JST+6+I)
   60 CONTINUE
      TOFG   = Q(JST+10)
      SAFETY = Q(JST+11)
      UPWGHT = Q(JST+12)
#ifdef USE_UPWGHT_AS_REPEAT_COUNT
*     print *, 'pop stacked track',istak,', copy',int(Q(JST+12)),
*    +         ', generation',IQ(JST+3)
      Q(JST+12) = Q(JST+12)-1
      UPWGHT = 1
#endif
*
* *** Prepare tracking parameters
*
      VECT(1) = VERT(1)
      VECT(2) = VERT(2)
      VECT(3) = VERT(3)
      PXD = PVERT(1)
      PYD = PVERT(2)
      PZD = PVERT(3)
      P2 = PXD**2+PYD**2+PZD**2
      IF(P2.GT.0.) THEN
         PP    = SQRT(P2)
         HNORM = ONE/PP
         VECT(4) = PVERT(1)*HNORM
         VECT(5) = PVERT(2)*HNORM
         VECT(6) = PVERT(3)*HNORM
         VECT(7) = PP
      ELSE
         VECT(4) = 0.
         VECT(5) = 0.
         VECT(6) = 1.
         VECT(7) = 0.
      ENDIF
*
*  ** Reload Particle characteristics, if needed
*
      IF (IPART.NE.IPAOLD) THEN
         JPA = LQ(JPART-IPART)
         DO 90 I = 1,5
            NAPART(I) = IQ(JPA+I)
   90    CONTINUE
         ITRTYP = Q(JPA+6)
         AMASS  = Q(JPA+7)
         CHARGE = Q(JPA+8)
         TLIFE  = Q(JPA+9)
         IUPD   = 0
         IPAOLD = IPART
      ENDIF
*
      DAMASS = AMASS
      GETOTD = SQRT(P2+DAMASS**2)
      GEKIND = GETOTD - DAMASS
      GETOT  = GETOTD
      GEKIN  = GEKIND
*
      IF (ITRTYP.EQ.7) THEN
*
* *** Cerenkov photon. Retrieve polarisation
         JPO = LQ(JSTAK-1)+(ISTAK-1)*3
         POLAR(1) = Q(JPO+1)
         POLAR(2) = Q(JPO+2)
         POLAR(3) = Q(JPO+3)
      ELSE
         CALL GEKBIN
      ENDIF
*
      SLENG  = 0.
      NSTEP  = 0
      NTMSTO = NTMSTO +1
      NTMULT = NTMSTO
#ifdef USE_UPWGHT_AS_REPEAT_COUNT
      ISTORY = IQ(JST+3)
#else
      ISTORY = 0
#endif
*
*  ** Initialize interaction probabilities
*
      IF (ITRTYP.EQ.1) THEN
*      Gammas
         CALL GRNDM(RNDM,5)
         ZINTPA = -LOG(RNDM(1))
         ZINTCO = -LOG(RNDM(2))
         ZINTPH = -LOG(RNDM(3))
         ZINTPF = -LOG(RNDM(4))
         ZINTRA = -LOG(RNDM(5))
      ELSE IF (ITRTYP.EQ.2) THEN
*       Electrons
         CALL GRNDM(RNDM,3)
         ZINTBR = -LOG(RNDM(1))
         ZINTDR = -LOG(RNDM(2))
         ZINTAN = -LOG(RNDM(3))
      ELSE IF (ITRTYP.EQ.3) THEN
*       Neutral hadrons
         CALL GRNDM(RNDM,2)
         SUMLIF = -CLIGHT*TLIFE*LOG(RNDM(1))
         ZINTHA = -LOG(RNDM(2))
      ELSE IF (ITRTYP.EQ.4) THEN
*       Charged hadrons
         CALL GRNDM(RNDM,3)
         SUMLIF = -CLIGHT*TLIFE*LOG(RNDM(1))
         ZINTHA = -LOG(RNDM(2))
         ZINTDR = -LOG(RNDM(3))
      ELSE IF (ITRTYP.EQ.5) THEN
*       Muons
         CALL GRNDM(RNDM,5)
         SUMLIF = -CLIGHT*TLIFE*LOG(RNDM(1))
         ZINTBR = -LOG(RNDM(2))
         ZINTPA = -LOG(RNDM(3))
         ZINTDR = -LOG(RNDM(4))
         ZINTMU = -LOG(RNDM(5))
      ELSE IF (ITRTYP.EQ.7) THEN
*       Cerenkov photons
         CALL GRNDM(RNDM,1)
         ZINTLA = -LOG(RNDM(1))
      ELSE IF (ITRTYP.EQ.8) THEN
*       Ions
         CALL GRNDM(RNDM,2)
         ZINTHA = -LOG(RNDM(1))
         ZINTDR = -LOG(RNDM(2))
      ENDIF
*
*   * Prepare common /GCVOLU/ and structure JGPAR, if needed
*
      IF (NJTMAX.LE.0) THEN
        IF (GONLY(NLEVEL).EQ.0.) NLEVEL=0
        CALL GMEDIA (VECT, NUMED)
      ENDIF
      INFROM = 0
*                                                             END GLTRAC
      END
 
