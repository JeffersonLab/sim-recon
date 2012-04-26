*
* $Id$
*
* $Log$
* Revision 1.14  2004/12/08 14:43:24  davidl
* Change argument 4 of second call to GSVERT from 0 to 0.0 to avoid compiler warnings
*
* Revision 1.13  2003/12/10 15:32:57  jonesrt
* -control.in : never mind [rtj]
* -gukine.F : fixed a bug in the setting of polarization ppol used
*        for polarization studies of the coherent bremsstrahlung source [rtj]
* -gustep.F : changed background studies facility to split data across
*        separate ntuples, one for each virtual detector [rtj]
* -hit*.F : modified behaviour from "quit" to "truncate" in the case where
*        the number of hits exceeds the maximum allowed for that counter [rtj]
*
* Revision 1.12  2003/07/28 15:42:58  jonesrt
* - gukine.F - added photon polarization as an additional attached info to
*              vertex [rtj]
*
* Revision 1.11  2003/07/28 15:31:33  jonesrt
* - gukine.F - added conditional BEAM_BOX_SIZE to enable simulations with
*   artificial electron beam motion superimposed on coherent bremsstrahlung [rtj]
*
* Revision 1.10  2003/01/02 23:49:33  jonesrt
* - included updates in gustep.F with conditional code for background
*   studies, by R. Jones and C. Gauthier [rtj]
* - moved the beam origin a meter upstream in gukine.F to make room for
*   additional shielding in the collimator cave [rtj]
*
* Revision 1.9  2002/07/10 14:57:18  jonesrt
* - fixed wierd problem with g77 compiler that wanted to interpret "slash star"
*   in a fortran comment line as a comment indicator a-la-c (complained about
*   unterminated comment) so I just removed the asterisk - rtj.
* - corrected the statistics printout from gelh_last() -rtj.
* - changed confusing use of VSCAN (card SCAP) to define the origin for single
*   particle generation; now gukine.F uses PKINE (card KINE) for both origin
*   and direction of single-particle generator, with the following format:
*      KINE kind energy theta phi vertex(1) vertex(2) vertex(3)
* - fixed gelh_outp() to remove the BaBar-dependent code so that it correctly
*   updates the photo-hadronic statistics that get reported at gelh_last() -rtj.
* - updated gelhad/Makefile to follow the above changes -rtj.
*
* Revision 1.8  2001/12/18 20:32:12  jonesrt
* I added the track="int" information to the output event, at the request of
* Dave Doughty.  Track numbers are assigned by Geant in the order of declaration,
* which is just the order they appear in the Reaction section, so it is not too
* difficult to figure out which track goes with which final-state product.
* However there is presently no internal identifier in the Reaction section that
* matches up to the track number. Even calling it a track is a bit of a stretch
* because it is assigned to neutrals as well as charged particles.  But that is
* the Geant nomenclature and it is simple to decode.
* I also added some comments to the control cards file control.in that might
* make it easier for a newbie to run his own simulations.
* -rtj-
*
* Revision 1.7  2001/10/30 11:52:36  jonesrt
* - fixed bug in gukine.F in coherent beam simulation
*   where variable spot was in meters but treated as if it were cm    -rtj-
*
* Revision 1.6  2001/10/29 17:39:23  jonesrt
* - added mc truth info to output event for internal track/photon generators
* - added special code for background studies, selected by the conditional
*   #define BACKGROUND_STUDIES   (in gustep.F)
* - added conditional code to disable normal event output for bg studies, using
*   #define DISABLE_OUTPUT    (in guout.F)
* Both of the above defines are disabled in the distribution code by default.
* -rtj-
*
* Revision 1.5  2001/08/02 03:08:05  jonesrt
* Now the BEAM data card is supported, with correct generation of
* coherent bremsstrahlung radiation. -rtj
*
* Revision 1.4  2001/07/27 21:04:09  jonesrt
* With this release, HDGeant version 1.0 is now in beta.  -rtj
*
* Revision 1.3  2001/07/24 05:37:16  jonesrt
* First working prototype of hits package -rtj
*
* Revision 1.2  2001/07/15 07:31:37  jonesrt
* HDGeant now supportskinematic input from Monte Carlo generators
* via the routines in hddmInput.c  -rtj
*
* Revision 1.1  2001/07/10 18:05:46  jonesrt
* imported several of the gu*.F user subroutines for Hall D customization -rtj
*
* Revision 1.1.1.1  1995/10/24 10:21:52  cernlib
* Geant
*
*

#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.35  by  S.Giani
*-- Author :
      SUBROUTINE GUKINE
*
************************************************************************
*                                                                      *
*             Generates Kinematics for primary tracks                  *
*                                                                      *
************************************************************************
*
#include "geant321/gcunit.inc"
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcscan.inc"
#include "geant321/gcomis.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcnum.inc"
#include "hdtrackparams.inc"
#include "backgrounds.inc"
#include "cobrems.inc"
*
      DIMENSION VERTEX(4),PLAB(5)
      DIMENSION RNDM(20)

      real tgen
      real unif01(100)
      integer i,j
      character*20 pname
      integer nubuf
      real ubuf(99)
      real pmin, pmax, thetamin, thetamax

*
*     -----------------------------------------------------------------
*
      UPWGHT = 1
      ISTORY = 0

      ev = IDEVT
      do i=1,10
        ev = ev/10.
        if (ev.lt.10) goto 2
      enddo
    2 if (int(ev).eq.ev) then
        write(LOUT,*) IDEVT," events simulated"
      endif

*     Get the current values of the random number seeds. Do this
*     here so we can overwrite them in the first "if" block below
*     for the case when we find seeds in the generated events file.
*     Not all generated events files will contain seeds so we want
*     to record the seed values as they are (i.e. specified in control.in)
*     If one of the built-in generators is used, the seed values
*     will also come from control.in. In all cases, the seeds that are actually
*     used are stored in the file via the storeseeds call at the end
*     of this subroutine.
*     
      call GRNDMQ(iseed1,iseed2,0,'G')

*
*
*              Try input from MonteCarlo generator first
*
      itry = nextInput()
      if (itry .eq. 0) then
        itry = loadInput()
*
*       Check for random number seeds in the input file. If they are
*       there, then the values of iseed1 and iseed2 will be overwritten
*       by this call. If they are not there, then they will be untouched
*       by the call.
*
        call getseeds(iseed1, iseed2)
        call GRNDMQ(iseed1,iseed2,0,'S')

*
*       Fake a tagger hit of the correct energy by adding up the energy
*       of all generated tracks and assigning them to the trigger time.
*
        PLAB(1) = 0
        PLAB(2) = 0
        PLAB(3) = 0
        PLAB(4) = -PMASS
        do nt=1,NTRACK
          call GFKINE(nt,VERT,PVERT,IPART,IVERT,ubuf,nubuf)
          call GFPART(IPART,pname,ITRTYP,AMASS,CHARGE,TLIFE,ubuf,nubuf)
          PLAB(1) = PLAB(1) + PVERT(1)
          PLAB(2) = PLAB(2) + PVERT(2)
          PLAB(3) = PLAB(3) + PVERT(3)
          PLAB(4) = PLAB(4)
     +            + sqrt(AMASS**2+PVERT(1)**2+PVERT(2)**2+PVERT(3)**2)
        enddo
        PLAB(5) = PLAB(4)
        VERTEX(1) = VERT(1)
        VERTEX(2) = VERT(2)
        VERTEX(3) = VERT(3)
        VERTEX(4) = TOFG
        call hitTagger(VERTEX,VERTEX,PLAB,PLAB,0.,nt,0,0)

        if (bgrate.gt.0) then
*
*        Superimpose background in the form of coherent bremsstrahlung
*        beam photons sent down the photon beam line.  They are generated
*       with a random time distribution over the duration of the gate
*       to simulate the actual conditions of a bremsstrahlung beam.
*
          ngen=0
          tgen=bggate(1)
          do i=1,99999
            call grndm(unif01,100)
            do j=1,100
              tgen=tgen-log(unif01(j))/bgrate
              if (tgen.gt.bggate(2)) goto 10
              call beamgen(tgen)
              ngen=ngen+1
            enddo
          enddo
   10     continue
c         print *, ngen,' background photons generated this event'
        endif
      elseif (itry .ne. 9) then
         ieorun = 1
         ieotri = 1
         return
*
*              Try coherent bremsstrahlung beam generation next
*
      elseif (E.gt.0) then
         call beamgen(0.)
         call storeInput(IDRUN,IDEVT,1);
*
*              If all else fails, do automatic single-track generation
*
      else
         TOFG=0
         VERTEX(1)=VSCAN(1)
         VERTEX(2)=VSCAN(2)
         VERTEX(3)=VSCAN(3)+(15.-30.*RNDM(1))
         IF (IKINE.GT.100) THEN
            IK=IKINE-100
            CALL GRNDM(RNDM,3)
*           
*           If the PLOG(TLOG) card is non-zero in control.in, then
*           distribute evenly in the log of total momentum(theta).
*           Otherwise, distribute evenly in total momentum(theta).
*           3/17/2009  DL
                            IF (plog_particle_gun.EQ.0) THEN
               PABS=PKINE(1)+PKINE(4)*(RNDM(1)-0.5)
                            ELSE
               pmin=PKINE(1)-0.5*PKINE(4)
               pmax=PKINE(1)+0.5*PKINE(4)
               IF (pmin.LE.0) pmin=0.100
               PABS=pmin*(pmax/pmin)**RNDM(1)
                            ENDIF
                            IF (tlog_particle_gun.EQ.0) THEN
               THETA=(PKINE(2)+PKINE(5)*(RNDM(2)-0.5))*DEGRAD
                            ELSE
               thetamin=PKINE(2)-0.5*PKINE(5)
               thetamax=PKINE(2)+0.5*PKINE(5)
               IF (thetamin.LE.0) thetamin=0.9
               THETA=(thetamin*(thetamax/thetamin)**RNDM(2))*DEGRAD
                            ENDIF
            PHI=(PKINE(3)+PKINE(6)*(RNDM(3)-0.5))*DEGRAD
         ELSE
            IK=IKINE
            CALL GRNDM(RNDM,2)
            PABS=PKINE(1)
            THETA=PI*RNDM(1)
            PHI=TWOPI*RNDM(2)
         ENDIF
         PLAB(1) = PABS*SIN(THETA)*COS(PHI)
         PLAB(2) = PABS*SIN(THETA)*SIN(PHI)
         PLAB(3) = PABS*COS(THETA)

         CALL GSVERT(VERTEX,0,0,0.0,0,NVERT)
         CALL GSKINE(PLAB,IK,NVERT,0,0,NT)

         call storeInput(IDRUN,IDEVT,NT);

      endif
*
*              Kinematic debug (controled by ISWIT(1))
*
      IF(IDEBUG.EQ.1.AND.ISWIT(1).EQ.1) THEN
        CALL GPRINT('VERT',0)
        CALL GPRINT('KINE',0)
      ENDIF
*
* If storing particle trajectories, clear the buffers
      if (storetraj.ne.0) then
        call cleartrajectories()
      endif

*     Store the random number seeds used for this event in the 
*     output file. The values for the seeds are determined
*     above since the coherent bremstrahlung generator uses the
*     random number generator.
*
      call storeseeds(iseed1, iseed2)

      END
