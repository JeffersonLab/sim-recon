*
* $Id$
*
* $Log$
* Revision 1.5  2002/07/10 19:53:08  jonesrt
* - moved open/close of hbook file from hdgeant.f to uginit.F/uglast.F so that
*   it gets called from hdgeant++ [rtj]
* - fixed Makefile to make hdgeant++ export its symbols to dynamic COMIS
*   functions [rtj]
*
* Revision 1.4  2002/07/10 14:57:18  jonesrt
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
* Revision 1.3  2002/06/28 19:01:03  jonesrt
* Major revision 1.1 -Richard Jones, Chris Gauthier, University of Connecticut
*
* 1. Added hadronic interactions for photons with the Gelhad package
*    http://www.slac.stanford.edu/BFROOT/www/Computing/Offline/Simulation/gelhad.html
*    Routines affected are:
*     - uginit.F : added new card GELH to set up gelhad parameters and
*                  call to gelh_vrfy() to print out their values.
*     - uglast.F : added call to gelh_last() to print out summary info.
*     - gtgama.F : Gelhad replacement for standard Geant routine that adds
*                  simulation of hadronic photoproduction processes.
*     - gelhad/  : contains a number of new functions (Fortran) and includes
*                  to support the hadronic photoproduction simulation.
*
* 2. Added muon-pair production by stealing every (Melectron/Mmuon)**2 pair
*    production events and trying to convert to muon pairs.  The deficit in
*    e+/e- events resulting from this theft is negligible.  The angular
*    distribution of muon pairs is generated using the general Geant method
*    in gpairg.F with the electron mass replaced by the muon mass.
*    Routines affected are:
*     - gpairg.F : added a switch to replace e+/e- with mu+/mu- in a small
*                  fraction of the pair-production vertices.
*
* Revision 1.2  2001/07/24 05:37:19  jonesrt
* First working prototype of hits package -rtj
*
* Revision 1.1  2001/07/08 06:24:35  jonesrt
* First release of the Geant3 geometry package for Hall D based on hdds.  -rtj
*
* Revision 1.1.1.1  1995/10/24 10:21:52  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.35  by  S.Giani
*-- Author :
      SUBROUTINE UGLAST
*
************************************************************************
*                                                                      *
*            Termination routine to print histograms and statistics    *
*                                                                      *
************************************************************************
#include "geant321/gcomis.inc"
*
*     -----------------------------------------------------------------
*
      call gelh_last()
      CALL GLAST
*
*             Close HIGZ
*
      CALL IGEND
*
      call HROUT(0,icycle,' ')
      call HREND('RZfile')
      call closeOutput()
      END
