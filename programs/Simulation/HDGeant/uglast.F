*
* $Id$
*
* $Log$
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
*     - gelhad/* : contains a number of new functions (Fortran) and includes
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
      call closeOutput()
      END
 
