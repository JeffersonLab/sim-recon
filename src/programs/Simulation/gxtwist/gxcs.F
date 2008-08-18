#define CERNLIB_COMIS true
*
* June 12, 2000 -rtj
* Modified calls to csext so that arg1 (string) is not too long
*
* $Id$
*
* $Log$
* Revision 1.1  2006/01/24 21:49:49  jonesrt
* Initial revision
*
* Revision 1.1  2001/07/08 06:24:33  jonesrt
* First release of the Geant3 geometry package for Hall D based on hdds.  -rtj
*
* Revision 1.3  2001/03/07 00:42:19  radphi
* Changes made by jonesrt@zeus
* several geometry fixes, LGD gain improvement -rtj
*
* Revision 1.2  1998/07/02 03:55:41  radphi
* Changes made by kurylov@jlabs4
* Small improvements to geometry, corrections to materials, hits definitions -AAK
*
* Revision 1.1.1.1  1995/10/24 10:21:49  cernlib
* Geant
*
*
#include "geant321/pilot.h"
#if defined(CERNLIB_COMIS)
*CMZ :          20/06/95  09.32.44  by  S.Ravndal
*-- Author :
      SUBROUTINE GXCS
C.
C.    ******************************************************************
C.    *                                                                *
C.    *        To initialize the COMIS package                         *
C.    *        To declare addresses of FORTRAN routines and COMMONs    *
C.    *        which may be invoked from COMIS routines                *
C.    *        (one can call CSOMAP instead)                           *
C.    *                                                                *
C.    ******************************************************************
#include "geant321/gcbank.inc"
#include "geant321/gcmate.inc"
#include "geant321/gctmed.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gcflag.inc"
#include "geant321/gctrak.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gccuts.inc"
#include "geant321/gclist.inc"
#include "geant321/gcnum.inc"
#include "geant321/gconst.inc"
#include "geant321/gcphys.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gcmulo.inc"
#include "geant321/gcomis.inc"
#include "geant321/gcsets.inc"

      DIMENSION P(1)
*
      EXTERNAL GINIT,GZINIT,GDINIT,GPRINT,GPSETS,GXCLOS
      EXTERNAL GSVERT,GSKINE,GSKING,GOPEN,GFIN,GCLOSE
      EXTERNAL GFOUT
      EXTERNAL GMATE,GSMATE,GSMIXT,GSTMED,GSTPAR,GPART,GPHYSI
      EXTERNAL GFMATE,GPIONS
      EXTERNAL GTRIG,GTRIGI,GTRIGC,GTREVE,GIDROP
      EXTERNAL GSVOLU,GSPOS,GSPOSP,GSDVN,GSDVS,GGCLOS,GOPTIM
      EXTERNAL GSROTM,GSORD,GSDET,GSDETH,GSDETV,GSATT
      EXTERNAL GPLMAT,GSAHIT,GSCHIT,GSDIGI,GSXYZ,GDEBUG
      EXTERNAL GPCXYZ,GDCXYZ,GDXYZ,GDAHIT,GDCHIT,GDHITS,GDHEAD
      EXTERNAL GDOPEN,GDCLOS,GDRAW,GDRAWC,GDSCAL,GDMAN,GDCOL
      EXTERNAL GDELET,GDAXIS,GDRAWT
      EXTERNAL GSCANK,GSCANU,GSCANO
      EXTERNAL UGLAST
*
C.
C.    ------------------------------------------------------------------
C.
      CALL PAWCS
*
      CALL CSCOM('GCLINK,GCBANK,GCCUTS,GCFLAG,GCKINE,GCLIST#'
     +,           JDIGI,NZEBRA,CUTGAM,IDEBUG,IKINE,NHSTA,P,P,P,P)
      CALL CSCOM('GCMATE,GCNUM,GCONST,GCPHYS,GCTMED,GCTRAK#'
     +,           NMAT,NMATE,PI,IPAIR,NUMED,VECT,P,P,P,P)
      CALL CSCOM('GCUNIT,GCVOLU,GCDRAW,GCKING,GCMULO#',LIN,NLEVEL,NUMNOD
     +,           KCASE,SINMUL,P,P,P,P,P)
*
      CALL CSEXT('GINIT,GZINIT,GDINIT,GPRINT,GPSETS,GXCLOS#'
     +,           GINIT,GZINIT,GDINIT,GPRINT,GPSETS,GXCLOS,P,P,P,P)
      CALL CSEXT(
     +'GSVERT,GSKINE,GSKING,GFIN,GOPEN,GCLOSE,GFOUT#',
     + GSVERT,GSKINE,GSKING,GFIN,GOPEN,GCLOSE,GFOUT,
     +            P,P,P)
      CALL CSEXT('GMATE,GSMATE,GFMATE,GSMIXT,GSTMED,GSTPAR,GPART#'
     +,           GMATE,GSMATE,GFMATE,GSMIXT,GSTMED,GSTPAR,GPART,P,P,P)
      CALL CSEXT('GPIONS,GPHYSI#'
     +,           GPIONS,GPHYSI,P,P,P,P,P,P,P,P)
      CALL CSEXT('GTRIG,GTRIGI,GTRIGC,GTREVE,GIDROP#'
     +,           GTRIG,GTRIGI,GTRIGC,GTREVE,GIDROP,P,P,P,P,P)
      CALL CSEXT('GSVOLU,GSPOS,GSPOSP,GSDVN,GSDVS,GGCLOS,GOPTIM#'
     +,           GSVOLU,GSPOS,GSPOSP,GSDVN,GSDVS,GGCLOS,GOPTIM
     +,           P,P,P)
      CALL CSEXT('GSROTM,GSORD,GSDET,GSDETH,GSDETV,GSATT#'
     +,           GSROTM,GSORD,GSDET,GSDETH,GSDETV,GSATT,P,P,P,P)
      CALL CSEXT('GPLMAT,GSAHIT,GSCHIT,GSDIGI,GSXYZ,GDEBUG#'
     +,           GPLMAT,GSAHIT,GSCHIT,GSDIGI,GSXYZ,GDEBUG,P,P,P,P)
      CALL CSEXT('GPCXYZ,GDCXYZ,GDXYZ,GDAHIT,GDCHIT,GDHITS,GDHEAD#'
     +,           GPCXYZ,GDCXYZ,GDXYZ,GDAHIT,GDCHIT,GDHITS,GDHEAD
     +,           P,P,P)
      CALL CSEXT('GDOPEN,GDCLOS,GDELET,GDRAW,GDRAWC,GDAXIS,GDSCAL#'
     +,           GDOPEN,GDCLOS,GDELET,GDRAW,GDRAWC,GDAXIS,GDSCAL
     +,           P,P,P)
      CALL CSEXT('GDMAN,GDCOL#'
     +,           GDMAN,GDCOL,P,P,P,P,P,P,P,P)
      CALL CSEXT('GDRAWT#',GDRAWT,P,P,P,P,P,P,P,P,P)
      CALL CSEXT('GSCANK,GSCANU,GSCANO,GBRSGE#',GSCANK,GSCANU,GSCANO,
     +            GBRSGE,P,P,P,P,P,P)
      CALL CSEXT('UGLAST#',UGLAST,P,P
     +,           P,P,P,P,P,P,P)
*
      CALL GUXCS
*
      END
 
#endif
