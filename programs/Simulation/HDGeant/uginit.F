*
* $Id$
*
* $Log$
* Revision 1.2  2001/07/15 07:31:38  jonesrt
* HDGeant now supportskinematic input from Monte Carlo generators
* via the routines in hddmInput.c  -rtj
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
      SUBROUTINE UGINIT
*
************************************************************************
*                                                                      *
*              To initialise GEANT3 program and read data cards        *
*                                                                      *
************************************************************************
*
#include "geant321/gckine.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcphys.inc"
*
      integer iskip
      integer infile(20)
      common /inputFile/ iskip,infile
      data iskip/0/
      data infile/20*0/
      integer openInput, skipInput
      external openInput, skipInput
*
*     -----------------------------------------------------------------
*
*             Initialize GEANT
C..geant..
      CALL GINIT
*
*             Prints version number
*
      WRITE(LOUT,1000)
*
*             IKINE  = particle type (default=1=gamma)
*             PKINE(1)=particle energy
*             IKINE and PKINE can be changed with the data card KINE
*
      PKINE(1)=10.
      PKINE(5)=4.
      IKINE=1
      ICOMP=1
      IPAIR=1
      IBREM=1
      IANNI=1
      IPHOT=1
      IHADR=1
*
*             Initialize GEANT/ZBOOK data structures
*
C..geant..
      CALL GZINIT
*
*             Define user FFREAD data cards (format free input)
*
*
*             Read the data cards
*
      OPEN(UNIT=4,FILE='control.in',STATUS='UNKNOWN')
      CALL FFSET('LINP',4)
      call ffkey('file',infile,20,'MIXED')
      call ffkey('skip',iskip,1,'INTEGER')
      CALL GFFGO
*
*             Open the input stream
*
      if (infile(1) .gt. 0) then
        ifail = openInput(infile)
        if (ifail .ne. 0) then
           write(lout,9000) infile
 9000      format('GUKINE ERROR - Could not open input stream ',20a4)
           stop
        endif
        if (iskip .gt. 0) then
          ifail = skipInput(iskip)
        endif
        if (nevent .eq. 0) then
          nevent = 999999999
        endif
      endif
*
*             Initialize graphics package
*
      CALL GDINIT
*
*             Initialize the Hall D geometry
*
      call HDDSgeant3
      call Goptimize
*
*             Initialize GEANT tracking structures
*
C..geant..
      CALL GGCLOS
      CALL GPART
      CALL GPHYSI
*
 1000 FORMAT(/,'  MODE VERSION 1.00 : ',/)
      END
