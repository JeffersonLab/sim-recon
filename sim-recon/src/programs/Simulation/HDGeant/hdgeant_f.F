C
C Note: In order to link with C++ code, "main" must defined from C++.
C Thus, we turn the fortran program "HDGEANT" into a subroutine and
C call it from the C++ main() which is in hdgeant.cc

      SUBROUTINE HDGEANT

      implicit none
      real q,hq
      integer ispace,hspace
      parameter (ispace=5000000)
      parameter (hspace=5000000)
      common /gcbank/ q(ispace)
      common /pawc/ hq(hspace)
      real secmax
      parameter (secmax=300000.)
c      integer istat,icycle

C---- Initialization of HBOOK, ZEBRA, clock
      call GZEBRA(ispace)
      call HLIMIT(-hspace)
      call TIMEST(secmax)

C---- Initialization of GEANT
      call HPLINT(0)
      call UGINIT

C---- Simulation
      call GRUN

C---- Termination ----
      CALL UGLAST
      WRITE(0,999)
  999 FORMAT(//,' ***** Normal exit from Hall D GEANT *****',//)
      END
