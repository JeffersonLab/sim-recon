      implicit none
      real q,hq
      integer ispace,hspace
      parameter (ispace=2000000)
      parameter (hspace=5000000)
      common /gcbank/ q(ispace)
      common /pawc/ hq(hspace)
      real secmax
      parameter (secmax=300000.)
      integer istat,icycle

C---- Initialization of HBOOK, ZEBRA, clock
      call GZEBRA(ispace)
      call HLIMIT(-hspace)
      call TIMEST(secmax)

C---- Initialization of GEANT
      call HPLINT(0)
      call UGINIT

C---- Simulation
      call HROPEN(50,'RZfile','geant.hbook','N',65536,istat)
        call GRUN
        call HROUT(0,icycle,' ')
      call HREND('RZfile')

C---- Termination ----
      CALL UGLAST
      WRITE(0,999)
  999 FORMAT(//,' ***** Normal exit from Hall D GEANT *****',//)
      END
