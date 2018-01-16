C
C---  A binary reaction definition: gamma+A-->B(recoil)+C
C
      INTEGER MXPNTR
      PARAMETER (MXPNTR=3000)
      COMMON/BG_REAC/  IPREAC(2),TSLREAC,ELREAC(2),NPXREAC
     +                ,XSREAC(MXPNTR),ESREAC(MXPNTR)
C
      INTEGER IPREAC      ! (1:2) - GEANT types of the secondary particles (=0 - photon) (F)
     +       ,NPXREAC     ! number of points in the array of energies/cross_sections (F)
      REAL    TSLREAC     ! t-slope 1/(GeV/C)**2, ds/dt~exp(-TSLREAC*t) (F)
     +       ,ELREAC      ! (1:2) min,max of the beam energy (F)
     +       ,ESREAC      ! energies ELREAC(1)+(ELREAC(2)-ELREAC(1))/(NPXREAC-1)*i 
     +       ,XSREAC      ! (1:NPXREAC) cross sections in the points ESREAC (F)
C
