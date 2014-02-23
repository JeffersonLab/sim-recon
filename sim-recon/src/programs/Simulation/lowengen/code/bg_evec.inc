C
C---  Simulated event
C
      INTEGER MXTRA
      PARAMETER (MXTRA=10)
      COMMON/BG_EVEC/ IEVPROC,NTRA 
     +               ,ITYPTR(MXTRA),AMTRA(MXTRA),PTRA(3,MXTRA)
     +               ,NDECTR(MXTRA)
      INTEGER IEVPROC  ! the process number
     +       ,NTRA     ! number of particles including the beam and the target
     +       ,ITYPTR   ! track type (GEANT)
     +       ,NDECTR   ! number of decay products
      REAL    AMTRA    ! mass
     +        ,PTRA    ! 3-momentum
C

 
