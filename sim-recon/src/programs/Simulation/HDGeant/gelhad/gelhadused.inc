c |gelhadused| is set to |.true.| when |gphad| is called
      integer   lmecgelh        ! GELHAD mechanism name
      parameter (lmecgelh = 212)
      common/gelhadused/gelhadused
     $  ,nevtot_gelh
     $  ,nevhad_gelh
     $  ,nevout_gelh
     $  ,ngelhperev
     $  ,ninthad_gelh
      logical gelhadused  ! there was an interaction in this event
      integer nevtot_gelh ! total events
      integer nevhad_gelh ! total hadronic interaction events
      integer nevout_gelh ! total events output
      integer ngelhperev  ! Number of gelhad interactions per event
      integer ninthad_gelh! total interactions
