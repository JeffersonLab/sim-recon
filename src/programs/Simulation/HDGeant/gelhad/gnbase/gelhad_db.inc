***** structure gelhad *****
      integer jphadr_gelhad !  Run:  (0=default)
      real ecut_gelhad !  Energy below which no interactions take place
      real scale_gelhad !  Cross section scale factor (1.0=nominal)
      integer mode_gelhad !  GPHAD Model control: (4=default)
      real ethresh_gelhad !  GPHAD Effective pion threshold for mode=4.
*
      common /gelhad_pc/
     + jphadr_gelhad,ecut_gelhad,scale_gelhad,mode_gelhad,
     + ethresh_gelhad
*
