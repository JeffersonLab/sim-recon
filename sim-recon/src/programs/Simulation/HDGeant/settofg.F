      subroutine settofg(vertex,time0)
      real vertex(3)        ! cm
      real time0            ! ns
*
* Sets the Geant variable TOFG which determines the start time of the
* tracking for subsequent particles placed on the primary stack.  The
* start time is determined assuming a beam photon is being generated.
* It is set so that the photon will cross the reference time plane at
* TOF=t0 if it makes it that far.
*
#include "geant321/gconst.inc"
#include "geant321/gctrak.inc"

      real reference_time_plane_z
      parameter (reference_time_plane_z=65.0)
      real beam_bucket_period_ns
      parameter (beam_bucket_period_ns=2.0)

      TOFG=2*int(time0/beam_bucket_period_ns+0.5)*1e-9
     +    +(vertex(3)-reference_time_plane_z)/CLIGHT
      end
