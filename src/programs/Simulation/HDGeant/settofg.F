      function settofg(vertex,time0)
      real settofg
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
#include "controlparams.inc"

      real reference_time_plane_z
      parameter (reference_time_plane_z=65.0)
      real beam_bucket_period_ns
      parameter (beam_bucket_period_ns=2.0)
      real xnormal(2)
      real t0

c smear the time0 value by the trigger time sigma
      call GRANOR(xnormal(1),xnormal(2))
      t0 = time0 + trigger_time_sigma_ns*xnormal(1)
c discretize the time according to the beam microstructure
      TOFG=beam_bucket_period_ns*int(t0/beam_bucket_period_ns+0.5)
c synchronize the time to the accelerator clock, with the phase
c set to zero when the bunch crosses the reference plane
      TOFG=TOFG*1e-9 + (vertex(3)-reference_time_plane_z)/CLIGHT
      settofg=TOFG
      end
