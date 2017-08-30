      function settofg(vertex,time0)
      real settofg
      real vertex(3)        ! cm
      real time0            ! ns
*
* Sets the Geant variable TOFG which determines the start time of the
* tracking for subsequent particles placed on the primary stack.  The
* start time is determined assuming a beam photon is being generated.
* It is set so that the beam photon will cross the reference plane at
* TOF=time0 if it makes it that far.
*
#include "geant321/gconst.inc"
#include "geant321/gctrak.inc"
#include "controlparams.inc"

      real t0
      real xnormal(2)
      real beam_period_ns
      real reference_time_plane_z
      beam_period_ns = get_beam_period()
      reference_time_plane_z = get_reference_plane()
      if (time0 .eq. 0) then
c smear the time0 value by the trigger time sigma
        call GRANOR(xnormal(1),xnormal(2))
        t0=trigger_time_sigma_ns*xnormal(1)
c discretize the time according to the beam microstructure
        t0=beam_period_ns*floor(t0/beam_period_ns+0.5)
c synchronize the time to the accelerator clock, with the phase
c set to zero when the bunch crosses the reference plane
        TOFG=t0*1e-9 + (vertex(3)-reference_time_plane_z)/CLIGHT
      else
        TOFG=time0*1e-9
      endif
      settofg=TOFG
      end

*
* Read beam_bucket_period_ns from calibdb
*
      function get_beam_period()
      real get_beam_period
      real beam_bucket_period_ns
      common /beam_bucket_saver/beam_bucket_period_ns
      data beam_bucket_period_ns/-1./
      character*80 dbpath
      integer ndata
      integer GetCalib
      external GetCalib
      if (beam_bucket_period_ns .lt. 0) then
        dbpath = "/PHOTON_BEAM/RF/beam_period"
        ndata = 1
        if (GetCalib(dbpath, ndata, beam_bucket_period_ns) .ne. 0) then
          write (6,*) "HDGeant error in settofg: ",
     *                "failed to read RF period ",
     *                "from calibdb, cannot continue."
          stop
        else
          write (6,*) "settofg: beam_bucket_period_ns set to ", 
     *                beam_bucket_period_ns, "ns"
        endif
      endif
      get_beam_period = beam_bucket_period_ns
      end

*
* Read reference_plane_z_cm from calibdb
*
      function get_reference_plane()
      real get_reference_plane
      real reference_plane_z_cm
      common /time_reference_plane_z_saver/reference_plane_z_cm
      data reference_plane_z_cm/1e30/
      character*80 dbpath
      integer ndata
      integer GetCalib
      external GetCalib
      if (reference_plane_z_cm .gt. 1e6) then
        dbpath = "/PHOTON_BEAM/RF/reference_plane_z"
        ndata = 1
        res = GetCalib(dbpath, ndata, reference_plane_z_cm)
        if (res .ne. 0 .or. reference_plane_z_cm .gt. 1e6) then
          reference_plane_z_cm = 65.
          write (6,*) "HDGeant warning in settofg: ",
     *                "failed to read RF reference_plane_z ",
     *                "from calibdb, using default z=65cm."
        else 
          write (6,*) "settofg: RF reference plane set to ", 
     *                reference_plane_z_cm, "cm"
        endif
      endif
      get_reference_plane = reference_plane_z_cm
      end
