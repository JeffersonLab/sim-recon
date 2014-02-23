      common/mode/mode
      character *80 mode
c mode=string describing decay mode ... same as *.newdat filename
c
      common/upid/probmin,qtof,qdedx,qcrid
      real *4 probmin
      logical qtof,qdedx,qcrid
c probmin=min probability of good pid
c qtof=use tof ?
c qdedx=use dedx ?
c qcrid=use crid ?
c
      common/cuts/dmk0
      real *4 dmk0
c dmk0=k0 mass cut
c
