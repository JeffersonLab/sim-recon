* mcgeom.f - materials and geometry definition for Geant simulator
*
* Updates:
*	March 14, 2001   original version       -rtj
*
* The routines in this file perform the mapping from MCFast detector
* objects to Geant.  The match is not perfect, and some decisions have
* to be made in this file about what Geant wants to see.  All information
* from the mcfast .db files are made available here, and one can use or
* not use it, as appropriate.  To replace the MCFast description with
* your own, you should comment

* The only check made by gpp during the
* conversion from mcfast is that the parameter list from the mcfast db
* matches the one in the common blocks below.  If not, it complains and
* quits.  Otherwise it assumes that the subroutine below is correct and
* will use it.  If it encounters a detector type in the db that is not
* defined below it includes a stub (dummy) subroutine in its output file
* mcfast.f and warns you that the new volume type is not yet functional.
* That way, things will still compile and run, but the new volume will not
* appear to Geant until some work is done.  The dummy subroutine had to be
* cut out of the mcfast.f file, pasted into this file and the appropriate
* snippet of Geant code filled in to define the new detector element.


      subroutine makedetector
      character*40 name
      character*40 geom_id
      common /detector/name,geom_id

* makedetector - declare the root volume HALL for Geant's geometry tree

	real origin(3),KGauss
	common /environ/origin,KGauss
common arguments for calls to gsmate
	integer imate
	character*20 chnama
	real aa(99),zz(99),dens,radl,absl,vbuf(99),nvbuf
	common /cgsmate/ imate,chnama,aa,zz,dens,radl,absl,vbuf,nvbuf
common arguments for calls to gstmed
	integer itmed,nmat,isvol,ifield,nwbuf
	real fieldm,tmaxfd,stemax,deemax,epsil,stmin,ubuf(99)
	character*20 natmed
	common /cgstmed/itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,
     +                  stemax,deemax,epsil,stmin,ubuf,nwbuf
common arguments for calls to gsvolu
	character*4 chname,chshap
	integer nmed,npar,ivolu
	real par(12)
	common /cgsvolu/chname,chshap,nmed,par,npar,ivolu
      imate = 100
      itmed = 1
      natmed = "atmosphere"
      nmat = 15 ! air
      isvol = 0
      ifield = 0
      fieldm = 0
      tmaxfd = 0
      stemax = 0
      deemax = 0
      epsil = 1e-3
      stmin = 0
      nwbuf = 0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     +            deemax,epsil,stmin,ubuf,nwbuf)
      chname = "HALL"
      chshap = "BOX "
      nmed = itmed
      par(1) = 10000
      par(2) = 10000
      par(3) = 10000
      npar = 3
      call GSVOLU(chname,chshap,nmed,par,npar)
      origin(1) = 0
      origin(2) = 0
      origin(3) = 0
      KGauss = 0
      end

      subroutine makeMaterial
      character*40 name
      real a
      real z
      real density
      real radlen
      real abslen
      real collen
      real dedx
      common /Material/name,a,z,density,radlen,abslen,collen,dedx

* makeMaterial - declare a new material to Geant

	real origin(3),KGauss
	common /environ/origin,KGauss
common arguments for calls to gsmate
	integer imate
	character*20 chnama
	real aa(99),zz(99),dens,radl,absl,vbuf(99),nvbuf
	common /cgsmate/ imate,chnama,aa,zz,dens,radl,absl,vbuf,nvbuf
      if (lenocc(name).gt.20) then
        write(6,*) "makeMaterial warning: material name '",
     +             name(1:lenocc(name)),
     +             "' longer than the Geant limit of 20 chars"
      endif
      imate = imate+1
      chnama = name
      aa(1) = a
      zz(1) = z
      dens = density
      radl = radlen
      absl = abslen
      vbuf(1) = collen
      vbuf(2) = dedx
      nvbuf = 2
      call GSMATE(imate,chnama,a,z,density,radlen,abslen,vbuf,nvbuf)
      end

      subroutine makeMixture
      character*40 name
      integer nmat
      character*40 matnames(5)
      real prop(5)
      common /Mixture/name,nmat,matnames,prop

* makeMixture - declare a new material as a mixture of known materials
*	Here I assume that the prop() above specifies mixtures in terms
*	of the proportion by weight, otherwise I don't know how to get
*	the density of the mixture for the general case.

	real origin(3),KGauss
	common /environ/origin,KGauss
common arguments for calls to gsmate
	integer imate
	character*20 chnama
	real aa(99),zz(99),dens,radl,absl,vbuf(99),nvbuf
	common /cgsmate/ imate,chnama,aa,zz,dens,radl,absl,vbuf,nvbuf
      real specvol
      character*20 matname
      integer m,im
      specvol = 0
      do m=1,nmat
        matname = matnames(m)
        do im=101,imate
          call GFMATE(im,chnama,aa(m),zz(m),dens,radl,absl,vbuf,nvbuf)
          if (matname.eq.chnama) goto 2
        end do
        write(6,*) "makeMixture error: undefined material ",
     +             matname(1:lenocc(matname))," in mixture!"
        STOP
    2   continue
        specvol = specvol + prop(m)/dens
      end do
      imate = imate+1
      if (lenocc(name).gt.20) then
        write(6,*) "makeMixture warning: material name '",
     +             name(1:lenocc(name)),
     +             "' longer than the Geant limit of 20 chars"
      endif
      chnama = name
      dens = 1/specvol
      call GSMIXT(imate,chnama,aa,zz,dens,nmat)
      end

      subroutine makeBPipe
      character*40 name
      real rmin
      real rmax
      real z0
      real zlen
      character*40 mat_fill
      real bndrthk(4)
      character*40 matrbnd(4)
      common /BPipe/name,rmin,rmax,z0,zlen,mat_fill,bndrthk,matrbnd

* makeBPipe - defines the beam pipe running through the apparatus

	real origin(3),KGauss
	common /environ/origin,KGauss
common arguments for calls to gsmate
	integer imate
	character*20 chnama
	real aa(99),zz(99),dens,radl,absl,vbuf(99),nvbuf
	common /cgsmate/ imate,chnama,aa,zz,dens,radl,absl,vbuf,nvbuf
common arguments for calls to gstmed
	integer itmed,nmat,isvol,ifield,nwbuf
	real fieldm,tmaxfd,stemax,deemax,epsil,stmin,ubuf(99)
	character*20 natmed
	common /cgstmed/itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,
     +                  stemax,deemax,epsil,stmin,ubuf,nwbuf
common arguments for calls to gsvolu
	character*4 chname,chshap
	integer nmed,npar,ivolu
	real par(12)
	common /cgsvolu/chname,chshap,nmed,par,npar,ivolu
common arguments for calls to gspos
	integer nr,irot
	character*4 chmoth,chonly
	real xc,yc,zc
	common /cgspos/nr,chmoth,xc,yc,zc,irot,chonly
      do i=1,4
        itmed = itmed+1
        natmed = "pipe-"//matrbnd(i)
        do nmat=101,imate
          call GFMATE(nmat,chnama,aa,zz,dens,radl,absl,vbuf,nvbuf)
          if (matrbnd(i).eq.chnama) goto 2
        end do
        write(6,*) "makeBPipe error: undefined material ",
     +             matrbnd(i)(1:lenocc(matrbnd(i)))," in BPipe!"
        STOP
    2   continue
        isvol = 0
        ifield = 3
        fieldm = KGauss
        tmaxfd = 10
        stemax = 0
        deemax = 0
        epsil = 1e-3
        stmin = 0
        nwbuf = 0
        call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     +              deemax,epsil,stmin,ubuf,nwbuf)
      enddo
      itmed = itmed+1
      natmed = "pipe-"//mat_fill
      do nmat=101,imate
        call GFMATE(nmat,chnama,aa,zz,dens,radl,absl,vbuf,nvbuf)
        if (mat_fill.eq.chnama) goto 4
      end do
      write(6,*) "makeBPipe error: undefined material ",
     +           mat_fill(1:lenocc(mat_fill))," in BPipe!"
      STOP
    4 continue
      isvol = 0
      ifield = 3
      fieldm = KGauss
      tmaxfd = 10
      stemax = 0
      deemax = 0
      epsil = 1e-3
      stmin = 0
      nwbuf = 0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     +            deemax,epsil,stmin,ubuf,nwbuf)
      if (lenocc(name).gt.4) then
        write(6,*) "makeBPipe warning: volume '",
     +             name(1:lenocc(name)),
     +             "' longer than the Geant limit of 4 chars"
      endif
      chname = name
      chshap = "TUBE"
      par(1) = rmin
      par(2) = rmax
      par(3) = zlen/2
      npar = 3
      nmed = itmed
      call GSVOLU(chname,chshap,nmed,par,npar)
      nr = 1
      chmoth = "HALL"
      xc = -origin(1)
      yc = -origin(2)
      zc = z0+par(3)-origin(3)
      irot = 1
      chonly = "ONLY"
      call GSPOS(chname,nr,chmoth,xc,yc,zc,irot,chonly)
      chmoth = chname
      if (bndrthk(1).gt.0) then
        nmed = itmed-4
        write(chname,"(a1,i3.3)") "B",nmed
        chshap = "TUBE"
        par(1) = rmin
        par(2) = rmin+bndrthk(1)
        call GSVOLU(chname,chshap,nmed,par,npar)
        call GSPOS(chname,nr,chmoth,xc,yc,zc,irot,chonly)
      endif
      if (bndrthk(2).gt.0) then
        nmed = itmed-3
        write(chname,"(a1,i3.3)") "B",nmed
        chshap = "TUBE"
        par(1) = rmax-bndrthk(2)
        par(2) = rmax
        call GSVOLU(chname,chshap,nmed,par,npar)
        xc = 0
        yc = 0
        zc = 0
        call GSPOS(chname,nr,chmoth,xc,yc,zc,irot,chonly)
      endif
      if (bndrthk(3).gt.0) then
        nmed = itmed-2
        write(chname,"(a1,i3.3)") "B",nmed
        par(1) = rmin
        par(2) = rmax
        par(3) = bndrthk(3)/2
        call GSVOLU(chname,chshap,nmed,par,npar)
        xc = 0
        yc = 0
        zc = z0+par(3)
        call GSPOS(chname,nr,chmoth,xc,yc,zc,irot,chonly)
      endif
      if (bndrthk(4).gt.0) then
        nmed = itmed-1
        write(chname,"(a1,i3.3)") "B",nmed
        par(1) = rmin
        par(2) = rmax
        par(3) = bndrthk(4)/2
        call GSVOLU(chname,chshap,nmed,par,npar)
        xc = 0
        yc = 0
        zc = z0+zlen-par(3)
        call GSPOS(chname,nr,chmoth,xc,yc,zc,irot,chonly)
      endif
      end

      subroutine makeSolenoid
      character*40 name
      real bfield
      real rmin
      real rmax
      real z0
      real zlen
      character*40 mat_fill
      real thick_boun(4)
      character*40 mat_boun(4)
      common /Solenoid/name,bfield,rmin,rmax,z0,zlen,mat_fill,
     + thick_boun,mat_boun

* makeSolenoid - define the solenoidal magnet to Geant

	real origin(3),KGauss
	common /environ/origin,KGauss
common arguments for calls to gsmate
	integer imate
	character*20 chnama
	real aa(99),zz(99),dens,radl,absl,vbuf(99),nvbuf
	common /cgsmate/ imate,chnama,aa,zz,dens,radl,absl,vbuf,nvbuf
common arguments for calls to gstmed
	integer itmed,nmat,isvol,ifield,nwbuf
	real fieldm,tmaxfd,stemax,deemax,epsil,stmin,ubuf(99)
	character*20 natmed
	common /cgstmed/itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,
     +                  stemax,deemax,epsil,stmin,ubuf,nwbuf
common arguments for calls to gsvolu
	character*4 chname,chshap
	integer nmed,npar,ivolu
	real par(12)
	common /cgsvolu/chname,chshap,nmed,par,npar,ivolu
common arguments for calls to gspos
	integer nr,irot
	character*4 chmoth,chonly
	real xc,yc,zc
	common /cgspos/nr,chmoth,xc,yc,zc,irot,chonly
      do i=1,4
        itmed = itmed+1
        natmed = "soln-"//mat_boun(i)
        do nmat=101,imate
          call GFMATE(nmat,chnama,aa,zz,dens,radl,absl,vbuf,nvbuf)
          if (mat_boun(i).eq.chnama) goto 2
        end do
        write(6,*) "makeSolenoid error: undefined material ",
     +             mat_boun(i)(1:lenocc(mat_boun(i)))," in Solenoid!"
        STOP
    2   continue
        isvol = 0
        ifield = 3
        fieldm = bfield*10
        tmaxfd = 10
        stemax = 0
        deemax = 0
        epsil = 1e-3
        stmin = 0
        nwbuf = 0
        call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     +              deemax,epsil,stmin,ubuf,nwbuf)
      enddo
      itmed = itmed+1
      natmed = "soln-"//mat_fill
      do nmat=101,imate
        call GFMATE(nmat,chnama,aa,zz,dens,radl,absl,vbuf,nvbuf)
        if (mat_fill.eq.chnama) goto 4
      end do
      write(6,*) "makeSolenoid error: undefined material ",
     +           mat_fill(1:lenocc(mat_fill))," in Solenoid!"
      STOP
    4 continue
      isvol = 0
      ifield = 3
      fieldm = bfield*10
      tmaxfd = 10
      stemax = 0
      deemax = 0
      epsil = 1e-3
      stmin = 0
      nwbuf = 0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     +            deemax,epsil,stmin,ubuf,nwbuf)
      if (lenocc(name).gt.4) then
        write(6,*) "makeSolenoid warning: volume name '",
     +             name(1:lenocc(name)),
     +             "' longer than the Geant limit of 4 chars"
      endif
      chname = name
      chshap = "TUBE"
      par(1) = rmin
      par(2) = rmax
      par(3) = zlen/2
      npar = 3
      nmed = itmed
      call GSVOLU(chname,chshap,nmed,par,npar)
      nr = 1
      chmoth = "HALL"
      xc = -origin(1)
      yc = -origin(2)
      zc = z0+par(3)-origin(3)
      irot = 1
      chonly = "ONLY"
      call GSPOS(chname,nr,chmoth,xc,yc,zc,irot,chonly)
      chmoth = chname
      if (thick_boun(1).gt.0) then
        nmed = itmed-4
        write(chname,"(a1,i3.3)") "S",nmed
        par(2) = rmin
        par(2) = rmin+thick_boun(1)
        call GSVOLU(chname,chshap,nmed,par,npar)
        xc = 0
        yc = 0
        zc = 0
        call GSPOS(chname,nr,chmoth,xc,yc,zc,irot,chonly)
      endif
      if (thick_boun(2).gt.0) then
        nmed = itmed-3
        write(chname,"(a1,i3.3)") "S",nmed
        par(1) = rmax-thick_boun(2)
        par(2) = rmax
        nmed = itmed-3
        call GSVOLU(chname,chshap,nmed,par,npar)
        xc = 0
        yc = 0
        zc = 0
        call GSPOS(chname,nr,chmoth,xc,yc,zc,irot,chonly)
      endif
      if (thick_boun(3).gt.0) then
        nmed = itmed-2
        write(chname,"(a1,i3.3)") "S",nmed
        par(1) = rmin
        par(2) = rmax
        par(3) = thick_boun(3)/2
        nmed = itmed-2
        call GSVOLU(chname,chshap,nmed,par,npar)
        xc = 0
        yc = 0
        zc = z0+par(3)
        call GSPOS(chname,nr,chmoth,xc,yc,zc,irot,chonly)
      endif
      if (thick_boun(4).gt.0) then
        nmed = itmed-1
        write(chname,"(a1,i3.3)") "S",nmed
        par(1) = rmin
        par(2) = rmax
        par(3) = thick_boun(4)/2
        nmed = itmed-1
        call GSVOLU(chname,chshap,nmed,par,npar)
        xc = 0
        yc = 0
        zc = z0+zlen-par(3)
        call GSPOS(chname,nr,chmoth,xc,yc,zc,irot,chonly)
      endif
      origin(1) = 0
      origin(2) = 0
      origin(3) = z0+zlen/2
      KGauss = fieldm
      end

      subroutine makeDrift
      integer num
      character*40 name
      integer num_anode
      integer num_cathode
      real rmin
      real rmax
      real z0
      real zlen
      character*40 material
      real thick_boun(4)
      character*40 mat_boun(4)
      common /Drift/num,name,num_anode,num_cathode,rmin,rmax,z0,zlen,
     + material,thick_boun,mat_boun

* makeDrift - declare a new cylindrical drift chamber

	real origin(3),KGauss
	common /environ/origin,KGauss
common arguments for calls to gsmate
	integer imate
	character*20 chnama
	real aa(99),zz(99),dens,radl,absl,vbuf(99),nvbuf
	common /cgsmate/ imate,chnama,aa,zz,dens,radl,absl,vbuf,nvbuf
common arguments for calls to gstmed
	integer itmed,nmat,isvol,ifield,nwbuf
	real fieldm,tmaxfd,stemax,deemax,epsil,stmin,ubuf(99)
	character*20 natmed
	common /cgstmed/itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,
     +                  stemax,deemax,epsil,stmin,ubuf,nwbuf
common arguments for calls to gsvolu
	character*4 chname,chshap
	integer nmed,npar,ivolu
	real par(12)
	common /cgsvolu/chname,chshap,nmed,par,npar,ivolu
common arguments for calls to gspos
	integer nr,irot
	character*4 chmoth,chonly
	real xc,yc,zc
	common /cgspos/nr,chmoth,xc,yc,zc,irot,chonly
      do i=1,4
        itmed = itmed+1
        natmed = "drift-"//mat_boun(i)
        do nmat=101,imate
          call GFMATE(nmat,chnama,aa,zz,dens,radl,absl,vbuf,nvbuf)
          if (mat_boun(i).eq.chnama) goto 2
        end do
        write(6,*) "makeDrift error: undefined material ",
     +             mat_boun(i)(1:lenocc(mat_boun(i)))," in Drift!"
        STOP
    2   continue
        isvol = 0
        ifield = 3
        fieldm = KGauss
        tmaxfd = 10
        stemax = 0
        deemax = 0
        epsil = 1e-3
        stmin = 0
        nwbuf = 0
        call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     +              deemax,epsil,stmin,ubuf,nwbuf)
      enddo
      itmed = itmed+1
      natmed = "drift-"//material
      do nmat=101,imate
        call GFMATE(nmat,chnama,aa,zz,dens,radl,absl,vbuf,nvbuf)
        if (material.eq.chnama) goto 4
      end do
      write(6,*) "makeDrift error: undefined material ",
     +           material(1:lenocc(material))," in Drift!"
      STOP
    4 continue
      isvol = 0
      ifield = 3
      fieldm = KGauss
      tmaxfd = 10
      stemax = 0
      deemax = 0
      epsil = 1e-3
      stmin = 0
      nwbuf = 0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     +            deemax,epsil,stmin,ubuf,nwbuf)
      write(chname,"(A2,i2.2)") "DC",num
      chshap = "TUBE"
      par(1) = rmin
      par(2) = rmax
      par(3) = zlen/2
      npar = 3
      nmed = itmed
      call GSVOLU(chname,chshap,nmed,par,npar)
      call GSATT(chname,"NODE",float(nmed))
      nr = 1
      chmoth = "HALL"
      xc = -origin(1)
      yc = -origin(2)
      zc = z0+par(3)-origin(3)
      irot = 1
      chonly = "ONLY"
      call GSPOS(chname,nr,chmoth,xc,yc,zc,irot,chonly)
      chmoth = chname
      if (thick_boun(1).gt.0) then
        nmed = itmed-4
        write(chname,"(a1,i3.3)") "D",nmed
        par(1) = rmin
        par(2) = rmin+thick_boun(1)
        call GSVOLU(chname,chshap,nmed,par,npar)
        par(2) = rmax
        xc = 0
        yc = 0
        zc = 0
        call GSPOS(chname,nr,chmoth,xc,yc,zc,irot,chonly)
      endif
      if (thick_boun(2).gt.0) then
        nmed = itmed-3
        write(chname,"(a1,i3.3)") "D",nmed
        par(1) = rmax-thick_boun(2)
        par(2) = rmax
        call GSVOLU(chname,chshap,nmed,par,npar)
        xc = 0
        yc = 0
        zc = 0
        call GSPOS(chname,nr,chmoth,xc,yc,zc,irot,chonly)
      endif
      if (thick_boun(3).gt.0) then
        nmed = itmed-2
        write(chname,"(a1,i3.3)") "D",nmed
        par(1) = rmin
        par(2) = rmax
        par(3) = thick_boun(3)/2
        call GSVOLU(chname,chshap,nmed,par,npar)
        xc = 0
        yc = 0
        zc = z0+par(3)
        call GSPOS(chname,nr,chmoth,xc,yc,zc,irot,chonly)
      endif
      if (thick_boun(4).gt.0) then
        nmed = itmed-1
        write(chname,"(a1,i3.3)") "D",nmed
        par(1) = rmin
        par(2) = rmax
        par(3) = thick_boun(4)/2
        call GSVOLU(chname,chshap,nmed,par,npar)
        xc = 0
        yc = 0
        zc = z0+zlen-par(3)
        call GSPOS(chname,nr,chmoth,xc,yc,zc,irot,chonly)
      endif
      end

      subroutine makeLayerDRFAno
      integer det
      integer lyr
      real radius
      real zlen
      real cell_height
      integer nwires
      integer ID_readout
      integer ID_cathode
      real phi0
      real stereo_tau
      real stereo_offset
      real eff_hit
      real eff_dedx
      real siga
      real sigb
      real sigc
      common /LayerDRFAno/det,lyr,radius,zlen,cell_height,nwires,
     + ID_readout,ID_cathode,phi0,stereo_tau,stereo_offset,eff_hit,
     + eff_dedx,siga,sigb,sigc

* makeLayerDRFAno - declare an anode layer for cylindrical drift chamber

	real origin(3),KGauss
	common /environ/origin,KGauss
common arguments for calls to gsmate
	integer imate
	character*20 chnama
	real aa(99),zz(99),dens,radl,absl,vbuf(99),nvbuf
	common /cgsmate/ imate,chnama,aa,zz,dens,radl,absl,vbuf,nvbuf
common arguments for calls to gstmed
	integer itmed,nmat,isvol,ifield,nwbuf
	real fieldm,tmaxfd,stemax,deemax,epsil,stmin,ubuf(99)
	character*20 natmed
	common /cgstmed/itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,
     +                  stemax,deemax,epsil,stmin,ubuf,nwbuf
common arguments for calls to gsvolu
	character*4 chname,chshap
	integer nmed,npar,ivolu
	real par(12)
	common /cgsvolu/chname,chshap,nmed,par,npar,ivolu
common arguments for calls to gspos
	integer nr,irot
	character*4 chmoth,chonly
	real xc,yc,zc
	common /cgspos/nr,chmoth,xc,yc,zc,irot,chonly
      write(chmoth,"(a2,i2.2)") "DC",det
      call GFATT(chmoth,"NODE",float(nmed))
      write(chname,"(a1,i3.3)") "A",det*10+lyr
      chshap = "TUBE"
      par(1) = radius-cell_height/2
      par(2) = radius+cell_height/2
      par(3) = zlen/2
      npar = 3
      call GSVOLU(chname,chshap,nmed,par,npar)
      nr = 1
      xc = 0
      yc = 0
      zc = 0
      chonly = "ONLY"
      call GSPOS(chname,nr,chmoth,xc,yc,zc,irot,chonly)
      chmoth = chname
      write(chname,"(a1,i3.3)") "L",det*10+lyr
      call GSDVN2(chname,chmoth,nwires,2,phi0,nmed)
      end

      subroutine makeOffsetDRFAno
      integer det
      integer lyr
      real cell_offset
      real sag
      real offset(3)
      real dircos(3)
      common /OffsetDRFAno/det,lyr,cell_offset,sag,offset,dircos

CC---> add the appropriate GEANT calls here

      end

      subroutine makeLayerDRFCatho
      integer det
      integer lyr
      real delta_r
      real zlen
      integer nstrips
      integer n_phi_segm
      integer ID_anode
      integer cell_offset
      real eff_hit
      real resa
      real resb
      real resc
      common /LayerDRFCatho/det,lyr,delta_r,zlen,nstrips,n_phi_segm,
     + ID_anode,cell_offset,eff_hit,resa,resb,resc

CC---> add the appropriate GEANT calls here

      end

      subroutine makeAbsorber
      character*40 name
      character*40 shape
      integer type
      real rmin(2)
      real rmax(2)
      real z0
      real zlen
      character*40 material
      common /Absorber/name,shape,type,rmin,rmax,z0,zlen,material

* makeAbsorber - declare a region containing passive material

	real origin(3),KGauss
	common /environ/origin,KGauss
common arguments for calls to gsmate
	integer imate
	character*20 chnama
	real aa(99),zz(99),dens,radl,absl,vbuf(99),nvbuf
	common /cgsmate/ imate,chnama,aa,zz,dens,radl,absl,vbuf,nvbuf
common arguments for calls to gstmed
	integer itmed,nmat,isvol,ifield,nwbuf
	real fieldm,tmaxfd,stemax,deemax,epsil,stmin,ubuf(99)
	character*20 natmed
	common /cgstmed/itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,
     +                  stemax,deemax,epsil,stmin,ubuf,nwbuf
common arguments for calls to gsvolu
	character*4 chname,chshap
	integer nmed,npar,ivolu
	real par(12)
	common /cgsvolu/chname,chshap,nmed,par,npar,ivolu
common arguments for calls to gspos
	integer nr,irot
	character*4 chmoth,chonly
	real xc,yc,zc
	common /cgspos/nr,chmoth,xc,yc,zc,irot,chonly
      itmed = itmed+1
      natmed = material
      do nmat=101,imate
        call GFMATE(nmat,chnama,aa,zz,dens,radl,absl,vbuf,nvbuf)
        if (natmed.eq.chnama) goto 2
      end do
      write(6,*) "makeAbsorber error: undefined material ",
     +           natmed(1:lenocc(natmed))," in Absorber!"
      STOP
    2 continue
      isvol = 0
      ifield = 3
      fieldm = KGauss
      tmaxfd = 10
      stemax = 0
      deemax = 0
      epsil = 1e-3
      stmin = 0
      nwbuf = 0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     +            deemax,epsil,stmin,ubuf,nwbuf)
      nmed = itmed
      if (lenocc(name).gt.4) then
        write(6,*) "makeAbsorber warning: volume name '",
     +             name(1:lenocc(name)),
     +             "' longer than the Geant limit of 4 chars"
      endif
      chname = name
      chshap = shape
      par(1) = rmin(1)
      par(2) = rmax(1)
      par(3) = zlen/2
      npar = 3
      call GSVOLU(chname,chshap,nmed,par,npar)
      nr = 1
      chmoth = "LASS"
      xc = -origin(1)
      yc = -origin(2)
      zc = z0+par(3)-origin(3)
      irot = 1
      chonly = "ONLY"
      call GSPOS(chname,nr,chmoth,xc,yc,zc,irot,chonly)
      end

      subroutine makeSiDisk
      integer num
      character*40 name
      integer nlyr
      real zpos
      common /SiDisk/num,name,nlyr,zpos

* makeSiDisk - declare a new disk-shaped detector (not necessarily silicon)

	real origin(3),KGauss
	common /environ/origin,KGauss
common arguments for calls to gsvolu
	character*4 chname,chshap
	integer nmed,npar,ivolu
	real par(12)
	common /cgsvolu/chname,chshap,nmed,par,npar,ivolu
common arguments for calls to gspos
	integer nr,irot
	character*4 chmoth,chonly
	real xc,yc,zc
	common /cgspos/nr,chmoth,xc,yc,zc,irot,chonly
      write(chname,"(a3,i1)") "FDC",num
      chshap = "TUBE"
      npar = 0
      nmed = 1   ! atmosphere
      call GSVOLU(chname,chshap,nmed,par,npar)
      call GSATT(chname,"NODE",zpos)
      end

      subroutine makeLayerSiDi
      integer det
      integer lyr
      character*40 mat
      integer nwed
      real z_local
      real thick
      real rmin
      real rmax
      real phi(2)
      real dphi
      integer type
      common /LayerSiDi/det,lyr,mat,nwed,z_local,thick,rmin,rmax,phi,
     + dphi,type

* makeLayerSiDi - declare a layer for disk-shaped volume (not only silicon)

	real origin(3),KGauss
	common /environ/origin,KGauss
common arguments for calls to gsmate
	integer imate
	character*20 chnama
	real aa(99),zz(99),dens,radl,absl,vbuf(99),nvbuf
	common /cgsmate/ imate,chnama,aa,zz,dens,radl,absl,vbuf,nvbuf
common arguments for calls to gstmed
	integer itmed,nmat,isvol,ifield,nwbuf
	real fieldm,tmaxfd,stemax,deemax,epsil,stmin,ubuf(99)
	character*20 natmed
	common /cgstmed/itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,
     +                  stemax,deemax,epsil,stmin,ubuf,nwbuf
common arguments for calls to gsvolu
	character*4 chname,chshap
	integer nmed,npar,ivolu
	real par(12)
	common /cgsvolu/chname,chshap,nmed,par,npar,ivolu
common arguments for calls to gspos
	integer nr,irot
	character*4 chmoth,chonly
	real xc,yc,zc
	common /cgspos/nr,chmoth,xc,yc,zc,irot,chonly
      write(chmoth,"(a3,i1)") "FDC",det
      call GFATT(chmoth,"NODE",zc)
      if (zc.ne.999.99) then
        call GSATT(chmoth,"NODE",999.99)
        nr = 1
        xc = -origin(1)
        yc = -origin(2)
        zc = zc-origin(3)
        irot = 1
        chonly = "ONLY"
        npar = 3
        call GSPOSP(chmoth,nr,"LASS",xc,yc,zc,irot,chonly,par,npar)
      endif
      itmed = itmed+1
      natmed = "sidi-"//mat
      do nmat=101,imate
        call GFMATE(nmat,chnama,aa,zz,dens,radl,absl,vbuf,nvbuf)
        if (mat.eq.chnama) goto 2
      end do
      write(6,*) "makeLayerSiDi error: undefined material ",
     +           mat(1:lenocc(natmed))," in LayerSiDi!"
      STOP
    2 continue
      isvol = 0
      ifield = 3
      fieldm = KGauss
      tmaxfd = 10
      stemax = 0
      deemax = 0
      epsil = 1e-3
      stmin = 0
      nwbuf = 0
      call GSTMED(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     +            deemax,epsil,stmin,ubuf,nwbuf)
      write(chname,"(a1,i3.3)") "P",det*10+lyr
      chshap = "TUBE"
      par(1) = rmin
      par(2) = rmax
      par(3) = thick/2
      npar = 3
      nmed = itmed
      call GSVOLU(chname,chshap,nmed,par,npar)
      nr = 1
      xc = 0
      yc = 0
      zc = z_local
      irot = 1
      chonly = "ONLY"
      call GSPOS(chname,nr,chmoth,xc,yc,zc,irot,chonly)
      chmoth = chname
      write(chname,"(a1,i3.3)") "W",nmed
      step = (phi(2)-phi(1))/nwed *180/3.1416
      phi0 = dphi *180/3.1416
      call GSDVT2(chname,chmoth,step,phi0,nmed,0)
      end

      subroutine makeWedge
      character*40 speci
      integer det
      integer lyr
      integer nwed
      integer nstrip
      real c0_r
      real c0_phi
      real pitch
      real stereo
      real eff_hit
      real siga
      real sigb
      real sigc
      common /Wedge/speci,det,lyr,nwed,nstrip,c0_r,c0_phi,pitch,stereo,
     + eff_hit,siga,sigb,sigc

CC---> add the appropriate GEANT calls here

      end

      subroutine makeHitsOnTrack
      integer all
      integer z
      integer svx
      common /HitsOnTrack/all,z,svx

CC---> add the appropriate GEANT calls here

      end

      subroutine makeBeamVrtx
      real xyz(3)
      real sig(3)
      common /BeamVrtx/xyz,sig

CC---> add the appropriate GEANT calls here

      end

      subroutine makeGeometry
      call detectorDef
      call MaterialDef
      call MixtureDef
      call SolenoidDef
      call BPipeDef
      call DriftDef
      call LayerDRFAnoDef
      call OffsetDRFAnoDef
      call LayerDRFCathoDef
      call AbsorberDef
      call SiDiskDef
      call LayerSiDiDef
      call WedgeDef
      call HitsOnTrackDef
      call BeamVrtxDef
      end
