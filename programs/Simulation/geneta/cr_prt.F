C
C   12/22/09  sross section with Lage/Sergey and 3 variable parameters
C
C     f77 laget_eta.f -o laget -L $CERN/2004/lib -lpacklib -lmathlib -lkernlib
C
      SUBROUTINE CRSEC6(egama,thpi0,sumall)
C
C      program-laget_eta ! dsigma/dt(microbarn/Gev^2)
c	Labels: o-omega contribution,r-rho meson,c-coulomb
      implicit none
      double precision meta, mn, ebeam
      parameter(meta = 0.54775d0, mn = 0.938272d0)
      real dgamma
      real egama,sumall
      double precision thpi0
      double precision s, s0, t,ao,ar,alpha,phi
      common/spar_com/ s, s0, t,ao,ar,alpha,phi
      double precision pi,co,cr,cc,w,wo,wr,wc
      double precision wdomeg
      double precision go,gr,fo,fr,fc
      double precision par(3),dwidth
C
	parameter(pi = 3.1416d0)
C                    PDG after 2004
Ct	parameter(dwidth = 0.51d0)
C                    PDG before 2004
	parameter(dwidth = 0.46d0)
C
       integer n
      double precision tmin, theta
C
C       !I put here the phase from Cornell, Sergey
Ct      phi=pi/3.d0     ! not correct, Ashot
Ct      phi=1.0d0
Ct      phi=-0.35d0
C
      phi=0.0d0
C
	par(1) = dwidth   ! magn. of Primakoff process
	par(2) = 1.0d0    ! magn. of hadr. amplitude
	par(3) = phi      ! diff. of interf. phase angle
C
	ebeam = egama
C
C       do n=0,100
C      theta=n*0.03d0
*      thetal = theta  / sqrt((mn+2.d0*ebeam)/mn)
C
      theta = thpi0
C
      s0 = 1.0d0
      s  = 2.d0*mn*ebeam + mn**2
      
      co=0.063d0     !0.29*sqrt(6.44/137)
      cr=0.0664d0    !0.81*sqrt(0.92/137)
	 
      tmin = meta**4/(4.d0*ebeam**2)
      t = -4.d0*ebeam**2*dsin(0.5d0*theta*pi/180.d0)**2-tmin
      ao=0.44d0+0.9d0*t
      ar=0.55d0+0.8d0*t
		
		if((ao.lt.0.0).OR.(ar.lt.0.0)) then
			sumall = 0
			RETURN
		endif

      go=dgamma(ao)      !Euler Gamma functions
      gr=dgamma(ar)      !Euler Gamma functions
C                           the omega amplitude
      fo= par(2)*20.0d0*co*ebeam**2*s**(ao-1.d0)*dsin(theta*pi/180.d0)*
     &(pi*0.9d0/(go*dsin(pi*ao)))/(1.414d0*meta)
       wo=fo**2
      fr= par(2)*20.0d0*cr*ebeam**2*s**(ar-1.d0)*dsin(theta*pi/180.d0)*
     &(pi*0.8d0/(gr*dsin(pi*ar)))/(1.414d0*meta)
      wr=fr**2

	fc=-dsqrt((8.0d0*par(1))/(137.0d0*meta**3))*ebeam**2*
     &dsin(theta*pi/180.0d0)/(t*50.0d0)
	wc=fc**2
	w=pi*(wo+wr+wc+2.0d0*dcos(pi*(ar-ao))*fo*fr+
     &2*fc*(fo*dcos(pi*ao+par(3))+fr*dcos(pi*ar+par(3))))/ebeam**2
C
C         !to change to 
c      dsigma/domega drop out pi at front and ebeam**2 at end of 
C      this expression.
C
*       from Ashot, calculation of d(sigma)/d(omega)
C
	wdomeg = (ebeam*ebeam*w)/pi
C
C	OPEN(2,FILE='laget_eta.dat')
C      write(2,1) -t,w
C	write(*,1) -t,w
C    1 format(x,7(E10.4,x))
C
C      enddo
C
        sumall = SNGL(wdomeg/1000.0)
C
      return
      end

     

      
