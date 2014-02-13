      block data halo
      include 'halo.inc'
      data hparX/2.95e1,4.19,-2.40e5,3.53e6,-1.5e-4,1.13e-3,
     +           1.93e2,-4.54e2,-1.59e6,3.53e6,-3.78e-4,1.13e-3,
     +           5.95e2,-1.92e3,-4.97e6,3.52e6,-6.75e-4,1.13e-3/
      data hparY/1.58e1,-1.06e3,-2.26e5,7.06e6,-1.06e-3,5.67e-4,
     +           5.27e1,-7.48e2,-6.21e5,7.06e6,-1.17e-3,5.66e-4,
     +           2.14e2,-6.97e2,-2.10e6,7.06e6,-1.33e-3,5.66e-4/
      end

      real function halox(x,ns)  ! halo fits from JLAB-TN-06-048
      real x
      integer ns
      include 'halo.inc'
      halox = (hparX(1,ns)+hparX(2,ns)*x+hparX(3,ns)*x**2)
     +        +hparX(4,ns)*exp(-0.5*((x-hparX(5,ns))/hparX(6,ns))**2)
      halox = max(halox,0.)
      end

      real function haloy(y,ns)  ! halo fits from JLAB-TN-06-048
      real y
      integer ns
      include 'halo.inc'
      haloy = (hparY(1,ns)+hparY(2,ns)*y+hparY(3,ns)*y**2)
     +        +hparY(4,ns)*exp(-0.5*((y-hparY(5,ns))/hparY(6,ns))**2)
      haloy = max(haloy,0.)
      end

      real function haloxy(x,y,ns)  ! 2D model by R.T. Jones
      real x,y
      integer ns
      real a,b                      ! Note regarding halo normalization:
      real rr0,p0                   !  52% of the integral of this halo
      parameter (a=1.1e-3)          !  intensity function lies outside the
      parameter (b=0.55e-3)         !  nominal 5-sigma ellipse that is often
      parameter (rr0=15.)           !  used to define the boundary of the halo.
      parameter (p0=1.e-3)          !  For example, to generate a 100ppm halo
      real c0(3),c1,c2              !  beyond the (5a,5b) ellipse boundary,
      data c0/2.5,15.,48./          !  sample one halo event for every 5.2e3
      parameter (c1=0)              !  events in the central gaussian.
      parameter (c2=-1/1.3e-4)
      real r,rr,theta
      real p,f
      r=sqrt(x**2+y**2)
      rr=sqrt((x/a)**2+(y/b)**2)
      theta=atan2(y/b,x/a)
      p=1+p0*(rr*exp(-rr/rr0))**6
      f=c0(ns)*(1+c1*r+c2*r**2)
      haloxy=((f/b)*(cos(theta)**2)**p+(f/(2*a))*(sin(theta)**2)**p)
     +      *sqrt(a**2+b**2) ! + 1e6*exp(-0.5*rr**2)
      haloxy = max(haloxy,0.)
      end

      real function haloc(x,y)  ! cut function occludes the central ellipse
      real x,y
      haloc=1
      if ((x/1.1e-3)**2+(y/0.55e-3)**2.lt.25) then
        haloc=0
      endif
      end
      
      subroutine hmake2(idin,idout,n)  ! 2D random point generator
      integer idin,idout,n
      real x,y
      do i=1,n
        call hrndm2(idin,x,y)
        call hfill(idout,x,y,1.)
      enddo
      end
