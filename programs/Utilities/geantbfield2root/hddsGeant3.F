      function getMap()
      implicit none
      integer getMap
      integer nlevel,names,number,lvolum
      common /gcvolu/nlevel,names(15),number(15),lvolum(15)
      integer i,istart(215)
      data (istart(i),i=1,100) /
     +     0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     +     0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     +     0,    0,    0,    0,    0,    0,    0,    0,    0,    1,
     +     0,    0,    0,    0,    0,    0,    0,    2,    0,    0,
     +     0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     +     0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     +     3,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     +     0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     +     0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     +     0,    0,    0,    0,    0,    0,    0,    0,    0,    0/
      data (istart(i),i=101,200) /
     +     0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     +     0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     +     0,    0,    0,    0,    0,    0,    0,    4,    0,    0,
     +     5,    0,    0,    6,    0,    0,    7,    0,    0,    8,
     +     0,    0,    9,    0,    0,   10,    0,    0,   11,    0,
     +     0,   12,    0,    0,   13,    0,    0,   14,    0,    0,
     +    15,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     +     0,    0,    0,    0,    0,    0,    0,   16,    0,    0,
     +     0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     +     0,    0,    0,    0,    0,    0,    0,    0,    0,    0/
      data (istart(i),i=201,215) /
     +     0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     +     0,    0,    0,    0,    0/
      integer lookup(16)
      data (lookup(i),i=1,16) /
     +     1,    2,    3,    4,    4,    4,    4,    4,    4,    4,
     +     4,    4,    4,    4,    4,    5/
      integer level,index
      integer map
      getMap = 0
      do level=1,nlevel
        index = istart(lvolum(level))
        if (index.gt.0) then
          map = lookup(index + number(level) - 1)
          if (map.gt.0) then
            getMap = map
          endif
        endif
      enddo
		
		getMap=3
		
      end

      subroutine gufld(r,B)
      implicit none
      real r(3),B(3)
      real rr(3),rs(3),BB(3)
      integer iregion
      integer getMap
      external getMap

      real orig1(3),rmat1(3,3)
      data orig1/0,0,-1975/
      data rmat1/1,0,0,
     +          0,6.12303e-17,-1,
     +          0,1,6.12303e-17/
      real orig2(3),rmat2(3,3)
      data orig2/0,0,-1345/
      data rmat2/1,0,0,
     +          0,6.12303e-17,-1,
     +          0,1,6.12303e-17/
      real orig3(3),rmat3(3,3)
      data orig3/0,0,0/
      data rmat3/1,0,0,
     +          0,1,0,
     +          0,0,1/
      real orig4(3),rmat4(3,3)
      data orig4/0,0,212/
      data rmat4/1,0,0,
     +          0,1,0,
     +          0,0,1/
      real orig5(3),rmat5(3,3)
      data orig5/0,0,622.8/
      data rmat5/1,0,0,
     +          0,1,0,
     +          0,0,1/

      iregion = getMap()
      if (iregion.eq.0) then
        B(1) = 0
        B(2) = 0
        B(3) = 0
      else if (iregion.eq.1) then
        B(1) = 0
        B(2) = -2
        B(3) = 1.22461e-16
      else if (iregion.eq.2) then
        B(1) = 0
        B(2) = -2
        B(3) = 1.22461e-16
      else if (iregion.eq.3) then
        rs(1) = r(1)-orig3(1)
        rs(2) = r(2)-orig3(2)
        rs(3) = r(3)-orig3(3)
        rr(1) = rs(1)*rmat3(1,1)+rs(2)*rmat3(1,2)+rs(3)*rmat3(1,3)
        rr(2) = rs(1)*rmat3(2,1)+rs(2)*rmat3(2,2)+rs(3)*rmat3(2,3)
        rr(3) = rs(1)*rmat3(3,1)+rs(2)*rmat3(3,2)+rs(3)*rmat3(3,3)
        call gufld1(rr,BB)
        B(1) = BB(1)*rmat3(1,1)+BB(2)*rmat3(2,1)+BB(3)*rmat3(3,1)
        B(2) = BB(1)*rmat3(1,2)+BB(2)*rmat3(2,2)+BB(3)*rmat3(3,2)
        B(3) = BB(1)*rmat3(1,3)+BB(2)*rmat3(2,3)+BB(3)*rmat3(3,3)
        B(1) = B(1)*10
        B(2) = B(2)*10
        B(3) = B(3)*10
      endif
      end

      subroutine gufld1(r,B)
      implicit none
      real r(3),B(3)
      real rho,phi,alpha
      real u(3)
      real twopi
      parameter (twopi=6.28318530717959)

      real bound0(3,2)
      data bound0/0,0,10.16,
     +            101.6,6.28319,645.16/
      integer reverse0(3)
      data reverse0/1,1,1/
      real Bmap(3,41,1,251)
      integer nsites(3)
      data nsites/41,1,251/
      logical loaded
      data loaded/.false./
      save Bmap,nsites,loaded
      integer i,i1,i2,i3

      if (.not.loaded) then
        open(unit=78,file='solenoid.map',status='old',err=7)
        read(unit=78,fmt=*,err=5,end=6)
     +      ((((Bmap(i,i1,i2,i3),i=1,3),
     +         i3=1,251),
     +        i1=1,41),
     +       i2=1,1)
        go to 8
    5   stop 'error reading magnetic field map, stop'
    6   stop 'EOF encountered reading magnetic field map, stop'
    7   stop 'error opening magnetic field map, stop'
    8   loaded=.true.
      endif

      rho = sqrt(r(1)**2+r(2)**2)
      phi = atan2(r(2),r(1))
      u(1) = (rho-bound0(1,1))/(bound0(1,2)-bound0(1,1))
      u(2) = (phi-bound0(2,1))/(bound0(2,2)-bound0(2,1))
      alpha = abs(twopi/(bound0(2,2)-bound0(2,1)))
      u(2) = u(2)-int(u(2)/alpha)*alpha
      if (u(2).lt.0) then
        u(2) = u(2)+alpha
      endif
      u(3) = (r(3)-bound0(3,1))/(bound0(3,2)-bound0(3,1))
      if ((u(1).ge.0.and.u(1).le.1).and.
     +    (u(2).ge.0.and.u(2).le.1).and.
     +    (u(3).ge.0.and.u(3).le.1)) then
        call interpol3(Bmap,nsites,u,B)
        B(1)=B(1)*reverse0(1)
        B(2)=B(2)*reverse0(2)
        B(3)=B(3)*reverse0(3)
        return
      endif

      B(1) = 0
      B(2) = 0
      B(3) = 0
      return
      end

      subroutine interpol3(Bmap,nsites,u,B)
      implicit none
      integer nsites(3)
      real Bmap(3,nsites(1),nsites(2),nsites(3))
      real u(3),B(3)
      integer ir(3),ir0(3),ir1(3)
      real ur(3),dur(3),ugrad(3,3)
      integer i
      do i=1,3
        ur(i)=u(i)*(nsites(i)-1)+1
        ir(i)=nint(ur(i))
        ir0(i)=max(ir(i)-1,1)
        ir1(i)=min(ir(i)+1,nsites(i))
        dur(i)=(ur(i)-ir(i))/(ir1(i)-ir0(i)+1e-20)
      enddo
      ugrad(1,1)=(Bmap(1,ir1(1),ir(2),ir(3))-Bmap(1,ir0(1),ir(2),ir(3)))
      ugrad(2,1)=(Bmap(2,ir1(1),ir(2),ir(3))-Bmap(2,ir0(1),ir(2),ir(3)))
      ugrad(3,1)=(Bmap(3,ir1(1),ir(2),ir(3))-Bmap(3,ir0(1),ir(2),ir(3)))
      ugrad(1,2)=(Bmap(1,ir(1),ir1(2),ir(3))-Bmap(1,ir(1),ir0(2),ir(3)))
      ugrad(2,2)=(Bmap(2,ir(1),ir1(2),ir(3))-Bmap(2,ir(1),ir0(2),ir(3)))
      ugrad(3,2)=(Bmap(3,ir(1),ir1(2),ir(3))-Bmap(3,ir(1),ir0(2),ir(3)))
      ugrad(1,3)=(Bmap(1,ir(1),ir(2),ir1(3))-Bmap(1,ir(1),ir(2),ir0(3)))
      ugrad(2,3)=(Bmap(2,ir(1),ir(2),ir1(3))-Bmap(2,ir(1),ir(2),ir0(3)))
      ugrad(3,3)=(Bmap(3,ir(1),ir(2),ir1(3))-Bmap(3,ir(1),ir(2),ir0(3)))
      B(1)=Bmap(1,ir(1),ir(2),ir(3))
     +       +ugrad(1,1)*dur(1)+ugrad(1,2)*dur(2)+ugrad(1,3)*dur(3)
      B(2)=Bmap(2,ir(1),ir(2),ir(3))
     +       +ugrad(2,1)*dur(1)+ugrad(2,2)*dur(2)+ugrad(2,3)*dur(3)
      B(3)=Bmap(3,ir(1),ir(2),ir(3))
     +       +ugrad(3,1)*dur(1)+ugrad(3,2)*dur(2)+ugrad(3,3)*dur(3)
      end
