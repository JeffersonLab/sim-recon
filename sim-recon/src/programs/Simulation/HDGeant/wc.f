      subroutine wc
      end

      subroutine wcprint(n,x,y,z)
      integer n
      real x(*),y(*),z(*)
c gcvolu.inc
      COMMON/GCVOLU/NLEVEL,NAMES(15),NUMBER(15),
     +LVOLUM(15),LINDEX(15),INFROM,NLEVMX,NLDEV(15),LINMX(15),
     +GTRAN(3,15),GRMAT(10,15),GONLY(15),GLX(3)
c gtvolu.inc
      INTEGER NLEVEL,NAMES,NUMBER,LVOLUM,LINDEX,INFROM,NLEVMX,
     +        NLDEV,LINMX
      REAL GTRAN,GRMAT,GONLY,GLX
      integer getrow,getcolumn,getlayer,getring
      integer getsector,getplane,getmodule,getcell
c      
      real xc(3),ubuf(99)
      real F(3)
      character*20 natmed
      do i=1,n
        xc(1)=x(i)
        xc(2)=y(i)
        xc(3)=z(i)
        call gmedia(xc,numed)
        call gftmed(numed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,
     +              stemax,deemax,epsil,stmin,ubuf,nwbuf)
        call gufld(xc,F)
        print 1000, i,xc,natmed,(names(l),number(l),l=1,nlevel)
        if (getrow().gt.0)    print 1002, 'row',getrow()
        if (getcolumn().gt.0) print 1002, 'column',getcolumn()
        if (getlayer().gt.0)  print 1002, 'layer',getlayer()
        if (getring().gt.0)   print 1002, 'ring',getring()
        if (getsector().gt.0) print 1002, 'sector',getsector()
        if (getplane().gt.0)  print 1002, 'plane',getplane()
        if (getmodule().gt.0) print 1002, 'module',getmodule()
c       if (getcell().gt.0)   print 1002, 'cell',getcell()
        if (ifield.eq.1) then
          print 1010, F
          if (sqrt(F(1)**2+F(2)**2+F(3)**2).gt.fieldm) then
            print *, 'WARNING: local magnetic field exceeds ',
     +               'upper bound of',fieldm,' specified for medium!'
          endif
        elseif (ifield.eq.2) then
          print 1020, F
          if (sqrt(F(1)**2+F(2)**2+F(3)**2).gt.fieldm) then
            print *, 'WARNING: local magnetic field exceeds ',
     +               'upper bound of',fieldm,' specified for medium!'
          endif
        elseif (ifield.eq.3) then
          print 1030, 0,0,fieldm
        endif
      enddo
 1000 format(' point',i3,':',3g15.8,1x,a20,99('/',a4,i3))
 1002 format(' ',a12,'=',i3)
 1010 format(10x,'strongly inhomogeneous field (',2(g12.5,','),
     +       g12.5,') kG')
 1020 format(10x,'inhomogeneous field (',2(g12.5,','),
     +       g12.5,') kG')
 1030 format(10x,'uniform field (',2(g12.5,','),
     +       g12.5,') kG')
      end

      subroutine wcpnorm(n,x,y,z)
      integer n
      real x(*),y(*),z(*)
      real xnorm(3,2),u(3,2),v(3,2)
      integer ierr
      do i=1,n
        xnorm(1,1)=x(i)
        xnorm(2,1)=y(i)
        xnorm(3,1)=z(i)
        call gmedia(xnorm,numed)
        call ggperp(xnorm,u,ierr)
        xnorm(1,2)=xnorm(1,1)+u(1,1)
        xnorm(2,2)=xnorm(2,1)+u(2,1)
        xnorm(3,2)=xnorm(3,1)+u(3,1)
        call GDFR3D(xnorm,2,u,v)
        call IPL(2,u,v)
      enddo
      end

      subroutine wc3dpline(n,x,y,z)
      integer n
      real x(*),y(*),z(*)
      real u(999),v(999)
      real x3d(3,999)
      if (n.gt.999) then
        print *, 'Warning from wc3dpline - cannot plot more than 999'
        print *, 'points in a single polyline, request ignored.'
        return
      endif
      do i=1,n
        x3d(1,i)=x(i)
        x3d(2,i)=y(i)
        x3d(3,i)=z(i)
      enddo
      call GDFR3D(x3d,n,u,v)
      call IPL(n,u,v)
      end
