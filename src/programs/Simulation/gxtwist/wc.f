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
 1000 format(' point',i3,':',3g12.5,1x,a20,99('/',a4,i3))
 1010 format(10x,'strongly inhomogeneous field (',2(g12.5,','),
     +       g12.5,') kG')
 1020 format(10x,'inhomogeneous field (',2(g12.5,','),
     +       g12.5,') kG')
 1030 format(10x,'uniform field (',2(g12.5,','),
     +       g12.5,') kG')
      end
