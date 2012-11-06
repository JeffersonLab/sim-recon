!Tue Nov 15 17:22:07 PST 1994 - version 1.01 - A. Snyder
!Fortran 90ism removed
!Mon Mar 17 10:57:46 PST 1997 - version 2.00  - A. Snyder
!Modified to |gelhad| version - add cos(theta) dependence
!Fri Apr  4 10:32:29 PST 1997 - version 2.01 - A. Snyder
!Make rho decay angular distribution default to sin^2(theta)
!
      subroutine geltwob(p,k,l)
      implicit none
      save
      real *4 p(5),k(5),l(5)
      real *4 mp,mk,ml,ek,el,pk,pl
      real *4 ctheta,stheta,cphi,sphi,phi,temp
      real *4 pi/3.141592/
      real *4 a/1.0/
      real *4 b/-1.0/
      real *4 value
      character *(*) command
c      integer *4 dum
       mp=p(5)
       mk=k(5)
       ml=l(5)
       ek=(mp**2+mk**2-ml**2)/(2.0*mp)
       el=(mp**2+ml**2-mk**2)/(2.0*mp)
       pk=sqrt(ek**2-mk**2)
       pl=pk
       call gelpicang(a,b,ctheta,stheta)
       call gelrndm(temp,1)
       phi=2.0*pi*temp
       cphi=cos(phi)
       sphi=sin(phi)
       k(4)=ek
       k(1)=pk*stheta*cphi
       k(2)=pk*stheta*sphi
       k(3)=pk*ctheta
       l(4)=el
       l(1)=-k(1)
       l(2)=-k(2)
       l(3)=-k(3)
       call gelboost(k,p,k,-1.0)
       call gelboost(l,p,l,-1.0)
      return
!
      entry geltwobdo(command,value)
      if(command.eq.'set:a') then
         a=value
      else if(command.eq.'set:b') then
         b=value
      else if(command.eq.'tell:a') then
         value=a
      else if(command.eq.'tell:b') then
         value=b
      endif
      return
      end
!
      subroutine gelpicang(a,b,c,s)
      real *4 a,b,c,s,temp(2),max,test
      max=a
      if((max+b).gt.max) max=max+b
 1000 continue
      call gelrndm(temp,2)
      c=-1.0+2.0*temp(1)
      test=(a+b*c**2)/max
      if(temp(2).lt.test) go to 1099
      go to 1000
 1099 continue
      s=sqrt(1.0-c**2)
      return
      end
