!Tue Nov 15 17:22:07 PST 1994 - version 1.01 - A. Snyder
!Remove fortran 90ism
!Mon Mar 17 10:57:46 PST 1997 - A. Snyder
!Modified to be |gelhad| version
!
      subroutine gelboost(pin,to,pout,sign)
      implicit none
      real *4 pin(5),to(5),pout(5),sign
      real *4 pto,work,pinp,ein,poutp,eout,pinx,piny,poutx,pouty
      real *4 x(3),y(3),z(3)
      real *4 beta,gamma
      logical ok
!
       pto=sqrt(to(1)**2+to(2)**2+to(3)**2)
       gamma=to(4)/to(5)
       beta=pto/to(4)
       z(1)=to(1)/pto
       z(2)=to(2)/pto
       z(3)=to(3)/pto
       work=sqrt(z(1)**2+z(2)**2)
       if(work.gt.0.1) then
          x(1)=z(1)*z(3)/work
          x(2)=z(2)*z(3)/work
          x(3)=-work
       else if(work.eq.0.0) then
          x(1)=1.0
          x(2)=0.0
          x(3)=0.0
       else
          x(1)=1.0
          x(2)=0.0
          x(3)=0.0
          call gelperp(x,z,x,ok)
          if(.not.ok) then
             x(1)=0.0
             x(2)=1.0
             x(3)=0.0
             call gelperp(x,z,x,ok)
          endif
       endif
       y(1)=z(2)*x(3)-z(3)*x(2)
       y(2)=z(3)*x(1)-z(1)*x(3)
       y(3)=z(1)*x(2)-z(2)*x(1)
       ein=pin(4)
       pinp=z(1)*pin(1)+z(2)*pin(2)+z(3)*pin(3)
       pinx=x(1)*pin(1)+x(2)*pin(2)+x(3)*pin(3)
       piny=y(1)*pin(1)+y(2)*pin(2)+y(3)*pin(3)
       poutx=pinx
       pouty=piny
       poutp=gamma*(pinp-sign*beta*ein)
       eout=gamma*(ein-sign*beta*pinp)
       pout(1)=poutp*z(1)+poutx*x(1)+pouty*y(1)
       pout(2)=poutp*z(2)+poutx*x(2)+pouty*y(2)
       pout(3)=poutp*z(3)+poutx*x(3)+pouty*y(3)
       pout(4)=eout
      return
      end
c
      subroutine geldot(u,v,dotp)
      implicit none
      real *4 u(3),v(3),dotp
      dotp=u(1)*v(1)+u(2)*v(2)+u(3)*v(3)
      return
      end
c
      subroutine gelperp(x,z,xp,ok)
      implicit none
      real *4 x(3),z(3),xp(3),xdotz,temp
      logical ok
      call geldot(x,z,xdotz)
      xp(1)=x(1)-xdotz*z(1)
      xp(2)=x(2)-xdotz*z(2)
      xp(3)=x(3)-xdotz*z(3)
      call geldot(xp,xp,temp)
      temp=sqrt(temp)
      ok=temp.gt.0.1
      return
      end
