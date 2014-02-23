!Tue Mar  4 11:21:49 PST 1997 - A. Snyder
!Boost produced particles from recoil frame pack to lab
      subroutine labframe(u,mrecoil,precoil,np,p,id)
      implicit none
      real *4 mrecoil,precoil,erecoil
      integer *4 np,id(*),i
      real *4 p(3,*),e,pznew,u(3),pz
      real *4 gamma,beta
c
      character *20 name
      real *4 temp(10)
      integer *4 ntemp,type
      real *4 mass,charge,lifetime
c
      erecoil=sqrt(precoil**2+mrecoil**2)
      beta=precoil/erecoil
      gamma=erecoil/mrecoil
c
      do i=1,np
       call gfpart(id(i),name,type,mass,charge,lifetime,temp,ntemp)
       pz=u(1)*p(1,i)+u(2)*p(2,i)+u(3)*p(3,i)
       e=sqrt(p(1,i)**2+p(2,i)**2+p(3,i)**2+mass**2)
       pznew=gamma*(pz+beta*e)
       p(1,i)=p(1,i)+(pznew-pz)*u(1)
       p(2,i)=p(2,i)+(pznew-pz)*u(2)
       p(3,i)=p(3,i)+(pznew-pz)*u(3)
      enddo
c      
      return
      end

