      subroutine sigmag(a,z,egamma,sigma)
********************************************************
* Routine to compute gN xsections
* Bothers nobody else
* Author: A. Snyder
* version 1.00 - Wed Nov 22 17:57:32 PST 1995
* version 1.01 - Wed Jan 24 16:34:07 PST 1996
* modified to give 0 xsection for A<1 (vaccuum)
*******************************************************
      implicit none
*
* input:
      real *4 a !atomic number (number of n+p)
      real *4 z !number of ps
      real *4 egamma !photon energy (gev)
*
* output:
      real *4 sigma !cross-section in mb
*
* internal:
      real *4 l !levinger factor
*
      sigma=0.0
      if(a.lt.0.99) return !no material
      sigma=0.2 !a fixed but fairly meaning less value
      if(egamma.lt.0.050) return
      if(egamma.lt.0.2) then !quasi-deutron region 
         call sigmatf(egamma,sigma) !thorlacius&fearing
         l=7
         if(a.lt.4) l=a
         sigma=sigma*l*(a-z)*z/a
         return
      else
         call sigmapdg(egamma,sigma)
         sigma=0.5*sigma*a
      endif                     !(egamma.lt.0.2)
*
      return
      end
*
      subroutine sigmatf(e,sigma) !thorlacius&fearing parameterization
      real *4 e,sigma
      real *4 fourpi/12.566371/
      real *4 c1/261.0/
      real *4 c2/-110.0/
      real *4 c3/24.6/
      real *4 c4/-17.1/
      real *4 c5/5.76/
      real *4 c6/-2.05/
      real *4 c7/0.267/
      real *4 c8/113.0/
      sigma=c1*exp(c2*e)+c3*exp(c4*e)+(c5+c6*e)/(1.0+c8*(e-c7)**2)
      sigma=sigma*fourpi/1000.0 !4pi and ->mb
      return
      end
*
      subroutine sigmapdg(eg,sigma) !D2 xsections from pdg
      real *4 eg,sigma
      real *4 e(14) /.20,.238,.258,.316,.341,.419,.543,.570,.706,.940
     >,.975,1.11,1.22,17.5/
      real *4 s(14)/0.24,0.40,0.82,0.90,0.90,0.56,0.35,0.4,0.5,0.34
     >,0.35,0.33,0.29,0.2/
      real *4 le(14) !log table energy
      integer *4 i,up,lo
      logical init/.false./

      save init,le
*
      if(.not.init) then
         do 10 i=1,14
            le(i)=alog(e(i))
 10      continue
         init=.true.
      endif                     !.not.init
*
* handle off table
      sigma=0.0
      if(eg.lt.e(1)) return
      sigma=s(14)
      if(eg.gt.e(14)) return
*
* look up position in table
      do 1000 i=2,14
         if(eg.gt.e(i)) go to 1000
         up=i
         goto 1099
 1000 continue
 1099 continue
*
* interpolate between values in table
      lo=up-1
      sigma=s(lo)+(s(up)-s(lo))*(alog(eg)-le(lo))/(le(up)-le(lo))
*
      return
      end



