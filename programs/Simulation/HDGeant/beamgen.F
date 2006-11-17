      subroutine beamgen(t0)
      real t0
*
* Generates a single beam photon according to the coherent bremsstrahlung
* model in cobrems.F using beam energy and primary coherent edge energies
* specified by the user.  The photon begins its lifetime just upstream of
* the primary collimator (WARNING: position hard-wired in the code below)
* and is tracked by the simulation from there forward.  The time of the
* photon is generated so that it will reach the midplane of the GlueX
* target at TOF=t0 if it makes it that far.
*
#include "geant321/gcunit.inc"
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcscan.inc"
#include "geant321/gcomis.inc"
#include "cobrems.inc"
 
      real vertex(3),plab(3)
      integer nvert,nt
      real rndm(20)
 
      integer coSwitch
      real freqMaximum
      common /coherentGen/coSwitch,freqMaximum
      data coSwitch,freqMaximum/0,0/
      save /coherentGen/
      real xMinimum,Theta02,CollimPos
      parameter (xMinimum=1e-2)
      parameter (Theta02=1.8)
      parameter (CollimPos=-2200.0)
      integer nubuf
      real ubuf(10)

      call GRNDM(rndm,7)
      phim = rndm(1)*TWOPI
      rhom = mospread*sqrt(-2*log(rndm(2)))
      thxMosaic = rhom*cos(phim)
      thyMosaic = rhom*sin(phim)
      phib = rndm(3)*TWOPI
      varBeam = (emit/spot)**2
      rhob = sqrt(-2*varBeam*log(rndm(4)))
      thxBeam = rhob*cos(phib)
      thyBeam = rhob*sin(phib)
      phis = rndm(5)*TWOPI
      varMS = sigma2MS(t*rndm(6))
      rhos = sqrt(-2*varMS*log(rndm(7)))
      thxMS = rhos*cos(phis)
      thyMS = rhos*sin(phis)
      cos45 = 1/sqrt(2d0)
      rotate(1,1) = 0
      rotate(1,2) = cos45			!point (1,0,0) along beam
      rotate(1,3) = -cos45			!point (0,1,1) vertically
      rotate(2,1) = 0
      rotate(2,2) = cos45
      rotate(2,3) = cos45
      rotate(3,1) = 1
      rotate(3,2) = 0
      rotate(3,3) = 0
      call rotmat(rotate,thxBeam+thxMS-thx-thxMosaic,0d0,0d0)
      call rotmat(rotate,0d0,thyBeam+thyMS-thy-thyMosaic,0d0)
      do i=1,10000
         call GRNDM(rndm,4)
         phi = rndm(1)*TWOPI
         x = xMinimum**rndm(2)
         theta2 = Theta02*rndm(3)/(1-rndm(3)+1e-30)
         coSwitch = coSwitch+1
         if (coSwitch/2*2.eq.coSwitch) then    !try coherent generation
            freq = dNcdxdp(x,phi)
            f = freq*rndm(3)
            do ip=1,q2points
               if (f.le.q2weight(ip)) then
                  theta2 = q2theta2(ip)
                  goto 5
               endif
            enddo
 5          continue
            ppol = polarization(x,theta2)
         else                                  !try incoherent generation
            freq = dNidxdt2(x,theta2)
            freq = freq*(Theta02+theta2)**2/(TWOPI*Theta02)
            ppol = 0
         endif
         freq = x*freq
         if (freq.gt.freqMaximum) then
            freqMaximum = freq*1.2
         elseif (freq.ge.freqMaximum*rndm(4)) then
            goto 50
         endif
      enddo
      write(5,*) 'Gukine: photon beam generator failed, giving up'

50    continue
      theta = sqrt(theta2)*(me/E)
      thetaX = thxMS+theta*cos(phi)
      thetaY = thyMS+theta*sin(phi)
      plab(3) = E*x
      plab(1) = plab(3)*thetaX
      plab(2) = plab(3)*thetaY
      call GRNDM(rndm,2)
      phic = rndm(1)*TWOPI
      rhoc = spot*sqrt(-2*log(rndm(2)))
      vertex(1) = (rhoc*cos(phic) + D*thetaX)*100
      vertex(2) = (rhoc*sin(phic) + D*thetaY)*100
      vertex(3) = CollimPos
      ubuf(1) = ppol
      nubuf = 1
#if defined BEAM_BOX_SIZE
      call GRNDM(rndm,2)
      ubuf(2) = rndm(1)*BEAM_BOX_SIZE
      ubuf(3) = rndm(2)*BEAM_BOX_SIZE
      vertex(1) = vertex(1) + ubuf(2)
      vertex(2) = vertex(2) + ubuf(3)
      nubuf = 3
#endif
      call settofg(vertex,t0)
      call GSVERT(vertex,0,0,ubuf,nubuf,nvert)
      call GSKINE(plab,1,nvert,0,0,nt)
      end

      subroutine settofg(vertex,time0)
      real vertex(3)
      real time0
#include "geant321/gconst.inc"
#include "geant321/gctrak.inc"
      real target_center
      parameter (target_center=65.0)
      TOFG=time0*1e-9+(vertex(3)-target_center)/CLIGHT
      end
