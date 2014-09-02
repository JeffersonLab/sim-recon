      subroutine beamgen(t0)
      real t0                   ! beam bucket, ns
*
* Generates a single beam photon according to the coherent bremsstrahlung
* model in cobrems.F using beam energy and primary coherent edge energies
* specified by the user.  The photon begins its lifetime just upstream of
* the primary collimator (WARNING: position hard-wired in the code below)
* and is tracked by the simulation from there forward, but its time t0 
* identifies its beam bucket, ie. the time the photon would reach the time-
* reference plane at the middle of the target.  Beam bucket t0=0 is the
* one that generates the event and defines the time origin for all hit times.
*
* To enable beam motion spreading, define the beam box size below (cm)
* #define BEAM_BOX_SIZE 5

#include "geant321/gcunit.inc"
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcscan.inc"
#include "geant321/gcomis.inc"
#include "geant321/gctrak.inc"
#include "cobrems.inc"
 
      real vertex(4),plab(5),pbeam
      real rhom,phim
      real rhop,phip
      real rhoc,phic
      integer nvert,nt
      real rndm(20)
 
c  freqMaximum = probability density cutoff for coherent/incoherent 
c  bremsstrahlung generator, defined on the measure [dx dphi dy] where
c    x = E_gamma/E_end_point
c    phi = azimuthal angle (radians)
c    y = a normalized polar angle parameter defined by the relation
c    dy = theta0^2 dtheta^2 / (theta0^2 + theta^2)^2 with 0<=y<=1.
c  The probability is for a single electron, so the scale is that of
c  the target thickness (radiation lengths) divided by 2pi.  A good
c  choice for freqMaximum is the target thickness in radiation lengths.
c  A warning is printed in the simulation output log each time a value
c  freq > freqMaximum is generated; a few ppm of these is no problem.
      real xMinimum,freqMaximum,beamStartZ,Theta02
      common /coherentGen/xMinimum,freqMaximum,beamStartZ,Theta02
      data xMinimum/0.01/
      data freqMaximum/9.0e-4/
      data beamStartZ/-2200.0/
      data Theta02/1.8/
      save /coherentGen/

      integer nProfileBins
      parameter (nProfileBins=500)
      real freqProfile(nProfileBins)
      real freqIntegral(nProfileBins)
      common /freqTables/freqProfile,freqIntegral
      data freqProfile/nProfileBins*0/
      data freqIntegral/nProfileBins*0/
      save /freqTables/
      real Wincoh
      parameter (Wincoh=0.1)

      integer nubuf
      real ubuf(10)

      logical hexist
      external hexist
      common /genstate/ppol,rndm
      save /genstate/

      call GRNDM(rndm,7)
      phim = rndm(1)*TWOPI
      rhom = mospread*sqrt(-2*log(rndm(2)))
      thxMosaic = rhom*cos(phim)
      thyMosaic = rhom*sin(phim)
      phib = rndm(3)*TWOPI
      rhob = sqrt(-2*log(rndm(4)))
      thxBeam = (emitx/spot)*rhob*cos(phib)
      thyBeam = (emity/spot)*rhob*sin(phib)
      phis = rndm(5)*TWOPI
      varMS = sigma2MS(t*rndm(6))
      rhos = sqrt(-2*varMS*log(rndm(7)))
      thxMS = rhos*cos(phis)
      thyMS = rhos*sin(phis)
      cos45 = 1/sqrt(2d0)
      rotate(1,1) = 0
      rotate(1,2) = cos45   !point (1,0,0) along beam
      rotate(1,3) = -cos45   !point (0,1,1) vertically
      rotate(2,1) = 0
      rotate(2,2) = cos45
      rotate(2,3) = cos45
      rotate(3,1) = 1
      rotate(3,2) = 0
      rotate(3,3) = 0
      call rotmat(rotate,thxBeam+thxMS-thx-thxMosaic,0d0,0d0)
      call rotmat(rotate,0d0,thyBeam+thyMS-thy-thyMosaic,0d0)
      if (freqIntegral(nProfileBins).eq.0) then
         x1 = 1
         x0 = xMinimum**(1./nProfileBins)
         freqProfile(1) = dNcdxdp((x0+x1)/2,TWOPI/4)
         freqIntegral(1) = freqProfile(1)*(x1-x0)
         do ip=2,nProfileBins
            x1 = x0
            x0 = xMinimum**(ip*1./nProfileBins)
            freqProfile(ip) = dNcdxdp((x0+x1)/2,TWOPI/4)
            freqIntegral(ip) = freqIntegral(ip-1)
     +                        +freqProfile(ip)*(x1-x0)
         enddo
      endif
      do i=1,1000000000
         call GRNDM(rndm,5)
         if (rndm(1).gt.1/(Wincoh+1)) then             !try coherent generation
            f = freqIntegral(nProfileBins)*rndm(2)
            do ip=1,nProfileBins
               if (f.le.freqIntegral(ip)) then
                  x0 = xMinimum**((ip-1.)/nProfileBins)
                  x1 = xMinimum**((ip*1.)/nProfileBins)
                  if (ip.gt.1) then
                     f0 = freqIntegral(ip-1)
                  else
                     f0 = 0
                  endif
                  f1 = freqIntegral(ip)
                  fp = freqProfile(ip)
                  x = (x0*(f1-f)+x1*(f-f0))/(f1-f0)
                  go to 4
               endif
            enddo
 4          continue
            phi = rndm(3)*TWOPI
            freq = dNcdxdp(x,phi)
            f = freq*rndm(4)
            do iq=1,q2points
               if (f.le.q2weight(iq)) then
                  theta2 = q2theta2(iq)
                  goto 5
               endif
            enddo
 5          continue
            freq = freq*freqIntegral(nProfileBins)/fp
            freq = freq*TWOPI
            ppol = polarization(x,theta2)
         else                                  !try incoherent generation
            x = xMinimum**rndm(2)
            phi = rndm(3)*TWOPI
            theta2 = Theta02*rndm(4)/(1-rndm(4)+1e-30)
            freq = dNidxdt2(x,theta2)
            freq = freq*(Theta02+theta2)**2/Theta02
            freq = freq*x*(-log(xMinimum))
            freq = freq*Wincoh
            ppol = 0
         endif
         if (freq.gt.freqMaximum) then
            print *, 'Warning from beamgen: freq=',freq,
     +               ' is greater than freqMaximum=',freqMaximum
         endif
         if (freq.ge.freqMaximum*rndm(5)) then
            goto 50
         endif
      enddo
      print *, 'Error in beamgen:',
     +         ' photon beam generator failed, giving up!'
      stop

50    continue

#if DEBUG_CB_BEAM_GENERATOR
      print *, 'success after',i,' attempts'
      if (.not.hexist(20)) then
        call hbnt(20,'coherent generator state',' ')
        call hbnt(21,'incoherent generator state',' ')
        call hbname(20,'genstate',ppol,'ppol:r')
        call hbname(20,'genstate',rndm(1),'varndm(5):r')
        call hbname(21,'genstate',ppol,'ppol:r')
        call hbname(21,'genstate',rndm(1),'varndm(5):r')
      endif
      if (ppol.eq.0) then
        call hfnt(21)
      else
        call hfnt(20)
      endif
#endif

      call GRNDM(rndm,2)
      phip = rndm(1)*TWOPI
      rhop = sqrt(-2*log(rndm(2)))
      pbeam = E+Erms*rhop*cos(phip)
      theta = sqrt(theta2)*(me/E)
      thetaX = thxBeam+thxMS+theta*cos(phi)
      thetaY = thyBeam+thyMS+theta*sin(phi)
      plab(5) = pbeam*x
      plab(1) = plab(5)*thetaX
      plab(2) = plab(5)*thetaY
      plab(3) = sqrt(plab(5)**2-plab(1)**2-plab(2)**2)
      plab(4) = plab(5)
      call GRNDM(rndm,2)
      phic = rndm(1)*TWOPI
      rhoc = spot*sqrt(-2*log(rndm(2)))
      vertex(1) = (rhoc*cos(phic)-D*thxBeam+D*thetaX)*100
      vertex(2) = (rhoc*sin(phic)-D*thyBeam+D*thetaY)*100
      vertex(3) = beamStartZ
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
      call GSKINE(plab,1,nvert,0,0,nt) ! push the beam photon on the stack
      vertex(4) = TOFG
      call hitTagger(vertex,vertex,plab,plab,0.,1,0,0)
      end
