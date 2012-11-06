      real function ymax()
      include 'nt.inc'
      vector ym(100)
c     theta = atan2(-pin(1),pin(3))*180/3.1416-9.944
      theta = -xin(3)
      ibin=pin(4)/12*100+1
      if (theta.gt.ym(ibin)) then
        ym(ibin)=theta
      endif
      ymax=ym(ibin)
      end
