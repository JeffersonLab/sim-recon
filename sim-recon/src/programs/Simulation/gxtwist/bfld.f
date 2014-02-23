      real function bfld(i,x,y,z)
      real x,y,z
      integer i
      real r(3),F(3)
      r(1) = x
      r(2) = y
      r(3) = z
      call gmedia(r,numed)
      call gufld(r,F)
      bfld = F(i)
      end
