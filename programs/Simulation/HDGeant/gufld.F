*CMZ :  3.21/02 29/03/94  15.41.25  by  S.Giani
*-- Author :
      SUBROUTINE GUFLD(X,F)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       User routine to compute the magnetic field F             *
C.    *       at space point X                                         *
C.    *                                                                *
C.    *    ==>Called by : GHELIX,GRKUTA                                *
C.    *                                                                *
C.    ******************************************************************
C.
      DIMENSION X(3),F(3)
C.
C.    ------------------------------------------------------------------
C.
      real Tesla2kG, cm2inch, zorigin
      parameter (Tesla2kG=10, cm2inch=0.3937, zorigin=26.0)

      real r(3),B(3)
      real rsampl(20000,3),Bsampl(20000,3)
      integer nsampl,nrows,ncolumns
      common /fieldMap/rsampl,Bsampl,nsampl,ncolumns,nrows
      data Bsampl/20000*0, 20000*0, 20000*0/
      data nsampl,ncolumns,nrows/0,0,0/
      save /fieldMap/
      if (nsampl.eq.0) then
        open(unit=78,file='dsolenoid.table',status='old',err=9)
        rlast = 0
        do 5 i=1,20000
          read(78,*,err=5,end=6) r,B
          if (B(3).ne.0) then
            nsampl = nsampl + 1
            rsampl(nsampl,1) = r(3)
            rsampl(nsampl,2) = r(1)
            Bsampl(nsampl,1) = B(1)
            Bsampl(nsampl,2) = B(2)
            Bsampl(nsampl,3) = B(3)
          endif
          if (r(1).ne.rlast) then
            nrows = nrows + 1
            rlast = r(1)
          endif
    5   enddo
    6   continue
        nrows = nrows + 1
        ncolumns = nsampl/nrows
        do i=1,nrows
          rsampl(ncolumns+i,1) = rsampl(i*ncolumns,2)
        enddo
        close(unit=78)
      endif

      r(1) = zorigin + X(3)*cm2inch
      r(2) = sqrt(X(1)**2 + X(2)**2)*cm2inch
      phi = atan2(X(2),X(1))
      F(3) = FINT(2,r,ncolumns,rsampl,Bsampl(1,3))*Tesla2kG
      F(1) = FINT(2,r,ncolumns,rsampl,Bsampl(1,1))*Tesla2kG
      F(2) = F(1)*sin(phi)
      F(1) = F(1)*cos(phi)
      return

    9 stop 'ERROR - Could not open input field map file dsolenoid.table'
      end
