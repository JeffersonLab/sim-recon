      REAL FUNCTION XSEC_JPSI(EB)
C---     J/psi photoproduction cross section (nb) near threshold: S.Brodsky et al PLB 498 (2001)
C---     The normalization is fit to the data at E>12 GeV (1 parameter)
C---     Calculates the cross_section * kineam_factor (phase space and tmin dependence)
C---     for 3 diagrams:
C         0) (s-s_thres)**0
C         1) (s-s_thres)**2
C         2) (s-s_thres)**4
C
      IMPLICIT NONE
      REAL EB
      REAL PHASE_V
C      VECTOR PKIN_PSI(15)
C      VECTOR CROSS_P(10,2)
C      
C
      INTEGER idiag
      REAL    e0,fkin,bsig,pow,s,s_thres,s0,am1,am2,am3,xx,xx0,scal
     +       ,tt,tslfac,xsec
C
      REAL xsecfac(3)   ! cross section factors for 3 diagrams (=cross section at E0=11 GeV)
C
      XSEC_JPSI=0.
      xsecfac(1)=0.    ! 3-gluon exchange (1-x)**0
      xsecfac(2)=0.36  ! 2-gluon exchange (1-x)**2
      xsecfac(3)=0.    ! 1-gluon exchange (1-x)**4

      e0=11.  ! Reference beam energy in GeV - this value was used in the fit (initially it was an arbitrary value) 
C               
      am1=0.938  ! target mass
      am2=3.1    ! J/psi mass
      am3=0.938  ! recoil mass
      tslfac=1.13  ! t-slope - the initial factor
C     
      s= am1**2+2.*am1*EB
      s0=am1**2+2.*am1*e0
      s_thres=(am2+am3)**2
      IF(s.LE.s_thres.OR.s0.LE.s_thres) GO TO 999
C
      xx= (s_thres-am1**2)/(s-am1**2)
      xx0=(s_thres-am1**2)/(s0-am1**2)
C
C-- Comment: Mark Srtikman defines x=1-(E_final-P_z_final)/M
C                            where E,P,M - parameters of the recoil
C
      xsec=0.
      DO idiag=1,3
         scal=xsecfac(idiag)
         IF(scal.GT.1.E-10) THEN
            pow=0.  
            tt=tslfac
            IF(idiag.EQ.2) THEN
               pow=2.
               tt=tslfac*9./4.
            ELSE IF(idiag.EQ.3) THEN 
               pow=4.
               tt=tslfac*9./1.
            ENDIF
            fkin=PHASE_V(EB,am1,am2,am3,tt)/PHASE_V(e0,am1,am2,am3,tt)
C
            bsig=1.
            IF(idiag.GT.1) THEN
               bsig=(1.-xx)**pow/(1.-xx0)**pow
            ENDIF
            bsig=bsig*(s-am1**2)**2/(s0-am1**2)**2 ! photon - ccbar coupling (from J.-M.Laget)
C
            xsec=xsec+fkin*bsig*scal
         ENDIF
      ENDDO
      XSEC_JPSI=xsec
C
 999  RETURN
      END
C
      REAL FUNCTION PHASE_V(EB,AM1,AM2,AM3,TT)
*       phase space
*  photoproduction: gamm+A-->Psi+A
*
      IMPLICIT NONE
      REAL EB,AM1,AM2,AM3,TT
      INTEGER ifl
C
      REAL s,ss,tmin,tmax,dsdt,b,pcm(2),ecm(2)
C
C      PHASE_V=1.
c      return
      s=am1**2+2.*am1*EB
C
      ss=SQRT(s)
      dsdt=1./16/3.1415/(s**2+am1**4-2.*s*am1**2)
      ecm(1)=(s-am1**2)/2./ss
      ecm(2)=(s+am2**2-am3**2)/2./ss
      pcm(1)=ecm(1)
      pcm(2)=SQRT(ecm(2)**2-am2**2)
      tmin=2.*(ecm(1)*ecm(2)-pcm(1)*pcm(2))-am2**2
      tmax=2.*(ecm(1)*ecm(2)+pcm(1)*pcm(2))-am2**2
      PHASE_V=dsdt/tt*(EXP(-tt*tmin)-EXP(-tt*tmax))
C
C      write(6,*) ss,dsdt,ecm,pcm,tmin,tmax,PHASE_V
      END
