      SUBROUTINE BGGEN_EVE(ID)
C
C---  Simulates low energy (<3 GeV) photoproduction by a coherent Bremsstrahlung beam
C     Reaction: gamma+p
C     ID - is the starting number of a set of predefined histograms with E,cos(th) distributions
C     
C     Processes:
C               1) p pi0
C               2) n pi+
C               3) p pi+ pi- non res
C               4) p rho0
C               5) Delta++ pi-
C               6) p pi0 pi0
C               7) n pi+ pi0
C               8) p eta
C               9) p pi+ pi- pi0
C              10) n pi+ pi+ pi-
C
      IMPLICIT NONE
      INTEGER ID
C
      INCLUDE 'bg_partc.inc'
      INCLUDE 'bg_evec.inc'
C
      REAL HRNDM1,RNDM,HI,GBRWIGN
      LOGICAL HEXIST
C
      INTEGER i,j,ip,np,npm,ibin,nproc,iproc,ityp,ok,ihi,ierr,ntry,ires
     +       ,id1,ifla,ityd,ntry1,ihel,np1
      REAL ebeam,xstot,xssum,xstmp,rnd,ecm,ecm2,bet(4),qq,ct,st,phi,wgt
     +    ,twopi
     +    ,amtot    ! sum of the masses 
     +    ,pcmm(4)  ! 4-mom of the mesons
     +    ,betm(4)  ! vel of CM as seen from the rest frame of the mesons
     +    ,ppf,epf1,epf2,ppi,tt,tmn,tmx,amdec,amd(6),xfac,pcms(4),par(6)
     +    ,wdm
C
      INTEGER mxpr,mxout
      PARAMETER (mxpr=10,mxout=6)
C 
      REAL ami(2),pcmi(4,2),plabi(4,2)
     +     ,am(mxout),pcm(4,mxout),plab(4,mxout)
     +     ,wgt4mx(mxpr)   ! max weight for the 4-body process (potentially, for each process)  
      INTEGER ity(mxout),ndec(mxout),kdec(3,mxout),kdectyp(mxout)
C
      INTEGER itypr(mxout,mxpr)
      DATA itypr/14, 7, 0, 0, 0, 0
     +          ,13, 8, 0, 0, 0, 0
     +          ,14, 8, 9, 0, 0, 0
     +          ,14,80, 0, 0, 0, 0
     +          ,82, 9, 0, 0, 0, 0
     +          ,14, 7, 7, 0, 0, 0
     +          ,13, 8, 7, 0, 0, 0
     +          ,14,17, 0, 0, 0, 0
     +          ,14, 8, 9, 7, 0, 0
     +          ,13, 8, 8, 9, 0, 0
     + /
      DATA wgt4mx/10*-1./
C
C     ------------------------------------------------------------------
C
      IEVPROC=0
      nproc=10  ! number of defined processes
C
C---    Beam energy
C
      ebeam=HRNDM1(ID+11)
      CALL HXI(ID+11,ebeam,ibin)  ! get ibin  - the bin number for this energy
C
      NTRA=2
      ITYPTR(1)=1
      AMTRA(1)=AM_PART(ITYPTR(1))
      PTRA(1,1)=0.
      PTRA(2,1)=0.
      PTRA(3,1)=ebeam
      ITYPTR(2)=14
      AMTRA(2)=AM_PART(ITYPTR(2))
      PTRA(1,2)=0.
      PTRA(2,2)=0.
      PTRA(3,2)=0.
      DO i=1,2
         NDECTR(i)=0
      ENDDO
C
C---      Initial state (beam goes along Z - no rotation applied)
C
      DO i=1,2
         ami(i)=AMTRA(i)
         DO j=1,3
            plabi(j,i)=PTRA(j,i)
         ENDDO
         qq=plabi(1,i)**2+plabi(2,i)**2+plabi(3,i)**2
         plabi(4,i)=SQRT(qq+ami(i)**2)
      ENDDO
      DO j=1,4
         pcms(j)=plabi(j,1)+plabi(j,2)
      ENDDO
C
      ecm2=ami(1)**2+ami(2)**2+2.*plabi(4,1)*plabi(4,2)
      ecm=SQRT(ecm2)
C
C---    Choose a process
C
      xstot=HI(ID+9,ibin)
      IF(xstot.LE.0.) GO TO 999  ! no simulation (low energy?) 
      xssum=HI(ID+15,ibin)       ! sum of all processes   
      IF(xssum.LE.0.) GO TO 999  
C
      xstmp=0.
      rnd=RNDM(xstmp)
      iproc=1
      DO i=1,nproc-1
         xstmp=xstmp+HI(ID+10000*iproc,ibin)/xssum
         IF(xstmp.GE.rnd) GO TO 20
         iproc=i+1
      ENDDO
 20   CONTINUE
C
      IEVPROC=iproc
C
      ntry=0
 30   np=0
      ntry=ntry+1
      amtot=0.
      ires=0
      DO ip=1,mxout
         ityp=itypr(ip,IEVPROC)
         IF(ityp.GT.0.AND.ityp.LE.MXPART) THEN
            np=np+1
            ity(np)=ityp
            am(np)=AM_PART(ityp)
            amdec=0.
            ndec(np)=0
            DO i=1,3
               ityd=KD_PART(i,ityp)
               IF(ityd.GT.0.AND.ityd.LE.MXPART) THEN
                  ndec(np)=ndec(np)+1
                  kdec(ndec(np),np)=ityd
                  amdec=amdec+AM_PART(ityd)
               ENDIF
            ENDDO
            kdectyp(np)=KD_PART(4,ityp)
            IF(WG_PART(ityp).GT.0.) THEN
               ires=1
               ntry1=0
 35            ntry1=ntry1+1
               wdm=WG_PART(ityp)*GBRWIGN(am)
C               write(6,*) am(np),wdm,amdec
               IF(am(np)+wdm.LT.amdec+0.01) THEN
                  IF(ntry1.LT.1000) GO TO 35
                  WRITE(6,*) ' *** BGGEN_EVE unsuff mass for decay '
     +            ,ityp,am(np),wdm,am(np)+wdm,amdec 
                  GO TO 999
               ENDIF
               am(np)=am(np)+wdm
            ENDIF
            amtot=amtot+am(np)
         ENDIF
      ENDDO
      IF(np.LT.1) GO TO 999
      IF(amtot.GE.ecm-0.01) THEN
         IF(ntry.LT.1000) GO TO 30
         GO TO 999
      ENDIF
C
      DO i=1,3
         bet(i)=(plabi(i,1)+plabi(i,2))/(plabi(4,1)+plabi(4,2))
      ENDDO
      bet(4)=(plabi(4,1)+plabi(4,2))/ecm
      DO i=1,2
         CALL GLOREN(bet,plabi(1,i),pcmi(1,i))
      ENDDO
      DO i=1,3
         bet(i)=-bet(i)
      ENDDO
C
C---  Treat the kinematics as 2-body one, in CM
C
      twopi=ACOS(0.)*4.
      ierr=1
      IF(np.EQ.2) THEN         
C        IF(IEVPROC.LE.2.OR.  ! SAID
C     +     IEVPROC.EQ.6      ! eta
C     +    ) THEN 
C---      In CM: momentum and energies of the particles
C
         epf1=(ecm2+am(1)**2-am(2)**2)/2./ecm
         ppf =SQRT(epf1**2-am(1)**2)        ! final momentum
         ppi=SQRT(pcmi(4,2)**2-ami(2)**2)  ! initial momentum
         IF(ppf.LE.0.) GO TO 999
C
         id1=ID+10000*IEVPROC
         ihi=0
         IF(HEXIST(id1+1)) THEN
            ihi=INT(HI(ID+10000*IEVPROC+1,ibin)+0.1)
            IF(ihi.GT.0) THEN 
               ct= HRNDM1(ID+10000*IEVPROC+10+ihi)
               ct=-ct           ! first particle is the recoil - invert the COS(TH)
               ierr=0
            ENDIF
         ENDIF
         IF(ierr.NE.0) THEN
            qq=ami(2)**2+am(1)**2-2.*epf1*pcmi(4,2)
            tmn=-(qq+2.*ppf*ppi)
            tmx=-(qq-2.*ppf*ppi)
            CALL GPXCOSTHR(IEVPROC,ebeam,tmn,tmx,ct,ierr) ! generated for the secondary baryon
         ENDIF
         IF(ierr.NE.0) THEN
            WRITE(6,*) ' *** Error in simulating COS(TH) ',ierr
     +                    ,' proc=',IEVPROC
         ENDIF
C
         st=SQRT(1.-ct**2)
         phi=twopi*RNDM(st)
C
C---    2-body
C
         pcm(4,1)=epf1
C
         pcm(1,1)=ppf*st*COS(phi)
         pcm(2,1)=ppf*st*SIN(phi)
         pcm(3,1)=ppf*ct
C
         DO i=1,3
            pcm(i,2)=-pcm(i,1)
         ENDDO
         pcm(4,2)=ecm-pcm(4,1)
C
C---       Boost to Lab
C
         DO i=1,2
            CALL GLOREN(bet,pcm(1,i),plab(1,i))
         ENDDO
C
C---       Decays?
C
         DO i=1,2
            IF(ndec(i).GT.0) THEN
               DO j=1,ndec(i)
                  amd(j)=AM_PART(kdec(j,i))
                  am (np+j)=amd(j)
                  ity(np+j)=kdec(j,i)
                  ndec(np+j)=0
               ENDDO
               IF(ndec(i).EQ.2) THEN  ! 2-body decay
                  ihel=kdectyp(i)  ! decay angle flag =0 - unoform, =1 - rho-like, =2 - j/psi-like
                  CALL OMDECA2(plab(1,i),amd(1),ihel,plab(1,np+1))
               ELSE IF(ndec(i).EQ.3) THEN
                  CALL OMDECA3(plab(1,i),amd(1),0.,plab(1,np+1))
               ENDIF
               np=np+ndec(i)
            ENDIF
         ENDDO
C
      ELSE IF(np.EQ.3) THEN         
C
         xfac=0.
         CALL OMDECA3(pcms(1),am(1),xfac,plab(1,1))
C
      ELSE IF(np.EQ.4) THEN         
C
C---      Phase space: 
C
         IF(wgt4mx(IEVPROC).LT.0.) THEN  ! initialize the max weight
            DO i=1,20000
               wgt=0.
               CALL GDECAN(np,ecm,am,wgt,pcm(1,1))
               wgt4mx(IEVPROC)=MAX(wgt4mx(IEVPROC),wgt)
            ENDDO
            wgt4mx(IEVPROC)=wgt4mx(IEVPROC)*1.2
         ENDIF
         wgt=wgt4mx(IEVPROC)
         CALL GDECAN(np,ecm,am,wgt,pcm(1,1))
         DO i=1,np
            CALL GLOREN(bet,pcm(1,i),plab(1,i))
         ENDDO
C
      ENDIF
C
      DO i=1,np
         DO j=1,3
            PTRA(j,i+2)=plab(j,i)
         ENDDO
         AMTRA(i+2)=am(i)
         ITYPTR(i+2)=ity(i)
         NDECTR(i+2)=ndec(i)
      ENDDO
      NTRA=2+np
C
 999  CONTINUE
C      write(6,*) ebeam,IEVPROC,ibin,xstot,xssum,NTRA
C
      END
C
C
      REAL FUNCTION RNDM(X)
      IMPLICIT NONE
      REAL X
      REAL a
C
      CALL RANMAR(a,1)
      RNDM=a
      RETURN
      END

