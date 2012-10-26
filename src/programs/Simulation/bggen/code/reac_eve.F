      SUBROUTINE REAC_EVE(IERR)
C
C---  Simulates 1 event - a single reaction
C     Reaction: gamma+p --> A+B , A - recoil, B - meson 
C               ISIMUL=1 A=p (14), B=J/psi (83)
C
      IMPLICIT NONE
      INTEGER IERR
C
      INCLUDE 'bg_ctrl.inc'
      INCLUDE 'bg_partc.inc'
      INCLUDE 'bg_evec.inc'
      INCLUDE 'bg_reac.inc'
C
      REAL RNDM,GBRWIGN
C
      INTEGER i,j,ip,np,ityp,ntry,ires
     +       ,ityd,ntry1,ihel
      REAL ebeam,ecm,ecm2,bet(4),qq,ct,st,phi
     +    ,twopi
     +    ,amtot    ! sum of the masses 
     +    ,ppf,epf1,ppi,tt,tmn,tmx,amdec,amd(6),pcms(4)
     +    ,wdm
C 
      INTEGER mxoutl
      PARAMETER (mxoutl=6)
      REAL ami(2),pcmi(4,2),plabi(4,2)
     +     ,am(mxoutl),pcm(4,mxoutl),plab(4,mxoutl)
      INTEGER ity(mxoutl),ndec(mxoutl),kdec(3,mxoutl),kdectyp(mxoutl)
     +       ,it1dec(mxoutl),itorig(mxoutl)
C
C     ------------------------------------------------------------------
C
      IERR=1
      NTRA=0
      IF(ISIMUL.LT.1.OR.ISIMUL.GT.1) GO TO 999
C
C---    Beam energy
C
      ebeam=PIN(3,1)
C
C---      Initial state (beam goes along Z - no rotation applied)
C
      DO i=1,2
         ami(i)=AMIN(i)
         DO j=1,3
            plabi(j,i)=PIN(j,i)
         ENDDO
         qq=plabi(1,i)**2+plabi(2,i)**2+plabi(3,i)**2
         plabi(4,i)=SQRT(qq+ami(i)**2)
      ENDDO
      DO j=1,4
         pcms(j)=plabi(j,1)+plabi(j,2)
      ENDDO
C
C      write(6,*) 'ami', ami,plab(4,1),plab(4,2)
      ecm2=ami(1)**2+ami(2)**2+2.*plabi(4,1)*plabi(4,2)
      ecm=SQRT(ecm2)
C
      ntry=0
 30   np=0
      ntry=ntry+1
      amtot=0.
      ires=0
      DO ip=1,2
         ityp=IPREAC(ip)
         IF(ityp.GT.0.AND.ityp.LE.MXPART) THEN
            np=np+1
            ity(np)=ityp
            am(np)=AM_PART(ityp)
            amdec=0.
            ndec(np)=0
            itorig(np)=0
            it1dec(np)=0
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
C      write(6,*) ' np..', np,amtot,ecm-0.01
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
      IF(np.EQ.2) THEN         
C
C---      In CM: momentum and energies of the particles
C
         epf1=(ecm2+am(1)**2-am(2)**2)/2./ecm
         ppf =SQRT(epf1**2-am(1)**2)        ! final momentum
         ppi=SQRT(pcmi(4,2)**2-ami(2)**2)  ! initial momentum
         IF(ppf.LE.0.) GO TO 999
C
         qq=ami(2)**2+am(1)**2-2.*epf1*pcmi(4,2)
         tmx=qq+2.*ppf*ppi
         tmn=qq-2.*ppf*ppi
         IF(TSLREAC.LT.0.001) THEN
            tt=tmn+(tmx-tmn)*RNDM(qq)
         ELSE
            tt=1./TSLREAC*ALOG(EXP(TSLREAC*tmn)
     +         +RNDM(qq)*(EXP(TSLREAC*tmx)-EXP(TSLREAC*tmn)))
         ENDIF
         ct=(tt-qq)/2./ppf/ppi
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
               it1dec(i)=np+1
               DO j=1,ndec(i)
                  amd(j)=AM_PART(kdec(j,i))
                  am (np+j)=amd(j)
                  ity(np+j)=kdec(j,i)
                  ndec(np+j)=0
                  itorig(np+j)=i
                  it1dec(np+j)=0
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
      ENDIF
C
      DO i=1,np
         DO j=1,3
            PTRA(j,i)=plab(j,i)
         ENDDO
         AMTRA(i)=am(i)
         ITPTRA(1,i)=ity(i)
         DO j=2,6
            ITPTRA(j,i)=0
         ENDDO
C         write(6,*) i,ity(i),MXPGEANT,IPLUND(ity(i)),itorig(i),it1dec(i)
         IF(ity(i).GT.0.AND.ity(i).LE.MXPGEANT) THEN
            ITPTRA(3,i)=IPLUND(ity(i))
         ENDIF
         ITPTRA(4,i)=itorig(i)
         ITPTRA(5,i)=it1dec(i)
         IF(it1dec(i).GT.0) ITPTRA(6,i)=it1dec(i)+ndec(i)-1
         ITPTRA(2,i)=1
         IF(it1dec(i).NE.0) ITPTRA(2,i)=10   ! indicates that this particle should not be used in GEANT 
      ENDDO
      NTRA=np
C
      IERR=0
 999  CONTINUE
C      write(6,*) ebeam,IEVPROC,ibin,xstot,xssum,NTRA
C
      END
C

