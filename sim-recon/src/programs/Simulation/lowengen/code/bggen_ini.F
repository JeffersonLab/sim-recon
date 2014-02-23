C
      SUBROUTINE BGGEN_INI(ID,E0,EP,ELIM,ZCOLLIM)
C
C---   Photoproduction by the coherent Brem. beam
C---   ID - histogram with the dN/dE per 1uA beam
C      E0 - e- energy
C      EP - coherent peak energy
C      ELIM - energy limits
C      ZCOLLIM - distance to the collimator
C
      IMPLICIT NONE
      INTEGER ID
      REAL E0,EP,ELIM(2),ZCOLLIM
C
      REAL DNIDX,DNCDX,GPXSECT,GPXSECP
      REAL SAIDE,SAIDXSECA
      LOGICAL HEXIST
C
      INTEGER i,nb,ipro,npro,id1,idt,ncth,icth,ihi
      REAL emn,emx,flx1,flx2,flx,xsec,dx,de,e,xstot,dcth,cth
     +    ,xlum,ecurr,xsth,targth,frate
C
      CHARACTER tit*132,cpro*2,cenr*4
C
C---   Initialize the coherent Bremsstrahlung 
C
      CALL COBREMS(E0,EP,ZCOLLIM/100.)  ! collimator distance in m
C
      targth=30.  ! target thickness
      ecurr=2.25E-6  ! beam current on a 20um crystal - the "high luminosity" setting
      xlum=ecurr/1.6E-19*targth*0.071*0.602 ! luminosity factor:  2.25uA on 20um (1.7e-4 RL) crystal, 30 cm LH2 (1/b) 
C
      emn=ELIM(1)
      emx=ELIM(2)
      de=0.005
      dx=de/E0
      nb=INT((emx-emn)/de+0.001)
C
      CALL HBOOK1(ID  ,'Beam flux dn/dE',nb,emn,emx,0.)
      CALL HBOOK1(ID+1,'Beam incoh flux dn/dE',nb,emn,emx,0.)
      CALL HBOOK1(ID+2,'Beam coher flux dn/dE',nb,emn,emx,0.)
      CALL HBOOK1(ID+5,'Beam incoh flux dn/dx ',nb,emn/E0,1.,0.)
      CALL HBOOK1(ID+9,'Total cross section in mb',nb,emn,emx,0.)
      CALL HBOOK1(ID+10,'Beam flux dn/dE',nb,emn,emx,0.)
C
      frate=0.
      DO i=1,nb
         e=ELIM(1)+(i-0.5)*de
         flx1=DNIDX(e/E0)
         flx2=DNCDX(e/E0)
         CALL HF1(ID+5,e/E0,flx1)
         flx1=flx1*dx/de
         flx2=flx2*dx/de
         flx=flx1+flx2
         CALL HF1(ID+1,e,flx1)
         CALL HF1(ID+2,e,flx2)
         CALL HF1(ID,e,flx)
         xstot=GPXSECT(e) ! x-sec in mb
         CALL HF1(ID+9,e,xstot)
         xstot=xstot*1.E-3 ! x-sec in b
         CALL HF1(ID+10,e,xstot*flx*xlum)
         frate=frate+xstot*flx*xlum*de
      ENDDO
      WRITE(6,FMT='(//10X,''Rates:'')')
      WRITE(6,1000) ecurr*1.E6,emn,emx
 1000 FORMAT(10X,'Beam: ',F8.2,' uA e-, gamma in ',2F6.2,' GeV') 
      WRITE(6,1005) targth
 1005 FORMAT(10X,'Target: ',F8.2,' cm LH2')
      WRITE(6,1010) frate/1000.
 1010 FORMAT(10X,'Interaction rate: ',F8.1,' kHz') 
C
      CALL HCOPY(ID+10,ID+11,'    ')  ! a copy of the final histogram to be used for HRNDM1
C
C---      Initialize the processes
C
      npro=10     ! number of defined processes 
      ncth=100    ! number of bins in the COS(th) distribution
      dcth=2./ncth
      DO ipro=1,npro
         id1=ID+10000*ipro
         WRITE(cpro,FMT='(I2)') ipro
         CALL HBOOK1(id1,'X-section for process '//cpro,nb,emn,emx,0.)
C
C---      Define the COS(TH) plots?
C
         IF(ipro.LE.2.OR.     ! SAID
     +      ipro.EQ.8) THEN   ! eta
            CALL HBOOK1(id1+1,'refer for COS(TH) for process '//cpro
     +                 ,nb,emn,emx,0.)
         ENDIF
C
         DO i=1,nb
            e=ELIM(1)+(i-0.5)*de
C
C---            Full x-section
C
            IF(ipro.LE.2) THEN
C
C---          SAID is used
C
               xsec=SAIDE(e,ipro,1)  ! SAID cross section, supressed above 2 GeV
               ihi=i
            ELSE
               xsec=GPXSECP(e,ipro) ! x-sec in mb
               ihi=0
C               CALL GPXCOSTH(e,ipro,0.,ihi,xsth)  ! check the COS(th) distribution
            ENDIF
C
            CALL HF1(id1,e,xsec)
            CALL HF1(id1+1,e,REAL(ihi))
C
C--- Get the cos(th) distributions
C
            idt=id1+10+ihi
            IF(ihi.NE.0.AND..NOT.HEXIST(idt)) THEN   ! fill the COS(th) distrib if not yet filled
               WRITE(cenr,FMT='(I4)') i
               CALL HBOOK1(idt,'COS(TH), proc '//cpro//' energy '//cenr
     +                     ,ncth,-1.,1.,0.)
               DO icth=1,ncth
                  cth=-1.+(icth-0.5)*dcth
                  IF(ipro.LE.2) THEN
                     xsth=SAIDXSECA(e,cth,ipro,1)
                  ELSE
C                     CALL GPXCOSTH(e,ipro,cth,ihi,xsth) ! get the COS(th) distribution
                  ENDIF
                  CALL HF1(idt,cth,xsth)
               ENDDO
            ENDIF
         ENDDO
C
         IF(ipro.EQ.1) THEN
            CALL HCOPY(id1,ID+15,'X-section for all process ')
         ELSE
            CALL HOPERA(id1,'+',ID+15,ID+15,1.,1.)
         ENDIF
C
      ENDDO
C
 999  RETURN
      END
