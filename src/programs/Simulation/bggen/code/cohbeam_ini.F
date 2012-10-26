C
      SUBROUTINE COHBEAM_INI(ID,E0,EP,ELIM,ZCOLLIM)
C
C---   Photoproduction by the coherent Brem. beam
C---   ID - histogram with the dN/dE*sigma(E), 
C              where dN/dE - coh. Brem., sigma(E) - total photoprod. on protons 
C      E0 - e- energy
C      EP - coherent peak energy (<0.01 - incoherent only)
C      ELIM - energy limits
C      ZCOLLIM - distance to the collimator
C
      IMPLICIT NONE
      INTEGER ID
      REAL E0,EP,ELIM(2),ZCOLLIM
C
      INCLUDE 'bg_ctrl.inc'
C
      REAL DNIDX,DNCDX,GPXSECT,GETXSEC
C
      INTEGER i,nb,ibrem
      REAL emn,emx,flx,dx,de,e,xstot,xlum,ecurr,targth,frate,vv
     +    ,cohfac,epk_tmp,fflux
C
      ibrem=1
      IF(ELIM(2).LT.ELIM(1)*1.006) ibrem=0
C
C---   Initialize the coherent Bremsstrahlung 
C
      targth=30.                ! target thickness
      ecurr=2.25E-6             ! beam current on a 20um crystal - the "high luminosity" setting
      xlum=ecurr/1.6E-19*targth*0.071*0.602 ! luminosity factor:  2.25uA on 20um (1.7e-4 RL) crystal, 30 cm LH2 (1/b) 
C
      epk_tmp=EP
      cohfac=1.
      IF(epk_tmp.LT.0.01) THEN  ! only use the incoherent radiation (a 15% discrepancy in the flux may occur)
         epk_tmp=0.75*E0
         cohfac=0.
      ENDIF
C
      IF(ibrem.NE.0) THEN
         CALL COBREMS(E0,epk_tmp,ZCOLLIM/100.) ! collimator distance in m
      ENDIF
C
      emn=ELIM(1)
      emx=ELIM(2)
      de=MIN(0.005,emx-emn)
      dx=de/E0
      nb=INT((emx-emn)/de+0.001)
      NHBEA=nb
C
      CALL HBOOK1(ID   ,'Beam flux dn/dE*sigma(E)' ,nb,emn,emx,0.)
      CALL HBOOK1(ID+10,'Total cross section in mb',nb,emn,emx,0.)
      CALL HBOOK1(ID+11,'Beam flux dn/dE'          ,nb,emn,emx,0.)
C
      fflux=0.
      frate=0.
      DO i=1,nb
         e=ELIM(1)+(i-0.5)*de
         IF(ISIMUL.EQ.0) THEN
            xstot=GPXSECT(e)*1.E-3 ! x-sec in b
         ELSE IF(ISIMUL.EQ.1) THEN
            xstot=GETXSEC(e)       ! x-sec in b
C            write(6,*) e,xstot
         ENDIF
         IF(ibrem.NE.0) THEN
            flx=(DNIDX(e/E0)+DNCDX(e/E0)*cohfac)*dx/de
            fflux=fflux+flx*de
            vv=xstot*flx*xlum
            frate=frate+vv*de
            CALL HF1(ID+11,e,flx)
         ELSE
            vv=1.
         ENDIF
         CALL HF1(ID   ,e,vv)
         CALL HF1(ID+10,e,xstot)
      ENDDO
C      CALL HPRINT(ID)
C
      RATESEC=0.
      IF(ibrem.NE.0) THEN
         RATESEC=frate
         WRITE(6,FMT='(//10X,''Rates:'')')
         WRITE(6,1000) ecurr*1.E6,emn,emx,fflux*ecurr/1.6e-19
 1000    FORMAT(10X,'Beam: ',F8.2,' uA e-, gamma in ',2F6.2,' GeV '
     +         ,'  flux ',E11.3,' Hz') 
         WRITE(6,1005) targth
 1005    FORMAT(10X,'Target: ',F8.2,' cm LH2')
         WRITE(6,1010) frate
 1010    FORMAT(10X,'Event rate : ',E11.3,' Hz'//)
      ENDIF 
C
      CALL HCOPY(ID,ID+1,'    ')  ! a copy of the final histogram to be used for HRNDM1
C
      RETURN
      END
