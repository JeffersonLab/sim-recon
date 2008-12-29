! Change decay status of certain particles in Lund
! For now, all we have is pi and K decay ...

	subroutine pythia6_setdecay(decaypik,decaylamks)

!	implicit none

!#include "py6dat3.inc"
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      INTEGER MDCY,MDME,KFDP
      DOUBLE PRECISION BRAT
      SAVE/PYDAT3/

	integer pycomp
	logical decaypik, decaylamks
	external pycomp

! Perhaps have Lund decay pi+/-, K+/-, and K0_long
! 	KF=211 --> charged pi
! 	KF=321 --> charged K
! 	KF=130 --> K_L0

	if (decayPiK) then
	  mdcy(pycomp(211),1)=1
	  mdcy(pycomp(321),1)=1
	  mdcy(pycomp(130),1)=1
	else
	  mdcy(pycomp(211),1)=0
	  mdcy(pycomp(321),1)=0
	  mdcy(pycomp(130),1)=0
	endif

! Perhaps have Lund decay Lambda (,Sigma,Xi,Omega) and K0_short
! 	KF=310  --> K_S0
! 	KF=3222 --> Sigma+ and anti-Sigma-
! 	KF=3112 --> Sigma- and anti-Sigma+
! 	KF=3322 --> Xi0 and anti-Xi0
! 	KF=3312 --> Xi- and anti-Xi+
! 	KF=3122 --> Lambda0 and anti-Lambda0
! 	KF=3334 --> Omega- and anti-Omega+

	if (decayLamKs) then
	  mdcy(pycomp(310),1)=1
	  mdcy(pycomp(3222),1)=1
	  mdcy(pycomp(3112),1)=1
	  mdcy(pycomp(3322),1)=1
	  mdcy(pycomp(3312),1)=1
	  mdcy(pycomp(3122),1)=1
	  mdcy(pycomp(3334),1)=1
	else
	  mdcy(pycomp(310),1)=0
	  mdcy(pycomp(3222),1)=0
	  mdcy(pycomp(3112),1)=0
	  mdcy(pycomp(3322),1)=0
	  mdcy(pycomp(3312),1)=0
	  mdcy(pycomp(3122),1)=0
	  mdcy(pycomp(3334),1)=0
	endif

	return
	end
