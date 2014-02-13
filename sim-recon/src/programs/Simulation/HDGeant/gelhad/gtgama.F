c mon oct  3 16:35:06 pdt 1994
c version implementing hadronic interactions via |gelhad|.
c a. snyder
c
c Mon Oct 30 17:38:39 PST 1995
c modifications of my version of |gtgama| for use as 
c "stanadard" part of |bbsim|. Note: i've not started from
c the geant version 3.21 code, but i've looked carefully and
c the code i use seems identical apart from the |gelhad| 
c hooks.
c A. Snyder
c
c Thu Jan 18 11:23:06 PST 1996
c - Change from FFREAD based to dbin based control.
c   Any FFREAD values get overwritten by dbin if bbgeom_read is called.
c   Thus overwriting dbin data is now the default operation in bbsim.
c   Add verify step to print values (leave behind a check on common blocks).
c - Change to implicit none found some mis-spellings and a bug with epcut.
c   Art Snyder fixed the bugs (|epcut| and |one| not being defined) by
c   looking at geant321.car version.
c F. Kral
c
c Fri Jan 19 17:08:03 PST 1996
c - Change to process/mechanism name LMECGELH which is not taken by others.
c F. Kral
c
c Wed Jan 24 16:34:07 PST 1996
c - Change to protect against interaction in vaccuum
c A. Snyder
c
*cmz :          05/08/94  13.55.57  by  unknown
*-- author :a
      subroutine gtgama
c.
c.    ******************************************************************
c.    *                                                                *
c.    *   photon track. computes step size and propagates particle     *
c.    *    through step.                                               *
c.    *                                                                *
c.    *   ==>called by : gtrack                                        *
c.    *      authors    r.brun, f.bruyant l.urban ********             *
c.    *                                                                *
c.    ******************************************************************
c.
      implicit none
#include "gelhad/gelhadused.inc"
#include "gnbase/gelhad_db.inc"
#include "geant321/gcbank.inc"
#include "geant321/gccuts.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcphys.inc"
#include "geant321/gcstak.inc"
#include "geant321/gctmed.inc"
#include "geant321/gcmulo.inc" 
#include "geant321/gctrak.inc"
#if  _debug_
#include "geant321/gcunit.inc"
#endif
      real epsmac
#if  ! _singl_
c double precison
      parameter (epsmac=1.e-6)
      double precision one,xcoef1,xcoef2,xcoef3,zero
      parameter (one=1.0D0,zero=0.0D0)
#endif
#if  _singl_
c single precision
      parameter (epsmac=1.e-11)
      real *4 one,xcoef1,xcoef2,xcoef3,zero
      parameter (one=1.0,zero=0.0)
#endif
      real *4 epcut
      parameter (epcut=1.022e-3)

c
c hadronic interaction variables - defaults here apply when
c call to |gtgamaff| entry is not activivated. Note that
c since |jphadr| default is 0, it is not possible to enable
c |gelhad| unless |_GELH_| flag has been defined in |Flags.h|.
c
c Kral - do not set default values here since expected to crash on AIX.
c      - still can relatively safely assume that jphadr wakes up equal to 0.
      integer *4 jphadr    !hadronic interaction on/off switch
c                          !any non-zero turns it one; 
      real *4 ecutphadr        !min energy for hadronic interactions
      real *4 scale           !scale factor for xsections
      integer *4 mode_gphad            !mode for /gphad/
      real *4 ethresh_gphad         !ethresh for /gphad/
c      integer *4 jphadr/0/ !hadronic interaction on/off switch
cc                          !any non-zero turns it one; 
c      real *4 ecutphadr/0.500/ !min energy for hadronic interactions
c      real *4 scale/1.0/      !scale factor for xsections
c      integer *4 mode_gphad/4/         !mode for /gphad/
c      real *4 ethresh_gphad/0.150/  !ethresh for /gphad/
c      real *4 sigma      !inverse interaction lengths (cm**-1)
c      real *4 stephadr   !step to next hadronic interaction
      real *4 sigma      !inverse interaction lengths (cm**-1)
      real *4 stephadr   !step to next hadronic interaction
c
c     force order of parameters for /ffkey/
      integer *4 iparameters(5)
C       you cannot mix SAVE and COMMON for iparameters - rgj 96/1/28
C      save iparameters
      equivalence(iparameters(1),jphadr)
      equivalence(iparameters(2),ecutphadr)
      equivalence(iparameters(3),scale)
      equivalence(iparameters(4),mode_gphad)
      equivalence(iparameters(5),ethresh_gphad)
c
c     backwards compatibility -- ffread to dbin
      equivalence (iparameters(1),jphadr_gelhad)
c
      real *4 temp(10)   !work space
c
c     setup for command
      integer *4 idumm
      integer *4 ival
      real *4 fval
      equivalence(ival,fval)
      character *(*) action
c
c     declarations needed after switched to implicit none
c
      integer  iproc
      real     gekrt1
      integer  ist
      integer  jst
      integer  i
      real     vectmp
      integer  i1
      real     stopmx
c
c     first call
c
c      logical epcutfirst/.true./
c      save epcutfirst
c
c     for gtgamavrfydeb
c
      logical vrfyprint/.true./
      save vrfyprint
      logical vrfyinit/.false./
      save vrfyinit
c.
c.    ------------------------------------------------------------------
*
*
*
* *** particle below energy threshold ?  short circuit
*
*
      if (gekin.le.cutgam) goto 998
*
* *** update local pointers if medium has changed
*
      if(iupd.eq.0)then
         iupd  = 1
         jphot = lq(jma-6)
         jcomp = lq(jma-8)
         jpair = lq(jma-10)
         jpfis = lq(jma-12)
         jrayl = lq(jma-13)
      endif
*
* *** compute current step size
*
      iproc  = 103
      step   = stemax
      gekrt1 = 1 .-gekrat
*
*  **   step limitation due to pair production ?
*
*
c      if (epcutfirst) then
c        epcutfirst = .false.
c        write (6, *) 'GTGAMA:  epcut = ', epcut
c      endif
      if (getot.gt.epcut) then
         if (ipair.gt.0) then
            steppa = gekrt1*q(jpair+iekbin) +gekrat*q(jpair+iekbin+1)
            spair  = steppa*zintpa
            if (spair.lt.step) then
               step  = spair
               iproc = 6
            endif
         endif
      endif
*
*  **   step limitation due to compton scattering ?
*
      if (icomp.gt.0) then
         stepco = gekrt1*q(jcomp+iekbin) +gekrat*q(jcomp+iekbin+1)
         scomp  = stepco*zintco
         if (scomp.lt.step) then
            step  = scomp
            iproc = 7
         endif
      endif
*
*  **   step limitation due to photo-electric effect ?
*
      if (gekin.lt.0.4) then
         if (iphot.gt.0) then
            stepph = gekrt1*q(jphot+iekbin) +gekrat*q(jphot+iekbin+1)
            sphot  = stepph*zintph
            if (sphot.lt.step) then
               step  = sphot
               iproc = 8
            endif
         endif
      endif
*
*  **   step limitation due to photo-fission ?
*
      if (jpfis.gt.0) then
         steppf = gekrt1*q(jpfis+iekbin) +gekrat*q(jpfis+iekbin+1)
         spfis  = steppf*zintpf
         if (spfis.lt.step) then
            step  = spfis
            iproc = 23
         endif
      endif
*
*  **   step limitation due to hadronic interactions ?
*
*     THIS IS GELHAD STUFF
      if (jphadr.gt.0) then
       if(getot.ge.ecutphadr) then !enough energy to bother with ?
        call gpsig(sigma)
        if(sigma.gt.0.0) then !not in vaccuum ?
         sigma=sigma*scale
         call grndm(temp,1)
         stephadr=-alog(temp(1))/sigma
         if(stephadr.lt.step) then
          step=stephadr
          iproc=lmecgelh
         endif !stephadr.lt.step
        endif !sigma.gt.0.0
       endif !getot.gt.ecutphadr
      endif !jphadr.gt.0
*
*  **   step limitation due to rayleigh scattering ?
*
      if (irayl.gt.0) then
         if (gekin.lt.0.01) then
            stepra = gekrt1*q(jrayl+iekbin) +gekrat*q(jrayl+iekbin+1)
            srayl  = stepra*zintra
            if (srayl.lt.step) then
               step  = srayl
               iproc = 25
            endif
         endif
      endif
*
      if (step.lt.0.) step = 0.
*
*  **   step limitation due to geometry ?
*
      if (step.ge.safety) then
         call gtnext
         if (ignext.ne.0) then
            step   = snext + prec
            inwvol= 2
            iproc = 0
            nmec =  1
            lmec(1)=1
         endif
*
*        update safety in stack companions, if any
         if (iq(jstak+3).ne.0) then
            do 10 ist = iq(jstak+3),iq(jstak+1)
               jst = jstak +3 +(ist-1)*nwstak
               q(jst+11) = safety
   10       continue
            iq(jstak+3) = 0
         endif
*
      else
         iq(jstak+3) = 0
      endif
*
* *** linear transport
*
      if (inwvol.eq.2) then
         do 20 i = 1,3
            vectmp  = vect(i) +step*vect(i+3)
            if(vectmp.eq.vect(i)) then
*
* *** correct for machine precision
*
               if(vect(i+3).ne.0.) then
                  vectmp = vect(i)+abs(vect(i))*sign(1.,vect(i+3))*
     +            epsmac
                  if(nmec.gt.0) then
                     if(lmec(nmec).eq.104) nmec=nmec-1
                  endif
                  nmec=nmec+1
                  lmec(nmec)=104
#if  _debug_
                  write(chmail, 10000)
                  call gmail(0,0)
                  write(chmail, 10100) gekin, numed, step, snext
                  call gmail(0,0)
10000 format(' boundary correction in GTGAMA: ',
     +       '    GEKIN      NUMED       STEP      SNEXT')
10100 format(31x,e10.3,1x,i10,1x,e10.3,1x,e10.3,1x)
#endif
               endif
            endif
            vect(i) = vectmp
   20    continue
      else
         do 30 i = 1,3
            vect(i)  = vect(i) +step*vect(i+3)
   30    continue
      endif
*
      sleng = sleng +step
*
* *** update time of flight
*
      tofg = tofg +step/clight
*
* *** update interaction probabilities
*
      if (getot.gt.epcut) then
         if (ipair.gt.0) zintpa = zintpa -step/steppa
      endif
      if (icomp.gt.0)    zintco = zintco -step/stepco
      if (gekin.lt.0.4) then
         if (iphot.gt.0) zintph = zintph -step/stepph
      endif
      if (jpfis.gt.0)    zintpf = zintpf -step/steppf
      if (irayl.gt.0) then
         if (gekin.lt.0.01) zintra = zintra -step/stepra
      endif
*
      if (iproc.eq.0) go to 999
      nmec = 1
      lmec(1) = iproc
*
*  ** pair production ?
*
      if (iproc.eq.6) then
         call gpairg
*
*  ** compton scattering ?
*
      else if (iproc.eq.7) then
         call gcomp
*
*  ** photo-electric effect ?
*
      else if (iproc.eq.8) then
*       calculate range of the photoelectron ( with kin. energy ephot)
*
        if(gekin.le.0.001)  then
            jcoef = lq(jma-17)
         if(gekrat.lt.0.7) then
            i1 = max(iekbin-1,1)
         else
            i1 = min(iekbin,nekbin-1)
         endif
         i1 = 3*(i1-1)+1
         xcoef1 = q(jcoef+i1)
         xcoef2 = q(jcoef+i1+1)
         xcoef3 = q(jcoef+i1+2)
         if(xcoef1.ne.0.) then
            stopmx = -xcoef2+sign(one,xcoef1)*sqrt(xcoef2**2 - (xcoef3-
     +      gekin/xcoef1))
         else
            stopmx = - (xcoef3-gekin)/xcoef2
         endif
*
*        do not call gphot if this (overestimated) range is smaller
*                than safety
*
         if (stopmx.le.safety) goto 998
        endif
 
         call gphot
*
*  ** rayleigh effect ?
*
      else if (iproc.eq.25) then
         call grayl
*
*  ** photo-fission ?
*
      else if (iproc.eq.23) then
         call gpfis
*
*  ** electro-production ?
*
*     AGAIN THIS IS GELHAD STUFF
      else if(iproc.eq.lmecgelh) then
       call gphad(mode_gphad,ethresh_gphad)
*
      endif
*
         goto 999
998      destep = gekin
         gekin  = 0.
         getot  = 0.
         vect(7)= 0.
         istop  = 2
         nmec   = 1
         lmec(1)= 30
999   continue
      return
c
      entry gtgamaff !read parameters from cards
c
      call ffkey('GELH',iparameters,5,'MIXED')
c
c set defaults
      jphadr=0
      ecutphadr=0.500           ! original version had this mis-spelled
      scale=1.0
      mode_gphad=4
      ethresh_gphad=0.150
c
c also start things off with use flag set false.
      gelhadused=.false.
c
      return
c
      entry gtgamac(action,idumm) !control programs behavior
c
      if(len(action).le.4) return
      if(action(1:3).eq.'put') then
       ival=idumm
      endif !action(1:3).eq.'put'
      if(action.eq.'get:jphadr') then
       ival=jphadr
      else if(action.eq.'put:jphadr') then
       jphadr=ival
      else if(action.eq.'get:ecutphadr') then
       fval=ecutphadr
      else if(action.eq.'put:ecutphadr') then
        ecutphadr=fval          ! original version had this mis-spelled
      else if(action.eq.'get:scale') then
       fval=scale
      else if(action.eq.'put:scale') then
       scale=fval
      else if(action.eq.'get:mode_gphad') then
       ival=mode_gphad
      else if(action.eq.'put:mode_gphad') then
       mode_gphad=ival
      else if(action.eq.'get:ethresh_gphad') then
       fval=ethresh_gphad
      else if(action.eq.'put:ethresh_gphad') then
       ethresh_gphad=fval
      endif !action.eq.'get:jphadr'
      if(action(1:3).eq.'get') then
       idumm=ival
      endif !action(1:3).eq.'get'


c
      return
c
      entry gtgamavrfydeb ! verify parameters from database/ffread cards
c
c  Kral 1/18/95 - Use this when the db common block is NOT used (old way).
c
      if (.not. vrfyprint) return
      if (.not. vrfyinit) then
        vrfyinit = .true.
        if (jphadr .ne. 0) then
          write (6, *)
          write (6, *) 'GELHAD parameter verification in gtgamavrfydeb'
          write (6, *) '  jphadr  = ', jphadr
          write (6, *) '  ecut    = ', ecutphadr
          write (6, *) '  scale   = ', scale
          write (6, *) '  mode    = ', mode_gphad
          write (6, *) '  ethresh = ', ethresh_gphad
          write (6, *)
        endif
        
      endif
c
      return
      end
