c Tue Oct 11 11:32:47 PDT 1994
c user special version of /gheish/ ... /gheishp/ to
c avoid interfering with particle types
c
c Mon Oct  3 16:35:06 PDT 1994
c /geant321/ version of /gelhad/. A. Snyder
c
c Mon Dec 11 15:52:37 PST 1995
c add flag to indicate |gphad| has been invoked

c Fri Jan 19 16:37:11 PST 1996
c Add count of interactions per event
c Set the name of the mechanism (only have a list of 30 to choose between)
c F. Kral
c
c Mon Mar  3 10:07:47 PST 1997
c Add momentum conserving models 5, 6 and 7
c A. Snyder
c
c Mon Mar 17 10:57:46 PST 1997
c Add vector dominance model 8
c A. Snyder
c
c Wed Apr  9 10:27:54 PDT 1997 
c add quasi-deutron model 9
c A. Snyder
c
c Fri Mar 12 15:31:37 PST 1999
c Fix A=0 problem by not letting mixtures that contain hydrogen interact
c

c
      subroutine gphad(mode,ethresh)
      implicit none
      save
c
c /geant321/ commons
c
#include "geant321/gcmate.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gcbank.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcunit.inc"
#include "gelhad/gelhadused.inc"
c
c     declarations when switch to implicit none
c
      real ethresh
      integer i
      logical ok
      logical goth
C
c
c /gelhad/ variables
      integer *4 mode      !/gelhad/ model
      integer *4 fail      !failure code
      logical interact/.true./
      integer *4 alias     !gamma's alias
      real *4 temp(10)     !work space
      real *4 mpi          !pi mass when needed
      real *4 ppipl(3)     !pi plus 3-mom
      real *4 mompl        !momenutum
      real *4 kinpl        !kinetic energy
      integer *4 ipiplus   !geant code
      real *4 ppimn(3)     !pi minus 3-mom
      real *4 mommn        !mometum
      real *4 kinmn        !kinetic energy 
      integer *4 ipimins   !geant code
      real *4 piplint/0.5/ !pi plus interaction prob.
      real *4 pimnint/0.5/ !pi minus interaction prob.
      integer *4 iprot     !proton code
      integer *4 ineut     !neutron code
      real *4 mn           !neutron mass
      real *4 mp           !proton mass
      real *4 pn(3)        !mom neutron from quasi deutron decay
      real *4 pp(3)        !mom proton from quasi deutron decay
      real *4 ed           !quasi-deuteron energy
      real *4 pd           !quasi-deuteron momentum
      real *4 momn         !mom neut
      real *4 momp         !mom prot
      real *4 kinn         !kinetic energy neut
      real *4 kinp         !kinetic energy prot
      real *4 pnint/0.5/   !neutron interaction prob.
      real *4 ppint/0.5/   !proton interaction prob.
      real *4 kinergy      !kinetic energy
      real *4 test         !accept/reject test number
      integer *4 ntemp     !words used in work space
      character *20 name
      integer *4 type
      real *4 mass,qchrg,lifetime !particle properties
      real *4 mtarget !target mass
      real *4 mrecoil !recoil mass
      real *4 precoil !lab momentum of recoil nucleon
c
      integer *4 mprod
      parameter (mprod=100)
      integer *4 nprod
      integer *4 iprod(mprod)
      real *4 tprod(mprod)
      real *4 pprod(3,mprod)
c
      integer *4 ncall/0/
c
      ncall=ncall+1
c
c set mechanism
      KCASE = NAMEC(12)         ! Hadronic mechanism (generic name)
c
      fail=1301
      gelhadused=.true.
      ngelhperev=ngelhperev + 1
      nprod=0
      if(.not.interact) return !off ?
      call gelrndmset(1) !set to use |geant| random numbers
c
      if(mode.eq.0) go to 1
      if(mode.eq.2) go to 2
      if(mode.eq.4) then
       if(getot.lt.ethresh) go to 1
       test=(1.0-ethresh/getot)**0.5
       call grndm(temp,1)
       if(temp(1).gt.test) go to 1
       go to 2
      endif !mode.eq.4
      if(mode.eq.5) go to 5
      if(mode.eq.6) go to 6
      if(mode.eq.7) then
         if(getot.lt.ethresh) go to 5
         test=(1.0-ethresh/getot)**0.5
         call grndm(temp,1)
         if(temp(1).gt.test) go to 5
         go to 6
      endif
      if(mode.eq.8) go to 8
      if(mode.eq.9) go to 9
      if(mode.eq.10) then
         if(getot.lt.ethresh) go to 9
         test=(1.0-ethresh/getot)**0.5
         call grndm(temp,1)
         if(temp(1).gt.test) go to 9
         go to 8
      endif
      go to 1313
c
1     continue !g->N model
c
c select proton or neutron
      call grndm(temp(1),1)
      alias=13 !neutron
      if(temp(1).gt.(a-z)/a) alias=14 !proton
      ntemp=10
      call gfpart(alias,name,type,mass,qchrg,lifetime,temp,ntemp)
c
c conserve energy and charge and all /gheish/
      kinergy=getot
      getot=getot+mass
      temp(7)=sqrt(getot**2-mass**2)
      call gelpcalc(temp(7),vect(4),temp)
      call gamate(1.0,1.0,goth)
      if(goth) then
*        call gelfill(alias,vect,mprod,nprod,iprod,tprod,pprod)
         call gelfill(alias,temp,mprod,nprod,iprod,tprod,pprod)
         go to 99
      endif
c
      call grmate(1.0,qchrg)       !remove nuc
      call gpgheip(alias,temp(7),kinergy)
      call gheishp(alias,temp,mprod,nprod,iprod,tprod,pprod)
      call grmate(-1.0,-qchrg)    !put nuc back
c
c
      go to 99
c
2     continue !g->pi model
c
c pick pi+ or pi0
      call grndm(temp(1),1)
      alias=8 !pi+
      if(temp(1).gt.z/a) alias=9 !pi-
      ntemp=10
      call gfpart(alias,name,type,mass,qchrg,lifetime,temp,ntemp)
c
c
c conserve energy and charge and call /gheish/
      temp(1)=getot-mass
      if(temp(1).lt.0.0) return !below threshold ?
      kinergy=temp(1)
      temp(7)=sqrt(getot**2-mass**2)
      call caspimset('SET:NFLFORCE',1)
      call caspipset('SET:NFLFORCE',2)
      temp(7)=sqrt(getot**2-mass**2)
      call gelpcalc(temp(7),vect(4),temp)
      call grmate(0.0,qchrg)
      call gpgheip(alias,temp(7),kinergy)
      call gheishp(alias,temp,mprod,nprod,iprod,tprod,pprod)
      call grmate(0.0,-qchrg)
      call caspimset('SET:NFLFORCE',0)
      call caspipset('SET:NFLFORCE',0)
c
      go to 99
c
 5    continue !g->N model with energy and momentum conservation
c
c select proton or neutron
      call grndm(temp(1),1)
      alias=13 !neutron
      if(temp(1).gt.(a-z)/a) alias=14 !proton
      ntemp=10
      call gfpart(alias,name,type,mass,qchrg,lifetime,temp,ntemp)
c
c conserve energy and charge and all /gheish/
      call gmmate(mtarget)

      call gamate(1.0,1.0,goth)      
      if(goth) then
*        call gelfill(alias,vect,mprod,nprod,iprod,tprod,pprod)
         call gelpcalc(vect(7),vect(4),temp)
         call gelfill(alias,temp,mprod,nprod,iprod,tprod,pprod)
         go to 99
      endif

      call grmate(1.0,qchrg)       !remove nuc
      call gmmate(mrecoil)
      call recoilframe(getot,mass,mtarget,mrecoil,temp(7),precoil,ok)
      if(.not.ok) go to 99
      kinergy=sqrt(temp(7)**2+mass**2)-mass
      call gelpcalc(temp(7),vect(4),temp)
      call gpgheip(alias,temp(7),kinergy)
      call gheishp(alias,temp,mprod,nprod,iprod,tprod,pprod)
      call labframe(vect(4),mrecoil,precoil,nprod,pprod,iprod)
      call grmate(-1.0,-qchrg)    !put nuc back
c
      go to 99
c
c
      
c
 6    continue !g->pi model with energy and momentum conservation
c
c pick pi+ or pi0
      call grndm(temp(1),1)
      alias=8 !pi+
      if(temp(1).gt.z/a) alias=9 !pi-
      ntemp=10
      call gfpart(alias,name,type,mass,qchrg,lifetime,temp,ntemp)
c
c
c conserve energy and charge and call /gheish/
      temp(1)=getot-mass
      if(temp(1).lt.0.0) return !below threshold ?
      call caspimset('SET:NFLFORCE',1)
      call caspipset('SET:NFLFORCE',2)
      call gmmate(mtarget)
      call grmate(0.0,qchrg)
      call gmmate(mrecoil)
      call recoilframe(getot,mass,mtarget,mrecoil,temp(7),precoil,ok)
      if(.not.ok) go to 99
      call gelpcalc(temp(7),vect(4),temp)
      kinergy=sqrt(temp(7)**2+mass**2)-mass
      call gpgheip(alias,temp(7),kinergy)
      call gheishp(alias,temp,mprod,nprod,iprod,tprod,pprod)
      call labframe(vect(4),mrecoil,precoil,nprod,pprod,iprod)
      call grmate(0.0,-qchrg)
      call caspimset('SET:NFLFORCE',0)
      call caspipset('SET:NFLFORCE',0)
c
      go to 99
c
 8    continue !vector dominance
c
c g->rho .. vector dominance model
      ipiplus=8
      ipimins=9
      call gfpart(ipiplus,name,type,mpi,qchrg,lifetime,temp,ntemp)
      temp(1)=getot-2.0*mpi
      if(temp(1).lt.0.0) return !below threshold ?
      call gelrhom(getot,mpi,mass)
      call gmmate(mtarget)
      mrecoil=mtarget
      call recoilframe(getot,mass,mtarget,mrecoil,temp(7),precoil,ok)
      if(.not.ok) go to 99
      call gelpcalc(temp(7),vect(4),temp)
      call gelrhodk(temp,mass,mpi,ppipl,ppimn)
      call grndm(temp(1),3)
      if(temp(1).lt.0.5) then
         if(temp(1).lt.piplint) then !pi+ interacts
            mompl=sqrt(ppipl(1)**2+ppipl(2)**2+ppipl(3)**2)
            kinpl=sqrt(mompl**2+mpi**2)-mpi
            call gpgheip(ipiplus,mompl,kinpl)
            call gheishp(ipiplus,ppipl,mprod,nprod,iprod,tprod,pprod)
         else
            call gelfill(ipiplus,ppipl,mprod,nprod,iprod,tprod,pprod)
         endif                  !(temp(1).lt.piplint)
         call gelfill(ipimins,ppimn,mprod,nprod,iprod,tprod,pprod)
      else
         if(temp(2).lt.pimnint) then !pi- interacts
            mommn=sqrt(ppimn(1)**2+ppimn(2)**2+ppimn(3)**2)
            kinmn=sqrt(mommn**2+mpi**2)-mpi
            i=nprod+1
            ntemp=0
            call gpgheip(ipimins,mommn,kinmn)
            call gheishp
     .      (ipimins,ppimn,mprod,ntemp,iprod(i),tprod(i),pprod(1,i))
            nprod=nprod+ntemp
         else
            call gelfill(ipimins,ppimn,mprod,nprod,iprod,tprod,pprod)
         endif                  !(temp(2).lt.pimnint)
         call gelfill(ipiplus,ppipl,mprod,nprod,iprod,tprod,pprod)
      endif                     !temp(3).lt.0.5
      call labframe(vect(4),mrecoil,precoil,nprod,pprod,iprod)
c
      go to 99
c
 9    continue !quasi-deutron model; g->D
c
      iprot=14
      ineut=13
      ntemp=10
      call gfpart(ineut,name,type,mn,qchrg,lifetime,temp,ntemp)
      call gfpart(iprot,name,type,mp,qchrg,lifetime,temp,ntemp)
      ed=mp+mn+getot
      pd=getot
      mass=sqrt(ed**2-pd**2)
c
      call gmmate(mtarget)

      call gamate(1.0,2.0,goth)      
      if(goth) then
*        call gelfill(iprot,vect,mprod,nprod,iprod,tprod,pprod)
*        call gelfill(ineut,vect,mprod,nprod,iprod,tprod,pprod)
         call gelpcalc(vect(7),vect(4),temp)
         call gelfill(iprot,temp,mprod,nprod,iprod,tprod,pprod)
         call gelfill(ineut,temp,mprod,nprod,iprod,tprod,pprod)
         go to 99
      endif

      call grmate(2.0,1.0)      !remove a deutron from target material
      call gmmate(mrecoil)
c
      call gelpcalc(pd,vect(4),temp)
      call geldeutdk(temp,mass,mn,mp,pn,pp)
c
      call grndm(temp,3)
      if(temp(3).lt.0.5) then
         if(temp(1).lt.pnint) then !neutron interacts
            momn=sqrt(pn(1)**2+pn(2)**2+pn(3)**2)
            kinn=sqrt(momn**2+mn**2)-mn
            call gpgheip(ineut,momn,kinn)
            call gheishp
     .           (ineut,pn,mprod,nprod,iprod,tprod,pprod)
         else
c            call call gelfill(ineut,pn,mprod,nprod,iprod,tprod,pprod)
            call gelfill(ineut,pn,mprod,nprod,iprod,tprod,pprod)
         endif
c         call call gelfill(iprot,pp,mprod,nprod,iprod,tprod,pprod)
         call gelfill(iprot,pp,mprod,nprod,iprod,tprod,pprod)
      else
         if(temp(2).lt.ppint) then !proton interacts
            momp=sqrt(pp(1)**2+pp(2)**2+pp(3)**2)
            kinp=sqrt(momp**2+mp**2)-mp
            i=nprod+1
            ntemp=0
            call gheishp
     .           (iprot,pp,mprod,ntemp,iprod(i),tprod(i),pprod(1,i))
            nprod=nprod+ntemp
         else
            call gelfill(iprot,pp,mprod,nprod,iprod,tprod,pprod)
         endif
         call gelfill(ineut,pn,mprod,nprod,iprod,tprod,pprod)
      endif                     !temp(3).lt.0.5
      call grmate(-2.0,-1.0) !put target material back to normal
c
      go to 99
c
1313  continue
       write(chmail,*) 'fail=',fail,'in /gelhad/'
      return
c
99    continue
c
c copy to /gckine/ stack
      do 1000 i=1,nprod
       ngkine=ngkine+1
       gkin(1,ngkine)=pprod(1,i)
       gkin(2,ngkine)=pprod(2,i)
       gkin(3,ngkine)=pprod(3,i)
       gkin(5,ngkine)=iprod(i)
       tofd(ngkine)=tofg+tprod(i)
       ntemp=10
       call gfpart(iprod(i),name,type,mass,qchrg,lifetime,temp,ntemp)
       temp(1)=mass**2+pprod(1,i)**2+pprod(2,i)**2+pprod(3,i)**2
       gkin(4,ngkine)=sqrt(temp(1))
       gpos(1,ngkine)=vect(1)
       gpos(2,ngkine)=vect(2)
       gpos(3,ngkine)=vect(3)
1000  continue
c
      istop=1 !end the gamma
      return
c
      end
c
      subroutine gelfill(ipart,pmom,mp,np,id,t,p)
      implicit none
      integer *4 ipart
      real *4 pmom(3)
      integer *4 mp
      integer *4 np
      integer *4 id(mp)
      real *4 t(mp)
      real  *4 p(3,mp)
c
      if(np.lt.mp) then
         np=np+1
         id(np)=ipart
         p(1,np)=pmom(1)
         p(2,np)=pmom(2)
         p(3,np)=pmom(3)
         t(np)=0.0
      else
c
         write(6,*) 
     .   'gelfill: |prod| array full; particle list truncated'
c         
      endif
c
      return
      end
c
      subroutine gelpcalc(mom,u,p)
      implicit none
      real *4 mom     !momentum
      real *4 u(3)    !unit vector
      real *4 p(3)    !momentum vector
      p(1)=u(1)*mom
      p(2)=u(2)*mom
      p(3)=u(3)*mom
      return
      end
c
      subroutine geldeutdk(pd,md,mn,mp,pn,pp)
      implicit none
      real *4 md,mn,mp
      real *4 pd(3),pn(3),pp(3)
      real *4 kd(5),kn(5),kp(5)
c
c ang distribution parms
      real *4 a/1.0/ 
      real *4 b/-1.0/
c     
      kd(1)=pd(1)
      kd(2)=pd(2)
      kd(3)=pd(3)
      kd(4)=sqrt(kd(1)**2+kd(2)**2+kd(3)**2+md**2)
      kd(5)=md
      kn(5)=mn
      kp(5)=mp
      call geltwobdo('set:a',a)
      call geltwobdo('set:b',b)
      call geltwob(kd,kn,kp)
      pn(1)=kn(1)
      pn(2)=kn(2)
      pn(3)=kn(3)
      pp(1)=kp(1)
      pp(2)=kp(2)
      pp(3)=kp(3)
c
      return
      end
c
      subroutine gelrhodk(prho,mrho,mpi,ppl,pmn)
      implicit none
      real *4 prho(3) !rho mom vec
      real *4 mrho    !rho mass
      real *4 mpi     !pi mass
      real *4 ppl(3)  !pi+ mom vec 
      real *4 pmn(3)  !pi- mom vec
c
c ang distribution parms
      real *4 a/1.0/
      real *4 b/-1.0/
      real *4 krho(5),kpipl(5),kpimn(5)
      krho(1)=prho(1)
      krho(2)=prho(2)
      krho(3)=prho(3)
      krho(4)=sqrt(krho(1)**2+krho(2)**2+krho(3)**2+mrho**2)
      krho(5)=mrho
      kpipl(5)=mpi
      kpimn(5)=mpi
      call geltwobdo('set:a',a)
      call geltwobdo('set:b',b)
      call geltwob(krho,kpipl,kpimn)
      ppl(1)=kpipl(1)
      ppl(2)=kpipl(2)
      ppl(3)=kpipl(3)
      pmn(1)=kpimn(1)
      pmn(2)=kpimn(2)
      pmn(3)=kpimn(3)
      return
      end
c
      subroutine gelrhom(energy,mpi,mrho)
      implicit none
      real *4 energy !available energy
      real *4 mpi    !pi mass
      real *4 mrho   !rho mass picked
      real *4 gamma/0.1512/  !width of rho
      real *4 mrho0/0.7699/  !central mass value
      real *4 temp(2),max,test
      call gelbw(mrho0,mrho0,gamma,max)
      if(energy.lt.mrho0) call gelbw(energy,mrho0,gamma,max)
 1000 continue
      call gelrndm(temp,2)
      mrho=2.0*mpi+(energy-2.0*mpi)*temp(1)
      call gelbw(mrho,mrho0,gamma,test)
      test=test/max
      if(temp(2).lt.test) go to 1099
      go to 1000
 1099 continue
      return
      end
c
      subroutine gelbw(e,e0,width,value)
      implicit none
      real *4 e,e0,width,value
      real *4 widthsq,de,desq
      integer *4 mode/0/
      character *(*) command
      integer *4 val
      de=e-e0
      if(mode.eq.0) then
         widthsq=width**2
         desq=de**2
         value=widthsq/(desq+widthsq/4.0)
      else
         write(6,*) 'gelpw: undefined BW model=',mode
      endif !mode.eq.0
      return
c
      entry gelbwset(command,val)
      if(command.eq.'set:mode') then
         mode=val
      endif
      return
c
      end
c
      subroutine gelrndm(a,n) !select random generator
      implicit none
      save
      real *4 a(*)              !array of random numbers
      integer *4 n              !number of random numbers wanted
      integer *4 mode/1/        !=1=>use |grndm| 2=>use |begran|
      real *4 begran            !|beget| random number generator
      integer *4 newmode        !new mode to be set
      integer *4 i              !dummy index
      if(mode.eq.1) then
         call grndm(a,n)
      else if(mode.eq.2) then
         do 1000 i=1,n
            a(i)=begran()
 1000    continue
      else
         write(6,*) 'gelrndm: mode=',mode,' undefined mode'
      endif !mode.eq.1
      return
c
      entry gelrndmset(newmode)
      mode=newmode
      return
c
      end
