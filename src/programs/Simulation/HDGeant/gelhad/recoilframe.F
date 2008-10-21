!Mon Mar  3 10:07:47 PST 1997
!Find momentum of stand-in particle in target rest from
!A. Snyder
!
      subroutine recoilframe
     .(egamma,mstandin,mtarget,mrecoil,pstandin,precoil,ok)
      implicit none
      real *4 egamma !energy of input gamma
      real *4 mstandin !mass of stand-in particle
      real *4 mtarget !mass of target
      real *4 pstandin !momentum of stand-in in recoil rest frame
      real *4 mrecoil !mass of recoil fragment
      real *4 precoil !momentum of recoil in lab
c
      real *4 mall,eall,pall,ercm,prcm,escm,pscm,beta,gamma
      logical ok
c
      ok=.false. !assume the worst
c
c cm system
      eall=egamma+mtarget
      pall=egamma
      if(eall.lt.pall) return
      mall=sqrt(eall**2-pall**2)
c
c do recoil
      ercm=(mall**2+mrecoil**2-mstandin**2)/(2.0*mall)
      if(ercm.lt.mrecoil) return
      prcm=sqrt(ercm**2-mrecoil**2)
      beta=pall/eall
      gamma=eall/mall
      precoil=gamma*(-prcm+beta*ercm)
c
c do projectile
      escm=(mall**2+mstandin**2-mrecoil**2)/(2.0*mall)
      if(escm.lt.mstandin) return
      pscm=sqrt(escm**2-mstandin**2)
      beta=-prcm/ercm
      gamma=ercm/mrecoil
      pstandin=gamma*(pscm-beta*escm)
c
      ok=.true. !success
c
      return
      end
