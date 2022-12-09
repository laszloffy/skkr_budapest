c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine dens(lmax,rs,dx,ns,rb,gz,fz,gj,fj,glz,flz,glj,flj,tau,
     >                zrhoeh,zrho)
c===================
c
c input:  lmax  - maximum of angular momentum index
c         rs,dx,ns - as in 'readpot'
c         rb - orientation of exchange field; either z (localmode) or real one
c         gz,fz - large and small component of regular radial solution * r
c         gj,fj - large and small component of irreg. radial solution * r
c         glz,flz,glj,flj - same for left-side solutions to the Dirac equation
c         tau - site-diagonal tau matrix
c output: zrho - radial density distribution
c         zrhosp - spin resolution of radial density distribution
c         zrhodsp - only l=2 subspace is considered
c         zmagrho - radial magnetic density distribution
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      dimension f1(nrad),f2(nrad),r(nrad)
c
      real*8     rb(3)
      complex*16 sum1,sum2
      complex*16 gz(dbogomaxp,dbogomaxp,nrad)
      complex*16 fz(dbogomaxp,dbogomaxp,nrad)
      complex*16 glz(dbogomaxp,dbogomaxp,nrad)
      complex*16 flz(dbogomaxp,dbogomaxp,nrad)
      complex*16 gz0(dbogomaxp,dbogomaxp)
      complex*16 fz0(dbogomaxp,dbogomaxp)
      complex*16 gj0(dbogomaxp,dbogomaxp)
      complex*16 fj0(dbogomaxp,dbogomaxp)
      complex*16 gj(dbogomaxp,dbogomaxp,nrad)
      complex*16 fj(dbogomaxp,dbogomaxp,nrad)
      complex*16 glj(dbogomaxp,dbogomaxp,nrad)
      complex*16 flj(dbogomaxp,dbogomaxp,nrad)
      complex*16 glzt(dbogomaxp,dbogomaxp)
      complex*16 flzt(dbogomaxp,dbogomaxp)
      complex*16 gljt(dbogomaxp,dbogomaxp)
      complex*16 fljt(dbogomaxp,dbogomaxp)
      complex*16 zrho(nrad),zrhosp(nrad,2),zrhodsp(nrad,2),zmagrho(nrad)
      complex*16 zrhoeh(nrad),zrhohe(nrad)
      complex*16 tau(dbogomaxp,dbogomaxp)
      complex*16 trmat(dbogomaxp,dbogomaxp)
      complex*16 gmat(dbogomaxp,dbogomaxp),fmat(dbogomaxp,dbogomaxp)
      complex*16 gmat1(dbogomaxp,dbogomaxp),fmat1(dbogomaxp,dbogomaxp)
      complex*16 rmat(dbogomaxp,dbogomaxp),rlms(kmymaxp),cf(nrad)
      complex*16 rlmseh(kmymaxp),rlmshe(kmymaxp)
      complex*16 deltakmy(kmymaxp,kmymaxp)
c
      complex*16 scoeff(kmymaxp,kmymaxp),sbcoeff(kmymaxp,kmymaxp)
      complex*16 sxcoeff(kmymaxp,kmymaxp),sxbcoeff(kmymaxp,kmymaxp)
      complex*16 sycoeff(kmymaxp,kmymaxp),sybcoeff(kmymaxp,kmymaxp)
      complex*16 szcoeff(kmymaxp,kmymaxp),szbcoeff(kmymaxp,kmymaxp)
      common/sigmat/sxcoeff,sxbcoeff,sycoeff,sybcoeff,szcoeff,szbcoeff
c
      data small/1.0d-12/
c
      pi=4.d0*datan(1.d0)
      fac=-1.d0/pi
      lmmax=(lmax+1)*(lmax+1)
      kmymax=2*lmmax
      dbogomax=2*kmymax
      zrho=(0.d0,0.d0)
      zrhosp=(0.d0,0.d0)
      zrhodsp=(0.d0,0.d0)
      zmagrho=(0.d0,0.d0)
c
      !handle optional global mode: sigma*rb operator needed for zmagrho
      scoeff = rb(1)*sxcoeff+rb(2)*sycoeff+rb(3)*szcoeff
      sbcoeff = rb(1)*sxbcoeff+rb(2)*sybcoeff+rb(3)*szbcoeff
c
      do i=1,ns
c
        do kmy=1,dbogomax
        do kmyp=1,dbogomax
          gz0(kmy,kmyp)=gz(kmy,kmyp,i)
          fz0(kmy,kmyp)=fz(kmy,kmyp,i)
          gj0(kmy,kmyp)=gj(kmy,kmyp,i)
          fj0(kmy,kmyp)=fj(kmy,kmyp,i)
          glzt(kmy,kmyp)=glz(kmyp,kmy,i)
          flzt(kmy,kmyp)=flz(kmyp,kmy,i)
          gljt(kmy,kmyp)=glj(kmyp,kmy,i)
          fljt(kmy,kmyp)=flj(kmyp,kmy,i)
        end do
        end do
c
        ! G ~ Z^R tau Z^L - J^R Z^L, permute to have (Z^L Z^R) * tau - Z^L J^R
        gmat(1:dbogomax,1:dbogomax)=
     >   matmul(glzt(1:dbogomax,1:dbogomax),gz0(1:dbogomax,1:dbogomax))
        fmat(1:dbogomax,1:dbogomax)=
     >   matmul(flzt(1:dbogomax,1:dbogomax),fz0(1:dbogomax,1:dbogomax))
        gmat(1:dbogomax,1:dbogomax)=
     >   matmul(gmat(1:dbogomax,1:dbogomax),tau(1:dbogomax,1:dbogomax))
        fmat(1:dbogomax,1:dbogomax)=
     >   matmul(fmat(1:dbogomax,1:dbogomax),tau(1:dbogomax,1:dbogomax))
        gmat(1:dbogomax,1:dbogomax)=gmat(1:dbogomax,1:dbogomax)-
     >   matmul(gljt(1:dbogomax,1:dbogomax),gz0(1:dbogomax,1:dbogomax))   ! old version had Z^R J^L...
!     >   matmul(glzt(1:dbogomax,1:dbogomax),gj0(1:dbogomax,1:dbogomax))  ! new version with J^R Z^L
        fmat(1:dbogomax,1:dbogomax)=fmat(1:dbogomax,1:dbogomax)-
     >   matmul(fljt(1:dbogomax,1:dbogomax),fz0(1:dbogomax,1:dbogomax))   ! old version had Z^R J^L...
!     >   matmul(flzt(1:dbogomax,1:dbogomax),fj0(1:dbogomax,1:dbogomax))  ! new version with J^R Z^L

        rmat(1:dbogomax,1:dbogomax)=
     >   gmat(1:dbogomax,1:dbogomax)+fmat(1:dbogomax,1:dbogomax)
c
c     Transform anomalous part
c
       call pairing(deltakmy)
c
       trmat(1:dbogomax,1:dbogomax) = rmat(1:dbogomax,1:dbogomax)
c
       call doubmt1(rmat(1:kmymax,kmymax+1:2*kmymax),
     >       deltakmy(1:kmymaxp,1:kmymaxp),
     >       trmat(1:kmymax,kmymax+1:2*kmymax),kmymax,kmymaxp)
c
       call doubmt1(rmat(kmymax+1:2*kmymax,1:kmymax),
     >       deltakmy(1:kmymaxp,1:kmymaxp),
     >       trmat(kmymax+1:2*kmymax,1:kmymax),kmymax,kmymaxp)
c
c
c        call replms(rlmseh,trmat(1:kmymax,kmymax+1:2*kmymax),lmax)
c        call replms(rlmshe,trmat(kmymax+1:2*kmymax,1:kmymax),lmax)
c
c anomalous densities
        sum1=(0.d0,0.d0)
        do kmy=1,kmymax
          sum1=sum1+trmat(kmy,kmy+kmymax)
        end do
        zrhoeh(i)=fac*sum1
c
        sum1=(0.d0,0.d0)
        do kmy=1,kmymax
          sum1=sum1+rmat(kmy+kmymax,kmy)
        end do
        zrho(i)=fac*sum1
c
c spin-projected charge density
c        sum1=(0.d0,0.d0)
c        sum2=(0.d0,0.d0)
c        do k=1,lmmax
c          sum1=sum1+rlms(k)
c          sum2=sum2+rlms(k+lmmax)
c        end do
c        zrhosp(i,1)=fac*sum1
c        zrhosp(i,2)=fac*sum2
c
c spin-projected d-like (l=2) charge density
c        sum1=(0.d0,0.d0)
c        sum2=(0.d0,0.d0)
c        do k=5,9
c          sum1=sum1+rlms(k)
c          sum2=sum2+rlms(k+lmmax)
c        end do
c        zrhodsp(i,1)=fac*sum1
c        zrhodsp(i,2)=fac*sum2
c
c spin-only magnetization density
c in case of local mode or z direction: szcoeff/szbcoeff only
c        call doubmt1(gmat,scoeff,gmat1,kmymax,kmymaxp)
c        call doubmt1(fmat,sbcoeff,fmat1,kmymax,kmymaxp)
c        call submat1(gmat1,fmat1,rmat,kmymax,kmymaxp)
c        sum1=(0.d0,0.d0)
c        do kmy=1,kmymax
c          sum1=sum1+rmat(kmy,kmy)
c        end do
c        zmagrho(i)=fac*sum1


!
! fix zrhosp for global version: up/down component is meaningless as-is
! either 1. rotate GF to local frame before computing zrhosp the old way,
!     or 2. recompute zrhosp from zmagrho ~ zrhosp_up - zrhosp_down, zrho ~ zrhosp_up + zrhosp_down
!
!   we'll go with number 2 for now, TODO number 1 for backwards compatibility/general well-workingness?
! !TODO remove original stuff above, a.k.a. clean-up
c        zrhosp(i,1) = (zrho(i)-zmagrho(i))/2.0  ! down component
c        zrhosp(i,2) = (zrho(i)+zmagrho(i))/2.0  ! up component
c
      end do
c
      return
      end
