c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine wafu(ce,cl,soc,lmax,idpot,v0,E_Fermi,
     >           singrat,urat,drat,delta,vr,br,rb,dx,rs,ns,
     >           tminv,gz,fz,gj,fj,glz,flz,glj,flj,iflag)
!====================
!
! input:  e - any complex energy
!         cl - light velocity
!         soc - spin-orbit scaling
!         lmax  - maximum of angular momentum index
!         v0 - vacuum potential
!         idpot,vr,br,dx,ns,rs - as in 'readpot'
!         iflag=0 only tminv is calculated
! output: tminv - inverse of single-site scattering matrix
!         gz,fz - large and small component of right-side regular radial solution*r
!         gj,fj - large and small component of right-side irregular radial solution*r
!         glz,flz - large and small component of left-side regular radial solution*r
!         glj,flj - large and small component of left-side irregular radial solution*r
!
      implicit none
!
      include '../param.h'
!
! I/O variables
      character*10, intent(in)  :: idpot
      integer,      intent(in)  :: lmax,iflag,ns
      real*8,       intent(in)  :: cl,soc,E_Fermi
      real*8,       intent(in)  :: vr(nrad),br(nrad),rb(3),v0
      real*8,       intent(in)  :: dx,rs
      real*8,       intent(in)  :: singrat,urat,drat
      complex*16,   intent(in)  :: ce,delta(nrad)
      complex*16,   intent(out) :: tminv(dbogomaxp,dbogomaxp)
      complex*16,   intent(out) :: gz(dbogomaxp,dbogomaxp,nrad)
      complex*16,   intent(out) :: fz(dbogomaxp,dbogomaxp,nrad)
      complex*16,   intent(out) :: gj(dbogomaxp,dbogomaxp,nrad)
      complex*16,   intent(out) :: fj(dbogomaxp,dbogomaxp,nrad)
      complex*16,   intent(out) :: glz(dbogomaxp,dbogomaxp,nrad)
      complex*16,   intent(out) :: flz(dbogomaxp,dbogomaxp,nrad)
      complex*16,   intent(out) :: glj(dbogomaxp,dbogomaxp,nrad)
      complex*16,   intent(out) :: flj(dbogomaxp,dbogomaxp,nrad)
!
c
      external rsimp
c
!
      integer lml(lmmaxp),lmm(lmmaxp)
      common/mome/lml,lmm
c
!
      write(6,*) 'Ask for the original version of this file from the 
     >owners of the Copyright, e.g. from laszloffy.andras(at)wigner.hu'
!
      return
      end subroutine wafu

