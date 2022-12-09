c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine gf(lmax,reg,dx,ns,rs,gz,fz,gj,fj,glz,flz,glj,flj,
     >              tau,ggrmat,fgrmat)    
c===================
c
c input:  lmax  - maximum of angular momentum index
c         reg - if .true. only regular part of the Green's function is considered
c         dx,ns,rs - as in 'readpot'
c         gz,fz   - large and small component of regular radial solution*r
c         gj,fj   - large and small component of irreg. radial solution*r
c         glz,flz,glj,flj - same for left-side solutions of Dirac equation
c         tau - site-diagonal tau matrix
c output: ggrmat, fgrmat - matrices of the radial integrals of the Green function 
c         tggrmat, tfgrmat - matrices of the radial integrals of the Green function with transformed electron-hole indices
c
c
      implicit none
c
      include '../param.h'
c
      logical reg
      integer lmax,ns
      real*8 rs,dx
      complex*16 gz(dbogomaxp,dbogomaxp,nrad)
      complex*16 fz(dbogomaxp,dbogomaxp,nrad)
      complex*16 gj(dbogomaxp,dbogomaxp,nrad)
      complex*16 fj(dbogomaxp,dbogomaxp,nrad)
      complex*16 glz(dbogomaxp,dbogomaxp,nrad)
      complex*16 flz(dbogomaxp,dbogomaxp,nrad)
      complex*16 glj(dbogomaxp,dbogomaxp,nrad)
      complex*16 flj(dbogomaxp,dbogomaxp,nrad)
      complex*16 tau(dbogomaxp,dbogomaxp)
      complex*16 ggrmat(dbogomaxp,dbogomaxp,0:lsup)
      complex*16 fgrmat(dbogomaxp,dbogomaxp,0:lsup)
c
      integer kmymax,i,kmy,kmyp,kmy1,kmyp1,lambda,dbogomax
      real*8 x,rad(nrad),rlam,rsimp
      real*8 rggr(nrad,dbogomaxp,dbogomaxp)
      real*8 iggr(nrad,dbogomaxp,dbogomaxp)
      real*8 rfgr(nrad,dbogomaxp,dbogomaxp)
      real*8 ifgr(nrad,dbogomaxp,dbogomaxp)
      real*8 rgg(nrad),igg(nrad),rfg(nrad),ifg(nrad)
      real*8 r1gg,i1gg,r1fg,i1fg
      complex*16 cggr,cfgr,irreg,detl
c      complex*16 ctrafo(dbogomaxp,dbogomaxp)
c      complex*16 ctrafo2(dbogomaxp,dbogomaxp)
c
      external rsimp
c
c      call conjinkmy(ctrafo)
c      ctrafo2=ctrafo
c      call gjinv(ctrafo2,dbogomaxp,dbogomaxp,detl)
c
      irreg=(1.0d0,0.0d0)
      if(reg) irreg=(0.0d0,0.0d0)
      ggrmat=(0.d0,0.d0)
      fgrmat=(0.d0,0.d0)
c
      kmymax=2*(lmax+1)*(lmax+1) 
      dbogomax=2*kmymax
      x=dlog(rs)-(ns-1)*dx
      do i=1,ns
         rad(i)=dexp(x)
         x=x+dx
      end do
c
c calculate matrices of radial integrals
c
      do i=1,ns
c
        do kmy=1,dbogomax
        do kmyp=1,dbogomax
c
          cggr=(0.d0,0.d0)
          cfgr=(0.d0,0.d0)
          do kmy1=1,dbogomax
          do kmyp1=1,dbogomax
            ! G ~ Z^R tau Z^L - J^R Z^L
            cggr=cggr+glz(kmy,kmy1,i)*gz(kmyp,kmyp1,i)*tau(kmyp1,kmy1)  ! Z^R tau Z^L
            cfgr=cfgr+flz(kmy,kmy1,i)*fz(kmyp,kmyp1,i)*tau(kmyp1,kmy1)  ! Z^R tau Z^L
          end do
          cggr=cggr-irreg*glz(kmy,kmy1,i)*gj(kmyp,kmy1,i)   ! - J^R Z^L 
          cfgr=cfgr-irreg*flz(kmy,kmy1,i)*fj(kmyp,kmy1,i)   ! - J^R Z^L 
          end do
c
          rggr(i,kmy,kmyp)=dreal(cggr)
          iggr(i,kmy,kmyp)=dimag(cggr)
          rfgr(i,kmy,kmyp)=dreal(cfgr)
          ifgr(i,kmy,kmyp)=dimag(cfgr)
c
        end do
        end do
c
      end do
c
      ggrmat=(0d0,0d0)
      fgrmat=(0d0,0d0)
      do kmy=1,dbogomax
      do kmyp=1,dbogomax
      do lambda=0,2*lmax
c
        do i=1,ns
          rlam=rad(i)**lambda
          rgg(i)=rlam*rggr(i,kmy,kmyp)
          igg(i)=rlam*iggr(i,kmy,kmyp)
          rfg(i)=rlam*rfgr(i,kmy,kmyp)
          ifg(i)=rlam*ifgr(i,kmy,kmyp)
        end do
        r1gg=rsimp(rgg,rad,ns,dx)
        i1gg=rsimp(igg,rad,ns,dx)
        r1fg=rsimp(rfg,rad,ns,dx)
        i1fg=rsimp(ifg,rad,ns,dx)
c
        ggrmat(kmyp,kmy,lambda)=dcmplx(r1gg,i1gg)
        fgrmat(kmyp,kmy,lambda)=dcmplx(r1fg,i1fg)
c         
      end do
      end do
      end do
c
c      do lambda=0,2*lmax
c         ggrmat(1:dbogomax,1:dbogomax,lambda)=
c     >        matmul(ggrmat(1:dbogomax,1:dbogomax,lambda),ctrafo)
c         fgrmat(1:dbogomax,1:dbogomax,lambda)=
c     >        matmul(fgrmat(1:dbogomax,1:dbogomax,lambda),ctrafo)
c         ggrmat(1:dbogomax,1:dbogomax,lambda)=
c     >        matmul(ctrafo2,ggrmat(1:dbogomax,1:dbogomax,lambda))
c         fgrmat(1:dbogomax,1:dbogomax,lambda)=
c     >        matmul(ctrafo2,fgrmat(1:dbogomax,1:dbogomax,lambda))
c      end do
c
      return
      end
c
c
c
      subroutine trafogeh(ggrmat,fgrmat,deltakmy,lmax,ggeh,fgeh)
c
      implicit none
c
      include '../param.h'
c
      integer lmax,kmymax,lambda
      complex*16 ggrmat(dbogomaxp,dbogomaxp,0:lsup)
      complex*16 fgrmat(dbogomaxp,dbogomaxp,0:lsup)
      complex*16 ggeh(kmymaxp,kmymaxp,0:lsup)
      complex*16 fgeh(kmymaxp,kmymaxp,0:lsup)
      complex*16 deltakmy(kmymaxp,kmymaxp)
c
      kmymax=2*(lmax+1)*(lmax+1) 
c
      do lambda=0,2*lmax
c
       call doubmt1(deltakmy,ggrmat(1:kmymax,kmymax+1:2*kmymax,lambda),
     >       ggeh(1:kmymax,1:kmymax,lambda),kmymax,kmymaxp)
       call doubmt1(deltakmy,fgrmat(1:kmymax,kmymax+1:2*kmymax,lambda),
     >       fgeh(1:kmymax,1:kmymax,lambda),kmymax,kmymaxp)
c
      end do
      return
      end
c
      subroutine trafoghe(ggrmat,fgrmat,deltakmy,lmax,gghe,fghe)
c
      implicit none
c
      include '../param.h'
c
      integer lmax,kmymax,lambda
      complex*16 ggrmat(dbogomaxp,dbogomaxp,0:lsup)
      complex*16 fgrmat(dbogomaxp,dbogomaxp,0:lsup)
      complex*16 gghe(kmymaxp,kmymaxp,0:lsup)
      complex*16 fghe(kmymaxp,kmymaxp,0:lsup)
      complex*16 deltakmy(kmymaxp,kmymaxp)
c
      kmymax=2*(lmax+1)*(lmax+1) 
c
      do lambda=0,2*lmax
c
       call doubmt1(deltakmy,ggrmat(kmymax+1:2*kmymax,1:kmymax,lambda),
     >       gghe(1:kmymax,1:kmymax,lambda),kmymax,kmymaxp)
       call doubmt1(deltakmy,fgrmat(1:kmymax,kmymax+1:2*kmymax,lambda),
     >       fghe(1:kmymax,1:kmymax,lambda),kmymax,kmymaxp)
c
      end do
      return
      end
c
