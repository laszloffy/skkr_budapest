c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine tata(conc,n,mc,ma,mb)
c====================
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
      parameter (ndim=lmmaxp)
c
      real*8 cc(4)
      complex*16 mc(ndim,ndim),ma(ndim,ndim,2),mb(ndim,ndim,2)
      complex*16 m(ndim,ndim,4),det
c
      cc(1)=0.5d0*conc
      cc(2)=0.5d0*conc
      cc(3)=0.5d0*(1.0d0-conc)
      cc(4)=0.5d0*(1.0d0-conc)
c
      call repl(m(1,1,1),ma(1,1,1),n,ndim)
      call repl(m(1,1,2),ma(1,1,2),n,ndim)
      call repl(m(1,1,3),mb(1,1,1),n,ndim)
      call repl(m(1,1,4),mb(1,1,2),n,ndim)
c
      do ic=1,4
        call gjinv(m(1,1,ic),n,ndim,det)
      end do
c
      do i=1,n
      do j=1,n
        mc(i,j)=(0.0d0,0.0d0)
        do ic=1,4
          mc(i,j)=mc(i,j)+cc(ic)*m(i,j,ic)
        end do
      end do
      end do
      call gjinv(mc,n,ndim,det)
c
      return
      end
      subroutine cpacor(conc,mc,ma,mb,tauc,taua,taub,n,nintfc,
     >itcpa,itcpam,cpatol,cpacon,concpa,cpatest,error)
c======================
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
      parameter (ndim=lmmaxp)
c
      logical cpacon,concpa,cpatest,cpalay
c
      dimension conc(mintfc)
      complex*16 mc(ndim,ndim,mintfc),ma(ndim,ndim,2,mintfc),
     >mb(ndim,ndim,2,mintfc),tauc(ndim,ndim,mintfc),
     >taua(ndim,ndim,2,mintfc),taub(ndim,ndim,2,mintfc)
c
      data tiny/1.0d-6/,tol/1.0d-10/
c
      error=0.d0
      do li=1,nintfc
        call compmat(ma(1,1,1,li),mc(1,1,li),n,ndim,tol,ica1)
        call compmat(ma(1,1,2,li),mc(1,1,li),n,ndim,tol,ica2)
        call compmat(mb(1,1,1,li),mc(1,1,li),n,ndim,tol,icb1)
        call compmat(mb(1,1,2,li),mc(1,1,li),n,ndim,tol,icb2)
        if((ica1+ica2+icb1+icb2.ge.1)) then
          call repl(taua(1,1,1,li),tauc(1,1,li),n,ndim)
          call repl(taua(1,1,2,li),tauc(1,1,li),n,ndim)
          call repl(taub(1,1,1,li),tauc(1,1,li),n,ndim)
          call repl(taub(1,1,2,li),tauc(1,1,li),n,ndim)
          goto 10
        end if
        call cpasol(conc(li),mc(1,1,li),ma(1,1,1,li),mb(1,1,1,li),
     >              tauc(1,1,li),taua(1,1,1,li),taub(1,1,1,li),
     >              n,error1,error2,error3)
        if(cpatest)
     >  write(6,100) itcpa,li,error1,error2,error3
c
        if(error3.ge.error) error=error3
  10    continue
      end do
      if(error.le.cpatol) concpa=.false.
      if(itcpa.ge.itcpam) cpacon=.false.
      itcpa=itcpa+1
c
  100 format(' Icpa',i2,'  L',i3,3x,1pd13.6,2x,1pd13.6,2x,1pd13.6)
c
      return
      end
      subroutine cpasol(conc,mc,ma,mb,tauc,taua,taub,n,
     >                  error1,error2,error3)
c======================
c
c update cpa t-matrix via b. ginatempo and j.b. staunton,
c                         j. phys. f 18, 1827 (1988)
c
      implicit real*8(a-h,o-z)
c
      include '../param.h'
      parameter (ndim=lmmaxp)
c
      complex*16 mc(ndim,ndim),ma(ndim,ndim,2),mb(ndim,ndim,2)
      complex*16 tauc(ndim,ndim),taua(ndim,ndim,2),taub(ndim,ndim,2)
c
      real*8 cc(4)
      complex*16 m(ndim,ndim,4),x(ndim,ndim,4),tau(ndim,ndim,4)
      complex*16 h(ndim,ndim),xc(ndim,ndim)
      complex*16 det
c
      data cone/(1.d0,0.d0)/
c
c #####
c
      cc(1)=0.5d0*conc
      cc(2)=0.5d0*conc
      cc(3)=0.5d0*(1.0d0-conc)
      cc(4)=0.5d0*(1.0d0-conc)
      call repl(m(1,1,1),ma(1,1,1),n,ndim)
      call repl(m(1,1,2),ma(1,1,2),n,ndim)
      call repl(m(1,1,3),mb(1,1,1),n,ndim)
      call repl(m(1,1,4),mb(1,1,2),n,ndim)
c
c excess scattering matrices: x
      do ic=1,4
        call repl(x(1,1,ic),mc,n,ndim)
        call submat(x(1,1,ic),m(1,1,ic),n,ndim)
        call gjinv(x(1,1,ic),n,ndim,det)
        call submat(x(1,1,ic),tauc,n,ndim)
        call gjinv(x(1,1,ic),n,ndim,det)
      end do
c
c average excess scattering matrix: xc
      error1=0.d0
      do i=1,n
      do j=1,n
        xc(i,j)=(0.0d0,0.0d0)
        do ic=1,4
          xc(i,j)=xc(i,j)+cc(ic)*x(i,j,ic)
        end do
        dummy=cdabs(xc(i,j))
        if(dummy.gt.error1) error1=dummy
      end do
      end do
c
c compute next guess for effective t**(-1)
      call repl(h,xc,n,ndim)
      call gjinv(h,n,ndim,det)
      call addmat(h,tauc,n,ndim)
      call gjinv(h,n,ndim,det)
      call submat(mc,h,n,ndim)
c
c error in conventional cpa conditions and 
c subsequent t**(-1) matrices
      d1=0.d0
      d2=0.d0
      do i=1,n
      do j=1,n
        d1=d1+h(i,j)*dconjg(h(i,j))
        d2=d2+mc(i,j)*dconjg(mc(i,j))
      end do
      end do
      error2=dsqrt(d2)
      error3=dsqrt(d1/d2)
c
c compute tau matrices with respect to species A and B
      do ic=1,4
        call repl(tau(1,1,ic),x(1,1,ic),n,ndim)
        call tripmt(tauc,tau(1,1,ic),tauc,n,n,ndim)
        call addmat(tau(1,1,ic),tauc,n,ndim)
      end do
c
      call repl(taua(1,1,1),tau(1,1,1),n,ndim)
      call repl(taua(1,1,2),tau(1,1,2),n,ndim)
      call repl(taub(1,1,1),tau(1,1,3),n,ndim)
      call repl(taub(1,1,2),tau(1,1,4),n,ndim)
c
      return
      end
