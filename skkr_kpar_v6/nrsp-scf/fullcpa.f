c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine tata(conc,n,mc,ma,mb)
c====================
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
      parameter (ndim=lmmaxp)
c
      complex*16 mc(ndim,ndim),ma(ndim,ndim),mb(ndim,ndim),det
c
      cca=conc
      ccb=1.d0-cca
      do i=1,n
      do j=1,n
        mc(i,j)=cca*mb(i,j)+ccb*ma(i,j)
      end do
      end do
      call gjinv(mc,n,ndim,det)
      call tripmt(mb,mc,ma,n,n,ndim)
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
      complex*16 mc(ndim,ndim,mintfc),ma(ndim,ndim,mintfc),
     >mb(ndim,ndim,mintfc),tauc(ndim,ndim,mintfc),
     >taua(ndim,ndim,mintfc),taub(ndim,ndim,mintfc)
c
      data tiny/1.0d-6/,tol/1.0d-10/
c
      error=0.d0
      do li=1,nintfc
        cpalay=1.d0-conc(li).gt.tiny
        call compmat(ma(1,1,li),mc(1,1,li),n,ndim,tol,ica)
        call compmat(mb(1,1,li),mc(1,1,li),n,ndim,tol,icb)
        if((.not.cpalay).or.(ica+icb.ge.1)) then
          call repl(taua(1,1,li),tauc(1,1,li),n,ndim)
          call repl(taub(1,1,li),tauc(1,1,li),n,ndim)
          goto 10
        end if
        call cpasol(conc(li),mc(1,1,li),ma(1,1,li),mb(1,1,li),
     >              tauc(1,1,li),taua(1,1,li),taub(1,1,li),
     >              n,error1,error2,error3)
        if(cpatest.and.cpalay)
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
      complex*16 mc(ndim,ndim),ma(ndim,ndim),mb(ndim,ndim)
      complex*16 tauc(ndim,ndim),taua(ndim,ndim),taub(ndim,ndim)
      complex*16 h(ndim,ndim),da(ndim,ndim),db(ndim,ndim)
      complex*16 xc(ndim,ndim),xa(ndim,ndim),xb(ndim,ndim),det
c
      data cone/(1.d0,0.d0)/
c
c #####
c
      cca=conc
      ccb=1.d0-cca
c
c delta t**(-1) matrices
      call repl(da,mc,n,ndim)
      call repl(db,mc,n,ndim)
      call submat(da,ma,n,ndim)
      call submat(db,mb,n,ndim)
c
c excess scattering matrices: xa, xb
      call repl(xa,da,n,ndim)
      call repl(xb,db,n,ndim)
      call gjinv(xa,n,ndim,det)
      call gjinv(xb,n,ndim,det)
      call submat(xa,tauc,n,ndim)
      call submat(xb,tauc,n,ndim)
      call gjinv(xa,n,ndim,det)
      call gjinv(xb,n,ndim,det)
c
      error1=0.d0
      do i=1,n
      do j=1,n
c average excess scattering matrix: xc
        xc(i,j)=cca*xa(i,j)+ccb*xb(i,j)
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
      call repl(xc,tauc,n,ndim)
      call tripmt(da,xc,db,n,n,ndim)
      error2=0.d0
      d1=0.d0
      d2=0.d0
      do i=1,n
      do j=1,n
        xc(i,j)=xc(i,j)-cca*da(i,j)-ccb*db(i,j)
        dummy=cdabs(xc(i,j))
        if(dummy.gt.error2) error2=dummy
        d1=d1+h(i,j)*dconjg(h(i,j))
        d2=d2+mc(i,j)*dconjg(mc(i,j))
      end do
      end do
      error3=dsqrt(d1/d2)
c
c compute tau matrices with respect to species A and B
      call repl(taua,xa,n,ndim)
      call repl(taub,xb,n,ndim)
      call tripmt(tauc,taua,tauc,n,n,ndim)
      call tripmt(tauc,taub,tauc,n,n,ndim)
      call addmat(taua,tauc,n,ndim)
      call addmat(taub,tauc,n,ndim)
c
      return
      end
