c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine repl(x,y,n,ndim)
c====================
c
c    x=y
c
      complex*16 x(ndim,ndim),y(ndim,ndim)
c
      x(1:n,1:n) = y(1:n,1:n)
      return
      end
      subroutine repldim(x,y,n,ndim,mdim)
c====================
c
c    x=y
c
      complex*16 x(ndim,ndim),y(mdim,mdim)
c
      x(1:n,1:n) = y(1:n,1:n)
      return
      end
      subroutine replt(x,y,n,ndim)
c====================
c
c    x=y^T
c
      complex*16 x(ndim,ndim),y(ndim,ndim)
c
      do 1 i=1,n
      do 1 j=1,n
    1 x(i,j)=y(j,i)
      return
      end
      subroutine replrel(x,y,n,ndim)
c=======================
c
c    x(k,k')=(-)**(l-l')* y(k',k)*
c
      include '../param.h'
      complex*16 x(ndim,ndim),y(ndim,ndim)
c
      if(n.gt.50) stop 'Increase index of ldex in reprel'
c
      do 1 i=1,n
        li=ldex(i)
      do 1 j=1,n
        lj=ldex(j)
        x(i,j)=dconjg(y(j,i))
        if(mod(iabs(li-lj),2).eq.1) x(i,j)=-x(i,j)
    1 continue
      return
      end
      subroutine compmat(x,y,n,ndim,tol,ic)
c=======================
c
c check whether  x=y
c
      implicit real*8 (a-h,o-z)
      complex*16 x(ndim,ndim),y(ndim,ndim),xx,dx
      data tol0/1.0d-12/
c
      if(tol.le.0.d0) tol=tol0
      ic=0
      xnorm=0.d0
      dxnorm=0.d0
      do i=1,n
      do j=1,n
        xx=x(i,j)
        dx=x(i,j)-y(i,j)
        xxre=dreal(xx)
        xxim=dimag(xx)
        dxre=dreal(dx)
        dxim=dimag(dx)
        xnorm=xnorm+xxre*xxre+xxim*xxim
        dxnorm=dxnorm+dxre*dxre+dxim*dxim
      end do
      end do
      xnorm=dsqrt(xnorm)
      dxnorm=dsqrt(dxnorm)
      if(dxnorm.lt.tol*xnorm) ic=1
      return
      end
      subroutine addmat(x,y,n,ndim)
c====================
c
c    x=x+y
c
      complex*16 x(ndim,ndim),y(ndim,ndim)
c
      x(1:n,1:n) = x(1:n,1:n) + y(1:n,1:n)
      return
      end
c
c
      subroutine addmatc(x,y,c,n,ndim)
c====================
c
c    x=x+c*y
c
      complex*16 x(ndim,ndim),y(ndim,ndim),c
c
      x(1:n,1:n) = x(1:n,1:n) + c*y(1:n,1:n)
      end subroutine addmatc
c
c
      subroutine addmat1(a,b,c,n,ndim)
c====================
c
c    c=a+b
c
      complex*16 a(ndim,ndim),b(ndim,ndim),c(ndim,ndim)
c
      c(1:n,1:n) = a(1:n,1:n) + b(1:n,1:n)
      return
      end
      subroutine submat(x,y,n,ndim)
c====================
c
c    x=x-y
c
      complex*16 x(ndim,ndim),y(ndim,ndim)
c
      x(1:n,1:n) = x(1:n,1:n) - y(1:n,1:n)
      return
      end
      subroutine submat1(a,b,c,n,ndim)
c====================
c
c    c=a-b
c
      complex*16 a(ndim,ndim),b(ndim,ndim),c(ndim,ndim)
c
      c(1:n,1:n) = a(1:n,1:n) - b(1:n,1:n)
      return
      end
      subroutine symmat(x,n,ndim)
c======================
c Symmetrize matrix x:
c    x -> 1/2*(x+xT)
c
      complex*16 x(ndim,ndim)
c
      do 1 i=1,n
      do 1 j=i+1,n
      x(i,j)=0.5d0*(x(i,j)+x(j,i))
    1 x(j,i)=x(i,j)
      return
      end
      subroutine doubmt(amat,bmat,ndim,ndimp)
c======================
c
c  amat=amat*bmat
c
      implicit real*8 (a-h,o-z)
c
      complex*16 amat(ndimp,ndimp),bmat(ndimp,ndimp)
c
      amat(1:ndim,1:ndim) = matmul(amat(1:ndim,1:ndim),
     >                             bmat(1:ndim,1:ndim))
c
      return
      end
      subroutine doubmt1(amat,bmat,cmat,ndim,ndimp)
c======================
c
c  cmat=amat*bmat
c
      implicit real*8 (a-h,o-z)
c
      complex*16 amat(ndimp,ndimp),bmat(ndimp,ndimp),cmat(ndimp,ndimp)
c
      cmat(1:ndim,1:ndim) = matmul(amat(1:ndim,1:ndim),
     >                             bmat(1:ndim,1:ndim))
c
      return
      end
      subroutine doubmt2(amat,bmat,cmat,n1,n2,n3,ndimp)
c======================
c
c  cmat=amat*bmat
c
c amat:   ndi1 x ndi2
c bmat:   ndi2 x ndi3
c
      implicit real*8 (a-h,o-z)
c
      complex*16 amat(ndimp,ndimp),bmat(ndimp,ndimp),cmat(ndimp,ndimp)
c
      cmat(1:n1,1:n3) = matmul(amat(1:n1,1:n2),bmat(1:n2,1:n3))
c
      return
      end
      function trdbmt(amat,bmat,ndim,ndimp)
c===============================
c
c  trdbmt=Tr(amat*bmat)
c
      integer i,j,ndim,ndimp
      complex*16 trdbmt,amat(ndimp,ndimp),bmat(ndimp,ndimp)
c
      trdbmt=dcmplx(0.d0,0.d0)
      do i=1,ndim
      do j=1,ndim
        trdbmt=trdbmt+bmat(j,i)*amat(i,j)
      end do
      end do
c
      return
      end
      subroutine tripmtr(u,b,ust,ndi1,ndi2,ndim)
c =====================
c
c Triple product of real rectangular matrices
c b= u * b * ust
c
c u:    ndi1 x ndi1
c b:    ndi1 x ndi2
c ust:  ndi2 x ndi2
c
      implicit real*8 (a-h,o-z)
      dimension u(ndim,ndim),ust(ndim,ndim),b(ndim,ndim)
c
      b(1:ndi1,1:ndi2) = matmul(u(1:ndi1,1:ndi1),b(1:ndi1,1:ndi2))
      b(1:ndi1,1:ndi2) = matmul(b(1:ndi1,1:ndi2),ust(1:ndi2,1:ndi2))
c
      return
      end
      subroutine tripmt(u,b,ust,ndi1,ndi2,ndim)
c =====================
c
c Triple product of complex rectangular matrices
c b= u * b * ust
c
c u:    ndi1 x ndi1
c b:    ndi1 x ndi2
c ust:  ndi2 x ndi2
c
      implicit real*8 (a-h,o-z)
      complex*16 u(ndim,ndim),ust(ndim,ndim),b(ndim,ndim)
c
      b(1:ndi1,1:ndi2) = matmul(u(1:ndi1,1:ndi1),b(1:ndi1,1:ndi2))
      b(1:ndi1,1:ndi2) = matmul(b(1:ndi1,1:ndi2),ust(1:ndi2,1:ndi2))
c
      return
      end
      subroutine tripmt1(u,b,ust,b1,ndi1,ndi2,ndim)
c =====================
c
c Triple product of rectangular matrices
c b1= u * b * ust
c
c u:    ndi1 x ndi1
c b:    ndi1 x ndi2
c ust:  ndi2 x ndi2
c
      implicit real*8 (a-h,o-z)
      complex*16 u,ust,b,b1
      dimension u(ndim,ndim),ust(ndim,ndim)
      dimension b(ndim,ndim),b1(ndim,ndim)
c
      b1(1:ndi1,1:ndi2) = matmul(u(1:ndi1,1:ndi1),b(1:ndi1,1:ndi2))
      b1(1:ndi1,1:ndi2) = matmul(b1(1:ndi1,1:ndi2),ust(1:ndi2,1:ndi2))
c
      return
      end
      subroutine tripmt2(u,b,v,b1,ndi1,ndi2,ndi3,ndi4,ndim,mdim)
c =====================
c
c Triple product of non-rectangular matrices
c
c b1= u * b * v
c
c u:  ndi1 x ndi2
c b:  ndi2 x ndi3
c v:  ndi3 x ndi4
c b1: ndi1 x ndi4
c
      implicit real*8 (a-h,o-z)
      complex*16 u,v,b,b1,temp
      dimension u(ndim,ndim),v(ndim,ndim)
      dimension b(ndim,ndim),b1(mdim,mdim)
      dimension temp(ndim,ndim)
c
      temp(1:ndi1,1:ndi3) = matmul(u(1:ndi1,1:ndi2),b(1:ndi2,1:ndi3))
      b1(1:ndi1,1:ndi4) = matmul(temp(1:ndi1,1:ndi3),v(1:ndi3,1:ndi4))
c
      return
      end
      subroutine outmat(mat,n,m,ndim,nper)
c=====================
      implicit real*8 (a-h,o-z)
      parameter(nndim=50)
c
      complex*16 mat(ndim,ndim)
      complex*16 mat1(nndim,nndim)
c
      if((nndim.lt.n).or.(nndim.lt.m))
     > stop 'Increase nndim in outmat'
c
      do i=1,n
       do j=1,m
        r1=dreal(mat(i,j))
        r2=dimag(mat(i,j))
        if(dabs(r1).lt.1.d-15) r1=0.d0
        if(dabs(r2).lt.1.d-15) r2=0.d0
        mat1(i,j)=dcmplx(r1,r2)
       end do
      end do
c
      write(nper,*) ' real part'
      do 1 i=1,n
   1  write(nper,10) (dreal(mat1(i,j)),j=1,m)
      write(nper,*) ' imaginary part'
      do 2 i=1,n
   2  write(nper,10) (dimag(mat1(i,j)),j=1,m)
c
  10  format(9(1pd14.6))
c
      return
      end
      subroutine outmat1(mat,n,m,ndim,tol,nper)
c=======================
      implicit real*8 (a-h,o-z)
      complex*16 mat(ndim,ndim)
c
      do j=1,m
      do 1 i=1,n
        r1=dreal(mat(i,j))
        r2=dimag(mat(i,j))
        if(dabs(r1).lt.tol.and.dabs(r2).lt.tol) goto 1 
        write(nper,'(2i4,5x,1pd20.10,1pd20.10)') i,j,r1,r2
   1  continue
      end do
c
      return
      end
c
c
      subroutine infinitynorm(norm,a,ndim,ndimp)
c=======================
c L infinity matrix norm
c
      integer*4,  intent(in)  :: ndim, ndimp
      complex*16, intent(in)  :: a(ndimp,ndimp)
      real*8,     intent(out) :: norm
      integer*4 :: i, j
c
      norm = -1.d0
      do j=1,ndim
        do i=1,ndim
          norm = max(norm,abs(a(i,j)))
        end do
      end do
c   Done!
      end subroutine infinitynorm
c
c
      subroutine frobeniusnorm(norm,a,ndim,ndimp)
c=======================
c Frobenius matrix norm
c
      integer*4,  intent(in)  :: ndim, ndimp
      complex*16, intent(in)  :: a(ndimp,ndimp)
      real*8,     intent(out) :: norm
      integer*4 :: i, j
c
      norm = -1.d0
      do j=1,ndim
        do i=1,ndim
          norm = max(norm,abs(a(i,j))**2)
        end do
      end do
      norm = sqrt(norm)
c   Done!
      end subroutine frobeniusnorm
c
