c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine matnr(tlm,lmax,dmat,dmat1)
c ====================
c
c R=(tvec,phi)
c dmat = D(R) and dmat1 = D(R)+ 
c
      implicit real*8 (a-h,o-z)
      include '../param.h'
      parameter(jdim=2*lsup+2)
c
      dimension tlm(4)
      complex*16 dmat(lmsup,lmsup),dmat1(lmsup,lmsup)
      complex*16 cmat(lmsup,lmsup)
      complex*16 d(jdim,jdim)
      common/test/itest
c
      lmaxs=2*lmax
      lmmaxs=(lmaxs+1)*(lmaxs+1)
c
c Set up matrix of rotation
c
      do i=1,lmmaxs
      do j=1,lmmaxs
        dmat(i,j)=(0.d0,0.d0)
      end do
      end do
      ist=0
      do l=0,lmaxs
        il=2*l+1
        call rotmat(d,trd,il,tlm,1)
        do m1=1,il
        do m2=1,il
c fill up matrix according to ascending m index
          dmat(ist+m1,ist+m2)=d(il+1-m1,il+1-m2)
        end do
        end do
        ist=ist+il
      end do
c
      do i=1,lmmaxs
      do j=1,lmmaxs
        dmat1(i,j)=dconjg(dmat(j,i))
      end do
      end do
c
      if(itest.lt.4) return
      write(6,'(/'' Matrix of rotation'')') 
      call outmat(dmat,lmmaxs,lmmaxs,lmsup,6)
      write(6,'(/'' Inverse'')') 
      call outmat(dmat1,lmmaxs,lmmaxs,lmsup,6)
      call repl(cmat,dmat,lmmaxs,lmsup)
      call doubmt(cmat,dmat1,lmmaxs,lmsup)
      write(6,'(/'' D(R) * D(R**-1) '')')
      call outmat(cmat,lmmaxs,lmmaxs,lmsup,6)
c
      return
      end
