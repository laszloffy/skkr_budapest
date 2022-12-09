c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine matr(tlm,lmax,dmat,dmat1)
c ===================
c
c R=(tvec,phi)
c dmat = D(R) and dmat1 = D(R)+ 
c
      implicit real*8 (a-h,o-z)
      include '../param.h'
      parameter(jdim=2*lsup+2)
c
      dimension tlm(4)
      complex*16 dmat(kmymaxp,kmymaxp),dmat1(kmymaxp,kmymaxp)
      complex*16 cmat(kmymaxp,kmymaxp)
      complex*16 d(jdim,jdim)
      common/test/itest
c
      kmax=2*lmax+1
      kmymax=2*(lmax+1)*(lmax+1)
c
      call czero(dmat,kmymaxp*kmymaxp)
      call czero(dmat1,kmymaxp*kmymaxp)
c
c Set up matrix of rotation
c
      ist=0
      do j2=1,2*lmax-1,2
        il=j2+1
        call rotmat(d,trd,il,tlm,1)
        do icase=1,2
          do m1=1,il
          do m2=1,il
c fill up matrix according to ascending m index
            dmat(ist+m1,ist+m2)=d(il+1-m1,il+1-m2)
          end do
          end do
          ist=ist+il
        end do
      end do
      j2=2*lmax+1
      il=j2+1
      call rotmat(d,trd,il,tlm,1)
      do m1=1,il
      do m2=1,il
c fill up matrix according to ascending m index
        dmat(ist+m1,ist+m2)=d(il+1-m1,il+1-m2)
      end do
      end do
c
      do i=1,kmymax
      do j=1,kmymax
        dmat1(i,j)=dconjg(dmat(j,i))
      end do
      end do
c
      if(itest.lt.4) return
      write(6,'(/'' Matrix of rotation'')') 
      call outmat(dmat,kmymax,kmymax,kmymaxp,6)
      write(6,'(/'' Inverse'')') 
      call outmat(dmat1,kmymax,kmymax,kmymaxp,6)
      call repl(cmat,dmat,kmymax,kmymaxp)
      call doubmt(cmat,dmat1,kmymax,kmymaxp)
      write(6,'(/'' D(R) * D(R**-1) '')')
      call outmat(cmat,kmymax,kmymax,kmymaxp,6)
c
      return
      end
