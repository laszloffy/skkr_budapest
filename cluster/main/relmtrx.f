c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine relmtrx(a,b,lmax)
c=======================
c
c *****************************************************
c * transformation from a non-relativistic matrix 'a' *
c *                to a relativistic matrix 'b'       *
c *****************************************************
c
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
      complex*16 a,b
      complex*16 a1(36,36)
c
      dimension a(lmmaxp,lmmaxp),b(kmymaxp,kmymaxp)
      dimension u1(50),ind1(50)
      dimension u2(50),ind2(50)
      common/cgc/u1,u2,ind1,ind2
c
      complex*16 fac
      common/relfac/fac
c
      lmmax=(lmax+1)*(lmax+1)
      kmymax=2*lmmax
      if(lmmax.gt.36) stop 'Increase dimension of a1 in relmtrx!'
c
c     call czero(a1,26*26)
      a1=(0.d0,0.d0)
      do i=1,lmmax
      do j=1,lmmax
        a1(i,j)=fac*a(i,j)
      end do
      end do
c
      do i=1,kmymax
        i1=ind1(i)
        i2=ind2(i)
        do j=1,kmymax
          j1=ind1(j)
          j2=ind2(j)
          b(i,j)=u1(i)*a1(i1,j1)*u1(j)+
     >           u2(i)*a1(i2,j2)*u2(j)
        end do
      end do
c
      return
      end
