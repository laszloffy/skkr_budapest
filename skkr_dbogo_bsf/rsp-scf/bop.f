c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine bop(rs,ns,dx,rhodsp,lz,bopr)
c=======================
c
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
      real*8 lz(2)
      dimension rhodsp(nrad,2),bopr(nrad,2)
      dimension help(nrad),r(nrad),b1(nrad),b2(nrad)
c
      pi=dacos(-1.0d0)
c     fac=8.0d0*pi/441.0d0
      fac=2.0d0/441.0d0
c
      xs=dlog(rs)
      do is=1,2
c
        x=xs-(ns-1)*dx
        do i=1,ns
          r(i)=dexp(x)
          help(i)=rhodsp(i,is)
          x=x+dx
        end do
        rnorm=rsimp(help,r,ns,dx)
        do i=1,ns
          help(i)=help(i)/rnorm
        end do
c
        call coulint(rs,ns,dx,2,help,b1)
        call coulint(rs,ns,dx,4,help,b2)
        do i=1,ns
          bopr(i,is)=-r(i)*fac*(9.0d0*b1(i)-5.0d0*b2(i))*lz(is)
        end do
c
      end do
c
      return
      end
