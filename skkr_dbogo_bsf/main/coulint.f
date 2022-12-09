c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine coulint(rs,ns,dx,lam,f,g)
c=======================
c
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
      dimension f(nrad),r(nrad),g(nrad)
      dimension f1(nrad),f2(nrad),g1(nrad),g2(nrad)
c
      xs=dlog(rs)
      x=xs-(ns-1)*dx
      do i=1,ns
        r(i)=dexp(x)
        rl=r(i)**lam
        f1(i)=f(i)*rl
        f2(i)=f(i)/(rl*r(i))
        x=x+dx
      end do
c
      do i=1,ns
        g1(i)=rsimp(f1,r,i,dx)
        g2(i)=rsimp(f2,r,i,dx)
      end do
c
      do i=1,ns
        rl=r(i)**lam
        g(i)=g1(i)/(rl*r(i))+rl*(g2(ns)-g2(i))
      end do
c
      return
      end
