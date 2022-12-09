c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
!info
!info calculate delta(n) for Kambe structure constant
!info 
      subroutine delkam(delta,n,x,zz)
      implicit none
!i/o
      integer n
      complex*16 delta(0:n),x,z,zz
!local
      integer i
      double precision pi,dc,d
      parameter (pi=3.1415926535897932384626d0) 
      complex*16 a,b,sx,c,ci,ex,am,ap,wm,wp,del0,del1,del2,cx,fc
      parameter (ci=(0.d0,1.d0))
!external
      complex*16 cw,cerf,zn,sxn
      external cw,cerf
!======================================================================
c     write(6,*)'input:n,x,zz',n,x,zz
      if(dreal(zz).lt.0.d0) then
         zn=ci*sqrt(abs(zz))
         sxn=sqrt(x)
      else
         sxn=-ci*sqrt(abs(x))
         zn=sqrt(zz)
      endif
      z=sqrt(zz)
      sx=sqrt(x)
      dc=dsqrt(pi)/2.d0
      b=2.d0/z
      a=z/2.d0/sx
      ex=exp(a*a-x)
      c=ci*sx
      am=-a+c
      ap= a+c
c     write(6,*)'am,ap',am,ap
      wm=cerf(am)
      wp=cerf(ap)
      c=ex*dc
      del0=wm+wp
      del0=del0*c
      del1=wm-wp
      del1=del1*c*b*ci
      delta(0)=del0
      delta(1)=del1
!===================== - beginning of recursion - =====================
      fc=sx
      b=b*b
      d=0.5d0
      do i=0,n-2
         fc=fc/x
         delta(i+2)=b*(-d*delta(i+1)-delta(i)+ex*fc)
         d=d+1.0d0
      enddo
      return
      end
