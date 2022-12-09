c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
!info
!info calculate delta(n) for Kambe structure constant
!info 
      subroutine ilkam(delta,l,alpha,t,kappa)
      implicit none
!i/o
      integer l
      double precision t
      complex*16 alpha,delta(0:l),kappa
!local
      integer i
      double precision pi,dc,d
      parameter (pi=3.1415926535897932384626d0) 
      complex*16 a,b,sx,c,ci,ex,am,ap,wm,wp,del,delm1,del2,cx,z,x
      parameter (ci=(0.d0,1.d0))
!external
      complex*16 cw,cerf
      external cw,cerf
!======================================================================
      z=kappa*t
      x=alpha
      sx=sqrt(x)
!==================
      b=2.d0/z
      dc=dsqrt(pi)/2.d0
      a=z/2.d0/sx
      ex=exp(x-a*a)
      c=ci*a
      ap= sx+c
      am=-sx+c
      wp=cerf(ap)
      wm=cerf(am)
      c=ex*dc
      del=wp+wm
      del=del*c*b
      delm1=wp-wm
      delm1=delm1*c/ci
      delta(0)=del
!===================== - beginning of recursion - ====================
      b=b*b
      cx=sx
      d=0.5d0
      i=0
 1    i=i+1
      if(i.gt.l) goto 10
      cx=cx/x
      del2=d*del-delm1+cx*ex
      del2=del2*b
      delta(i)=del2
      delm1=del
      del=del2
      d=d+1.d0
      goto 1
 10   return
      end
