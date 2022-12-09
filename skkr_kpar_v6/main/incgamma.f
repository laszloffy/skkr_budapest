c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
!info
!info incomplete gamma function gam(1/2-n,z)
!info                                        by A.Ernst 25.01.98 Vienna
! also for n < 0, array gam(1/2,z),gam(3/2,z),...  Robert Hammerling, 3.7.2002

      subroutine incgamma(gam,n,z)
!=============================================================================
      implicit real*8 (a-h,o-z)
!i/o
      include '../param.h'
      integer n
      complex*16 z, gam(0:lshape+1)
!local
      integer i
      double precision pi,b
      parameter (pi=3.1415926535897932384626d0) 
      complex*16 a,ci,ex,zb,sz,s 
      parameter (ci=(0.d0,1.d0))
!external
      complex*16 cerf
      external cerf
!==============================================
      ex=exp(-z)
      sz=sqrt(z)
      a=sz*exp(ci*pi/2.d0)
      s=sqrt(pi)*ex*cerf(a)
! literature: gam(1/2,z)  
      gam(0)=s
      if (n .lt. 0) then
         do i=1,-n,1
            gam(i)=z**(dfloat(i)-0.5d0)*ex + (dfloat(i)-0.5d0)*gam(i-1)
         end do
      else
         b=1.d0/2.d0-1.d0
         zb=sz
         i=0
 1       i=i+1
         if(i.gt.n) goto 10
         zb=zb/z
         s=(s-zb*ex)/b
         gam(i)=s
         b=b-1.d0
         goto 1
 10      return
      end if

      end
