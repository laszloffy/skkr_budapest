c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
!info
!info YLMC calculate complex spherical harmonics up to l=lmax
!info
      subroutine ylmc(x,y,z,lmax,ylm)
!=======================================================================
! Input:
!  x,y,z :
!  lmax  :
! Output:
!  ylm   : 
!=======================================================================
      implicit none
!input/output
      integer lmax
      double precision x,y,z
      complex*16 ylm(*)
!local
      integer i,k,l,lmx,m,mb,n,nit,lmmx
c     data nit/1/
      parameter (lmx=12)
      parameter (lmmx=(lmx+1)*(lmx+1))
      double precision b(0:lmx,0:lmx,0:lmx),fact(0:2*lmx),tiny,fa1,!pi4,
     &     fa2,r,cth,d,p(0:lmx),px,t,a(0:lmx,0:lmx),f(0:lmx,0:lmx)
      complex*16 c,ci,cf(lmmx),ephi
      double precision pi,pi4 ! 01/03/2021 Laszloffy
c     data ci/(0.d0,1.d0)/,pi4/12.56637061d0/,tiny/1.d-30/
      data ci/(0.d0,1.d0)/,tiny/1.d-30/
c
      pi=dacos(-1.0d0)
      pi4=4.d0*pi ! 01/03/2021 Laszloffy
c     if(nit.eq.1) then
         do l=0,lmx
            do n=0,lmx
               a(l,n)=0.d0
               f(l,n)=0.d0
               do k=0,lmx
                  b(l,n,k)=0.d0
               enddo ! k
            enddo ! n
         enddo ! l
         a(0,0)= 1.0d0
         a(1,1)= 1.0d0
         a(2,0)=-0.5d0
         a(2,2)= 1.5d0
         
         do l=3,lmx
            fa1= float(2*l-1)/float(l)
            fa2=-float(l-1)/float(l)
            a(l,0)=fa2*a(l-2,0)
            do n=1,l-2
               a(l,n)=fa1*a(l-1,n-1)+fa2*a(l-2,n)
            enddo ! n
            a(l,l-1)=fa1*a(l-1,l-2)
            a(l,l)=fa1*a(l-1,l-1)
         enddo ! l
         
         do n=0,lmx
            f(n,0)=1.d0
            do m=0,n-1
               f(n,m+1)=f(n,m)*float(n-m)
            enddo ! m
         enddo ! n
         do l=0,lmx
            do m=0,l
               do k=0,l-m
                  b(l,m,k)=a(l,k+m)*f(k+m,m)
               enddo ! k
            enddo ! m
         enddo ! l
! calculate factorials
         fact(0)=1.d0
         do n=1,2*lmx
            fact(n)=fact(n-1)*float(n)
         enddo ! n
! calculate coficients for spherical harmonics
         i=0
         do l=0,lmx
            do m=-l,l
               i=i+1
               mb=abs(m)
               cf(i)=sqrt(float(2*l+1)*fact(l-mb)/(pi4*fact(l+mb)))*  
     &              (ci**(m+mb))
            enddo
         enddo
c        nit=0
c     endif

      r=dsqrt(x*x+y*y+z*z)
      if(r.lt.1.d-16) then
         r=0.d0
         cth=0.d0
         ephi=cmplx(1.d0,0.d0)
      else
         cth=z/r
         d=dsqrt(x*x+y*y)
         if(d.lt.1.d-6) then
            ephi=cmplx(1.d0,0.d0)
         else
            ephi=(x+ci*y)/d
         endif
      endif
      px=sqrt(1.d0-cth*cth)+tiny

      p(0)=1.d0
      do k=1,lmax
         p(k)=p(k-1)*cth
      enddo ! k
      i=0
      do l=0,lmax
         do m=-l,l
            i=i+1
            mb=abs(m)
            c=cf(i)
            t=0.d0
            do k=0,l
               t=t+b(l,mb,k)*p(k)
            enddo ! k
            ylm(i)=c*(px**mb)*t*(ephi**m)
         enddo ! m
      enddo ! l
      return
      end
