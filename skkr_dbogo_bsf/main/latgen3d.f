c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
C***latgen.spg  processed by SPAG 5.11R  at 12:38 on 17 Feb 2000
      subroutine latgen3d(n1,n2,n3,m,r,mmax,dist,a1,a2,a3,overflow)
c --------------------------------------------------------------------**
c generate lattice vectors with basis a1, a2, a3
c ------------------------- dummy arguments --------------------------**
c                                                                     **
c  input  - n1                                                        **
c  input  - n2                                                        **
c  input  - n3                                                        **
c  input  - m                                                         **
c  output - r                                                         **
c  output - mmax                                                       **
c  input  - dist                                                      **
c  input  - a1                                                        **
c  input  - a2                                                        **
c  input  - a3                                                        **
c  output - overflow
c                                                                     **
c ------------------------- common variables -------------------------**
c                                                                     **
c  modifies    ** nothing **                                          **
c  uses value  ** nothing **                                          **
c                                                                     **
c ----------------------- external subprograms -----------------------**
c                                                                     **
c  calls       ** nothing **                                          **
c  called by   vector                                                 **
c                                                                     **
c --------------------------------------------------------------------**
c
      implicit none
C*** Start of declarations inserted by SPAG
      real*8 a1,a2,a3,dist,dist2,r,r2,x,y,z
      integer i,j,k,m,mmax,n1,n2,n3,nr
      logical overflow
C*** End of declarations inserted by SPAG
      dimension r(3,m),a1(3),a2(3),a3(3)
c
      overflow=.false.
      nr=0
      dist2=dist*dist
      do i=-n1,n1
        do j=-n2,n2
          do k=-n3,n3
            x=a1(1)*i+a2(1)*j+a3(1)*k
            y=a1(2)*i+a2(2)*j+a3(2)*k
            z=a1(3)*i+a2(3)*j+a3(3)*k
            r2=x*x+y*y+z*z
            if (r2.le.dist2) then
              nr=nr+1
              if (nr.gt.m) then
                 overflow=.true.
                 return
              endif
              r(1,nr)=x
              r(2,nr)=y
              r(3,nr)=z
            endif
          enddo
        enddo
      enddo
c
      mmax=nr
      return
      end
