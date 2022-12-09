c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine latgen2d(n1,n2,m,r,mmax,dist,a1,a2,overflow)
c========================
c
c generate lattice vectors in 2D  with basis a1, a2
c
      implicit real*8 (a-h,o-z)
      logical overflow
      dimension r(2,m),a1(2),a2(2)
c
      overflow=.false.
      nr=0
      dist2=dist*dist
      do i=-n1,n1
      do j=-n2,n2
        x=a1(1)*i+a2(1)*j
        y=a1(2)*i+a2(2)*j
        r2=x*x+y*y
        if(r2.le.dist2) then      
          nr=nr+1
          if(nr.gt.m) then
            overflow=.true.
            return
          end if
          r(1,nr)=x
          r(2,nr)=y
        end if
      end do    
      end do    
c
      mmax=nr
      return
      end
