c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine rzero(f,n)
      real*8 f(*)
      do i=1,n
        f(i)=0.d0
      end do
      return
      end
      subroutine czero(f,n)
      complex*16 f(*)
      do i=1,n
        f(i)=(0.d0,0.0d0)
      end do
      return
      end
