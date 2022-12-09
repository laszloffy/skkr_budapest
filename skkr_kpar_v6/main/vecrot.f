c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine vecrot(w,phi)
c     Rotate real vector by an angle of phi
c
      implicit real*8(a-h,o-z)
      dimension w(2),v(2),a(2,2)
c
      a(1,1) = dcos(phi)
      a(1,2) =-dsin(phi)
      a(2,1) = dsin(phi)
      a(2,2) = dcos(phi) 
c
      do i=1,2
        v(i) = a(i,1)*w(1) + a(i,2)*w(2)
      enddo
c
      do i=1,2
        w(i) = v(i)
      enddo
c
      return
      end
