c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine findangles(x,y,z,theta,phi)
c     =====================
c Find theta and phi angles for a unit vector (x,y,z)
c
      implicit real*8(a-h,o-z)
c
      data tol/1.0d-08/
      pi=dacos(-1.0d0)
c
      rxy=dsqrt(x*x+y*y)
      r  =dsqrt(x*x+y*y+z*z)
      if(r.lt.tol) then
         theta=0.0d0
      else
         theta=dacos(z/r)
      endif
      if(rxy.lt.tol) then
         phi=0.d0
      else
         phi=dacos(x/rxy)
         if(y.lt.0.d0) phi=2.0d0*pi-phi 
      endif
c
      return
      end
