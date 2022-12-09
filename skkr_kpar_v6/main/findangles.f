c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine findangles(x,y,z,theta,phi)
c     =====================
c Find theta and phi angles for a unit vector (x,y,z)
c
      implicit real*8(a-h,o-z)
c
      data tol/1.0d-05/
      pi=dacos(-1.0d0)
c
      rxy=dsqrt(x*x+y*y)
      r  =dsqrt(x*x+y*y+z*z)
      dx =dabs(x)
      dy =dabs(y)
      dz =dabs(z)
c
      if(r.lt.tol) then
         theta=0.0d0
         phi=0.0d0
         return
      end if
c
      if(dz.gt.tol) then
         theta=dacos(dz/r)
         if(z.lt.0.d0) theta=pi-theta
      else
         theta=0.5d0*pi
      endif
c
      if(dx.lt.tol.and.dy.lt.tol) then
         phi=0.d0
      else
         if(dy.lt.tol) then
           if(x.gt.0.d0) phi=0.d0
           if(x.lt.0.d0) phi=pi
           return
         end if
         if(dx.lt.tol) then
           if(y.gt.0.d0) phi=0.5d0*pi
           if(y.lt.0.d0) phi=1.5d0*pi
           return
         end if
         phi=dacos(dx/rxy)
         if(x.lt.0.d0.and.y.gt.0.d0) phi=pi-phi 
         if(x.lt.0.d0.and.y.lt.0.d0) phi=pi+phi 
         if(x.gt.0.d0.and.y.lt.0.d0) phi=2.0d0*pi-phi 
      endif
c
      return
      end
