c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine setphase(z0,z1)
      complex*16 z0,z1
      real*8     r0,r1,c0,c1
      integer    iround
      real*8     x,pi
c
      pi = 4.d0*datan(1.d0)
c
      r0 = dreal(z0)      
      r1 = dreal(z1)      
      c0 = dimag(z0)
      c1 = dimag(z1)
c
      x = (c1 - c0)/pi
      c1 = c1 - iround(x)*pi
      z1 = dcmplx(r1,c1)
      return
      end
      integer function iround(x)
      real*8 x
c
      ix = int(x)      
      if(x.gt.0.d0) then
       if(x-ix.gt.0.5d0) ix = ix + 1
      else
       if(dabs(x-ix).gt.0.5d0) ix = ix - 1
      endif 
      iround = ix
      return
      end
