c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine dens(tau,y1,y2,rho,imax,lmax,para)
c====================
c
c calculate partial differential charge distribution
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      logical para
c
      complex*16 tau(lmmaxp,lmmaxp),sum,tsum
      complex*16 y1(nrad,0:lmaxp),y2(nrad,0:lmaxp),rho(nrad)
      data pi/3.1415926535897932d0/
c
      fac=-1.d0/pi
      if(para) fac=2.0d0*fac
      do i=1,imax
        sum=(0.d0,0.d0)
        lm=0
        do l=0,lmax
          tsum=(0.0d0,0.0d0)
          tlp1=dfloat(2*l+1)
          do m=-l,l
             lm=lm+1
             tsum=tsum+tau(lm,lm)
          end do
          sum=sum+y1(i,l)*y1(i,l)*tsum-tlp1*y1(i,l)*y2(i,l)
        end do
        rho(i)=fac*sum
      end do
c
      return
      end
