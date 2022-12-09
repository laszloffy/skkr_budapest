c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine spherh(lmax,x,y,z,cylm)
c======================
c
c calculate spherical harmonics with condon-shortley phase convention
c as given by altmann: quaternions; subroutine init has to be called
c first to calculate normalization factors ylnorm.
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      parameter(small=1.d-12,nmax=lshape+30)
c
      dimension ylm(2,-lshape:lshape,0:lshape)
      complex*16 cylm(-lshape:lshape,0:lshape)
c
      common/ssphnrm/ylnorm(0:lshape,0:lshape),al(0:nmax),
     >              alinv(nmax),tlp1(-nmax:nmax)
c
      rxy=dsqrt(x*x+y*y)
      r  =dsqrt(x*x+y*y+z*z)
      if(r.lt.small) then
      costh=1.d0
      sinth=0.d0
         else
      rinv=1.d0/r
      costh=  z*rinv
      sinth=rxy*rinv
         endif
      if(rxy.lt.small) then
      cosfi=1.d0
      sinfi=0.d0
         else
      rinv=1.d0/rxy
      cosfi=x*rinv
      sinfi=y*rinv
         endif
c
      cosmfi=1.d0
      sinmfi=0.d0
      plm=1.d0
c
      do 2 m=0,lmax-1
      sig=1.d0-2.d0*mod(m,2)
      plp1m  = tlp1(m)*plm*costh
      plp1mp1=-tlp1(m)*plm*sinth
      yr=plm*ylnorm(m,m)*cosmfi
      yi=plm*ylnorm(m,m)*sinmfi
c
      ylm(1,-m,m)= yr*sig
      ylm(2,-m,m)=-yi*sig
      ylm(1, m,m)= yr
      ylm(2, m,m)= yi
c
      yr=plp1m*ylnorm(m,m+1)*cosmfi
      yi=plp1m*ylnorm(m,m+1)*sinmfi
c
      ylm(1,-m,m+1)= yr*sig
      ylm(2,-m,m+1)=-yi*sig
      ylm(1, m,m+1)= yr
      ylm(2, m,m+1)= yi
c
      do 1 l=m+1,lmax-1
      plm1m=plm
      plm=plp1m
      plp1m=(tlp1(l)*plm*costh-al(l+m)*plm1m)*alinv(l-m+1)
      yr=plp1m*ylnorm(m,l+1)*cosmfi
      yi=plp1m*ylnorm(m,l+1)*sinmfi
c
      ylm(1,-m,l+1)= yr*sig
      ylm(2,-m,l+1)=-yi*sig
      ylm(1, m,l+1)= yr
      ylm(2, m,l+1)= yi
c
    1 continue
c
      plm=plp1mp1
      r=cosmfi*cosfi-sinmfi*sinfi
      sinmfi=sinmfi*cosfi+cosmfi*sinfi
      cosmfi=r
    2 continue
c
      yr=plm*ylnorm(lmax,lmax)*cosmfi
      yi=plm*ylnorm(lmax,lmax)*sinmfi
c
      sig=1.d0-2.d0*mod(lmax,2)
      ylm(1,-lmax,lmax)= yr*sig
      ylm(2,-lmax,lmax)=-yi*sig
      ylm(1, lmax,lmax)= yr
      ylm(2, lmax,lmax)= yi
c
      do l=0,lmax
      do m=-l,l
        cylm(m,l)=dcmplx(ylm(1,m,l),ylm(2,m,l))
      end do
      end do
c
      return
      end
