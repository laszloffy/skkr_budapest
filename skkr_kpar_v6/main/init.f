c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine init
c====================
c
c normalization of spherical harmonics
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      parameter(nmax=lshape+30)
      parameter(small=1.d-14)
c
      complex*16 ylm(-lshape:lshape,0:lshape)
c
      common/cgaunt/yg(lshape,lmshape)
      common/ssphnrm/ylnorm(0:lshape,0:lshape),al(0:nmax),
     >               alinv(nmax),tlp1(-nmax:nmax)
c
      pi=4.0d0*datan(1.0d0)
      do 1 l=0,lshape
    1 ylnorm(0,l)=dsqrt((2*l+1)/(4.d0*pi))
c
      do 2 m=1,lshape
      do 2 l=m,lshape
    2 ylnorm(m,l)=ylnorm(m-1,l)/dsqrt((l-m+1.d0)*(l+m))
c
      al(0)=0.d0
      do 6 l=1,nmax
      al(l)=l
    6 alinv(l)=1.d0/al(l)
c
      do 7 l=-nmax,nmax
    7 tlp1(l)=2*l+1
c
c associated legendre polynomials and corresponding
c weights to be used in function gaunt
c
      ng=2*lshape
c
      do 5 i=1,lshape
      phi=(4*i-1)*pi/(4*ng+2.d0)
      z=dcos(phi+1.d0/(8*ng*ng*dtan(phi)))
      ind=0
c
c find zeros of legendre polynomials to be used in gaussian integration:
c newton method with starting values according to abramowitz & stegun
c
    3 call legpol(z,ng,pl,plm1,dpl)
      ind=ind+1
      if(ind.gt.100) stop 'init:1'
      dz=pl/dpl
      z=z-dz
      if(dabs(dz).gt.small) goto 3
      x=1.d0-z*z
      w=2.d0*x/(ng*ng*plm1*plm1)
c
      fac=(4.d0*pi*w)**(1.d0/3.d0)
      y=0.d0
      x=dsqrt(x)
      call spherh(lshape,x,y,z,ylm)
      ind=0
      do 4 l=0,lshape
      do 4 m=-l,l
      ind=ind+1
    4 yg(i,ind)=dreal(ylm(m,l))*fac
    5 continue
c
      return
      end
