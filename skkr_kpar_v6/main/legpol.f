c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine legpol(x,l,pl,plm1,dpl)
c======================
c
c calculate legendre polynomials of degree l and l-1
c and derivative of degree l
c
      implicit real*8 (a-h,o-z)
c
      pnm1=1.d0
      pn  =x
c
      do 1 n=1,l-1
      pnp1=((n+n+1)*x*pn-n*pnm1)/(n+1.d0)
      pnm1=pn
    1 pn  =pnp1
c
      pl  =pn
      plm1=pnm1
      dpl =l*(x*pl-plm1)/(x*x-1.d0)
      return
      end
