c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine legzero(ng,xx,wx)
c=======================
c
c find zero's of legendre polinomial of order n and corresponding
c weighting factors to gaussian quadrature
c newton method with starting values according to abramowitz & stegun
c
      implicit real*8 (a-h,o-z)
c
      parameter(pi=3.1415926535897932d0,small=1.d-14)
c
      dimension xx(*),wx(*)
c
      n=ng/2
      do i=1,n
        phi=(4*i-1)*pi/(4*ng+2.d0)
        x=dcos(phi+1.d0/(8*ng*ng*dtan(phi)))
        ind=0
    1   call legpol(x,ng,pl,plm1,dpl)
        ind=ind+1
        if(ind.gt.100) stop 'legzero'
        dx=pl/dpl
        x=x-dx
        if(dabs(dx).gt.small) goto 1
        w=2.d0*(1.d0-x*x)/(ng*ng*plm1*plm1)
        xx(ng+1-i)=x
        xx(i)=-x
        wx(ng+1-i)=w
        wx(i)=w
      end do
c
      if(2*n.ne.ng) then
        x=0.d0
        call legpol(x,ng,pl,plm1,dpl)
        w=2.d0/(ng*ng*plm1*plm1)
        xx(n+1)=x
        wx(n+1)=w
      end if
c
      return
      end
