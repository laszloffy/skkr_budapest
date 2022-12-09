c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine ytrafo(u,up,lmax)
c
c ********************************************************
c *                                                      *
c * transformation matrix from real spherical harmonics  *
c * to complex spherical harmonics                       *
c *                                                      *
c * YC(L) = sum_{L'} YR(L') U(L',L)                      *
c * YC(L)* = sum_{L'} UP(L,L') YR(L')                    *
c *                                                      *
c ********************************************************
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      complex*16 u(lmmaxp,lmmaxp),up(lmmaxp,lmmaxp)
      complex*16 sum
c
c ========================================================
c
      lmmax=(lmax+1)*(lmax+1)
      fac=1.0d0/dsqrt(2.0d0)
      call czero(u,lmmaxp*lmmaxp)
      call czero(up,lmmaxp*lmmaxp)
c
      do l=0,lmax
      do m=-l,l
        lm=l*l+l+m+1
        lmp=l*l+l-m+1
c
        if(m.lt.0) then
          u(lm,lm)=dcmplx(0.0d0,-fac)
          u(lmp,lm)=fac
        end if
        if(m.eq.0) u(lm,lm)=1.0d0
        if(m.gt.0) then
          c=1.0d0
          n=mod(m,2)
          if(n.eq.1) c=-1.0d0
          u(lm,lm)=c*fac
          u(lmp,lm)=dcmplx(0.0d0,c*fac)
        end if
c
      end do
      end do
c
      do lm=1,lmmax
      do lmp=1,lmmax
        up(lm,lmp)=dconjg(u(lmp,lm))
      end do
      end do
c
      return
      end
