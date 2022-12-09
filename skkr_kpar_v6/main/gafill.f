c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine gafill
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      complex*16 gacoeff(lmmaxp,lmmaxp,lmsup)
      complex*16 rgacoeff(kmymaxp,kmymaxp,lmsup)
c
      common/ggaunt/gacoeff
      common/rggaunt/rgacoeff
c
      sqfpi=dsqrt(16.d0*datan(1.d0))
      i=0
      do lam=0,lsup
      fac=sqfpi/(2*lam+1)
      do nu=-lam,lam
        i=i+1
c
        call czero(gacoeff(1,1,i),lmmaxp*lmmaxp)
        call czero(rgacoeff(1,1,i),kmymaxp*kmymaxp)
c
        do l=0,lmaxp
        do m=-l,l
          lm=l*l+l+m+1
          mp=m-nu
          lpmin=iabs(l-lam)
          lpmax=l+lam
          do 1 lp=lpmin,lpmax,2
            if(lp.gt.lmaxp) goto 1
            if(iabs(mp).gt.lp) goto 1
            lmp=lp*lp+lp+mp+1
            gacoeff(lm,lmp,i)=fac*gaunt(l,lam,lp,m,nu,mp)
    1     continue
        end do
        end do
c       call outmat(gacoeff(1,1,i),lmmaxp,lmmaxp,lmmaxp,6)
c       call relmtrx(gacoeff(1,1,i),rgacoeff(1,1,i),lmaxp)
        call relmtrx(gacoeff(1,1,i),rgacoeff(1,1,i),(1.d0,0.d0),lmaxp) ! 01/03/2021 Laszloffy
c       call outmat(rgacoeff(1,1,i),kmymaxp,kmymaxp,kmymaxp,6)
c
      end do
      end do
c
      return
      end
