c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine initzero1(qvpa,qva,vmadid,vmadiq,
     & spin_magvpa,spin_magva,
     & orb_magvpa,orb_magva,
     & enba,denba,enorba,qmoma,
     & rhova,rhospa,rhodspa,rhomaga,
     & enbifc,qvifc,omifc,enpba)
c
      implicit real*8 (a-h,o-z) 
      include '../param.h'
c
      integer mmax
      parameter (mmax=mimp)
c
      complex*16 qmoma(lmsup,mmax)
c
      real*8 qva(mmax)
      real*8 qvpa(kmymaxp,mmax)
      real*8 spin_magvpa(kmymaxp,mmax,3)
      real*8 spin_magva(mmax,3)
      real*8 orb_magvpa(kmymaxp,mmax,3)
      real*8 orb_magva(mmax,3)
      real*8 enba(mmax)
      real*8 enpba(kmymaxp,mmax)
      real*8 denba(mmax)
      real*8 enorba(mmax)
      real*8 rhova(nrad,mmax)
      real*8 rhospa(nrad,2,mmax)
      real*8 rhodspa(nrad,2,mmax)
      real*8 rhomaga(nrad,mmax)
c
      real*8 vmadid(mmax)
      real*8 vmadiq(mmax)
c
      real*8 enbifc
      real*8 qvifc
      real*8 omifc
c
      call rzero(qvpa,kmymaxp*mmax)
      call rzero(qva,mmax)
c
      call rzero(spin_magvpa,3*kmymaxp*mmax)
      call rzero(spin_magva,3*mmax)
c
      call rzero(orb_magvpa,3*kmymaxp*mmax)
      call rzero(orb_magva,3*mmax)
c
      call rzero(enba,mmax)
      call rzero(enpba,kmymaxp*mmax)
      call rzero(denba,mmax)
      call rzero(enorba,mmax)
c
      call rzero(rhova,nrad*mmax)
      call rzero(rhospa,2*nrad*mmax)
      call rzero(rhodspa,2*nrad*mmax)
      call rzero(rhomaga,nrad*mmax)
c
      call rzero(vmadid,mmax)
      call rzero(vmadiq,mmax)
c
      call czero(qmoma,lmsup*mmax)
c
      enbifc=0.d0
      qvifc=0.d0
      omifc=0.d0
c
      return
      end
