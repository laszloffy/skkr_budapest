c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine tmat(ce,lmax,idpot,v0,vr,dx,ns,rs,rel,tminv,alphakkr)
c====================
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      logical rel
      character*10 idpot
      dimension vr(nrad)
      dimension lml(lmmaxp),lmm(lmmaxp)
      complex*16 ce,detl,alphakkr(0:lmaxp)
      complex*16 tminv(lmmaxp,lmmaxp),tml(0:lmaxp)
      complex*16 fint(0:lmaxp),gint(0:lmaxp),fpint(lmaxp)
      complex*16 yl1(nrad,0:lmaxp),yl2(nrad,0:lmaxp)
c
      common/test/itest
      common/mome/lml,lmm
c
      nl2=(lmax+1)*(lmax+1)
c     -------------------------------------------------------
      call czero(tminv,lmmaxp*lmmaxp)
      call wafu(ce,lmax,idpot,v0,vr,dx,ns,rs,0,rel,
     >          tml,fint,gint,fpint,yl1,yl2)
c     -------------------------------------------------------
c
c - Apply screening-transformation on tminv.
c   Supposed that tminv is diagonal!!!!!!!.
c
      do lm=1,nl2
         l=lml(lm)
         detl=tml(l)/(1.d0-tml(l)*alphakkr(l))
         tminv(lm,lm)=detl
      end do
c
      if(itest.lt.3) return
         write(6,'(''alphakkr:'',2x,10f12.8)') (alphakkr(l),l=0,lmax)
         write(6,'(2x,a,''      e='',2f10.6)') idpot,ce
         write(6,*) '  unscreened inverse t - matrix'
         write(6,'(8e13.5)') (tml(l),l=0,lmax)
         write(6,*) '  screened inverse t - matrix'
         write(6,'(8e13.5)') (tminv(l*l+1,l*l+1),l=0,lmax)
c
      return
      end
