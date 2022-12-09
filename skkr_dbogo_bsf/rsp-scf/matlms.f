c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine matlms(alms,akmy,lmax)
c======================
c
c akmy: (kappa,my) representation
c alms: (l,m,s,eh) representation
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      complex*16 alms(lmmaxp,lmmaxp,2,2,2,2),akmy(dbogomaxp,dbogomaxp)
      dimension u1(72),ind1(72)
      dimension u2(72),ind2(72)
      common/cgc/u1,u2,ind1,ind2
c
      lmmax=(lmax+1)*(lmax+1)
      kmymax=2*lmmax
c
      alms=(0.d0,0.d0)
      call czero(alms,16*lmmaxp*lmmaxp)
      do lm1=1,lmmax
      do lm2=1,lmmax
        do kmy=1,kmymax
        do kmyp=1,kmymax
c
c-------->electron-electron block
c
          if(ind1(kmy).eq.lm1.and.ind1(kmyp).eq.lm2) 
     >    alms(lm1,lm2,1,1,1,1)=alms(lm1,lm2,1,1,1,1)
     >                      +u1(kmy)*akmy(kmy,kmyp)*u1(kmyp)
          if(ind1(kmy).eq.lm1.and.ind2(kmyp).eq.lm2) 
     >    alms(lm1,lm2,1,2,1,1)=alms(lm1,lm2,1,2,1,1)
     >                      +u1(kmy)*akmy(kmy,kmyp)*u2(kmyp)
          if(ind2(kmy).eq.lm1.and.ind1(kmyp).eq.lm2) 
     >    alms(lm1,lm2,2,1,1,1)=alms(lm1,lm2,2,1,1,1)
     >                      +u2(kmy)*akmy(kmy,kmyp)*u1(kmyp)
          if(ind2(kmy).eq.lm1.and.ind2(kmyp).eq.lm2) 
     >    alms(lm1,lm2,2,2,1,1)=alms(lm1,lm2,2,2,1,1)
     >                      +u2(kmy)*akmy(kmy,kmyp)*u2(kmyp)
c
c-------->electron-hole block
c
          if(ind1(kmy).eq.lm1.and.ind1(kmyp).eq.lm2) 
     >    alms(lm1,lm2,1,1,1,2)=alms(lm1,lm2,1,1,1,2)
     >                      +u1(kmy)*akmy(kmy,kmyp+kmymax)*u1(kmyp)
          if(ind1(kmy).eq.lm1.and.ind2(kmyp).eq.lm2) 
     >    alms(lm1,lm2,1,2,1,2)=alms(lm1,lm2,1,2,1,2)
     >                      +u1(kmy)*akmy(kmy,kmyp+kmymax)*u2(kmyp)
          if(ind2(kmy).eq.lm1.and.ind1(kmyp).eq.lm2) 
     >    alms(lm1,lm2,2,1,1,2)=alms(lm1,lm2,2,1,1,2)
     >                      +u2(kmy)*akmy(kmy,kmyp+kmymax)*u1(kmyp)
          if(ind2(kmy).eq.lm1.and.ind2(kmyp).eq.lm2) 
     >    alms(lm1,lm2,2,2,1,2)=alms(lm1,lm2,2,2,1,2)
     >                      +u2(kmy)*akmy(kmy,kmyp+kmymax)*u2(kmyp)
c
c-------->hole-electron block
c
          if(ind1(kmy).eq.lm1.and.ind1(kmyp).eq.lm2) 
     >    alms(lm1,lm2,1,1,2,1)=alms(lm1,lm2,1,1,2,1)
     >                      +u1(kmy)*akmy(kmy+kmymax,kmyp)*u1(kmyp)
          if(ind1(kmy).eq.lm1.and.ind2(kmyp).eq.lm2) 
     >    alms(lm1,lm2,1,2,2,1)=alms(lm1,lm2,1,2,2,1)
     >                      +u1(kmy)*akmy(kmy+kmymax,kmyp)*u2(kmyp)
          if(ind2(kmy).eq.lm1.and.ind1(kmyp).eq.lm2) 
     >    alms(lm1,lm2,2,1,2,1)=alms(lm1,lm2,2,1,2,1)
     >                      +u2(kmy)*akmy(kmy+kmymax,kmyp)*u1(kmyp)
          if(ind2(kmy).eq.lm1.and.ind2(kmyp).eq.lm2) 
     >    alms(lm1,lm2,2,2,2,1)=alms(lm1,lm2,2,2,2,1)
     >                      +u2(kmy)*akmy(kmy+kmymax,kmyp)*u2(kmyp)
c
c-------->hole-hole block
c
          if(ind1(kmy).eq.lm1.and.ind1(kmyp).eq.lm2) 
     >    alms(lm1,lm2,1,1,2,2)=alms(lm1,lm2,1,1,2,2)
     >                   +u1(kmy)*akmy(kmy+kmymax,kmyp+kmymax)*u1(kmyp)
          if(ind1(kmy).eq.lm1.and.ind2(kmyp).eq.lm2) 
     >    alms(lm1,lm2,1,2,2,2)=alms(lm1,lm2,1,2,2,2)
     >                   +u1(kmy)*akmy(kmy+kmymax,kmyp+kmymax)*u2(kmyp)
          if(ind2(kmy).eq.lm1.and.ind1(kmyp).eq.lm2) 
     >    alms(lm1,lm2,2,1,2,2)=alms(lm1,lm2,2,1,2,2)
     >                   +u2(kmy)*akmy(kmy+kmymax,kmyp+kmymax)*u1(kmyp)
          if(ind2(kmy).eq.lm1.and.ind2(kmyp).eq.lm2) 
     >    alms(lm1,lm2,2,2,2,2)=alms(lm1,lm2,2,2,2,2)
     >                   +u2(kmy)*akmy(kmy+kmymax,kmyp+kmymax)*u2(kmyp)
        end do
        end do
      end do
      end do
c
      return
      end
      subroutine matkmy(akmy,alms,lmax)
c======================
c
c akmy: (kappa,my) representation
c alms: (l,m,s) representation
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      complex*16 alms(lmmaxp,lmmaxp,2,2),akmy(kmymaxp,kmymaxp)
      dimension u1(72),ind1(72)
      dimension u2(72),ind2(72)
      common/cgc/u1,u2,ind1,ind2
c
      lmmax=(lmax+1)*(lmax+1)
      kmymax=2*lmmax
c
      call czero(akmy,4*lmmaxp*lmmaxp)
      do k1=1,kmymax
        i1k1=ind1(k1)
        i2k1=ind2(k1)
      do k2=1,kmymax
        i1k2=ind1(k2)
        i2k2=ind2(k2)
          akmy(k1,k2)=u1(k1)*alms(i1k1,i1k2,1,1)*u1(k2)
     >               +u1(k1)*alms(i1k1,i2k2,1,2)*u2(k2)
     >               +u2(k1)*alms(i2k1,i1k2,2,1)*u1(k2)
     >               +u2(k1)*alms(i2k1,i2k2,2,2)*u2(k2)
      end do
      end do
c
      return
      end
c
c
