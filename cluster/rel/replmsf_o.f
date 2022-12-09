c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine replmsf(alms,akmy,lmax)
c=======================
c
c find matrix in (l,m,s) representation
c to be given in (kappa,my) representation
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      complex*16 alms(kmymaxp,kmymaxp),akmy(kmymaxp,kmymaxp)
      dimension u1(50),ind1(50)
      dimension u2(50),ind2(50)
      common/cgc/u1,u2,ind1,ind2
c
      lmmax=(lmax+1)*(lmax+1)
      kmymax=2*lmmax
c
      lm1=0
      do l1=0,lmax
      do m1=-l1,l1
        lm1=lm1+1
        lmm1=lm1
        lmp1=lm1+lmmax
      lm2=0
      do l2=0,lmax
      do m2=-l2,l2
        lm2=lm2+1
        lmm2=lm2
        lmp2=lm2+lmmax
c
        alms(lmm1,lmm2)=(0.d0,0.d0)
        alms(lmm1,lmp2)=(0.d0,0.d0)
        alms(lmp1,lmm2)=(0.d0,0.d0)
        alms(lmp1,lmp2)=(0.d0,0.d0)
c
        do kmy=1,kmymax
        do kmyp=1,kmymax
c
          if(ind1(kmy).eq.lm1.and.ind1(kmyp).eq.lm2) 
     >    alms(lmm1,lmm2)=alms(lmm1,lmm2)+
     >                    u1(kmy)*akmy(kmy,kmyp)*u1(kmyp)
c
          if(ind1(kmy).eq.lm1.and.ind2(kmyp).eq.lm2) 
     >    alms(lmm1,lmp2)=alms(lmm1,lmp2)+
     >                    u1(kmy)*akmy(kmy,kmyp)*u2(kmyp)
c
          if(ind2(kmy).eq.lm1.and.ind1(kmyp).eq.lm2) 
     >    alms(lmp1,lmm2)=alms(lmp1,lmm2)+
     >                    u2(kmy)*akmy(kmy,kmyp)*u1(kmyp)
c
          if(ind2(kmy).eq.lm1.and.ind2(kmyp).eq.lm2) 
     >    alms(lmp1,lmp2)=alms(lmp1,lmp2)+
     >                    u2(kmy)*akmy(kmy,kmyp)*u2(kmyp)
c
        end do
        end do
c
      end do
      end do
      end do
      end do
c
      return
      end
