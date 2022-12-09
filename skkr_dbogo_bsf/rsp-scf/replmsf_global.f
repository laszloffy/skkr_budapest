c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine replmsf_global(alms,akmy,lmax,rb)
c======================
c
c find diagonal elements of a matrix in (l,m,s) representation
c to be given in (kappa,my) representation
c involves a spinor basis change according to the eigenstates of rb*S
c
      implicit none
c
      include '../param.h'
c
      complex*16,  intent(in) :: akmy(kmymaxp,kmymaxp)
      integer,     intent(in) :: lmax
      real*8,      intent(in) :: rb(3)
      complex*16, intent(out) :: alms(kmymaxp,kmymaxp)
c
      complex*16 :: Smat(2,2),Smatinv(2,2)
      integer    :: lmmax,kmymax,kmy,kmyp,lm1,lm2,lmm,lmp,s1,s2
      integer    :: l1,m1,lmm1,lmp1,l2,m2,lmm2,lmp2

      real*8     :: utmp(2,72),indtmp(2,72)
c
      real*8     :: u1(72),u2(72)
      integer    :: ind1(72),ind2(72)
      common/cgc/u1,u2,ind1,ind2
c
      lmmax = (lmax+1)*(lmax+1)
      kmymax = 2*lmmax
c
      ! construct spin-dependent clebsch and lm(kmy,s) arrays for looping
      utmp(1,:) = u1
      utmp(2,:) = u2
      indtmp(1,:) = ind1
      indtmp(2,:) = ind2
c
      ! get SU2 transformation matrix
      call SU2mat(rb,Smat,Smatinv)
c
      alms = (0.d0,0.d0)
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
          do s1=1,2
          do s2=1,2
            do kmy=1,kmymax
            do kmyp=1,kmymax
c
          if(indtmp(s1,kmy).eq.lm1.and.indtmp(s2,kmyp).eq.lm2) then
          alms(lmm1,lmm2)=alms(lmm1,lmm2)+
     >                    Smat(1,s1)*
     >                    utmp(s1,kmy)*akmy(kmy,kmyp)*utmp(s2,kmyp)*
     >                    Smatinv(s2,1)
c
c         if(ind1(kmy).eq.lm1.and.ind2(kmyp).eq.lm2) 
          alms(lmm1,lmp2)=alms(lmm1,lmp2)+
     >                    Smat(1,s1)*
     >                    utmp(s1,kmy)*akmy(kmy,kmyp)*utmp(s2,kmyp)*
     >                    Smatinv(s2,2)
c
c         if(ind2(kmy).eq.lm1.and.ind1(kmyp).eq.lm2) 
          alms(lmp1,lmm2)=alms(lmp1,lmm2)+
     >                    Smat(2,s1)*
     >                    utmp(s1,kmy)*akmy(kmy,kmyp)*utmp(s2,kmyp)*
     >                    Smatinv(s2,1)
c
c         if(ind2(kmy).eq.lm1.and.ind2(kmyp).eq.lm2) 
          alms(lmp1,lmp2)=alms(lmp1,lmp2)+
     >                    Smat(2,s1)*
     >                    utmp(s1,kmy)*akmy(kmy,kmyp)*utmp(s2,kmyp)*
     >                    Smatinv(s2,2)
          end if
c
            end do
            end do
          end do
          end do
        end do
        end do
      end do
      end do
c
c     do s1=1,2
c     do s2=1,2
c       do kmy=1,kmymax
c         lm1 = indtmp(s1,kmy)
c         if (lm1 > lmmax) cycle
c       do kmyp=1,kmymax
c         lm2 = indtmp(s2,kmyp)
c         if (lm2 > lmmax) cycle
c
c         if (lm1 == lm2) then
c           lmm = lm1
c           lmp = lm1 + lmmax
c           alms(lmm) = alms(lmm) +
c    +                     Smat(1,s1)*utmp(s1,kmy)*
c    +                     akmy(kmy,kmyp)*
c    +                     utmp(s2,kmyp)*Smatinv(s2,1)
c           alms(lmp) = alms(lmp) +
c    +                     Smat(2,s1)*utmp(s1,kmy)*
c    +                     akmy(kmy,kmyp)*
c    +                     utmp(s2,kmyp)*Smatinv(s2,2)
c         end if
c       end do
c       end do
c     end do
c     end do
c
      return
      end subroutine replmsf_global
