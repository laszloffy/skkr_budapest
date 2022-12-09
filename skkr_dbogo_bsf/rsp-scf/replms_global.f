c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine replms_global(alms,akmy,lmax,rb)
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
      complex*16, intent(out) :: alms(kmymaxp)
c
      complex*16 :: Smat(2,2),Smatinv(2,2)
      integer    :: lmmax,kmymax,kmy,kmyp,lm1,lm2,lmm,lmp,s1,s2
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
      do s1=1,2
      do s2=1,2
        do kmy=1,kmymax
          lm1 = indtmp(s1,kmy)
          if (lm1 > lmmax) cycle
        do kmyp=1,kmymax
          lm2 = indtmp(s2,kmyp)
          if (lm2 > lmmax) cycle

          if (lm1 == lm2) then
            lmm = lm1
            lmp = lm1 + lmmax
            alms(lmm) = alms(lmm) +
     +                     Smat(1,s1)*utmp(s1,kmy)*
     +                     akmy(kmy,kmyp)*
     +                     utmp(s2,kmyp)*Smatinv(s2,1)
            alms(lmp) = alms(lmp) +
     +                     Smat(2,s1)*utmp(s1,kmy)*
     +                     akmy(kmy,kmyp)*
     +                     utmp(s2,kmyp)*Smatinv(s2,2)
          end if
        end do
        end do
      end do
      end do
c
      return
      end subroutine replms_global
