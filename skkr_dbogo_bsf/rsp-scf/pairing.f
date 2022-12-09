c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine genpairing(singrat,urat,drat,deltakmy)
c======================
c
c  i*sigma_y*conjg operator in (kappa,my) representation
c  combined with uu and dd triplet pairing
c
      implicit none
c
      include '../param.h'
c
      complex*16 xlms(lmmaxp,lmmaxp,2,2)
      complex*16 deltakmy(kmymaxp,kmymaxp)
c
      real*8 singrat,urat,drat
      integer lm,lmp,m,mp,ll,llp
      integer lml(lmmaxp),lmm(lmmaxp)
      common/mome/lml,lmm
c
      xlms(:,:,:,:)=(0.d0,0.d0)
      deltakmy=(0.d0,0.d0)
      do lm=1,lmmaxp
      do lmp=1,lmmaxp
        ll=lml(lm)
        llp=lml(lmp)
        m=lmm(lm)
        mp=lmm(lmp)
        if(ll.eq.(llp)) then
        if(m.eq.(-mp)) then
          xlms(lm,lmp,1,2)=singrat*(-1)**m
          xlms(lm,lmp,2,1)=-singrat*(-1)**m
          xlms(lm,lmp,1,1)=urat*(-1)**m
          xlms(lm,lmp,2,2)=drat*(-1)**m
        end if
        end if
      end do
      end do
c
      call matkmy(deltakmy,xlms,lmaxp)
c
      return
      end
c
c
c
      subroutine pairing(deltakmy)
c======================
c
c  i*sigma_y*conjg operator in (kappa,my) representation
c
      implicit none
c
      include '../param.h'
c
      complex*16 xlms(lmmaxp,lmmaxp,2,2)
      complex*16 deltakmy(kmymaxp,kmymaxp)
c
      integer lm,lmp,m,mp,ll,llp
      integer lml(lmmaxp),lmm(lmmaxp)
      common/mome/lml,lmm
c
      xlms(:,:,:,:)=(0.d0,0.d0)
      deltakmy=(0.d0,0.d0)
      do lm=1,lmmaxp
      do lmp=1,lmmaxp
        ll=lml(lm)
        llp=lml(lmp)
        m=lmm(lm)
        mp=lmm(lmp)
        if(ll.eq.(llp)) then
        if(m.eq.(-mp)) then
          xlms(lm,lmp,1,2)=0.5d0*(-1)**m
          xlms(lm,lmp,2,1)=-0.5d0*(-1)**m
        end if
        end if
      end do
      end do
c
      call matkmy(deltakmy,xlms,lmaxp)
c
      return
      end
c
c
c
      subroutine t0pairing(deltakmy)
c======================
c
c    triplet with unequal spin 
c
      implicit none
c
      include '../param.h'
c
      complex*16 xlms(lmmaxp,lmmaxp,2,2)
      complex*16 deltakmy(kmymaxp,kmymaxp)
c
      integer lm,lmp,m,mp,ll,llp
      integer lml(lmmaxp),lmm(lmmaxp)
      common/mome/lml,lmm
c
      xlms(:,:,:,:)=(0.d0,0.d0)
      deltakmy=(0.d0,0.d0)
      do lm=1,lmmaxp
      do lmp=1,lmmaxp
        ll=lml(lm)
        llp=lml(lmp)
        m=lmm(lm)
        mp=lmm(lmp)
        if(ll.eq.(llp)) then
        if(m.eq.(-mp)) then
           xlms(lm,lmp,1,2)=0.5d0*(-1)**mp
           xlms(lm,lmp,2,1)=0.5d0*(-1)**mp
        end if
        end if
      end do
      end do
c
      call matkmy(deltakmy,xlms,lmaxp)
c
      return
      end
c
c
      subroutine tuppairing(deltakmy)
c======================
c
c
      implicit none
c
      include '../param.h'
c
      complex*16 xlms(lmmaxp,lmmaxp,2,2)
      complex*16 deltakmy(kmymaxp,kmymaxp)
c
      integer lm,lmp,m,mp,ll,llp
      integer lml(lmmaxp),lmm(lmmaxp)
      common/mome/lml,lmm
c
      xlms(:,:,:,:)=(0.d0,0.d0)
      deltakmy=(0.d0,0.d0)
      do lm=1,lmmaxp
      do lmp=1,lmmaxp
        ll=lml(lm)
        llp=lml(lmp)
        m=lmm(lm)
        mp=lmm(lmp)
        if(ll.eq.(llp)) then
        if(m.eq.(-mp)) then
c          xlms(lm,lmp,1,1)=((-1)**mp+(-1)**m)/2 !G.Csire version
          xlms(lm,lmp,1,1)=0.5d0*(-1)**mp
          xlms(lm,lmp,2,2)=-0.5d0*(-1)**mp
        end if
        end if
      end do
      end do
c
      call matkmy(deltakmy,xlms,lmaxp)
c
      return
      end
c
c
      subroutine tdownpairing(deltakmy)
c======================
c
c
      implicit none
c
      include '../param.h'
c
      complex*16 xlms(lmmaxp,lmmaxp,2,2)
      complex*16 deltakmy(kmymaxp,kmymaxp)
c
      integer lm,lmp,m,mp,ll,llp
      integer lml(lmmaxp),lmm(lmmaxp)
      common/mome/lml,lmm
c
      xlms(:,:,:,:)=(0.d0,0.d0)
      deltakmy=(0.d0,0.d0)
      do lm=1,lmmaxp
      do lmp=1,lmmaxp
        ll=lml(lm)
        llp=lml(lmp)
        m=lmm(lm)
        mp=lmm(lmp)
        if(ll.eq.(llp)) then
        if(m.eq.(-mp)) then
c          xlms(lm,lmp,2,2)=-((-1)**mp+(-1)**m)/2 !G.Csire version
          xlms(lm,lmp,1,1)=-0.5d0*(-1)**mp*CMPLX(0.0d0,1.0d0)
          xlms(lm,lmp,2,2)=-0.5d0*(-1)**mp*CMPLX(0.0d0,1.0d0)
        end if
        end if
      end do
      end do
c
      call matkmy(deltakmy,xlms,lmaxp)
c
      return
      end
