c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine relmtrx(a,b,lmax)
c=======================
c
c ***************************************************************
c * transform structure constant from a non-relativistic matrix *
c *                to a relativistic matrix                     *
c ***************************************************************
c
      implicit real*8(a-h,o-z)
      include '../param.h'
c
      integer lmax,lmmax,kmymax,i,j,i1,i2,j1,j2
      complex*16 a(bogomaxp,bogomaxp),b(dbogomaxp,dbogomaxp)
      complex*16 a1(36,36)
c
      integer ind1(72)
      integer ind2(72)
      real*8 u1(72)
      real*8 u2(72)
      common/cgc/u1,u2,ind1,ind2
c
      complex*16 fac,fach
      common/relfac/fac,fach
c
      lmmax=(lmax+1)*(lmax+1)
      kmymax=2*lmmax
      if(lmmax.gt.36) stop 'Increase dimension of a1 in relmtrx!'
c
c     write(6,'('' relfac'',2d15.6)') fac
      a1=(0.d0,0.d0)
      b=(0.d0,0.d0)
      do i=1,(2*lmmax)
      do j=1,(2*lmmax)
        a1(i,j)=fac*a(i,j)
        if(j.gt.lmmax) then
           a1(i,j)=fach*a(i,j)
        end if
      end do
      end do
c
      do i=1,kmymax
        i1=ind1(i)
        i2=ind2(i)
        do j=1,kmymax
          j1=ind1(j)
          j2=ind2(j)
          b(i,j)=u1(i)*a1(i1,j1)*u1(j)+
     >           u2(i)*a1(i2,j2)*u2(j)
          b(i+kmymax,j+kmymax)=u1(i)*a1(i1+lmmax,j1+lmmax)*u1(j)+
     >           u2(i)*a1(i2+lmmax,j2+lmmax)*u2(j)
        end do
      end do
c
      return
      end
c
c======================================================================
c
      subroutine relmtrxg(a,b,lmax)
c=======================
c
c *****************************************************
c * transformation from a non-relativistic matrix 'a' *
c *                to a relativistic matrix 'b'       *
c *****************************************************
c   (applied for gacoeff)
c
c
      implicit real*8(a-h,o-z)
      include '../param.h'
c
      integer lmax,lmmax,kmymax,i,j,i1,i2,j1,j2
c
      complex*16 a(lmmaxp,lmmaxp),b(kmymaxp,kmymaxp)
c
      integer ind1(72)
      integer ind2(72)
      real*8 u1(72)
      real*8 u2(72)
      common/cgc/u1,u2,ind1,ind2
c
      lmmax=(lmax+1)*(lmax+1)
      kmymax=2*lmmax
      if(lmmax.gt.36) stop 'Increase dimension of a1 in relmtrx!'
c
      do i=1,kmymax
        i1=ind1(i)
        i2=ind2(i)
        do j=1,kmymax
          j1=ind1(j)
          j2=ind2(j)
          b(i,j)=u1(i)*a(i1,j1)*u1(j)+
     >           u2(i)*a(i2,j2)*u2(j)
        end do
      end do
c
      return
      end
c
c======================================================================
c
      subroutine lmtolmconj(a)
c=======================
c
c *****************************************************
c * transform a matrix from |lm> basis                *
c *                to |lm>* basis                     *
c *****************************************************
c
c
      implicit none
      include '../param.h'
c
      complex*16 a(lmmaxp,lmmaxp)
      complex*16 trafo(lmmaxp,lmmaxp)
      integer lm,lmp,m,mp,ll,llp
      complex*16 detl
      integer lml(lmmaxp),lmm(lmmaxp)
      common/mome/lml,lmm
c
c
      trafo(:,:)=0.0d0
      do lm=1,lmmaxp
      do lmp=1,lmmaxp
        ll=lml(lm)
        llp=lml(lmp)
        m=lmm(lm)
        mp=lmm(lmp)
        if(ll.eq.(llp)) then
        if(m.eq.(-mp)) then
          trafo(lm,lmp)=(-1)**m
        end if
        end if
      end do
      end do 
c
      a=matmul(trafo,a)
      call gjinv(trafo,lmmaxp,lmmaxp,detl)
      a=matmul(a,trafo)
c
      return
      end
c
c
      subroutine conjinkmy(ctrafo)
c======================
c
c     conjg operator in (kappa,my) representation
c
      implicit none
c
      include '../param.h'
c
      complex*16 xlms(lmmaxp,lmmaxp,2,2)
      complex*16 conjkmy(kmymaxp,kmymaxp)
      complex*16 ctrafo(dbogomaxp,dbogomaxp)
c
      integer lm,lmp,m,mp,ll,llp,kmy
      integer lml(lmmaxp),lmm(lmmaxp)
      common/mome/lml,lmm
c
      xlms(:,:,:,:)=(0.d0,0.d0)
      conjkmy=(0.d0,0.d0)
      do lm=1,lmmaxp
      do lmp=1,lmmaxp
        ll=lml(lm)
        llp=lml(lmp)
        m=lmm(lm)
        mp=lmm(lmp)
        if(ll.eq.(llp)) then
        if(m.eq.(-mp)) then
          xlms(lm,lmp,1,1)=(-1)**m
          xlms(lm,lmp,2,2)=(-1)**m
        end if
        end if
      end do
      end do
c
      call matkmy(conjkmy,xlms,lmaxp)
c
      ctrafo=(0.d0,0.d0)
      do kmy=1,kmymaxp
        ctrafo(kmy,kmy)=(1.d0,0.d0)
      end do
      ctrafo((kmymaxp+1):2*kmymaxp,(kmymaxp+1):2*kmymaxp)=conjkmy  
c
c
      return
      end

