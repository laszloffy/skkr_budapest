c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine genlatt2d
c
c=====================
c
c input:  a1,a2 in common/brav2d/ - basis translations in direct space
c         b1,b2 in common/brill2d/ - basis translations in rec. space
c work:   rho - non-primitive lattice  vector in 2D
c output: rvec,distr,maxr in common/dirvec2d/ - explanation in vector2d
c         indr,nshr,maxshr in common/vindr2d/   -
c               maxshr : number of shells around rho
c               nshr(i): number of vectors in shell i
c               indr(i): order index of real space vector: if i is the
c                        ordinal number then rvec(2,indr(i)) is the 
c                        corresponding vector  
c
c         the same variables with k instead of r are used in recip.
c         space
c ====================
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      dimension rho(2),ak(2)
      common/brav2d/a1(2),a2(2)
      common/brill2d/b1(2),b2(2)
c
      common/dirvec2d/rvec(2,mdir),distr,maxr
      common/vindr2d/indr(mdir),nshr(mr),maxshr
c
      common/recvec2d/veck(2,mrec),distk,maxk
      common/vindk2d/indk(mrec),nshk(mk),maxshk
c
      common/test/itest
c
c generate direct and reciprocal lattices
c
      call vector2d(a1,a2,b1,b2)
c
c sort real space vectors around origin 
c
      rho(1)=0.d0
      rho(2)=0.d0
      call vecsrt2d(maxr,mr,rvec,rho,distr,indr,nshr,maxshr)
c
c sort reciprocal space vectors around ak
c
      ak(1)=0.d0 
      ak(2)=0.d0 
      call vecsrt2d(maxk,mk,veck,ak,distk,indk,nshk,maxshk)
c
      if(itest.ge.1) then
      write(6,'(/2x,''routine GENLATT2D>''/)')
      write(6,'('' distr : '',f12.7)') distr
      write(6,'('' maxshr: '',i4)') maxshr
      write(6,'('' maxr: '',i7)') maxr
      write(6,'('' distk : '',f12.7)') distk
      write(6,'('' maxshk: '',i4)') maxshk
      write(6,'('' maxk: '',i7)') maxk
      end if
c
      return
      end
