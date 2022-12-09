c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine genlatt3d
c
c=========================
c
c input:  a1l,a2l,a3l in common/brav3dl/ - basis translations in direct space
c                                          in the left semi-infinite region
c         a1r,a2r,a3r in common/brav3dr/ - basis translations in direct space
c                                          in the right semi-infinite region
c
c work:   rho,ak - non-primitive lattice vectors in 3D
c
c output: b1l,b2l,b3l in common/brill3dl/ - basis translations in rec. space
c                                          in the left semi-infinite region
c         b1r,b2r,b3r in common/brill3dr/ - basis translations in rec. space
c                                          in the right semi-infinite region
c         rvecl,volwsl,distrl,maxrl in common/dirvec3dl/ 
c         rvecr,volwsr,distrr,maxrr in common/dirvec3dr/ - 
c                                          explanation in vector3d
c         indrl,nshrl,maxshrl in common/vindr3dl/    
c         indrr,nshrr,maxshrr in common/vindr3dr/   -
c               maxshr : number of shells around rho
c               nshr(i): number of vectors in shell i
c               indr(i): order index of real space vector: if i is the
c                        ordinal number then rvec(2,indr(i)) is the 
c                        corresponding vector  
c
c         the same variables with k instead of r are used in recip. space
c ====================
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      dimension rho(3),ak(3)
c
      common/brav3dl/a1l(3),a2l(3),a3l(3)
      common/brill3dl/b1l(3),b2l(3),b3l(3)
      common/dirvec3dl/rvecl(3,mdir),volwsl,distrl,maxrl
      common/vindr3dl/indrl(mdir),nshrl(mr),maxshrl
      common/recvec3dl/veckl(3,mrec),volbzl,distkl,maxkl
      common/vindk3dl/indkl(mrec),nshkl(mk),maxshkl
c
      common/brav3dr/a1r(3),a2r(3),a3r(3)
      common/brill3dr/b1r(3),b2r(3),b3r(3)
      common/dirvec3dr/rvecr(3,mdir),volwsr,distrr,maxrr
      common/vindr3dr/indrr(mdir),nshrr(mr),maxshrr
      common/recvec3dr/veckr(3,mrec),volbzr,distkr,maxkr
      common/vindk3dr/indkr(mrec),nshkr(mk),maxshkr
c
      common/lay2d/cvec(mtotal,3),nextra,nbulkl,nbulkr,
     &             nprc,ninprc(0:mprc+1)
c
      common/test/itest
c
      if(itest.ge.0) then
      write(6,'(/2x,''routine GENLATT3D>''/)')
      write(6,'(2x,''Left semi-infinite bulk:'')')
      write(6,'(2x,''------------------------'')')
      end if
c
c generate direct and reciprocal lattices
c
      call vector3d(a1l,a2l,a3l,nbulkl,b1l,b2l,b3l,
     >              rvecl,volwsl,distrl,maxrl,
     >              veckl,volbzl,distkl,maxkl)
      call flush(6)
c
c sort real space vectors around origin 
c
      rho(1)=0.d0
      rho(2)=0.d0
      rho(3)=0.d0
      call vecsrt3d(maxrl,mr,rvecl,rho,distrl,indrl,nshrl,maxshrl)
c
c sort reciprocal space vectors around ak
c
      ak(1)=0.d0 
      ak(2)=0.d0 
      ak(3)=0.d0 
      call vecsrt3d(maxkl,mk,veckl,ak,distkl,indkl,nshkl,maxshkl)
c
      if(itest.ge.0) then
      write(6,'(/2x,''Right semi-infinite bulk:'')')
      write(6,'(2x,''-------------------------'')')
      end if
c
c generate direct and reciprocal lattices
c
      call vector3d(a1r,a2r,a3r,nbulkr,b1r,b2r,b3r,
     >              rvecr,volwsr,distrr,maxrr,
     >              veckr,volbzr,distkr,maxkr)
c
c sort real space vectors around origin 
c
      rho(1)=0.d0
      rho(2)=0.d0
      rho(3)=0.d0
      call vecsrt3d(maxrr,mr,rvecr,rho,distrr,indrr,nshrr,maxshrr)
      call flush(6)
c
c sort reciprocal space vectors around ak
c
      ak(1)=0.d0 
      ak(2)=0.d0 
      ak(3)=0.d0 
      call vecsrt3d(maxkr,mk,veckr,ak,distkr,indkr,nshkr,maxshkr)
c
      return
      end
