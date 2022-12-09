c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine directg (park,sep,ecomp,lmax,gp,gm)
c     =================
c
c----- generates the bloch green functn for propagation between the
c      layers.  since z > 0, we dont need to use kambe split, but
c      can sum in reciprocal space
c ---------------------------------------------------------------------
c     implicit automatic (a-h,o-z)
      implicit real*8 (a-h,o-z)
c
      include  '../param.h'
      parameter (dthr=1.0d-12)
c
      complex*16 gp,gm,ylm,ylm1,ylm2,ylm3,ylm4,ecomp,il,kappa,
     > kappa2,f1,fp,fm,ctermp,ctermm,cth,sth,phi,gkz,gkz2,cylm
      complex*16 fac
c
      common/test/itest
      common/recvec2d/veck(2,mrec),distk,maxk
      common/vindk2d/indk(mrec),nshk(mk),maxshk

      common/crystl2d/vol,volbz
c  -------------------------------------------------------------
c  veck(2,mrec): lattice vectors of the 2d reciprocal lattice; 
c  indk(mrec):   containes the index of the vectors from sorting (not
c                used in this routine)
c  nshk(mk):     gives the number of rec. space vectors in shell mr
c  maxshk:       gives the number of shells in reciprocal space
c  vol,volbz:    volume of the elementary cell of the 2d direct and rec.
c                lattice
c  mrec,mk:      parameters set in 'param.h'
c  ---------------------------------------------------------------

      dimension sep(3),park(2)
      dimension  gp(lmmaxp,lmmaxp),gm(lmmaxp,lmmaxp)
      dimension ylm(lmmaxp),ylm1(lmmaxp),ylm2(lmmaxp),ylm3(lmmaxp),
     > ylm4(lmmaxp),cylm(-lsup:lsup,0:lsup)
c
c----- initialise
c
      area=vol
      akx=park(1)
      aky=park(2)
c
      pi=4.0d0*atan(1.0d0)
      lmx=lmax+1
      lmx2=lmx*lmx
      kappa2=ecomp
      kappa=sqrt(kappa2)
      if(dimag(kappa).lt.0.d0) kappa=-kappa
      f1=8.0d0*pi*pi/(area)
c
      do j=1,lmx2
        do i=1,lmx2
          gp(i,j)=0.0
          gm(i,j)=0.0
        enddo
      enddo
c
c----- loop over reciprocal lattice vectors
c
      if(itest.gt.5) write(6,*) 'DIRECTG> sep :',(sep(i),i=1,3)
c     write(6,*) 'DIRECTG> :'
c     write(6,'(3f10.5)') (sep(i),i=1,3)
c     write(6,'(2f10.5)') akx,aky
      num=0
      do i1=1,maxshk
        term=0.0d0
        do j1=1,nshk(i1)
          num=num+1
          ig=indk(num)
          if(itest.gt.5) write(6,*) 'reciprocal vec:',(veck(i,ig),i=1,2)
          gkx=akx+veck(1,ig)
          gky=aky+veck(2,ig)
          gk2=gkx*gkx+gky*gky
          gk=sqrt(gk2)
          gkz2=kappa2-gk2
          gkz=sqrt(gkz2)
          phi=1.0d0
          if (gk.gt.1.0e-06) phi=dcmplx(gkx,gky)/gk
          cth=gkz/kappa
          sth=gk/kappa
c
c----- generate the spherical harmonics 
c
          call sphrm (ylm,lmx2,cth,sth,phi,lmax)
c 
          lm=0
          il=(0.0d0,-1.0d0)
          sl=-1.0d0
          sm=-1.0d0
          do l=0,lmax
            sl=-sl
            il=il*(0.0d0,1.0d0)
            do m=-l,l
              sm=-sm
              lm=lm+1
              lmmm=lm-m-m
              ylm1(lm)=sl*il*ylm(lm)
              ylm2(lm)=il*sm*ylm(lmmm)
              ylm3(lm)=il*sm*ylm(lm)
              ylm4(lm)=sl*il*ylm(lmmm)
            enddo
          enddo
c
          sum=0.0d0
          fp=exp((0.0d0,1.0d0)*( gkx*sep(1)+gky*sep(2)+gkz*sep(3)))
          fm=exp((0.0d0,1.0d0)*(-gkx*sep(1)-gky*sep(2)+gkz*sep(3)))
          fp=fp*f1/gkz
          fm=fm*f1/gkz
          fac=(0.d0,-1.d0)
          do lmi=1,lmx2
            do lmj=1,lmx2
              ctermp=ylm1(lmi)*fp*ylm2(lmj)
              ctermm=ylm3(lmi)*fm*ylm4(lmj)
              gp(lmj,lmi)=gp(lmj,lmi)+ctermp*fac
              gm(lmj,lmi)=gm(lmj,lmi)+ctermm*fac
              term=term+abs(ctermp)+abs(ctermm)
              sum=sum+abs(gp(lmi,lmj))+abs(gm(lmi,lmj))
            enddo
          enddo
        enddo
c
c----- check convergence in g
c
c       if (sum.gt.1.0e-6) then
          if((itest.gt.5).and.(i1.gt.1)) 
     &     write(6,'(''directg-convergence'',2e14.5)') term,sum
          if (term/sum.lt.dthr.and.i1.gt.1) go to 90
c       endif
c
      enddo
      write (6,*) ' error: bloch g failed to converge'
      write (6,'(" term = ",2e14.5)') term,sum
      stop
c
   90 continue
c
      if (itest.gt.3) then
        write(6,*) ' recip g(+) g(-)' 
        do i=1,lmx2
          do j=1,lmx2
            if(abs(gp(i,j))+abs(gm(i,j)).gt.1.0e-5) then
              write (6,'(2i5,2e14.5,2x,2e14.5)') i,j,gp(i,j),gm(i,j)
            endif
          enddo
        enddo
      endif
c
      return
      end
