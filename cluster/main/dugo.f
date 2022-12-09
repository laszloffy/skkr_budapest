c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
c**********************************************************************
c INTRODUCE A MAX.NUMBER OF ITERATIONS IN surfgf BASED ON A
c LIMITING NUMBER OF INTERACTING LAYERS (decay of structure constants).
c
c (bulk) case: use the 3D repetition unit (see new madelung) to stop 
c number of iterations
c**********************************************************************
      subroutine dugo(lmax,lkmax,tminvb,lkmaxp,g2d,g2di,xy,mdimp,
     &                side,bulk,irel,sgfconv)
c----------------------------------------------------------------
c  _cpl  = current P.L (0 for LEFT, nprc+1 for RIGHT)
c  _cpli = adjacent P.L of Intermediate-Region (either 1, nprc)
c
c----------------------------------------------------------------
c     (nonrel):  lkmax=(lmax+1)**2
c     (rel):     lkmax=kmymax=2*lkmax
c----------------------------------------------------------------
c
c
      implicit real*8 (a-h,o-z)
      include '../param.h'
      parameter(mdim=mdimr)
c
      character*1 side
      logical bulk,sgfconv
      complex*16 g2d(lmmaxp,lmmaxp,-minprc:minprc,minprc)
      complex*16 g2di(lmmaxp,lmmaxp,-minprc:minprc,minprc)
      complex*16 tminvb(lkmaxp,lkmaxp,minprc)
      complex*16 m00(mdim,mdim),m01(mdim,mdim),m10(mdim,mdim)
      complex*16 mt(mdim,mdim),detl
      complex*16 xy(mdimp,mdimp),xystart(mdim,mdim)
c
      common/test/itest
      common/iterpar/errmax,itermaxl,itermaxr,ichk
      common/lay2d/cvec(mtotal,3),nextra,nbulkl,nbulkr,
     &              nprc,ninprc(0:mprc+1)
      data tol/1.0d-08/
c
c
c  First, obtain mt-matrix for tau2d (dummy 'm00,m01/10').
c  Then, the largest m-matrices to be used in surfgfd and tau2d.
c
      if(.not.bulk) then
        if(side.eq.'L') then
          ncpl  = ninprc(0)
          ncpli = ninprc(1)
          nip = max(ncpl,ncpli)
          nim = ncpl
          call packer(ncpli,ncpl,ninprc(2),g2di,m00,mt,m01,
     &                lmax,lkmax,mdim,irel)
        endif
c
        if(side.eq.'R') then
          ncpl  = ninprc(nprc+1)
          ncpli = ninprc(nprc)
          nim = max(ncpl,ncpli)
          nip = ncpl
          call packer(ncpli,ninprc(nprc-1),ncpl,g2di,m00,m10,mt,
     &                lmax,lkmax,mdim,irel)
        endif
c
        call packer(ncpl,nim,nip,g2d,m00,m10,m01,
     &              lmax,lkmax,mdim,irel)
c
c-- bulk case
c
      else
        ncpl = ninprc(1)
        ncpli = ncpl
        if(side.eq.'L') call packer(ncpl,ncpl,ncpl,g2d,m00,m10,mt,
     &                              lmax,lkmax,mdim,irel)
        if(side.eq.'R') call packer(ncpl,ncpl,ncpl,g2d,m00,mt,m01,
     &                              lmax,lkmax,mdim,irel)
        call packer(ncpl,ncpl,ncpl,g2di,m00,m10,m01,
     &              lmax,lkmax,mdim,irel)
      endif
c
c
      nn = ncpl * lkmax
      ni = ncpli * lkmax
      do inprc=0,ncpl-1
         do kmy=1,lkmax
         do kmyp=1,lkmax
            i=inprc*lkmax+kmy
            j=inprc*lkmax+kmyp
            m00(i,j)=tminvb(kmy,kmyp,inprc+1)+m00(i,j)
         end do
         end do
      end do
c
      if(ichk.eq.1) write(6,'('' Side: '',a1)') side
c     ----------------------------------------------
      if(side.eq.'L') then
        call surfgfd(m10,m00,m01,xystart,nn,itermaxl,errmax,ichk,
     >               sgfconv)
        if(bulk) then
          call tripmt2(m10,xystart,mt,xy,ni,nn,nn,ni,mdim,mdimp)
        else
          call tripmt2(mt,xystart,m01,xy,ni,nn,nn,ni,mdim,mdimp)
        endif
      else 
        call surfgfd(m01,m00,m10,xystart,nn,itermaxr,errmax,ichk,
     >               sgfconv)
        if(bulk) then
          call tripmt2(m01,xystart,mt,xy,ni,nn,nn,ni,mdim,mdimp)
        else
          call tripmt2(mt,xystart,m10,xy,ni,nn,nn,ni,mdim,mdimp)
        endif
      end if
c     ----------------------------------------------
ccccccccccc
      if(itest.ge.3) then
      WRITE(64,*) 'DUGO: XSAVE  ',ni,nn
      CALL OUTMAT1(xy,ni,ni,mdimp,tol,64)
      endif
ccccccccccc
c
      return
      end
