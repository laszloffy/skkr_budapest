c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine tau2d(irel,wgeff,bulkgeo,lmax,lkmax,lkmaxp,nintfc,
     >                 tminv,g2d,x,y,tau,mdimp)
c==================================================================
c
c  Calculates the principal-layer diagonal elements of tau(park;E)
c  using the procedure described by Godfrin
c  (J.Phys. Condens. Matt. 3 7843 (1991).
c  Only diagonal blocks are calculated!!!
c
c-------------------------------------------------
c     (nonrel): mdimp=mdimnr, lkmax=(lmax+1)**2
c     (rel):    mdimp=mdimr,  lkmax=kmymax=2*lkmax
c-------------------------------------------------    
c
      implicit real*8 (a-h,o-z)
      include '../param.h'
      parameter(mdim=mdimr)
c
      logical bulkgeo
      complex*16 tminv(lkmaxp,lkmaxp,mintfc)
      complex*16 g2d(lmmaxp,lmmaxp,-minprc:minprc,minprc,0:mprc1+1)
      complex*16 x(mdimp,mdimp),y(mdimp,mdimp)
      complex*16 xx(mdim,mdim,mprc),yy(mdim,mdim,mprc)
      complex*16 tau00(mdim,mdim,mprc)
      complex*16 tau(lkmaxp,lkmaxp,mintfc),detl
      complex*16 m00(mdim,mdim),m01(mdim,mdim),m10(mdim,mdim)
      complex*16 help(mdim,mdim)
c
      common/test/itest
      common/lay2d/cvec(mtotal,3),nextra,nbulkl,nbulkr,
     &       nprc,ninprc(0:mprc+1)
      data tol/1.0d-10/
c
c  Calculate first the XX and YY matrices required in the inversion.
c  Since X and Y are the SSPO-s, the procedure starts from there.
c
      ndim=lkmax*ninprc(1)  
      do 11 i=1,ndim
      do 11 j=1,ndim
        yy(i,j,1)=y(i,j)
 11   continue
c
      ndim=lkmax*ninprc(nprc) 
      do 22 i=1,ndim
      do 22 j=1,ndim
        xx(i,j,nprc)=x(i,j)
 22   continue
c
c
c  Y matrices
c
      ilay = 0
      do ipl=1,nprc
        iplg2d=ipl
        if(bulkgeo) iplg2d=1
        ndim  = lkmax*ninprc(ipl)
        call packer(ninprc(ipl),ninprc(ipl-1),ninprc(ipl+1),
     &              g2d(1,1,-minprc,1,iplg2d),m00,m10,m01,
     &              lmax,lkmax,mdim,irel)
c
c     Store m00 for further use
        do 20 li=1,ninprc(ipl)
          ilay = ilay+1
          do 20 k1=1,lkmax
            kk1 = (li-1)*lkmax + k1
            do 20 k2=1,lkmax
               kk2 = (li-1)*lkmax + k2
               m00(kk1,kk2) = tminv(k1,k2,ilay) + m00(kk1,kk2)
 20     continue
        call repl(tau00(1,1,ipl),m00,ndim,mdim)
c
        if(ipl.gt.1) then
          ndim0 = lkmax*ninprc(ipl-1)
          call doubmt2(m10,yy(1,1,ipl),help,ndim,ndim0,ndim,mdim)
          call repl(yy(1,1,ipl),help,ndim,mdim)
        endif
c
        if(ipl.lt.nprc) then
          ndim1 = lkmax*ninprc(ipl+1)
          call submat(m00,yy(1,1,ipl),ndim,mdim)
          call gjinv(m00,ndim,mdim,detl)
          call doubmt2(m00,m01,yy(1,1,ipl+1),ndim,ndim,ndim1,mdim)
        endif
      enddo
c
c  X matrices       (m00 is dummy now)
c
      do ipl=nprc,1,-1
        iplg2d=ipl
        if(bulkgeo) iplg2d=1
        ndim  = lkmax*ninprc(ipl)
        call packer(ninprc(ipl),ninprc(ipl-1),ninprc(ipl+1),
     &              g2d(1,1,-minprc,1,iplg2d),m00,m10,m01,
     &              lmax,lkmax,mdim,irel)
c
        if(ipl.lt.nprc) then
          ndim1 = lkmax*ninprc(ipl+1)
          call doubmt2(m01,xx(1,1,ipl),m00,ndim,ndim1,ndim,mdim)
          call repl(xx(1,1,ipl),m00,ndim,mdim)
        endif
c
        if(ipl.gt.1) then
          ndim0 = lkmax*ninprc(ipl-1)
          call submat1(tau00(1,1,ipl),xx(1,1,ipl),m00,ndim,mdim)  
          call gjinv(m00,ndim,mdim,detl)
          call doubmt2(m00,m10,xx(1,1,ipl-1),ndim,ndim,ndim0,mdim)
        endif
      enddo
c
c
c  Diagonal blocks of tau-matrix
c
      ilay=0
      do iiprc=1,nprc
        ndim=lkmax*ninprc(iiprc)
c
        do i=1,ndim
        do j=1,ndim
           tau00(i,j,iiprc)=tau00(i,j,iiprc)-xx(i,j,iiprc)-yy(i,j,iiprc)
        end do
        end do
c
c       -------------------------------------------
        call gjinv(tau00(1,1,iiprc),ndim,mdim,detl)
c       -------------------------------------------
c       ----------------------------------------
c       call hmatinv(mdim,ndim,tau00(1,1,iiprc))
c       ----------------------------------------
c
        il0=0
        do 1 inprc=1,ninprc(iiprc)
          ilay=ilay+1
c
          do l=1,lkmax
          do lp=1,lkmax
            tau(l,lp,ilay)=wgeff*tau00(il0+l,il0+lp,iiprc)
          end do
          end do
c
          il0=il0+lkmax
c
          if(itest.le.2) goto 1
          write(6,'('' tau-matrix for layer: '',i5)') ilay
          call outmat1(tau(1,1,ilay),lkmax,lkmax,lkmaxp,tol,6)
c
   1    continue
c
      end do  
c
      return
      end
