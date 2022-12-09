c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine tau2dgod(irel,wgeff,bulkgeo,lmax,lkmax,lkmaxp,nintfc,
     >                 tminv,g2d,x,y,tau,mdimp)
c==================================================================
c
c  Calculates the full matrix of tau(park;E)
c  using the procedure described by Godfrin
c  (J.Phys. Condens. Matt. 3 7843 (1991).
c  Only diagonal blocks are calculated!!!
c
c-------------------------------------------------
c     (nonrel): mdimp=mdimnr, lkmax=(lmax+1)**2
c     (rel):    mdimp=mdimr,  lkmax=kmymax=2*lkmax
c-------------------------------------------------    
c
      implicit none
c     implicit real*8 (a-h,o-z)
      include '../param.h'
      integer mdim
      parameter(mdim=mdimr)
c
      INTEGER*4 NEXTRA ,NBULKL ,NBULKR ,NPRC ,NINPRC ,ITEST
      INTEGER*4 IIPRC ,IPL ,NDIM ,I ,J ,ILAY
      INTEGER*4 IPLG2D ,LI ,K1 ,KK1 ,K2 ,KK2
      INTEGER*4 NDIM0 ,NDIM1 ,II0 ,IL0 ,JLAY ,JL0
      INTEGER*4 IL ,JL ,JJPRC ,JDIM ,JDIM0 ,NN
      INTEGER*4 JJ0 ,LLAY ,LL0 ,L ,LP ,JDIM1
      INTEGER*4 IREL ,LMAX ,LKMAX ,LKMAXP ,NINTFC ,MDIMP
      REAL*8 WGEFF
      REAL*8 CVEC ,TOL
c      
      logical bulkgeo
      complex*16 tminv(lkmaxp,lkmaxp,mintfc)
      complex*16 g2d(lmmaxp,lmmaxp,-minprc:minprc,minprc,0:mprc1+1)
      complex*16 x(mdimp,mdimp),y(mdimp,mdimp)
      complex*16 xx(mdim,mdim,mprc),yy(mdim,mdim,mprc)
      complex*16 tau00(mdim,mdim,mprc)
      complex*16 tau(lkmaxp,lkmaxp,mintfc,mintfc),detl
      complex*16 m00(mdim,mdim),m01(mdim,mdim),m10(mdim,mdim)
      complex*16 help(mdim,mdim),cc(mdim,mdim,mprc),dd(mdim,mdim,mprc)
      complex*16 help1(mdim,mdim)
c
      integer    nnprc(mprc),nlay
      common/test/itest
      common/lay2d/cvec(mtotal,3),nextra,nbulkl,nbulkr,
     &       nprc,ninprc(0:mprc+1)
      data tol/1.0d-10/
c
c  fill nnprc array
      iiprc = 0
      do ipl = 1,nprc
       nnprc(ipl) = iiprc
       iiprc = iiprc + ninprc(ipl)
      end do 
c
c  Calculate first the XX and YY matrices required in the inversion.
c  Since X and Y are the SSPO-s, the procedure starts from there.
c
      ndim=lkmax*ninprc(1)  
      do  i=1,ndim
      do  j=1,ndim
        yy(i,j,1)=y(i,j)
      end do
      end do
c
      ndim=lkmax*ninprc(nprc) 
      do  i=1,ndim
      do  j=1,ndim
        xx(i,j,nprc)=x(i,j)
      end do
      end do
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
        do  li=1,ninprc(ipl)
          ilay = ilay+1
          do  k1=1,lkmax
            kk1 = (li-1)*lkmax + k1
            do  k2=1,lkmax
               kk2 = (li-1)*lkmax + k2
               m00(kk1,kk2) = tminv(k1,k2,ilay) + m00(kk1,kk2)
            end do
          end do
        end do
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
c store dd (only for the full matrix inversion)          
          call repl(dd(1,1,ipl),yy(1,1,ipl+1),ndim1,mdim)
          do i = 1,ndim1
           do j = 1,ndim
            dd(j,i,ipl) = -dd(j,i,ipl)
           end do 
          end do 
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
c store cc (only for the full matrix inversion)          
          call repl(cc(1,1,ipl),xx(1,1,ipl-1),ndim0,mdim)
          do i = 1,ndim0
           do j = 1,ndim
            cc(j,i,ipl) = -cc(j,i,ipl)
           end do 
          end do 
        endif
      enddo
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
c --------------- original version ----------------------
c       il0=0
c       do  inprc=1,ninprc(iiprc)
c         ilay=ilay+1
c
c         do l=1,lkmax
c         do lp=1,lkmax
c tau matrix is multiplied by weight factor          
c           tau(l,lp,ilay,ilay)=wgeff*tau00(il0+l,il0+lp,iiprc)
c           tau(l,lp,ilay,ilay)=tau00(il0+l,il0+lp,iiprc)
c         end do
c         end do
c
c         il0=il0+lkmax
c
c         if(itest.gt.2) then
c          write(6,'('' tau-matrix for layer: '',i5)') ilay
c          call outmat1(tau(1,1,ilay,ilay),lkmax,lkmax,lkmaxp,tol,6)
c         endif 
c
c       end do
c --------------- end of original version ----------------------
        ii0 = nnprc(iiprc)
        nlay = ninprc(iiprc)
        do ilay = 1,nlay
         il0 = (ilay - 1)*lkmax
         do jlay = 1,nlay
          jl0 = (jlay - 1)*lkmax
          do il = 1,lkmax
           do jl = 1,lkmax
            tau(il,jl,ilay+ii0,jlay+ii0) = tau00(il+il0,jl+jl0,iiprc)
           end do 
          end do 
         end do 
c
          if(itest.gt.2) then
           write(6,'('' tau-matrix for layer: '',i5)') ilay
           call outmat1(tau(1,1,ilay,ilay),lkmax,lkmax,lkmaxp,tol,6)
          endif 
c
         end do 
c
c calculates off-diagonal elements of the tau matrix
c
        call repl(help,tau00(1,1,iiprc),ndim,mdim)
        do jjprc = iiprc+1,nprc
         jdim = lkmax*ninprc(jjprc)
         jdim0 = lkmax*ninprc(jjprc-1)
         nn = max(jdim,jdim0)
         call doubmt2(cc(1,1,jjprc),help,help1,jdim,jdim0,jdim0,mdim)
         call repl(help,help1,nn,mdim)
c
         ii0 = nnprc(iiprc)
         jj0 = nnprc(jjprc)
         do llay = 1,ninprc(iiprc)
          ll0 = (llay-1)*lkmax
          do jlay = 1,ninprc(jjprc)
           jl0 = (jlay-1)*lkmax
           do l = 1,lkmax
            il = ll0 + l
            do lp = 1,lkmax
             jl = jl0 + lp
             tau(lp,l,jlay+jj0,llay+ii0) = help(jl,il)
            end do 
           end do 
          end do 
         end do 
        end do 
c          
        call repl(help,tau00(1,1,iiprc),ndim,mdim)
        do jjprc = iiprc-1,1,-1
         jdim = lkmax*ninprc(jjprc)
         jdim1 = lkmax*ninprc(jjprc+1)
         nn = max(jdim,jdim1)
         call doubmt2(dd(1,1,jjprc),help,help1,jdim,jdim1,jdim1,mdim)
         call repl(help,help1,nn,mdim)
c
         ii0 = nnprc(iiprc)
         jj0 = nnprc(jjprc)
         do llay = 1,ninprc(iiprc)
          ll0 = (llay-1)*lkmax
          do jlay = 1,ninprc(jjprc)
           jl0 = (jlay-1)*lkmax
           do l = 1,lkmax
            il = ll0 + l
            do lp = 1,lkmax
             jl = jl0 + lp
             tau(lp,l,jlay+jj0,llay+ii0) = help(jl,il)
            end do 
           end do 
          end do 
         end do 
        end do 
c          
      end do  
c
      return
      end
