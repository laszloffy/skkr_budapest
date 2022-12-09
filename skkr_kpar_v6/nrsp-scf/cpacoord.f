c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine cpacoord(
     > nsp,itscf,ie,ce,we,lmax,nintfc,eta,rightm,
     > iesublatt,bulk,slab,bulkgeo,kset,xk,wk,nk,intbz,iek,
     > conc,itcpam,cpatol,cpatest,
     > tminvl,tminvr,tminv,tminva,tminvb,taua,taub,taukdiag)
c
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'

      include 'mpif.h'
      parameter(mdim=mdimnr)
      parameter(mek=meknr)
c
      logical bulk,slab,bulkgeo,sgf,sgfx,sgfy,sgfconv
      logical cpacon,concpa,cpatest,cpalay
c
c MPI declaration      
      integer root, myrank, nprocs, ierr
      common/mpi/root,myrank,nprocs
      complex*16 taujoin(lmmaxp,lmmaxp,mintfc)

      character*1 rightm
c
      dimension park(2),xk(mkpar,2),wk(mkpar)
      dimension conc(mintfc)
      dimension invg(melem),iesublatt(mintfc)
c
      complex*16 tminv(lmmaxp,lmmaxp,mintfc)
      complex*16 tminva(lmmaxp,lmmaxp,mintfc)
      complex*16 tminvb(lmmaxp,lmmaxp,mintfc)
      complex*16 tminvl(lmmaxp,lmmaxp,minprc)
      complex*16 tminvr(lmmaxp,lmmaxp,minprc)
      complex*16 tau(lmmaxp,lmmaxp,mintfc)
      complex*16 tau00(lmmaxp,lmmaxp,mintfc)
      complex*16 taua(lmmaxp,lmmaxp,mintfc)
      complex*16 taub(lmmaxp,lmmaxp,mintfc)
      complex*16 dg(lmmaxp,lmmaxp,melem)
      complex*16 dg1(lmmaxp,lmmaxp,melem)
      complex*16 tau1(lmmaxp,lmmaxp)
      complex*16 taukdiag(mkpar,0:lmaxp,mintfc)
c
      complex*16 ce,we,cwgeff
c
      complex*16 xx(mdim,mdim),yy(mdim,mdim),
     > gg2d(lmmaxp,lmmaxp,-minprc:minprc,minprc,0:mprc1+1),
     > gg2dbulk(lmmaxp,lmmaxp,-minprc:minprc,minprc,0:mprc1+1) 
c    > gg2dbulk(lmmaxp,lmmaxp,-minprc:minprc,minprc,0:2) 
c
      complex*16 xsave(mdim,mdim,mek+1,2),
     >           ysave(mdim,mdim,mek+1,2),
     >g2d(lmmaxp,lmmaxp,-minprc:minprc,minprc,0:mprc1+1,mek+1)
c
      complex*16 alphalkkr(0:lmaxp,minprc)
      complex*16 alpharkkr(0:lmaxp,minprc)
      complex*16 alphaintkkr(0:lmaxp,mintfc)
      common/scrpar/alphalkkr,alpharkkr,alphaintkkr  
c
      common/test/itest
      common/lay2d/cvec(mtotal,3),nextra,nbulkl,nbulkr,
     &             nprc,ninprc(0:mprc+1)
      common/iterpar/errmax,itermaxl,itermaxr,ichk
c
      data outtol/1.0d-12/
      data tol/1.0d-8/
      data tiny/1.0d-6/
c
      save xsave,ysave,g2d 
c
      call flush(6)
c
c ********************
c initialize constants
c ********************
c 
      irel=0
c
      nl=lmax+1
      nl2=nl*nl
      l2=2*lmax
c
      ninprcl=ninprc(0)
      ninprcr=ninprc(nprc+1)
      ndiml=ninprcl*nl2
      ndimr=ninprcr*nl2
c
      nprc1=nprc
      if(bulkgeo) nprc1=1
      if(nprc1.gt.mprc1) stop ' increase mprc1 !!!'
c
c ****************************************************************
c initialize irreducible representations of point group operations
c ****************************************************************
c
c     ---------------------------
      call sphbas(dg,dg1,ng,invg)
c     ---------------------------
      if(intbz.eq.0.or.intbz.eq.2) ngeff=1
      if(intbz.eq.1) ngeff=ng
c
c
c                        ***************************
c                        * starting CPA iterations *
c                        ***************************
c
      itcpa=1
      concpa=.true.        ! controls CPA tolerance
      cpacon=.true.        ! controls no. of CPA iterations 
c
      iekstart=iek
  100 continue
      iek=iekstart
c
c *******************************
c adjust t-matrices for bulk case
c *******************************
c
      if(bulk) then
c The I region is supposed to comprise one PL only !
        if(kset.ge.1) then
c Force L and R t-matrices to be identical with the
c corresponding t-matrices in I
          do li=1,ninprcl
            call repl(tminvl(1,1,li),tminv(1,1,li),nl2,lmmaxp)
            call repl(tminvr(1,1,li),tminv(1,1,li),nl2,lmmaxp)
          enddo
        else
c This is for only spectral-DOS calculation: force I t-matrices
c to be identical with the corresponding self-consistent L t-matrices
          do li=1,nintfc
            call repl(tminv(1,1,li),tminvl(1,1,li),nl2,lmmaxp)
          end do
        end if
      end if
c
c
c ****************************
c loop over k points in 2D IBZ
c ****************************
c 
c     -------------------------------------
      call czero(tau,lmmaxp*lmmaxp*mintfc)
c     -------------------------------------
      nk2=0
      do ik=1,nk
c
      if(modulo(ik-1,nprocs) == myrank) then 
         park(1)=xk(ik,1)
         park(2)=xk(ik,2)
c        write(6,'('' k='',2f13.5)') park
c
c *************************************************
c initialize saving of structure constants and SSPO
c *************************************************
c
         iek=iek+1
         if(iek.le.mek) then
           iekuse=iek
         else
           iekuse=mek+1
         end if
         if(iek.eq.mek+1) 
     >   write(6,'(/''WARNING: iek='',i3,'' greater than mek='',i3/
     >              '' causes dramatic slow-down !'')') iek,mek
c
         if(bulk.or.slab) then
c    calculate g2d only for the first CPA iteration and spin channel
           nocalcstr=itcpa+nsp-2
         else
c    calculate g2d only for the 1st. scf, CPA iteration and spin channel
           nocalcstr=itscf+itcpa+nsp-3
         end if
c    do not store g2d for iek.gt.mek
         if(iek.gt.mek) nocalcstr=0

         if(bulk) then
c    always calculate SSPO 
           sgf=.true.
         else
c    calculate SSPO only for the first scf and CPA iteration
           sgf=itcpa.eq.1.and.itscf.eq.1
c    do not store SSPO for iek.gt.mek           
           sgf=sgf.or.(iek.gt.mek)
         end if
         sgfy=sgf
         sgfx=sgf
c    if vacuum, recalculate R-SSPO 
         if(rightm.eq.'V') sgfx=itcpa.eq.1
         if(slab) then
            sgfx=itcpa.eq.1
            sgfy=itcpa.eq.1
         end if
c
c ********************************************
c k-resolved layer-indexed structure constants
c ********************************************
c
         if(nocalcstr.ne.0) then
           call getg2d(g2d(1,1,-minprc,1,0,iekuse),gg2d,
     >                 bulkgeo,lmax,nintfc)
         else
c          -------------------------------------------------
           call gstore(park,ce,eta,bulkgeo,'L',lmax,nintfc,gg2d)
c          -------------------------------------------------
           call getg2d(gg2d,g2d(1,1,-minprc,1,0,iekuse),
     >                 bulkgeo,lmax,nintfc)
         end if
c
c **********************
c surface Green function
c **********************
c
         if(sgfy) then
            if(.not.bulkgeo) then
            call gstore(park,ce,eta,.true.,'L',lmax,nintfc,gg2dbulk)
c           ------------------------------------------
            call dugo(lmax,nl2,tminvl,lmmaxp,
     >                gg2dbulk(1,1,-minprc,1,0),
     >                gg2dbulk(1,1,-minprc,1,1),
     >                yy,mdim,'L',.true.,irel,sgfconv)
c           ------------------------------------------
            else
c           ------------------------------------------
            call dugo(lmax,nl2,tminvl,lmmaxp,
     >                gg2d(1,1,-minprc,1,0),
     >                gg2d(1,1,-minprc,1,1),
     >                yy,mdim,'L',.true.,irel,sgfconv)
c           ------------------------------------------
            end if
c
            if(.not.sgfconv) then
             write(6,'('' Left SGF'')')
             write(6,'('' k='',2f15.6)') park
             call flush(6)
             stop ' stopping due to SGF'
            end if
c
            call repl(ysave(1,1,iekuse,nsp),yy,ndiml,mdim)
         else
            call repl(yy,ysave(1,1,iekuse,nsp),ndiml,mdim)
         end if
c
         if(sgfx) then
c           if(rightm.eq.'V') then
c             call czero(xx,mdim*mdim)
c           else
            if(.not.bulkgeo) then
            call gstore(park,ce,eta,.true.,'R',lmax,nintfc,gg2dbulk)
c           --------------------------------------------
            call dugo(lmax,nl2,tminvr,lmmaxp,
     >                gg2dbulk(1,1,-minprc,1,2),
     >                gg2dbulk(1,1,-minprc,1,1),
     >                xx,mdim,'R',.true.,irel,sgfconv)
c           --------------------------------------------
            else
c           --------------------------------------------
            call dugo(lmax,nl2,tminvr,lmmaxp,
     >                gg2d(1,1,-minprc,1,2),
     >                gg2d(1,1,-minprc,1,1),
     >                xx,mdim,'R',.true.,irel,sgfconv)
c           --------------------------------------------
            end if
            if(.not.sgfconv) then
             write(6,'('' Right SGF'')')
             write(6,'('' k='',2f15.6)') park
             call flush(6)
             stop ' stopping due to SGF'
            end if
c
c           end if
            call repl(xsave(1,1,iekuse,nsp),xx,ndimr,mdim)
         else
            call repl(xx,xsave(1,1,iekuse,nsp),ndimr,mdim)
         end if 
c
c ***********************************
c k-resolved tau matrix for interface
c ***********************************
c
          wgeff=wk(ik)/ngeff
          cwgeff=dcmplx(wgeff,0.d0)
          w0=1.0d0
c        -------------------------------------------------
         call tau2d(irel,w0,bulkgeo,lmax,nl2,lmmaxp,nintfc,
     >              tminv,gg2d,xx,yy,tau00,mdim)
c        -------------------------------------------------
c
c **************************************
c BZ sum for layer diagonal tau matrices 
c **************************************
c
         do ig=1,ngeff
           ig1=invg(ig)
           do li=1,nintfc
             call repl(tau1,tau00(1,1,li),nl2,lmmaxp)
             call tripmt(dg(1,1,ig),tau1,dg1(1,1,ig),nl2,nl2,lmmaxp)
             call addmatc(tau(1,1,li),tau1,cwgeff,nl2,lmmaxp)
             call phystau(tau1,tminv(1,1,li),tminv(1,1,li),
     >                    alphaintkkr(0,li),lmax,1)
             do l=0,lmax
               lm1=l*l+1
               lm2=(l+1)*(l+1)
               taukdiag(ik,l,li)=(0.d0,0.d0)
               do lm=lm1,lm2
                 taukdiag(ik,l,li)=taukdiag(ik,l,li)+tau1(lm,lm)
               end do
             enddo
           end do
         enddo
c 33     continue
      end if
! >>> MPI if END
!   **********
c
      end do
c
c *************************
c end of loop over k points
c ************************* 
c
      call sublattmat(nintfc,iesublatt,nl2,lmmaxp,tau)
c
      taujoin = (0.d0,0.d0)
      call mpi_allreduce(tau,taujoin,lmmaxp*lmmaxp*mintfc,
     >mpi_double_complex,mpi_sum,mpi_comm_world,ierr)
      tau = taujoin
c        
c
c **********
c CPA solver
c **********
c     --------------------------------------------------------------
      call cpacor(conc,tminv,tminva,tminvb,tau,taua,taub,nl2,nintfc,
     >            itcpa,itcpam,cpatol,cpacon,concpa,cpatest,cpaerr)
c     --------------------------------------------------------------
c
      if(itest.ge.2) then
      do li=1,nintfc
        cpalay=(1.d0-conc(li)).gt.tiny
        if(cpalay) then
         write(6,*) ' Layer:',li
         write(6,*) ' M_A'
         call outmat1(tminva(1,1,li),nl2,nl2,lmmaxp,outtol,6)
         write(6,*) ' M_B'
         call outmat1(tminvb(1,1,li),nl2,nl2,lmmaxp,outtol,6)
         write(6,*) ' M_C'
         call outmat1(tminv(1,1,li),nl2,nl2,lmmaxp,outtol,6)
         write(6,*) ' TAU_C'
         call outmat1(tau(1,1,li),nl2,nl2,lmmaxp,outtol,6)
         write(6,*) ' TAU_A'
         call outmat1(taua(1,1,li),nl2,nl2,lmmaxp,outtol,6)
        end if
      end do
      end if
c
      if(concpa.and.cpacon) goto 100
      if(itest.ge.1) write(6,'('' CPA iterations:'',i3,2x,
     >'' ERROR:'',d15.6)') itcpa-1,cpaerr
c
c                        *************************
c                        * ending CPA iterations *
c                        *************************
c
c *************************************************************
c * loop over layers to transform layer diagonal tau matrices *
c * into physical representation                              *
c *************************************************************
c                                                                 
c
      do li=1,nintfc
         cpalay=(1.d0-conc(li)).gt.tiny
c       
         call phystau(taua(1,1,li),tminva(1,1,li),tminva(1,1,li),
     >                alphaintkkr(0,li),lmax,1)
c
         if(itest.ge.2) then
           write(6,'('' Tau -A: Layer'',i2)') li
           call outmat1(taua(1,1,li),nl2,nl2,lmmaxp,outtol,6)
         end if
c
         if(cpalay) then
c+------------+
c+ BIG CPA IF +
c+------------+
c
         call phystau(taub(1,1,li),tminvb(1,1,li),tminvb(1,1,li),
     >                alphaintkkr(0,li),lmax,1)
c
         if(itest.ge.2) then
           write(6,'('' Tau -B: Layer'',i2)') li
           call outmat1(taub(1,1,li),nl2,nl2,lmmaxp,outtol,6)
         end if
c
         end if
c+----------------+
c+ END BIG CPA IF +
c+----------------+ 
c
      end do
c * end loop over layers *
c
c
      return
      end
