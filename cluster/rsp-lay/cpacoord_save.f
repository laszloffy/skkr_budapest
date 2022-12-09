c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
c Store E- and k- for magnetic bulk
      subroutine cpacoord(
c     ====================
     > itscf,ie,ce,lmax,nintfc,eta,rightm,bulk,bulkgeo,wrel,
     > kset,xk,wk,nk,intbz,iek,
     > conc,itcpam,cpatol,cpatest,
     > dmata,dmatb,dmatpa,dmatpb,
     > tminvl,tminvr,tminv,tminva,tminvb,taua,taub,gtaua,gtaub,
c --------- cluster --------------------)
     > kpair,npair,tau_ii,tau_ij,tau_ji,taupath,
     > calcoff,offready)
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
#ifdef MPIP
      include 'mpif.h'
#endif
      parameter(mdim=mdimr) 
      parameter(mek=mekr) 
c
      logical tautest
      logical wrel,bulk,bulkgeo,sgf,sgfx,sgfy,sgfconv
      logical cpacon,concpa,cpatest,cpalay
c
      character*1 rightm
c
      real*8 park(2)
      real*8 xk(mkpar,2)
      real*8 wk(mkpar)
      real*8 wgeff
c
      integer invg(melem)
      integer ieqg(melem)
      integer ieqgl(melem)
      integer ieqgr(melem)
      integer igordl(melem)
      integer igordr(melem)
c
      real*8 conc(mintfc)
c
      complex*16 ce
      complex*16 psq
      complex*16 fac
c
      complex*16 tminv(kmymaxp,kmymaxp,mintfc)
      complex*16 tminva(kmymaxp,kmymaxp,mintfc)
      complex*16 tminvb(kmymaxp,kmymaxp,mintfc)
      complex*16 tminvl(kmymaxp,kmymaxp,minprc)
      complex*16 tminvr(kmymaxp,kmymaxp,minprc)
      complex*16 ttmp(kmymaxp,kmymaxp,mintfc,melem)
      complex*16 ttmpl(kmymaxp,kmymaxp,minprc,melem)
      complex*16 ttmpr(kmymaxp,kmymaxp,minprc,melem)
c
      complex*16 xx(mdim,mdim),yy(mdim,mdim),
     > gg2d(lmmaxp,lmmaxp,-minprc:minprc,minprc,0:mprc1+1),
     > gg2dbulk(lmmaxp,lmmaxp,-minprc:minprc,minprc,0:2)
c
      complex*16 xsave(mdim,mdim,mek+1,melem1),
     >           ysave(mdim,mdim,mek+1,melem1),
     > g2d(lmmaxp,lmmaxp,-minprc:minprc,minprc,0:mprc1+1,mek+1)
c
      complex*16 tau(kmymaxp,kmymaxp,mintfc)
      complex*16 taua(kmymaxp,kmymaxp,mintfc)
      complex*16 taub(kmymaxp,kmymaxp,mintfc)
      complex*16 gtaua(kmymaxp,kmymaxp,mintfc)
      complex*16 gtaub(kmymaxp,kmymaxp,mintfc)
      complex*16 tau2(kmymaxp,kmymaxp)
c
c --- For the cluster code ---
c
      character*15 taupath
c
      logical cpaok
      logical offdiag
      logical calcoff
      logical offready
c
      integer kpair(4,mpair)
      integer npair
      integer ipair
      integer kk
      integer n0
c
      real*8 amat
      real*8 a0mat
      integer mg
      integer mg0
      common/rrotmat/amat(2,2,melem),a0mat(2,2,melem),mg,mg0
c
      real*8 a1,a2
      common/brav2d/a1(2),a2(2)
c
      real*8 dvec(2)
      real*8 bvec(2)
      real*8 tbvec(2)
c
      real*8 rvec(2)
      real*8 trvec(2)
      real*8 dcpar(2)
      real*8 kx
c
      complex*16 tmusch(kmymaxp,kmymaxp,mintfc)
      complex*16 tau1(kmymaxp,kmymaxp,mintfc,mintfc,melem)
      complex*16 tauij(kmymaxp,kmymaxp,mintfc,mintfc)
      complex*16 tau_ii(kmymaxp,kmymaxp,mintfc)
      complex*16 tau_ij(kmymaxp,kmymaxp,mpair)
      complex*16 tau_ji(kmymaxp,kmymaxp,mpair)
      complex*16 taudiag(kmymaxp)
      complex*16 ffacx
      complex*16 ffacxp
c ----------------------------
c
      complex*16 dmata(kmymaxp,kmymaxp,mintfc)
      complex*16 dmatpa(kmymaxp,kmymaxp,mintfc)
      complex*16 dmatb(kmymaxp,kmymaxp,mintfc)
      complex*16 dmatpb(kmymaxp,kmymaxp,mintfc)
c
      complex*16 rmat(kmymaxp,kmymaxp,melem)
      complex*16 rmatp(kmymaxp,kmymaxp,melem)
c
      complex*16 alphalkkr(0:lmaxp,minprc)
      complex*16 alpharkkr(0:lmaxp,minprc)
      complex*16 alphaintkkr(0:lmaxp,mintfc)
      common/scrpar/alphalkkr,alpharkkr,alphaintkkr
c
      common/test/itest
      common/relfac/fac
      common/lay2d/cvec(mtotal,3),nextra,nbulkl,nbulkr,
     &             nprc,ninprc(0:mprc+1)
c
      data tol/1.0d-8/ 
      data tiny/1.0d-6/
c
      save xsave,ysave,g2d
c
c MPI
      integer root, myrank, nprocs, ierr
      common/mpi/root,myrank,nprocs
#ifdef MPIP
      complex*16 taujoin(kmymaxp,kmymaxp,mintfc)
      complex*16, allocatable :: tau_ijsum(:,:,:) ! offdiag
      integer*4 AllocateStatus
#endif
c ********************
c initialize constants
c ********************
c
c---> c in rydberg units:
      c=274.072d0
      if(.not.wrel) then
        psq=ce+ce*ce/(c*c)
      else
        psq=ce
      end if
      fac = psq/ce
      tautest=.false.
      if(tautest) then
      write(6,*) '<cpacoord> : psq,ce= ',psq,ce
      write(6,*) '<cpacoord> : fac= ',fac
      end if
c
      irel=1
c
      nl=lmax+1
      kmymax=2*nl*nl
c
      ninprcl = ninprc(0)
      ninprcr = ninprc(nprc+1)
      ndiml=ninprcl*kmymax
      ndimr=ninprcr*kmymax
c
      nprc1=nprc
      if(bulkgeo) nprc1=1
      if(nprc1.gt.mprc1) stop ' increase mprc1 !!!'
c
c ****************************************************************
c initialize irreducible representations of point group operations
c ****************************************************************
c
c     --------------------------------------
      call sphbas(rmat,rmatp,ng,invg)
c     --------------------------------------
      if(intbz.eq.0.or.intbz.eq.2) ngeff=1
      if(intbz.eq.1) ngeff=ng
      if(ngeff.gt.melem1) then
        write(6,'(/'' Parameter melem1 is too small!!!'')')
        stop
      end if
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
c
      cpaok=.false.        ! is the cpa calculation convergent (+1 iteration)
      offdiag=.false.      ! are the off-diag element calculated?
c
  100 continue
      iek=iekstart
c
c     ----------------------------------------
      call czero(tau,kmymaxp*kmymaxp*mintfc)
      call czero(tau_ii,kmymaxp*kmymaxp*mintfc)
      call czero(tau_ij,kmymaxp*kmymaxp*mpair)
      call czero(tau_ji,kmymaxp*kmymaxp*mpair)
c     ----------------------------------------
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
            call repl(tminvl(1,1,li),tminv(1,1,li),kmymax,kmymaxp)
            call repl(tminvr(1,1,li),tminv(1,1,li),kmymax,kmymaxp)
          enddo
        else
c This is for only spectral-DOS calculation: force I t-matrices
c to be identical with the corresponding self-consistent L t-matrices
          do li=1,nintfc
            call repl(tminv(1,1,li),tminvl(1,1,li),kmymax,kmymaxp)
          end do
        end if
      end if
c
c *******************************************************************
c transform inverse t-matrices with respect to point group operations
c and find degeneracies
c *******************************************************************
c
c     ---------------------------------------------------
      call sorttmat(
     > lmax,nintfc,ninprcl,ninprcr,rmat,rmatp,ngeff,invg,
     > tminvl,tminvr,tminv,ttmpl,ttmpr,ttmp,
     > ieqgl,ieqgr,ieqg,igordl,igordr,0)
c     ---------------------------------------------------
c
c ****************************
c loop over k points in 2D IBZ
c ****************************
c
      do ik=1,nk
#ifdef MPIP
       if(modulo(ik-1,nprocs).eq.myrank) then
#endif
        park(1)=xk(ik,1)
        park(2)=xk(ik,2)
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
     >  write(6,'(/''WARNING: iek='',i3,'' greater than mek='',i3/
     >             ''         causes dramatic slow-down !'')') iek,mek
c
        if(bulk) then
c    calculate g2d only for the first CPA iteration
          nocalcstr=itcpa-1
        else
c    calculate g2d only for the first scf and CPA iteration
          nocalcstr=itscf+itcpa-2
        end if
c    do not store g2d for iek.gt.mek
        if(iek.gt.mek) nocalcstr=0
c
        if(bulk) then
c    always calculate SSPO
          sgf=.true.
        else
c    calculate SSPO only for the first scf and CPA iteration
c    or for iek.gt.mek
          sgf=(itcpa.eq.1.and.itscf.eq.1).or.(iek.gt.mek)
        end if
        sgfy=sgf
        sgfx=sgf
c    if vacuum, recalculate R SSPO after every tenth step
c        if(rightm.eq.'V')
c     >   sgfx=((mod(itscf,10).eq.1).and.(itcpa.eq.1)).or.(iek.gt.mek)
c    if vacuum, calculate R SSPO in every step
        if(rightm.eq.'V')
     >   sgfx=(itcpa.eq.1).or.(iek.gt.mek)
c
c ********************************************
c k-resolved layer-indexed structure constants
c ********************************************
c
        if(nocalcstr.ne.0) then
          call getg2d(g2d(1,1,-minprc,1,0,iekuse),gg2d,
     >                bulkgeo,lmax,nintfc)
        else
c         -------------------------------------------------
          call gstore(park,psq,eta,bulkgeo,'L',lmax,nintfc,gg2d) ! ce changed to psq laszloffy
c         -------------------------------------------------
          call getg2d(gg2d,g2d(1,1,-minprc,1,0,iekuse),
     >                bulkgeo,lmax,nintfc)
        end if
c
c ********************************
c loop over point group operations
c ********************************
c
        do ig=1,ngeff
          ige=ieqg(ig)
          igel=ieqgl(ig)
          iger=ieqgr(ig)
          ig1=invg(ig)
c
c ---------------------------------------------------------------------
c ----- Surface Green function for left side --------------------------
c ---------------------------------------------------------------------
          if((igel.eq.ig).and.sgfy) then
            if(.not.bulkgeo) then
            call gstore(park,psq,eta,.true.,'L',lmax,nintfc,gg2dbulk)   
            call dugo(lmax,kmymax,ttmpl(1,1,1,igel),kmymaxp,
     >                gg2dbulk(1,1,-minprc,1,0),
     >                gg2dbulk(1,1,-minprc,1,1),
     >                yy,mdim,'L',.true.,irel,sgfconv)
            else
            call dugo(lmax,kmymax,ttmpl(1,1,1,igel),kmymaxp,
     >                gg2d(1,1,-minprc,1,0),
     >                gg2d(1,1,-minprc,1,1),
     >                yy,mdim,'L',.true.,irel,sgfconv)
            end if
            if(.not.sgfconv) then
             write(6,'('' Left SGF'')')
             write(6,'('' k='',2f15.6)') park
             stop
            end if
c
           call repl(ysave(1,1,iekuse,igordl(igel)),yy,ndiml,mdim)
          else
           call repl(yy,ysave(1,1,iekuse,igordl(igel)),ndiml,mdim)
          end if 
c ---------------------------------------------------------------------
c ----- Surface Green function for right side --------------------------
c ---------------------------------------------------------------------
          if((iger.eq.ig).and.sgfx) then
            if(.not.bulkgeo) then
            call gstore(park,psq,eta,.true.,'R',lmax,nintfc,gg2dbulk)   
            call dugo(lmax,kmymax,ttmpr(1,1,1,iger),kmymaxp,
     >                gg2dbulk(1,1,-minprc,1,2),
     >                gg2dbulk(1,1,-minprc,1,1),
     >                xx,mdim,'R',.true.,irel,sgfconv)
            else 
            call dugo(lmax,kmymax,ttmpr(1,1,1,iger),kmymaxp,
     >                gg2d(1,1,-minprc,1,2),
     >                gg2d(1,1,-minprc,1,1),
     >                xx,mdim,'R',.true.,irel,sgfconv)
            end if
            if(.not.sgfconv) then
             write(6,'('' Right SGF'')')
             write(6,'('' k='',2f15.6)') park
             stop
            end if
c
            call repl(xsave(1,1,iekuse,igordr(iger)),xx,ndimr,mdim)
          else
            call repl(xx,xsave(1,1,iekuse,igordr(iger)),ndimr,mdim)
          end if 
c
          if(ige.eq.ig) then
c
c ***********************************
c k-resolved tau matrix for interface
c ***********************************
c
            wgeff = wk(ik)/ngeff
c           ----------------------------------------------------
c           call tau2d(irel,wgeff,bulkgeo,lmax,kmymax,kmymaxp,nintfc,
c    >                 ttmp(1,1,1,ig),gg2d,xx,yy,tau11(1,1,1,ig),mdim)
c           ----------------------------------------------------
       tautest=.false.
       if(tautest) then
       write(6,*) '<cpacoord> : before tau2dgod'
       write(6,*) '<cpacoord> : irel = ',irel
       write(6,*) '<cpacoord> : wgeff = ',wgeff
       write(6,*) '<cpacoord> : bulkgeo = ',bulkgeo
       write(6,*) '<cpacoord> : lmax = ',lmax
       write(6,*) '<cpacoord> : kmymax = ',kmymax
       write(6,*) '<cpacoord> : kmymaxp = ',kmymaxp
       write(6,*) '<cpacoord> : nintfc = ',nintfc
       write(6,*) '<cpacoord> : ttmp(1,1,1,ig) = ',ttmp(1,1,1,ig)
       write(6,*) '<cpacoord> : gg2d = ',gg2d(1:nl*nl,1:nl*nl,0,1,0)
       write(6,*) '<cpacoord> : xx = ',xx(1:kmymax,1:kmymax)
       write(6,*) '<cpacoord> : yy = ',yy(1:kmymax,1:kmymax)
       write(6,*) '<cpacoord> : tau1(1,1,1,1,ig) = ',
     >            wgeff*tau1(1:kmymax,1:kmymax,1,1,ig)
       write(6,*) '<cpacoord> : mdim = ',mdim
       end if
c           ----------------------------------------------------
            call tau2dgod(irel,wgeff,bulkgeo,lmax,kmymax,kmymaxp,nintfc,
     >                  ttmp(1,1,1,ig),gg2d,xx,yy,tau1(1,1,1,1,ig),mdim)
c           ----------------------------------------------------
       if(tautest) then
       write(6,*) '<cpacoord> : after tau2dgod'
       write(6,*) '<cpacoord> : irel = ',irel
       write(6,*) '<cpacoord> : wgeff = ',wgeff
       write(6,*) '<cpacoord> : bulkgeo = ',bulkgeo
       write(6,*) '<cpacoord> : lmax = ',lmax
       write(6,*) '<cpacoord> : kmymax = ',kmymax
       write(6,*) '<cpacoord> : kmymaxp = ',kmymaxp
       write(6,*) '<cpacoord> : nintfc = ',nintfc
       write(6,*) '<cpacoord> : ttmp(1,1,1,ig) = ',ttmp(1,1,1,ig)
       write(6,*) '<cpacoord> : gg2d = ',gg2d(1:nl*nl,1:nl*nl,0,1,0)
       write(6,*) '<cpacoord> : xx = ',xx(1:kmymax,1:kmymax)
       write(6,*) '<cpacoord> : yy = ',yy(1:kmymax,1:kmymax)
       write(6,*) '<cpacoord> : tau1(1,1,1,1,ig) = ',
     >            wgeff*tau1(1:kmymax,1:kmymax,1,1,ig)
       write(6,*) '<cpacoord> : mdim = ',mdim
       end if
       tautest=.false.
       if(tautest) then
         write(6,*) ' <cpacoord> : tautest, tauii ig=',ig
         call outmat1(wgeff*tau1(1:kmymax,1:kmymax,1,1,ig),
     >                kmymax,kmymax,kmymax,
     >                tol,6)
         write(6,*) ' <cpacoord> : tautest, tauij ig=',ig
         call outmat1(wgeff*tau1(1:kmymax,1:kmymax,1,2,ig),
     >                kmymax,kmymax,kmymax,
     >                tol,6)
       end if
c
            if(itest.gt.3) then
              write(6,*) '<cpacoord>: after tau2d_god' 
              do li=1,nintfc
                write(6,*) '<cpacoord>: TAU_before BZ sum layer=',li
                call replms(taudiag,tau1(1,1,li,li,ige),lmax)
                do kk=1,kmymax
                   write(6,'(i3,2d20.10)') kk,taudiag(kk)
                end do
              end do
            end if

c DEBUG BEGINN
c           write(6,*) '<cpacoord>: after tau2d_god' 
c           do li=3,3
c             write(6,*) '<cpacoord>: TAU_before BZ sum layer=',li
c             do kk=1,kmymax
c                write(6,'(i3,2d20.10)') kk,wgeff*tau1(kk,kk,li,li,ige)
c                write(6,'(i3,2d20.10)') kk,tau11(kk,kk,li,ige)
c             end do
c           end do
c DEBUG BEGINN
          end if
c
c **************************************
c BZ sum for layer diagonal tau matrices
c **************************************
c -old version beginn
c         do li=1,nintfc
c           --------------------------------------------------
c           call repl(tau2,tau11(1,1,li,ige),kmymax,kmymaxp)
c           --------------------------------------------------
c           call tripmt(rmatp(1,1,ig1),tau2,rmat(1,1,ig1),
c    >                  kmymax,kmymax,kmymaxp)        
c           -----------------------------------------------
c           call addmat(tau(1,1,li),tau2,kmymax,kmymaxp)
c           -----------------------------------------------
c         end do
c -old version end
          do li=1,nintfc
c           --------------------------------------------------
            call repl(tau2,tau1(1,1,li,li,ige),kmymax,kmymaxp)
c           --------------------------------------------------
            call tripmt(rmatp(1,1,ig1),tau2,rmat(1,1,ig1),
     >                  kmymax,kmymax,kmymaxp)        
c           -----------------------------------------------
            call addmata(tau(1,1,li),tau2,wgeff,kmymax,kmymaxp)
c           -----------------------------------------------
          end do
c
c ======================================================================
c ----------------- off diagonal elements of tau matrix ----------------
c
          if((cpaok.AND.calcoff).AND.(npair.gt.0)) then
c           write(6,*) '<cpacoord>: Calculation of off-diag elments' 
            tautest=.false.
            do li = 1,nintfc
             do lj = 1,nintfc
              if(tautest) then
               write(6,*) '<cpacoord> : tautest tauij1(:,:,li,lj,ige)
     > before rotation',
     >                   li,lj,ige
               call outmat1(wgeff*tau1(1:kmymax,1:kmymax,li,lj,ige),
     >                kmymax,kmymax,kmymax,tol,6)
              end if
              call repl(tauij(1,1,li,lj),tau1(1,1,li,lj,ige),
     >                    kmymax,kmymaxp)
              call tripmt(rmatp(1,1,ig1),tauij(1,1,li,lj),rmat(1,1,ig1),
     >                       kmymax,kmymax,kmymaxp)
              if(tautest) then
               write(6,*) '<cpacoord> : tautest tauij(:,:,li,lj,ige)
     > after rotation',
     >                   li,lj,ige
               call outmat1(wgeff*tauij(1:kmymax,1:kmymax,li,lj),
     >                 kmymax,kmymax,kmymax,tol,6)
              end if
             end do
            end do
c
            n0 = ninprc(0)*(nextra+1)
            do ipair = 1,npair
              li = kpair(3,ipair)
              lj = kpair(4,ipair)
c
              bvec(1)=cvec(li+n0,1)-cvec(lj+n0,1)
              bvec(2)=cvec(li+n0,2)-cvec(lj+n0,2)
              tbvec(1)=a0mat(1,1,ig1)*bvec(1)+a0mat(1,2,ig1)*bvec(2)
              tbvec(2)=a0mat(2,1,ig1)*bvec(1)+a0mat(2,2,ig1)*bvec(2)
c
              rvec(1)=kpair(1,ipair)*a1(1) + kpair(2,ipair)*a2(1)
              rvec(2)=kpair(1,ipair)*a1(2) + kpair(2,ipair)*a2(2)
c
              trvec(1)=a0mat(1,1,ig1)*rvec(1)+a0mat(1,2,ig1)*rvec(2)
              trvec(2)=a0mat(2,1,ig1)*rvec(1)+a0mat(2,2,ig1)*rvec(2)
c
              dcpar(1)=-trvec(1)-tbvec(1)+bvec(1)
              dcpar(2)=-trvec(2)-tbvec(2)+bvec(2)
c
              kx=park(1)*dcpar(1)+park(2)*dcpar(2)
              ffacx = dcmplx(dcos(kx),-dsin(kx))
              ffacxp= dcmplx(dcos(kx), dsin(kx))
              if(tautest) then
                write(6,*) '<cpacoord> : tautest ffacx, ffacxp'
                write(6,*) ffacx,ffacxp
              end if
c ======================================================================
c ORIGINAL
c             bvec(1) = kpair(1,ipair)*a1(1) + kpair(2,ipair)*a2(1) +
c    >                  (cvec(lj+n0,1) - cvec(li+n0,1))
c             bvec(2) = kpair(1,ipair)*a1(2) + kpair(2,ipair)*a2(2) +
c    >                  (cvec(lj+n0,2) - cvec(li+n0,2))
c             dvec(1) = a0mat(1,1,ig1)*bvec(1) + a0mat(1,2,ig1)*bvec(2)-
c    >                  (cvec(lj+n0,1)-cvec(li+n0,1))
c             dvec(2) = a0mat(2,1,ig1)*bvec(1) + a0mat(2,2,ig1)*bvec(2)-
c    >                  (cvec(lj+n0,2)-cvec(li+n0,2))
c
c             kx = park(1)*dvec(1) + park(2)*dvec(2)
c             ffacx = dcmplx(dcos(kx),-dsin(kx))
c             ffacxp= dcmplx(dcos(kx),dsin(kx))
c ORIGINAL
c ======================================================================
c
c
              do k1 = 1,kmymax
                do k2 = 1,kmymax
                  tau_ij(k1,k2,ipair) = tau_ij(k1,k2,ipair) +
     >                  ffacx*wgeff*tauij(k1,k2,li,lj)
                  tau_ji(k1,k2,ipair) = tau_ji(k1,k2,ipair) +
     >                  ffacxp*wgeff*tauij(k1,k2,lj,li)
                end do
              end do
            end do
            offdiag=.true.
            offready=.true.
c           write(6,*) '<cpacoord>: Off-diag elments ready' 
          end if
          tautest=.false.
c ----------------------------------------------------------------------
c ======================================================================
c
        end do
c *********************************
c end of loop over point operations
c *********************************
c
#ifdef MPIP
      end if
#endif
c MPI if END
      end do
c *************************
c end of loop over k points
c *************************
c
#ifdef MPIP
       taujoin = (0.d0,0.d0)
       call mpi_allreduce(tau,taujoin,kmymaxp*kmymaxp*mintfc,
     >     mpi_double_complex,mpi_sum,mpi_comm_world,ierr)
c
       tau = taujoin
c
       if(calcoff.and.(npair.gt.0)) then
         ALLOCATE ( tau_ijsum(kmymax,kmymax,npair),
     >                   STAT = AllocateStatus)
         call alloccheck(AllocateStatus,'tau_ijsum in cpacoord')
         tau_ijsum = (0.d0,0.d0)
         call mpi_allreduce(tau_ij(1:kmymax,1:kmymax,1:npair),
     >       tau_ijsum,kmymax*kmymax*npair,
     >       mpi_double_complex,mpi_sum,mpi_comm_world,ierr)
c
         tau_ij(1:kmymax,1:kmymax,1:npair) = tau_ijsum
c
         tau_ijsum = (0.d0,0.d0)
         call mpi_allreduce(tau_ji(1:kmymax,1:kmymax,1:npair),
     >       tau_ijsum,kmymax*kmymax*npair,
     >       mpi_double_complex,mpi_sum,mpi_comm_world,ierr)
c
         tau_ji(1:kmymax,1:kmymax,1:npair) = tau_ijsum
         deallocate(tau_ijsum,stat=AllocateStatus)
         call alloccheck(AllocateStatus,'tau_ijsum dealloc in cpacoord')
       end if
c
#endif
c **********
c CPA solver
c **********
c     -----------------------------------------------------------------
      call cpacor(conc,tminv,tminva,tminvb,tau,taua,taub,kmymax,nintfc,
     >            itcpa,itcpam,cpatol,cpacon,concpa,cpatest,cpaerr)
c     -----------------------------------------------------------------
c
c     if(concpa.and.cpacon) goto 100
      if((.NOT.cpaok).AND.concpa.and.cpacon) goto 100
c
      if(itest.ge.1.AND.(.NOT.cpaok)) write(6,'('' CPA iterations:'',
     >   i3,2x,'' ERROR:'',d15.6)') itcpa-1,cpaerr
c
      cpaok=.true.
      if((.NOT.offdiag).AND.calcoff) goto 100
c
      if(cpaok.AND.concpa.and.cpacon) then
        write(6,'('' <cpacoord>: WARNING cpa convergence again worse 
     >                          in the +1 cpa iteration'',
     >        i3,2x,'' ERROR:'',d15.6)') itcpa-1,cpaerr
      end if
c
      if(itest.ge.1) write(6,'('' +1 CPA iterations:'',i3,2x,
     >'' ERROR:'',d15.6)') itcpa-1,cpaerr
c
c                        *************************
c                        * ending CPA iterations *
c                        *************************
c
c ======================================================================
      tautest=.false.
      ipair=1
      if(tautest) then
         write(6,*) ' <cpacoord> : tautest, tauij before phystau2,
     > electron-electron, ipair=',ipair
         call outmat1(tau_ij(1,1,ipair),kmymax,kmymax,kmymaxp,
     >                tol,6)
         write(6,*) ' <cpacoord> : tautest, tauji before phystau2,
     > electron-electron, ipair=',ipair
         call outmat1(tau_ji(1,1,ipair),kmymax,kmymax,kmymaxp,
     >                tol,6)
      end if
c - the off-diagonal elements of tau
      if(calcoff) then
       if(npair.gt.0) then
        do ipair = 1,npair
          li = kpair(3,ipair)
          lj = kpair(4,ipair)
          idiag = 0
c         if(li.eq.lj) idiag=1
c         call phystau2(tau_ij,tminva(1,1,li),tminva(1,1,lj),
c    >                  alphaintkkr(0,li),alphaintkkr(0,lj),lmax,idiag)
          call phystau2(tau_ij(1,1,ipair),tminv(1,1,li),tminv(1,1,lj),
     >                  alphaintkkr(0,li),alphaintkkr(0,lj),lmax,idiag)
          call phystau2(tau_ji(1,1,ipair),tminv(1,1,lj),tminv(1,1,li),
     >                  alphaintkkr(0,lj),alphaintkkr(0,li),lmax,idiag)
        end do
       end if
c unscreening the tmatices
c
        do li=1,nintfc
          call unsct(tminv(1,1,li),tmusch(1,1,li),kmymax,
     >               alphaintkkr(0,li),li)
        end do
c
c - the diagonal elements of tau
        do li=1,nintfc
          call repl(tau_ii(1,1,li),tau(1,1,li),kmymax,kmymaxp)
c         --------------------------------------------------------
          call phystau(tau_ii(1,1,li),tminv(1,1,li),tminv(1,1,li),
     >                 alphaintkkr(0,li),lmax,1)
c         --------------------------------------------------------
        end do
c
c       write(6,*) '<cpacoord>:',taupath
        call tauhout(ie,ce,npair,nintfc,nk,kmymax,
     >               tau_ij,tau_ji,tau_ii,taupath)
        call thout(ie,ce,nintfc,kmymax,tmusch,taupath)
      end if
c
c
      tautest=.false.
      ipair=1
      if(tautest) then
         write(6,*) ' <cpacoord> : tautest, tauij after phystau2,
     > electron-electron, ipair=',ipair
         call outmat1(tau_ij(1,1,ipair),kmymax,kmymax,kmymaxp,
     >                tol,6)
         write(6,*) ' <cpacoord> : tautest, tauji after phystau2,
     > electron-electron, ipair=',ipair
         call outmat1(tau_ji(1,1,ipair),kmymax,kmymax,kmymaxp,
     >                tol,6)
      end if
c ======================================================================
c
c
c *************************************************************
c * loop over layers to transform layer diagonal tau matrices *
c * into physical representation and rotate to the local      *
c * frame of reference                                        *
c *************************************************************
c 
      if(ie.eq.1) then
      tautest=.false.
      li=1
      if(itest.gt.2.or.tautest) then
         write(6,'(/'' tau -A '',i2)') li
         call outmat1(taua(1,1,li),kmymax,kmymax,kmymaxp,tol,6)
      end if
      end if
      do li=1,nintfc
         cpalay=(1.d0-conc(li)).gt.tiny
c
c        --------------------------------------------------------
         call phystau(taua(1,1,li),tminva(1,1,li),tminva(1,1,li),
     >                alphaintkkr(0,li),lmax,1)
c        --------------------------------------------------------
c
c rotate tau-matrix to local frame of reference
c
         call repl(gtaua(1,1,li),taua(1,1,li),kmymax,kmymaxp)
c        ------------------------------------------------------
         call tripmt(dmatpa(1,1,li),taua(1,1,li),dmata(1,1,li),
     >               kmymax,kmymax,kmymaxp)
c        ------------------------------------------------------
c
         tautest=.false.
         if(itest.gt.2.or.tautest) then
            write(6,'(/'' tau -A '',i2)') li
            call outmat1(taua(1,1,li),kmymax,kmymax,kmymaxp,tol,6)
         end if
c
         if(cpalay) then
c+------------+
c+ BIG CPA IF +
c+------------+
c
c        --------------------------------------------------------
         call phystau(taub(1,1,li),tminvb(1,1,li),tminvb(1,1,li),
     >                alphaintkkr(0,li),lmax,1)
c        --------------------------------------------------------
c
c rotate tau-matrix to local frame of reference
c
         call repl(gtaub(1,1,li),taub(1,1,li),kmymax,kmymaxp)
c        ---------------------------------------------------
         call tripmt(dmatpb(1,1,li),taub(1,1,li),dmatb(1,1,li),
     >               kmymax,kmymax,kmymaxp)
c        ---------------------------------------------------
c
         if(itest.gt.2) then
            write(6,'(/'' tau -B '',i2)') li
            call outmat1(taub(1,1,li),kmymax,kmymax,kmymaxp,tol,6)
         end if
c
         end if
c+----------------+
c+ END BIG CPA IF +
c+----------------+
      end do
c *** end loop over layers ***
c
      return
      end
