c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
c Store k- for magnetic bulk
      subroutine cpacoord(
c     ====================
     > itscf,ie,ce,E_Fermi,lmax,nbogomax,
     > nintfc,eta,rightm,bulk,bulkgeo,wrel,
     > kset,xk,wk,nk,intbz,iek,
     > conc,itcpam,cpatol,cpatest,
     > dmata,dmatb,dmatpa,dmatpb,
     > tminvl,tminvr,tminv,tminva,tminvb,
     > tau,taua,taub,gtaua,gtaub,c_light,
c --------- cluster --------------------)
     > kpair,npair,npair0,npair1,tau_ii,tau_ij,tau_ji,taupath,
     > clucalc,hostprint)
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
#ifdef MPIP
      include 'mpif.h'
#endif
      parameter(mdim=mdimdbogo) 
c
      logical wrel,bulk,bulkgeo,sgf,sgfx,sgfy,sgfconv
      logical cpacon,concpa,cpatest,cpalay
      logical clucalc, hostprint, calcoff
c
      logical tautest
c
      character*1 rightm
c
      dimension park(2),xk(mkpar,2),wk(mkpar)
c
      dimension invg(melem),ieqg(melem),ieqgl(melem),ieqgr(melem)
      dimension igordl(melem),igordr(melem)
c
      dimension conc(mintfc)
c
      complex*16 ce,psq,fac,celoc,psqh,fach,celoch,detl
      complex*16 help(dbogomaxp,dbogomaxp),help1(dbogomaxp,dbogomaxp)
      complex*16 taudiag(kmymaxp)
c
      complex*16 tminv(dbogomaxp,dbogomaxp,mintfc)
      complex*16 tminva(dbogomaxp,dbogomaxp,mintfc)
      complex*16 tminvb(dbogomaxp,dbogomaxp,mintfc)
      complex*16 tminvl(dbogomaxp,dbogomaxp,minprc)
      complex*16 tminvr(dbogomaxp,dbogomaxp,minprc)
      complex*16 ttmp(dbogomaxp,dbogomaxp,mintfc,melem)
      complex*16 ttmpl(dbogomaxp,dbogomaxp,minprc,melem)
      complex*16 ttmpr(dbogomaxp,dbogomaxp,minprc,melem)
      complex*16 ctrafo(dbogomaxp,dbogomaxp)
      complex*16 ctrafo2(dbogomaxp,dbogomaxp)
c
      complex*16 xx(mdim,mdim),yy(mdim,mdim),
     > gg2d(bogomaxp,bogomaxp,-minprc:minprc,minprc,0:mprc1+1),
     > gg2dbulk(bogomaxp,bogomaxp,-minprc:minprc,minprc,0:2),
     > gg2de(lmmaxp,lmmaxp,-minprc:minprc,minprc,0:mprc1+1),
     > gg2dbulke(lmmaxp,lmmaxp,-minprc:minprc,minprc,0:2),
     > gg2dh(lmmaxp,lmmaxp,-minprc:minprc,minprc,0:mprc1+1),
     > gg2dbulkh(lmmaxp,lmmaxp,-minprc:minprc,minprc,0:2)
c
      complex*16 xsave(mdim,mdim,mkpar1+1,melem1),
     >           ysave(mdim,mdim,mkpar1+1,melem1),
     > g2de(lmmaxp,lmmaxp,-minprc:minprc,minprc,0:mprc1+1,mkpar1+1),
     > g2dh(lmmaxp,lmmaxp,-minprc:minprc,minprc,0:mprc1+1,mkpar1+1)
c
      complex*16 tau(dbogomaxp,dbogomaxp,mintfc)
      complex*16 taua(dbogomaxp,dbogomaxp,mintfc)
      complex*16 taub(dbogomaxp,dbogomaxp,mintfc)
      complex*16 gtaua(dbogomaxp,dbogomaxp,mintfc)
      complex*16 gtaub(dbogomaxp,dbogomaxp,mintfc)
      complex*16 tau1(dbogomaxp,dbogomaxp,mintfc,melem)
      complex*16 tau2(dbogomaxp,dbogomaxp)
c
c for cluster code
c
      character*300 taupath
c
c      logical cpaok
c      logical offdiag
c      logical calcoff
c      logical offready
c
      integer kpair(4,npair0)
      integer npair,npair0,npair1
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
c      real*8 dvec(2)
      real*8 bvec(2)
      real*8 tbvec(2)
c
      real*8 rvec(2)
      real*8 trvec(2)
      real*8 dcpar(2)
      real*8 kx
c cluster tau matrices
      integer*4  AllocateStatus,ie
      complex*16, allocatable :: tauij1(:,:,:,:,:) ! offdiag
      complex*16, allocatable :: tauij(:,:,:,:) ! offdiag
      complex*16, allocatable :: tau_ijsum(:,:,:) ! offdiag
      complex*16, allocatable :: tmusch(:,:,:)
      complex*16 tau_ii(nbogomax,nbogomax,nintfc)
      complex*16 tau_ij(nbogomax,nbogomax,npair1)
c     complex*16 tau_ijsum(dbogomaxp,dbogomaxp,mpair)
      complex*16 tau_ji(nbogomax,nbogomax,npair1)
      complex*16 ffacx
      complex*16 ffacxp
c
      complex*16 dmata(kmymaxp,kmymaxp,mintfc)
      complex*16 dmatpa(kmymaxp,kmymaxp,mintfc)
      complex*16 dmatb(kmymaxp,kmymaxp,mintfc)
      complex*16 dmatpb(kmymaxp,kmymaxp,mintfc)
c
      complex*16 rmat(dbogomaxp,dbogomaxp,melem)
      complex*16 rmatp(dbogomaxp,dbogomaxp,melem)
c
      complex*16 alphalkkr(0:lmaxp,minprc)
      complex*16 alpharkkr(0:lmaxp,minprc)
      complex*16 alphaintkkr(0:lmaxp,mintfc)
      complex*16 alphalkkrh(0:lmaxp,minprc)
      complex*16 alpharkkrh(0:lmaxp,minprc)
      complex*16 alphaintkkrh(0:lmaxp,mintfc)
c
      common/scrpar/alphalkkr,alpharkkr,alphaintkkr
      common/scrparh/alphalkkrh,alpharkkrh,alphaintkkrh
c
      common/test/itest
      common/lay2d/cvec(mtotal,3),nextra,nbulkl,nbulkr,
     &             nprc,ninprc(0:mprc+1)
      common/iterpar/errmax,itermaxl,itermaxr,ichk
c
      common/relfac/fac,fach
c
      data tol/1.0d-8/ 
      data outtol/1.0d-12/ 
      data tiny/1.0d-6/
c
      save xsave,ysave,g2d
c MPI
      integer root, myrank, nprocs, ierr
      common/mpiparam/root,myrank,nprocs
      complex*16 taujoin(dbogomaxp,dbogomaxp,mintfc)
c ********************
c initialize constants
c ********************
c
c           write(6,*) ' <cpacoord> : begin ie=',ie
c---> c in rydberg units:
c      c=274.072d0 -- from input
      c=c_light
c---> energy for holes (conjugation is componsated in gstore.f)
      celoch=E_Fermi-conjg(ce)
      psqh=celoch+celoch*celoch/(c*c)
      fach=conjg(psqh)/conjg(celoch)
c---> energy for electrons 
      celoc=E_Fermi+ce
      psq=celoc+celoc*celoc/(c*c)
      fac=psq/celoc
c     fac=(1d0,0d0) ! TODO: only for tests laszloffy
c
      irel=1
c
      nl=lmax+1
      lmmax=nl*nl
      kmymax=2*nl*nl
c
      ninprcl = ninprc(0)
      ninprcr = ninprc(nprc+1)
      ndiml=ninprcl*kmymax*2
      ndimr=ninprcr*kmymax*2
c
      call conjinkmy(ctrafo)
      ctrafo2=ctrafo
      call gjinv(ctrafo2,dbogomaxp,dbogomaxp,detl)
c
      nprc1=nprc
      if(bulkgeo) nprc1=1
      if(nprc1.gt.mprc1) stop ' increase mprc1 !!!'
c
c ****************************************************************
c initialize irreducible representations of point group operations
c ****************************************************************
c
c---> Extended for holes
c     --------------------------------------
      call sphbas(rmat,rmatp,ng,invg,kmymax)
c     --------------------------------------
      if(intbz.eq.0.or.intbz.eq.2) ngeff=1
      if(intbz.eq.1) ngeff=ng
c      write(6,*) 'ngeff in cpacord: ', ngeff
      if(ngeff.gt.melem1) then
        write(6,'(/'' Parameter melem1 is too small!!!'')')
        stop
      end if
c
c ****************************************************************
c allocate arrays for cluster calculations
c ****************************************************************
      calcoff=clucalc.or.hostprint
      if(calcoff) then
        ALLOCATE ( tauij1(nbogomax,nbogomax,nintfc,nintfc,ngeff),
     >                  STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >        'tauij1 in cpacoord                                ')
        ALLOCATE ( tauij(nbogomax,nbogomax,nintfc,nintfc),
     >                  STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >        'tauij in cpacoord                                 ')
      end if 
c        
c           write(6,*) ' <cpacoord> : alloc ie=',ie
c
c                        ***************************
c                        * starting CPA iterations *
c                        ***************************
c
      itcpa=1
      concpa=.true.        ! controls CPA tolerance
      cpacon=.true.        ! controls no. of CPA iterations
      tautest=.false.
c
      iekstart=iek
  100 continue
      iek=iekstart
c
c     ----------------------------------------
      call czero(tau,dbogomaxp*dbogomaxp*mintfc)
      if(calcoff) then
        call czero(tau_ii,nbogomax*nbogomax*nintfc)
        call czero(tau_ij,nbogomax*nbogomax*npair1)
        call czero(tau_ji,nbogomax*nbogomax*npair1)
      end if
c     ----------------------------------------
c
c *******************************
c adjust t-matrices for bulk case
c *******************************
c
      if(bulk) then
c The I region is supposed to comprise one PL only !
c        if(kset.ge.1) then
c Force L and R t-matrices to be identical with the
c corresponding t-matrices in I
          do li=1,ninprcl
            call repl(tminvl(1,1,li),tminv(1,1,li),nbogomax,dbogomaxp)
            call repl(tminvr(1,1,li),tminv(1,1,li),nbogomax,dbogomaxp)
          enddo
c TODO: bnyari,check!! Changed in tcpa_main, interfaces matrices are read
c        else
c This is for only spectral-DOS calculation: force I t-matrices
c to be identical with the corresponding self-consistent L t-matrices
c          do li=1,nintfc
c            call repl(tminv(1,1,li),tminvl(1,1,li),2*kmymax,dbogomaxp)
c          end do
c        end if
      end if
c
c *******************************************************************
c transform inverse t-matrices with respect to point group operations
c and find degeneracies
c *******************************************************************
c
c     do li=1,nintfc
c       write(6,*) '  screened tm1 matrix'
c       call outmat1(tminv(1,1,li),kmymax,kmymax,kmymaxp,tol,6)
c     end do
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
       if(modulo(ik-1,nprocs) == myrank) then
#endif
        park(1)=xk(ik,1)
        park(2)=xk(ik,2)
c        write(6,*) 'k point in cpacoord:', park
c
c *************************************************
c initialize saving of structure constants and SSPO
c *************************************************
c
        iek=iek+1
        if(ik.le.mkpar1) then
          ikuse=ik
        else
          ikuse=mkpar1+1
        end if
c
c    calculate g2d only for the first CPA iteration
        nocalcstr=itcpa-1
c    store g2d only for the first mk1 points
        if(ik.gt.mkpar1) nocalcstr=0
c
        if(bulk) then
c    always calculate SSPO
          sgf=.true.
        else
c    calculate SSPO only for the first scf and CPA iteration
c    or for iek.gt.mek
c         sgf=(itcpa.eq.1.and.itscf.eq.1).or.(ik.gt.mkpar1)
          sgf=(itcpa.eq.1).or.(ik.gt.mkpar1)
        end if
        sgfy=sgf
        sgfx=sgf
c    if vacuum, recalculate R SSPO after every tenth step
        if(rightm.eq.'V')
     >   sgfx=((mod(itscf,10).eq.1).and.(itcpa.eq.1)).or.(ik.gt.mkpar1)
c
c ********************************************
c k-resolved layer-indexed structure constants
c ********************************************
c           write(6,*) ' <cpacoord> before getg2d : ie=',ie
c
        if(nocalcstr.ne.0) then
          call getg2d(g2de(1,1,-minprc,1,0,ikuse),gg2de,
     >                bulkgeo,lmax,nintfc)
          call getg2d(g2dh(1,1,-minprc,1,0,ikuse),gg2dh,
     >                bulkgeo,lmax,nintfc)
c         --------------------------------------------------
          call g2d_pack(gg2d(:,:,:,:,:),gg2de(:,:,:,:,:),
     >                   gg2dh(:,:,:,:,:),lmmax,mprc1)
c         --------------------------------------------------
        else
c         --------------------------------------------------
          call gstore(park,psq,eta,bulkgeo,'L',lmax,nintfc,
     >                gg2de,.false.,c_light,mprc1)
          call gstore(park,psqh,eta,bulkgeo,'L',lmax,nintfc,
     >                gg2dh,.true.,c_light,mprc1)
c         --------------------------------------------------
          call getg2d(gg2de,g2de(1,1,-minprc,1,0,ikuse),
     >                bulkgeo,lmax,nintfc)
          call getg2d(gg2dh,g2dh(1,1,-minprc,1,0,ikuse),
     >                 bulkgeo,lmax,nintfc)
c          -------------------------------------------------
          call g2d_pack(gg2d(:,:,:,:,:),gg2de(:,:,:,:,:),
     >                   gg2dh(:,:,:,:,:),lmmax,mprc1)
        endif
      tautest=.false.
      if(tautest.and.myrank.eq.root) then
       write(6,*) '<cpacoord> : psq = ',psq
       write(6,*) '<cpacoord> : gg2d pack = '
         call outmat1(gg2d(1:nl*nl,1:nl*nl,0,1,0),nl*nl,nl*nl,nl*nl,
     >     tol,6)
       write(6,*) '<cpacoord> : xx = '
         call outmat1(xx(1:kmymax,1:kmymax),kmymax,kmymax,kmymax,tol,6)
       write(6,*) '<cpacoord> : yy = '
         call outmat1(yy(1:kmymax,1:kmymax),kmymax,kmymax,kmymax,tol,6)
      end if
c
c           write(6,*) ' <cpacoord> : after getg2d ie=',ie
c ********************************
c loop over point group operations
c ********************************
c
        do ig=1,ngeff
          ige=ieqg(ig)
          igel=ieqgl(ig)
          iger=ieqgr(ig)
          ig1=invg(ig)
          tautest=.false.
          if(tautest.and.myrank.eq.root) then
            write(6,*) '<cpacoord> : ig= ',ig
            write(6,*) '<cpacoord> : ige= ',ige
            write(6,*) '<cpacoord> : igel= ',igel
            write(6,*) '<cpacoord> : iger= ',iger
            write(6,*) '<cpacoord> : ig1= ',ig1
          end if
c
c **********************
c surface Green function
c **********************
c
          if((igel.eq.ig).and.sgfy) then
          if(ichk.gt.0) then
            write(6,'(''L-SGF Point group'',i3,'' K-point'',i5,2f11.5)')
     >         ig,ik,park(1),park(2)
          endif
          if(.not.bulkgeo) then
          call gstore(park,psq,eta,.true.,'L',lmax,nintfc,
     >                gg2dbulke,.false.,c_light,1)
          call gstore(park,psqh,eta,.true.,'L',lmax,nintfc,
     >                gg2dbulkh,.true.,c_light,1)
c           ------------------------------------------
          call g2d_pack(gg2dbulk(:,:,:,:,:),gg2dbulke(:,:,:,:,:),
     >                  gg2dbulkh(:,:,:,:,:),lmmax,1)
c         ------------------------------------------------
          tautest=.false.
c         if(tautest.and.myrank.eq.root) then
          if(tautest.and.myrank.eq.10) then
            write(6,*) '<cpacoord> : lmax= ',lmax
            write(6,*) '<cpacoord> : nbogomax= ',nbogomax
            write(6,*) '<cpacoord> : ttmpl='
            call outmat1(ttmpl(1:nbogomax,1:nbogomax,1,igel),
     >                   nbogomax,nbogomax,nbogomax,tol,6)
            write(6,*) '<cpacoord> : dbogomaxp= ',dbogomaxp
            write(6,*) '<cpacoord> : gg2dbulk0= '
            call outmat1(gg2dbulk(1:bogomaxp,1:bogomaxp,-1,1,0),
     >                   bogomaxp,bogomaxp,bogomaxp,tol,6)
            write(6,*) '<cpacoord> : gg2dbulk1= '
            call outmat1(gg2dbulk(1:bogomaxp,1:bogomaxp,-1,1,1),
     >                   bogomaxp,bogomaxp,bogomaxp,tol,6)
            write(6,*) '<cpacoord> : yy = '
          call outmat1(yy(1:kmymax,1:kmymax),kmymax,kmymax,kmymax,tol,6)
            write(6,*) '<cpacoord> : mdim= ',mdim
            write(6,*) '<cpacoord> : irel= ',irel
c           write(6,*) '<cpacoord> : sgfconv= ',sgfconv
c          else
c            call sleep(10)
          end if
          call dugo(lmax,nbogomax,ttmpl(1,1,1,igel),dbogomaxp,
     >              gg2dbulk(1,1,-minprc,1,0),
     >              gg2dbulk(1,1,-minprc,1,1),
     >              yy,mdim,'L',.true.,irel,sgfconv)
c         ------------------------------------------------
c         if(tautest.and.myrank.eq.root) then
          if(tautest.and.myrank.eq.10) then
            write(6,*) '<cpacoord> : lmax= ',lmax
            write(6,*) '<cpacoord> : nbogomax= ',nbogomax
            write(6,*) '<cpacoord> : ttmpl='
            call outmat1(ttmpl(1:nbogomax,1:nbogomax,1,igel),
     >                   nbogomax,nbogomax,nbogomax,tol,6)
            write(6,*) '<cpacoord> : dbogomaxp= ',dbogomaxp
            write(6,*) '<cpacoord> : gg2dbulk0= '
            call outmat1(gg2dbulk(1:bogomaxp,1:bogomaxp,-1,1,0),
     >                   bogomaxp,bogomaxp,bogomaxp,tol,6)
            write(6,*) '<cpacoord> : gg2dbulk1= '
            call outmat1(gg2dbulk(1:bogomaxp,1:bogomaxp,-1,1,1),
     >                   bogomaxp,bogomaxp,bogomaxp,tol,6)
            write(6,*) '<cpacoord> : yy = '
          call outmat1(yy(1:kmymax,1:kmymax),kmymax,kmymax,kmymax,tol,6)
            write(6,*) '<cpacoord> : mdim= ',mdim
            write(6,*) '<cpacoord> : irel= ',irel
            write(6,*) '<cpacoord> : sgfconv= ',sgfconv
          end if
          else
c         ------------------------------------------------
          call dugo(lmax,nbogomax,ttmpl(1,1,1,igel),dbogomaxp,
     >              gg2d(1,1,-minprc,1,0),
     >              gg2d(1,1,-minprc,1,1),
     >              yy,mdim,'L',.true.,irel,sgfconv)
c         ------------------------------------------------
          end if
          if(.not.sgfconv) then
             write(6,'('' Left SGF'')')
             write(6,'('' k='',2f15.6)') park
c#ifdef MPIP
c             call flush(6)
c             call mpi_finalize(ierror)
c#endif
             stop
          end if
c
          call repl(ysave(1,1,ikuse,igordl(igel)),yy,
     >              ndiml,mdim)
          else
            call repl(yy,ysave(1,1,ikuse,igordl(igel)),
     >                ndiml,mdim)
          endif
      tautest=.false.
      if(tautest.and.myrank.eq.root) then
       write(6,*) '<cpacoord> : psq = ',psq
       write(6,*) '<cpacoord> : gg2d dugo 1 = '
         call outmat1(gg2d(1:nl*nl,1:nl*nl,0,1,0),nl*nl,nl*nl,nl*nl,
     >     tol,6)
       write(6,*) '<cpacoord> : xx = '
         call outmat1(xx(1:kmymax,1:kmymax),kmymax,kmymax,kmymax,tol,6)
       write(6,*) '<cpacoord> : yy = '
         call outmat1(yy(1:kmymax,1:kmymax),kmymax,kmymax,kmymax,tol,6)
      end if
c
          if((iger.eq.ig).and.sgfx) then
          if(ichk.gt.0) then
          write(6,'(''R-SGF Point group'',i3,'' K-point'',i5,2f11.5)')
     >    ig,ik,park(1),park(2)
          endif
          if(.not.bulkgeo) then
          call gstore(park,psq,eta,.true.,'R',lmax,nintfc,
     >                gg2dbulke,.false.,c_light,1)
          call gstore(park,psqh,eta,.true.,'R',lmax,nintfc,
     >                gg2dbulkh,.true.,c_light,1)
c           ------------------------------------------
          call g2d_pack(gg2dbulk(:,:,:,:,:),gg2dbulke(:,:,:,:,:),
     >                  gg2dbulkh(:,:,:,:,:),lmmax,1)
c         ------------------------------------------------
          call dugo(lmax,nbogomax,ttmpr(1,1,1,iger),dbogomaxp,
     >              gg2dbulk(1,1,-minprc,1,2),
     >              gg2dbulk(1,1,-minprc,1,1),
     >              xx,mdim,'R',.true.,irel,sgfconv)
c         ------------------------------------------------
          else
c         ------------------------------------------------
          call dugo(lmax,nbogomax,ttmpr(1,1,1,iger),dbogomaxp,
     >              gg2d(1,1,-minprc,1,2),
     >              gg2d(1,1,-minprc,1,1),
     >              xx,mdim,'R',.true.,irel,sgfconv)
c         ------------------------------------------------
          end if
          if(.not.sgfconv) then
             write(6,'('' Right SGF'')')
             write(6,'('' k='',2f15.6)') park
             stop
          end if
c           write(6,*) ' <cpacoord> : dugo ie=',ie
      tautest=.false.
      if(tautest.and.myrank.eq.root) then
       write(6,*) '<cpacoord> : psq = ',psq
       write(6,*) '<cpacoord> : gg2d dugo 2 = '
         call outmat1(gg2d(1:nl*nl,1:nl*nl,0,1,0),nl*nl,nl*nl,nl*nl,
     >     tol,6)
       write(6,*) '<cpacoord> : xx = '
         call outmat1(xx(1:kmymax,1:kmymax),kmymax,kmymax,kmymax,tol,6)
       write(6,*) '<cpacoord> : yy = '
         call outmat1(yy(1:kmymax,1:kmymax),kmymax,kmymax,kmymax,tol,6)
      end if
      tautest=.false.
      if(tautest.and.myrank.eq.root) then
        write(6,*) '<cpacoord> : xsave, ikuse=',ikuse
        write(6,*) '<cpacoord> : xsave, igordr(iger)=',igordr(iger)
        write(6,*) '<cpacoord> : xsave, ndimr=',ndimr
        write(6,*) '<cpacoord> : xsave, mdir=',mdir
        write(6,*) '<cpacoord> : xsave, mkpar1+1=',mkpar1+1
        write(6,*) '<cpacoord> : xsave, melem1=',melem1
      end if
c
c         call repl(xsave(1,1,ikuse,igordr(iger)),xx,
c    >              ndimr,mdim)
c         else
c         call repl(xx,xsave(1,1,ikuse,igordr(iger)),
c    >              ndimr,mdim)
          endif
c
      tautest=.false.
      if(tautest.and.myrank.eq.root) then
        write(6,*) '<cpacoord> : xsave, ikuse=',ikuse
        write(6,*) '<cpacoord> : xsave, igordr(iger)=',igordr(iger)
        write(6,*) '<cpacoord> : xsave, ndimr=',ndimr
        write(6,*) '<cpacoord> : xsave, mdir=',mdir
        write(6,*) '<cpacoord> : xsave, mkpar1+1=',mkpar1+1
        write(6,*) '<cpacoord> : xsave, melem1=',melem1
      end if
      tautest=.false.
      if(tautest.and.myrank.eq.root) then
       write(6,*) '<cpacoord> : psq = ',psq
       write(6,*) '<cpacoord> : gg2d xpack = '
         call outmat1(gg2d(1:nl*nl,1:nl*nl,0,1,0),nl*nl,nl*nl,nl*nl,
     >     tol,6)
       write(6,*) '<cpacoord> : xx = '
         call outmat1(xx(1:kmymax,1:kmymax),kmymax,kmymax,kmymax,tol,6)
       write(6,*) '<cpacoord> : yy = '
         call outmat1(yy(1:kmymax,1:kmymax),kmymax,kmymax,kmymax,tol,6)
      end if
c
          if(ige.eq.ig) then
c
c ***********************************
c k-resolved tau matrix for interface
c ***********************************
c
c          if(itest.ge.4) then
c          write(6,*) ' YY'
c          do il=1,ninprcl
c          do jl=1,ninprcl
c            i0=(il-1)*(2*kmymax)
c            j0=(jl-1)*(2*kmymax)
c            do i=1,(2*kmymax)
c            do j=1,(2*kmymax)
c              help1(i,j)=yy(i0+i,j0+j)
c            end do
c            end do
c            write(6,*) il,jl
c            call replmsf(help,help1,lmax)
c            call outmat1(help,2*kmymax,2*kmymax,dbogomaxp,outtol,6)
c          end do
c          end do
c          write(6,*) ' XX'
c          do il=1,ninprcr
c          do jl=1,ninprcr
c            i0=(il-1)*(2*kmymax)
c            j0=(jl-1)*(2*kmymax)
c            do i=1,(2*kmymax)
c            do j=1,(2*kmymax)
c              help1(i,j)=xx(i0+i,j0+j)
c            end do
c            end do
c            write(6,*) il,jl
c            call replmsf(help,help1,lmax)
c            call outmat1(help,2*kmymax,2*kmymax,dbogomaxp,outtol,6)
c          end do
c          end do
c          end if
c
      tautest=.false.
       if(tautest.and.myrank.eq.root) then
       write(6,*) '<cpacoord> : psq = ',psq
       write(6,*) '<cpacoord> : gg2d ige = '
         call outmat1(gg2d(1:nl*nl,1:nl*nl,0,1,0),nl*nl,nl*nl,nl*nl,
     >     tol,6)
       write(6,*) '<cpacoord> : xx = '
         call outmat1(xx(1:kmymax,1:kmymax),kmymax,kmymax,kmymax,tol,6)
       write(6,*) '<cpacoord> : yy = '
         call outmat1(yy(1:kmymax,1:kmymax),kmymax,kmymax,kmymax,tol,6)
      end if
          wgeff = wk(ik)/ngeff
      tautest=.false.
       if(tautest.and.myrank.eq.root) then
       write(6,*) '<cpacoord> : psq = ',psq
       write(6,*) '<cpacoord> : gg2d wgeff = '
         call outmat1(gg2d(1:nl*nl,1:nl*nl,0,1,0),nl*nl,nl*nl,nl*nl,
     >     tol,6)
       write(6,*) '<cpacoord> : xx = '
         call outmat1(xx(1:kmymax,1:kmymax),kmymax,kmymax,kmymax,tol,6)
       write(6,*) '<cpacoord> : yy = '
         call outmat1(yy(1:kmymax,1:kmymax),kmymax,kmymax,kmymax,tol,6)
      end if
          if(calcoff) then
       tautest=.false.
       if(tautest.and.myrank.eq.root) then
       write(6,*) '<cpacoord> : before tau2dgod'
       write(6,*) '<cpacoord> : irel = ',irel
       write(6,*) '<cpacoord> : wgeff = ',wgeff
       write(6,*) '<cpacoord> : bulkgeo = ',bulkgeo
       write(6,*) '<cpacoord> : lmax = ',lmax
       write(6,*) '<cpacoord> : kmymax = ',nbogomax
       write(6,*) '<cpacoord> : kmymaxp = ',dbogomaxp
       write(6,*) '<cpacoord> : nintfc = ',nintfc
       write(6,*) '<cpacoord> : ttmp(1,1,1,ig) = ',ttmp(1,1,1,ig)
       write(6,*) '<cpacoord> : gg2d = '
         call outmat1(gg2d(1:nl*nl,1:nl*nl,0,1,0),nl*nl,nl*nl,nl*nl,
     >     tol,6)
       write(6,*) '<cpacoord> : xx = '
         call outmat1(xx(1:kmymax,1:kmymax),kmymax,kmymax,kmymax,tol,6)
       write(6,*) '<cpacoord> : yy = '
         call outmat1(yy(1:kmymax,1:kmymax),kmymax,kmymax,kmymax,tol,6)
       write(6,*) '<cpacoord> : tau1(1,1,1,1,ig) = ',tauij1(1,1,1,1,ig)
       write(6,*) '<cpacoord> : mdim = ',mdim
       end if
c           ----------------------------------------------------
            call tau2dgod(irel,wgeff,bulkgeo,lmax,nbogomax,dbogomaxp,
     >        nintfc,ttmp(1,1,1,ig),gg2d,xx,yy,tauij1(1,1,1,1,ig),mdim)
c           ----------------------------------------------------
       if(tautest.and.myrank.eq.root) then
       write(6,*) '<cpacoord> : after tau2dgod'
       write(6,*) '<cpacoord> : irel = ',irel
       write(6,*) '<cpacoord> : wgeff = ',wgeff
       write(6,*) '<cpacoord> : bulkgeo = ',bulkgeo
       write(6,*) '<cpacoord> : lmax = ',lmax
       write(6,*) '<cpacoord> : kmymax = ',nbogomax
       write(6,*) '<cpacoord> : kmymaxp = ',dbogomaxp
       write(6,*) '<cpacoord> : nintfc = ',nintfc
       write(6,*) '<cpacoord> : ttmp(1,1,1,ig) = ',ttmp(1,1,1,ig)
       write(6,*) '<cpacoord> : gg2d = ',gg2d(1:nl*nl,1:nl*nl,0,1,0)
       write(6,*) '<cpacoord> : xx = '
         call outmat1(xx(1:kmymax,1:kmymax),kmymax,kmymax,kmymax,tol,6)
       write(6,*) '<cpacoord> : yy = '
         call outmat1(yy(1:kmymax,1:kmymax),kmymax,kmymax,kmymax,tol,6)
       write(6,*) '<cpacoord> : tau1(1,1,1,1,ig) = ',tauij1(1,1,1,1,ig)
       write(6,*) '<cpacoord> : mdim = ',mdim
       end if
       tautest=.false.
       if(tautest.and.myrank.eq.root) then
         write(6,*) ' <cpacoord> : tautest, tauii ig=',ig
         call outmat1(tauij1(1,1,1,1,ig),nbogomax,nbogomax,nbogomax,
     >                tol,6)
         write(6,*) ' <cpacoord> : tautest, tauij ig=',ig
         call outmat1(tauij1(1,1,1,2,ig),nbogomax,nbogomax,nbogomax,
     >                tol,6)
       end if
          else
      tautest=.false.
      if(tautest.and.myrank.eq.root) then
       write(6,*) '<cpacoord> : psq = ',psq
       write(6,*) '<cpacoord> : gg2d tau2d = '
         call outmat1(gg2d(1:nl*nl,1:nl*nl,0,1,0),nl*nl,nl*nl,nl*nl,
     >     tol,6)
      end if
       tautest=.false.
       if(tautest.and.myrank.eq.root) then
c      write(6,*) '<cpacoord> : psq = ',psq
c      write(6,*) '<cpacoord> : before tau2d'
c      write(6,*) '<cpacoord> : irel = ',irel
c      write(6,*) '<cpacoord> : wgeff = ',wgeff
c      write(6,*) '<cpacoord> : bulkgeo = ',bulkgeo
c      write(6,*) '<cpacoord> : lmax = ',lmax
c      write(6,*) '<cpacoord> : kmymax = ',nbogomax
c      write(6,*) '<cpacoord> : kmymaxp = ',dbogomaxp
c      write(6,*) '<cpacoord> : nintfc = ',nintfc
       write(6,*) '<cpacoord> : gg2d = '
         call outmat1(gg2d(1:nl*nl,1:nl*nl,0,1,0),nl*nl,nl*nl,nl*nl,
     >     tol,6)
c      write(6,*) '<cpacoord> : ttmp(1,1,1,ig) = ',ttmp(1,1,1,ig)
       write(6,*) '<cpacoord> : xx = '
         call outmat1(xx(1:kmymax,1:kmymax),kmymax,kmymax,kmymax,tol,6)
       write(6,*) '<cpacoord> : yy = '
         call outmat1(yy(1:kmymax,1:kmymax),kmymax,kmymax,kmymax,tol,6)
c      write(6,*) '<cpacoord> : mdim = ',mdim
       end if
c           ----------------------------------------------------
            call tau2d(irel,wgeff,bulkgeo,lmax,nbogomax,dbogomaxp,
     >        nintfc,ttmp(1,1,1,ig),gg2d,xx,yy,tau1(1,1,1,ig),mdim)
c           ----------------------------------------------------
       tautest=.false.
       if(tautest.and.myrank.eq.root) then
       write(6,*) '<cpacoord> : before tau2d'
       write(6,*) '<cpacoord> : irel = ',irel
       write(6,*) '<cpacoord> : wgeff = ',wgeff
       write(6,*) '<cpacoord> : bulkgeo = ',bulkgeo
       write(6,*) '<cpacoord> : lmax = ',lmax
       write(6,*) '<cpacoord> : kmymax = ',nbogomax
       write(6,*) '<cpacoord> : kmymaxp = ',dbogomaxp
       write(6,*) '<cpacoord> : nintfc = ',nintfc
       write(6,*) '<cpacoord> : ttmp(1,1,1,ig) = ',ttmp(1,1,1,ig)
       write(6,*) '<cpacoord> : gg2d = '
         call outmat1(gg2d(1:nl*nl,1:nl*nl,0,1,0),nl*nl,nl*nl,nl*nl,
     >     tol,6)
       write(6,*) '<cpacoord> : xx = '
         call outmat1(xx(1:kmymax,1:kmymax),kmymax,kmymax,kmymax,tol,6)
       write(6,*) '<cpacoord> : yy = '
         call outmat1(yy(1:kmymax,1:kmymax),kmymax,kmymax,kmymax,tol,6)
       write(6,*) '<cpacoord> : mdim = ',mdim
       end if
       tautest=.false.
       if(tautest) then
         write(6,*) ' <cpacoord> : tautest, tauii ig=',ig
         call outmat1(tau1(1,1,1,ig),nbogomax,nbogomax,dbogomaxp,
     >                tol,6)
       end if
          end if 
c           write(6,*) ' <cpacoord> : tau2d ie=',ie
       tautest=.false.
c
c          do li=1,nintfc
c             if(itest.gt.2) then
c                write(6,*) ' TAU_before BZ sum'
c                call replms(taudiag,tau1(1,1,li,ige),lmax)
c                do k=1,kmymax
c                   write(6,'(i3,2d20.10)') k,taudiag(k)
c                end do
c             end if
c          end do
c
          end if
c
c **************************************
c BZ sum for layer diagonal tau matrices
c **************************************
c
          do li=1,nintfc
c           -----------------------------------------------
            if(calcoff) then
              call repldim(tau2,tauij1(1,1,li,li,ige),nbogomax,
     >                     dbogomaxp,nbogomax)
            else
              call repl(tau2,tau1(1,1,li,ige),nbogomax,dbogomaxp)
c              call repl(tau2,tau1(1,1,li,ige),kmymax,kmymaxp)
            end if
c           -----------------------------------------------
            call tripmt(rmatp(1,1,ig1),tau2,rmat(1,1,ig1),
     >                  nbogomax,nbogomax,dbogomaxp)        
c           -----------------------------------------------
            call addmat(tau(1,1,li),tau2,nbogomax,dbogomaxp)
c           -----------------------------------------------
          end do
c
c ======================================================================
c ----------------- off diagonal elements of tau matrix ----------------
c
          if((calcoff).AND.(npair.gt.0)) then
c           write(6,*) '<cpacoord>: Calculation of off-diag elments' 
            tautest=.false.
            do li = 1,nintfc
             do lj = 1,nintfc
              if(tautest) then
               write(6,*) '<cpacoord> : tautest tauij1(:,:,li,lj,ige)
     > before rotation',
     >                   li,lj,ige
               call outmat1(tauij1(1:kmymax,1:kmymax,li,lj,ige),
     >                kmymax,kmymax,kmymax,tol,6)
              end if
              call repldim(help,tauij1(1,1,li,lj,ige),
     >                  nbogomax,dbogomaxp,nbogomax)
              call tripmt(rmatp(1,1,ig1),help,rmat(1,1,ig1),
     >                       nbogomax,nbogomax,dbogomaxp)
              call repldim(tauij(1,1,li,lj),help,
     >                  nbogomax,nbogomax,dbogomaxp)
              if(tautest) then
               write(6,*) '<cpacoord> : tautest tauij(:,:,li,lj,ige)
     > after rotation',
     >                   li,lj,ige
               call outmat1(tauij(1:kmymax,1:kmymax,li,lj),
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
c
              if(tautest) then
                write(6,*) '<cpacoord> : tautest ffacx, ffacxp'
                write(6,*) ffacx,ffacxp
              end if
c
              do k1 = 1,nbogomax
                do k2 = 1,nbogomax
                  tau_ij(k1,k2,ipair) = tau_ij(k1,k2,ipair) +
     >                  ffacx*tauij(k1,k2,li,lj)
                  tau_ji(k1,k2,ipair) = tau_ji(k1,k2,ipair) +
     >                  ffacxp*tauij(k1,k2,lj,li)
                end do
              end do
            end do
c           write(6,*) '<cpacoord>: Off-diag elments ready' 
          end if
          tautest=.false.
c ----------------------------------------------------------------------
c ======================================================================
        end do
c           write(6,*) ' <cpacoord> : offdiag ie=',ie
c *********************************
c end of loop over point operations
c *********************************
c
#ifdef MPIP
      end if
#endif
c MPI if END
      end do
      if(calcoff) then
        deallocate(tauij1,stat=AllocateStatus)
        call alloccheck( AllocateStatus,
     >        'tauij1 dealloc in cpacoord                        ')
        deallocate(tauij,stat=AllocateStatus)
        call alloccheck( AllocateStatus,
     >        'tauij dealloc in cpacoord                         ')
      end if
c *************************
c end of loop over k points
c *************************
c
#ifdef MPIP
       taujoin = (0.d0,0.d0)
       call mpi_allreduce(tau,taujoin,dbogomaxp*dbogomaxp*mintfc,
     >     mpi_double_complex,mpi_sum,mpi_comm_world,ierr)
c
       tau = taujoin
c
       tautest=.false.
       if(calcoff.and.npair.gt.0) then
         if(tautest) then
           write(6,*) ' <cpacoord> : tau_ij/ji(n,n,npair) before 
     >mpi_allreduce:'
           write(6,*) tau_ij(nbogomax,nbogomax,npair)
           write(6,*) tau_ji(nbogomax,nbogomax,npair)
         end if
         ALLOCATE ( tau_ijsum(nbogomax,nbogomax,npair),
     >                   STAT = AllocateStatus)
         call alloccheck(AllocateStatus,
     >        'tau_ijsum in cpacoord                             ')
         tau_ijsum = (0.d0,0.d0)
         call mpi_allreduce(tau_ij,tau_ijsum,nbogomax*nbogomax*npair,
     >       mpi_double_complex,mpi_sum,mpi_comm_world,ierr)
c
         tau_ij = tau_ijsum
c
         tau_ijsum = (0.d0,0.d0)
         call mpi_allreduce(tau_ji,tau_ijsum,nbogomax*nbogomax*npair,
     >       mpi_double_complex,mpi_sum,mpi_comm_world,ierr)
c
         tau_ji = tau_ijsum
         deallocate(tau_ijsum,stat=AllocateStatus)
         call alloccheck(AllocateStatus,
     >        'tau_ijsum dealloc in cpacoord                     ')
         if(tautest) then
           write(6,*) ' <cpacoord> : tau_ij/ji(n,n,npair) after 
     >mpi_allreduce:'
           write(6,*) tau_ij(nbogomax,nbogomax,npair)
           write(6,*) tau_ji(nbogomax,nbogomax,npair)
         end if
       end if
c
c           write(6,*) ' <cpacoord> : allreduce ie=',ie
#endif
c **********
c CPA solver
c **********
      if (itcpam.gt.0) then 
c     -----------------------------------------------------------------
      call cpacor(conc,tminv,tminva,tminvb,tau,taua,taub,nbogomax,
     >            nintfc,
     >            itcpa,itcpam,cpatol,cpacon,concpa,cpatest,cpaerr)
c     -----------------------------------------------------------------
c
      if(concpa.and.cpacon) goto 100
      if(itest.ge.1.and.myrank.eq.root) 
     >   write(6,'('' CPA iterations:'',i3,2x,
     >  '' ERROR:'',d15.6)') itcpa-1,cpaerr
      else 
      if(itest.ge.2) write(6,*) 'cpacor is skipped'
      end if 
c
c                        *************************
c                        * ending CPA iterations *
c                        *************************
c
c ======================================================================
c - the off-diagonal elements of tau
c
c
      if(calcoff) then
      tautest=.false.
      ipair=1
      if(tautest) then
         write(6,*) ' <cpacoord> : tautest, tauij before phystau2,
     > electron-electron, ipair=',ipair
         call outmat1(tau_ij(1:kmymax,1:kmymax,ipair),
     >                kmymax,kmymax,kmymax,
     >                tol,6)
         write(6,*) ' <cpacoord> : tautest, tauij before phystau2,
     > electron-hole, ipair=',ipair
         call outmat1(tau_ij(1:kmymax,kmymax+1:nbogomax,ipair),
     >                kmymax,kmymax,kmymax,
     >                tol,6)
         write(6,*) ' <cpacoord> : tautest, tauij before phystau2,
     > hole-electron, ipair=',ipair
         call outmat1(tau_ij(kmymax+1:nbogomax,1:kmymax,ipair),
     >                kmymax,kmymax,kmymax,
     >                tol,6)
         write(6,*) ' <cpacoord> : tautest, tauij before phystau2,
     > hole-hole, ipair=',ipair
         call outmat1(tau_ij(kmymax+1:nbogomax,kmymax+1:nbogomax,ipair),
     >                kmymax,kmymax,kmymax,
     >                tol,6)
      end if
       if(npair.gt.0) then
        do ipair = 1,npair
          li = kpair(3,ipair)
          lj = kpair(4,ipair)
          idiag = 0
c         if(li.eq.lj) idiag=1
c         call phystau2(tau_ij,tminva(1,1,li),tminva(1,1,lj),
c    >                  alphaintkkr(0,li),alphaintkkr(0,lj),lmax,idiag)
          call phystau2(tau_ij(1,1,ipair),nbogomax,tminv(1,1,li),
     >                  tminv(1,1,lj),
     >                  alphaintkkr(0,li),alphaintkkrh(0,li),
     >                  alphaintkkr(0,lj),alphaintkkrh(0,lj),lmax,idiag)
          call phystau2(tau_ji(1,1,ipair),nbogomax,tminv(1,1,lj),
     >                  tminv(1,1,li),
     >                  alphaintkkr(0,lj),alphaintkkrh(0,lj),
     >                  alphaintkkr(0,li),alphaintkkrh(0,li),lmax,idiag)
        end do
       end if
       ALLOCATE ( tmusch(nbogomax,nbogomax,nintfc),
     >                 STAT = AllocateStatus)
       call alloccheck( AllocateStatus,
     >        'tmusch in cpacoord                                ' )
c unscreening the tmatices
c
c       do li=1,nintfc
c           write(6,*) ' <cpacoord> : unsct ie=',ie
c         call unsct(tminv(1,1,li),help,nbogomax,
c    >               alphaintkkr(0,li),alphaintkkrh(0,li),li)
c           write(6,*) ' <cpacoord> : repldim ie=',ie
c         call repldim(tmusch(1,1,li),help,nbogomax,nbogomax,dbogomaxp)
c       end do
c
c - the diagonal elements of tau
        do li=1,nintfc
c           write(6,*) ' <cpacoord> : repldim ie=',ie
          call repldim(tau_ii(1,1,li),tau(1,1,li),nbogomax,nbogomax,
     >                 dbogomaxp)
c         --------------------------------------------------------
c           write(6,*) ' <cpacoord> : phystau ie=',ie
          call phystau(tau_ii(1,1,li),nbogomax,tminv(1,1,li),
     >                 tminv(1,1,li),alphaintkkr(0,li),
     >                 alphaintkkrh(0,li),lmax,1)
c         --------------------------------------------------------
        end do
c
c       write(6,*) '<cpacoord>:',taupath
         if (hostprint.and.myrank.eq.root) then
c          call thout(ie,ce,nintfc,nbogomax,tmusch,taupath)
          tautest=.false.
          if(tautest) then
            write(6,*) ' <cpacoord> : ie=',ie
            write(6,*) 'ce=',ce
            write(6,*) 'npair=',npair
            write(6,*) 'nintfc=',nintfc
            write(6,*) 'nk=',nk
            write(6,*) 'nbogomax=',nbogomax
            write(6,*) 'tau_ij(1,1,1)=',tau_ij(1,1,1)
            write(6,*) 'tau_ji(1,1,1)=',tau_ji(1,1,1)
            write(6,*) 'tau_ii(1,1,1)=',tau_ii(1,1,1)
            write(6,*) 'taupath=',taupath
          end if
          call tauhout(ie,ce,npair,npair1,kpair,nintfc,nk,nbogomax,
     >               tau_ij,tau_ji,tau_ii,taupath)
         end if
c
      tautest=.false.
      ipair=1
      if(tautest) then
         write(6,*) ' <cpacoord> : tautest, tauij after phystau2,
     > electron-electron, ipair=',ipair
         call outmat1(tau_ij(1:kmymax,1:kmymax,ipair),
     >                kmymax,kmymax,kmymax,
     >                tol,6)
         write(6,*) ' <cpacoord> : tautest, tauij after phystau2,
     > electron-hole, ipair=',ipair
         call outmat1(tau_ij(1:kmymax,kmymax+1:nbogomax,ipair),
     >                kmymax,kmymax,kmymax,
     >                tol,6)
         write(6,*) ' <cpacoord> : tautest, tauij after phystau2,
     > hole-electron, ipair=',ipair
         call outmat1(tau_ij(kmymax+1:nbogomax,1:kmymax,ipair),
     >                kmymax,kmymax,kmymax,
     >                tol,6)
         write(6,*) ' <cpacoord> : tautest, tauij after phystau2,
     > hole-hole, ipair=',ipair
         call outmat1(tau_ij(kmymax+1:nbogomax,kmymax+1:nbogomax,ipair),
     >                kmymax,kmymax,kmymax,
     >                tol,6)
      end if
       if(tautest) then
        do li=1,nintfc
         write(6,*) ' <cpacoord> : tautest, tm li=',li
         call outmat1(tmusch(1,1,li),nbogomax,nbogomax,dbogomaxp,
     >                tol,6)
         write(6,*) ' <cpacoord> : tautest, tauii li=',li
         call outmat1(tau_ii(1,1,li),nbogomax,nbogomax,dbogomaxp,
     >                tol,6)
        end do
        do ipair = 1,npair
         write(6,*) ' <cpacoord> : tautest, tauji ipair=',ipair
         call outmat1(tau_ji(1,1,ipair),nbogomax,nbogomax,dbogomaxp,
     >                tol,6)
         write(6,*) ' <cpacoord> : tautest, tauij ipair=',ipair
         call outmat1(tau_ij(1,1,ipair),nbogomax,nbogomax,dbogomaxp,
     >                tol,6)
        end do
       end if
       deallocate(tmusch,stat=AllocateStatus)
       call alloccheck( AllocateStatus,
     >        'tmusch dealloc in cpacoord                        ' )
      end if
c ======================================================================
c
c *************************************************************
c * loop over layers to transform layer diagonal tau matrices *
c * into physical representation and rotate to the local      *
c * frame of reference                                        *
c *************************************************************
c 
      do li=1,nintfc
         cpalay=(1.d0-conc(li)).gt.tiny
c
      if (.not.cpalay) then
        if (itest.gt.2) then
         write(6,*) 'taua is replaced with tau, beacause we have no cpa'
        end if
         taua(:,:,li)=tau(:,:,li)
      end if 
c
c        --------------------------------------------------------
         call phystau(taua(1,1,li),dbogomaxp,tminva(1,1,li),
     >                tminva(1,1,li),
     >                alphaintkkr(0,li),alphaintkkrh(0,li),lmax,1)    
c         write(6,*) 'tau_a after phystau in cpacord for li = ', li
c         call outmat1(taua(1,1,li),2*kmymax,2*kmymax,dbogomaxp,outtol,6)
c         write(6,*) 'tau in cpacord for li = ', li
c         call outmat1(tau(1,1,li),2*kmymax,2*kmymax,dbogomaxp,outtol,6)
         call phystau(tau(1,1,li),dbogomaxp,tminv(1,1,li),tminv(1,1,li),
     >                alphaintkkr(0,li),alphaintkkrh(0,li),lmax,1)    
c         write(6,*) 'tau after phystau in cpacord for li = ', li
c         call outmat1(tau(1,1,li),2*kmymax,2*kmymax,dbogomaxp,outtol,6)

c         taua(1:2*kmymax,1:2*kmymax,li)=
c     >        matmul(taua(1:2*kmymax,1:2*kmymax,li),ctrafo)
c         taua(1:2*kmymax,1:2*kmymax,li)=
c     >        matmul(ctrafo2,taua(1:2*kmymax,1:2*kmymax,li))
        
c         if(itest.ge.2) then
c         write(6,'(/'' tau -A '',i2)') li
c         call outmat1(taua(1,1,li),2*kmymax,2*kmymax,dbogomaxp,tol,6)
c         end if
c        --------------------------------------------------------
c
c rotate tau-matrix to local frame of reference if necessary
c
         call repl(gtaua(1,1,li),taua(1,1,li),nbogomax,dbogomaxp)
c        ------------------------------------------------------
c         if (localmode) then    ! localmode defined in param.h!
c            call tripmt(dmatpa(1,1,li),taua(1,1,li),dmata(1,1,li),
c     >                  2*kmymax,2*kmymax,dbogomaxp)
c         end if
c        ------------------------------------------------------
c
         tautest=.false.
         if(itest.ge.2.or.tautest) then
            write(6,'(/'' tau -A '',i2)') li
c           call replms(taudiag,taua(1,1,li),lmax)
c           do k=1,kmymax
c             write(6,'(i3,2d20.10)') k,taudiag(k)
c           end do
            call outmat1(taua(1,1,li),nbogomax,nbogomax,dbogomaxp,tol,6)
            write(6,'(/'' tau    '',i2)') li
            call outmat1(tau(1,1,li),nbogomax,nbogomax,dbogomaxp,tol,6)
         end if
c
         if(cpalay) then
c+------------+
c+ BIG CPA IF +
c+------------+
c
c        --------------------------------------------------------
         call phystau(taub(1,1,li),dbogomaxp,tminvb(1,1,li),
     >                tminvb(1,1,li),
     >                alphaintkkr(0,li),alphaintkkrh(0,li),lmax,1)
c        --------------------------------------------------------
c         write(6,*) 'tau_b after phystau in cpacord for li = ', li
c         call outmat1(taub(1,1,li),2*kmymax,2*kmymax,dbogomaxp,outtol,6)
c
c rotate tau-matrix to local frame of reference if necessary
c
         call repl(gtaub(1,1,li),taub(1,1,li),nbogomax,dbogomaxp)
c        ---------------------------------------------------
c         if (localmode) then    ! localmode defined in param.h!
c            call tripmt(dmatpb(1,1,li),taub(1,1,li),dmatb(1,1,li),
c     >                  kmymax,kmymax,kmymaxp)
c         end if
c        ---------------------------------------------------
c
         if(itest.gt.2) then
            write(6,'(/'' tau -B '',i2)') li
            call outmat1(taub(1,1,li),nbogomax,nbogomax,dbogomaxp,tol,6)
         end if
c
         else 
            taub(:,:,li) = taua(:,:,li)
         end if
c+----------------+
c+ END BIG CPA IF +
c+----------------+
      end do
c *** end loop over layers ***
c
      return
      end

