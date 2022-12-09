c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      program main
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
#ifdef MPIP
      include 'mpif.h'
#endif
c
      logical bulk,dos,wrel,lms,bulkgeo,linbw,calcder
      logical orbpol,opotl,opotr,opot
      logical cpatest,cpamatin,cpamatinl,cpamatinr,cpain
      logical singlesite,reg
      logical bsfcalc
      logical skiploq
c
      character*10 idpota(mintfc),idpotla(minprc),idpotra(minprc)
      character*10 idpotb(mintfc),idpotlb(minprc),idpotrb(minprc)
      character*30 printout,potout,dosout,pdosout,dosoutimp,tmatout
      character*34 pdosoutimp,tripletout
      character*30 momout
      character*30 leftpot,rightpot,laypot,laycore
      character*30 leftdelta,rightdelta,laydelta
      character*30 leftmat,rightmat,laymat
      character*30 leftmom,rightmom
      character*1 rightm
c
      real*8 lza(2,mintfc),lzb(2,mintfc)
      real*8 tstart,tfinish,tarrive
      real*8 entotifc
      real*8 deltaprint
c
      dimension park(2),kset(me)
      dimension xk(mkpar,2),wk(mkpar)
      dimension nepanel(5),ebottom(5),etop(5),eps(5)
      dimension iesublatt(mintfc)
c
      complex*16 deltala(nrad,minprc),deltara(nrad,minprc)
      complex*16 deltaa(nrad,mintfc)
      complex*16 deltaold(nrad)
      complex*16 deltalb(nrad,minprc),deltarb(nrad,minprc)
      complex*16 deltab(nrad,mintfc)
c
c      complex*16 deltal(nrad,minprc),deltar(nrad,minprc)
c      complex*16 delta(nrad,mintfc)
c      complex*16 delta0(nrad,mintfc)
c      dimension redelta(nrad,mintfc)
c      dimension imdelta(nrad,mintfc)
      dimension vra(nrad,mintfc),bra(nrad,mintfc)
      dimension bopra(nrad,2,mintfc)
      dimension vrb(nrad,mintfc),brb(nrad,mintfc)
      dimension boprb(nrad,2,mintfc)
      dimension rs(mintfc),dx(mintfc),ns(mintfc),rmt(mintfc)
      dimension vrla(nrad,minprc),brla(nrad,minprc)
      dimension vrlb(nrad,minprc),brlb(nrad,minprc)
      dimension boprla(nrad,2,minprc),boprlb(nrad,2,minprc)
      dimension rsl(minprc),dxl(minprc),nsl(minprc),rmtl(minprc)
      dimension vrra(nrad,minprc),brra(nrad,minprc)
      dimension vrrb(nrad,minprc),brrb(nrad,minprc)
      dimension boprra(nrad,2,minprc),boprrb(nrad,2,minprc)
      dimension rsr(minprc),dxr(minprc),nsr(minprc),rmtr(minprc)
c
      dimension za(mintfc),zb(mintfc),zla(minprc),zlb(minprc),
     &          zra(minprc),zrb(minprc)
      dimension conc(mintfc),sxa(mintfc),sxb(mintfc)
      dimension concl(minprc),sxl(minprc)
      dimension concr(minprc),sxr(minprc)
      dimension qca(mintfc),qcb(mintfc)
      dimension qva(mintfc),qvb(mintfc),qvha(mintfc),qvhb(mintfc)
      dimension qvpa(kmymaxp,mintfc),qvpb(kmymaxp,mintfc)
c
      dimension rbl(3,minprc),rbr(3,minprc)
      dimension rba(3,mintfc),rbb(3,mintfc)
      dimension vecna(3,mintfc),phia(mintfc)
      dimension vecnb(3,mintfc),phib(mintfc)
c
      complex*16 rhoveha(nrad,mintfc)
      complex*16 rhovhea(nrad,mintfc)
      complex*16 rhovehb(nrad,mintfc)
      complex*16 rhovheb(nrad,mintfc)
      complex*16 qveha(mintfc),qvhea(mintfc)
      complex*16 qvehb(mintfc),qvheb(mintfc)
      complex*16 qvteha(mintfc),qvthea(mintfc)
      complex*16 qvtehb(mintfc),qvtheb(mintfc)
c
      complex*16 dmatl(kmymaxp,kmymaxp,minprc)
      complex*16 dmatlp(kmymaxp,kmymaxp,minprc)
      complex*16 dmatr(kmymaxp,kmymaxp,minprc)
      complex*16 dmatrp(kmymaxp,kmymaxp,minprc)
      complex*16 dmata(kmymaxp,kmymaxp,mintfc)
      complex*16 dmatpa(kmymaxp,kmymaxp,mintfc)
      complex*16 dmatb(kmymaxp,kmymaxp,mintfc)
      complex*16 dmatpb(kmymaxp,kmymaxp,mintfc)
      complex*16 rmata(lmsup,lmsup,mintfc),rmatpa(lmsup,lmsup,mintfc)
      complex*16 rmatb(lmsup,lmsup,mintfc),rmatpb(lmsup,lmsup,mintfc)
c
      complex*16 ddpha(kmymaxp,kmymaxp,mintfc)
      complex*16 ddtha(kmymaxp,kmymaxp,mintfc)
      complex*16 d2dpha(kmymaxp,kmymaxp,mintfc)
      complex*16 d2dtha(kmymaxp,kmymaxp,mintfc)
      complex*16 d2dthpha(kmymaxp,kmymaxp,mintfc)
      complex*16 ddphpa(kmymaxp,kmymaxp,mintfc)
      complex*16 ddthpa(kmymaxp,kmymaxp,mintfc)
      complex*16 d2dphpa(kmymaxp,kmymaxp,mintfc)
      complex*16 d2dthpa(kmymaxp,kmymaxp,mintfc)
      complex*16 d2dthphpa(kmymaxp,kmymaxp,mintfc)
c
      complex*16 ddphb(kmymaxp,kmymaxp,mintfc)
      complex*16 ddthb(kmymaxp,kmymaxp,mintfc)
      complex*16 d2dphb(kmymaxp,kmymaxp,mintfc)
      complex*16 d2dthb(kmymaxp,kmymaxp,mintfc)
      complex*16 d2dthphb(kmymaxp,kmymaxp,mintfc)
      complex*16 ddphpb(kmymaxp,kmymaxp,mintfc)
      complex*16 ddthpb(kmymaxp,kmymaxp,mintfc)
      complex*16 d2dphpb(kmymaxp,kmymaxp,mintfc)
      complex*16 d2dthpb(kmymaxp,kmymaxp,mintfc)
      complex*16 d2dthphpb(kmymaxp,kmymaxp,mintfc)
c
      dimension dosa(kmymaxp,mintfc,me),dosb(kmymaxp,mintfc,me)
      dimension dosmaga(kmymaxp,mintfc,me),dosmagb(kmymaxp,mintfc,me)
      dimension enbdiffa(mintfc,me),enbdiffb(mintfc,me)
      dimension qvdiffa(mintfc,me),qvdiffb(mintfc,me)
c
      dimension spin_magvpa(kmymaxp,mintfc,3)
      dimension spin_magvpb(kmymaxp,mintfc,3)
      dimension spin_magva(mintfc,3),spin_magvb(mintfc,3)
      dimension orb_magvpa(kmymaxp,mintfc,3)
      dimension orb_magvpb(kmymaxp,mintfc,3)
      dimension orb_magva(mintfc,3),orb_magvb(mintfc,3)
      dimension enba(mintfc),enbb(mintfc)
      dimension ddphenba(mintfc),ddphenbb(mintfc)
      dimension ddthenba(mintfc),ddthenbb(mintfc)
      dimension d2dphenba(mintfc),d2dphenbb(mintfc)
      dimension d2dthenba(mintfc),d2dthenbb(mintfc)
      dimension d2dthphenba(mintfc),d2dthphenbb(mintfc)
      dimension enca(mintfc),encb(mintfc)
      dimension enela(mintfc),enelb(mintfc)
      dimension enxca(mintfc),enxcb(mintfc)
      dimension enpota(mintfc),enpotb(mintfc)
      dimension enmaga(mintfc),enmagb(mintfc)
      dimension enorba(mintfc),enorbb(mintfc)
      dimension entota(mintfc),entotb(mintfc),entot(mintfc)
      dimension rhoca(nrad,mintfc),rhocb(nrad,mintfc)
      dimension rhova(nrad,mintfc),rhovb(nrad,mintfc)
      dimension rhospa(nrad,2,mintfc),rhospb(nrad,2,mintfc)
      dimension rhodspa(nrad,2,mintfc),rhodspb(nrad,2,mintfc)
      dimension rhomaga(nrad,mintfc),rhomagb(nrad,mintfc)
c
      dimension spinmoma(mintfc),orbmoma(mintfc)
      dimension th0a(mintfc),th1a(mintfc)
      dimension ph0a(mintfc),ph1a(mintfc)
      dimension spinmomb(mintfc),orbmomb(mintfc)
      dimension th0b(mintfc),th1b(mintfc)
      dimension ph0b(mintfc),ph1b(mintfc)
c
      dimension dosha(kmymaxp,mintfc,me),doshb(kmymaxp,mintfc,me)
      dimension doseha(kmymaxp,mintfc,me),doshea(kmymaxp,mintfc,me)
      dimension dosehb(kmymaxp,mintfc,me),dosheb(kmymaxp,mintfc,me)
      dimension dosteha(kmymaxp,mintfc,me),dosthea(kmymaxp,mintfc,me)
      dimension dostehb(kmymaxp,mintfc,me),dostheb(kmymaxp,mintfc,me)
      dimension dosmagha(kmymaxp,mintfc,me)
      dimension dosmaghb(kmymaxp,mintfc,me)
c
      dimension rlambdala(minprc),rlambdara(minprc)
      dimension rlambdaa(mintfc)
      dimension rlambdalb(minprc),rlambdarb(minprc)
      dimension rlambdab(mintfc)
c
      complex*16 cear(me),we(me),ce,ze
      complex*16 tminvl(dbogomaxp,dbogomaxp,minprc)
      complex*16 tminvr(dbogomaxp,dbogomaxp,minprc)
      complex*16 tminv(dbogomaxp,dbogomaxp,mintfc)
      complex*16 tminva(dbogomaxp,dbogomaxp,mintfc)
      complex*16 tminvb(dbogomaxp,dbogomaxp,mintfc)
      complex*16 ptminva(dbogomaxp,dbogomaxp,mintfc)
      complex*16 ptminvb(dbogomaxp,dbogomaxp,mintfc)
      complex*16 tau(dbogomaxp,dbogomaxp,mintfc)
      complex*16 taua(dbogomaxp,dbogomaxp,mintfc)
      complex*16 taub(dbogomaxp,dbogomaxp,mintfc)
      complex*16 tau_kint(dbogomaxp,dbogomaxp,mintfc)
      complex*16 taua_kint(dbogomaxp,dbogomaxp,mintfc)
      complex*16 taub_kint(dbogomaxp,dbogomaxp,mintfc)
      complex*16 gtaua(dbogomaxp,dbogomaxp,mintfc)
      complex*16 gtaub(dbogomaxp,dbogomaxp,mintfc)
      complex*16 tempmat(dbogomaxp,dbogomaxp)
      complex*16 qmoma(lmsup,mintfc),qmomb(lmsup,mintfc)
      complex*16 qmomla(lmsup,minprc),qmomlb(lmsup,minprc)
      complex*16 qmomra(lmsup,minprc),qmomrb(lmsup,minprc)
      complex*16 xlms(lmmaxp,lmmaxp,2,2,2,2)
      complex*16 detl
      complex*16 tautmp(dbogomaxp,dbogomaxp)
c cluster
      integer*4  AllocateStatus,ie
      complex*16, allocatable :: tminvh(:,:,:,:)
      complex*16, allocatable :: tau_ii1(:,:,:)
      complex*16, allocatable :: tau_ij1(:,:,:)
      complex*16, allocatable :: tau_ji1(:,:,:)
      complex*16, allocatable :: tau_ii(:,:,:,:)
      complex*16, allocatable :: tau_ij(:,:,:,:)
      complex*16, allocatable :: tau_ji(:,:,:,:)
      complex*16, allocatable :: taucl(:,:,:)
c     complex*16, allocatable :: gtaucl(:,:,:)
      complex*16, allocatable :: tminvcl(:,:,:)
      complex*16, allocatable :: ptminvcl(:,:,:)
c
      common/test/itest
      common/broypot/cmix,wbr,nbrmax(100)
      common/lay2d/cvec(mtotal,3),nextra,nbulkl,nbulkr,
     &             nprc,ninprc(0:mprc+1) 
c
      real*8 tol
      data tol/1.0d-8/ 
      data tiny/1.0d-6/
      data rnul/0.d0/
      data sqrtm1/(0.d0,1.d0)/
c     data tolerr,tolef,tolen,tolv0/1.0d-10,1.0d-10,1.0d-9,1.0d-08/
c MPI
#ifdef MPIP
      integer root, myrank, nprocs, ierror
      integer sendto,senddata,tag,sendcount,from
      integer status(MPI_STATUS_SIZE) 
      common/mpiparam/root,myrank,nprocs
#endif
c Version control
      logical clucalc,hostprint
      logical tautest
c
c -----------immpurities----------
      character*30 imppot
      character*300 taupath
c
      real*8 b0
c
      integer ib0
      integer impb0
      integer impb0f
      integer, allocatable ::  isigb(:)
c  pothandle imp
      character*10, allocatable :: idpotimp(:)
      real*8, allocatable ::  vrimp(:,:)
      real*8, allocatable :: brimp(:,:)
      real*8, allocatable :: boprimp(:,:,:)
      real*8, allocatable :: dximp(:)
      real*8, allocatable :: rsimp(:)
      real*8, allocatable :: zimp(:)
      complex*16, allocatable :: deltaimp(:,:)
c     dimension, allocatable :: rlambdaimp(:)
      real*8, allocatable :: rlambdaimp(:)
      real*8, allocatable :: singratimp(:)
      real*8, allocatable :: uratimp(:)
      real*8, allocatable :: dratimp(:)
      integer, allocatable ::  nsimp(:)
c  readclu
c     integer iscreencl
      integer npair,npair0,npair1
      integer madmax
      integer, allocatable ::  nposimp(:,:)
      integer, allocatable ::  kpairind(:,:)
      integer, allocatable ::  kpair(:,:)
      real*8, allocatable ::  rbcl(:,:)
      real*8, allocatable ::  bthcl(:)
      real*8, allocatable ::  bphcl(:)
      real*8, allocatable ::  sxcl(:)
c locquant2
      real*8, allocatable ::  dosimp(:,:,:)
      real*8, allocatable ::  doshimp(:,:,:)
      real*8, allocatable ::  dosehimp(:,:,:)
      real*8, allocatable ::  dosheimp(:,:,:)
      real*8, allocatable ::  dostehimp(:,:,:)
      real*8, allocatable ::  dostheimp(:,:,:)
      real*8, allocatable ::  qvimp(:),qvpimp(:,:)
      real*8, allocatable ::  qvhimp(:),qvhpimp(:,:)
      real*8, allocatable ::  qvehimp(:),qvheimp(:)
      real*8, allocatable ::  qvtehimp(:),qvtheimp(:)
      complex*16, allocatable ::  qmomimpa(:,:)
      complex*16, allocatable ::  qmomimpha(:,:)
      complex*16, allocatable ::  qmomimpeha(:)
      complex*16, allocatable ::  qmomimphea(:)
      complex*16, allocatable ::  qmomimpteha(:)
      complex*16, allocatable ::  qmomimpthea(:)
      real*8, allocatable ::  vmadih(:)
      real*8, allocatable ::  vmadich(:)
      real*8, allocatable ::  vmadid(:)
      real*8, allocatable ::  vmadiq(:)
      real*8, allocatable ::  enbimp(:),enbhimp(:)
      real*8, allocatable ::  enbdiffimp(:,:),enbhdiffimp(:,:)
      real*8, allocatable ::  rhovimp(:,:)
      real*8, allocatable ::  rhospimp(:,:,:)
      real*8, allocatable ::  rhodspimp(:,:,:)
      real*8, allocatable ::  rhomagimp(:,:)
      real*8, allocatable ::  dosmagimp(:,:,:),dosmaghimp(:,:,:)
      real*8, allocatable ::  spin_magvpimp(:,:,:),spin_magvimp(:,:)
      real*8, allocatable ::  spin_magvhpimp(:,:,:),spin_magvhimp(:,:)
      real*8, allocatable ::  orb_magvpimp(:,:,:),orb_magvimp(:,:)
      real*8, allocatable ::  orb_magvhpimp(:,:,:),orb_magvhimp(:,:)
      real*8, allocatable ::  tripletmat(:,:,:,:,:)
c
      integer rnet,cnet !reduced net, current net
c
c     ***************************************************************
c                       INITIALIZE LAYER CALCULATION
c     ***************************************************************
c
      call cpu_time(tstart)
c MPI begin
#ifdef MPIP
      root = 0
      call mpi_init(ierror)
      call mpi_comm_rank(mpi_comm_world,myrank,ierror)
      call mpi_comm_size(mpi_comm_world,nprocs,ierror)
      write(*,*)  'Process ', myrank, ' of ', nprocs, ' is alive'
#else 
      root=0
      myrank=0
      nprocs=1
#endif
c   MPI END
c
      open(unit=5,file='input_rsp.in',status='old')
c     ------------------------------------------------------------------
      call readini(
     > imesh,ne,npanel,nepanel,ne1,ne2,ne3,ebottom,etop,eps,lmax,eta,
     > sigma,park,kset,ksetmax,itcpam,cpatol,cpatest,cpain,
     > printout,potout,dosout,pdosout,dosoutimp,pdosoutimp,tripletout,
     > tmatout,momout,
     > leftpot,rightpot,laypot,laycore,
     > leftmat,rightmat,laymat,leftmom,rightmom,
     > wrel,lms,bulk,dos,rightm,newvac,ivacpot,v0,v00,vrsh,
     > vscreen,iscreen,itscfmax,tolerr,tolef,tolen,intbz,kunning,
     > orbpol,opotl,opotr,opot,lliter,E_Fermi,
c     > leftdelta,rightdelta,laydelta,singlesite,reg,c_light,
     > singlesite,reg,c_light,
     > clucalc,hostprint,nimp,npair,npair0,npair1)
c     ------------------------------------------------------------------
      close(5)
c
      if (myrank.eq.root) then
      write(6,'(/" Running on",i3," parallel processes."/)') nprocs
      end if
c
      open(unit=5,file='input_geo.in',status='old')
c     ------------------------------------------------------------------
      call readgeom(bulk,bulkgeo,rightm,intbz,nintfc,arot,
     > rsl,concl,rbl,sxl,
     > rs,conc,rba,rbb,sxa,sxb,
     > rsr,concr,rbr,sxr,
     > lfix1,lfix2,linbw,igraph,iesublatt)    
c     ------------------------------------------------------------------
      close(5)
c     ------------------------------------------------------------------
c     Allocate cluster related variables and read cluster geo 
      ALLOCATE ( kpair(4,npair0),STAT = AllocateStatus)
      call alloccheck( AllocateStatus,
     >     'kpair in skkr                                     ' )
      if(clucalc.or.hostprint) then
        ALLOCATE ( nposimp(3,nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'nposimp in skkr                                   ' )
        ALLOCATE ( kpairind(nimp,nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'kpairimp in skkr                                  ' )
        ALLOCATE ( rbcl(3,nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'rbcl in skkr                                      ' )
        ALLOCATE ( bthcl(nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'bthcl in skkr                                     ' )
        ALLOCATE ( bphcl(nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'bphcl in skkr                                     ' )
        ALLOCATE ( sxcl(nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'sxcl in skkr                                      ' )
        ALLOCATE ( isigb(nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'isigb in skkr                                     ' )
c      
        call readclu(nintfc,nimp,nposimp,npair,npair0,npair1,kpair,
     >             kpairind,
c    >             iscreencl,
     >             bthcl,bphcl,rbcl,sxcl,taupath,
     >             imppot,madmax,lfix1,lfix2,
     >             ib0,b0,impb0,impb0f,isigb)
      end if
c
c
      write(6,'(/2x,''routine MAIN>''/)')
c
      nl=lmax+1
      lmmax=nl*nl
      kmymax=2*lmmax
      nbogomax=2*kmymax
      ntotal=nintfc+ninprc(0)*(nextra+1)+ninprc(nprc+1)*(nextra+1)
c
      call flush(6)
      write(6,*) ' Dynamical arrays sizes:'
      write(6,*) '   nintfc=',nintfc
      write(6,*) '   nimp=',nimp
      write(6,*) '   npair=',npair
      write(6,*) '   npair0=',npair0
      write(6,*) '   npair1=',npair1
      call flush(6)
c -allocate cluster t and tau matrices
c     if(clucalc.neqv.hostprint) then
        ALLOCATE ( tau_ii1(nbogomax,nbogomax,nintfc), ! only for one energy point
     >                  STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'tau_ii in skkr                                    ' )
        ALLOCATE ( tau_ij1(nbogomax,nbogomax,npair1), ! only for one energy point
     >                  STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'tau_ij in skkr                                    ' )
        ALLOCATE ( tau_ji1(nbogomax,nbogomax,npair1), ! only for one energy point
     >                  STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'tau_ji in skkr                                    ' )
c     else
c       ALLOCATE ( tau_ii1(1,1,1), ! unused, but need to be passed
c    >                  STAT = AllocateStatus)
c       call alloccheck( AllocateStatus,'tau_ii in skkr' )
c       ALLOCATE ( tau_ij1(1,1,1), ! unused, but need to be passed
c    >                  STAT = AllocateStatus)
c       call alloccheck( AllocateStatus,'tau_ij in skkr' )
c       ALLOCATE ( tau_ji1(1,1,1), ! unused, but need to be passed
c    >                  STAT = AllocateStatus)
c       call alloccheck( AllocateStatus,'tau_ji in skkr' )
c     end if
      if(clucalc) then
        rnet=ne/nprocs
        if(myrank.lt.mod(ne,nprocs)) then
          rnet=rnet+1
        endif
        if(myrank.eq.root) then
          write(6,*) 'rnet',rnet
        endif
        ALLOCATE ( tminvh(nbogomax,nbogomax,nintfc,rnet),
     >                  STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'tminvh in skkr                                    ' )
        ALLOCATE ( taucl(nbogomax,nbogomax,nimp),
     >                  STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'taucl in skkr                                     ' )
c       ALLOCATE ( gtaucl(nbogomax,nbogomax,nimp), ! is it used?
c    >                  STAT = AllocateStatus)
c       call alloccheck( AllocateStatus,'gtaucl in skkr' )
        ALLOCATE ( tau_ii(nbogomax,nbogomax,nintfc,rnet),
     >                    STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'tau_ii in skkr                                    ' )
        ALLOCATE ( tau_ij(nbogomax,nbogomax,npair1,rnet),
     >                    STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'tau_ij in skkr                                    ' )
        ALLOCATE ( tau_ji(nbogomax,nbogomax,npair1,rnet),
     >                    STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'tau_ji in skkr                                    ' )
      end if ! end of clucalc if 
c
      if(clucalc.and.hostprint) goto 410 
c -read Left,Right and Layer potentials
c
c     ------------------------------------------------------------------
! layer potential needed
        call pothandle(
     > bulk,linbw,rightm,lmax,nintfc,ninprc(0),ninprc(nprc+1),
     > vrsbulk,v00,ivacpot,opotl,opotr,opot,
     > leftpot,leftmom,concl,qmomla,qmomlb,dxl,nsl,rsl,
     > idpotla,vrla,brla,boprla,zla,idpotlb,vrlb,brlb,boprlb,zlb,
     > rightpot,rightmom,concr,qmomra,qmomrb,dxr,nsr,rsr,
     > idpotra,vrra,brra,boprra,zra,idpotrb,vrrb,boprrb,brrb,zrb,
     > laypot,conc,dx,ns,rs,
     > idpota,vra,bra,bopra,za,idpotb,vrb,brb,boprb,zb,igraph,
c     > leftdelta,rightdelta,laydelta,
     > deltala,deltara,deltaa,rlambdala,rlambdara,rlambdaa,
     > deltalb,deltarb,deltab,rlambdalb,rlambdarb,rlambdab)
c     ------------------------------------------------------------------
c
c -read in data from update file if exists
c    -------------------------------------------------------------------
      call restart_in(itscf0,etop(npanel),v0,conc,opot,
     >                vra,bra,bopra,rba,vrb,brb,boprb,rbb,
     >                bulk,vrsbulk,nintfc,idpota,idpotb,za,zb,
     >                laycore,printout)
c    -------------------------------------------------------------------
      if (myrank.eq.root) then
      write(6,*)
      write(6,'(2x,''etop= '',f17.13)') etop(npanel)
      if(.not.bulk) write(6,*) 'v0=  ',v0
      end if
c
c -fix manipulations with files containing inverse cpa t-matrices
c     ------------------------------------------------------------------
      call tcpa_main(
     > leftmat,rightmat,laymat,bulk,dos,ksetmax,nintfc,ninprc(0),
     > ninprc(nprc+1),conc,concl,concr,
     > cpain,cpamatin,cpamatinl,cpamatinr)
c     ------------------------------------------------------------------
c
c -initialize common blocks
c
  410 continue
c     ------------------------------------------------------------------
      call initia
c     ------------------------------------------------------------------
      call cpu_time(tfinish)
      if (myrank.eq.root) then
      write(6,'('' Time ellapsed with init:'',t35,f8.1,'' s'')')
     > tfinish-tstart
      end if
c
c     *************************************************************
c               START SELFCONSISTENT ITERATIONS 
c     *************************************************************
c
      itscf=0
      istop=0
      efermi=etop(npanel)
  100 continue
#ifdef MPIP
      call mpi_barrier(mpi_comm_world,ierr)
#endif
      call cpu_time(tstart)
      itscf=itscf+1
      itscfcont=itscf0+itscf
      if (myrank.eq.root )
     >  write(6,'(/'' Selfconsistent iteration: '',i5/)') itscf
      efermi0=efermi
c
      if(clucalc.and.hostprint) goto 420 
c     ------------------------------------------------------------------
      call sublattpot(nintfc,ns,vra,iesublatt)
      call sublattpot(nintfc,ns,vrb,iesublatt)
      call sublattpot(nintfc,ns,bra,iesublatt)
      call sublattpot(nintfc,ns,brb,iesublatt)
      call sublattpot(nintfc,ns,bopra,iesublatt)
      call sublattpot(nintfc,ns,boprb,iesublatt)
c     call sublattpot(nintfc,ns,deltaa,iesublatt) !TODO: deltaa is complex, line commented by Laszloffy 07/07/2021
c     call sublattpot(nintfc,ns,deltab,iesublatt)
c     ------------------------------------------------------------------
      if(bulk) then
c     ------------------------------------------------------------------
         call averagebpot(nsl,nintfc,nbulkl,vra,vrla,vrra)
         call averagebpot(nsl,nintfc,nbulkl,vrb,vrlb,vrrb)
         call averagebpot(nsl,nintfc,nbulkl,bra,brla,brra)
         call averagebpot(nsl,nintfc,nbulkl,brb,brlb,brrb)
         do j=1,minprc               !!! if bulk, deltal and deltar are empty
           deltala(:,j)=deltaa(:,j)
           deltara(:,j)=deltaa(:,j)
           deltalb(:,j)=deltab(:,j)
           deltarb(:,j)=deltab(:,j)
         end do
c        bnyari: added average delta for scf calculations
c         call averagebpotc(nsl,nintfc,nbulkl,deltaa,detala,deltara)
c         call averagebpotc(nsl,nintfc,nbulkl,deltab,detalb,deltarb)
c     ------------------------------------------------------------------
         if(orbpol) then
c     ------------------------------------------------------------------
         call averagebpotsp(nsl,nintfc,nbulkl,bopra,boprla,boprra)
         call averagebpotsp(nsl,nintfc,nbulkl,boprb,boprlb,boprrb)
c     ------------------------------------------------------------------
         end if
      end if
c
c -solve Diraq equation for core states
c  (non spin-polarized!)
c
c     ------------------------------------------------------------------
      call chcore(itscfcont,nintfc,conc,idpota,vra,za,
     >idpotb,vrb,zb,dx,ns,rs,qca,qcb,rhoca,rhocb,enca,encb)
c     ------------------------------------------------------------------
c
c -generate new energy mesh in the case of bulk iterations
  420 continue
c
      if(itscf.eq.1.or.bulk) then
c     ------------------------------------------------------------------
         call zmesh(imesh,npanel,nepanel,ne1,ne2,ne3,ebottom,etop,
     &              eps,cear,we)
c     ------------------------------------------------------------------
      end if
      if(itest.ge.2.and.itscf.eq.1.and.myrank.eq.root) then 
        write(6,*)
        do ie=1,ne
          write(6,'('' ie='',i2,'' e='',2f12.8,'' w='',2f12.8,   
     >    ''  kset='',i5)') ie,cear(ie),we(ie),kset(ie)
        end do
        write(6,*)
      end if
c
      if(clucalc.and.hostprint) goto 400 
c
c -open binary files containing (inverse of) effective t-matrices
c     ------------------------------------------------------------------
      call tcpa_open(
     > tmatout,leftmat,rightmat,laymat,cpamatin,cpamatinl,cpamatinr)
c     ------------------------------------------------------------------
c
c - set quantities for contour integration to zero
c     ------------------------------------------------------------------
      call initzero(qvpa,qvpb,qva,qvha,qveha,qvhea,qvteha,qvthea,
     & qvb,qvhb,qvehb,qvheb,qvtehb,qvtheb,
     & spin_magvpa,spin_magvpb,spin_magva,spin_magvb,
     & orb_magvpa,orb_magvpb,orb_magva,orb_magvb,
     & enba,enbb,enorba,enorbb,qmoma,qmomb,
     & ddphenba,ddthenba,d2dphenba,d2dthenba,d2dthphenba,
     & ddphenbb,ddthenbb,d2dphenbb,d2dthenbb,d2dthphenbb,
     & rhova,rhoveha,rhovhea,rhovb,rhospa,rhospb,rhodspa,rhodspb,
     & rhomaga,rhomagb,enbifc,qvifc,omifc) 
c     ------------------------------------------------------------------
c
c Rotation matrices between local and global frames of reference
c
      calcder=.false.
      do li=1,nintfc
        if(dabs(rba(3,li)-1.0d0).lt.1.0d-06) then
          calcder=.false.
          goto 300
        end if
        if(dabs(rbb(3,li)-1.0d0).lt.1.0d-06) then
          calcder=.false.
          goto 300
        end if
      end do
 300  continue                       
      if(calcder) then
c     write(6,*) 'calcder',calcder
c     ------------------------------------------------------------------
      call localrot2(
     > lmax,nintfc,rbl,rbr,rba,rbb,vecna,vecnb,phia,phib,
     > dmatl,dmatlp,dmatr,dmatrp,dmata,dmatb,dmatpa,dmatpb,
     > ddpha,ddtha,d2dpha,d2dtha,d2dthpha,
     > ddphpa,ddthpa,d2dphpa,d2dthpa,d2dthphpa,
     > ddphb,ddthb,d2dphb,d2dthb,d2dthphb,
     > ddphpb,ddthpb,d2dphpb,d2dthpb,d2dthphpb,
     > rmata,rmatpa,rmatb,rmatpb)
c     ------------------------------------------------------------------
c      write(6,*) 'ddtha'
c      call outmat1(ddtha(1,1,1),kmymax,kmymax,kmymaxp,1.0d-10,6)
c      write(6,*) 'ddthpa'
c      call outmat1(ddthpa(1,1,1),kmymax,kmymax,kmymaxp,1.0d-10,6)
      else
c     write(6,*) 'calcder',calcder
      call czero(ddpha,kmymaxp*kmymaxp*mintfc)
      call czero(ddtha,kmymaxp*kmymaxp*mintfc)
      call czero(d2dpha,kmymaxp*kmymaxp*mintfc)
      call czero(d2dtha,kmymaxp*kmymaxp*mintfc)
      call czero(d2dthpha,kmymaxp*kmymaxp*mintfc)
      call czero(ddphpa,kmymaxp*kmymaxp*mintfc)
      call czero(ddthpa,kmymaxp*kmymaxp*mintfc)
      call czero(d2dphpa,kmymaxp*kmymaxp*mintfc)
      call czero(d2dthpa,kmymaxp*kmymaxp*mintfc)
      call czero(d2dthphpa,kmymaxp*kmymaxp*mintfc)
      call czero(ddphb,kmymaxp*kmymaxp*mintfc)
      call czero(ddthb,kmymaxp*kmymaxp*mintfc)
      call czero(d2dphb,kmymaxp*kmymaxp*mintfc)
      call czero(d2dthb,kmymaxp*kmymaxp*mintfc)
      call czero(d2dthphb,kmymaxp*kmymaxp*mintfc)
      call czero(ddphpb,kmymaxp*kmymaxp*mintfc)
      call czero(ddthpb,kmymaxp*kmymaxp*mintfc)
      call czero(d2dphpb,kmymaxp*kmymaxp*mintfc)
      call czero(d2dthpb,kmymaxp*kmymaxp*mintfc)
      call czero(d2dthphpb,kmymaxp*kmymaxp*mintfc)
c     ------------------------------------------------------------------
      call localrot(
     > lmax,nintfc,rbl,rbr,rba,rbb,vecna,vecnb,phia,phib,
     > dmatl,dmatlp,dmatr,dmatrp,dmata,dmatb,dmatpa,dmatpb,
     > rmata,rmatpa,rmatb,rmatpb)
c     ------------------------------------------------------------------
      end if
      call cpu_time(tfinish)
      if (myrank.eq.root) then
      write(6,'('' Time ellapsed with int scf:'',t35,f8.1,'' s'')')
     > tfinish-tstart
      end if
c
c                 ********************************
c                 *** Start layer energy loop ***
c                 ********************************
  400 continue 
c
      ie0=1
      iek=0
      cnet=0
c
      do 10 ie=1,ne
c
        if(kset(ie).eq.0) then
          itcpam=0
          intbz=0
          bsfcalc=.true.
        else 
          bsfcalc=.false.
        end if
c
      ce=cear(ie)
c
      if(clucalc.and.hostprint) then
        if (myrank.eq.root.and.ie.eq.1) then
          write(6,*) ' Skipping energy loop for host system.'
          write(6,*) ' Host t and tau matrices will be read from files.'
        end if
      else
c -read in (inverse of) effective t-matrices if necessary
c     ------------------------------------------------------------------
      call cpu_time(tstart)
      call tcpa_in(leftmat,rightmat,laymat,ce,nintfc,
     > ninprc(0),ninprc(nprc+1),conc,concl,concr,
     > cpamatin,cpamatinl,cpamatinr,2*kmymax,
     > tminv,tminvl,tminvr,
     > tau_kint,taua_kint,taub_kint,dbogomaxp)
      call cpu_time(tfinish)
c      write(6,'('' Time ellapsed in tcpa_in:'',t35,f8.1,'' s'')')
c     > tfinish-tstart
c     ------------------------------------------------------------------
c
c  generate k-mesh in irreducible Brillouin-zone
c
c     ------------------------------------------------------------------
      call cpu_time(tstart)
      call kmesh2d(intbz,kunning,kset(ie),park,xk(1,1),xk(1,2),wk,nk,
     &             mkpar)
      call cpu_time(tfinish)
c      write(6,'('' Time ellapsed in kmesh2d:'',t35,f8.1,'' s'')')
c     > tfinish-tstart
c     ------------------------------------------------------------------ 
c
      if(itest.ge.1 .and. myrank.eq.root) then
        write(6,'('' ie='',i3,''   e='',2f12.8,''  nk='',i5)') ie,ce,nk
        call flush(6)
      end if
      if(itest.ge.2.and. myrank.eq.root) then
        if(ie.gt.1) ie0=ie-1
        if((kset(ie).ne.kset(ie0)).or.(ie.eq.1)) then
          do ik=1,nk
            write(6,'(i5,3f14.10)') ik,xk(ik,1),xk(ik,2),wk(ik)
          end do
        end if
      end if
c
      call cpu_time(tstart)
c  t-matrices, ATA for CPA
c     ------------------------------------------------------------------ 
      call tmatini(
     > ce,lmax,nintfc,iscreen,vscreen,v0,E_Fermi,
     > bulk,wrel,sxl,sxr,sxa,sxb,
     > concl,concr,conc,cpamatin,cpamatinl,cpamatinr,
     > idpotla,vrla,brla,idpotlb,vrlb,brlb,rsl,dxl,nsl,
     > idpotra,vrra,brra,idpotrb,vrrb,brrb,rsr,dxr,nsr,
     > idpota,vra,bra,idpotb,vrb,brb,dx,ns,rs,
     > dmatl,dmatlp,dmatr,dmatrp,dmata,dmatb,dmatpa,dmatpb,
     > tminvl,tminvr,tminv,tminva,tminvb,ptminva,ptminvb,
     > rbl,rbr,rba,rbb,
     > deltala,deltara,deltaa,
     > deltalb,deltalb,deltab,c_light)
c     ------------------------------------------------------------------ 
      if (hostprint.and.myrank.eq.root) then
         call thout(ie,ce,nintfc,nbogomax,
     >              ptminva(1:nbogomax,1:nbogomax,1:nintfc),taupath)
      end if
      call cpu_time(tfinish)
      if (myrank.eq.root) 
     >  write(6,'('' Time ellapsed in tmatini:'',t35,f8.1,'' s'')')
     > tfinish-tstart
c
       tautest=.false.
       if(itest.ge.2.or.tautest.and.myrank.eq.root) then
        write(6,'(''Inverse screened t-matrix in lms representation'')')
        write(6,'(a)') idpotla(1)
       do li=1,nintfc
c
       write(6,'(i3,2x,a)') li,idpota(li)
c
       xlms=(0.0d0,0.0d0)
c
       call matlms(xlms,tminva(1,1,li),lmax)
       write(6,'(''Electron-electron Spin: 1 1'')') 
       call outmat1(xlms(:,:,1,1,1,1),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Electron-electron Spin: 2 2'')') 
       call outmat1(xlms(:,:,2,2,1,1),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Electron-electron Spin: 1 2'')') 
       call outmat1(xlms(:,:,1,2,1,1),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Electron-electron Spin: 2 1'')') 
       call outmat1(xlms(:,:,2,1,1,1),lmmax,lmmax,lmmaxp,tolerr,6)
c
       write(6,'(''Electron-hole Spin: 1 1'')') 
       call outmat1(xlms(:,:,1,1,1,2),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Electron-hole Spin: 2 2'')') 
       call outmat1(xlms(:,:,2,2,1,2),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Electron-hole Spin: 1 2'')') 
       call outmat1(xlms(:,:,1,2,1,2),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Electron-hole Spin: 2 1'')') 
       call outmat1(xlms(:,:,2,1,1,2),lmmax,lmmax,lmmaxp,tolerr,6)
c
       write(6,'(''Hole-electron Spin: 1 1'')') 
       call outmat1(xlms(:,:,1,1,2,1),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Hole-electron Spin: 2 2'')') 
       call outmat1(xlms(:,:,2,2,2,1),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Hole-electron Spin: 1 2'')') 
       call outmat1(xlms(:,:,1,2,2,1),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Hole-electron Spin: 2 1'')') 
       call outmat1(xlms(:,:,2,1,2,1),lmmax,lmmax,lmmaxp,tolerr,6)
c
       write(6,'(''Hole-hole Spin: 1 1'')') 
       call outmat1(xlms(:,:,1,1,2,2),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Hole-hole Spin: 2 2'')') 
       call outmat1(xlms(:,:,2,2,2,2),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Hole-hole Spin: 1 2'')') 
       call outmat1(xlms(:,:,1,2,2,2),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Hole-hole Spin: 2 1'')') 
       call outmat1(xlms(:,:,2,1,2,2),lmmax,lmmax,lmmaxp,tolerr,6)
c
       end do
       end if
c
c
c
      call cpu_time(tstart)
c  perform CPA, BZ integration & calculate tau-matices
c     ------------------------------------------------------------------ 
       if (.not. singlesite) then
c        if (cpamatin) then 
c          write(6,*) 'tau calculation for a single k point'
c          write(6,*) xk(nk,1:2)
c        end if
        call cpacoord(
     >   itscf,ie,ce,E_Fermi,lmax,nbogomax,
     >   nintfc,eta,rightm,bulk,bulkgeo,wrel,
     >   kset(ie),xk,wk,nk,intbz,iek,
     >   conc,itcpam,cpatol,cpatest,
     >   dmata,dmatb,dmatpa,dmatpb,
     >   tminvl,tminvr,tminv,tminva,tminvb,
     >   tau,taua,taub,gtaua,gtaub,c_light,
c  --------------- cluster -----------------
     >   kpair,npair,npair0,npair1,tau_ii1,tau_ij1,tau_ji1,taupath,
     >   clucalc,hostprint)
c
      else
        write(6,*) 'Skipping cpacoord, tmat used instead of tau, ie = ',
     >  ie
        do li=1,nintfc
c           ! physical tminv^{-1} => physical tau
            tempmat = ptminva(:,:,li)
            call gjinv(tempmat,2*kmymax,dbogomaxp,detl)
            taua(:,:,li) = tempmat
            write(6,*) 't-matrix in kappamy representation'
            call outmat1(tempmat(:,:),
     >            2*kmymax,2*kmymax,dbogomaxp,1.0d-8,6)
         end do
       end if
c     ------------------------------------------------------------------ 
c
      call cpu_time(tfinish)
      if (myrank.eq.root) then
       write(6,'('' Time ellapsed in cpacoord:'',t35,f8.1,'' s'')')
     > tfinish-tstart
      end if
      tautest=.false.
      if(tautest.and.myrank.eq.root) then
       do li=1,nintfc
        write(6,*) ' <skkr> : tautest, ptminv li=',li
        call outmat1(ptminva(1,1,li),nbogomax,nbogomax,dbogomaxp,
     >               tol,6)
        write(6,*) ' <skkr> : tautest, tminv li=',li
        call outmat1(tminva(1,1,li),nbogomax,nbogomax,dbogomaxp,
     >               tol,6)
       end do
      end if
c
       tautest=.false.
       if(itest.ge.2.or.tautest) then
       if (.not. singlesite) then
          write(6,'(''Physical Tau in lms representation'')')
       else
          write(6,'(''t-matrix in lms representation'')')
       end if
c
       do li=1,nintfc
c
       write(6,'(i3,2x,a)') li,idpota(li)
c
       xlms=(0.0d0,0.0d0)
c
       call matlms(xlms,taua(1,1,li),lmax)
       write(6,'(''Electron-electron Spin: 1 1'')') 
       call outmat1(xlms(:,:,1,1,1,1),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Electron-electron Spin: 2 2'')') 
       call outmat1(xlms(:,:,2,2,1,1),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Electron-electron Spin: 1 2'')') 
       call outmat1(xlms(:,:,1,2,1,1),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Electron-electron Spin: 2 1'')') 
       call outmat1(xlms(:,:,2,1,1,1),lmmax,lmmax,lmmaxp,tolerr,6)
c
       write(6,'(''Electron-hole Spin: 1 1'')') 
       call outmat1(xlms(:,:,1,1,1,2),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Electron-hole Spin: 2 2'')') 
       call outmat1(xlms(:,:,2,2,1,2),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Electron-hole Spin: 1 2'')') 
       call outmat1(xlms(:,:,1,2,1,2),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Electron-hole Spin: 2 1'')') 
       call outmat1(xlms(:,:,2,1,1,2),lmmax,lmmax,lmmaxp,tolerr,6)
c
       write(6,'(''Hole-electron Spin: 1 1'')') 
       call outmat1(xlms(:,:,1,1,2,1),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Hole-electron Spin: 2 2'')') 
       call outmat1(xlms(:,:,2,2,2,1),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Hole-electron Spin: 1 2'')') 
       call outmat1(xlms(:,:,1,2,2,1),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Hole-electron Spin: 2 1'')') 
       call outmat1(xlms(:,:,2,1,2,1),lmmax,lmmax,lmmaxp,tolerr,6)
c
       write(6,'(''Hole-hole Spin: 1 1'')') 
       call outmat1(xlms(:,:,1,1,2,2),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Hole-hole Spin: 2 2'')') 
       call outmat1(xlms(:,:,2,2,2,2),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Hole-hole Spin: 1 2'')') 
       call outmat1(xlms(:,:,1,2,2,2),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Hole-hole Spin: 2 1'')') 
       call outmat1(xlms(:,:,2,1,2,2),lmmax,lmmax,lmmaxp,tolerr,6)
c
       end do
       write(6,'(''Physical Tau ij offdiagonal block
     > in lms representation'')')
       do li=1,npair
c
       write(6,'(i3)') li
c
       xlms=(0.0d0,0.0d0)
c
       call matlms(xlms,tau_ij1(1,1,li),lmax)
       write(6,'(''Electron-electron Spin: 1 1'')') 
       call outmat1(xlms(:,:,1,1,1,1),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Electron-electron Spin: 2 2'')') 
       call outmat1(xlms(:,:,2,2,1,1),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Electron-electron Spin: 1 2'')') 
       call outmat1(xlms(:,:,1,2,1,1),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Electron-electron Spin: 2 1'')') 
       call outmat1(xlms(:,:,2,1,1,1),lmmax,lmmax,lmmaxp,tolerr,6)
c
       write(6,'(''Electron-hole Spin: 1 1'')') 
       call outmat1(xlms(:,:,1,1,1,2),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Electron-hole Spin: 2 2'')') 
       call outmat1(xlms(:,:,2,2,1,2),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Electron-hole Spin: 1 2'')') 
       call outmat1(xlms(:,:,1,2,1,2),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Electron-hole Spin: 2 1'')') 
       call outmat1(xlms(:,:,2,1,1,2),lmmax,lmmax,lmmaxp,tolerr,6)
c
       write(6,'(''Hole-electron Spin: 1 1'')') 
       call outmat1(xlms(:,:,1,1,2,1),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Hole-electron Spin: 2 2'')') 
       call outmat1(xlms(:,:,2,2,2,1),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Hole-electron Spin: 1 2'')') 
       call outmat1(xlms(:,:,1,2,2,1),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Hole-electron Spin: 2 1'')') 
       call outmat1(xlms(:,:,2,1,2,1),lmmax,lmmax,lmmaxp,tolerr,6)
c
       write(6,'(''Hole-hole Spin: 1 1'')') 
       call outmat1(xlms(:,:,1,1,2,2),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Hole-hole Spin: 2 2'')') 
       call outmat1(xlms(:,:,2,2,2,2),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Hole-hole Spin: 1 2'')') 
       call outmat1(xlms(:,:,1,2,2,2),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Hole-hole Spin: 2 1'')') 
       call outmat1(xlms(:,:,2,1,2,2),lmmax,lmmax,lmmaxp,tolerr,6)
c
       end do
       write(6,'(''Physical Tau ji offdiagonal block
     > in lms representation'')')
       do li=1,npair
c
       write(6,'(i3)') li
c
       xlms=(0.0d0,0.0d0)
c
       call matlms(xlms,tau_ji1(1,1,li),lmax)
       write(6,'(''Electron-electron Spin: 1 1'')') 
       call outmat1(xlms(:,:,1,1,1,1),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Electron-electron Spin: 2 2'')') 
       call outmat1(xlms(:,:,2,2,1,1),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Electron-electron Spin: 1 2'')') 
       call outmat1(xlms(:,:,1,2,1,1),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Electron-electron Spin: 2 1'')') 
       call outmat1(xlms(:,:,2,1,1,1),lmmax,lmmax,lmmaxp,tolerr,6)
c
       write(6,'(''Electron-hole Spin: 1 1'')') 
       call outmat1(xlms(:,:,1,1,1,2),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Electron-hole Spin: 2 2'')') 
       call outmat1(xlms(:,:,2,2,1,2),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Electron-hole Spin: 1 2'')') 
       call outmat1(xlms(:,:,1,2,1,2),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Electron-hole Spin: 2 1'')') 
       call outmat1(xlms(:,:,2,1,1,2),lmmax,lmmax,lmmaxp,tolerr,6)
c
       write(6,'(''Hole-electron Spin: 1 1'')') 
       call outmat1(xlms(:,:,1,1,2,1),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Hole-electron Spin: 2 2'')') 
       call outmat1(xlms(:,:,2,2,2,1),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Hole-electron Spin: 1 2'')') 
       call outmat1(xlms(:,:,1,2,2,1),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Hole-electron Spin: 2 1'')') 
       call outmat1(xlms(:,:,2,1,2,1),lmmax,lmmax,lmmaxp,tolerr,6)
c
       write(6,'(''Hole-hole Spin: 1 1'')') 
       call outmat1(xlms(:,:,1,1,2,2),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Hole-hole Spin: 2 2'')') 
       call outmat1(xlms(:,:,2,2,2,2),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Hole-hole Spin: 1 2'')') 
       call outmat1(xlms(:,:,1,2,2,2),lmmax,lmmax,lmmaxp,tolerr,6)
       write(6,'(''Hole-hole Spin: 2 1'')') 
       call outmat1(xlms(:,:,2,1,2,2),lmmax,lmmax,lmmaxp,tolerr,6)
c
       end do
       end if
c
c

       call cpu_time(tstart)
       if (myrank.eq.root) then
       if (itcpam.gt.0) then
c -write out (inverse of) effective t-matrices
c      ------------------------------------------------------------------ 
          call tcpa_out(ce,nintfc,conc,2*kmymax,
     >    tminv,tau,taua,taub,dbogomaxp)
c      ------------------------------------------------------------------ 
c          stop 'after tcpa_out'
       end if
       end if 
       call cpu_time(tfinish)
c       write(6,'('' Time ellapsed in tcpa_out:'',t35,f8.1,'' s'')')
c     > tfinish-tstart
c
c
c      stop 'The code is stopped before locquant'
c
c
c
c  calculate local physical quantities
c  integrate with respect to energy
      skiploq=.false.
c     ------------------------------------------------------------------
      if (.not.bsfcalc) then
      call cpu_time(tstart)
      if (.not.skiploq) then
      call locquant(
     > ie,ce,we(ie),dos,lmax,nintfc,wrel,lms,sxa,sxb,conc,v0,
     > E_Fermi,deltaa,deltab,
     > idpota,vra,bra,rba,idpotb,vrb,brb,rbb,dx,ns,rs,
c     > rmata,rmatpa,rmatb,rmatpb,
     > taua,taub,gtaua,gtaub,
c     > ptminva,ptminvb,tminv,
     > dosa(1,1,ie),dosha(1,1,ie),doseha(1,1,ie),doshea(1,1,ie),
     > dosteha(1,1,ie),dosthea(1,1,ie),
     > qvpa,qva,qvha,qveha,qvhea,qvteha,qvthea,
     > qvdiffa(1,ie),enba,enbdiffa(1,ie),qmoma,
c     > ddphenba,ddthenba,d2dphenba,d2dthenba,d2dthphenba,
     > rhova,rhoveha,rhovhea,rhospa,rhodspa,rhomaga,
     > dosmaga(1,1,ie),dosmagha(1,1,ie),
     > spin_magvpa,spin_magva,orb_magvpa,orb_magva,
c quantites for B component
     > dosb(1,1,ie),doshb(1,1,ie),dosehb(1,1,ie),dosheb(1,1,ie),
     > dostehb(1,1,ie),dostheb(1,1,ie),
     > qvpb,qvb,qvhb,qvehb,qvheb,qvtehb,qvtheb,
     > qvdiffb(1,ie),enbb,enbdiffb(1,ie),qmomb,
     > rhovb,rhovehb,rhovheb,rhospb,rhodspb,rhomagb,
     > dosmagb(1,1,ie),dosmaghb(1,1,ie),
     > spin_magvpb,spin_magvb,orb_magvpb,orb_magvb,
     > lliter,linbw,efermi,enbifc,qvifc,omifc,reg,c_light)
      end if
      call cpu_time(tfinish)
      if (myrank.eq.root) then
      write(6,'('' Time ellapsed in locquant:'',t35,f8.1,'' s'')')
     > tfinish-tstart
      end if
      else 
      call cpu_time(tstart)
      call bsf(
     > ie,ce,we(ie),dos,lmax,nintfc,wrel,lms,sxa,sxb,conc,v0,
     > E_Fermi,deltaa,deltab,
     > idpota,vra,bra,rba,idpotb,vrb,brb,rbb,dx,ns,rs,
     > tminv,tau,tau_kint,taua_kint,taub_kint,gtaua,gtaub,
c we store the bsf values in the dos variables
     > dosa(1,:,ie),dosha(1,:,ie),doseha(1,:,ie),doshea(1,:,ie),
     > dosteha(1,:,ie),dosthea(1,:,ie),
     > lliter,linbw,efermi,enbifc,qvifc,omifc,reg,c_light)
      call cpu_time(tfinish)
      if (myrank.eq.root) then
      write(6,'('' Time ellapsed in bsf:'',t35,f8.1,'' s'')')
     > tfinish-tstart
      end if
      end if 
      end if ! end of if(clucalc.and.hostprint)
c     ------------------------------------------------------------------
c
c      stop 'The code is stopped after locquant'
c                                    
       if(clucalc) then
         if(mod(ie-1,nprocs).eq.myrank) then
           cnet=cnet+1
           ze=cear(ie)
           if(hostprint) then
c call with 2*kmymax instead of kmymax
           call tauhin(ie,ze,npair,npair1,kpair,nimp,2*kmymax,nintfc,
     >                 tau_ij(1,1,1,cnet),
     >                 tau_ji(1,1,1,cnet),
     >                 tau_ii(1,1,1,cnet),taupath)!,imesh)
           call thin(ie,ze,nimp,2*kmymax,nintfc,tminvh(1,1,1,cnet),
     >               taupath)
           write(*,*) 'after read tauhin',myrank,ie
           else
             tminvh(:,:,:,cnet)=ptminva(1:nbogomax,1:nbogomax,1:nintfc)
             tau_ij(:,:,:,cnet)=tau_ij1
             tau_ji(:,:,:,cnet)=tau_ji1
             tau_ii(:,:,:,cnet)=tau_ii1
c TODO: if cluster calculations, copy tau_ii1, tau_ij1 and tau_ji1 to the
c       corresponding arrays (tau_ii[...cnet]=tau_ii1 etc.).
c       Unscreened host tminv is also needed from cpacoord.
           end if
         endif
       end if
 10   continue
c
c                 **************************
c                 *** End of energy loop ***
c                 **************************
c
      if (myrank.eq.root) then
      if (bsfcalc) then
      write(6,*) 'The result of the BSF calculation'
      write(6,*) 'li  Re(E)  Im(E)  kx  ky  bsf e,h,s,t0,t1,t-1'
      do ie=1,ne
      do li=1,nintfc
        write(6,'(i3,2f12.8,f12.5,f12.5,6d14.6)')
     >       li,cear(ie),xk(nk,1:2),dosa(1,li,ie),
     >        dosha(1,li,ie),doseha(1,li,ie),
     >        doshea(1,li,ie),dosteha(1,li,ie),dosthea(1,li,ie)
c        write(6,'(''e '')') ce, xk(nk,1:2), dosa(1,li,ie), dosha(1,li,ie)
      end do
      end do
      call flush(6)
      stop 'bsf calculation is finished'
      end if
      end if
c 
c -close binary files containing (inverse of) effective t-matrices
c     ------------------------------------------------------------------
      call tcpa_close(
     >leftmat,rightmat,laymat,cpamatin,cpamatinl,cpamatinr)
c     ------------------------------------------------------------------
c
c  update spin-quantization axes
c     ------------------------------------------------------------------
      call newdir(
     > itscf,itscfcont,nintfc,lfix1,lfix2,
     > rba,rbb,vecna,phia,vecnb,phib,lmax,
     > spin_magvpa,spin_magvpb,spin_magva,spin_magvb,
     > orb_magvpa,orb_magvpb,orb_magva,orb_magvb,lza,lzb,
     > spinmoma,orbmoma,th0a,th1a,ph0a,ph1a,
     > spinmomb,orbmomb,th0b,th1b,ph0b,ph1b)
c     ------------------------------------------------------------------
c
c  write out results
c     ------------------------------------------------------------------
      if(myrank.eq.root) then
      call printscf(
     > linbw,nintfc,itscfcont,conc,enbifc,qvifc,omifc,
     > qva,qvha,qveha,qvhea,qvteha,qvthea,
     > qvb,qvhb,qvehb,qvheb,qvtehb,qvtheb,
     > qca,za,qcb,zb,enba,enbb,efermi,
     > ddphenba,ddthenba,d2dphenba,d2dthenba,d2dthphenba,
     > ddphenbb,ddthenbb,d2dphenbb,d2dthenbb,d2dthphenbb,
     > arot,rba,rbb,spin_magva,spin_magvb,orb_magva,orb_magvb,
     > spinmoma,orbmoma,th0a,th1a,ph0a,ph1a,
     > spinmomb,orbmomb,th0b,th1b,ph0b,ph1b)
       write(6,*)
       end if
c bnyari May 2022 delta update
       do li=1,nintfc
           deltaold = (0.0d0,0.0d0)
           deltaold(:)=deltaa(:,li)
           do j=1,nrad
           deltaa(j,li) = 2*dreal(qveha(li))*rlambdaa(li)*cmix
     >                  +(1-cmix)*deltaold(j)
           end do
           if (1-conc(li).gt.tiny) then
               deltaold = (0.0d0,0.0d0)
               deltaold(:)=deltaa(:,li)
               do j=1,nrad
               deltab(j,li) = dreal(qvehb(li))*rlambdab(li)*cmix
     >                      +(1-cmix)*deltaold(j)
               end do
           end if
       end do 
      if (myrank.eq.root) then 
          do li=1,nintfc
            if (1-conc(li).lt.tiny) then
              deltaprint=dreal(deltaa(1,li))*13605.662285137
          write(6,'('' I'',i4,''  L'',i4,'' Delta'',e18.10,'' meV'')')
     >               itscfcont,li,deltaprint 
            else
              deltaprint=dreal(deltaa(1,li))*13605.662285137
          write(6,'('' I'',i4,''  L'',i4,'' DeltaA'',e18.10,'' meV'')')
     >               itscfcont,li,deltaprint
              deltaprint=dreal(deltab(1,li))*13605.662285137
          write(6,'('' I'',i4,''  L'',i4,'' DeltaB'',e18.10,'' meV'')')
     >               itscfcont,li,deltaprint
            end if
          end do
          write(6,*) 
      end if
      call flush(6)
c     ------------------------------------------------------------------
c
      if(linbw) goto 200
c for DOS calculation print out results and stop
      if(dos) goto 99
c
c -find new Fermi level --> fixed by normal state calculation
c
c      if(bulk) then
c        kmy0=1
c     ------------------------------------------------------------------
c        call newfl(itscfcont,kmy0,kmymax,kmymaxp,ne,nbulkl,conc,
c     >             dosa,dosb,qca,qcb,qva,qvb,za,zb,efermi,defermi)
c     ------------------------------------------------------------------
c      end if
c
      if (myrank.eq.root) then
c      call deltamix(itscf,itscfcont,nintfc,lmax,iesublatt,
c     > dx,ns,rs,ferr1,rlambda,delta,qveha,qvhea)
c
c -calculate and print out energies
      entotifc =0.d0
c      write(6,*) 'before enprint', entotifc
c     ------------------------------------------------------------------
      call enprint(itscfcont,nintfc,conc,enca,encb,enba,enbb,
     >  enela,enelb,enxca,enxcb,enpota,enpotb,enmaga,enmagb,
     >  enorba,enorbb,entota,entotb,entot,entotifc)
c     ------------------------------------------------------------------
      end if
c      call delta_bcast(ferr1,rlambda,delta,qveha,qvhea,entotifc)
c
      if(itscf.eq.1) then
        dentotifc=10.0d0*tolen
      else
        dentotifc=entotifc-entotifc0
      end if
      entotifc0=entotifc
c
c -update output necessary for restart
c     ------------------------------------------------------------------
c      call restart_out(itscfcont,efermi,v0,conc,orbpol,
c     >                 vra,bra,bopra,rba,vrb,brb,boprb,rbb,
c     >                 bulk,vrsbulk,nintfc,printout)
c     ------------------------------------------------------------------
c      if(bulk) etop(npanel)=efermi
c
       istop=0
c       if(ferr1.lt.tolerr) istop=1
c      if(bulk.and.(dabs(defermi).gt.tolef)) istop=0
c      if(bulk.and.(dabs(dentotifc).gt.nintfc*tolen)) istop=0
c     if(rightm.eq.'V'.and.dv0.gt.tolv0) istop=0
c
c -Print results
c
   99 continue
      call flush(6)
c     ------------------------------------------------------------------
      if (myrank.eq.root) then
      call flush(6)
      call printres(
     >     dosout,pdosout,potout,momout,
     >     itscfcont,lmax,nintfc,ne,cear,efermi0,
     >     conc,dosa,dosb,dosha,doshb,
     >     doseha,doshea,dosteha,dosthea,
     >     dosehb,dosheb,dostehb,dostheb,
     >     qvpa,qvpb,qva,qvb,dosmaga,dosmagb,dosmagha,
     >     spin_magvpa,spin_magvpb,spin_magva,spin_magvb,
     >     orb_magvpa,orb_magvpb,orb_magva,orb_magvb,
     >     arot,qmoma,qmomb,za,zb,qca,qcb,
     >     qvdiffa,qvdiffb,enbdiffa,enbdiffb,
     >     entota,entotb,entot,entotifc,
     >     idpota,idpotb,vra,bra,bopra,vrb,brb,boprb,
     >     rs,dx,ns,vrsbulk,lms,orbpol,bulk,
     >     qveha,qvehb,rlambdaa,rlambdab)
      end if
c
c     ------------------------------------------------------------------
c
c
      if(clucalc.or.hostprint) then
        if(itscf.lt.itscfmax.and.istop.eq.0) goto 100
        if(istop.eq.1) write(6,'(/" SCF cycle converged!")')
      else 
        if(itscf.lt.itscfmax.and.istop.eq.0) goto 100
        if(istop.eq.1 .and. myrank .eq. root) 
     >          write(6,'(/" SCF cycle converged!")')
      end if
c      if(itscf.lt.itscfmax) goto 100
c
c     ****************************************************************
c                    END OF SELFCONSISTENT ITERATIONS
c     ****************************************************************
c
cc TODO: deallocate arrays for layer calculations (e.g. tau_ii1 if exsists)
      deallocate(tau_ii1,stat=AllocateStatus)
      call alloccheck( AllocateStatus,
     >     'tauii1 dealloc in skkr                            ' )
      deallocate(tau_ij1,stat=AllocateStatus)
      call alloccheck( AllocateStatus,
     >     'tauij1 dealloc in skkr                            ' )
      deallocate(tau_ji1,stat=AllocateStatus)
      call alloccheck( AllocateStatus,
     >     'tauji1 dealloc in skkr                            ' )
c
c     ***************************************************************
c                       INITIALIZE CLUSTER CALCULATION
c     ***************************************************************
c     ------------------------------------------------------------------
      if(clucalc) then ! cluster case
        if (myrank.eq.root) write(6,*) '<skkr> initialize cluster
     > calculation'
       call flush(6)
c     Allocate potential related variables and read cluster potentials 
        ALLOCATE ( idpotimp(nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'idpotimp in skkr                                  ' )
        ALLOCATE ( vrimp(nrad,nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'vrimp in skkr                                     ' )
        ALLOCATE ( brimp(nrad,nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'brimp in skkr                                     ' )
        ALLOCATE ( boprimp(nrad,2,nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'boprimp in skkr                                   ' )
        ALLOCATE ( deltaimp(nrad,nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'deltaimp in skkr                                  ' )
        ALLOCATE ( rlambdaimp(nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'rlambdaimp in skkr                                ' )
        ALLOCATE ( singratimp(nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'singratimp in skkr                                ' )
        ALLOCATE ( uratimp(nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'uratimp in skkr                                   ' )
        ALLOCATE ( dratimp(nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'dratimp in skkr                                   ' )
        ALLOCATE ( dximp(nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'dximp in skkr                                     ' )
        ALLOCATE ( nsimp(nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'nsimp in skkr                                     ' )
        ALLOCATE ( rsimp(nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'rsimp in skkr                                     ' )
        ALLOCATE ( zimp(nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'zimp in skkr                                      ' )
c  
        do iimp=1,nimp
          rsimp(iimp) = rs(nposimp(3,iimp))
        end do
        rlambdaimp = (0.0d0,0.0d0) 
        call pothandleimp(
     >  lmax,nimp,v00,ivacpot,opot,imppot,dximp,nsimp,rsimp,idpotimp,
     >  vrimp,brimp,boprimp,zimp,ib0,b0,impb0,impb0f,isigb,
     >  deltaimp,rlambdaimp,singratimp,uratimp,dratimp)
c     ------------------------------------------------------------------
c allocate cluster scattering matrices
        ALLOCATE ( tminvcl(nbogomax,nbogomax,nimp),
     >                  STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'tminvcl in skkr                                   ' )
        ALLOCATE ( ptminvcl(nbogomax,nbogomax,nimp),
     >                  STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'ptminvcl in skkr                                  ' )
c     ------------------------------------------------------------------
c allocate local quantites
        ALLOCATE ( dosimp(kmymaxp,nimp,me),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'dosimp in skkr                                    ' )
        ALLOCATE ( doshimp(kmymaxp,nimp,me),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'doshimp in skkr                                   ' )
        ALLOCATE ( dosehimp(kmymaxp,nimp,me),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'dosehimp in skkr                                  ' )
        ALLOCATE ( dosheimp(kmymaxp,nimp,me),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'dosheimp in skkr                                  ' )
        ALLOCATE ( dostehimp(kmymaxp,nimp,me),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'dostehimp in skkr                                 ' )
        ALLOCATE ( dostheimp(kmymaxp,nimp,me),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'dostheimp in skkr                                 ' )
        ALLOCATE ( qvimp(nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'qvimp in skkr                                     ' )
        ALLOCATE ( qvpimp(kmymaxp,nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'qvpimp in skkr                                    ' )
        ALLOCATE ( qvhpimp(kmymaxp,nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'qvhpimp in skkr                                   ' )
        ALLOCATE ( qvhimp(nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'qvhimp in skkr                                    ' )
        ALLOCATE ( qvehimp(nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'qvehimp in skkr                                   ' )
        ALLOCATE ( qvheimp(nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'qvheimp in skkr                                   ' )
        ALLOCATE ( qvtehimp(nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'qvtehimp in skkr                                  ' )
        ALLOCATE ( qvtheimp(nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'qvtheimp in skkr                                  ' )
        ALLOCATE ( qmomimpa(lmsup,nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'qmomimpa in skkr                                  ' )
        ALLOCATE ( qmomimpha(lmsup,nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'qmomimpha in skkr                                 ' )
        ALLOCATE ( qmomimpeha(nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'qmomimpeha in skkr                                ' )
        ALLOCATE ( qmomimphea(nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'qmomimphea in skkr                                ' )
        ALLOCATE ( qmomimpteha(nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'qmomimpteha in skkr                               ' )
        ALLOCATE ( qmomimpthea(nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'qmomimpthea in skkr                               ' )
        ALLOCATE ( vmadih(nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'vmadih in skkr                                    ' )
        ALLOCATE ( vmadich(nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'vmadich in skkr                                   ' )
        ALLOCATE ( vmadid(nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'vmadid in skkr                                    ' )
        ALLOCATE ( vmadiq(nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'vmadiq in skkr                                    ' )
        ALLOCATE ( enbimp(nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'enbimp in skkr                                    ' )
        ALLOCATE ( enbhimp(nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'enbhimp in skkr                                   ' )
        ALLOCATE ( enbdiffimp(nimp,me),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'enbdiffimp in skkr                                ' )
        ALLOCATE ( enbhdiffimp(nimp,me),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'enbhdiffimp in skkr                               ' )
        ALLOCATE ( rhovimp(nrad,nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'rhovimp in skkr                                   ' )
        ALLOCATE ( rhospimp(nrad,2,nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'rhospimp in skkr                                  ' )
        ALLOCATE ( rhodspimp(nrad,2,nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'rhodspimp in skkr                                 ' )
        ALLOCATE ( rhomagimp(nrad,nimp),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'rhomagimp in skkr                                 ' )
        ALLOCATE ( dosmagimp(kmymaxp,nimp,me),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'dosmagimp in skkr                                 ' )
        ALLOCATE ( dosmaghimp(kmymaxp,nimp,me),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'dosmaghimp in skkr                                ' )
        ALLOCATE ( spin_magvpimp(kmymaxp,nimp,3),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'spin_magvpimp in skkr                             ' )
        ALLOCATE ( spin_magvimp(nimp,3),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'spin_magvimp in skkr                              ' )
        ALLOCATE (spin_magvhpimp(kmymaxp,nimp,3),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'spin_magvhpimp in skkr                            ' )
        ALLOCATE ( spin_magvhimp(nimp,3),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'spin_magvhimp in skkr                             ' )
        ALLOCATE ( orb_magvpimp(kmymaxp,nimp,3),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'orb_magvpimp in skkr                              ' )
        ALLOCATE ( orb_magvimp(nimp,3),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'orb_magvimp in skkr                               ' )
        ALLOCATE ( orb_magvhpimp(kmymaxp,nimp,3),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'orb_magvhpimp in skkr                             ' )
        ALLOCATE ( orb_magvhimp(nimp,3),STAT = AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'orb_magvhimp in skkr                              ' )
        ALLOCATE ( tripletmat(lmmaxp,lmmaxp,3,nimp,me),STAT = 
     >  AllocateStatus)
        call alloccheck( AllocateStatus,
     >       'tripletmat in skkr                                ' )
      end if
c TODO: initzero for local quantities
c     ------------------------------------------------------------------
cc =====================================================================
cc
c
c            *************************************
c            *** Start energy loop for cluster ***
c            *************************************
c
      if(clucalc) then
      itscf=0
  101 continue
      itscf=itscf+1
      itscfcont=itscf0+itscf
      write(6,'(/'' Selfconsistent iteration: '',i3/)') itscf
      ie0=1
      iek=0
      cnet=0
c
c Initialize locquant quantities
c due to collect_data, some variables (e.g. dosimp) will not set directly
c
      call initzero1(kmymax,nimp,ne,! TODO: why is qva complex and qvimp
! real?
     > dosimp,doshimp,dosehimp,
     > dosheimp,dostehimp,dostheimp,
     > qvpimp,qvhpimp,qvimp,qvhimp,qvehimp,qvheimp,
     > qvtehimp,qvtheimp,
     > qmomimpha,qmomimpeha,qmomimphea,qmomimpteha,qmomimpthea)
c
      if (myrank.eq.root) write(6,*) 'Start energy loop for cluster'
      do 11 ie=1,ne
      if(mod(ie-1,nprocs).eq.myrank) then
      cnet=cnet+1
      ce=cear(ie)
      if (myrank.eq.root) then
        write(6,*)
        do iie=1,ne
          write(6,'('' ie='',i3,'' e='',2f12.8,'' w='',2f12.8,   
     >    ''  kset='',i5)') iie,cear(iie),we(iie),kset(iie)
        end do
        write(6,*)
      end if
      call flush(6)
c      call mpi_barrier(mpi_comm_world,ierr)
      write(*, '(''energy point  '',i4,'' on process '',i2,'' cnet= '',
     > i2,'' rnet= '',i2)') ie,myrank,cnet,rnet
      if (myrank.eq. root) then
      write(6,*) 'energy point',ie,'on process',myrank,'cnet',cnet
      call flush(6)
      end if
c
c check host t and tau matrices
c minor change in tauhin fixed some 0 blocks in tauii
c nimp -> nintfc changed in the routine
c
      tautest=.false.
      if(tautest) then
       do li=1,nintfc
        write(6,*) ' <skkr> : tautest, tm li=',li
        call outmat1(tminvh(1,1,li,cnet),nbogomax,nbogomax,dbogomaxp,
     >               tol,6)
        write(6,*) ' <skkr> : tautest, tauii li=',li
        call outmat1(tau_ii(1,1,li,cnet),nbogomax,nbogomax,dbogomaxp,
     >               tol,6)
       end do
       do li = 1,npair
        write(6,*) ' <skkr> : tautest, tauji ipair=',li
        call outmat1(tau_ji(1,1,li,cnet),nbogomax,nbogomax,dbogomaxp,
     >               tol,6)
        write(6,*) ' <skkr> : tautest, tauij ipair=',li
        call outmat1(tau_ij(1,1,li,cnet),nbogomax,nbogomax,dbogomaxp,
     >               tol,6)
       end do
      end if
c
      tautest=.false.
      if(tautest) then
       write(6,*) "ce=",ce
       write(6,*) "wrel=",wrel
       write(6,*) "sxcl=",sxcl
       write(6,*) "idpotimp=",idpotimp
       write(6,*) "vrimp=",vrimp
       write(6,*) "brimp=",brimp
       write(6,*) "boprimp=",boprimp
       write(6,*) "dximp=",dximp
       write(6,*) "nsimp=",nsimp
       write(6,*) "rsimp=",rsimp
       write(6,*) "rbcl=",rbcl
       write(6,*) "deltaimp=",deltaimp
       write(6,*) "E_Fermi=",E_Fermi
       write(6,*) "c_light=",c_light
      end if
c
c
cc TODO: implement methods below
cc  calculate cluster t matrices
cc =====================================================================
      write(6,*) 'before tmatiniimp'    
      call tmatiniimp(
     > ce,lmax,nbogomax,nintfc,nimp,nposimp,0,vscreen,v0,
     > wrel,sxcl,idpotimp,vrimp,brimp,boprimp,dximp,nsimp,rsimp,
c    > dmata,dmatpa,
     > tminvcl,ptminvcl,!ptminv0,
     > rbcl,deltaimp,E_Fermi,c_light,singratimp,uratimp,dratimp)
      print *, 'after tmatiniimp',myrank   
      write(6,*) 'after tmatiniimp'
c
c
      tautest=.false.
      if(tautest) then
       do li=1,nimp
        write(6,*) ' <skkr> : tautest, tm cl li=',li
        call outmat1(tminvcl(1,1,li),nbogomax,nbogomax,nbogomax,
     >               tol,6)
        write(6,*) ' <skkr> : tautest, ptm cl li=',li
        call outmat1(ptminvcl(1,1,li),nbogomax,nbogomax,nbogomax,
     >               tol,6)
       end do
      end if
c
      write(6,*) 'before ecoreimp'    
c
      tautest=.false.
      if(tautest) then
       do li=1,nintfc
        write(6,*) ' <skkr> : tautest, tau_ii li=',li
        call outmat1(tau_ii(1:kmymax,1:kmymax,li,cnet),
     >               kmymax,kmymax,kmymax,
     >               tol,6)
       end do
      end if
c
cc  calculate cluster tau-matrices
cc =====================================================================
      call ecoreimp(nimp,npair,npair1,nintfc,0,! No screening for cluster calculations
     > tau_ij(1,1,1,cnet),tau_ji(1,1,1,cnet),tau_ii(1,1,1,cnet),
     > tminvcl,tminvh(1,1,1,cnet),lmax,nbogomax,nposimp,kpairind,
     > taucl,detl)!,cnet,rnet)
c
      print *, 'after ecoreimp',myrank    
      write(6,*) 'after ecoreimp'    
c
      tautest=.false.
      if(tautest) then
       do li=1,nimp
        write(6,*) ' <skkr> : tautest, taucl iimp=',li
        call outmat1(taucl(1:kmymax,1:kmymax,li),kmymax,kmymax,kmymax,
     >               tol,6)
       end do
      end if
c
cc  calculate local physical quantities
cc  integrate with respect to energy
cc =====================================================================
c remove spin orbital moments, band energy, charge moments, Madelung 
c     write(6,*) 'before locquant2'    
       call locquant2(
     > ie,ce,we(ie),lmax,madmax,nimp,wrel,lms,sxcl,v0,
     > deltaimp,singratimp,uratimp,dratimp,
     > idpotimp,vrimp,brimp,rbcl,boprimp,dximp,nsimp,rsimp,nimp,
     > nposimp,
     > taucl,!gtaucl,
c    > ptminva,
     > dosimp(1,1,ie),doshimp(1,1,ie),dosehimp(1,1,ie),
     > dosheimp(1,1,ie),dostehimp(1,1,ie),dostheimp(1,1,ie),
     > tripletmat(1,1,1,1,ie),
     > qvpimp,qvhpimp,
     > qvimp,qvhimp,qvehimp,qvheimp,qvtehimp,qvtheimp,
     > qmomimpha,qmomimpeha,qmomimphea,qmomimpteha,qmomimpthea,! 17/11/2020 Laszloffy added here to handle correctly dynamical arrays
!    > qvdiffimp(1,ie),nqva,  ! TODO : commented on 27/8/2020 Laszloffy
     > vmadid,vmadiq,
     > enbimp,enbhimp,enbdiffimp(1,ie),enbhdiffimp(1,ie),
!    > denba,!nenba, ! commented on 07/07/2021 Laszloffy
!     > enorba,
     > qmomimpa,
     > rhovimp,rhospimp,rhodspimp,rhomagimp,
     > dosmagimp(1,1,ie),dosmaghimp(1,1,ie),
     > spin_magvpimp,spin_magvhpimp,spin_magvimp,spin_magvhimp,
     > orb_magvpimp,orb_magvhpimp,orb_magvimp,orb_magvhimp,
     > lliter,enbifc,enbhifc,qvifc,qvhifc,omifc,iwtau,nwtau,
     > reg,c_light,E_Fermi,tripletout)
c     write(6,*) 'after locquant2'    
cc
c
      tautest=.false.
      if(tautest) then
      write(6,*) ' <skkr> : dosimp locquant'
      do li=1,nimp
        write(6,*) 'li ie',li,ie
        write(6,*) dosimp(:,li,ie)
      end do
      end if
c
c
      end if
 11   continue
c
c            **************************************
c            *** End of energy loop for cluster ***
c            **************************************
cc
#ifdef MPIP
cc  add and broadcast data to all nodes
cc =====================================================================
       write(6,*) 'End of energy loop for cluster, collecting data'
       call cpu_time(tarrive)
       print *,'<skkr> node', myrank,'arrived at',tarrive
       call mpi_barrier(mpi_comm_world,ierr)
       print *,'<skkr> after barrier, before collect_data'
       call collect_data(kmymax,nimp,ne,
     > dosimp,doshimp,dosehimp,
     > dosheimp,dostehimp,dostheimp,
     > qvpimp,qvhpimp,qvimp,qvhimp,qvehimp,qvheimp,
     > qvtehimp,qvtheimp,tripletmat)
c      > nqva,
c      > vmadid,vmadiq,enba,enbdiffa,denba,nenba,enorba,qmoma,
c      > rhova,rhospa,rhodspa,rhomaga,
c      > dosmaga,spin_magvpa,spin_magva,orb_magvpa,orb_magva,
c      > vij,gradph,gradth,
c      > ebl,detl2)
#endif
cc TODO
c
      tautest=.false.
      if(tautest) then
      write(6,*) ' <skkr> : dosimp collect'
      do ie=1,ne
      do li=1,nimp
        write(6,*) 'li ie',li,ie
        write(6,*) dosimp(:,li,ie)
      end do
      end do
      end if
c
cc =====================================================================
       if(myrank.eq.root) then
       call printscf1(
     > nimp,1,
     > qvimp,qvhimp,qvehimp,qvheimp,qvtehimp,qvtheimp)
       end if
c =====================================================================
c Delta update for cluster       
       do li=1,nimp
           deltaold = (0.0d0,0.0d0)
           deltaold(:)=deltaimp(:,li)
           do j=1,nrad
             if (rlambdaimp(li).gt.1.0d-8) then 
               deltaimp(j,li) = 2*qvehimp(li)*rlambdaimp(li)*cmix
     >                  +(1-cmix)*deltaold(j)
             else 
               deltaimp(j,li) = 0.0d0
             end if
           end do
       end do 
      if (myrank.eq.root) then 
          write(6,*) 
          do li=1,nimp
              deltaprint=deltaimp(1,li)*13605.662285137
         write(6,'('' I'',i4,''  iimp'',i4,'' Delta'',e18.10,'' meV'')')
     >               itscfcont,li,deltaprint 
          end do
          write(6,*) 
          do li=1,nimp
          deltaprint=qvehimp(li)*2
         write(6,'('' I'',i4,''  iimp'',i4,'' chi'',e18.10)')
     >               itscfcont,li,deltaprint 
          end do
          write(6,*) 
          do li=1,nimp
          deltaprint=rlambdaimp(li)*13605.662285137
        write(6,'('' I'',i4,''  iimp'',i4,'' lambda'',e18.10,'' meV'')')
     >               itscfcont,li,deltaprint 
          end do
          write(6,*) 
      end if
      call flush(6)
      if (myrank.eq.root) then
        write(6,*)
        do iie=1,ne
          write(6,'('' ie='',i3,'' e='',2f12.8,'' w='',2f12.8,   
     >    ''  kset='',i5)') iie,cear(iie),we(iie),kset(iie)
        end do
        write(6,*)
      end if
      call flush(6)
cc =====================================================================
cc  print results
cc =====================================================================
      if (myrank.eq.root) then
       call printres1(
     >     dosoutimp,pdosoutimp,itscf,lmax,nimp,ne,cear,efermi,
     >     dosimp,doshimp,dosehimp,dosheimp,dostehimp,dostheimp,
     >     qvpimp,qvhpimp,
     >     qvimp,qvhimp,qvehimp,qvheimp,qvtehimp,qvtheimp,
     >     dos,lms)
c================
      call tripletprint(
     >     tripletmat,nimp,ne,cear,tripletout)
      end if
c
      tautest=.false.
      if(tautest) then
      write(6,*) ' <skkr> : dosimp printres'
      do ie=1,ne
      do li=1,nimp
        write(6,*) 'li ie',li,ie
        write(6,*) dosimp(:,li,ie)
      end do
      end do
      end if
      end if ! End of cluster calculations
cc
c
#ifdef MPIP
      call cpu_time(start)
      call mpi_barrier(mpi_comm_world,ierror)
      call cpu_time(finish)
      call flush(6)
      if (ierror /= MPI_SUCCESS) then
          write(6,'("Error in MPI_barrier! (err=",i2)') ierror
          call flush(6)
          call mpi_finalize(ierror)
          stop
      end if
      write(6,'(/" MPI sync idle time =",f10.3," s")') finish-start
      if(itscf.lt.itscfmax) goto 101
#endif
  200 close(6)  
#ifdef MPIP
      call mpi_barrier(mpi_comm_world,ierror)
      call mpi_finalize(ierror)
#endif
      end program main 
