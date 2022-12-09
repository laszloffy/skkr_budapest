c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      program main
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
      include 'mpif.h'
c
      logical bulk,dos,wrel,lms,bulkgeo,linbw,calcder
      logical orbpol,opotl,opotr,opot
      logical cpatest,cpamatin,cpamatinl,cpamatinr,cpain
      logical singlesite,reg
      logical bsfcalc
c
      character*10 idpota(mintfc),idpotla(minprc),idpotra(minprc)
      character*10 idpotb(mintfc),idpotlb(minprc),idpotrb(minprc)
      character*30 printout,potout,dosout,pdosout,tmatout,momout
      character*30 leftpot,rightpot,laypot,laycore
      character*30 leftdelta,rightdelta,laydelta
      character*30 leftmat,rightmat,laymat
      character*30 leftmom,rightmom
      character*1 rightm
c
      real*8 lza(2,mintfc),lzb(2,mintfc)
      real*8 tstart,tfinish
      real*8 entotifc
c
      dimension park(2),kset(me)
      dimension xk(mkpar,2),wk(mkpar)
      dimension nepanel(5),ebottom(5),etop(5),eps(5)
      dimension iesublatt(mintfc)
c
      complex*16 deltala(nrad,minprc),deltara(nrad,minprc)
      complex*16 deltaa(nrad,mintfc)
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
      complex*16 cear(me),we(me),ce
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
c
      common/test/itest
      common/lay2d/cvec(mtotal,3),nextra,nbulkl,nbulkr,
     &             nprc,ninprc(0:mprc+1) 
c
      data tiny/1.0d-6/
      data rnul/0.d0/
      data sqrtm1/(0.d0,1.d0)/
c     data tolerr,tolef,tolen,tolv0/1.0d-10,1.0d-10,1.0d-9,1.0d-08/
c MPI
      integer root, myrank, nprocs, ierror
      common/mpi/root,myrank,nprocs
c
c     ***************************************************************
c                            INITIALIZE
c     ***************************************************************
c
c MPI begin
      root = 0
      call mpi_init(ierror)
      call mpi_comm_rank(mpi_comm_world,myrank,ierror)
      call mpi_comm_size(mpi_comm_world,nprocs,ierror)
c   MPI END
c
      open(unit=5,file='input_rsp.in',status='old')
c     ------------------------------------------------------------------
      call readini(
     > imesh,ne,npanel,nepanel,ne1,ne2,ne3,ebottom,etop,eps,lmax,eta,
     > sigma,park,kset,ksetmax,itcpam,cpatol,cpatest,cpain,
     > printout,potout,dosout,pdosout,tmatout,momout,
     > leftpot,rightpot,laypot,laycore,
     > leftmat,rightmat,laymat,leftmom,rightmom,
     > wrel,lms,bulk,dos,rightm,newvac,ivacpot,v0,v00,vrsh,
     > vscreen,iscreen,itscfmax,tolerr,tolef,tolen,intbz,kunning,
     > orbpol,opotl,opotr,opot,lliter,E_Fermi,
     > leftdelta,rightdelta,laydelta,singlesite,reg,c_light)
c     ------------------------------------------------------------------
      close(5)
c
      write(6,'(/" Running on",i3," parallel processes."/)') nprocs
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
c
      write(6,'(/2x,''routine MAIN>''/)')
c
      nl=lmax+1
      lmmax=nl*nl
      kmymax=2*lmmax
      ntotal=nintfc+ninprc(0)*(nextra+1)+ninprc(nprc+1)*(nextra+1)
c
c -read Left,Right and Layer potentials
c
c     ------------------------------------------------------------------
      call pothandle(
     > bulk,linbw,rightm,lmax,nintfc,ninprc(0),ninprc(nprc+1),
     > vrsbulk,v00,ivacpot,opotl,opotr,opot,
     > leftpot,leftmom,concl,qmomla,qmomlb,dxl,nsl,rsl,
     > idpotla,vrla,brla,boprla,zla,idpotlb,vrlb,brlb,boprlb,zlb,
     > rightpot,rightmom,concr,qmomra,qmomrb,dxr,nsr,rsr,
     > idpotra,vrra,brra,boprra,zra,idpotrb,vrrb,boprrb,brrb,zrb,
     > laypot,conc,dx,ns,rs,
     > idpota,vra,bra,bopra,za,idpotb,vrb,brb,boprb,zb,igraph,
     > leftdelta,rightdelta,laydelta,
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
      write(6,*)
      write(6,'(2x,''etop= '',f17.13)') etop(npanel)
      if(.not.bulk) write(6,*) 'v0=  ',v0
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
c     ------------------------------------------------------------------
      call initia
c     ------------------------------------------------------------------
c
c     *************************************************************
c               START SELFCONSISTENT ITERATIONS 
c     *************************************************************
c
      itscf=0
      istop=0
      efermi=etop(npanel)
  100 continue
      itscf=itscf+1
      itscfcont=itscf0+itscf
      write(6,'(/'' Selfconsistent iteration: '',i5/)') itscf
      efermi0=efermi
c
c     ------------------------------------------------------------------
      call sublattpot(nintfc,ns,vra,iesublatt)
      call sublattpot(nintfc,ns,vrb,iesublatt)
      call sublattpot(nintfc,ns,bra,iesublatt)
      call sublattpot(nintfc,ns,brb,iesublatt)
      call sublattpot(nintfc,ns,bopra,iesublatt)
      call sublattpot(nintfc,ns,boprb,iesublatt)
      call sublattpot(nintfc,ns,deltaa,iesublatt)
      call sublattpot(nintfc,ns,deltab,iesublatt)
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
c
      if(itscf.eq.1.or.bulk) then
c     ------------------------------------------------------------------
         call zmesh(imesh,npanel,nepanel,ne1,ne2,ne3,ebottom,etop,
     &              eps,cear,we)
c     ------------------------------------------------------------------
      end if
      if(itest.ge.2.and.itscf.eq.1) then 
        write(6,*)
        do ie=1,ne
          write(6,'('' ie='',i2,'' e='',2f12.8,'' w='',2f12.8,   
     >    ''  kset='',i5)') ie,cear(ie),we(ie),kset(ie)
        end do
        write(6,*)
      end if
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
c
c                 *************************
c                 *** Start energy loop ***
c                 *************************
c
      ie0=1
      iek=0
      do 10 ie=1,ne
c
      ce=cear(ie)
c
c -read in (inverse of) effective t-matrices if necessary
c     ------------------------------------------------------------------
      call tcpa_in(leftmat,rightmat,laymat,ce,nintfc,
     > ninprc(0),ninprc(nprc+1),conc,concl,concr,
     > cpamatin,cpamatinl,cpamatinr,2*kmymax,
     > tminv,tminvl,tminvr,
     > tau_kint,taua_kint,taub_kint,dbogomaxp)
c      if (cpamatin) then
c      write(6,*) 'tminv from tcpa'
c      call outmat1(tminv(:,:,1),2*kmymax,2*kmymax,dbogomaxp,tolerr,6)
c      write(6,*) 'tau from tcpa'
c      call outmat1(tau_kint(:,:,1),2*kmymax,2*kmymax,dbogomaxp,tolerr,6)
c      write(6,*) 'taua from tcpa'
c      call outmat1(taua_kint(:,:,1),2*kmymax,2*kmymax,dbogomaxp,tolerr,6)
c      write(6,*) 'taub from tcpa'
c      call outmat1(taub_kint(:,:,1),2*kmymax,2*kmymax,dbogomaxp,tolerr,6)
c      end if 
c     ------------------------------------------------------------------
c
c  generate k-mesh in irreducible Brillouin-zone
c
      if(kset(ie).eq.0) then
        itcpam=0
        intbz=0
        bsfcalc=.true.
      else 
        bsfcalc=.false.
      end if
c     ------------------------------------------------------------------
      call kmesh2d(intbz,kunning,kset(ie),park,xk(1,1),xk(1,2),wk,nk,
     &             mkpar)
c     ------------------------------------------------------------------ 
c
      if(itest.ge.1) then
        write(6,'('' ie='',i3,''   e='',2f12.8,''  nk='',i5)') ie,ce,nk
        call flush(6)
      end if
      if(itest.ge.2) then
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
c      write(6,*) 'tminva from tmatini for li=1'
c      call outmat1(tminva(:,:,1),2*kmymax,2*kmymax,dbogomaxp,tolerr,6)
c      write(6,*) 'tminvb from tmatini for li=1'
c      call outmat1(tminvb(:,:,1),2*kmymax,2*kmymax,dbogomaxp,tolerr,6)
c      write(6,*) 'tminv from tmatini for li=1'
c      call outmat1(tminv(:,:,1),2*kmymax,2*kmymax,dbogomaxp,tolerr,6)
c     ------------------------------------------------------------------ 
      call cpu_time(tfinish)
      write(6,'('' Time ellapsed in tmatini:'',t35,f8.1,'' s'')')
     > tfinish-tstart
c
       if(itest.ge.2) then
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
        if (cpamatin) then 
          write(6,*) 'tau calculation for a single k point'
          write(6,*) xk(nk,1:2)
        end if
        call cpacoord(
     >   itscf,ie,ce,E_Fermi,lmax,nintfc,eta,rightm,bulk,bulkgeo,wrel,
     >   kset(ie),xk,wk,nk,intbz,iek,
     >   conc,itcpam,cpatol,cpatest,
     >   dmata,dmatb,dmatpa,dmatpb,
     >   tminvl,tminvr,tminv,tminva,tminvb,
     >   tau,taua,taub,gtaua,gtaub,c_light)
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
      write(6,'('' Time ellapsed in cpacoord:'',t35,f8.1,'' s'')')
     > tfinish-tstart
c
       if(itest.ge.2) then
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
       end if
c
c
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
c
c
c      stop 'The code is stopped before locquant'
c
c
c
c  calculate local physical quantities
c  integrate with respect to energy
c     ------------------------------------------------------------------
      if (.not.bsfcalc) then
      call cpu_time(tstart)
      call locquant(
     > ie,ce,we(ie),dos,lmax,nintfc,wrel,lms,sxa,sxb,conc,v0,
     > E_Fermi,deltaa,deltab,
     > idpota,vra,bra,rba,idsotb,vrb,brb,rbb,dx,ns,rs,
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
      call cpu_time(tfinish)
      write(6,'('' Time ellapsed in locquant:'',t35,f8.1,'' s'')')
     > tfinish-tstart
      else 
      call cpu_time(tstart)
      write(6,*) '>- bsf calculation' 
      call bsf(
     > ie,ce,we(ie),dos,lmax,nintfc,wrel,lms,sxa,sxb,conc,v0,
     > E_Fermi,deltaa,deltab,
     > idpota,vra,bra,rba,idsotb,vrb,brb,rbb,dx,ns,rs,
     > tminv,tau,tau_kint,taua_kint,taub_kint,gtaua,gtaub,
c we store the bsf values in the dos variables
     > dosa(1,:,ie),dosha(1,:,ie),doseha(1,:,ie),doshea(1,:,ie),
     > dosteha(1,:,ie),dosthea(1,:,ie),
     > lliter,linbw,efermi,enbifc,qvifc,omifc,reg,c_light)
      call cpu_time(tfinish)
      write(6,'('' Time ellapsed in bsf:'',t35,f8.1,'' s'')')
     > tfinish-tstart
      end if 
c     ------------------------------------------------------------------
c
c      stop 'The code is stopped after locquant'
c                                    
 10   continue
c
c                 **************************
c                 *** End of energy loop ***
c                 **************************

      if (myrank.eq.root) then
      if (bsfcalc) then
      write(6,*) 'The result of the BSF calculation'
      do ie=1,ne
      do li=1,nintfc
        write(6,'('' li='',i3,''   en='',2f12.8,''  k=('',f12.5,'' 
     > '',f12.5'' )    bsf='',d14.6''  bsfh='',d14.6'' bsf s='',d14.6)')
     >       li,cear(ie),xk(nk,1:2),dosa(1,li,ie),
     >       dosha(1,li,ie),doseha(1,li,ie)
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
       end if 
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
      write(6,*) 'before enprint', entotifc
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
c     ------------------------------------------------------------------
      if (myrank.eq.root) then
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
     >     rs,dx,ns,vrsbulk,lms,orbpol,bulk)
      end if
c     ------------------------------------------------------------------
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
c

      if(itscf.lt.itscfmax.and.istop.eq.0) goto 100
      if(istop.eq.1) write(6,'(/" SCF cycle converged!")')
c      if(itscf.lt.itscfmax) goto 100
c
c     ****************************************************************
c                    END OF SELFCONSISTENT ITERATIONS 
c     ****************************************************************
c   
      call mpi_finalize(ierror)
  200 if (myrank.eq.root) close(6)
      stop
      end
