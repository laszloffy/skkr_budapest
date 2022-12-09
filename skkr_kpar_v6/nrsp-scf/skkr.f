c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      program main
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      include 'mpif.h'
c MPI parameters (esimon)
      integer root, myrank, nprocs, ierror
      common/mpi/root,myrank,nprocs
c      
      logical rel,bulk,dos,bulkgeo
      logical cpatest,cpamatin,cpamatinl,cpamatinr,cpain
      logical dlm,dlml,dlmr
      logical nspheres
c
      character*10 idpota(mintfc),idpotla(minprc),idpotra(minprc)
      character*10 idpotb(mintfc),idpotlb(minprc),idpotrb(minprc)
      character*30 printout,potout,dosout,pdosout,tmatout,momout
      character*30 leftpot,rightpot,laypot,laycore
      character*30 leftmat,rightmat,laymat
      character*30 leftmom,rightmom
      character*1 rightm
c
      dimension park(2),kset(me)
      dimension xk(mkpar,2),wk(mkpar)
      dimension nepanel(5),ebottom(5),etop(5),eps(5)
c
      dimension iesublatt(mintfc),lsublatt(mintfc,melem)
      dimension vra(nrad,mintfc),vrb(nrad,mintfc)
      dimension bra(nrad,mintfc),brb(nrad,mintfc)
      dimension wra(nrad,mintfc),wrb(nrad,mintfc)
      dimension dx(mintfc),ns(mintfc),rs(mintfc),rmt(mintfc)
      dimension vrla(nrad,minprc),vrlb(nrad,minprc)
      dimension brla(nrad,minprc),brlb(nrad,minprc)
      dimension wrla(nrad,minprc),wrlb(nrad,minprc)
      dimension rsl(minprc),dxl(minprc),nsl(minprc)
      dimension vrra(nrad,minprc),vrrb(nrad,minprc)
      dimension brra(nrad,minprc),brrb(nrad,minprc)
      dimension wrra(nrad,minprc),wrrb(nrad,minprc)
      dimension rsr(minprc),dxr(minprc),nsr(minprc)
c
      dimension za(mintfc),zb(mintfc),zla(minprc),zlb(minprc),
     &          zra(minprc),zrb(minprc) 
      dimension conc(mintfc),sxa(mintfc),sxb(mintfc)
      dimension concl(minprc),sxl(minprc)
      dimension concr(minprc),sxr(minprc) 
      dimension jspl(minprc),jspr(minprc),jspa(mintfc),jspb(mintfc)
      dimension qca(mintfc),qcb(mintfc)
      dimension qva(mintfc),qvb(mintfc)
      dimension qmagva(mintfc),qmagvb(mintfc)
      dimension qwa(mintfc,2),qwb(mintfc,2)
      dimension qvpa(0:lmaxp,mintfc),qvpb(0:lmaxp,mintfc)
      dimension qmagvpa(0:lmaxp,mintfc),qmagvpb(0:lmaxp,mintfc)
      dimension qwpa(0:lmaxp,mintfc,2),qwpb(0:lmaxp,mintfc,2)
c
      dimension dosa(0:lmaxp,mintfc,me,2),dosb(0:lmaxp,mintfc,me,2)
      dimension doslma(lmmaxp,mintfc,me,2),doslmb(lmmaxp,mintfc,me,2)
      dimension dosla(0:lmaxp,mintfc,me),doslb(0:lmaxp,mintfc,me)
      dimension spdosw(mkpar,0:lmaxp,mintfc,2)
      dimension sdosw(0:lmaxp,mintfc,2)
      dimension rhoca(nrad,mintfc),rhocb(nrad,mintfc)
      dimension rhova(nrad,mintfc),rhovb(nrad,mintfc)
      dimension rhomagva(nrad,mintfc),rhomagvb(nrad,mintfc)
      dimension rhowa(nrad,mintfc,2),rhowb(nrad,mintfc,2)
      dimension enbwa(mintfc,2),enbwb(mintfc,2)
      dimension enba(mintfc),enbb(mintfc)
      dimension enca(mintfc),encb(mintfc)
      dimension enela(mintfc),enelb(mintfc)
      dimension enxca(mintfc),enxcb(mintfc)
      dimension enpota(mintfc),enpotb(mintfc)
      dimension enkina(mintfc),enkinb(mintfc)
      dimension enmaga(mintfc),enmagb(mintfc)
c
      complex*16 cear(me),we(me),ce
      complex*16 tminvl(lmmaxp,lmmaxp,minprc)
      complex*16 tminvr(lmmaxp,lmmaxp,minprc)
      complex*16 tminv(lmmaxp,lmmaxp,mintfc)
      complex*16 tminva(lmmaxp,lmmaxp,mintfc)
      complex*16 tminvb(lmmaxp,lmmaxp,mintfc)
      complex*16 taua(lmmaxp,lmmaxp,mintfc)
      complex*16 taub(lmmaxp,lmmaxp,mintfc)
      complex*16 taukdiag(mkpar,0:lmaxp,mintfc)
      complex*16 qmoma(lmsup,mintfc),qmomb(lmsup,mintfc)
      complex*16 qmomwa(lmsup,mintfc,2),qmomwb(lmsup,mintfc,2)
      complex*16 qmomla(lmsup,minprc),qmomlb(lmsup,minprc)
      complex*16 qmomra(lmsup,minprc),qmomrb(lmsup,minprc)
c
      common/test/itest
      common/lay2d/cvec(mtotal,3),nextra,nbulkl,nbulkr,
     &             nprc,ninprc(0:mprc+1)
      common/pot/potshift,b0,layshift1,layshift2,layb0,layb0f,ib0
c
      data tiny/1.0d-6/
c     data tolerr,tolef,tolv0/1.0d-11,1.0d-08,1.0d-08/
c
c
c     ***************************************************************
c                            INITIALIZE
c     ***************************************************************
c
!   MPI BEGIN
      root = 0
      call mpi_init(ierror)
      call mpi_comm_rank(mpi_comm_world,myrank,ierror)
      call mpi_comm_size(mpi_comm_world,nprocs,ierror)
!   MPI END
      open(unit=5,file='input_nrsp.in',status='old')
c     ------------------------------------------------------------------
      call readini(imesh,ne,npanel,nepanel,ne1,ne2,ne3,ebottom,etop,
     >eps,lmax,eta,sigma,park,kset,itcpam,cpatol,cpatest,cpain,
     >printout,potout,dosout,pdosout,tmatout,momout,
     >leftpot,rightpot,laypot,laycore,
     >leftmat,rightmat,laymat,leftmom,rightmom,
     >bulk,nspheres,rel,dos,rightm,newvac,ivacpot,
     >v0,v00,vrsh,iscreen,vscreen,itscfmax,tolerr,tolef,intbz,kunning,
     >dlm,dlml,dlmr)
c     ------------------------------------------------------------------
      close(5)
c
cMPI
      write(6,'(/" Running on",i3," parallel processes."/)') nprocs
c
      open(unit=5,file='input_geo.in',status='old')
c     ------------------------------------------------------------------
      call readgeom(bulk,bulkgeo,rightm,intbz,nintfc,arot,
     > rsl,concl,jspl,rs,conc,jspa,jspb,rsr,concr,jspr,
     > igraph,iesublatt,lsublatt,volint)
c     ------------------------------------------------------------------
      close(5)
c
      write(6,'(/2x,''routine MAIN>''/)')
c
      nl=lmax+1
      nl2=nl*nl
      l2=2*lmax
      lmmaxs=(l2+1)*(l2+1)
      ninprcl=ninprc(0)
      ninprcr=ninprc(nprc+1)
c
c -read Left,Right and Layer potentials
c     ------------------------------------------------------------------
      call pothandle(
     > bulk,rightm,lmax,nintfc,ninprcl,ninprcr,vrsbulk,v00,ivacpot,
     > leftpot,leftmom,concl,qmomla,qmomlb,dxl,nsl,rsl,jspl,
     > idpotla,vrla,brla,zla,idpotlb,vrlb,brlb,zlb,
     > rightpot,rightmom,concr,qmomra,qmomrb,dxr,nsr,rsr,jspr,
     > idpotra,vrra,brra,zra,idpotrb,vrrb,brrb,zrb,
     > laypot,conc,dx,ns,rs,idpota,vra,bra,za,jspa,
     > idpotb,vrb,brb,zb,jspb,dlm,dlml,dlmr,igraph)
c     ------------------------------------------------------------------
c
c -read in data from update file if exists
c     ------------------------------------------------------------------
      call restart_in(itscf0,etop(npanel),v0,conc,vra,vrb,bra,brb,
     >     bulk,vrsbulk,nintfc,idpota,idpotb,za,zb,laycore,printout)
c     ------------------------------------------------------------------
c
c Neutral spheres block
      if(nspheres) then
        write(6,'(/'' WS radii:'')')
        do li=1,nintfc
          write(6,'('' I'', i3,''  Sws'',f17.8)') li,rs(li)
        end do
      end if
c
c
c
c--------------------------------------------------------------------------
      write(6,*)
      write(6,*) 'etop=',etop(npanel)
      if(.not.bulk) write(6,*) 'v0=  ',v0
c
c -fix manipulations with files containing inverse cpa t-matrices
c     ------------------------------------------------------------------
      call tcpa_main(
     > leftmat,rightmat,laymat,bulk,dos,kset(1),nintfc,ninprcl,ninprcr,
     > conc,concl,concr,cpain,cpamatin,cpamatinl,cpamatinr)
c     ------------------------------------------------------------------
c
c -initialize common block 'mome' and spherical harmonics
c     ------------------------------------------------------------------
      call initia
c     ------------------------------------------------------------------
c
c
c    **************************************************************
c                   START SELFCONSISTENT ITERATIONS 
c    **************************************************************
c
      itscf=0
      istop=0
      efermi=etop(npanel)
  100 continue
      itscf=itscf+1
      itscfcont=itscf0+itscf
c
      write(6,'(/'' Selfconsistent iteration: '',i4/)') itscf
      call flush(6)
c
      call sublattpot(nintfc,ns,vra,iesublatt)
      call sublattpot(nintfc,ns,vrb,iesublatt)
      call sublattpot(nintfc,ns,bra,iesublatt)
      call sublattpot(nintfc,ns,brb,iesublatt)
      if(bulk) then
c     ------------------------------------------------------------------
         call averagebpot(nsl,nintfc,nbulkl,vra,vrla,vrra)
         call averagebpot(nsl,nintfc,nbulkl,vrb,vrlb,vrrb)
         call averagebpot(nsl,nintfc,nbulkl,bra,brla,brra)
         call averagebpot(nsl,nintfc,nbulkl,brb,brlb,brrb)
c     ------------------------------------------------------------------
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
c -generate energy mesh if required
c
      if(itscf.eq.1.or.bulk) then 
c     ------------------------------------------------------------------
         call zmesh(imesh,npanel,nepanel,ne1,ne2,ne3,ebottom,etop,
     &              eps,cear,we)
c     ------------------------------------------------------------------
c        do ie=1,ne
c          ce=cear(ie)
c          write(6,'('' e='',2f12.6,5x,2f12.6)') ce,we(ie)
c        end do
c        call flush(6)
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
      call initzero(qwpa(0,1,1),qwpb(0,1,1),qwa(1,1),qwb(1,1),
     >              enbwa(1,1),enbwb(1,1),qmomwa(1,1,1),qmomwb(1,1,1),
     >              rhowa(1,1,1),rhowb(1,1,1))
      call initzero(qwpa(0,1,2),qwpb(0,1,2),qwa(1,2),qwb(1,2),
     >              enbwa(1,2),enbwb(1,2),qmomwa(1,1,2),qmomwb(1,1,2),
     >              rhowa(1,1,2),rhowb(1,1,2))
c     ------------------------------------------------------------------
c
       if(dos.and.ne.eq.1) open(9,file=dosout,status='unknown')
c
c                 ******************************
c                 *** Start spin chanel loop ***
c                 ******************************
c
      do nsp=1,2
c
c -set up 'spin-down' or 'spin-up' potentials
c
      if(nsp.eq.1) isp=-1
      if(nsp.eq.2) isp= 1
      do li=1,ninprcl
        do i=1,nsl(li)
          wrla(i,li)=vrla(i,li)+isp*brla(i,li)
          wrlb(i,li)=vrlb(i,li)+isp*brlb(i,li)
        end do
      end do
      do li=1,nintfc
        do i=1,ns(li)
          wra(i,li)=vra(i,li)+isp*bra(i,li)
          wrb(i,li)=vrb(i,li)+isp*brb(i,li)
        end do
      end do
      do li=1,ninprcr
        do i=1,nsr(li)
          wrra(i,li)=vrra(i,li)+isp*brra(i,li)
          wrrb(i,li)=vrrb(i,li)+isp*brrb(i,li)
        end do
      end do
      if(itest.ge.1) write(6,'(/'' Spin channel '',i2)') nsp
c
c                   *************************
c                   *** Start energy loop ***
c                   *************************
c
      iek=0
      do 10 ie=1,ne
c
      ce=cear(ie)
c
c -read in (inverse of) effective t-matrices if necessary
c     ------------------------------------------------------------------
      call tcpa_in(leftmat,rightmat,laymat,ce,nintfc,
     > ninprcl,ninprcr,conc,concl,concr,cpamatin,cpamatinl,
     > cpamatinr,nl2,tminv,tminvl,tminvr,lmmaxp)
c     ------------------------------------------------------------------
c
c  generate k-mesh in irreducible Brillouin-zone
c
      if(kset(ie).eq.0) itcpam=0
c     ------------------------------------------------------------------
      call kmesh2d(intbz,kunning,kset(ie),park,xk(1,1),xk(1,2),wk,nk,
     &             mkpar)
c     ------------------------------------------------------------------
c
c     if(dos.and.ne.eq.1) then
c       ik=0
c       krow=2*kset(1)+1
c       xxk=dabs(park(1))
c       yyk=dabs(park(2))
c       dxk=xxk/kset(1)
c       dyk=yyk/kset(1)
c       do i1=1,krow
c       do i2=1,krow
c         ik=ik+1
c         xk(ik,1)=-xxk+dxk*(i1-1)
c         xk(ik,2)=-yyk+dyk*(i2-1)
c         wk(ik)=1.0d0
c       enddo
c       enddo
c       nk=ik
c     endif 
c
      write(6,'('' ie='',i3,''   e='',2f12.8,''  nk='',i5)') ie,ce,nk
      call flush(6)
      if(itest.ge.2.or.kset(ie).eq.0) then
        write(6,'('' no. of k-points in IBZ:'',i5)') nk
        do ik=1,nk
           write(6,'(i5,3f14.10)') ik,xk(ik,1),xk(ik,2),wk(ik)
        end do
      endif
      call flush(6)
c 
 
c  t-matrices to start CPA
c     ------------------------------------------------------------------
      call tmatini(
     > nintfc,ce,lmax,rel,v0,bulk,iscreen,vscreen,
     > concl,concr,conc,cpamatin,cpamatinl,cpamatinr,
     > idpotla,wrla,idpotlb,wrlb,dxl,nsl,rsl,
     > idpotra,wrra,idpotrb,wrrb,dxr,nsr,rsr,
     > idpota,wra,idpotb,wrb,dx,ns,rs,
     > tminvl,tminvr,tminv,tminva,tminvb)
c     ------------------------------------------------------------------
c
c  perform CPA, BZ integration & calculate tau-matices   
c     ------------------------------------------------------------------
      call cpacoord(
     > nsp,itscf,ie,ce,we(ie),lmax,nintfc,eta,rightm,
     > iesublatt,bulk,'false',bulkgeo,kset(ie),xk,wk,nk,intbz,iek,
     > conc,itcpam,cpatol,cpatest,
     > tminvl,tminvr,tminv,tminva,tminvb,taua,taub,taukdiag)
c     ------------------------------------------------------------------
c
c -write out (inverse of) effective t-matrices
c     ------------------------------------------------------------------
       if (myrank.eq.root) then
      call tcpa_out(ce,nintfc,conc,nl2,tminv,lmmaxp)
        end if
c     ------------------------------------------------------------------
c
c  calculate local physical quantities
c  integrate with respect to energy   
c     ------------------------------------------------------------------
       call locquant(
     > .false.,ce,we(ie),lmax,nintfc,rel,conc,v0,
     > idpota,wra,idpotb,wrb,dx,ns,rs,taua,taub,taukdiag,nk,wk,
     > dosa(0,1,ie,nsp),doslma(1,1,ie,nsp),qwpa(0,1,nsp),qwa(1,nsp),
     > enbwa(1,nsp),qmomwa(1,1,nsp),rhowa(1,1,nsp),
     > dosb(0,1,ie,nsp),doslmb(1,1,ie,nsp),qwpb(0,1,nsp),qwb(1,nsp),
     > enbwb(1,nsp),qmomwb(1,1,nsp),rhowb(1,1,nsp),
     > spdosw(1,0,1,nsp),sdosw(0,1,nsp))
c     ------------------------------------------------------------------
c
 10    continue
c
c                 ************************** 
c                 *** End of energy loop ***
c                 ************************** 
c
c      if(dos.and.ne.eq.1) then
c        write(9,'(''# energy '',2f13.5)') ce 
c        write(9,'(''# Spin '',i3)') nsp
c        do ik=1,nk
c          write(9,'(50f13.6)') xk(ik,1),xk(ik,2),
c    >     ((spdosw(ik,l,li,nsp),l=0,lmax),li=layshift1,layshift2)
c        end do
c        ik=0
c        do i1=1,krow
c        do i2=1,krow
c          ik=ik+1
c          write(9,'(50f13.6)') xk(ik,1),xk(ik,2),
c    >     ((spdosw(ik,l,li,nsp),l=0,lmax),li=layshift1,layshift2)
c        end do
c        write(9,*)
c        end do
c        write(9,'(50f13.6)') 
c    >   ((sdosw(l,li,nsp),l=0,lmax),li=layshift1,layshift2)
c        write(9,'(50f13.6)') 
c    >   ((dosa(l,li,1,nsp),l=0,lmax),li=layshift1,layshift2)
c      end if
c
      end do
c
c               *******************************
c               *** End of spin chanel loop ***
c               *******************************
c
c     if(dos.and.ne.eq.1) then
c        close(9)
c        stop
c     end if
c
c
c -recalculate quantities
c
      do li=1,nintfc
        qva(li)=qwa(li,1)+qwa(li,2)
        qvb(li)=qwb(li,1)+qwb(li,2)
        qmagva(li)=qwa(li,2)-qwa(li,1)
        qmagvb(li)=qwb(li,2)-qwb(li,1)
        do l=0,lmax
          qvpa(l,li)=qwpa(l,li,1)+qwpa(l,li,2)
          qvpb(l,li)=qwpb(l,li,1)+qwpb(l,li,2)
          qmagvpa(l,li)=qwpa(l,li,2)-qwpa(l,li,1)
          qmagvpb(l,li)=qwpb(l,li,2)-qwpb(l,li,1)
          dosla(l,li,ne)=dosa(l,li,ne,1)+dosa(l,li,ne,2)
          doslb(l,li,ne)=dosb(l,li,ne,1)+dosb(l,li,ne,2)
        end do
        do i=1,lmmaxs
          qmoma(i,li)=qmomwa(i,li,1)+qmomwa(i,li,2)
          qmomb(i,li)=qmomwb(i,li,1)+qmomwb(i,li,2)
        end do
        enba(li)=enbwa(li,1)+enbwa(li,2)
        enbb(li)=enbwb(li,1)+enbwb(li,2)
        do i=1,ns(li)
          rhova(i,li)=rhowa(i,li,1)+rhowa(i,li,2)
          rhovb(i,li)=rhowb(i,li,1)+rhowb(i,li,2)
          rhomagva(i,li)=rhowa(i,li,2)-rhowa(i,li,1)
          rhomagvb(i,li)=rhowb(i,li,2)-rhowb(i,li,1)
        end do
      end do
c
c -close binary files containing (inverse of) effective t-matrices
c     ------------------------------------------------------------------
      call tcpa_close(
     >leftmat,rightmat,laymat,cpamatin,cpamatinl,cpamatinr)
c     ------------------------------------------------------------------
c
c for DOS calculation print out results and stop
      if(dos) goto 99
c
c write out results
c     ------------------------------------------------------------------
      call printscf(nintfc,itscfcont,conc,qva,qvb,qmagva,qmagvb)
c     ------------------------------------------------------------------
c
c -find new Fermi level      
c
      if(bulk) then 
        nl0=0
c     ------------------------------------------------------------------
        call newfl(itscfcont,nl0,lmax,lmaxp,ne,nbulkl,conc,dosla,doslb,
     &             qca,qcb,qva,qvb,za,zb,efermi,defermi)
c     ------------------------------------------------------------------
      end if
c
      if (myrank.eq.root) then
c
c -Generate new layer potentials 
c     ------------------------------------------------------------------
      call vgen(
     > itscf,itscfcont,lmax,nintfc,sigma,iesublatt,
     > rightm,bulk,vrsbulk,vrsh,nspheres,volint,
     > conc,concl,concr,rsl,rsr,
     > za,qca,qva,qmoma,qmomla,qmomra,
     > zb,qcb,qvb,qmomb,qmomlb,qmomrb,
     > rhoca,rhova,rhocb,rhovb,rhomagva,rhomagvb,dx,ns,rs,
     > efermi,defermi,vra,vrb,bra,brb,v0,dv0,ferr1,ferr2,newvac,
     > enpota,enela,enxca,enpotb,enelb,enxcb,enmaga,enmagb,dlm)
c     ------------------------------------------------------------------
c
c -calculate and print out energies
c     ------------------------------------------------------------------
      call enprint(itscfcont,nintfc,conc,enca,encb,enba,enbb,
     >enela,enelb,enxca,enxcb,enpota,enpotb,enkina,enkinb,
     >enmaga,enmagb)
c     ------------------------------------------------------------------
c
c -update output necessary for restart
c     ------------------------------------------------------------------
      call restart_out(itscfcont,efermi,v0,conc,vra,vrb,bra,brb,
     >                 bulk,vrsbulk,nintfc,printout)
c     ------------------------------------------------------------------
c
      end if
c
      call vgen_bcast(efermi,defermi,vra,vrb,bra,brb,rs,rsl,rsr,
     >  v0,dv0,ferr1,ferr2,entotifc)
c
      if(bulk) etop(npanel)=efermi
c
      if(ferr1.lt.tolerr) istop=1
      if(bulk.and.dabs(defermi).gt.tolef) istop=0
c     if(rightm.eq.'V'.and.dv0.gt.tolv0) istop=0
c
c -Print results
c
   99 continue
c     ------------------------------------------------------------------
      if (myrank.eq.root) then
      call printres(
     >potout,dosout,pdosout,momout,
     >itscfcont,lmax,nintfc,ne,cear,efermi,conc,dlm,bulk,
     >dosa,dosb,doslma,doslmb,qvpa,qvpb,qva,qvb,
     >qmagvpa,qmagvpb,qmagva,qmagvb,qmoma,qmomb,za,zb,qca,qcb,
     >enba,enbb,enkina,enkinb,enela,enelb,enxca,enxcb,
     >idpota,idpotb,vra,vrb,bra,brb,dx,ns,rs,vrsbulk)
      end if
c     ------------------------------------------------------------------
      call cpu_time(start)
      call mpi_barrier(mpi_comm_world,ierror)
      call cpu_time(finish)
      call flush(6)
      if (ierror /= MPI_SUCCESS) then
          write(6,'("Error in MPI_barrier! (err=",i2)') ierror
          call flush(6)
          stop
      end if
      write(6,'(/" MPI sync idle time =",f10.3," s")') finish-start
c
c
      if(itscf.lt.itscfmax.and.istop.eq.0) goto 100
      if(istop.eq.1) write(6,'(/" SCF cycle converged!")')
c
c     ********************************************************
c               END OF SELFCONSISTENT ITERATIONS
c     ********************************************************
c
      call mpi_finalize(ierror)

  200 if (myrank.eq.root) close(6)
      stop
      end
