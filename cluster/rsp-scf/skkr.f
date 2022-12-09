c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      program main
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      logical bulk,dos,wrel,lms,bulkgeo,linbw
      logical orbpol,opotl,opotr,opot
      logical cpatest,cpamatin,cpamatinl,cpamatinr
c
      character*10 idpota(mintfc),idpotla(minprc),idpotra(minprc)
      character*10 idpotb(mintfc),idpotlb(minprc),idpotrb(minprc)
      character*30 for006,for007,for008,for009,for010
      character*30 leftpot,rightpot,laypot,laycore
      character*30 leftmat,rightmat,laymat
      character*1 rightm
c
      real*8 lza(2,mintfc),lzb(2,mintfc)
c
      dimension park(2),kset(me)
      dimension xk(mkpar,2),wk(mkpar)
      dimension nepanel(5),ebottom(5),etop(5),eps(5)
      dimension iesubl(mtotal),iesublatt(mintfc),lsublatt(mintfc,melem)
c
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
      dimension qva(mintfc),qvb(mintfc)
      dimension qvpa(kmymaxp,mintfc),qvpb(kmymaxp,mintfc)
c
      dimension rbl(3,minprc),rbr(3,minprc)
      dimension rba(3,mintfc),rbb(3,mintfc)
      dimension vecna(3,mintfc),phia(mintfc)
      dimension vecnb(3,mintfc),phib(mintfc)
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
      complex*16 ddphb(kmymaxp,kmymaxp,mintfc)
      complex*16 ddphpa(kmymaxp,kmymaxp,mintfc)
      complex*16 ddphpb(kmymaxp,kmymaxp,mintfc)
      complex*16 ddtha(kmymaxp,kmymaxp,mintfc)
      complex*16 ddthb(kmymaxp,kmymaxp,mintfc)
      complex*16 ddthpa(kmymaxp,kmymaxp,mintfc)
      complex*16 ddthpb(kmymaxp,kmymaxp,mintfc)
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
      dimension denba(mintfc),denbb(mintfc)
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
      complex*16 cear(me),we(me),ce
      complex*16 tminvl(kmymaxp,kmymaxp,minprc)
      complex*16 tminvr(kmymaxp,kmymaxp,minprc)
      complex*16 tminv(kmymaxp,kmymaxp,mintfc)
      complex*16 tminva(kmymaxp,kmymaxp,mintfc)
      complex*16 tminvb(kmymaxp,kmymaxp,mintfc)
      complex*16 ptminva(kmymaxp,kmymaxp,mintfc)
      complex*16 ptminvb(kmymaxp,kmymaxp,mintfc)
      complex*16 taua(kmymaxp,kmymaxp,mintfc)
      complex*16 taub(kmymaxp,kmymaxp,mintfc)
      complex*16 gtaua(kmymaxp,kmymaxp,mintfc)
      complex*16 gtaub(kmymaxp,kmymaxp,mintfc)
      complex*16 qmoma(lmsup,mintfc),qmomb(lmsup,mintfc)
      complex*16 qmomla(lmsup,minprc),qmomlb(lmsup,minprc)
      complex*16 qmomra(lmsup,minprc),qmomrb(lmsup,minprc)
c
      common/test/itest
      common/lay2d/cvec(mtotal,3),nextra,nbulkl,nbulkr,
     &             nprc,ninprc(0:mprc+1) 
c
c
      data tiny/1.0d-6/
      data rnul/0.d0/
c     data tolerr,tolef,tolen,tolv0/1.0d-10,1.0d-10,1.0d-9,1.0d-08/
c
c
c     ***************************************************************
c                            INITIALIZE
c     ***************************************************************
c
      open(unit=5,file='input_rsp.in',status='old')
c     ---------------------------------------------------------------
      call readini(
     > imesh,ne,npanel,nepanel,ne1,ne2,ne3,ebottom,etop,eps,lmax,eta,
     > sigma,park,kset,ksetmax,itcpam,cpatol,cpatest,
     > for006,for007,for008,for009,for010,iwrite,
     > leftpot,rightpot,laypot,laycore,leftmat,rightmat,laymat,
     > wrel,lms,bulk,dos,rightm,newvac,ivacpot,v0,v00,
     > vscreen,iscreen,itscfmax,tolerr,tolef,tolen,intbz,kunning,
     > orbpol,opotl,opotr,opot,lliter)
c     ---------------------------------------------------------------
      close(5)
c
      open(unit=5,file='input_geo.in',status='old')
c     ---------------------------------------------------------------
      call readgeom(bulk,bulkgeo,rightm,intbz,nintfc,arot,
     > rsl,concl,rbl,sxl,
     > rs,conc,rba,rbb,sxa,sxb,
     > rsr,concr,rbr,sxr,
     > lfix1,lfix2,linbw,igraph,iesublatt,lsublatt)    
c     ---------------------------------------------------------------
      close(5)
c
c     ---------------------------------------------------------------
      call readclu(nintfc,nimp,nposimp,bthcl,bphcl,rbcl,sxcl)
c     ---------------------------------------------------------------
c
      write(6,'(/2x,''routine MAIN>''/)')
c
      nl=lmax+1
      kmymax=2*nl*nl
      ntotal=nintfc+ninprc(0)*(nextra+1)+ninprc(nprc+1)*(nextra+1)
c
c -read Left,Right and Layer potentials
c
c     ---------------------------------------------------------------
      call pothandle(
     > bulk,linbw,rightm,lmax,nintfc,ninprc(0),ninprc(nprc+1),
     > v00,ivacpot,for006,opotl,opotr,opot,
     > leftpot,concl,qmomla,qmomlb,dxl,nsl,rsl,
     > idpotla,vrla,brla,boprla,zla,idpotlb,vrlb,brlb,boprlb,zlb,
     > rightpot,concr,qmomra,qmomrb,dxr,nsr,rsr,
     > idpotra,vrra,brra,boprra,zra,idpotrb,vrrb,boprrb,brrb,zrb,
     > laypot,conc,dx,ns,rs,
     > idpota,vra,bra,bopra,za,idpotb,vrb,brb,boprb,zb,igraph)
c     ---------------------------------------------------------------
c
c -read in data from update file if exists
c    ----------------------------------------------------------------
      call restart_in(itscf0,etop(npanel),v0,conc,opot,
     >                vra,bra,bopra,rba,vrb,brb,boprb,rbb,
     >                bulk,nintfc,idpota,idpotb,za,zb,laycore,for006)
c    ----------------------------------------------------------------
      write(6,*)
      write(6,'(2x,''etop= '',f17.13)') etop(npanel)
      if(.not.bulk) write(6,*) 'v0=  ',v0
c
c -fix manipulations with files containing inverse cpa t-matrices
c     --------------------------------------------------------------
      call tcpa_main(
     > leftmat,rightmat,laymat,bulk,dos,ksetmax,nintfc,ninprc(0),
     > ninprc(nprc+1),conc,concl,concr,cpamatin,cpamatinl,cpamatinr)
c     --------------------------------------------------------------
c
c -initialize common blocks
c
      call initia
c
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
      write(6,'(/'' Selfconsistent iteration: '',i3/)') itscf
c
      if(bulk) then
c        ---------------------------------------------------------
         call averagebpot(nsl,nintfc,nbulkl,vra,vrla,vrra)
         call averagebpot(nsl,nintfc,nbulkl,vrb,vrlb,vrrb)
         call averagebpot(nsl,nintfc,nbulkl,bra,brla,brra)
         call averagebpot(nsl,nintfc,nbulkl,brb,brlb,brrb)
         if(orbpol) then
         call averagebpotsp(nsl,nintfc,nbulkl,bopra,boprla,boprra)
         call averagebpotsp(nsl,nintfc,nbulkl,boprb,boprlb,boprrb)
         end if
c        ---------------------------------------------------------
      end if
c
c -solve Diraq equation for core states
c  (non spin-polarized!)
c
c     ------------------------------------------------------------
      call chcore(itscfcont,nintfc,conc,idpota,vra,za,
     >idpotb,vrb,zb,dx,ns,rs,qca,qcb,rhoca,rhocb,enca,encb,for006)
c     ------------------------------------------------------------
c
c -generate new energy mesh in the case of bulk iterations
c
      if(itscf.eq.1.or.bulk) then
c        ---------------------------------------------------------
         call zmesh(imesh,npanel,nepanel,ne1,ne2,ne3,ebottom,etop,
     &              eps,cear,we)
c        ---------------------------------------------------------
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
c     -------------------------------------------------------------
      call tcpa_open(
     > for010,leftmat,rightmat,laymat,cpamatin,cpamatinl,cpamatinr)
c     -------------------------------------------------------------
c
c - set quantities for contour integration to zero
c
c     -----------------------------------------------------------
      call initzero(qvpa,qvpb,qva,qvb,
     & spin_magvpa,spin_magvpb,spin_magva,spin_magvb,
     & orb_magvpa,orb_magvpb,orb_magva,orb_magvb,
     & enba,enbb,denba,denbb,enorba,enorbb,qmoma,qmomb,
     & rhova,rhovb,rhospa,rhospb,rhodspa,rhodspb,rhomaga,rhomagb,
     & enbifc,qvifc,omifc) 
c     -----------------------------------------------------------
c
c Rotation matrices between local and global frames of reference
c
c     -------------------------------------------------------------
      call localrot2(
     > lmax,nintfc,rbl,rbr,rba,rbb,vecna,vecnb,phia,phib,
     > dmatl,dmatlp,dmatr,dmatrp,dmata,dmatb,dmatpa,dmatpb,
     > ddpha,ddphb,ddphpa,ddphpb,ddtha,ddthb,ddthpa,ddthpb,
     > rmata,rmatpa,rmatb,rmatpb)
c     -------------------------------------------------------------
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
c     -----------------------------------------------
      call tcpa_in(leftmat,rightmat,laymat,ce,nintfc,
     & ninprc(0),ninprc(nprc+1),conc,concl,concr,
     > cpamatin,cpamatinl,cpamatinr,kmymax,
     & tminv,tminvl,tminvr,kmymaxp)
c     -----------------------------------------------
c
c  generate k-mesh in irreducible Brillouin-zone
c
      if(kset(ie).eq.0) itcpam=0
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
c  t-matrices to start CPA
c     ------------------------------------------------------
      call tmatini(
     > ce,lmax,nintfc,iscreen,vscreen,v0,
     > bulk,wrel,sxl,sxr,sxa,sxb,
     > concl,concr,conc,cpamatin,cpamatinl,cpamatinr,
     > idpotla,vrla,brla,boprla,
     > idpotlb,vrlb,brlb,boprlb,rsl,dxl,nsl,
     > idpotra,vrra,brra,boprra,
     > idpotrb,vrrb,brrb,boprrb,rsr,dxr,nsr,
     > idpota,vra,bra,bopra,idpotb,vrb,brb,boprb,dx,ns,rs,
     > dmatl,dmatlp,dmatr,dmatrp,dmata,dmatb,dmatpa,dmatpb,
     > tminvl,tminvr,tminv,tminva,tminvb,ptminva,ptminvb)
c     ------------------------------------------------------
c       write(6,'('' tmat'')')
c       call flush(6)
c
c  perform CPA, BZ integration & calculate tau-matices
c     ---------------------------------------------
      call cpacoord(
     > iwrite,itscf,ie,ce,lmax,nintfc,eta,rightm,bulk,bulkgeo,wrel,
     > kset(ie),xk,wk,nk,intbz,iek,
     > conc,itcpam,cpatol,cpatest,
     > dmata,dmatb,dmatpa,dmatpb,
     > tminvl,tminvr,tminv,tminva,tminvb,taua,taub,gtaua,gtaub)
c     ---------------------------------------------                        
c       write(6,'('' cpa'')')
c       call flush(6)
c
c -write out (inverse of) effective t-matrices
c     --------------------------------------------------
      call tcpa_out(ce,nintfc,conc,kmymax,tminv,kmymaxp)
c     --------------------------------------------------
c
c  calculate local physical quantities
c  integrate with respect to energy
c     --------------------------------------------------------------
      call locquant(
     > ie,ce,we(ie),lmax,nintfc,wrel,lms,sxa,sxb,conc,v0,
c
     > idpota,vra,bra,bopra,idpotb,vrb,brb,boprb,dx,ns,rs,
c
     > dmata,dmatb,dmatpa,dmatpb,
     > ddpha,ddphb,ddphpa,ddphpb,
     > ddtha,ddthb,ddthpa,ddthpb,
     > rmata,rmatpa,rmatb,rmatpb,
c
     > taua,taub,gtaua,gtaub,ptminva,ptminvb,tminv,
c
     > dosa(1,1,ie),qvpa,qva,qvdiffa(1,ie),
     > enba,enbdiffa(1,ie),denba,enorba,qmoma,
     > rhova,rhospa,rhodspa,rhomaga,
     > dosmaga(1,1,ie),spin_magvpa,spin_magva,orb_magvpa,orb_magva,
c
     > dosb(1,1,ie),qvpb,qvb,qvdiffb(1,ie),
     > enbb,enbdiffb(1,ie),denbb,enorbb,qmomb,
     > rhovb,rhospb,rhodspb,rhomagb,
     > dosmagb(1,1,ie),spin_magvpb,spin_magvb,orb_magvpb,orb_magvb,
c
     > lliter,linbw,efermi,enbifc,qvifc,omifc)
c     --------------------------------------------------------------  
c                                    
 10   continue
c
c                 **************************
c                 *** End of energy loop ***
c                 **************************
c 
c -close binary files containing (inverse of) effective t-matrices
c     -----------------------------------------------------
      call tcpa_close(
     >leftmat,rightmat,laymat,cpamatin,cpamatinl,cpamatinr)
c     -----------------------------------------------------
c
c for DOS calculation print out results and stop
      if(dos) goto 99
c
c  update spin-quantization axes
c     -----------------------------------------------
      call newdir(
     > itscf,itscfcont,nintfc,lfix1,lfix2,
     > rba,rbb,vecna,phia,vecnb,phib,lmax,
     > spin_magvpa,spin_magvpb,spin_magva,spin_magvb,
     > orb_magvpa,orb_magvpb,orb_magva,orb_magvb,lza,lzb,
     > spinmoma,orbmoma,th0a,th1a,ph0a,ph1a,
     > spinmomb,orbmomb,th0b,th1b,ph0b,ph1b)
c     -----------------------------------------------    
c
c  write out results
c     ---------------------------------------------------
      call printscf(
     > linbw,nintfc,itscfcont,conc,enbifc,qvifc,omifc,
     > qva,qca,za,qvb,qcb,zb,enba,enbb,denba,denbb,efermi,
     > arot,rba,rbb,spin_magva,spin_magvb,orb_magva,orb_magvb,
     > spinmoma,orbmoma,th0a,th1a,ph0a,ph1a,
     > spinmomb,orbmomb,th0b,th1b,ph0b,ph1b)
c     ---------------------------------------------------
      if(linbw) goto 200
c
c -find new Fermi level
c
      if(bulk) then
c       -----------------------------------------------------------
        kmy0=1
        call newfl(itscfcont,kmy0,kmymax,kmymaxp,ne,nbulkl,conc,
     >             dosa,dosb,qca,qcb,qva,qvb,za,zb,efermi,defermi)
c       -----------------------------------------------------------
      end if
c
c -Generate new layer potentials 
c
c     -----------------------------------------------------------
      call vgen(
     > itscf,itscfcont,nintfc,lmax,sigma,rightm,bulk,
     > orbpol,opot,conc,concl,concr,
     > za,qca,qva,qmoma,qmomla,qmomra,
     > zb,qcb,qvb,qmomb,qmomlb,qmomrb,
     > rhoca,rhova,rhospa,rhodspa,rhocb,rhovb,rhospb,rhodspb,dx,ns,rs,
     > efermi,defermi,vra,vrb,bra,brb,bopra,boprb,lza,lzb,
     > v0,dv0,ferr1,newvac,
     > enpota,enela,enxca,enpotb,enelb,enxcb,enmaga,enmagb)
c     -----------------------------------------------------------
c
c
c -calculate and print out energies
c     -------------------------------------------------------
      call enprint(itscfcont,nintfc,conc,enca,encb,enba,enbb,
     >  enela,enelb,enxca,enxcb,enpota,enpotb,enmaga,enmagb,
     >  enorba,enorbb,entota,entotb,entot,entotifc)
c     -------------------------------------------------------
      if(itscf.eq.1) then
        dentotifc=10.0d0*tolen
      else
        dentotifc=entotifc-entotifc0
      end if
      entotifc0=entotifc
c
c -update output necessary for restart
c     --------------------------------------------------------------
      call restart_out(itscfcont,efermi,v0,conc,orbpol,
     >                 vra,bra,bopra,rba,vrb,brb,boprb,rbb,
     >                 bulk,nintfc,for006)
c     --------------------------------------------------------------
      if(bulk) etop(npanel)=efermi
c
      if(ferr1.lt.tolerr) istop=1
      if(bulk.and.(dabs(defermi).gt.tolef)) istop=0
      if(bulk.and.(dabs(dentotifc).gt.nintfc*tolen)) istop=0
c     if(rightm.eq.'V'.and.dv0.gt.tolv0) istop=0
c
c
c -Print results
c
   99 continue
c     ---------------------------------------------------------------
      call printres(
     >     for007,for008,itscfcont,lmax,nintfc,ne,cear,efermi,
     >     conc,dosa,dosb,qvpa,qvpb,qva,qvb,dosmaga,dosmagb,
     >     spin_magvpa,spin_magvpb,spin_magva,spin_magvb,
     >     orb_magvpa,orb_magvpb,orb_magva,orb_magvb,
     >     arot,qmoma,qmomb,za,zb,qca,qcb,
     >     qvdiffa,qvdiffb,enbdiffa,enbdiffb,
     >     entota,entotb,entot,entotifc,
     >     idpota,idpotb,vra,bra,bopra,vrb,brb,boprb,
     >     rs,dx,ns,lms,orbpol,bulk)
c     ---------------------------------------------------------------
      if(itscf.lt.itscfmax.and.istop.eq.0) goto 100
c
c
c     ****************************************************************
c                    END OF SELFCONSISTENT ITERATIONS 
c     ****************************************************************
c   
  200 close(6)
      stop
      end
