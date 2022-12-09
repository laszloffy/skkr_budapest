c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      program main
c     implicit real*8 (a-h,o-z)
      implicit none
c
      include '../param.h'
c
      logical dos
      logical wrel
      logical lms
      logical orbpol
      logical opot
c
      character*10 idpota(mimp)
      character*30 for006
      character*30 for007
      character*30 for008
      character*30 for009
      character*30 for010
      character*30 laycore
c
      real*8 TINY,RNUL,ETA,SIGMA,TOLERR
      real*8 TOLEF,TOLEN,AROT,V00,V0
      real*8 EFERMI,ENBIFC,QVIFC,FERR1
      real*8 CVEC

      integer*4 IMESH,NE,NPANEL,NE1,NE2
      integer*4 NE3,LMAX,ISCREEN,ITSCFMAX,LLITER
      integer*4 NINTFC,NIMP,LFIX1,LFIX2
      integer*4 IIMP,NL,KMYMAX,NTOTAL,NINPRC
      integer*4 NEXTRA,NPRC,IVACPOT,ITSCF0,ITEST
      integer*4 IE,ITSCF,ISTOP
      integer*4 ITSCFCONT,OMIFC,IE0,IEK
      integer*4 NBULKL
      integer*4 NBULKR
c
c
      integer nepanel(5)
      integer ns(mimp)
c
      real*8 ebottom(5)
      real*8 etop(5)
      real*8 eps(5)
c
      real*8 vra(nrad,mimp)
      real*8 bra(nrad,mimp)
      real*8 bopra(nrad,2,mimp)
      real*8 rs(mimp)
      real*8 rs0(mintfc)
      real*8 dx(mimp)
c
      real*8 za(mimp)
      real*8 sxa(mimp)
      real*8 qca(mimp)
      real*8 qva(mimp)
      real*8 qvpa(kmymaxp,mimp)
c
c
      real*8 vecna(3,mimp)
      real*8 rba(3,mimp)
      real*8 phia(mimp)
c
      complex*16 dmata(kmymaxp,kmymaxp,mimp)
      complex*16 dmatpa(kmymaxp,kmymaxp,mimp)
      complex*16 rmata(lmsup,lmsup,mimp)
      complex*16 rmatpa(lmsup,lmsup,mimp)
c
      complex*16 ddpha(kmymaxp,kmymaxp,mimp)
      complex*16 ddphpa(kmymaxp,kmymaxp,mimp)
      complex*16 ddtha(kmymaxp,kmymaxp,mimp)
      complex*16 ddthpa(kmymaxp,kmymaxp,mimp)
c
      real*8 dosa(kmymaxp,mimp)
c
      real*8 spin_magvpa(kmymaxp,mimp,3)
      real*8 spin_magva(mimp,3)
      real*8 orb_magvpa(kmymaxp,mimp,3)
      real*8 orb_magva(mimp,3)
      real*8 enba(mimp)
      real*8 denba(mimp)
      real*8 enca(mimp)
      real*8 enorba(mimp)
      real*8 rhoca(nrad,mimp)
      real*8 rhova(nrad,mimp)
      real*8 rhospa(nrad,2,mimp)
      real*8 rhodspa(nrad,2,mimp)
      real*8 rhomaga(nrad,mimp)
c
      complex*16 cear(me)
      complex*16 we(me)
      complex*16 ce
c
      complex*16 tminva(kmymaxp,kmymaxp,mimp)
      complex*16 ptminva(kmymaxp,kmymaxp,mimp)
      complex*16 ptminv0(kmymaxp,kmymaxp,mimp)
c
c     complex*16 taua(kmymaxp,kmymaxp,mimp)
c     complex*16 gtaua(kmymaxp,kmymaxp,mimp)
c
      complex*16 qmoma(lmsup,mimp)
c
c -----------immpurities----------
      complex*16 detl
      character*30 imppot
      character*15 taupath
c
      real*8 b0
c
      integer ib0
      integer impb0
      integer impb0f
      integer isigb(mimp)
c
      integer iscreencl
      integer nposimp(3,mimp)
      integer kpairind(mimp,mimp)
      integer kpair(4,mpair)
      integer npair
      integer madmax
c
      real*8 bthcl(mimp)
      real*8 bphcl(mimp)
      real*8 sxcl(mimp)
c
      real*8 vmadid(mimp)
      real*8 vmadiq(mimp)
c
      complex*16 tau_ii(kmymaxp,kmymaxp,mintfc)
      complex*16 tminvh(kmymaxp,kmymaxp,mintfc)
c
      complex*16 tau_ij(kmymaxp,kmymaxp,mpair)
      complex*16 tau_ji(kmymaxp,kmymaxp,mpair)
      complex*16 taucl(kmymaxp,kmymaxp,mimp)
      complex*16 gtaucl(kmymaxp,kmymaxp,mimp)
c
      complex*16 tauclij(kmymaxp,kmymaxp,mimp,mimp)      
c
      complex*16 ze
c
      real*8 vscreen
c variables for the derivatives of the free energy
      complex*16 d2dph(kmymaxp,kmymaxp,mimp)
      complex*16 d2dphp(kmymaxp,kmymaxp,mimp)
      complex*16 d2dth(kmymaxp,kmymaxp,mimp)
      complex*16 d2dthp(kmymaxp,kmymaxp,mimp)
      complex*16 d2dthph(kmymaxp,kmymaxp,mimp)
      complex*16 d2dthphp(kmymaxp,kmymaxp,mimp)
      real*8     theta0,phi0
c
      real*8     vij(4,mimp,mimp),gradth(mimp),gradph(mimp)      
c---------------------------------
      integer    iwtau, nwtau
      integer    momsite(mimp)
c---------------------------------
c Laszloffy 06/04/17
      logical isscf,rotind,spinorb
      real*8  dmmin
c---------------------------------
c
      common/test/itest
      common/lay2d/cvec(mtotal,3),nextra,nbulkl,nbulkr,
     &             nprc,ninprc(0:mprc+1) 
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
     > imesh,ne,npanel,nepanel,ne1,ne2,ne3,ebottom,etop,eps,lmax,
     > eta,sigma,!park,kset,ksetmax,itcpam,cpatol,cpatest,
     > for006,for007,for008,for009,!for010,iwrite,
c    > leftpot,rightpot,laypot,
     > laycore,!leftmat,rightmat,laymat,
     > wrel,lms,!bulk,
     > dos,!rightm,newvac,ivacpot,v0,v00,
     > vscreen,iscreen,
     > itscfmax,tolerr,tolef,tolen,!intbz,kunning,
     > orbpol,!opotl,opotr,
     > opot,lliter,
     > isscf,rotind,spinorb,dmmin)
c     ---------------------------------------------------------------
      close(5)
c
      open(unit=5,file='input_geo.in',status='old')
c     ---------------------------------------------------------------
      call readgeom(nintfc,arot,rs0)
c     call readgeom(bulk,bulkgeo,rightm,intbz,nintfc,arot,
c    > rsl,concl,rbl,sxl,
c    > rs,conc,rba,rbb,sxa,sxb,
c    > rsr,concr,rbr,sxr,
c    > lfix1,lfix2,linbw,igraph,iesublatt,lsublatt)    
c     ---------------------------------------------------------------
      close(5)
c
c     ---------------------------------------------------------------
      call readclu(nintfc,nimp,nposimp,npair,kpair,kpairind,
     >             iscreencl,bthcl,bphcl,rba,sxcl,taupath,
     >             imppot,madmax,lfix1,lfix2,
     >             ib0,b0,impb0,impb0f,isigb,
     >             theta0,phi0, iwtau, nwtau,
     >             momsite)
c     ---------------------------------------------------------------
      do iimp=1,nimp
        rs(iimp)=rs0(nposimp(3,iimp))
      end do
c
      write(6,'(/2x,''routine MAIN>''/)')
c
      nl=lmax+1
      kmymax=2*nl*nl
      write(6,*) "nprc=",nprc
      ntotal=nintfc+ninprc(0)*(nextra+1)+ninprc(nprc+1)*(nextra+1)
c
c -read Impuity potentials
c
c     ---------------------------------------------------------------
      call pothandleimp(
     > lmax,nimp,v00,ivacpot,for006,opot,imppot,dx,ns,rs,idpota,
     > vra,bra,bopra,za,ib0,b0,impb0,impb0f,isigb)
c     ---------------------------------------------------------------
      do iimp=1,nimp
        if (momsite(iimp)/=1) then
          !ignore induced Weiss-fields at non-magnetic sites
          bra(:,iimp) = 0.d0
          write(6,'(''<skkr> : Effective field set to 0
     >           on non-magnetic site '',i3)') iimp
        end if
      end do
!     do imp=1,nimp
!       if(itest.ge.2) then
!        write(6,*) '<skkr>: after pothandleimp'
!         write(6,*) '<skkr>: imp=',imp
!         write(6,*) '<skkr>: idpota=',idpota(imp)
!         write(6,*) '<skkr>: v0=',v0
!         write(6,*) '<skkr>: dx=',dx(imp)
!         write(6,*) '<skkr>: ns=',ns(imp)
!         write(6,*) '<skkr>: rs=',rs(imp)
!         write(6,*) '<skkr>: sxcl=',sxcl(imp)
!       end if
!     end do
c
c
c -read in data from update file if exists
c    ----------------------------------------------------------------
      call restart_in1(itscf0,etop(npanel),v0,opot,
     >                vra,bra,bopra,rba,
     >                nimp,idpota,za,laycore,for006)
c    ----------------------------------------------------------------
c
c -initialize common blocks
c
      call initia
c
c -generate energy mesh (fix energy contour)
c
c     ---------------------------------------------------------
      call zmesh(imesh,npanel,nepanel,ne1,ne2,ne3,ebottom,etop,
     >           eps,cear,we)
c     ---------------------------------------------------------
      if(itest.ge.2) then 
        write(6,*)
        do ie=1,ne
          write(6,'('' ie='',i2,'' e='',2f12.8,'' w='',2f12.8)')
     >          ie,cear(ie),we(ie)
        end do
        write(6,*) 'CROSSCHECK RSCL!!!!!!!!!'
        write(6,*) 'GIVE MADMAX!!!!!!!!!!!!!'
        write(6,*) 'EMAD!!!!!!!!!!!!!'
        write(6,*) 'direction mixing!!!!!!!!!!!!!'
!       write(6,*) 'READ V0!!!!!!!!!!!!'
!       write(6,*) 'SCRPAR!!!!!!!!!!!!!'
!       write(6,*) 'RESTART!!!!!!!!!!!!!'
      end if
c
c =====================================================================
c
c     call vhin(nimp,vmadih,v0,iscreencl,vscreen,kpairind,taupath)
c =====================================================================
      write(6,*)
      write(6,'(2x,''etop= '',f17.13)') etop(npanel)
      write(6,*) 'v0=  ',v0
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
c -solve Diraq equation for core states
c  (non spin-polarized!)
c
c     ------------------------------------------------------------
      call chcore(itscfcont,nimp,idpota,vra,za,
     >            dx,ns,rs,qca,rhoca,enca,for006)
c     ------------------------------------------------------------
c
c - set quantities for contour integration to zero
c     -----------------------------------------------------------
      call initzero1(qvpa,qva,vmadid,vmadiq,
     & spin_magvpa,spin_magva,
     & orb_magvpa,orb_magva,
     & enba,denba,enorba,qmoma,
     & rhova,rhospa,rhodspa,rhomaga,
     & enbifc,qvifc,omifc) 
c set vij, gradth, gradph to zero      
      call rzero(vij,4*mimp*mimp)
      call rzero(gradth,mimp)
      call rzero(gradph,mimp)
c     -----------------------------------------------------------
c
c Rotation matrices between local and global frames of reference
c
c     -------------------------------------------------------------
      call localrot2imp(
     > lmax,nimp,rba,vecna,phia,dmata,dmatpa,
     > ddpha,ddphpa,ddtha,ddthpa,rmata,rmatpa,
     > theta0,phi0,
     > d2dph, d2dphp, d2dth, d2dthp, d2dthph, d2dthphp)  
c     -------------------------------------------------------------
c                 *************************
c                 *** Start energy loop ***
c                 *************************
c
      ie0=1
      iek=0
      imesh = 0
      open(8,file=for008,form='formatted')
      do 10 ie=1,ne
c
      ce=cear(ie)
        call tauhin(ie,ce,npair,nintfc,kmymax,
     >                  tau_ij(1,1,1),
     >                  tau_ji(1,1,1),
     >                  tau_ii(1,1,1),taupath,imesh)
        call thin(ie,ce,nintfc,kmymax,tminvh(1,1,1),taupath)
c
c  t-matrices 
c     ------------------------------------------------------
      call tmatiniimp(
     > ce,lmax,nintfc,nimp,nposimp,iscreencl,vscreen,v0,
     > wrel,sxcl,idpota,vra,bra,bopra,dx,ns,rs,
     > dmata,dmatpa,tminva,ptminva,ptminv0)
c     ------------------------------------------------------
c       write(6,'('' tmat'')')
c       call flush(6)
c  calculate derivatives of the physical t^{-1} matrices
c =====================================================================
c     call derivt(nimp,kmymax,ptminv0,
c    > dmata,dmatpa,
c    > ddpha,ddphpa,d2dph,d2dphp,
c    > ddtha,ddthpa,d2dth,d2dthp,
c    > d2dthph,d2dthphp,
c    > dtmph, dtmth, d2tmph, d2tmth, d2tmthph)
c
c  calculate tau-matices
c =====================================================================
      call ecoreimp(nimp,iscreencl,tau_ij(1,1,1),tau_ji(1,1,1),
     >              tau_ii(1,1,1),ptminva,
     >              tminvh(1,1,1),lmax,nposimp,kpairind,dmata,dmatpa,
     >              gtaucl,taucl,tauclij,detl)
c =====================================================================
c       write(6,'('' cpa'')')
c       call flush(6)
c
c  calculate first and second derivatives of the band energy
c       call jij(nimp,kmymax,we(ie),tauclij,
c    >               dtmph,dtmth,d2tmph,d2tmth,d2tmthph,
c    >               vij,gradph,gradth)
c =====================================================================
c
c  calculate local physical quantities
c  integrate with respect to energy
c     --------------------------------------------------------------
      call locquant_dos(
     > ce,lmax,nimp,wrel,lms,sxcl,v0,
     > idpota,vra,bra,bopra,dx,ns,rs,
     > taucl, dosa)
c
c     call locquant1(
c    > ie,ce,we(ie),lmax,madmax,nimp,wrel,lms,sxcl,v0,
c    > idpota,vra,bra,bopra,dx,ns,rs,nimp,nposimp,
c    > dmata,dmatpa,
c    > ddpha,ddphpa,
c    > ddtha,ddthpa,
c    > rmata,rmatpa,
c    > taucl,gtaucl,ptminva,
c    > dosa(1,1,ie),qvpa,qva,qvdiffa(1,ie),
c    > vmadid,vmadiq,enba,enbdiffa(1,ie),denba,enorba,qmoma,
c    > rhova,rhospa,rhodspa,rhomaga,
c    > dosmaga(1,1,ie),spin_magvpa,spin_magva,orb_magvpa,orb_magva,
c    > lliter,efermi,enbifc,qvifc,omifc, iwtau, nwtau)
c     --------------------------------------------------------------  
c                                    
 10   continue
      close(8)
c
c                 **************************
c                 *** End of energy loop ***
c                 **************************
c     ---------------------------------------------------------------
c     call printres1(
c    >     for007,for008,itscfcont,lmax,nimp,ne,cear,efermi,
c    >     dosa,qvpa,qva,dosmaga,
c    >     spin_magvpa,spin_magva,
c    >     orb_magvpa,orb_magva,
c    >     arot,qmoma,za,qca,
c    >     qvdiffa,enbdiffa,
c    >     entota,entot,entotifc,
c    >     idpota,vra,bra,bopra,
c    >     rs,dx,ns,lms,orbpol)
c     ---------------------------------------------------------------
c
c     ****************************************************************
c                    END OF SELFCONSISTENT ITERATIONS 
c     ****************************************************************
c   
! 200 close(6)
      stop
      end
