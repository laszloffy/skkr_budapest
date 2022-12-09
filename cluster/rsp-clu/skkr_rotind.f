c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      program main
c     implicit real*8 (a-h,o-z)
      implicit none
c
      include '../param.h'
c
      logical bulk
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
      real*8 EFERMI,ENBIFC,QVIFC,FERR1,ENTOTIFC
      real*8 DENTOTIFC,ENTOTIFC0,CVEC

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
      real*8 lza(2,mimp)
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
      real*8 rmt(mimp)
c
      real*8 za(mimp)
      real*8 sxa(mimp)
      real*8 qca(mimp)
      real*8 qva(mimp)
      real*8 qvpa(kmymaxp,mimp)
c
      real*8 zb(mimp)
      real*8 qcb(mimp)
      real*8 encb(mimp)
      real*8 rhocb(nrad,mimp)
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
      real*8 dosa(kmymaxp,mimp,me)
      real*8 dosmaga(kmymaxp,mimp,me)
      real*8 enbdiffa(mimp,me)
      real*8 qvdiffa(mimp,me)
c
      real*8 spin_magvpa(kmymaxp,mimp,3)
      real*8 spin_magva(mimp,3)
      real*8 orb_magvpa(kmymaxp,mimp,3)
      real*8 orb_magva(mimp,3)
      real*8 enba(mimp)
      real*8 denba(mimp)
      real*8 enca(mimp)
      real*8 enela(mimp)
      real*8 enxca(mimp)
      real*8 enpota(mimp)
      real*8 enmaga(mimp)
      real*8 enorba(mimp)
      real*8 entota(mimp)
      real*8 entot(mimp)
      real*8 rhoca(nrad,mimp)
      real*8 rhova(nrad,mimp)
      real*8 rhospa(nrad,2,mimp)
      real*8 rhodspa(nrad,2,mimp)
      real*8 rhomaga(nrad,mimp)
c
      real*8 spinmoma(mimp)
      real*8 orbmoma(mimp)
      real*8 th0a(mimp)
      real*8 th1a(mimp)
      real*8 ph0a(mimp)
      real*8 ph1a(mimp)
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
      complex*16 qmomla(lmsup,minprc)
      complex*16 qmomra(lmsup,minprc)
c
c -----------immpurities----------
      complex*16 detl,det0(me)
      real*8     ebl
c
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
      real*8 rbcl(3,mimp)
      real*8 bthcl(mimp)
      real*8 bphcl(mimp)
      real*8 sxcl(mimp)
c
      real*8 vmadih(mimp)
      real*8 vmadich(mimp)
      real*8 vmadid(mimp)
      real*8 vmadiq(mimp)
c
      complex*16 tau_ii(kmymaxp,kmymaxp,mintfc,net)
      complex*16 tminvh(kmymaxp,kmymaxp,mintfc,net)
c
      complex*16 tau_ij(kmymaxp,kmymaxp,mpair,net)
      complex*16 tau_ji(kmymaxp,kmymaxp,mpair,net)
      complex*16 taucl(kmymaxp,kmymaxp,mimp)
      complex*16 gtaucl(kmymaxp,kmymaxp,mimp)
c
      complex*16 tauclij(kmymaxp,kmymaxp,mimp,mimp)      
c
      complex*16 ze
c
      real*8 vscreen
      real*8 vscreen0
c variables for the derivatives of the free energy
      complex*16 d2dph(kmymaxp,kmymaxp,mimp)
      complex*16 d2dphp(kmymaxp,kmymaxp,mimp)
      complex*16 d2dth(kmymaxp,kmymaxp,mimp)
      complex*16 d2dthp(kmymaxp,kmymaxp,mimp)
      complex*16 d2dthph(kmymaxp,kmymaxp,mimp)
      complex*16 d2dthphp(kmymaxp,kmymaxp,mimp)
      real*8     theta0,phi0
c
      complex*16 dtmph(kmymaxp,kmymaxp,mimp)
      complex*16 dtmth(kmymaxp,kmymaxp,mimp)
      complex*16 d2tmph(kmymaxp,kmymaxp,mimp)
      complex*16 d2tmth(kmymaxp,kmymaxp,mimp)
      complex*16 d2tmthph(kmymaxp,kmymaxp,mimp)
c
      real*8     vij(4,mimp,mimp),gradth(mimp),gradph(mimp)      
      real*8     dedmx(mimp),dedmy(mimp),dedmz(mimp)
c---------------------------------
      integer    iwtau, nwtau
c---------------------------------
c Laszloffy 06/04/17
      logical isscf,rotind,spinorb
      real*8  dmmin
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
     >             theta0,phi0, iwtau, nwtau)
c     ---------------------------------------------------------------
      do iimp=1,nimp
        rs(iimp)=rs0(nposimp(3,iimp))
      end do
c
      write(6,'(/2x,''routine MAIN>''/)')
c
      nl=lmax+1
      kmymax=2*nl*nl
      ntotal=nintfc+ninprc(0)*(nextra+1)+ninprc(nprc+1)*(nextra+1)
c
c -read Impuity potentials
c
c     ---------------------------------------------------------------
      call pothandleimp(
     > lmax,nimp,v00,ivacpot,for006,opot,imppot,dx,ns,rs,idpota,
     > vra,bra,bopra,za,ib0,b0,impb0,impb0f,isigb)
c     ---------------------------------------------------------------
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
c - read in the host t and tau matrices
      do ie=1,ne
        ze=cear(ie)
        call tauhin(ie,ze,npair,nintfc,kmymax,
     >                  tau_ij(1,1,1,ie),
     >                  tau_ji(1,1,1,ie),
     >                  tau_ii(1,1,1,ie),taupath,imesh)
        call thin(ie,ze,nintfc,kmymax,tminvh(1,1,1,ie),taupath)
      end do
c
      call vhin(nimp,vmadih,v0,iscreencl,vscreen,kpairind,taupath)
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
c set det0 array
      do ie = 1,ne
       det0(ie) = (0.d0, 0.d0)
      end do
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
      ebl = 0.d0
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
c test localrot
c     write(6,*) 'ddph'
c     call testrot(dmata,dmatpa,ddpha,ddphpa,nimp,kmymax)
c     write(6,*) 'ddth'
c     call testrot(dmata,dmatpa,ddtha,ddthpa,nimp,kmymax)
c test localrot vege 
       if(iwtau.eq.1) then
        open(unit=16,file='adatok.bin',form='unformatted') 
        write(16) kmymax,ne
       endif  

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
      call derivt(nimp,kmymax,ptminv0,
     > dmata,dmatpa,
     > ddpha,ddphpa,d2dph,d2dphp,
     > ddtha,ddthpa,d2dth,d2dthp,
     > d2dthph,d2dthphp,
     > dtmph, dtmth, d2tmph, d2tmth, d2tmthph)
c
c  calculate tau-matices
c =====================================================================
      call ecoreimp(nimp,iscreencl,tau_ij(1,1,1,ie),tau_ji(1,1,1,ie),
     >              tau_ii(1,1,1,ie),ptminva,
     >              tminvh(1,1,1,ie),lmax,nposimp,kpairind,dmata,dmatpa,
     >              gtaucl,taucl,tauclij,detl)
c =====================================================================
c      if(itest.ge.2) call timit(tindex,6) 
c    band energy from Lloyd's formula
        write(6,*) 'det:',ie,ce,detl
        call setphase(det0(ie),detl)
        ebl = ebl + dimag(we(ie)*detl) 
c
c       write(6,'('' cpa'')')
c       call flush(6)
c
c  calculate first and second derivatives of the band energy
        call jij(nimp,kmymax,we(ie),tauclij,
     >               dtmph,dtmth,d2tmph,d2tmth,d2tmthph,
     >               vij,gradph,gradth)
c =====================================================================
c
c  calculate local physical quantities
c  integrate with respect to energy
c     --------------------------------------------------------------
      call locquant2(
     > ie,ce,we(ie),lmax,madmax,nimp,wrel,lms,sxcl,v0,
     > idpota,vra,bra,bopra,dx,ns,rs,nimp,nposimp,
     > dmata,dmatpa,
     > ddpha,ddphpa,
     > ddtha,ddthpa,
     > rmata,rmatpa,
     > taucl,gtaucl,ptminva,
     > dosa(1,1,ie),qvpa,qva,qvdiffa(1,ie),
     > vmadid,vmadiq,enba,enbdiffa(1,ie),denba,enorba,qmoma,
     > rhova,rhospa,rhodspa,rhomaga,
     > dosmaga(1,1,ie),spin_magvpa,spin_magva,orb_magvpa,orb_magva,
     > lliter,efermi,enbifc,qvifc,omifc, iwtau, nwtau)
c     --------------------------------------------------------------  
c                                    
 10   continue
      if(iwtau.eq.1) close(16)
c
c                 **************************
c                 *** End of energy loop ***
c                 **************************
c 
      call gradm(rba,nimp,gradth,gradph,theta0,phi0,
     >                 dedmx,dedmy,dedmz)
      call printjij(vij,gradth,gradph,dedmx,dedmy,dedmz,nimp)
c        

c for DOS calculation print out results and stop
      if(dos) goto 99
c
c Write band energy from Lloyd's formula
      ebl = ebl/3.141592653589793d0
      write(6,*) 'Band energy:',ebl
c  update spin-quantization axes
c     -----------------------------------------------
c     call newdir1(
c    > itscf,itscfcont,nimp,lfix1,lfix2,rba,vecna,phia,
c    > lmax,spin_magvpa,spin_magva,orb_magvpa,orb_magva,
c    > lza,spinmoma,orbmoma,th0a,th1a,ph0a,ph1a,
c    > rotind,spinorb,dmmin)
c
c  Kikommentezve milyen eredmenyt ad?
c
      call newdir1(
     > nimp,rba,vecna,phia,lmax,spin_magva,orb_magva,
     > spinmoma,orbmoma,
     > rotind,spinorb,dmmin)
c 
c
c
c     -----------------------------------------------    
c
c  write out results
c     ---------------------------------------------------
      call printscf1(
     > nimp,itscfcont,enbifc,qvifc,omifc,
     > qva,qca,za,enba,denba,efermi,
     > rba,spin_magva,orb_magva,
     > spinmoma,orbmoma,th0a,th1a,ph0a,ph1a)
c     ---------------------------------------------------
c -Generate new layer potentials 
c
c     -----------------------------------------------------------
c
c Laszloffy
c only if scf calculations performed
c
      if(isscf) then
      call vgen1(
     > itscf,itscfcont,nimp,nposimp,
     > orbpol,opot,madmax,
     > za,qca,qva,
     > rhoca,rhova,rhospa,rhodspa,dx,ns,rs,
     > vra,bra,bopra,lza,
     > ferr1,
     > enpota,enela,enxca,enmaga,
     > vmadih,vmadich,vmadid,vmadiq)
      end if
c     -----------------------------------------------------------
c
c
c -calculate and print out energies
c     -------------------------------------------------------
      call enprint1(itscfcont,nimp,enca,enba,
     >  enela,enxca,enpota,enmaga,
     >  enorba,entota,entot,entotifc)
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
      call restart_out1(itscfcont,efermi,v0,orbpol,
     >                 vra,bra,bopra,rba,
     >                 nimp,for006)
c     --------------------------------------------------------------
c     if(bulk) etop(npanel)=efermi
c
c Laszloffy
c
      if(isscf.and.(ferr1.lt.tolerr)) istop=1
c
c -Print results
c
   99 continue
c     ---------------------------------------------------------------
      call printres1(
     >     for007,for008,itscfcont,lmax,nimp,ne,cear,efermi,
     >     dosa,qvpa,qva,dosmaga,
     >     spin_magvpa,spin_magva,
     >     orb_magvpa,orb_magva,
     >     arot,qmoma,za,qca,
     >     qvdiffa,enbdiffa,
     >     entota,entot,entotifc,
     >     idpota,vra,bra,bopra,
     >     rs,dx,ns,lms,orbpol)
c     ---------------------------------------------------------------
      if(itscf.lt.itscfmax.and.istop.eq.0) goto 100
c
c
c     ****************************************************************
c                    END OF SELFCONSISTENT ITERATIONS 
c     ****************************************************************
c   
! 200 close(6)
      stop
      end
