c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
c locquant for cluster calculation
      subroutine locquant2(
c     ===================
     > ie,ce,we,lmax,madmax,nintfc,wrel,lms,sxa,v0,
     > deltaa,singratimp,dratimp,uratimp,     
     > idpota,vra,bra,rba,bopra,dx,ns,rs,nimp,nposimp,
     > taua,!gtaua,ptminva,
     > dosa,dosha,doseha,doshea,dosteha,dosthea,
     > tripletmat,
     > qvpa,qvhpa,
     > qva,qvha,qveha,qvhea,qvteha,qvthea,
     > qmomha,qmomeha,qmomhea,qmomteha,qmomthea,
     > vmadid,vmadiq,
     > enba,enbha,enbdiffa,enbhdiffa,
!    > denba,!nenba, ! commented on 07/07/2021 Laszloffy
!     > enorba, !commented on 30/06/2022 bnyari
     > qmoma,
     > rhova,rhospa,rhodspa,rhomaga,
     > dosmaga,dosmagha,
     > spin_magvpa, spin_magvhpa,spin_magva,spin_magvha,
     > orb_magvpa,orb_magvhpa,orb_magva,orb_magvha,
     > lliter,enbifc,enbhifc,qvifc,qvhifc,omifcint, iwtau, nwtau,
     > reg,c_light,E_Fermi,tripletout)
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      logical wrel,lms,reg
      logical tautest
c
      character*10 idpota(nimp)
      character*30 for006
c
      integer nuz(dbogomaxp)
      integer indz(nuzp,dbogomaxp)
      integer mpar(-lsup:lsup)
      integer ns(nimp)
c
      real*8 vra(nrad,nimp)
      real*8 bra(nrad,nimp)
      real*8 rba(3,nintfc)
      complex*16 deltaa(nrad,nimp)
      dimension singratimp(nimp),uratimp(nimp),dratimp(nimp)
      real*8 bopra(nrad,2,nimp)
      real*8 rs(nimp)
      real*8 dx(nimp)
c
      real*8 rtmp(3)
      real*8 rtmph(3)
c
      real*8 enba(nimp)
      real*8 enbha(nimp)
c     real*8 denba(nimp) ! commented on 07/07/2021 Laszloffy
c     complex*16 nenba(nimp,me)
!      real*8 enorba(nimp)
c
      real*8 dosa(kmymaxp,nimp)
      real*8 dosha(kmymaxp,nimp)
      real*8 doseha(kmymaxp,nimp),dosteha(kmymaxp,nimp)
      real*8 doshea(kmymaxp,nimp),dosthea(kmymaxp,nimp)
      real*8 dosmaga(kmymaxp,nimp)
      real*8 dosmagha(kmymaxp,nimp)
      real*8 qva(nimp), qvpa(kmymaxp,nimp)
      real*8 qvha(nimp),qvhpa(kmymaxp,nimp)
      real*8 qveha(nimp),qvhea(nimp)
      real*8 qvteha(nimp),qvthea(nimp)
      real*8 tripletmat(lmmaxp,lmmaxp,3,nimp)
      dimension spin_magvpa(kmymaxp,nimp,3)
      dimension spin_magva(nimp,3)
      dimension orb_magvpa(kmymaxp,nimp,3)
      dimension orb_magva(nimp,3)
      dimension spin_magvhpa(kmymaxp,nimp,3)
      dimension spin_magvha(nimp,3)
      dimension orb_magvhpa(kmymaxp,nimp,3)
      dimension orb_magvha(nimp,3)
c
      real*8 enbdiffa(nimp)
      real*8 enbhdiffa(nimp)
      real*8 qvdiffa(nimp)
      real*8 qvhdiffa(nimp)
      real*8 enbifc,enbhifc
      real*8 qvifc,qvhifc
c
      real*8 rhova(nrad,nimp)
      real*8 rhospa(nrad,2,nimp)
      real*8 rhodspa(nrad,2,nimp)
      real*8 rhomaga(nrad,nimp)
c      real*8 rhov(nrad)
c      real*8 rad(nrad)
c
      real*8 sxa(nimp)
c
      real*8 vmadid(nimp)
      real*8 vmadiq(nimp)
c
      complex*16 ce
      complex*16 we
      complex*16 psq
      complex*16 cih
      complex*16 temp
c
c ======================================================================
      complex*16 mlmij(lmsup,nimp)
c
      integer nposimp(3,nimp)
      integer madmax
      integer nimp
c
c ======================================================================
c
      complex*16 taua(dbogomaxp,dbogomaxp,nimp)
c     complex*16 gtaua(kmymaxp,kmymaxp,nimp)
c     complex*16 ptminva(kmymaxp,kmymaxp,nimp)
      complex*16 tm(dbogomaxp,dbogomaxp)
c
      complex*16 gz(dbogomaxp,dbogomaxp,nrad)
      complex*16 fz(dbogomaxp,dbogomaxp,nrad)
      complex*16 gj(dbogomaxp,dbogomaxp,nrad)
      complex*16 fj(dbogomaxp,dbogomaxp,nrad)
      complex*16 glz(dbogomaxp,dbogomaxp,nrad)
      complex*16 flz(dbogomaxp,dbogomaxp,nrad)
      complex*16 glj(dbogomaxp,dbogomaxp,nrad)
      complex*16 flj(dbogomaxp,dbogomaxp,nrad)
c
      complex*16 ggrmat(dbogomaxp,dbogomaxp,0:lsup)
      complex*16 fgrmat(dbogomaxp,dbogomaxp,0:lsup)
      complex*16 tggrmat(kmymaxp,kmymaxp,0:lsup)
      complex*16 tfgrmat(kmymaxp,kmymaxp,0:lsup)
c
      complex*16 qmoma(lmsup,nimp)
c     complex*16 zdos(dbogomaxp)
      complex*16 zrho(nrad),zmagrho(nrad)
      complex*16 zrhosp(nrad,2),zrhodsp(nrad,2)
      complex*16 zrhoeh(nrad),zrhohe(nrad)
      complex*16 zdmag(dbogomaxp,3),zdmor(dbogomaxp,3)
      complex*16 zqmom(dbogomaxp,lmsup)
      complex*16 zscmom(dbogomaxp,lmsup),zsctmom(dbogomaxp,lmsup)
      complex*16 zscdos(dbogomaxp),zsctdos(dbogomaxp)
      complex*16 zenorb
      complex*16 qmom(lmsup),qmomh(lmsup)
      complex*16 qmomrot(lmsup),qmomroth(lmsup)
      complex*16 qmomeh,qmomhe,qmomthe,qmomteh
c
      complex*16 qmomha(lmsup,nimp)
      complex*16 qmomeha(nimp),qmomhea(nimp)
      complex*16 qmomteha(nimp),qmomthea(nimp)
cc
      complex*16 deltakmy(kmymaxp,kmymaxp)
      complex*16 rgacoeff(kmymaxp,kmymaxp,lmsup)
      common/rggaunt/rgacoeff
c
      complex*16 scoeff(kmymaxp,kmymaxp,3)
      complex*16 sbcoeff(kmymaxp,kmymaxp,3)
      complex*16 sxcoeff(kmymaxp,kmymaxp)
      complex*16 sxbcoeff(kmymaxp,kmymaxp)
      complex*16 sycoeff(kmymaxp,kmymaxp)
      complex*16 sybcoeff(kmymaxp,kmymaxp)
      complex*16 szcoeff(kmymaxp,kmymaxp)
      complex*16 szbcoeff(kmymaxp,kmymaxp)
      common/sigmat/sxcoeff,sxbcoeff,sycoeff,sybcoeff,szcoeff,szbcoeff
c
      complex*16 lcoeff(kmymaxp,kmymaxp,3)
      complex*16 lbcoeff(kmymaxp,kmymaxp,3)
      complex*16 lxcoeff(kmymaxp,kmymaxp)
      complex*16 lxbcoeff(kmymaxp,kmymaxp)
      complex*16 lycoeff(kmymaxp,kmymaxp)
      complex*16 lybcoeff(kmymaxp,kmymaxp)
      complex*16 lzcoeff(kmymaxp,kmymaxp)
      complex*16 lzbcoeff(kmymaxp,kmymaxp)
c
      character*34 tripletout
      integer    tindex
      real       t1,t2
      common/lmat/lxcoeff,lxbcoeff,lycoeff,lybcoeff,lzcoeff,lzbcoeff
c
c     complex*16 alphalkkr(0:lmaxp,minprc)
c     complex*16 alpharkkr(0:lmaxp,minprc)
c     complex*16 alphaintkkr(0:lmaxp,mintfc)
c     common/scrpar/alphaintkkr
c
      common/test/itest
c
      data tol/1.0d-8/,cih/(0.0d0,-0.5d0)/
      data tiny/1.0d-6/
c
c ********************
c initialize constants
c ********************
c
c---> c in rydberg units:
c      c=274.072d0
      write(6,*) '<locquant2> Ef=',E_Fermi
      c=c_light
      psq=ce+ce*ce/(c*c)
c      if(.not.wrel) then
c        psq=ce+ce*ce/(c*c)
c      else
c        psq=ce
c      end if
c
      nl=lmax+1
      nl2=nl*nl
      kmax=2*lmax+1
      kmymax=2*nl2
      lmaxs=2*lmax
      lmmaxs=(lmaxs+1)*(lmaxs+1)
c
      mpar(0)=1
      do m=1,lmaxs
        mpar(m)=-mpar(m-1)
        mpar(-m)=mpar(m)
      end do
c
c      open(31,file=tripletout,status='unknown')
c
      tautest=.false.
      if(tautest) then
        write(6,*) "<locquant2> : I/O ie"
        write(6,*) ie
        write(6,*) "<locquant2> : I/O ce"
        write(6,*) ce
        write(6,*) "<locquant2> : I/O we"
        write(6,*) we
        write(6,*) "<locquant2> : I/O lmax"
        write(6,*) lmax
        write(6,*) "<locquant2> : I/O madmax"
        write(6,*) madmax
        write(6,*) "<locquant2> : I/O nintfc"
        write(6,*) nintfc
        write(6,*) "<locquant2> : I/O wrel"
        write(6,*) wrel
        write(6,*) "<locquant2> : I/O lms"
        write(6,*) lms
        write(6,*) "<locquant2> : I/O sxa"
        write(6,*) sxa
        write(6,*) "<locquant2> : I/O v0"
        write(6,*) v0
        write(6,*) "<locquant2> : I/O deltaa"
        write(6,*) deltaa
        write(6,*) "<locquant2> : I/O idpota"
        write(6,*) idpota
        write(6,*) "<locquant2> : I/O vra"
        write(6,*) vra
        write(6,*) "<locquant2> : I/O bra"
        write(6,*) bra
        write(6,*) "<locquant2> : I/O rba"
        write(6,*) rba
        write(6,*) "<locquant2> : I/O bopra"
        write(6,*) bopra
        write(6,*) "<locquant2> : I/O dx"
        write(6,*) dx
        write(6,*) "<locquant2> : I/O ns"
        write(6,*) ns
        write(6,*) "<locquant2> : I/O rs"
        write(6,*) rs
        write(6,*) "<locquant2> : I/O nimp"
        write(6,*) nimp
        write(6,*) "<locquant2> : I/O nposimp"
        write(6,*) nposimp
        write(6,*) "<locquant2> : I/O taua"
        write(6,*) taua
        write(6,*) "<locquant2> : I/O dosa"
        write(6,*) dosa
        write(6,*) "<locquant2> : I/O dosha"
        write(6,*) dosha
        write(6,*) "<locquant2> : I/O doseha"
        write(6,*) doseha
        write(6,*) "<locquant2> : I/O doshea"
        write(6,*) doshea
        write(6,*) "<locquant2> : I/O dosteha"
        write(6,*) dosteha
        write(6,*) "<locquant2> : I/O dosthea"
        write(6,*) dosthea
        write(6,*) "<locquant2> : I/O qvpa"
        write(6,*) qvpa
        write(6,*) "<locquant2> : I/O qvhpa"
        write(6,*) qvhpa
        write(6,*) "<locquant2> : I/O qva"
        write(6,*) qva
        write(6,*) "<locquant2> : I/O qvha"
        write(6,*) qvha
        write(6,*) "<locquant2> : I/O qveha"
        write(6,*) qveha
        write(6,*) "<locquant2> : I/O qvhea"
        write(6,*) qvhea
        write(6,*) "<locquant2> : I/O qvteha"
        write(6,*) qvteha
        write(6,*) "<locquant2> : I/O qvthea"
        write(6,*) qvthea
        write(6,*) "<locquant2> : I/O vmadid"
        write(6,*) vmadid
        write(6,*) "<locquant2> : I/O vmadiq"
        write(6,*) vmadiq
        write(6,*) "<locquant2> : I/O enba"
        write(6,*) enba
        write(6,*) "<locquant2> : I/O enbha"
        write(6,*) enbha
        write(6,*) "<locquant2> : I/O enbdiffa"
        write(6,*) enbdiffa
        write(6,*) "<locquant2> : I/O enbhdiffa"
        write(6,*) enbhdiffa
c       write(6,*) "<locquant2> : I/O denba"
c       write(6,*) denba
!        write(6,*) "<locquant2> : I/O enorba"
!        write(6,*) enorba
!        write(6,*) "<locquant2> : I/O qmoma"
!        write(6,*) qmoma
        write(6,*) "<locquant2> : I/O rhova"
        write(6,*) rhova
        write(6,*) "<locquant2> : I/O rhospa"
        write(6,*) rhospa
        write(6,*) "<locquant2> : I/O rhodspa"
        write(6,*) rhodspa
        write(6,*) "<locquant2> : I/O rhomaga"
        write(6,*) rhomaga
        write(6,*) "<locquant2> : I/O dosmaga"
        write(6,*) dosmaga
        write(6,*) "<locquant2> : I/O dosmagha"
        write(6,*) dosmagha
        write(6,*) "<locquant2> : I/O spin_magvpa"
        write(6,*) spin_magvpa
        write(6,*) "<locquant2> : I/O spin_magvhpa"
        write(6,*) spin_magvhpa
        write(6,*) "<locquant2> : I/O spin_magva"
        write(6,*) spin_magva
        write(6,*) "<locquant2> : I/O spin_magvha"
        write(6,*) spin_magvha
        write(6,*) "<locquant2> : I/O orb_magvpa"
        write(6,*) orb_magvpa
        write(6,*) "<locquant2> : I/O orb_magvhpa"
        write(6,*) orb_magvhpa
        write(6,*) "<locquant2> : I/O orb_magva"
        write(6,*) orb_magva
        write(6,*) "<locquant2> : I/O orb_magvha"
        write(6,*) orb_magvha
        write(6,*) "<locquant2> : I/O lliter"
        write(6,*) lliter
        write(6,*) "<locquant2> : I/O enbifc"
        write(6,*) enbifc
        write(6,*) "<locquant2> : I/O enbhifc"
        write(6,*) enbhifc
        write(6,*) "<locquant2> : I/O qvifc"
        write(6,*) qvifc
        write(6,*) "<locquant2> : I/O qvhifc"
        write(6,*) qvhifc
        write(6,*) "<locquant2> : I/O omifcint"
        write(6,*) omifcint
        write(6,*) "<locquant2> : I/O iwtau"
        write(6,*) iwtau
        write(6,*) "<locquant2> : I/O nwtau"
        write(6,*) nwtau
        write(6,*) "<locquant2> : I/O reg"
        write(6,*) reg
        write(6,*) "<locquant2> : I/O c_light"
        write(6,*) c_light
      end if
c
c
c
c
c ***************************************************
c * loop over layers to compute physical quantities *
c ***************************************************
c 
      qvifc=0.d0
      enbifc=0.d0
      do li=1,nintfc
         write(6,*) '<locquant> li=',li
!        cpalay=(1.d0-conc(li)).gt.tiny
c
c Compute scattering solutions 
c        --------------------------------------------------------
c        call cpu_time(t1)
c         call wafu(ce,psq,lmax,idpota(li),v0,vra(1,li),bra(1,li),
c     >             bopra(1,1,li),dx(li),ns(li),rs(li),
c     >             tm,gz,fz,gj,fj,nuz,indz,1,sxa(li))
        rtmp = rba(:,li) ! TODO : set initial value of rtmph
        rtmph(1) = rtmp(1) ! TODO : set initial value of rtmph
        rtmph(2) = -rtmp(2) ! TODO : set initial value of rtmph
        rtmph(3) = rtmp(3) ! TODO : set initial value of rtmph
c       ! TODO : set initial value for eh singlet, triplet
       call wafu(ce,c,sxa(li),lmax,idpota(li),v0,E_Fermi,
     >            singratimp(li),uratimp(li),dratimp(li),
     >            deltaa(1,li),vra(1,li),bra(1,li),
     >            rtmp,dx(li),rs(li),ns(li),tm,
     >            gz,fz,gj,fj,glz,flz,glj,flj,1)
c        call cpu_time(t2)
c        write(6,*) 'wafu:',t2-t1
c        --------------------------------------------------------
        call gf(lmax,reg,dx(li),ns(li),rs(li),
     >          gz,fz,gj,fj,glz,flz,glj,flj,taua(1,1,li),
     >          ggrmat,fgrmat)
c Density of multipole moments
         sig=1.d0
         i=0
c         write(31,*) 'Geh',ie,ce,li
         call tripletcalc(lmax,lms,rgacoeff(:,:,1),rgacoeff(:,:,1),sig,
     >                ggrmat(1:kmymax,kmymax+1:2*kmymax,0),
     >                fgrmat(1:kmymax,kmymax+1:2*kmymax,0),rtmp,
     >                tripletmat(1,1,1,li))
         do lam=0,lmaxs
         do mu=-lam,lam
           i=i+1
c        ---------------------------------------------------------- e-e
          call moment(lmax,lms,rgacoeff(:,:,i),rgacoeff(:,:,i),sig,
     >                ggrmat(1:kmymax,1:kmymax,lam),
     >                fgrmat(1:kmymax,1:kmymax,lam),
     >                zqmom(1:kmymax,i),rtmp)
c         --------------------------------------------------------- h-h
          call moment(lmax,lms,rgacoeff(:,:,i),rgacoeff(:,:,i),sig,
     >                ggrmat(kmymax+1:2*kmymax,kmymax+1:2*kmymax,lam),
     >                fgrmat(kmymax+1:2*kmymax,kmymax+1:2*kmymax,lam),
     >                zqmom(kmymax+1:2*kmymax,i),rtmph)
c
      tautest=.false.
      if(tautest) then
      write(6,*) " <locquant> : moment output lmax"
      write(6,*) lmax
      write(6,*) " <locquant> : moment output lms"
      write(6,*) lms
      write(6,*) " <locquant> : moment output rgacoeff(:,:,i)"
      write(6,*) rgacoeff(:,:,i)
      write(6,*) " <locquant> : moment output rgacoeff(:,:,i)"
      write(6,*) rgacoeff(:,:,i)
      write(6,*) " <locquant> : moment output sig"
      write(6,*) sig
      write(6,*) " <locquant> : moment output
     > ggrmat(kmymax+1:2*kmymax,kmymax+1:2*kmymax,lam)"
      write(6,*) ggrmat(kmymax+1:2*kmymax,kmymax+1:2*kmymax,lam)
      write(6,*) " <locquant> : moment output
     > fgrmat(kmymax+1:2*kmymax,kmymax+1:2*kmymax,lam)"
      write(6,*) fgrmat(kmymax+1:2*kmymax,kmymax+1:2*kmymax,lam)
      write(6,*) " <locquant> : moment output
     > zqmom(kmymax+1:2*kmymax,i)"
      write(6,*) zqmom(kmymax+1:2*kmymax,i)
      write(6,*) " <locquant> : moment output rtmph"
      write(6,*) rtmph
      end if
c---> Transformation applied -> in gf <- for the offdiag blocks to obtain the anomalous DOS 
c
c         --------------------------------------------------------- ehS           
          call pairing(deltakmy)
          call trafogeh(ggrmat,fgrmat,deltakmy,lmax,tggrmat,tfgrmat) 
          call moment(lmax,lms,rgacoeff(:,:,i),rgacoeff(:,:,i),sig,
     >                tggrmat(1:kmymax,1:kmymax,lam),
     >                tfgrmat(1:kmymax,1:kmymax,lam),
     >                zscmom(1:kmymax,i),rtmp)
c         --------------------------------------------------- ehT0 -> dz
          call t0pairing(deltakmy)
          call trafogeh(ggrmat,fgrmat,deltakmy,lmax,tggrmat,tfgrmat) 
          call moment(lmax,lms,rgacoeff(:,:,i),rgacoeff(:,:,i),sig,
     >                tggrmat(1:kmymax,1:kmymax,lam),
     >                tfgrmat(1:kmymax,1:kmymax,lam),
     >                zscmom(kmymax+1:2*kmymax,i),rtmp)
c         ---------------------------------------------------------
c --> Non-Unitary triplet components in the anomalous part
c         --------------------------------------------------- ehTU -> dx           
          call tuppairing(deltakmy)
          call trafogeh(ggrmat,fgrmat,deltakmy,lmax,tggrmat,tfgrmat) 
          call moment(lmax,lms,rgacoeff(:,:,i),rgacoeff(:,:,i),sig,
     >                tggrmat(1:kmymax,1:kmymax,lam),
     >                tfgrmat(1:kmymax,1:kmymax,lam),
     >                zsctmom(1:kmymax,i),rtmp)
c         --------------------------------------------------- ehTD -> dy
          call tdownpairing(deltakmy)
          call trafogeh(ggrmat,fgrmat,deltakmy,lmax,tggrmat,tfgrmat) 
          call moment(lmax,lms,rgacoeff(:,:,i),rgacoeff(:,:,i),sig,
     >                tggrmat(1:kmymax,1:kmymax,lam),
     >                tfgrmat(1:kmymax,1:kmymax,lam),
     >                zsctmom(kmymax+1:2*kmymax,i),rtmp)
c         --------------------------------------------------------- 

         end do
         end do
! commented by bnyari, we dont use the zdos quantities instead we use
! zqmom(:,1)
!         do kmy=1,2*kmymax
!            zdos(kmy)=zqmom(kmy,1)
!            zcsdos(kmy) = zscmom(kmy,1)
!            zsctdos(kmy) = zsctmom(kmy,1)
!         end do
c print out tau matrix and vawe function on unit 16
         if((li.eq.nwtau).and.(iwtau.eq.1))then
          write(16) we,ce
          write(16) rs(li),dx(li),ns(li)
c  write tau matrix
          write(16) ((taua(i,j,li),i=1,kmymax),j=1,kmymax)          
c  write big component of regular solution
          do ik = 1, kmymax
           do i = 1,2
            write(16) (gz(j,i,ik),j=1,ns(i))
           end do
          end do
c  write small component of regular solution
          do ik = 1, kmymax
           do i = 1,2
            write(16) (fz(j,i,ik),j=1,ns(i))
           end do
          end do
c  write big component of irregular solution
          do ik = 1, kmymax
           do i = 1,2
            write(16) (gj(j,i,ik),j=1,ns(i))
           end do
          end do
c  write small component of irregular solution
          do ik = 1, kmymax
           do i = 1,2
            write(16) (fj(j,i,ik),j=1,ns(i))
           end do
          end do
         endif
c
c Magnetic moments 
        scoeff(1:kmymaxp,1:kmymaxp,1) = sxcoeff(1:kmymaxp,1:kmymaxp)
        scoeff(1:kmymaxp,1:kmymaxp,2) = sycoeff(1:kmymaxp,1:kmymaxp)
        scoeff(1:kmymaxp,1:kmymaxp,3) = szcoeff(1:kmymaxp,1:kmymaxp)
        sbcoeff(1:kmymaxp,1:kmymaxp,1) = sxbcoeff(1:kmymaxp,1:kmymaxp)
        sbcoeff(1:kmymaxp,1:kmymaxp,2) = sybcoeff(1:kmymaxp,1:kmymaxp)
        sbcoeff(1:kmymaxp,1:kmymaxp,3) = szbcoeff(1:kmymaxp,1:kmymaxp)
        lcoeff(1:kmymaxp,1:kmymaxp,1) = lxcoeff(1:kmymaxp,1:kmymaxp)
        lcoeff(1:kmymaxp,1:kmymaxp,2) = lycoeff(1:kmymaxp,1:kmymaxp)
        lcoeff(1:kmymaxp,1:kmymaxp,3) = lzcoeff(1:kmymaxp,1:kmymaxp)
        lbcoeff(1:kmymaxp,1:kmymaxp,1) = lxbcoeff(1:kmymaxp,1:kmymaxp)
        lbcoeff(1:kmymaxp,1:kmymaxp,2) = lybcoeff(1:kmymaxp,1:kmymaxp)
        lbcoeff(1:kmymaxp,1:kmymaxp,3) = lzbcoeff(1:kmymaxp,1:kmymaxp)
        sig=-1.d0
        do i=1,3
c       ------------------------------------------------- e-e 
        call moment(lmax,lms,scoeff(1:kmymaxp,1:kmymaxp,i),
     >              sbcoeff(1:kmymaxp,1:kmymaxp,i),sig,
     >              ggrmat(1:kmymax,1:kmymax,0),
     >              fgrmat(1:kmymax,1:kmymax,0),zdmag(1:kmymax,i),rtmp)
        call moment(lmax,.true.,lcoeff(1:kmymaxp,1:kmymaxp,i),
     >              lbcoeff(1:kmymaxp,1:kmymaxp,i),sig,
     >              ggrmat(1:kmymax,1:kmymax,0),
     >              fgrmat(1:kmymax,1:kmymax,0),zdmor(1:kmymax,i),rtmp)
c        ------------------------------------------------- h-h
        call moment(lmax,lms,scoeff(1:kmymaxp,1:kmymaxp,i),
     >              sbcoeff(1:kmymaxp,1:kmymaxp,i),sig,
     >              ggrmat(kmymax+1:2*kmymax,kmymax+1:2*kmymax,0),
     >              fgrmat(kmymax+1:2*kmymax,kmymax+1:2*kmymax,0),
     >              zdmag(kmymax+1:2*kmymax,i),rtmp)
        call moment(lmax,.true.,lcoeff(1:kmymaxp,1:kmymaxp,i),
     >              lbcoeff(1:kmymaxp,1:kmymaxp,i),sig,
     >              ggrmat(kmymax+1:2*kmymax,kmymax+1:2*kmymax,0),
     >              fgrmat(kmymax+1:2*kmymax,kmymax+1:2*kmymax,0),
     >              zdmor(kmymax+1:2*kmymax,i),rtmp)
c        -------------------------------------------------
        end do

c Magnetic density of states - old version 
c        -------------------------------------------------
c         call magnet(lmax,rs(li),dx(li),ns(li),
c     >               gz,fz,gj,fj,nuz,indz,taua(1,1,li),
c     >               sxcoeff,sxbcoeff,zdmag(1,1),lms,1)
c         call magnet(lmax,rs(li),dx(li),ns(li),
c     >               gz,fz,gj,fj,nuz,indz,taua(1,1,li),
c     >               sycoeff,sybcoeff,zdmag(1,2),lms,1)
c         call magnet(lmax,rs(li),dx(li),ns(li),
c     >               gz,fz,gj,fj,nuz,indz,taua(1,1,li),
c     >               szcoeff,szbcoeff,zdmag(1,3),lms,1)
c         call magnet(lmax,rs(li),dx(li),ns(li),
c     >               gz,fz,gj,fj,nuz,indz,taua(1,1,li),
c     >               lxcoeff,lxbcoeff,zdmor(1,1),.true.,1)
c         call magnet(lmax,rs(li),dx(li),ns(li),
c     >               gz,fz,gj,fj,nuz,indz,taua(1,1,li),
c     >               lycoeff,lybcoeff,zdmor(1,2),.true.,1)
c         call magnet(lmax,rs(li),dx(li),ns(li),
c     >               gz,fz,gj,fj,nuz,indz,taua(1,1,li),
c     >               lzcoeff,lzbcoeff,zdmor(1,3),.true.,1)
c        -------------------------------------------------
        
        do k=1,kmymax
          dosa(k,li)=dimag(zqmom(k,1))                  ! electron part
          dosha(k,li)=dimag(zqmom(k+kmymax,1))          ! hole part
          doseha(k,li)=dimag(zscmom(k,1))               ! electron-hole part 
          doshea(k,li)=dimag(zscmom(k+kmymax,1))        ! triplet part 
          dosteha(k,li)=dimag(zsctmom(k,1))               ! triplet up part 
          dosthea(k,li)=dimag(zsctmom(k+kmymax,1))        ! triplet down part 
          dosmaga(k,li)=dimag(rtmp(1)*zdmag(k,1)     ! electron part of magnetic dos
     >                       +rtmp(2)*zdmag(k,2)
     >                       +rtmp(3)*zdmag(k,3))  ! rtmp either rz or rba
          dosmagha(k,li)=dimag(rtmph(1)*zdmag(k+kmymax,1)    ! hole part of magnetic dos
     >                       +rtmph(2)*zdmag(k+kmymax,2)
     >                       +rtmph(3)*zdmag(k+kmymax,3)) ! rtmph for holes
          tautest=.false.
          if(itest.ge.2.or.tautest) then
            write(6,*) 'DOS electron kmy layeri'
            write(6,*) dosa(k,li),k,li
            write(6,*) 'DOS hole kmy layeri'
            write(6,*) dosha(k,li),k,li
c     >    write(6,'(2d13.5,5x,2d13.5)') dosa(k,li),dosmaga(k,li)
          end if
        end do
c
c         do k=1,kmymax
c           dosa(k,li)=dimag(zdos(k))
c           dosmaga(k,li)=dimag(zdmag(k,3))
c           if(itest.ge.2)
c     >     write(6,'(2d13.5,5x,2d13.5)') dosa(k,li),dosmaga(k,li)
c         end do
c
c Radial distribution of charge and magnetic density
c        ------------------------------------------------
c        call cpu_time(t1)
c         call dens(lmax,rs(li),dx(li),ns(li),gz,fz,gj,fj,
c     >             nuz,indz,taua(1,1,li),bopra(1,1,li),
c     >             zrho,zrhosp,zrhodsp,zmagrho,zenorb)
         call dens(lmax,rs(li),dx(li),ns(li),rtmp,
     >             gz,fz,gj,fj,glz,flz,glj,flj,
     >             taua(1,1,li),zrhoeh,zrho)
c        call cpu_time(t2)
c        write(6,*) 'dens:',t2-t1
c        ------------------------------------------------
c add contribution to contour integral
c
         qvdiffa(li)=0.d0
         qvhdiffa(li)=0.d0
!        qvehdiffa(li)=0.d0 ! TODO : commented on 27/8/2020 Laszloffy
!        qvhediffa(li)=0.d0 ! TODO : commented on 27/8/2020 Laszloffy
!        qvtehdiffa(li)=0.d0 ! TODO : commented on 27/8/2020 Laszloffy
!        qvthediffa(li)=0.d0 ! TODO : commented on 27/8/2020 Laszloffy
c
         enbdiffa(li)=0.d0
         enbhdiffa(li)=0.d0
c            
c TODO check init zero before energy loop! 
c         qv=0.0d0
c         qvh=0.0d0
c         qveh=(0.0d0,0.0d0)
c         qvhe=(0.0d0,0.0d0)
c         qvteh=(0.0d0,0.0d0)
c         qvthe=(0.0d0,0.0d0)
c
c         qvdiffa(li)=0.d0
c         enbdiffa(li)=0.d0
c         qv=0.0d0
         do k=1,kmymax
c            qv=qv+dimag(we*zdos(k))
c            qvdiffa(li)=qvdiffa(li)+dimag(we*zdos(k))
c            enbdiffa(li)=enbdiffa(li)+dimag(we*ce*zdos(k))
c            qvpa(k,li)=qvpa(k,li)+dimag(we*zdos(k))
c            qva(li)=qva(li)+dimag(we*zdos(k))
c            nqva(li,ie)=nqva(li,ie)+zdos(k)
c            qv=qv+dimag(we*zdos(k))               ! electron charge
c            qvh=qvh+dimag(we*zdos(k+kmymax))      ! hole charge
c            qveh=qveh+we*zscdos(k)                ! e-h singlet charge
c            qvhe=qvhe+we*zscdos(k+kmymax)  ! t0 charge
c            qvteh=qvteh+we*zsctdos(k)         ! t up charge
c            qvthe=qvthe+we*zsctdos(k+kmymax)  ! t down charge
c
c            qvdiffa(li)=qvdiffa(li)+dimag(we*zdos(k))
c            qvhdiffa(li)=qvhdiffa(li)+dimag(we*zdos(k+kmymax))
!           qvehdiffa(li)=qvehdiffa(li)+dimag(we*zscdos(k))           ! TODO : commented on 27/8/2020 Laszloffy
!           qvhediffa(li)=qvhediffa(li)+dimag(we*zscdos(k+kmymax)) ! TODO : commented on 27/8/2020 Laszloffy
!           qvtehdiffa(li)=qvtehdiffa(li)+dimag(we*zsctdos(k)) ! TODO : commented on 27/8/2020 Laszloffy
!           qvthediffa(li)=qvthediffa(li)+dimag(we*zsctdos(k+kmymax)) ! TODO : commented on 27/8/2020 Laszloffy
c
c            enbdiffa(li)=enbdiffa(li)+dimag(we*ce*zdos(k))
c            enbhdiffa(li)=enbhdiffa(li)+dimag(we*ce*zdos(k+kmymax))
c            enbehdiffa(li)=enbehdiffa(li)+dimag(we*ce*zscdos(k))
c            enbhediffa(li)=enbhediffa(li)+dimag(we*ce*zscdos(k+kmymax))
c            enbtehdiffa(li)=enbtehdiffa(li)+dimag(we*ce*zsctdos(k))
c            enbthediffa(li)=enbthediffa(li)+dimag(we*ce*zsctdos(k+kmymax))
c
            qvpa(k,li)=qvpa(k,li)+dimag(we*zqmom(k,1))
            qvhpa(k,li)=qvhpa(k,li)+dimag(we*zqmom(k+kmymax,1))
!           qvehpa(k,li)=qvehpa(k,li)+dimag(we*zscdos(k)) ! TODO : commented on 27/8/2020 Laszloffy
!           qvhepa(k,li)=qvhepa(k,li)+dimag(we*zscdos(k+kmymax)) ! TODO : commented on 27/8/2020 Laszloffy
!           qvtehpa(k,li)=qvtehpa(k,li)+dimag(we*zsctdos(k)) ! TODO : commented on 27/8/2020 Laszloffy
!           qvthepa(k,li)=qvthepa(k,li)+dimag(we*zsctdos(k+kmymax)) ! TODO : commented on 27/8/2020 Laszloffy
c
            qva(li)=qva(li)+dimag(we*zqmom(k,1))
            qvha(li)=qvha(li)+dimag(we*zqmom(k+kmymax,1))
            qveha(li)=qveha(li)+dimag(we*zscmom(k,1))
            qvhea(li)=qvhea(li)+dimag(we*zscmom(k+kmymax,1))
            qvteha(li)=qvteha(li)+dimag(we*zsctmom(k,1))
            qvthea(li)=qvthea(li)+dimag(we*zsctmom(k+kmymax,1))
c
            do i=1,3
              spin_magvpa(k,li,i)=spin_magvpa(k,li,i)+
     >                            dimag(we*zdmag(k,i))
              spin_magvhpa(k,li,i)=spin_magvhpa(k,li,i)+
     >                            dimag(we*zdmag(k+kmymax,i))
c
              spin_magva(li,i)=spin_magva(li,i)+
     >                         dimag(we*zdmag(k,i))
              spin_magvha(li,i)=spin_magvha(li,i)+
     >                         dimag(we*zdmag(k+kmymax,i))
c
              orb_magvpa(k,li,i)=orb_magvpa(k,li,i)+
     >                           dimag(we*zdmor(k,i))
              orb_magvhpa(k,li,i)=orb_magvhpa(k,li,i)+
     >                           dimag(we*zdmor(k+kmymax,i))
c
              orb_magva(li,i)=orb_magva(li,i)+
     >                        dimag(we*zdmor(k,i))
              orb_magvha(li,i)=orb_magvha(li,i)+
     >                        dimag(we*zdmor(k+kmymax,i))
            end do
            enba(li)=enba(li)+dimag(we*ce*zqmom(k,1))
            enbha(li)=enbha(li)+dimag(we*ce*zqmom(k+kmymax,1))
c           nenba(li,ie)=nenba(li,ie)+ce*zdos(k) ! commented on 01/10/2020 Laszloffy
!           nenbha(li,ie)=nenbha(li,ie)+ce*zdos(k+kmymax) ! TODO : commented on 27/8/2020 Laszloffy
         end do
!         enorba(li)=enorba(li)+dimag(we*zenorb)
c
c        call enblloyd(enbll,we,gtaua(1,1,li),ptminva(1,1,li),
c    >                 tminv(1,1,li),alphaintkkr(0,li),lmax,lliter)
c        enba(li)=enba(li)+enbll
c        write(6,'(2d15.6)') enbll,enba(li)
c        
c todo check lloyd formula! 
!        call denblloyd(denbll,we,lmax, ! TODO : commented on 27/8/2020 Laszloffy
!    >                  gtaua(1,1,li),ptminva(1,1,li), ! TODO : commented on 27/8/2020 Laszloffy
!    >                  dmata(1,1,li),dmatpa(1,1,li), ! TODO : commented on 27/8/2020 Laszloffy
!    >                  ddtha(1,1,li),ddthpa(1,1,li)) ! TODO : commented on 27/8/2020 Laszloffy
c    >                  ddpha(1,1,li),ddphpa(1,1,li))
!        denba(li)=denba(li)+denbll ! TODO : commented on 27/8/2020 Laszloffy
c
         do irad=1,ns(li)
!           rhov(irad)=dimag(we*zrho(irad))
            rhova(irad,li)=rhova(irad,li)+dimag(we*zrho(irad))
c           rhoveha(irad,li)=rhoveha(irad,li)+we*zrhoeh(irad) ! TODO : commented on 27/8/2020 Laszloffy
            rhospa(irad,1,li)=rhospa(irad,1,li)+dimag(we*zrhosp(irad,1))
            rhospa(irad,2,li)=rhospa(irad,2,li)+dimag(we*zrhosp(irad,2))
            rhodspa(irad,1,li)=rhodspa(irad,1,li)+
     >                         dimag(we*zrhodsp(irad,1))
            rhodspa(irad,2,li)=rhodspa(irad,2,li)+
     >                         dimag(we*zrhodsp(irad,2))
            rhomaga(irad,li)=rhomaga(irad,li)+dimag(we*zmagrho(irad))
         end do
c
         do lam=0,lmaxs
         do mu=-lam,lam
            i=lam*(lam+1)+mu+1
            ii=lam*(lam+1)-mu+1
            qmom(i)=(0.0d0,0.0d0)
            do k=1,kmymax
              qmom(i)=qmom(i)+cih*(we*zqmom(k,i)-
     >                             mpar(mu)*dconjg(we*zqmom(k,ii)))
              qmomh(i)=qmomh(i)+cih*(we*zqmom(k+kmymax,i)-
     >                         mpar(mu)*dconjg(we*zqmom(k+kmymax,ii)))
            end do
         end do
         end do
c
c Rotate moments of charge density to the global frame of reference
c and add to the integral
         qmomeh=(0.0d0,0.0d0)
         qmomhe=(0.0d0,0.0d0)
         qmomteh=(0.0d0,0.0d0)
         qmomthe=(0.0d0,0.0d0)
         do k=1,kmymax
              qmomeh=qmomeh+dimag(we*zscdos(k))
              qmomhe=qmomhe+dimag(we*zscdos(k+kmymax))
              qmomteh=qmomteh+dimag(we*zsctdos(k))        
              qmomthe=qmomthe+dimag(we*zsctdos(k+kmymax))
         end do
c
c Add moments of charge density to the integral
c
!          qmoma(1:lmmaxs,li) = qmoma(1:lmmaxs,li) + qmom(1:lmmaxs)
          qmomha(1:lmmaxs,li) = qmomha(1:lmmaxs,li) + qmomh(1:lmmaxs)
          qmomeha(li) = qmomeha(li) + qmomeh
          qmomhea(li) = qmomhea(li) + qmomhe
          qmomteha(li) = qmomteha(li) + qmomteh
          qmomthea(li) = qmomthea(li) + qmomthe          
c TODO bnyari: how to calculate qmomrot?
c         do i=1,lmmaxs
c           qmomrot(i)=(0.d0,0.d0)
c           do j=1,lmmaxs
c             qmoma(i,li)=qmoma(i,li)+qmom(j)*rmatpa(j,i,li)
c             qmomrot(i)=qmomrot(i)+qmom(j)*rmatpa(j,i,li)
c           end do
c         end do
c
c ======================================================================
c *******************************************************
c * Madelung potential due to the impurities wihtout 00 *
c *******************************************************
!        call madcoimp(madmax,li,nimp,nposimp,mlmij) ! TODO : commented on 27/8/2020 Laszloffy
c        do jimp=1,nintfc
c          if(jimp.ne.li) then
c  from dipole momentum
c
c            temp=(0.d0,0.d0) !DEBUG
c            do ii=2,4
c              vmadid(jimp)=vmadid(jimp)+
c    >                      dreal(mlmij(ii,jimp)*qmomrot(ii))
c              temp=temp+mlmij(ii,jimp)*qmomrot(ii)
c            end do
c            if(dimag(temp).gt.1.0d-12) then
c               write(6,*) '<locquant1>: *********************'
c               write(6,*) '<locquant1>: IMAGINARY PART OF THE'
c               write(6,*) '<locquant1>: MADELUNG DIP .NE. 0  '
c               write(6,*) '<locquant1>: temp=',jimp,temp !DEBUG
c               write(6,*) '<locquant1>: STOP DEBUG'
c               write(6,*) '<locquant1>: *********************'
c            end if
c
c  from quadrupol momentum
c
c            temp=(0.d0,0.d0) !DEBUG
c            do ii=5,9
c              vmadiq(jimp)=vmadiq(jimp)+
c    >                      dreal(mlmij(ii,jimp)*qmomrot(ii))
c              temp=temp+mlmij(ii,jimp)*qmomrot(ii)
c            end do
c            if(dimag(temp).gt.1.0d-12) then
c               write(6,*) '<locquant1>: *********************'
c               write(6,*) '<locquant1>: IMAGINARY PART OF THE'
c               write(6,*) '<locquant1>: MADELUNG QUAD .NE. 0 '
c               write(6,*) '<locquant1>: temp=',jimp,temp !DEBUG
c               write(6,*) '<locquant1>: STOP DEBUG'
c               write(6,*) '<locquant1>: *********************'
c            end if
c
c          end if
c        end do
c ======================================================================
c Check charges
c        x=dlog(rs(li))-(ns(li)-1)*dx(li)
c        do inad=1,ns(li)
c          rad(irad)=dexp(x)
c          x=x+dx(li)
c        end do
c        qvp=rsimp(rhov,rad,ns(li),dx(li))
c        write(6,'(''Layer'',i5,''    qv,qvp:'',2d20.10)') li,qv,qvp
c
         if(itest.ge.3) then
           write(6,'('' Contribution to contour integral so far'')')
           write(6,'(''Q    '',d13.5)') qva(li)
           do i=1,lmmaxs
             write(6,'(''Q'',i2,2x,2d13.5)') i,qmoma(i,li)
           end do
           write(6,'(''Eband'',d13.5)') enba(li)
!           write(6,'(''Eorb '',d13.5)')  enorba(li)
           write(6,'(''Sx,Lx'',d13.5,5x,d13.5)') 
     >               spin_magva(li,1),orb_magva(li,1)
           write(6,'(''Sy,Ly'',d13.5,5x,d13.5)') 
     >               spin_magva(li,2),orb_magva(li,2)
           write(6,'(''Sz,Lz'',d13.5,5x,d13.5)') 
     >               spin_magva(li,3),orb_magva(li,3)
           if(itest.gt.3) then
             write(6,'(''Loop nrad - rho,rhomag'')')
             do irad=1,ns(li)
              write(6,'(d13.5,5x,d13.5)') 
     >              rhova(irad,li),rhomaga(irad,li)
              write(6,'(d13.5,5x,d13.5)') 
     >              rhospa(irad,1,li),rhospa(irad,2,li)
              write(6,'(d13.5,5x,d13.5)') 
     >              rhodspa(irad,1,li),rhodspa(irad,2,li)
             enddo
           endif
         endif
c
c
         enbifc=enbifc+enbdiffa(li)
         enbhifc=enbhifc+enbhdiffa(li)
c        enbifc=enbifc+conc(li)*enbdiffa(li)
c    >                +(1.0d0-conc(li))*enbdiffb(li)
         qvifc=qvifc+qvdiffa(li)
         qvhifc=qvifc+qvhdiffa(li)
c        qvifc=qvifc+conc(li)*qvdiffa(li)
c    >              +(1.0d0-conc(li))*qvdiffb(li)
c
      end do
c ************************
c * end loop over layers *
c ************************
c     if(linbw) then
c       omifc=enbifc-efermi*qvifc
c       enbifcint=enbifcint+enbifc
c       qvifcint=qvifcint+qvifc
c       omifcint=omifcint+omifc
c       write(6,'('' Enb'',i3,2f10.5,3d17.8)') 
c    >  ie,ce,enbifc,qvifc,omifc
c     end if
c
      return
      end
