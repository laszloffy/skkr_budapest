c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine locquant(
c     ===================
     > ie,ce,we,dos,lmax,nintfc,wrel,lms,sxa,sxb,conc,v0,
     > E_Fermi,deltaa,deltab,
c
     > idpota,vra,bra,rba,idpotb,vrb,brb,rbb,dx,ns,rs,
c     > dmata,dmatb,dmatpa,dmatpb,
c     > ddpha,ddtha,d2dpha,d2dtha,d2dthpha,
c     > ddphpa,ddthpa,d2dphpa,d2dthpa,d2dthphpa,
c     > ddphb,ddthb,d2dphb,d2dthb,d2dthphb,
c     > ddphpb,ddthpb,d2dphpb,d2dthpb,d2dthphpb,
c     > rmata,rmatpa,rmatb,rmatpb,
     > taua,taub,gtaua,gtaub,
c     > ptminva,ptminvb,tminv,
     > dosa,dosha,doseha,doshea,dosteha,dosthea,
     > qvpa,qva,qvha,qveha,qvhea,qvteha,qvthea,
     > qvdiffa,enba,enbdiffa,qmoma,
c     > ddphenba,ddthenba,d2dphenba,d2dthenba,d2dthphenba,
     > rhova,rhoveha,rhovhea,rhospa,rhodspa,rhomaga,
     > dosmaga,dosmagha,
     > spin_magvpa,spin_magva,orb_magvpa,orb_magva,
     > dosb,doshb,dosehb,dosheb,dostehb,dostheb,
     > qvpb,qvb,qvhb,qvehb,qvheb,qvtehb,qvtheb,
     > qvdiffb,enbb,enbdiffb,qmomb,
c     > ddphenbb,ddthenbb,d2dphenbb,d2dthenbb,d2dthphenbb,
     > rhovb,rhovehb,rhovheb,rhospb,rhodspb,rhomagb,
     > dosmagb,dosmaghb,
     > spin_magvpb,spin_magvb,orb_magvpb,orb_magvb,
c
     > lliter,linbw,efermi,enbifcint,qvifcint,omifcint,reg,c_light)
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      logical wrel,lms,cpalay,linbw,dos,reg
c
      character*10 idpota(mintfc),idpotb(mintfc)
c
      dimension mpar(-lsup:lsup)
c
      dimension vra(nrad,mintfc),bra(nrad,mintfc)
      dimension vrb(nrad,mintfc),brb(nrad,mintfc)
      dimension rs(mintfc),dx(mintfc),ns(mintfc)
      dimension rba(3,mintfc),rbb(3,mintfc)
      complex*16 deltaa(nrad,mintfc),deltab(nrad,mintfc)
c
      dimension enba(mintfc),enbb(mintfc)
      dimension enbha(mintfc),enbhb(mintfc)
c      dimension ddphenba(mintfc),ddphenbb(mintfc)
c      dimension ddthenba(mintfc),ddthenbb(mintfc)
c      dimension d2dphenba(mintfc),d2dphenbb(mintfc)
c      dimension d2dthenba(mintfc),d2dthenbb(mintfc)
c      dimension d2dthphenba(mintfc),d2dthphenbb(mintfc)
c
      dimension dosa(kmymaxp,mintfc),dosb(kmymaxp,mintfc)
      dimension dosha(kmymaxp,mintfc),doshb(kmymaxp,mintfc)
      dimension doseha(kmymaxp,mintfc),dosteha(kmymaxp,mintfc)
      dimension doshea(kmymaxp,mintfc),dosthea(kmymaxp,mintfc)
      dimension dosehb(kmymaxp,mintfc),dostehb(kmymaxp,mintfc)
      dimension dosheb(kmymaxp,mintfc),dostheb(kmymaxp,mintfc)
      dimension dosmaga(kmymaxp,mintfc),dosmagb(kmymaxp,mintfc)
      dimension dosmagha(kmymaxp,mintfc),dosmaghb(kmymaxp,mintfc)
      dimension qva(mintfc),qvb(mintfc)
      dimension qvpa(kmymaxp,mintfc),qvpb(kmymaxp,mintfc)
      dimension qvha(mintfc),qvhb(mintfc)
      dimension qvhpa(kmymaxp,mintfc),qvhpb(kmymaxp,mintfc)
      dimension spin_magvpa(kmymaxp,mintfc,3)
      dimension spin_magvpb(kmymaxp,mintfc,3)
      dimension spin_magva(mintfc,3),spin_magvb(mintfc,3)
      dimension orb_magvpa(kmymaxp,mintfc,3)
      dimension orb_magvpb(kmymaxp,mintfc,3)
      dimension orb_magva(mintfc,3),orb_magvb(mintfc,3)
      dimension spin_magvhpa(kmymaxp,mintfc,3)
      dimension spin_magvhpb(kmymaxp,mintfc,3)
      dimension spin_magvha(mintfc,3),spin_magvhb(mintfc,3)
      dimension orb_magvhpa(kmymaxp,mintfc,3)
      dimension orb_magvhpb(kmymaxp,mintfc,3)
      dimension orb_magvha(mintfc,3)
      dimension orb_magvhb(mintfc,3)
c
      dimension enbdiffa(mintfc),enbdiffb(mintfc)
      dimension enbhdiffa(mintfc),enbhdiffb(mintfc)
      dimension qvdiffa(mintfc),qvdiffb(mintfc)
      dimension qvhdiffa(mintfc),qvhdiffb(mintfc)
c
      dimension rhova(nrad,mintfc),rhovb(nrad,mintfc)
      dimension rhospa(nrad,2,mintfc),rhospb(nrad,2,mintfc)
      dimension rhodspa(nrad,2,mintfc),rhodspb(nrad,2,mintfc)
      dimension rhomaga(nrad,mintfc),rhomagb(nrad,mintfc)
c      dimension rhov(nrad)
c
      dimension rad(nrad)
c
      dimension conc(mintfc),sxa(mintfc),sxb(mintfc)
c
c      complex*16 qveha,qvhea,qvteha,qvthea
c      complex*16 qvehb,qvheb,qvtehb,qvtheb
c
c      complex*16 rhoveh(nrad)
c      complex*16 rhovhe(nrad)
      complex*16 rhoveha(nrad,mintfc)
      complex*16 rhovhea(nrad,mintfc)
      complex*16 rhovehb(nrad,mintfc)
      complex*16 rhovheb(nrad,mintfc)
      complex*16 qveha(mintfc)
      complex*16 qvehb(mintfc)
      complex*16 qvehpa(kmymaxp,mintfc)
      complex*16 qvehpb(kmymaxp,mintfc)
      complex*16 qvteha(mintfc)
      complex*16 qvtehb(mintfc)
      complex*16 qvtehpa(kmymaxp,mintfc)
      complex*16 qvtehpb(kmymaxp,mintfc)
      complex*16 qvhea(mintfc)
      complex*16 qvheb(mintfc)
      complex*16 qvhepa(kmymaxp,mintfc)
      complex*16 qvhepb(kmymaxp,mintfc)
      complex*16 qvthea(mintfc)
      complex*16 qvtheb(mintfc)
      complex*16 qvthepa(kmymaxp,mintfc)
      complex*16 qvthepb(kmymaxp,mintfc)
      complex*16 qvehdiffa(mintfc)
      complex*16 qvehdiffb(mintfc)
      complex*16 qvhediffa(mintfc)
      complex*16 qvhediffb(mintfc)
      complex*16 qvtehdiffa(mintfc)
      complex*16 qvtehdiffb(mintfc)
      complex*16 qvthediffa(mintfc)
      complex*16 qvthediffb(mintfc)
c
      complex*16 ce,we,psq,cih
c
c      complex*16 dmata(kmymaxp,kmymaxp,mintfc)
c      complex*16 dmatpa(kmymaxp,kmymaxp,mintfc)
c      complex*16 dmatb(kmymaxp,kmymaxp,mintfc)
c      complex*16 dmatpb(kmymaxp,kmymaxp,mintfc)
c
c      complex*16 ddpha(kmymaxp,kmymaxp,mintfc)
c      complex*16 ddphpa(kmymaxp,kmymaxp,mintfc)
c      complex*16 ddphb(kmymaxp,kmymaxp,mintfc)
c      complex*16 ddphpb(kmymaxp,kmymaxp,mintfc)
c      complex*16 ddtha(kmymaxp,kmymaxp,mintfc)
c      complex*16 ddthpa(kmymaxp,kmymaxp,mintfc)
c      complex*16 ddthb(kmymaxp,kmymaxp,mintfc)
c      complex*16 ddthpb(kmymaxp,kmymaxp,mintfc)
c      complex*16 d2dtha(kmymaxp,kmymaxp,mintfc)
c      complex*16 d2dthpa(kmymaxp,kmymaxp,mintfc)
c      complex*16 d2dthb(kmymaxp,kmymaxp,mintfc)
c      complex*16 d2dthpb(kmymaxp,kmymaxp,mintfc)
c      complex*16 d2dpha(kmymaxp,kmymaxp,mintfc)
c      complex*16 d2dphpa(kmymaxp,kmymaxp,mintfc)
c      complex*16 d2dphb(kmymaxp,kmymaxp,mintfc)
c      complex*16 d2dphpb(kmymaxp,kmymaxp,mintfc)
c      complex*16 d2dthpha(kmymaxp,kmymaxp,mintfc)
c      complex*16 d2dthphpa(kmymaxp,kmymaxp,mintfc)
c      complex*16 d2dthphb(kmymaxp,kmymaxp,mintfc)
c      complex*16 d2dthphpb(kmymaxp,kmymaxp,mintfc)
c
c      complex*16 rmata(lmsup,lmsup,mintfc),rmatpa(lmsup,lmsup,mintfc)
c      complex*16 rmatb(lmsup,lmsup,mintfc),rmatpb(lmsup,lmsup,mintfc)
c
      complex*16 taua(dbogomaxp,dbogomaxp,mintfc)
      complex*16 taub(dbogomaxp,dbogomaxp,mintfc)
      complex*16 gtaua(dbogomaxp,kmymaxp,mintfc)
      complex*16 gtaub(dbogomaxp,dbogomaxp,mintfc)
c      complex*16 tminv(dbogomaxp,dbogomaxp,mintfc)
c      complex*16 ptminva(dbogomaxp,dbogomaxp,mintfc)
c      complex*16 ptminvb(dbogomaxp,dbogomaxp,mintfc)
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
      complex*16 qmoma(lmsup,mintfc),qmomb(lmsup,mintfc)
      complex*16 zdos(dbogomaxp),zrho(nrad),zmagrho(nrad)
      complex*16 zrhosp(nrad,2),zrhodsp(nrad,2)
      complex*16 zrhoeh(nrad),zrhohe(nrad)
      complex*16 zdmag(dbogomaxp,3),zdmor(dbogomaxp,3)
      complex*16 zqmom(dbogomaxp,lmsup)
      complex*16 zscmom(dbogomaxp,lmsup),zsctmom(dbogomaxp,lmsup)
      complex*16 zscdos(dbogomaxp),zsctdos(dbogomaxp)
      complex*16 qmom(lmsup),qmomh(lmsup)
      complex*16 qmomeh,qmomhe,qmomthe,qmomteh
c
      complex*16 qmomha(lmsup,mintfc)
      complex*16 qmomhb(lmsup,mintfc)
      complex*16 qmomeha(mintfc),qmomhea(mintfc)
      complex*16 qmomehb(mintfc),qmomheb(mintfc)
      complex*16 qmomteha(mintfc),qmomthea(mintfc)
      complex*16 qmomtehb(mintfc),qmomtheb(mintfc)
c
      real*8     rz(3),rtmp(3),rtmph(3)
c
      complex*16 rgacoeff(kmymaxp,kmymaxp,lmsup)
      common/rggaunt/rgacoeff
c
      complex*16 sxcoeff(kmymaxp,kmymaxp),sxbcoeff(kmymaxp,kmymaxp)
      complex*16 sycoeff(kmymaxp,kmymaxp),sybcoeff(kmymaxp,kmymaxp)
      complex*16 szcoeff(kmymaxp,kmymaxp),szbcoeff(kmymaxp,kmymaxp)
      common/sigmat/sxcoeff,sxbcoeff,sycoeff,sybcoeff,szcoeff,szbcoeff
c
      complex*16 scoeff(kmymaxp,kmymaxp,3)
      complex*16 sbcoeff(kmymaxp,kmymaxp,3)
      complex*16 lxcoeff(kmymaxp,kmymaxp),lxbcoeff(kmymaxp,kmymaxp)
      complex*16 lycoeff(kmymaxp,kmymaxp),lybcoeff(kmymaxp,kmymaxp)
      complex*16 lzcoeff(kmymaxp,kmymaxp),lzbcoeff(kmymaxp,kmymaxp)
      complex*16 lcoeff(kmymaxp,kmymaxp,3)
      complex*16 lbcoeff(kmymaxp,kmymaxp,3)
      complex*16 deltakmy(kmymaxp,kmymaxp)
      common/lmat/lxcoeff,lxbcoeff,lycoeff,lybcoeff,lzcoeff,lzbcoeff
c
c      complex*16 alphalkkr(0:lmaxp,minprc)
c      complex*16 alpharkkr(0:lmaxp,minprc)
c      complex*16 alphaintkkr(0:lmaxp,mintfc)
c      common/scrpar/alphalkkr,alpharkkr,alphaintkkr
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
c      c=274.072d0 -- from input
      c=c_light
      psq=ce+ce*ce/(c*c)
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
      rz = 0.d0
      rz(3) = 1.d0
c
c *******************************************************************
c * loop over layers to compute physical quantities in global frame *
c *******************************************************************
c 
      qvifc=0.d0
      enbifc=0.d0
      do li=1,nintfc
c       write(6,'(a)') idpota(li)
        cpalay=(1.d0-conc(li)).gt.tiny
c
c       call cpu_time(tstart)
c Compute scattering solutions, call wafu for local frame (pass rb=(0,0,1)!) if necessary
c       --------------------------------------------------------
        if (localmode) then
           rtmp = rz
        else
           rtmp = rba(:,li)
           rtmph(1)=rtmp(1)
           rtmph(2)=-rtmp(2) !TODO: TEST !!! 
           rtmph(3)=rtmp(3)
        end if
        call wafu(ce,c,sxa(li),lmax,idpota(li),v0,E_Fermi,
     >            deltaa(1,li),vra(1,li),bra(1,li),
     >            rtmp,dx(li),rs(li),ns(li),tm,
     >            gz,fz,gj,fj,glz,flz,glj,flj,1)
c       --------------------------------------------------------
c       call cpu_time(tfinish)
c       write(6,'('' Time ellapsed in wafu:'',t35,f8.1,'' s'')')
c    >  tfinish-tstart
c       write(6,*) ' tau from locquant '
c       call outmat1(taua(:,:,li),2*kmymax,2*kmymax,dbogomaxp,tol,6)
c       write(6,*) ' gz(S)'
c       call outmat1(gz(:,:,ns),2*kmymax,2*kmymax,dbogomaxp,tol,6)
c       write(6,*) ' fz(S)'
c       call outmat1(fz(:,:,ns),2*kmymax,2*kmymax,dbogomaxp,tol,6)
c       write(6,*) ' gj(1)'
c       call outmat1(gj(:,:,1),2*kmymax,2*kmymax,dbogomaxp,tol,6)
c       write(6,*) ' fj(1)'
c       call outmat1(fj(:,:,1),2*kmymax,2*kmymax,dbogomaxp,tol,6)
c
c       call cpu_time(tstart)
c Green function
        call gf(lmax,reg,dx(li),ns(li),rs(li),
     >          gz,fz,gj,fj,glz,flz,glj,flj,taua(1,1,li),
     >          ggrmat,fgrmat)
c       call cpu_time(tfinish)
c       write(6,'('' Time ellapsed in gf:'',t35,f8.1,'' s'')')
c     >  tfinish-tstart
c
        if(itest.gt.2) then
          write(6,*) ' ggrmat(0)'
          call outmat1(ggrmat(:,:,0),2*kmymax,2*kmymax,dbogomaxp,tol,6)
          write(6,*) ' fgrmat(0)'
          call outmat1(fgrmat(:,:,0),2*kmymax,2*kmymax,dbogomaxp,tol,6)
        end if
c
c       call cpu_time(tstart)
c Multipole moments
        sig=1.d0
        i=0
        do lam=0,lmaxs
        do mu=-lam,lam
          i=i+1
c         --------------------------------------------------------- e-e
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
c---> Transformation applied -> in gf <- for the offdiag blocks to obtain the anomalous DOS 
c
c         --------------------------------------------------------- ehS           
          call pairing(deltakmy)
          call trafogeh(ggrmat,fgrmat,deltakmy,lmax,tggrmat,tfgrmat) 
          call moment(lmax,lms,rgacoeff(:,:,i),rgacoeff(:,:,i),sig,
     >                tggrmat(1:kmymax,1:kmymax,lam),
     >                tfgrmat(1:kmymax,1:kmymax,lam),
     >                zscmom(1:kmymax,i),rtmp)
c         --------------------------------------------------------- ehT0
          call t0pairing(deltakmy)
          call trafogeh(ggrmat,fgrmat,deltakmy,lmax,tggrmat,tfgrmat) 
          call moment(lmax,lms,rgacoeff(:,:,i),rgacoeff(:,:,i),sig,
     >                tggrmat(1:kmymax,1:kmymax,lam),
     >                tfgrmat(1:kmymax,1:kmymax,lam),
     >                zscmom(kmymax+1:2*kmymax,i),rtmp)
c         ---------------------------------------------------------
c --> Non-Unitary triplet components in the anomalous part
c         --------------------------------------------------------- ehTU            
          call tuppairing(deltakmy)
          call trafogeh(ggrmat,fgrmat,deltakmy,lmax,tggrmat,tfgrmat) 
          call moment(lmax,lms,rgacoeff(:,:,i),rgacoeff(:,:,i),sig,
     >                tggrmat(1:kmymax,1:kmymax,lam),
     >                tfgrmat(1:kmymax,1:kmymax,lam),
     >                zsctmom(1:kmymax,i),rtmp)
c         --------------------------------------------------------- ehTD
          call tdownpairing(deltakmy)
          call trafogeh(ggrmat,fgrmat,deltakmy,lmax,tggrmat,tfgrmat) 
          call moment(lmax,lms,rgacoeff(:,:,i),rgacoeff(:,:,i),sig,
     >                tggrmat(1:kmymax,1:kmymax,lam),
     >                tfgrmat(1:kmymax,1:kmymax,lam),
     >                zsctmom(kmymax+1:2*kmymax,i),rtmp)
c         --------------------------------------------------------- 
        end do
        end do
        do kmy=1,2*kmymax
          zdos(kmy)=zqmom(kmy,1)
          zscdos(kmy)=zscmom(kmy,1)
          zsctdos(kmy)=zsctmom(kmy,1)  
        enddo
c
c Magnetic moments TODO: do it in a for loop!
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
c        call moment(lmax,lms,sycoeff,sybcoeff,sig,
c     >              ggrmat(1:kmymax,1:kmymax,0),
c     >              fgrmat(1:kmymax,1:kmymax,0),zdmag(1:kmymax,2),rtmp)
c        call moment(lmax,lms,szcoeff,szbcoeff,sig,
c     >              ggrmat(1:kmymax,1:kmymax,0),
c     >              fgrmat(1:kmymax,1:kmymax,0),zdmag(1:kmymax,3),rtmp)
        call moment(lmax,.true.,lcoeff(1:kmymaxp,1:kmymaxp,i),
     >              lbcoeff(1:kmymaxp,1:kmymaxp,i),sig,
     >              ggrmat(1:kmymax,1:kmymax,0),
     >              fgrmat(1:kmymax,1:kmymax,0),zdmor(1:kmymax,i),rtmp)
c        call moment(lmax,.true.,lycoeff,lybcoeff,sig,
c     >              ggrmat(1:kmymax,1:kmymax,0),
c     >              fgrmat(1:kmymax,1:kmymax,0),zdmor(1:kmymax,2),rtmp)
c        call moment(lmax,.true.,lzcoeff,lzbcoeff,sig,
c     >              ggrmat(1:kmymax,1:kmymax,0),
c     >              fgrmat(1:kmymax,1:kmymax,0),zdmor(1:kmymax,3),rtmp)
c        ------------------------------------------------- h-h
        call moment(lmax,lms,scoeff(1:kmymaxp,1:kmymaxp,i),
     >              sbcoeff(1:kmymaxp,1:kmymaxp,i),sig,
     >              ggrmat(kmymax+1:2*kmymax,kmymax+1:2*kmymax,0),
     >              fgrmat(kmymax+1:2*kmymax,kmymax+1:2*kmymax,0),
     >              zdmag(kmymax+1:2*kmymax,i),rtmp)
c        call moment(lmax,lms,sycoeff,sybcoeff,sig,
c     >              ggrmat(kmymax+1:2*kmymax,kmymax+1:2*kmymax,0),
c     >              fgrmat(kmymax+1:2*kmymax,kmymax+1:2*kmymax,0),
c     >              zdmag(kmymax+1:2*kmymax,2),rtmp)
c        call moment(lmax,lms,szcoeff,szbcoeff,sig,
c     >              ggrmat(kmymax+1:2*kmymax,kmymax+1:2*kmymax,0),
c     >              fgrmat(kmymax+1:2*kmymax,kmymax+1:2*kmymax,0),
c     >              zdmag(kmymax+1:2*kmymax,3),rtmp)
        call moment(lmax,.true.,lcoeff(1:kmymaxp,1:kmymaxp,i),
     >              lbcoeff(1:kmymaxp,1:kmymaxp,i),sig,
     >              ggrmat(kmymax+1:2*kmymax,kmymax+1:2*kmymax,0),
     >              fgrmat(kmymax+1:2*kmymax,kmymax+1:2*kmymax,0),
     >              zdmor(kmymax+1:2*kmymax,i),rtmp)
c        call moment(lmax,.true.,lycoeff,lybcoeff,sig,
c     >              ggrmat(kmymax+1:2*kmymax,kmymax+1:2*kmymax,0),
c     >              fgrmat(kmymax+1:2*kmymax,kmymax+1:2*kmymax,0),
c     >              zdmor(kmymax+1:2*kmymax,2),rtmp)
c        call moment(lmax,.true.,lzcoeff,lzbcoeff,sig,
c     >              ggrmat(kmymax+1:2*kmymax,kmymax+1:2*kmymax,0),
c     >              fgrmat(kmymax+1:2*kmymax,kmymax+1:2*kmymax,0),
c     >              zdmor(kmymax+1:2*kmymax,3),rtmp)
c        -------------------------------------------------
        end do
c       call cpu_time(tfinish)
c       write(6,'('' Time ellapsed in moment:'',t35,f8.1,'' s'')')
c    >  tfinish-tstart
c
        do k=1,kmymax
          dosa(k,li)=dimag(zdos(k))                  ! electron part
          dosha(k,li)=dimag(zdos(k+kmymax))          ! hole part
          doseha(k,li)=dimag(zscdos(k))               ! electron-hole part 
          doshea(k,li)=dimag(zscdos(k+kmymax))        ! triplet part 
          dosteha(k,li)=dimag(zsctdos(k))               ! triplet up part 
          dosthea(k,li)=dimag(zsctdos(k+kmymax))        ! triplet down part 
          dosmaga(k,li)=dimag(rtmp(1)*zdmag(k,1)     ! electron part of magnetic dos
     >                       +rtmp(2)*zdmag(k,2)
     >                       +rtmp(3)*zdmag(k,3))  ! rtmp either rz or rba
          dosmagha(k,li)=dimag(rtmph(1)*zdmag(k+kmymax,1)    ! hole part of magnetic dos
     >                       +rtmph(2)*zdmag(k+kmymax,2)
     >                       +rtmph(3)*zdmag(k+kmymax,3)) ! rtmph for holes
          if(itest.ge.2) then
            write(6,*) 'DOS electron kmy layeri'
            write(6,*) dosa(k,li),k,li
            write(6,*) 'DOS hole kmy layeri'
            write(6,*) dosha(k,li),k,li
c     >    write(6,'(2d13.5,5x,2d13.5)') dosa(k,li),dosmaga(k,li)
          end if
        end do
c
c     call cpu_time(tstart)
c Radial distribution of charge and magnetic density;
c rtmp defined before wafu based on localmode switch
c        ------------------------------------------------ anomalous density for the e-h part
         call dens(lmax,rs(li),dx(li),ns(li),rtmp,
     >             gz,fz,gj,fj,glz,flz,glj,flj,
     >             taua(1,1,li),zrhoeh,zrho)
c        ------------------------------------------------
c     call cpu_time(tfinish)
c     write(6,'('' Time ellapsed in dens:'',t35,f8.1,'' s'')')
c    > tfinish-tstart
c
c add contribution to contour integral
c
         qvdiffa(li)=0.d0
         qvhdiffa(li)=0.d0
         qvehdiffa(li)=0.d0
         qvhediffa(li)=0.d0
         qvtehdiffa(li)=0.d0
         qvthediffa(li)=0.d0
c
         enbdiffa(li)=0.d0
         enbhdiffa(li)=0.d0
c            
         qv=0.0d0
         qvh=0.0d0
         qveh=(0.0d0,0.0d0)
         qvhe=(0.0d0,0.0d0)
         qvteh=(0.0d0,0.0d0)
         qvthe=(0.0d0,0.0d0)
c
         do k=1,kmymax
c                                                                           TODO:  INITZERO in the main routine (skkr.f)
            qv=qv+dimag(we*zdos(k))               ! electron charge
            qvh=qvh+dimag(we*zdos(k+kmymax))      ! hole charge
            qveh=qveh+we*zscdos(k)                ! e-h singlet charge
            qvhe=qvhe+we*zscdos(k+kmymax)  ! t0 charge
            qvteh=qvteh+we*zsctdos(k)         ! t up charge
            qvthe=qvthe+we*zsctdos(k+kmymax)  ! t down charge
c
            qvdiffa(li)=qvdiffa(li)+dimag(we*zdos(k))                          
            qvhdiffa(li)=qvhdiffa(li)+dimag(we*zdos(k+kmymax))
            qvehdiffa(li)=qvehdiffa(li)+dimag(we*zscdos(k))          
            qvhediffa(li)=qvhediffa(li)+dimag(we*zscdos(k+kmymax))
            qvtehdiffa(li)=qvtehdiffa(li)+dimag(we*zsctdos(k))
            qvthediffa(li)=qvthediffa(li)+dimag(we*zsctdos(k+kmymax))
c
            enbdiffa(li)=enbdiffa(li)+dimag(we*ce*zdos(k))
            enbhdiffa(li)=enbhdiffa(li)+dimag(we*ce*zdos(k+kmymax))
c            enbehdiffa(li)=enbehdiffa(li)+dimag(we*ce*zscdos(k))
c            enbhediffa(li)=enbhediffa(li)+dimag(we*ce*zscdos(k+kmymax))
c            enbtehdiffa(li)=enbtehdiffa(li)+dimag(we*ce*zsctdos(k))
c            enbthediffa(li)=enbthediffa(li)+dimag(we*ce*zsctdos(k+kmymax))
c
            qvpa(k,li)=qvpa(k,li)+dimag(we*zdos(k))
            qvhpa(k,li)=qvhpa(k,li)+dimag(we*zdos(k+kmymax))
            qvehpa(k,li)=qvehpa(k,li)+dimag(we*zscdos(k))
            qvhepa(k,li)=qvhepa(k,li)+dimag(we*zscdos(k+kmymax))
            qvtehpa(k,li)=qvtehpa(k,li)+dimag(we*zsctdos(k))
            qvthepa(k,li)=qvthepa(k,li)+dimag(we*zsctdos(k+kmymax))
c
            qva(li)=qva(li)+dimag(we*zdos(k))
            qvha(li)=qvha(li)+dimag(we*zdos(k+kmymax))
            qveha(li)=qveha(li)+dimag(we*zscdos(k))
            qvhea(li)=qvhea(li)+dimag(we*zscdos(k+kmymax))
            qvteha(li)=qvteha(li)+dimag(we*zsctdos(k))
            qvthea(li)=qvthea(li)+dimag(we*zsctdos(k+kmymax))
c
c            write(6,*) li,'qveha',qveha(li) 
            do i=1,3
c
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
c
            end do
c
            enba(li)=enba(li)+dimag(we*ce*zdos(k))                 !!! TODO: IMPLEMENT THE LLOYD FORMULA (+enbll)
            enbha(li)=enbha(li)+dimag(we*ce*zdos(k+kmymax))
c            enbeha(li)=enbeha(li)+dimag(we*ce*zscdos(k))
c            enbhea(li)=enbhea(li)+dimag(we*ce*zscdos(k+kmymax))
c            enbteha(li)=enbteha(li)+dimag(we*ce*zsctdos(k))
c            enbthea(li)=enbthea(li)+dimag(we*ce*zsctdos(k+kmymax))
c
         end do
c                                                                            
         do irad=1,ns(li)                                        
c            rhov(irad)=dimag(we*zrho(irad))
            rhova(irad,li)=rhova(irad,li)+dimag(we*zrho(irad))
c            rhoveh(irad)=we*zrhoeh(irad)
            rhoveha(irad,li)=rhoveha(irad,li)+we*zrhoeh(irad)
         end do
c
         do lam=0,lmaxs
         do mu=-lam,lam
            i=lam*(lam+1)+mu+1
            ii=lam*(lam+1)-mu+1
            qmom(i)=(0.0d0,0.0d0)
            qmomh(i)=(0.0d0,0.0d0)
            do k=1,kmymax
              qmom(i)=qmom(i)+cih*(we*zqmom(k,i)-
     >                             mpar(mu)*dconjg(we*zqmom(k,ii)))
              qmomh(i)=qmomh(i)+cih*(we*zqmom(k+kmymax,i)-
     >                         mpar(mu)*dconjg(we*zqmom(k+kmymax,ii)))
            end do
          end do
          end do
c
          qmomeh=(0.0d0,0.0d0)
          qmomhe=(0.0d0,0.0d0)
          qmomteh=(0.0d0,0.0d0)
          qmomthe=(0.0d0,0.0d0)
c
          do k=1,kmymax
              qmomeh=qmomeh+dimag(we*zscdos(k))
              qmomhe=qmomhe+dimag(we*zscdos(k+kmymax))
              qmomteh=qmomteh+dimag(we*zsctdos(k))        
              qmomthe=qmomthe+dimag(we*zsctdos(k+kmymax))
          end do
c
c Add moments of charge density to the integral
c
          qmoma(1:lmmaxs,li) = qmoma(1:lmmaxs,li) + qmom(1:lmmaxs)
          qmomha(1:lmmaxs,li) = qmomha(1:lmmaxs,li) + qmomh(1:lmmaxs)
          qmomeha(li) = qmomeha(li) + qmomeh
          qmomhea(li) = qmomhea(li) + qmomhe
          qmomteha(li) = qmomteha(li) + qmomteh
          qmomthea(li) = qmomthea(li) + qmomthe          
c
c Check charges
c
c         x=dlog(rs(li))-(ns(li)-1)*dx(li)
c         do irad=1,ns(li)
c           rad(irad)=dexp(x)
c           x=x+dx(li)
c         end do
c         qvp=rsimp(dimag(rhovhe),rad,ns(li),dx(li))
c         write(6,'(''Layer'',i5,''  qvhe,qvphe:'',3d20.10)') li,qvhe,qvp
c
         if(itest.ge.2) then
           write(6,'('' Contribution to contour integral so far'')')
c           write(6,'(''Q    '',d13.5)') qva(li)
           do i=1,lmmaxs
             write(6,'(''Qmom - electron '',i2,2x,2d13.5)') i,qmoma(i,li)
             write(6,'(''Qmom - hole     '',i2,2x,2d13.5)') i,qmomha(i,li)
           end do
           write(6,'(''Q eh spin offdiag '',2d13.5)') qmomeha(li)
           write(6,'(''Q T0 spin offdiag '',2d13.5)') qmomhea(li)
           write(6,'(''Q T up      '',2d13.5)') qmomteha(li)
           write(6,'(''Q T down     '',2d13.5)') qmomthea(li)
           write(6,'(''Eband -> electron '',d13.5)') enba(li)
           write(6,'(''Eband -> hole     '',d13.5)') enbha(li)
           write(6,'(''Sx,Lx -> electron '',d13.5,5x,d13.5)') 
     >               spin_magva(li,1),orb_magva(li,1)
           write(6,'(''Sx,Lx -> hole     '',d13.5,5x,d13.5)') 
     >               spin_magvha(li,1),orb_magvha(li,1)
           write(6,'(''Sy,Ly -> electron '',d13.5,5x,d13.5)') 
     >               spin_magva(li,2),orb_magva(li,2)
           write(6,'(''Sy,Ly -> hole     '',d13.5,5x,d13.5)') 
     >               spin_magvha(li,2),orb_magvha(li,2)
           write(6,'(''Sz,Lz -> electron '',d13.5,5x,d13.5)') 
     >               spin_magva(li,3),orb_magva(li,3)
           write(6,'(''Sz,Lz -> hole     '',d13.5,5x,d13.5)') 
     >               spin_magvha(li,3),orb_magvha(li,3)
           if(itest.gt.2) then
            write(6,'(''Loop: rad - rhoeh (re,im) and rho '')')
            do irad=1,ns(li)
              write(6,'(d13.5,5x,d13.5,5x,d13.5,5x,d13.5,5x,d13.5)') 
     >        rad(irad),dreal(rhoveha(irad,li)),dimag(rhoveha(irad,li)),
     >                  rhova(irad,li)
            enddo
           endif
         endif
c
         if(cpalay) then
c+------------+
c+ BIG CPA IF +
c+------------+
c       call cpu_time(tstart)
c Compute scattering solutions, call wafu for local frame (pass rb=(0,0,1)!) if necessary
c       --------------------------------------------------------
        if (localmode) then
           rtmp = rz
        else
           rtmp = rbb(:,li)
           rtmph(1)=rtmp(1)
           rtmph(2)=-rtmp(2) !TODO: TEST !!! 
           rtmph(3)=rtmp(3)
        end if
        call wafu(ce,c,sxb(li),lmax,idpotb(li),v0,E_Fermi,
     >            deltab(1,li),vrb(1,li),brb(1,li),
     >            rtmp,dx(li),rs(li),ns(li),tm,
     >            gz,fz,gj,fj,glz,flz,glj,flj,1)
c        write(6,*) 'delta B for li=',li
c        write(6,*) deltab(:,li)
c       --------------------------------------------------------
c       call cpu_time(tfinish)
c       write(6,'('' Time ellapsed in wafu:'',t35,f8.1,'' s'')')
c    >  tfinish-tstart
c       write(6,*) ' tau from locquant '
c       call outmat1(taua(:,:,li),2*kmymax,2*kmymax,dbogomaxp,tol,6)
c       write(6,*) ' gz(S)'
c       call outmat1(gz(:,:,ns),2*kmymax,2*kmymax,dbogomaxp,tol,6)
c       write(6,*) ' fz(S)'
c       call outmat1(fz(:,:,ns),2*kmymax,2*kmymax,dbogomaxp,tol,6)
c       write(6,*) ' gj(1)'
c       call outmat1(gj(:,:,1),2*kmymax,2*kmymax,dbogomaxp,tol,6)
c       write(6,*) ' fj(1)'
c       call outmat1(fj(:,:,1),2*kmymax,2*kmymax,dbogomaxp,tol,6)
c
c       call cpu_time(tstart)
c Green function
        call gf(lmax,reg,dx(li),ns(li),rs(li),
     >          gz,fz,gj,fj,glz,flz,glj,flj,taub(1,1,li),
     >          ggrmat,fgrmat)
c       call cpu_time(tfinish)
c       write(6,'('' Time ellapsed in gf:'',t35,f8.1,'' s'')')
c     >  tfinish-tstart
c
c We calculated the scattering solution and tha GF for the B type, GF and wavefunctions stored in the same variable
        if(itest.gt.2) then
          write(6,*) ' ggrmat(0)'
          call outmat1(ggrmat(:,:,0),2*kmymax,2*kmymax,dbogomaxp,tol,6)
          write(6,*) ' fgrmat(0)'
          call outmat1(fgrmat(:,:,0),2*kmymax,2*kmymax,dbogomaxp,tol,6)
        end if
c
c
c       call cpu_time(tstart)
c Multipole moments
        sig=1.d0
        i=0
        do lam=0,lmaxs
        do mu=-lam,lam
          i=i+1
c         --------------------------------------------------------- e-e
          call moment(lmax,lms,rgacoeff(:,:,i),rgacoeff(:,:,i),sig,
     >                ggrmat(1:kmymax,1:kmymax,lam),
     >            fgrmat(1:kmymax,1:kmymax,lam),zqmom(1:kmymax,i),rtmp)
c         --------------------------------------------------------- h-h
          call moment(lmax,lms,rgacoeff(:,:,i),rgacoeff(:,:,i),sig,
     >                ggrmat(kmymax+1:2*kmymax,kmymax+1:2*kmymax,lam),
     >                fgrmat(kmymax+1:2*kmymax,kmymax+1:2*kmymax,lam),
     >                zqmom(kmymax+1:2*kmymax,i),rtmph)
c
c---> Transformation applied -> in gf <- for the offdiag blocks to obtain the anomalous DOS 
c
c         --------------------------------------------------------- ehS           
          call pairing(deltakmy)
          call trafogeh(ggrmat,fgrmat,deltakmy,lmax,tggrmat,tfgrmat) 
          call moment(lmax,lms,rgacoeff(:,:,i),rgacoeff(:,:,i),sig,
     >                tggrmat(1:kmymax,1:kmymax,lam),
     >                tfgrmat(1:kmymax,1:kmymax,lam),
     >                zscmom(1:kmymax,i),rtmp)
c         --------------------------------------------------------- ehT0
          call t0pairing(deltakmy)
          call trafogeh(ggrmat,fgrmat,deltakmy,lmax,tggrmat,tfgrmat) 
          call moment(lmax,lms,rgacoeff(:,:,i),rgacoeff(:,:,i),sig,
     >                tggrmat(1:kmymax,1:kmymax,lam),
     >                tfgrmat(1:kmymax,1:kmymax,lam),
     >                zscmom(kmymax+1:2*kmymax,i),rtmp)
c         ---------------------------------------------------------
c --> Non-Unitary triplet components in the anomalous part
c         --------------------------------------------------------- ehTU            
          call tuppairing(deltakmy)
          call trafogeh(ggrmat,fgrmat,deltakmy,lmax,tggrmat,tfgrmat) 
          call moment(lmax,lms,rgacoeff(:,:,i),rgacoeff(:,:,i),sig,
     >                tggrmat(1:kmymax,1:kmymax,lam),
     >                tfgrmat(1:kmymax,1:kmymax,lam),
     >                zsctmom(1:kmymax,i),rtmp)
c         --------------------------------------------------------- ehTD
          call tdownpairing(deltakmy)
          call trafogeh(ggrmat,fgrmat,deltakmy,lmax,tggrmat,tfgrmat) 
          call moment(lmax,lms,rgacoeff(:,:,i),rgacoeff(:,:,i),sig,
     >                tggrmat(1:kmymax,1:kmymax,lam),
     >                tfgrmat(1:kmymax,1:kmymax,lam),
     >                zsctmom(kmymax+1:2*kmymax,i),rtmp)
c         --------------------------------------------------------- 
        end do
        end do
        do kmy=1,2*kmymax
          zdos(kmy)=zqmom(kmy,1)
          zscdos(kmy)=zscmom(kmy,1)
          zsctdos(kmy)=zsctmom(kmy,1)  
        end do
c
c Magnetic moments TODO: do it in a for loop 
        sig=-1.d0
        do i=1,3
c       ------------------------------------------------- e-e 
        call moment(lmax,lms,scoeff(1:kmymaxp,1:kmymaxp,i),
     >              sbcoeff(1:kmymaxp,1:kmymaxp,i),sig,
     >              ggrmat(1:kmymax,1:kmymax,0),
     >              fgrmat(1:kmymax,1:kmymax,0),zdmag(1:kmymax,i),rtmp)
c        call moment(lmax,lms,sycoeff,sybcoeff,sig,
c     >              ggrmat(1:kmymax,1:kmymax,0),
c     >              fgrmat(1:kmymax,1:kmymax,0),zdmag(1:kmymax,2),rtmp)
c        call moment(lmax,lms,szcoeff,szbcoeff,sig,
c     >              ggrmat(1:kmymax,1:kmymax,0),
c     >              fgrmat(1:kmymax,1:kmymax,0),zdmag(1:kmymax,3),rtmp)
        call moment(lmax,.true.,lcoeff(1:kmymaxp,1:kmymaxp,i),
     >              lbcoeff(1:kmymaxp,1:kmymaxp,i),sig,
     >              ggrmat(1:kmymax,1:kmymax,0),
     >              fgrmat(1:kmymax,1:kmymax,0),zdmor(1:kmymax,i),rtmp)
c        call moment(lmax,.true.,lycoeff,lybcoeff,sig,
c     >              ggrmat(1:kmymax,1:kmymax,0),
c     >              fgrmat(1:kmymax,1:kmymax,0),zdmor(1:kmymax,2),rtmp)
c        call moment(lmax,.true.,lzcoeff,lzbcoeff,sig,
c     >              ggrmat(1:kmymax,1:kmymax,0),
c     >              fgrmat(1:kmymax,1:kmymax,0),zdmor(1:kmymax,3),rtmp)
c        ------------------------------------------------- h-h
        call moment(lmax,lms,scoeff(1:kmymaxp,1:kmymaxp,i),
     >              sbcoeff(1:kmymaxp,1:kmymaxp,i),sig,
     >              ggrmat(kmymax+1:2*kmymax,kmymax+1:2*kmymax,0),
     >              fgrmat(kmymax+1:2*kmymax,kmymax+1:2*kmymax,0),
     >              zdmag(kmymax+1:2*kmymax,i),rtmp)
c        call moment(lmax,lms,sycoeff,sybcoeff,sig,
c     >              ggrmat(kmymax+1:2*kmymax,kmymax+1:2*kmymax,0),
c     >              fgrmat(kmymax+1:2*kmymax,kmymax+1:2*kmymax,0),
c     >              zdmag(kmymax+1:2*kmymax,2),rtmp)
c        call moment(lmax,lms,szcoeff,szbcoeff,sig,
c     >              ggrmat(kmymax+1:2*kmymax,kmymax+1:2*kmymax,0),
c     >              fgrmat(kmymax+1:2*kmymax,kmymax+1:2*kmymax,0),
c     >              zdmag(kmymax+1:2*kmymax,3),rtmp)
        call moment(lmax,.true.,lcoeff(1:kmymaxp,1:kmymaxp,i),
     >              lbcoeff(1:kmymaxp,1:kmymaxp,i),sig,
     >              ggrmat(kmymax+1:2*kmymax,kmymax+1:2*kmymax,0),
     >              fgrmat(kmymax+1:2*kmymax,kmymax+1:2*kmymax,0),
     >              zdmor(kmymax+1:2*kmymax,i),rtmp)
c        call moment(lmax,.true.,lycoeff,lybcoeff,sig,
c     >              ggrmat(kmymax+1:2*kmymax,kmymax+1:2*kmymax,0),
c     >              fgrmat(kmymax+1:2*kmymax,kmymax+1:2*kmymax,0),
c     >              zdmor(kmymax+1:2*kmymax,2),rtmp)
c        call moment(lmax,.true.,lzcoeff,lzbcoeff,sig,
c     >              ggrmat(kmymax+1:2*kmymax,kmymax+1:2*kmymax,0),
c     >              fgrmat(kmymax+1:2*kmymax,kmymax+1:2*kmymax,0),
c     >              zdmor(kmymax+1:2*kmymax,3),rtmp)
c        -------------------------------------------------
        end do
c        -------------------------------------------------
c       call cpu_time(tfinish)
c       write(6,'('' Time ellapsed in moment:'',t35,f8.1,'' s'')')
c    >  tfinish-tstart
c
        do k=1,kmymax
          dosb(k,li)=dimag(zdos(k))                  ! electron part
          doshb(k,li)=dimag(zdos(k+kmymax))          ! hole part
          dosehb(k,li)=dimag(zscdos(k))               ! electron-hole part 
          dosheb(k,li)=dimag(zscdos(k+kmymax))        ! triplet part 
          dostehb(k,li)=dimag(zsctdos(k))               ! triplet up part 
          dostheb(k,li)=dimag(zsctdos(k+kmymax))        ! triplet down part 
          dosmagb(k,li)=dimag(rtmp(1)*zdmag(k,1)     ! electron part of magnetic dos
     >                       +rtmp(2)*zdmag(k,2)
     >                       +rtmp(3)*zdmag(k,3))  ! rtmp either rz or rba
          dosmaghb(k,li)=dimag(rtmph(1)*zdmag(k+kmymax,1)    ! hole part of magnetic dos
     >                       +rtmph(2)*zdmag(k+kmymax,2)
     >                       +rtmph(3)*zdmag(k+kmymax,3)) ! rtmph for holes
          if(itest.ge.2) then
            write(6,*) 'DOS electron kmy layeri'
            write(6,*) dosb(k,li),k,li
            write(6,*) 'DOS hole kmy layeri'
            write(6,*) doshb(k,li),k,li
c     >    write(6,'(2d13.5,5x,2d13.5)') dosa(k,li),dosmaga(k,li)
          end if
        end do
c
c     call cpu_time(tstart)
c Radial distribution of charge and magnetic density;
c rtmp defined before wafu based on localmode switch
c        ------------------------------------------------ anomalous density for the e-h part
         call dens(lmax,rs(li),dx(li),ns(li),rtmp,
     >             gz,fz,gj,fj,glz,flz,glj,flj,
     >             taub(1,1,li),zrhoeh,zrho)
c        ------------------------------------------------
c     call cpu_time(tfinish)
c     write(6,'('' Time ellapsed in dens:'',t35,f8.1,'' s'')')
c    > tfinish-tstart
c
c add contribution to contour integral
c
         qvdiffb(li)=0.d0
         qvhdiffb(li)=0.d0
         qvehdiffb(li)=0.d0
         qvhediffb(li)=0.d0
         qvtehdiffb(li)=0.d0
         qvthediffb(li)=0.d0
c
         enbdiffb(li)=0.d0
         enbhdiffb(li)=0.d0
c            
         do k=1,kmymax
c                                                                           TODO:  INITZERO in the main routine (skkr.f)
            qvdiffb(li)=qvdiffb(li)+dimag(we*zdos(k))                          
            qvhdiffb(li)=qvhdiffb(li)+dimag(we*zdos(k+kmymax))
            qvehdiffb(li)=qvehdiffb(li)+dimag(we*zscdos(k))          
            qvhediffb(li)=qvhediffb(li)+dimag(we*zscdos(k+kmymax))
            qvtehdiffb(li)=qvtehdiffb(li)+dimag(we*zsctdos(k))
            qvthediffb(li)=qvthediffb(li)+dimag(we*zsctdos(k+kmymax))
c
            enbdiffb(li)=enbdiffb(li)+dimag(we*ce*zdos(k))
            enbhdiffb(li)=enbhdiffb(li)+dimag(we*ce*zdos(k+kmymax))
c            enbehdiffa(li)=enbehdiffa(li)+dimag(we*ce*zscdos(k))
c            enbhediffa(li)=enbhediffa(li)+dimag(we*ce*zscdos(k+kmymax))
c            enbtehdiffa(li)=enbtehdiffa(li)+dimag(we*ce*zsctdos(k))
c            enbthediffa(li)=enbthediffa(li)+dimag(we*ce*zsctdos(k+kmymax))
c
            qvpb(k,li)=qvpb(k,li)+dimag(we*zdos(k))
            qvhpb(k,li)=qvhpb(k,li)+dimag(we*zdos(k+kmymax))
            qvehpb(k,li)=qvehpb(k,li)+dimag(we*zscdos(k))
            qvhepb(k,li)=qvhepb(k,li)+dimag(we*zscdos(k+kmymax))
            qvtehpb(k,li)=qvtehpb(k,li)+dimag(we*zsctdos(k))
            qvthepb(k,li)=qvthepb(k,li)+dimag(we*zsctdos(k+kmymax))
c
            qvb(li)=qvb(li)+dimag(we*zdos(k))
            qvhb(li)=qvhb(li)+dimag(we*zdos(k+kmymax))
            qvehb(li)=qvehb(li)+dimag(we*zscdos(k))
            qvheb(li)=qvheb(li)+dimag(we*zscdos(k+kmymax))
            qvtehb(li)=qvtehb(li)+dimag(we*zsctdos(k))
            qvtheb(li)=qvtheb(li)+dimag(we*zsctdos(k+kmymax))
c
c            write(6,*) li,'qvehb',qvehb(li) 
            do i=1,3
c
              spin_magvpb(k,li,i)=spin_magvpb(k,li,i)+
     >                            dimag(we*zdmag(k,i))
              spin_magvhpb(k,li,i)=spin_magvhpb(k,li,i)+
     >                            dimag(we*zdmag(k+kmymax,i))
c
              spin_magvb(li,i)=spin_magvb(li,i)+
     >                         dimag(we*zdmag(k,i))
              spin_magvhb(li,i)=spin_magvhb(li,i)+
     >                         dimag(we*zdmag(k+kmymax,i))
c
              orb_magvpb(k,li,i)=orb_magvpb(k,li,i)+
     >                           dimag(we*zdmor(k,i))
              orb_magvhpb(k,li,i)=orb_magvhpb(k,li,i)+
     >                           dimag(we*zdmor(k+kmymax,i))
c
              orb_magvb(li,i)=orb_magvb(li,i)+
     >                        dimag(we*zdmor(k,i))
              orb_magvhb(li,i)=orb_magvhb(li,i)+
     >                        dimag(we*zdmor(k+kmymax,i))
c
            end do
c
            enbb(li)=enbb(li)+dimag(we*ce*zdos(k))                 !!! TODO: IMPLEMENT THE LLOYD FORMULA (+enbll)
            enbhb(li)=enbhb(li)+dimag(we*ce*zdos(k+kmymax))
c            enbeha(li)=enbeha(li)+dimag(we*ce*zscdos(k))
c            enbhea(li)=enbhea(li)+dimag(we*ce*zscdos(k+kmymax))
c            enbteha(li)=enbteha(li)+dimag(we*ce*zsctdos(k))
c            enbthea(li)=enbthea(li)+dimag(we*ce*zsctdos(k+kmymax))
c
         end do
c                                                                            
         do irad=1,ns(li)                                        
            rhovb(irad,li)=rhovb(irad,li)+dimag(we*zrho(irad))
            rhovehb(irad,li)=rhovehb(irad,li)+we*zrhoeh(irad)
c            rhov(irad)=dimag(we*zrho(irad))
c            rhoveh(irad)=we*zrhoeh(irad)
         end do
c
         do lam=0,lmaxs
         do mu=-lam,lam
            i=lam*(lam+1)+mu+1
            ii=lam*(lam+1)-mu+1
            qmom(i)=(0.0d0,0.0d0)
            qmomh(i)=(0.0d0,0.0d0)
            do k=1,kmymax
              qmom(i)=qmom(i)+cih*(we*zqmom(k,i)-
     >                             mpar(mu)*dconjg(we*zqmom(k,ii)))
              qmomh(i)=qmomh(i)+cih*(we*zqmom(k+kmymax,i)-
     >                         mpar(mu)*dconjg(we*zqmom(k+kmymax,ii)))
            end do
          end do
          end do
c
          qmomeh=(0.0d0,0.0d0)
          qmomhe=(0.0d0,0.0d0)
          qmomteh=(0.0d0,0.0d0)
          qmomthe=(0.0d0,0.0d0)
c
          do k=1,kmymax
              qmomeh=qmomeh+dimag(we*zscdos(k))
              qmomhe=qmomhe+dimag(we*zscdos(k+kmymax))
              qmomteh=qmomteh+dimag(we*zsctdos(k))        
              qmomthe=qmomthe+dimag(we*zsctdos(k+kmymax))
          end do
c
c Add moments of charge density to the integral
c
          qmomb(1:lmmaxs,li) = qmomb(1:lmmaxs,li) + qmom(1:lmmaxs)
          qmomhb(1:lmmaxs,li) = qmomhb(1:lmmaxs,li) + qmomh(1:lmmaxs)
          qmomehb(li) = qmomehb(li) + qmomeh
          qmomheb(li) = qmomheb(li) + qmomhe
          qmomtehb(li) = qmomtehb(li) + qmomteh
          qmomtheb(li) = qmomtheb(li) + qmomthe          
c
c Check charges
c
c         x=dlog(rs(li))-(ns(li)-1)*dx(li)
c         do irad=1,ns(li)
c           rad(irad)=dexp(x)
c           x=x+dx(li)
c         end do
c         qvp=rsimp(dimag(rhovhe),rad,ns(li),dx(li))
c         write(6,'(''Layer'',i5,''  qvhe,qvphe:'',3d20.10)') li,qvhe,qvp
c
         if(itest.ge.2) then
           write(6,'('' Contribution to contour integral so far'')')
c           write(6,'(''Q    '',d13.5)') qvb(li)
           do i=1,lmmaxs
             write(6,'(''Qmom - electron '',i2,2x,2d13.5)') i,qmomb(i,li)
             write(6,'(''Qmom - hole     '',i2,2x,2d13.5)') i,qmomhb(i,li)
           end do
           write(6,'(''Q eh spin offdiag '',2d13.5)') qmomehb(li)
           write(6,'(''Q T0 spin offdiag '',2d13.5)') qmomheb(li)
           write(6,'(''Q T up      '',2d13.5)') qmomtehb(li)
           write(6,'(''Q T down     '',2d13.5)') qmomtheb(li)
           write(6,'(''Eband -> electron '',d13.5)') enbb(li)
           write(6,'(''Eband -> hole     '',d13.5)') enbhb(li)
           write(6,'(''Sx,Lx -> electron '',d13.5,5x,d13.5)') 
     >               spin_magvb(li,1),orb_magvb(li,1)
           write(6,'(''Sx,Lx -> hole     '',d13.5,5x,d13.5)') 
     >               spin_magvhb(li,1),orb_magvhb(li,1)
           write(6,'(''Sy,Ly -> electron '',d13.5,5x,d13.5)') 
     >               spin_magvb(li,2),orb_magvb(li,2)
           write(6,'(''Sy,Ly -> hole     '',d13.5,5x,d13.5)') 
     >               spin_magvhb(li,2),orb_magvhb(li,2)
           write(6,'(''Sz,Lz -> electron '',d13.5,5x,d13.5)') 
     >               spin_magvb(li,3),orb_magvb(li,3)
           write(6,'(''Sz,Lz -> hole     '',d13.5,5x,d13.5)') 
     >               spin_magvhb(li,3),orb_magvhb(li,3)
           if(itest.gt.2) then
            write(6,'(''Loop: rad - rhoeh (re,im) and rho '')')
            do irad=1,ns(li)
              write(6,'(d13.5,5x,d13.5,5x,d13.5,5x,d13.5,5x,d13.5)') 
     >        rad(irad),dreal(rhovehb(irad,li)),dimag(rhovehb(irad,li)),
     >                  rhovb(irad,li)
            enddo
           endif
         endif
c
         else
c If we have no CPA then B component equals to A component  
           do k=1,kmymax
             qvpb(k,li)=qvpa(k,li)
             qvhpb(k,li)=qvhpa(k,li)
             qvehpb(k,li)=qvehpa(k,li)
             qvhepb(k,li)=qvhepa(k,li)
             qvtehpb(k,li)=qvtehpa(k,li)
             qvthepb(k,li)=qvthepa(k,li)
c
             do i=1,3
               spin_magvpb(k,li,i)=spin_magvpa(k,li,i)
               orb_magvpb(k,li,i)=orb_magvpa(k,li,i)
             end do
           end do
           do i=1,lmmaxs
             qmomb(i,li)=qmoma(i,li)
             qmomhb(i,li)=qmomha(i,li)
           end do
           qmomehb(li) = qmomeha(li)
           qmomheb(li) = qmomhea(li)
           qmomtehb(li) = qmomteha(li)
           qmomtheb(li) = qmomthea(li)          
           do i=1,3
             spin_magvpb(k,li,i)=spin_magvpa(k,li,i)
             spin_magvhpb(k,li,i)=spin_magvhpa(k,li,i)
             spin_magvb(li,i)=spin_magva(li,i)
             spin_magvhb(li,i)=spin_magvha(li,i)
             orb_magvpb(k,li,i)=orb_magvpa(k,li,i)
             orb_magvhpb(k,li,i)=orb_magvhpa(k,li,i)
             orb_magvb(li,i)=orb_magva(li,i)
             orb_magvhb(li,i)=orb_magvha(li,i)
c
           end do
           qvdiffb(li)=qvdiffa(li)
           enbdiffb(li)=enbdiffa(li)
           qvb(li)=qva(li)
           qvhb(li)=qvha(li)
           qvehb(li)=qveha(li)
           qvheb(li)=qvhea(li)
           qvtehb(li)=qvteha(li)
           qvtheb(li)=qvthea(li)

           enbb(li)=enba(li)
           enbhb(li)=enbha(li)
c           enorbb(li)=enorba(li)
           do irad=1,ns(li)
             rhovb(irad,li)=rhova(irad,li)
             rhovehb(irad,li)=rhoveha(irad,li)
c             rhospb(irad,1,li)=rhospa(irad,1,li)
c             rhospb(irad,2,li)=rhospa(irad,2,li)
c             rhodspb(irad,1,li)=rhodspa(irad,1,li)
c             rhodspb(irad,2,li)=rhodspa(irad,2,li)
c             rhomagb(irad,li)=rhomaga(irad,li)
           end do
c

c
         end if
c+----------------+
c+ END BIG CPA IF +
c+----------------+
c
c         enbifc=enbifc+conc(li)*enbdiffa(li)
c     >                +(1.0d0-conc(li))*enbdiffb(li)
c         qvifc=qvifc+conc(li)*qvdiffa(li)
c     >              +(1.0d0-conc(li))*qvdiffb(li)
c
      end do
c ************************
c * end loop over layers *
c ************************

c      if(linbw) then
c        omifc=enbifc-efermi*qvifc
c        enbifcint=enbifcint+enbifc
c        qvifcint=qvifcint+qvifc
c        omifcint=omifcint+omifc
c        write(6,'('' Enb'',i3,2f10.5,3d17.8)') 
c     >  ie,ce,enbifc,qvifc,omifc
c      end if
c
      return
      end
