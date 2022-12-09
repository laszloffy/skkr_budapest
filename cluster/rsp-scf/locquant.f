c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine locquant(
c     ===================
     > ie,ce,we,lmax,nintfc,wrel,lms,sxa,sxb,conc,v0,
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
     > dosa,qvpa,qva,qvdiffa,
     > enba,enbdiffa,denba,enorba,qmoma,
     > rhova,rhospa,rhodspa,rhomaga,
     > dosmaga,spin_magvpa,spin_magva,orb_magvpa,orb_magva,
c
     > dosb,qvpb,qvb,qvdiffb,
     > enbb,enbdiffb,denbb,enorbb,qmomb,
     > rhovb,rhospb,rhodspb,rhomagb,
     > dosmagb,spin_magvpb,spin_magvb,orb_magvpb,orb_magvb,
c
     > lliter,linbw,efermi,enbifcint,qvifcint,omifcint)
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      logical wrel,lms,cpalay,linbw
c
      character*10 idpota(mintfc),idpotb(mintfc)
      character*30 for006
c
      dimension nuz(kmymaxp),indz(nuzp,kmymaxp)
      dimension mpar(-lsup:lsup)
c
      dimension vra(nrad,mintfc),bra(nrad,mintfc)
      dimension bopra(nrad,2,mintfc)
      dimension vrb(nrad,mintfc),brb(nrad,mintfc)
      dimension boprb(nrad,2,mintfc)
      dimension rs(mintfc),dx(mintfc),ns(mintfc)
c
      dimension enba(mintfc),enbb(mintfc)
      dimension denba(mintfc),denbb(mintfc)
      dimension enorba(mintfc),enorbb(mintfc)
c
      dimension dosa(kmymaxp,mintfc),dosb(kmymaxp,mintfc)
      dimension dosmaga(kmymaxp,mintfc),dosmagb(kmymaxp,mintfc)
      dimension qva(mintfc),qvb(mintfc)
      dimension qvpa(kmymaxp,mintfc),qvpb(kmymaxp,mintfc)
      dimension spin_magvpa(kmymaxp,mintfc,3)
      dimension spin_magvpb(kmymaxp,mintfc,3)
      dimension spin_magva(mintfc,3),spin_magvb(mintfc,3)
      dimension orb_magvpa(kmymaxp,mintfc,3)
      dimension orb_magvpb(kmymaxp,mintfc,3)
      dimension orb_magva(mintfc,3),orb_magvb(mintfc,3)
c
      dimension enbdiffa(mintfc),enbdiffb(mintfc)
      dimension qvdiffa(mintfc),qvdiffb(mintfc)
c
      dimension rhova(nrad,mintfc),rhovb(nrad,mintfc)
      dimension rhospa(nrad,2,mintfc),rhospb(nrad,2,mintfc)
      dimension rhodspa(nrad,2,mintfc),rhodspb(nrad,2,mintfc)
      dimension rhomaga(nrad,mintfc),rhomagb(nrad,mintfc)
      dimension rhov(nrad),rad(nrad)
c
      dimension conc(mintfc),sxa(mintfc),sxb(mintfc)
c
      complex*16 ce,we,psq,cih
c
      complex*16 dmata(kmymaxp,kmymaxp,mintfc)
      complex*16 dmatpa(kmymaxp,kmymaxp,mintfc)
      complex*16 dmatb(kmymaxp,kmymaxp,mintfc)
      complex*16 dmatpb(kmymaxp,kmymaxp,mintfc)
      complex*16 ddpha(kmymaxp,kmymaxp,mintfc)
      complex*16 ddphb(kmymaxp,kmymaxp,mintfc)
      complex*16 ddphpa(kmymaxp,kmymaxp,mintfc)
      complex*16 ddphpb(kmymaxp,kmymaxp,mintfc)
      complex*16 ddtha(kmymaxp,kmymaxp,mintfc)
      complex*16 ddthb(kmymaxp,kmymaxp,mintfc)
      complex*16 ddthpa(kmymaxp,kmymaxp,mintfc)
      complex*16 ddthpb(kmymaxp,kmymaxp,mintfc)
      complex*16 rmata(lmsup,lmsup,mintfc),rmatpa(lmsup,lmsup,mintfc)
      complex*16 rmatb(lmsup,lmsup,mintfc),rmatpb(lmsup,lmsup,mintfc)
c
      complex*16 taua(kmymaxp,kmymaxp,mintfc)
      complex*16 taub(kmymaxp,kmymaxp,mintfc)
      complex*16 gtaua(kmymaxp,kmymaxp,mintfc)
      complex*16 gtaub(kmymaxp,kmymaxp,mintfc)
      complex*16 tminv(kmymaxp,kmymaxp,mintfc)
      complex*16 ptminva(kmymaxp,kmymaxp,mintfc)
      complex*16 ptminvb(kmymaxp,kmymaxp,mintfc)
      complex*16 tm(kmymaxp,kmymaxp)
c
      complex*16 gz(nrad,nuzp,kmymaxp),fz(nrad,nuzp,kmymaxp)
      complex*16 gj(nrad,nuzp,kmymaxp),fj(nrad,nuzp,kmymaxp)
c
      complex*16 qmoma(lmsup,mintfc),qmomb(lmsup,mintfc)
      complex*16 zdos(kmymaxp),zrho(nrad),zmagrho(nrad)
      complex*16 zrhosp(nrad,2),zrhodsp(nrad,2)
      complex*16 zdmag(kmymaxp,3),zdmor(kmymaxp,3)
      complex*16 zqmom(kmymaxp,lmsup),zenorb
      complex*16 qmom(lmsup)
c
      complex*16 rgacoeff(kmymaxp,kmymaxp,lmsup)
      common/rggaunt/rgacoeff
c
      complex*16 sxcoeff(kmymaxp,kmymaxp),sxbcoeff(kmymaxp,kmymaxp)
      complex*16 sycoeff(kmymaxp,kmymaxp),sybcoeff(kmymaxp,kmymaxp)
      complex*16 szcoeff(kmymaxp,kmymaxp),szbcoeff(kmymaxp,kmymaxp)
      common/sigmat/sxcoeff,sxbcoeff,sycoeff,sybcoeff,szcoeff,szbcoeff
c
      complex*16 lxcoeff(kmymaxp,kmymaxp),lxbcoeff(kmymaxp,kmymaxp)
      complex*16 lycoeff(kmymaxp,kmymaxp),lybcoeff(kmymaxp,kmymaxp)
      complex*16 lzcoeff(kmymaxp,kmymaxp),lzbcoeff(kmymaxp,kmymaxp)
      common/lmat/lxcoeff,lxbcoeff,lycoeff,lybcoeff,lzcoeff,lzbcoeff
c
      complex*16 alphalkkr(0:lmaxp,minprc)
      complex*16 alpharkkr(0:lmaxp,minprc)
      complex*16 alphaintkkr(0:lmaxp,mintfc)
      common/scrpar/alphalkkr,alpharkkr,alphaintkkr
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
      c=274.072d0
      if(.not.wrel) then
        psq=ce+ce*ce/(c*c)
      else
        psq=ce
      end if
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
c ***************************************************
c * loop over layers to compute physical quantities *
c ***************************************************
c 
      qvifc=0.d0
      enbifc=0.d0
      do li=1,nintfc
         cpalay=(1.d0-conc(li)).gt.tiny
c
c Compute scattering solutions 
c        --------------------------------------------------------
         call wafu(ce,psq,lmax,idpota(li),v0,vra(1,li),bra(1,li),
     >             bopra(1,1,li),dx(li),ns(li),rs(li),
     >             tm,gz,fz,gj,fj,nuz,indz,1,sxa(li))
c        --------------------------------------------------------
c
c Density of multipole moments
         i=0
         do lam=0,lmaxs
         do mu=-lam,lam
           i=i+1
c        --------------------------------------------
           call moment(lmax,rs(li),dx(li),ns(li),
     >                 gz,fz,gj,fj,nuz,indz,
     >                 rgacoeff(1,1,i),lam,taua(1,1,li),
     >                 zqmom(1,i),lms,1)
c        --------------------------------------------
         end do
         end do
         do kmy=1,kmymax
            zdos(kmy)=zqmom(kmy,1)
         enddo
c
c Magnetic density of states
c        -------------------------------------------------
         call magnet(lmax,rs(li),dx(li),ns(li),
     >               gz,fz,gj,fj,nuz,indz,taua(1,1,li),
     >               sxcoeff,sxbcoeff,zdmag(1,1),lms,1)
         call magnet(lmax,rs(li),dx(li),ns(li),
     >               gz,fz,gj,fj,nuz,indz,taua(1,1,li),
     >               sycoeff,sybcoeff,zdmag(1,2),lms,1)
         call magnet(lmax,rs(li),dx(li),ns(li),
     >               gz,fz,gj,fj,nuz,indz,taua(1,1,li),
     >               szcoeff,szbcoeff,zdmag(1,3),lms,1)
         call magnet(lmax,rs(li),dx(li),ns(li),
     >               gz,fz,gj,fj,nuz,indz,taua(1,1,li),
     >               lxcoeff,lxbcoeff,zdmor(1,1),.true.,1)
         call magnet(lmax,rs(li),dx(li),ns(li),
     >               gz,fz,gj,fj,nuz,indz,taua(1,1,li),
     >               lycoeff,lybcoeff,zdmor(1,2),.true.,1)
         call magnet(lmax,rs(li),dx(li),ns(li),
     >               gz,fz,gj,fj,nuz,indz,taua(1,1,li),
     >               lzcoeff,lzbcoeff,zdmor(1,3),.true.,1)
c        -------------------------------------------------
c
         do k=1,kmymax
           dosa(k,li)=dimag(zdos(k))
           dosmaga(k,li)=dimag(zdmag(k,3))
           if(itest.ge.2)
     >     write(6,'(2d13.5,5x,2d13.5)') dosa(k,li),dosmaga(k,li)
         end do
c
c Radial distribution of charge and magnetic density
c        ------------------------------------------------
         call dens(lmax,rs(li),dx(li),ns(li),gz,fz,gj,fj,
     >             nuz,indz,taua(1,1,li),bopra(1,1,li),
     >             zrho,zrhosp,zrhodsp,zmagrho,zenorb)
c        ------------------------------------------------
c
c add contribution to contour integral
c
         qvdiffa(li)=0.d0
         enbdiffa(li)=0.d0
         qv=0.0d0
         do k=1,kmymax
            qv=qv+dimag(we*zdos(k))
            qvdiffa(li)=qvdiffa(li)+dimag(we*zdos(k))
            enbdiffa(li)=enbdiffa(li)+dimag(we*ce*zdos(k))
            qvpa(k,li)=qvpa(k,li)+dimag(we*zdos(k))
            qva(li)=qva(li)+dimag(we*zdos(k))
            do i=1,3
              spin_magvpa(k,li,i)=spin_magvpa(k,li,i)+
     >                            dimag(we*zdmag(k,i))
              spin_magva(li,i)=spin_magva(li,i)+dimag(we*zdmag(k,i))
              orb_magvpa(k,li,i)=orb_magvpa(k,li,i)+
     >                           dimag(we*zdmor(k,i))
              orb_magva(li,i)=orb_magva(li,i)+dimag(we*zdmor(k,i))
            end do
            enba(li)=enba(li)+dimag(we*ce*zdos(k))
         end do
         enorba(li)=enorba(li)+dimag(we*zenorb)
c
c        call enblloyd(enbll,we,gtaua(1,1,li),ptminva(1,1,li),
c    >                 tminv(1,1,li),alphaintkkr(0,li),lmax,lliter)
c        enba(li)=enba(li)+enbll
c        write(6,'(2d15.6)') enbll,enba(li)
c        
         call denblloyd(denbll,we,lmax,
     >                  gtaua(1,1,li),ptminva(1,1,li),
     >                  dmata(1,1,li),dmatpa(1,1,li),
     >                  ddpha(1,1,li),ddphpa(1,1,li))
         denba(li)=denba(li)+denbll
c
         do irad=1,ns(li)
            rhov(irad)=dimag(we*zrho(irad))
            rhova(irad,li)=rhova(irad,li)+dimag(we*zrho(irad))
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
            end do
         end do
         end do
c
c Rotate moments of charge density to the global frame of reference
c and add to the integral
         do i=1,lmmaxs
           do j=1,lmmaxs
             qmoma(i,li)=qmoma(i,li)+qmom(j)*rmatpa(j,i,li)
           end do
         end do
c
c Check charges
c        x=dlog(rs(li))-(ns(li)-1)*dx(li)
c        do irad=1,ns(li)
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
           write(6,'(''Eorb '',d13.5)')  enorba(li)
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
         if(cpalay) then
c+------------+
c+ BIG CPA IF +
c+------------+
c
c Compute scattering solutions 
c        --------------------------------------------------------
         call wafu(ce,psq,lmax,idpotb(li),v0,vrb(1,li),brb(1,li),
     >             boprb(1,1,li),dx(li),ns(li),rs(li),
     >             tm,gz,fz,gj,fj,nuz,indz,1,sxb(li))
c        --------------------------------------------------------
c
c Density of multipole moments
         i=0
         do lam=0,lmaxs
         do mu=-lam,lam
           i=i+1
c          -------------------------------------------
           call moment(lmax,rs(li),dx(li),ns(li),
     >                 gz,fz,gj,fj,nuz,indz,
     >                 rgacoeff(1,1,i),lam,taub(1,1,li),
     >                 zqmom(1,i),lms,1)
c          -------------------------------------------
         end do
         end do
         do kmy=1,kmymax
            zdos(kmy)=zqmom(kmy,1)
         enddo
c
c Magnetic density of states
c        -------------------------------------------------
         call magnet(lmax,rs(li),dx(li),ns(li),
     >               gz,fz,gj,fj,nuz,indz,taub(1,1,li),
     >               sxcoeff,sxbcoeff,zdmag(1,1),lms,1)
         call magnet(lmax,rs(li),dx(li),ns(li),
     >               gz,fz,gj,fj,nuz,indz,taub(1,1,li),
     >               sycoeff,sybcoeff,zdmag(1,2),lms,1)
         call magnet(lmax,rs(li),dx(li),ns(li),
     >               gz,fz,gj,fj,nuz,indz,taub(1,1,li),
     >               szcoeff,szbcoeff,zdmag(1,3),lms,1)
         call magnet(lmax,rs(li),dx(li),ns(li),
     >               gz,fz,gj,fj,nuz,indz,taub(1,1,li),
     >               lxcoeff,lxbcoeff,zdmor(1,1),.true.,1)
         call magnet(lmax,rs(li),dx(li),ns(li),
     >               gz,fz,gj,fj,nuz,indz,taub(1,1,li),
     >               lycoeff,lybcoeff,zdmor(1,2),.true.,1)
         call magnet(lmax,rs(li),dx(li),ns(li),
     >               gz,fz,gj,fj,nuz,indz,taub(1,1,li),
     >               lzcoeff,lzbcoeff,zdmor(1,3),.true.,1)
c        -------------------------------------------------
c
         do k=1,kmymax
           dosb(k,li)=dimag(zdos(k))
           dosmagb(k,li)=dimag(zdmag(k,3))
           if(itest.ge.2)
     >     write(6,'(2d13.5,5x,2d13.5)') dosb(k,li),dosmagb(k,li)
         end do
c
c Radial distribution of charge and magnetic density
c        ------------------------------------------------
         call dens(lmax,rs(li),dx(li),ns(li),gz,fz,gj,fj,
     >             nuz,indz,taub(1,1,li),boprb(1,1,li),
     >             zrho,zrhosp,zrhodsp,zmagrho,zenorb)
c        ------------------------------------------------
c
c add contribution to contour integral
c
         qvdiffb(li)=0.d0
         enbdiffb(li)=0.d0
         do k=1,kmymax
            qvdiffb(li)=qvdiffb(li)+dimag(we*zdos(k))
            enbdiffb(li)=enbdiffb(li)+dimag(we*ce*zdos(k))
            qvpb(k,li)=qvpb(k,li)+dimag(we*zdos(k))
            qvb(li)=qvb(li)+dimag(we*zdos(k))
            do i=1,3
              spin_magvpb(k,li,i)=spin_magvpb(k,li,i)+
     >                            dimag(we*zdmag(k,i))
              spin_magvb(li,i)=spin_magvb(li,i)+dimag(we*zdmag(k,i))
              orb_magvpb(k,li,i)=orb_magvpb(k,li,i)+
     >                           dimag(we*zdmor(k,i))
              orb_magvb(li,i)=orb_magvb(li,i)+dimag(we*zdmor(k,i))
            end do
            enbb(li)=enbb(li)+dimag(we*ce*zdos(k))
         end do
         enorbb(li)=enorbb(li)+dimag(we*zenorb)
c
c        call enblloyd(enbll,we,gtaub(1,1,li),ptminvb(1,1,li),
c    >                 tminv(1,1,li),alphaintkkr(0,li),lmax,lliter)
c        enbb(li)=enbb(li)+enbll
c        write(6,'(2d15.6)') enbll,enbb(li)
c        
         call denblloyd(denbll,we,lmax,
     >                  gtaub(1,1,li),ptminvb(1,1,li),
     >                  dmatb(1,1,li),dmatpb(1,1,li),
     >                  ddphb(1,1,li),ddphpb(1,1,li))
         denbb(li)=denbb(li)+denbll
c
         do irad=1,ns(li)
            rhovb(irad,li)=rhovb(irad,li)+dimag(we*zrho(irad))
            rhospb(irad,1,li)=rhospb(irad,1,li)+dimag(we*zrhosp(irad,1))
            rhospb(irad,2,li)=rhospb(irad,2,li)+dimag(we*zrhosp(irad,2))
            rhodspb(irad,1,li)=rhodspb(irad,1,li)+
     >                         dimag(we*zrhodsp(irad,1))
            rhodspb(irad,2,li)=rhodspb(irad,2,li)+
     >                         dimag(we*zrhodsp(irad,2))
            rhomagb(irad,li)=rhomagb(irad,li)+dimag(we*zmagrho(irad))
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
            end do
         end do
         end do
c
c Rotate moments of charge density to the global frame of reference
c and add to the integral
         do i=1,lmmaxs
           do j=1,lmmaxs
             qmomb(i,li)=qmomb(i,li)+qmom(j)*rmatpb(j,i,li)
           end do
         end do
c
         else
c
           do k=1,kmymax
             qvpb(k,li)=qvpa(k,li)
             do i=1,3
               spin_magvpb(k,li,i)=spin_magvpa(k,li,i)
               orb_magvpb(k,li,i)=orb_magvpa(k,li,i)
             end do
           end do
           do i=1,lmmaxs
             qmomb(i,li)=qmoma(i,li)
           end do
           do i=1,3
             spin_magvb(li,i)=spin_magva(li,i)
             orb_magvb(li,i)=orb_magva(li,i)
           end do
           qvdiffb(li)=qvdiffa(li)
           enbdiffb(li)=enbdiffa(li)
           qvb(li)=qva(li)
           enbb(li)=enba(li)
           enorbb(li)=enorba(li)
           do irad=1,ns(li)
             rhovb(irad,li)=rhova(irad,li)
             rhospb(irad,1,li)=rhospa(irad,1,li)
             rhospb(irad,2,li)=rhospa(irad,2,li)
             rhodspb(irad,1,li)=rhodspa(irad,1,li)
             rhodspb(irad,2,li)=rhodspa(irad,2,li)
             rhomagb(irad,li)=rhomaga(irad,li)
           end do
c
         end if
c
         if(itest.ge.3) then
           write(6,'('' Contribution to contour integral'')')
           write(6,'(''Q    '',d13.5)') qvb(li)
           do i=1,lmmaxs
             write(6,'(''Q'',i2,2x,2d13.5)') i,qmomb(i,li)
           end do
           write(6,'(''Eband'',d13.5)') enbb(li)
           write(6,'(''Eorb '',d13.5)')  enorbb(li)
           write(6,'(''Sx,Lx'',d13.5,5x,d13.5)') 
     >               spin_magvb(li,1),orb_magvb(li,1)
           write(6,'(''Sy,Ly'',d13.5,5x,d13.5)') 
     >               spin_magvb(li,2),orb_magvb(li,2)
           write(6,'(''Sz,Lz'',d13.5,5x,d13.5)') 
     >               spin_magvb(li,3),orb_magvb(li,3)
           if(itest.gt.3) then
             write(6,'(''Loop nrad - rho,rhomag'')')
             do irad=1,ns(li)
              write(6,'(d13.5,5x,d13.5)') 
     >              rhovb(irad,li),rhomagb(irad,li)
              write(6,'(d13.5,5x,d13.5)') 
     >              rhospb(irad,1,li),rhospb(irad,2,li)
              write(6,'(d13.5,5x,d13.5)') 
     >              rhodspb(irad,1,li),rhodspb(irad,2,li)
             enddo
           endif
         endif
c
c+----------------+
c+ END BIG CPA IF +
c+----------------+
c
         enbifc=enbifc+conc(li)*enbdiffa(li)
     >                +(1.0d0-conc(li))*enbdiffb(li)
         qvifc=qvifc+conc(li)*qvdiffa(li)
     >              +(1.0d0-conc(li))*qvdiffb(li)
c
      end do
c ************************
c * end loop over layers *
c ************************

      if(linbw) then
        omifc=enbifc-efermi*qvifc
        enbifcint=enbifcint+enbifc
        qvifcint=qvifcint+qvifc
        omifcint=omifcint+omifc
        write(6,'('' Enb'',i3,2f10.5,3d17.8)') 
     >  ie,ce,enbifc,qvifc,omifc
      end if
c
      return
      end
