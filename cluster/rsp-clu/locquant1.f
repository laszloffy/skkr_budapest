c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine locquant1(
c     ===================
     > ie,ce,we,lmax,madmax,nintfc,wrel,lms,sxa,v0,
c
     > idpota,vra,bra,bopra,dx,ns,rs,nimp,nposimp,
c
     > dmata,dmatpa,
     > ddpha,ddphpa,
     > ddtha,ddthpa,
     > rmata,rmatpa,
c
     > taua,gtaua,ptminva,
c
     > dosa,qvpa,qva,qvdiffa,
     > vmadid,vmadiq,
     > enba,enbdiffa,denba,enorba,qmoma,
     > rhova,rhospa,rhodspa,rhomaga,
     > dosmaga,spin_magvpa,spin_magva,orb_magvpa,orb_magva,
c
     > lliter,efermi,enbifcint,qvifcint,omifcint)
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      integer mmax
      parameter (mmax=mimp)
c
      logical tautest
      logical wrel,lms
c
      character*10 idpota(mmax)
      character*30 for006
c
      integer nuz(kmymaxp)
      integer indz(nuzp,kmymaxp)
      integer mpar(-lsup:lsup)
      integer ns(mmax)
c
      real*8 vra(nrad,mmax)
      real*8 bra(nrad,mmax)
      real*8 bopra(nrad,2,mmax)
      real*8 rs(mmax)
      real*8 dx(mmax)
c
      real*8 enba(mmax)
      real*8 denba(mmax)
      real*8 enorba(mmax)
c
      real*8 dosa(kmymaxp,mmax)
      real*8 dosmaga(kmymaxp,mmax)
      real*8 qva(mmax)
      real*8 qvpa(kmymaxp,mmax)
      real*8 spin_magvpa(kmymaxp,mmax,3)
      real*8 spin_magva(mmax,3)
      real*8 orb_magvpa(kmymaxp,mmax,3)
      real*8 orb_magva(mmax,3)
c
      real*8 enbdiffa(mmax)
      real*8 qvdiffa(mmax)
c
      real*8 rhova(nrad,mmax)
      real*8 rhospa(nrad,2,mmax)
      real*8 rhodspa(nrad,2,mmax)
      real*8 rhomaga(nrad,mmax)
      real*8 rhov(nrad)
      real*8 rad(nrad)
c
      real*8 sxa(mmax)
c
      real*8 vmadid(mmax)
      real*8 vmadiq(mmax)
c
      complex*16 ce
      complex*16 we
      complex*16 psq
      complex*16 cih
      complex*16 temp
c
c ======================================================================
      complex*16 mlmij(lmsup,mimp)
c
      integer nposimp(3,mimp)
      integer madmax
      integer nimp
c
c ======================================================================
c
      complex*16 dmata(kmymaxp,kmymaxp,mmax)
      complex*16 dmatpa(kmymaxp,kmymaxp,mmax)
      complex*16 ddpha(kmymaxp,kmymaxp,mmax)
      complex*16 ddphpa(kmymaxp,kmymaxp,mmax)
      complex*16 ddtha(kmymaxp,kmymaxp,mmax)
      complex*16 ddthpa(kmymaxp,kmymaxp,mmax)
      complex*16 rmata(lmsup,lmsup,mmax)
      complex*16 rmatpa(lmsup,lmsup,mmax)
c
      complex*16 taua(kmymaxp,kmymaxp,mmax)
      complex*16 gtaua(kmymaxp,kmymaxp,mmax)
      complex*16 ptminva(kmymaxp,kmymaxp,mmax)
      complex*16 tm(kmymaxp,kmymaxp)
c
      complex*16 gz(nrad,nuzp,kmymaxp),fz(nrad,nuzp,kmymaxp)
      complex*16 gj(nrad,nuzp,kmymaxp),fj(nrad,nuzp,kmymaxp)
c
      complex*16 qmoma(lmsup,mmax)
      complex*16 zdos(kmymaxp)
      complex*16 zrho(nrad)
      complex*16 zmagrho(nrad)
      complex*16 zrhosp(nrad,2)
      complex*16 zrhodsp(nrad,2)
      complex*16 zdmag(kmymaxp,3)
      complex*16 zdmor(kmymaxp,3)
      complex*16 zqmom(kmymaxp,lmsup)
      complex*16 zenorb
      complex*16 qmom(lmsup)
      complex*16 qmomrot(lmsup)
c
      complex*16 rgacoeff(kmymaxp,kmymaxp,lmsup)
      common/rggaunt/rgacoeff
c
      complex*16 sxcoeff(kmymaxp,kmymaxp)
      complex*16 sxbcoeff(kmymaxp,kmymaxp)
      complex*16 sycoeff(kmymaxp,kmymaxp)
      complex*16 sybcoeff(kmymaxp,kmymaxp)
      complex*16 szcoeff(kmymaxp,kmymaxp)
      complex*16 szbcoeff(kmymaxp,kmymaxp)
      common/sigmat/sxcoeff,sxbcoeff,sycoeff,sybcoeff,szcoeff,szbcoeff
c
      complex*16 lxcoeff(kmymaxp,kmymaxp)
      complex*16 lxbcoeff(kmymaxp,kmymaxp)
      complex*16 lycoeff(kmymaxp,kmymaxp)
      complex*16 lybcoeff(kmymaxp,kmymaxp)
      complex*16 lzcoeff(kmymaxp,kmymaxp)
      complex*16 lzbcoeff(kmymaxp,kmymaxp)
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
!        cpalay=(1.d0-conc(li)).gt.tiny
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
         tautest=.false.
         write(6,*) '<locquant> : tautest dosa li=',li
         do k=1,kmymax
           dosa(k,li)=dimag(zdos(k))
           dosmaga(k,li)=dimag(zdmag(k,3))
           if(itest.ge.2.or.tautest) then
             write(6,*) 'DOS electron kmy layeri'
             write(6,*) dosa(k,li),k,li
           end if
c    >     write(6,'(2d13.5,5x,2d13.5)') dosa(k,li),dosmaga(k,li)
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
!           rhov(irad)=dimag(we*zrho(irad))
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
           qmomrot(i)=(0.d0,0.d0)
           do j=1,lmmaxs
             qmoma(i,li)=qmoma(i,li)+qmom(j)*rmatpa(j,i,li)
             qmomrot(i)=qmomrot(i)+qmom(j)*rmatpa(j,i,li)
           end do
         end do
c
c ======================================================================
c *******************************************************
c * Madelung potential due to the impurities wihtout 00 *
c *******************************************************
         call madcoimp(madmax,li,nimp,nposimp,mlmij)
         do jimp=1,nintfc
           if(jimp.ne.li) then
c  from dipole momentum
c
             temp=(0.d0,0.d0) !DEBUG
             do ii=2,4
               vmadid(jimp)=vmadid(jimp)+
     >                      dreal(mlmij(ii,jimp)*qmomrot(ii))
               temp=temp+mlmij(ii,jimp)*qmomrot(ii)
             end do
             if(dimag(temp).gt.1.0d-12) then
                write(6,*) '<locquant1>: *********************'
                write(6,*) '<locquant1>: IMAGINARY PART OF THE'
                write(6,*) '<locquant1>: MADELUNG DIP .NE. 0  '
                write(6,*) '<locquant1>: temp=',jimp,temp !DEBUG
                write(6,*) '<locquant1>: STOP DEBUG'
                write(6,*) '<locquant1>: *********************'
             end if
c
c  from quadrupol momentum
c
             temp=(0.d0,0.d0) !DEBUG
             do ii=5,9
               vmadiq(jimp)=vmadiq(jimp)+
     >                      dreal(mlmij(ii,jimp)*qmomrot(ii))
               temp=temp+mlmij(ii,jimp)*qmomrot(ii)
             end do
             if(dimag(temp).gt.1.0d-12) then
                write(6,*) '<locquant1>: *********************'
                write(6,*) '<locquant1>: IMAGINARY PART OF THE'
                write(6,*) '<locquant1>: MADELUNG QUAD .NE. 0 '
                write(6,*) '<locquant1>: temp=',jimp,temp !DEBUG
                write(6,*) '<locquant1>: STOP DEBUG'
                write(6,*) '<locquant1>: *********************'
             end if
c
           end if
         end do
c ======================================================================
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
c        if(cpalay) then
c+------------+
c+ BIG CPA IF +
c+------------+
c
c Compute scattering solutions 
c        --------------------------------------------------------
c        call wafu(ce,psq,lmax,idpotb(li),v0,vrb(1,li),brb(1,li),
c    >             boprb(1,1,li),dx(li),ns(li),rs(li),
c    >             tm,gz,fz,gj,fj,nuz,indz,1,sxb(li))
c        --------------------------------------------------------
c
c Density of multipole moments
c        i=0
c        do lam=0,lmaxs
c        do mu=-lam,lam
c          i=i+1
c          -------------------------------------------
c          call moment(lmax,rs(li),dx(li),ns(li),
c    >                 gz,fz,gj,fj,nuz,indz,
c    >                 rgacoeff(1,1,i),lam,taub(1,1,li),
c    >                 zqmom(1,i),lms,1)
c          -------------------------------------------
c        end do
c        end do
c        do kmy=1,kmymax
c           zdos(kmy)=zqmom(kmy,1)
c        enddo
c
c Magnetic density of states
c        -------------------------------------------------
c        call magnet(lmax,rs(li),dx(li),ns(li),
c    >               gz,fz,gj,fj,nuz,indz,taub(1,1,li),
c    >               sxcoeff,sxbcoeff,zdmag(1,1),lms,1)
c        call magnet(lmax,rs(li),dx(li),ns(li),
c    >               gz,fz,gj,fj,nuz,indz,taub(1,1,li),
c    >               sycoeff,sybcoeff,zdmag(1,2),lms,1)
c        call magnet(lmax,rs(li),dx(li),ns(li),
c    >               gz,fz,gj,fj,nuz,indz,taub(1,1,li),
c    >               szcoeff,szbcoeff,zdmag(1,3),lms,1)
c        call magnet(lmax,rs(li),dx(li),ns(li),
c    >               gz,fz,gj,fj,nuz,indz,taub(1,1,li),
c    >               lxcoeff,lxbcoeff,zdmor(1,1),.true.,1)
c        call magnet(lmax,rs(li),dx(li),ns(li),
c    >               gz,fz,gj,fj,nuz,indz,taub(1,1,li),
c    >               lycoeff,lybcoeff,zdmor(1,2),.true.,1)
c        call magnet(lmax,rs(li),dx(li),ns(li),
c    >               gz,fz,gj,fj,nuz,indz,taub(1,1,li),
c    >               lzcoeff,lzbcoeff,zdmor(1,3),.true.,1)
c        -------------------------------------------------
c
c        do k=1,kmymax
c          dosb(k,li)=dimag(zdos(k))
c          dosmagb(k,li)=dimag(zdmag(k,3))
c          if(itest.ge.2)
c    >     write(6,'(2d13.5,5x,2d13.5)') dosb(k,li),dosmagb(k,li)
c        end do
c
c Radial distribution of charge and magnetic density
c        ------------------------------------------------
c        call dens(lmax,rs(li),dx(li),ns(li),gz,fz,gj,fj,
c    >             nuz,indz,taub(1,1,li),boprb(1,1,li),
c    >             zrho,zrhosp,zrhodsp,zmagrho,zenorb)
c        ------------------------------------------------
c
c add contribution to contour integral
c
c        qvdiffb(li)=0.d0
c        enbdiffb(li)=0.d0
c        do k=1,kmymax
c           qvdiffb(li)=qvdiffb(li)+dimag(we*zdos(k))
c           enbdiffb(li)=enbdiffb(li)+dimag(we*ce*zdos(k))
c           qvpb(k,li)=qvpb(k,li)+dimag(we*zdos(k))
c           qvb(li)=qvb(li)+dimag(we*zdos(k))
c           do i=1,3
c             spin_magvpb(k,li,i)=spin_magvpb(k,li,i)+
c    >                            dimag(we*zdmag(k,i))
c             spin_magvb(li,i)=spin_magvb(li,i)+dimag(we*zdmag(k,i))
c             orb_magvpb(k,li,i)=orb_magvpb(k,li,i)+
c    >                           dimag(we*zdmor(k,i))
c             orb_magvb(li,i)=orb_magvb(li,i)+dimag(we*zdmor(k,i))
c           end do
c           enbb(li)=enbb(li)+dimag(we*ce*zdos(k))
c        end do
c        enorbb(li)=enorbb(li)+dimag(we*zenorb)
c
c        call enblloyd(enbll,we,gtaub(1,1,li),ptminvb(1,1,li),
c    >                 tminv(1,1,li),alphaintkkr(0,li),lmax,lliter)
c        enbb(li)=enbb(li)+enbll
c        write(6,'(2d15.6)') enbll,enbb(li)
c        
c        call denblloyd(denbll,we,lmax,
c    >                  gtaub(1,1,li),ptminvb(1,1,li),
c    >                  dmatb(1,1,li),dmatpb(1,1,li),
c    >                  ddphb(1,1,li),ddphpb(1,1,li))
c        denbb(li)=denbb(li)+denbll
c
c        do irad=1,ns(li)
c           rhovb(irad,li)=rhovb(irad,li)+dimag(we*zrho(irad))
c           rhospb(irad,1,li)=rhospb(irad,1,li)+dimag(we*zrhosp(irad,1))
c           rhospb(irad,2,li)=rhospb(irad,2,li)+dimag(we*zrhosp(irad,2))
c           rhodspb(irad,1,li)=rhodspb(irad,1,li)+
c    >                         dimag(we*zrhodsp(irad,1))
c           rhodspb(irad,2,li)=rhodspb(irad,2,li)+
c    >                         dimag(we*zrhodsp(irad,2))
c           rhomagb(irad,li)=rhomagb(irad,li)+dimag(we*zmagrho(irad))
c        end do
c
c        do lam=0,lmaxs
c        do mu=-lam,lam
c           i=lam*(lam+1)+mu+1
c           ii=lam*(lam+1)-mu+1
c           qmom(i)=(0.0d0,0.0d0)
c           do k=1,kmymax
c             qmom(i)=qmom(i)+cih*(we*zqmom(k,i)-
c    >                             mpar(mu)*dconjg(we*zqmom(k,ii)))
c           end do
c        end do
c        end do
c
c Rotate moments of charge density to the global frame of reference
c and add to the integral
c        do i=1,lmmaxs
c          do j=1,lmmaxs
c            qmomb(i,li)=qmomb(i,li)+qmom(j)*rmatpb(j,i,li)
c          end do
c        end do
c
c        else
c
c          do k=1,kmymax
c            qvpb(k,li)=qvpa(k,li)
c            do i=1,3
c              spin_magvpb(k,li,i)=spin_magvpa(k,li,i)
c              orb_magvpb(k,li,i)=orb_magvpa(k,li,i)
c            end do
c          end do
c          do i=1,lmmaxs
c            qmomb(i,li)=qmoma(i,li)
c          end do
c          do i=1,3
c            spin_magvb(li,i)=spin_magva(li,i)
c            orb_magvb(li,i)=orb_magva(li,i)
c          end do
c          qvdiffb(li)=qvdiffa(li)
c          enbdiffb(li)=enbdiffa(li)
c          qvb(li)=qva(li)
c          enbb(li)=enba(li)
c          enorbb(li)=enorba(li)
c          do irad=1,ns(li)
c            rhovb(irad,li)=rhova(irad,li)
c            rhospb(irad,1,li)=rhospa(irad,1,li)
c            rhospb(irad,2,li)=rhospa(irad,2,li)
c            rhodspb(irad,1,li)=rhodspa(irad,1,li)
c            rhodspb(irad,2,li)=rhodspa(irad,2,li)
c            rhomagb(irad,li)=rhomaga(irad,li)
c          end do
c
c        end if
c
c        if(itest.ge.3) then
c          write(6,'('' Contribution to contour integral'')')
c          write(6,'(''Q    '',d13.5)') qvb(li)
c          do i=1,lmmaxs
c            write(6,'(''Q'',i2,2x,2d13.5)') i,qmomb(i,li)
c          end do
c          write(6,'(''Eband'',d13.5)') enbb(li)
c          write(6,'(''Eorb '',d13.5)')  enorbb(li)
c          write(6,'(''Sx,Lx'',d13.5,5x,d13.5)') 
c    >               spin_magvb(li,1),orb_magvb(li,1)
c          write(6,'(''Sy,Ly'',d13.5,5x,d13.5)') 
c    >               spin_magvb(li,2),orb_magvb(li,2)
c          write(6,'(''Sz,Lz'',d13.5,5x,d13.5)') 
c    >               spin_magvb(li,3),orb_magvb(li,3)
c          if(itest.gt.3) then
c            write(6,'(''Loop nrad - rho,rhomag'')')
c            do irad=1,ns(li)
c             write(6,'(d13.5,5x,d13.5)') 
c    >              rhovb(irad,li),rhomagb(irad,li)
c             write(6,'(d13.5,5x,d13.5)') 
c    >              rhospb(irad,1,li),rhospb(irad,2,li)
c             write(6,'(d13.5,5x,d13.5)') 
c    >              rhodspb(irad,1,li),rhodspb(irad,2,li)
c            enddo
c          endif
c        endif
c
c+----------------+
c+ END BIG CPA IF +
c+----------------+
c
         enbifc=enbifc+enbdiffa(li)
c        enbifc=enbifc+conc(li)*enbdiffa(li)
c    >                +(1.0d0-conc(li))*enbdiffb(li)
         qvifc=qvifc+qvdiffa(li)
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
