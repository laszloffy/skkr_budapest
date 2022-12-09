c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine bsf(
     > ie,ce,we,dos,lmax,nintfc,wrel,lms,sxa,sxb,conc,v0,
     > E_Fermi,deltaa,deltab,
     > idpota,vra,bra,rba,idpotb,vrb,brb,rbb,dx,ns,rs,
     > tminv,tau_k,tau_kint,taua_kint,taub_kint,gtaua,gtaub,
     > bsfe,bsfh,bsfeh,bsfhe,
     > bsfteh,bshthe,
     > lliter,linbw,efermi,enbifcint,qvifcint,omifcint,reg,c_light)
c
      implicit real*8 (a-h,o-z)
      include '../param.h'
      logical wrel,lms,cpalay,linbw,dos,reg
      character*10 idpota(mintfc),idpotb(mintfc)
      dimension mpar(-lsup:lsup)
      dimension vra(nrad,mintfc),bra(nrad,mintfc)
      dimension vrb(nrad,mintfc),brb(nrad,mintfc)
      dimension rs(mintfc),dx(mintfc),ns(mintfc)
      dimension rba(3,mintfc),rbb(3,mintfc)
      complex*16 deltaa(nrad,mintfc),deltab(nrad,mintfc)
c
      dimension dosa(kmymaxp,mintfc),dosb(kmymaxp,mintfc)
      dimension dosha(kmymaxp,mintfc),doshb(kmymaxp,mintfc)
      dimension doseha(kmymaxp,mintfc),dosteha(kmymaxp,mintfc)
      dimension doshea(kmymaxp,mintfc),dosthea(kmymaxp,mintfc)
      dimension dosehb(kmymaxp,mintfc),dostehb(kmymaxp,mintfc)
      dimension dosheb(kmymaxp,mintfc),dostheb(kmymaxp,mintfc)
c
      dimension bsfe(mintfc),bsfh(mintfc),bsfeh(mintfc)
      dimension bsfhe(mintfc),bsfteh(mintfc),bsfthe(mintfc)
c
      dimension rad(nrad)
c
      dimension conc(mintfc),sxa(mintfc),sxb(mintfc)
c
      complex*16 ce,we,psq,cih
      complex*16 detl
c
c
      complex*16 tau_k(dbogomaxp,dbogomaxp,mintfc) ! tau_c(k)
      complex*16 tau_kint(dbogomaxp,dbogomaxp,mintfc) ! tau_c
      complex*16 tau_diff(dbogomaxp,dbogomaxp) ! tau_c(k)-tau_c
      complex*16 tauinv(dbogomaxp,dbogomaxp) ! tau_c^-1
      complex*16 taua_kint(dbogomaxp,dbogomaxp,mintfc) ! tau_a
      complex*16 taub_kint(dbogomaxp,dbogomaxp,mintfc) ! tau_b
      complex*16 taui(dbogomaxp,dbogomaxp,2) ! tau_alpha
      complex*16 dmat(dbogomaxp,dbogomaxp,2) ! tau_alpha*tauc^-1
      complex*16 dmatp(dbogomaxp,dbogomaxp,2) ! tauc^-1*tau_alpha
      complex*16 gtaua(dbogomaxp,kmymaxp,mintfc)
      complex*16 gtaub(dbogomaxp,dbogomaxp,mintfc)
      complex*16 tminv(dbogomaxp,dbogomaxp,mintfc)
      complex*16 tm(dbogomaxp,dbogomaxp)
      complex*16 tempmat(dbogomaxp,dbogomaxp)
c wave functions 
      complex*16 gz(dbogomaxp,dbogomaxp,nrad,2)
      complex*16 fz(dbogomaxp,dbogomaxp,nrad,2)
      complex*16 gj(dbogomaxp,dbogomaxp,nrad,2)
      complex*16 fj(dbogomaxp,dbogomaxp,nrad,2)
      complex*16 glz(dbogomaxp,dbogomaxp,nrad,2)
      complex*16 flz(dbogomaxp,dbogomaxp,nrad,2)
      complex*16 glj(dbogomaxp,dbogomaxp,nrad,2)
      complex*16 flj(dbogomaxp,dbogomaxp,nrad,2)
c Green's function matricies
      complex*16 ggrmat(dbogomaxp,dbogomaxp,0:lsup)
      complex*16 fgrmat(dbogomaxp,dbogomaxp,0:lsup)
      complex*16 ggrmat_tmp(dbogomaxp,dbogomaxp,0:lsup) ! only for lam=0
      complex*16 fgrmat_tmp(dbogomaxp,dbogomaxp,0:lsup) ! only for lam=0

      complex*16 tggrmat(kmymaxp,kmymaxp,0:lsup)
      complex*16 tfgrmat(kmymaxp,kmymaxp,0:lsup)
c
      complex*16 zdos(dbogomaxp),zrho(nrad),zmagrho(nrad)
      complex*16 zrhosp(nrad,2),zrhodsp(nrad,2)
      complex*16 zrhoeh(nrad),zrhohe(nrad)
      complex*16 zdmag(dbogomaxp,3),zdmor(dbogomaxp,3)
      complex*16 zqmom(dbogomaxp,lmsup)
      complex*16 zscmom(dbogomaxp,lmsup),zsctmom(dbogomaxp,lmsup)
      complex*16 zscdos(dbogomaxp),zsctdos(dbogomaxp)
c
c
      real*8     rz(3),rtmp(3),rtmph(3),cab(2),cij
      integer kmyx2
c
      complex*16 rgacoeff(kmymaxp,kmymaxp,lmsup)
      common/rggaunt/rgacoeff
c
      common/sigmat/sxcoeff,sxbcoeff,sycoeff,sybcoeff,szcoeff,szbcoeff
c
      complex*16 deltakmy(kmymaxp,kmymaxp)
      common/lmat/lxcoeff,lxbcoeff,lycoeff,lybcoeff,lzcoeff,lzbcoeff
c
c
      common/test/itest
c
      data tol/1.0d-8/,cih/(0.0d0,-0.5d0)/
      data tiny/1.0d-6/
c
c      write(6,*) '> Subroutine bsf'
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
      kmyx2=2*kmymax
      do li=1,nintfc
        cpalay=(1.d0-conc(li)).gt.tiny
c This way we don't sum up the contributions of different layers!       
        tempmat = (0.d0,0.d0)
        ggrmat = (0.d0,0.d0)
        fgrmat = (0.d0,0.d0)
c
      taui=(0.d0,0.d0)
      if (cpalay) then
        cab(1) = conc(li)
        cab(2) = 1.d0-conc(li)
        taui(1:2*kmymax,1:2*kmymax,1) = 
     >  taua_kint(1:2*kmymax,1:2*kmymax,li) 
        taui(1:2*kmymax,1:2*kmymax,2) = 
     >  taub_kint(1:2*kmymax,1:2*kmymax,li) 
c        write(6,*) 'taui is uploaded with taua/b'
      else
        write(6,*) 'we have no cpa in bsf'
        cab(1) = 1.d0
        cab(2) = 0.d0
      end if
       call cpu_time(tstart)
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
     >            gz(:,:,:,1),fz(:,:,:,1),gj(:,:,:,1),fj(:,:,:,1),
     >            glz(:,:,:,1),flz(:,:,:,1),glj(:,:,:,1),flj(:,:,:,1),1)

        if (cpalay) then 
        call wafu(ce,c,sxb(li),lmax,idpotb(li),v0,E_Fermi,
     >            deltab(1,li),vrb(1,li),brb(1,li),
     >            rtmp,dx(li),rs(li),ns(li),tm,
     >            gz(:,:,:,2),fz(:,:,:,2),gj(:,:,:,2),fj(:,:,:,2),
     >            glz(:,:,:,2),flz(:,:,:,2),glj(:,:,:,2),flj(:,:,:,2),1)
c        write(6,*) 'delta B for li=',li
c        write(6,*) deltab(:,li)
        end if 
c       --------------------------------------------------------
       call cpu_time(tfinish)
       write(6,'('' Time ellapsed in wafu:'',t35,f8.1,'' s'')')
     >  tfinish-tstart

c       call cpu_time(tstart)
c
c=================================
c--- Calculation of D matrices ---
c=================================

        call cpu_time(tstart)
c we need tauinv and D matrices only for cpa
        if (cpalay) then
        dmat = (0.d0,0.d0)
        dmatp = (0.d0,0.d0)
        tauinv = (0.d0,0.d0)
        tauinv(1:2*kmymax,1:2*kmymax)=tau_kint(1:2*kmymax,1:2*kmymax,li)
        call gjinv(tauinv(:,:),2*kmymax,dbogomaxp,detl)
        call doubmt1(taua_kint(:,:,li),tauinv,dmat(:,:,1),
     >       2*kmymax,dbogomaxp)
        call doubmt1(tauinv,taua_kint(:,:,li),dmatp(:,:,1),
     >       2*kmymax,dbogomaxp) 
c        write(6,*) 'dmat and dmatp for component A are calculated'
c        call outmat1(dmat(:,:,1),2*kmymax,2*kmymax,dbogomaxp,tol,6)
c        call outmat1(dmatp(:,:,1),2*kmymax,2*kmymax,dbogomaxp,tol,6)
c  We calculate for the B component only in the case of cpa 
        call doubmt1(taub_kint(:,:,li),tauinv,dmat(:,:,2),
     >       2*kmymax,dbogomaxp)
        call doubmt1(tauinv,taub_kint(:,:,li),dmatp(:,:,2),
     >       2*kmymax,dbogomaxp)
c        write(6,*) 'dmat and dmatp for component B are calculated'
c        call outmat1(dmat(:,:,2),2*kmymax,2*kmymax,dbogomaxp,tol,6)
c        call outmat1(dmatp(:,:,2),2*kmymax,2*kmymax,dbogomaxp,tol,6)
        end if 
c
c========================================
c--- Calculation of ggrmat and fgrmat ---
c========================================
c       
        tau_diff=(0.d0,0.d0)
c If we have cpa we need taudiff, else taudiff is simply tau_c(k)
        if (cpalay) then
          call submat1(tau_k(:,:,li),tau_kint(:,:,li),tau_diff,
     >         2*kmymax,dbogomaxp)
        else 
          tau_diff(:,:)=tau_k(:,:,li)
        end if
c
c       write(6,*) 'tau_c(k) in bsf for li= '
c       call outmat1(tau_k(:,:,li),2*kmymax,2*kmymax,dbogomaxp,tol,6)
c       write(6,*) 'tau_kint in bsf for li= '
c       call outmat1(tau_kint(:,:,li),2*kmymax,2*kmymax,dbogomaxp,tol,6)
c       write(6,*) 'tau_diff  in bsf for li= '
c       call outmat1(tau_diff(:,:),2*kmymax,2*kmymax,dbogomaxp,tol,6)
       do icpa=1,2
c      -------------------------
c      --- Calculation of A0 ---
c      -------------------------
         if (cpalay) then
c         write(6,*) 'we have single sum contribution for', icpa
            call gf(lmax,reg,dx(li),ns(li),rs(li),
     >      gz(:,:,:,icpa),fz(:,:,:,icpa),gj(:,:,:,icpa),fj(:,:,:,icpa),
     >      glz(:,:,:,icpa),flz(:,:,:,icpa),glj(:,:,:,icpa),flj(:,:,:,icpa),
     >      taui(:,:,icpa),ggrmat_tmp,fgrmat_tmp)
c           write(6,*) 'After gf for A0'
            call addmatc(ggrmat(:,:,0),ggrmat_tmp(:,:,0),cab(icpa),
     >           2*kmymax,dbogomaxp)
            call addmatc(fgrmat(:,:,0),fgrmat_tmp(:,:,0),cab(icpa),
     >           2*kmymax,dbogomaxp)
c           write(6,*) 'After addmatc for A0'
         end if 
       do jcpa=1,2
c      -------------------------
c      --- Calculation of A1 ---
c      -------------------------
       cij = cab(icpa)*cab(jcpa)
c      we do something only if we have nonzero contribution
         if (cij.gt.tiny) then
c         write(6,*)'we have double sum contribution for ',icpa,jcpa 
c      tempmat = D_alpha * (tau_c(k)-tau_c) * D_beta
          if (cpalay) then
          call tripmt1(dmat(:,:,icpa),tau_diff,dmatp(:,:,jcpa),tempmat,
     >        2*kmymax,2*kmymax,dbogomaxp)
          else
              tempmat(:,:) = tau_diff(:,:)
          end if
c       write(6,*) 'D_alpha * (tau_c(k)-tau_c) * D_beta calculated'
c
          call gf(lmax,reg,dx(li),ns(li),rs(li),
     >        gz(:,:,:,icpa),fz(:,:,:,icpa),gj(:,:,:,icpa),fj(:,:,:,icpa),
     >        glz(:,:,:,jcpa),flz(:,:,:,jcpa),glj(:,:,:,jcpa),flj(:,:,:,jcpa),
     >        tempmat,ggrmat_tmp,fgrmat_tmp)
c        write(6,*) 'After gf for A1'
c        write(6,*) 'cij for i,j',icpa,jcpa,cij
c        write(6,*) 'ggrmat_tmp'
c        call outmat1(ggrmat_tmp(:,:,0),2*kmymax,2*kmymax,dbogomaxp,tol,6)
c        write(6,*) 'fgrmat_tmp'
c        call outmat1(fgrmat_tmp(:,:,0),2*kmymax,2*kmymax,dbogomaxp,tol,6)

          call addmatc(ggrmat(:,:,0),ggrmat_tmp(:,:,0),cij,
     >         2*kmymax,dbogomaxp)
          call addmatc(fgrmat(:,:,0),fgrmat_tmp(:,:,0),cij,
     >         2*kmymax,dbogomaxp)
         end if
       end do
       end do
c===============================
c===  End of double CPA loop ===
c===============================
       call cpu_time(tfinish)
c       write(6,'('' Time ellapsed in bsf cpa loops:'',
c     >              t35,f8.1,'' s'')')
c     >  tfinish-tstart
c       write(6,*) 'ggrmat in bsf for li= '
c       call outmat1(ggrmat(:,:,0),2*kmymax,2*kmymax,dbogomaxp,tol,6)
c       write(6,*) 'fgrmat in bsf for li= '
c       call outmat1(fgrmat(:,:,0),2*kmymax,2*kmymax,dbogomaxp,tol,6)
c============================================
c
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
        i=1
        lam=0
c        do lam=0,lmaxs
c        do mu=-lam,lam
c          i=i+1
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
c        end do
c        end do
        do kmy=1,2*kmymax
          zdos(kmy)=zqmom(kmy,1)
          zscdos(kmy)=zscmom(kmy,1)
          zsctdos(kmy)=zsctmom(kmy,1)  
        enddo
c
c
c we sum up over the kappa my index
        do k=1,kmymax
          dosa(k,li)=dimag(zdos(k))                  ! electron part
          bsfe(li) = bsfe(li) + dosa(k,li)
          dosha(k,li)=dimag(zdos(k+kmymax))          ! hole part
          bsfh(li) = bsfh(li) + dosha(k,li)
          doseha(k,li)=dimag(zscdos(k))               ! electron-hole part 
          bsfeh(li) = bsfeh(li) + doseha(k,li)
          doshea(k,li)=dimag(zscdos(k+kmymax))        ! triplet part 
          bsfhe(li) = bsfhe(li) + doshea(k,li)
          dosteha(k,li)=dimag(zsctdos(k))               ! triplet up part 
          bsfteh(li) = bsfteh(li) + dosteha(k,li)
          dosthea(k,li)=dimag(zsctdos(k+kmymax))        ! triplet down part 
          bsfthe(li) = bsfthe(li) + dosthea(k,li)
          if(itest.ge.2) then
            write(6,*) 'DOS electron kmy layeri'
            write(6,*) dosa(k,li),k,li
            write(6,*) 'DOS hole kmy layeri'
            write(6,*) dosha(k,li),k,li
          end if
        end do
        write(6,'(''   li='',i3,''  bsf='',d14.6,''  bsfh='',d14.6,''
     >             bfs s='',d14.6)' ) 
c        write(6,*)
     >        li,bsfe(li),bsfh(li),bsfeh(li)
c
c
      end do
c ************************
c * end loop over layers *
c ************************
      return
      end
