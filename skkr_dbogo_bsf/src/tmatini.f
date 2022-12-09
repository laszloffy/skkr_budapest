c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine tmatini(
c     ==================
     > ce,lmax,nintfc,iscreen,vscreen,v0,E_Fermi,
     > bulk,wrel,sxl,sxr,sxa,sxb,
     > concl,concr,conc,cpamatin,cpamatinl,cpamatinr,
c
     > idpotla,vrla,brla,idpotlb,vrlb,brlb,rsl,dxl,nsl,
     > idpotra,vrra,brra,idpotrb,vrrb,brrb,rsr,dxr,nsr,
     > idpota,vra,bra,idpotb,vrb,brb,dx,ns,rs,
c
     > dmatl,dmatlp,dmatr,dmatrp,dmata,dmatb,dmatpa,dmatpb,
c
     > tminvl,tminvr,tminv,tminva,tminvb,ptminva,ptminvb,
     > rbl,rbr,rba,rbb,
     > deltala,deltara,deltaa,
     > deltalb,deltarb,deltab,c_light)
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      logical wrel,bulk,cpalay
      logical cpamatin,cpamatinl,cpamatinr
c
      character*10 idpota(mintfc),idpotla(minprc),idpotra(minprc)
      character*10 idpotb(mintfc),idpotlb(minprc),idpotrb(minprc)
c
      complex*16 deltala(nrad,minprc),deltara(nrad,minprc)
      complex*16 deltaa(nrad,mintfc)
      complex*16 deltalb(nrad,minprc),deltarb(nrad,minprc)
      complex*16 deltab(nrad,mintfc)
      dimension vra(nrad,mintfc),bra(nrad,mintfc)
      dimension vrb(nrad,mintfc),brb(nrad,mintfc)
      dimension rs(mintfc),dx(mintfc),ns(mintfc)
      dimension vrla(nrad,minprc),brla(nrad,minprc)
      dimension vrlb(nrad,minprc),brlb(nrad,minprc)
      dimension rsl(minprc),dxl(minprc),nsl(minprc)
      dimension vrra(nrad,minprc),brra(nrad,minprc)
      dimension vrrb(nrad,minprc),brrb(nrad,minprc)
      dimension rsr(minprc),dxr(minprc),nsr(minprc)
c
      dimension conc(mintfc),sxa(mintfc),sxb(mintfc)
      dimension concl(minprc),sxl(minprc)
      dimension concr(minprc),sxr(minprc)
c
      complex*16 ce,we,psq
!
      real*8     rbl(3,minprc), rbr(3,minprc)
      real*8     rba(3,mintfc), rbb(3,mintfc)
c
      complex*16 tx(dbogomaxp,dbogomaxp)
      complex*16 tminv(dbogomaxp,dbogomaxp,mintfc)
      complex*16 tminva(dbogomaxp,dbogomaxp,mintfc)
      complex*16 tminvb(dbogomaxp,dbogomaxp,mintfc)
      complex*16 ptminva(dbogomaxp,dbogomaxp,mintfc)
      complex*16 ptminvb(dbogomaxp,dbogomaxp,mintfc)
      complex*16 tminvl(dbogomaxp,dbogomaxp,minprc)
      complex*16 tminvr(dbogomaxp,dbogomaxp,minprc)
c
      complex*16 dmatl(kmymaxp,kmymaxp,minprc)
      complex*16 dmatlp(kmymaxp,kmymaxp,minprc)
      complex*16 dmatr(kmymaxp,kmymaxp,minprc)
      complex*16 dmatrp(kmymaxp,kmymaxp,minprc)
      complex*16 dmata(kmymaxp,kmymaxp,mintfc)
      complex*16 dmatpa(kmymaxp,kmymaxp,mintfc)
      complex*16 dmatb(kmymaxp,kmymaxp,mintfc)
      complex*16 dmatpb(kmymaxp,kmymaxp,mintfc)
c
      complex*16 alphalkkr(0:lmaxp,minprc)
      complex*16 alpharkkr(0:lmaxp,minprc)
      complex*16 alphaintkkr(0:lmaxp,mintfc)
      complex*16 alphalkkrh(0:lmaxp,minprc)
      complex*16 alpharkkrh(0:lmaxp,minprc)
      complex*16 alphaintkkrh(0:lmaxp,mintfc)
c
      real*8  rbtmp(3), rz(3)
c
      common/scrpar/alphalkkr,alpharkkr,alphaintkkr
      common/scrparh/alphalkkrh,alpharkkrh,alphaintkkrh
c
      common/test/itest
      common/lay2d/cvec(mtotal,3),nextra,nbulkl,nbulkr,
     &             nprc,ninprc(0:mprc+1)
c
      data tol/1.0d-8/ 
      data tiny/1.0d-6/
c
c ********************
c initialize constants
c ********************
c
c---> c in rydberg units:
c      c=274.072d0 -- from input
      c=c_light
         
c      psq=ce+ce*ce/(c*c) !!!!!!!!!!!!! attettem az alphamatba
      psq=ce
c
      nl=lmax+1
      nl2=nl*nl
      kmax=2*lmax+1
      kmymax=2*nl2
c
      ninprcl = ninprc(0)
      ninprcr = ninprc(nprc+1)
c
c  determine whether we solve the Dirac equation in global frame
c  if we do: tminv, ptminv are in global frame -> rotate back ptminv to local
c  if we don't: tminv, ptminv are in local frame -> rotate tminv to global 
      !globalmode = .true.  ! <-- this has been moved to param.h as .not.localmode!!!
      rz = 0.d0
      rz(3) = 1.d0
c
c ********************
c screening parameters
c ********************
c
c zero.f is real*8 => double arrays size!
      call czero(alphalkkr,minprc*(lmaxp+1))
      call czero(alpharkkr,minprc*(lmaxp+1))
      call czero(alphaintkkr,mintfc*(lmaxp+1))
      call czero(alphalkkrh,minprc*(lmaxp+1))
      call czero(alpharkkrh,minprc*(lmaxp+1))
      call czero(alphaintkkrh,mintfc*(lmaxp+1))
      if(iscreen.ge.1) then
c Use square well t-matrix
c        --------------------------------------------------------------
         call alphamat(psq,iscreen,vscreen,lmax,E_Fermi,
     >                 nintfc,ninprcl,ninprcr,c)
c        --------------------------------------------------------------
      end if
c
c **********************************************************
c t-matrix for bulk and interface: screened and rotated to
c                                  global frame of reference
c **********************************************************
c
      do li=1,ninprcl
        cpalay=(1.0d0-dabs(concl(li))).gt.tiny
        if(.not.cpalay) then
c         -----------------------------------------------------------
          if (.not.localmode) then
             rbtmp = rbl(:,li)
          else 
             rbtmp = rz
          end if
          call tmat(
     >    ce,c,sxl(li),lmax,idpotla(li),v0,E_Fermi,deltala(1,li),
     >    vrla(1,li),brla(1,li),
     >    dxl(li),rsl(li),nsl(li),tx,tminvl(1,1,li),alphalkkr(0,li),
     >    alphalkkrh(0,li),rbtmp)
c         -----------------------------------------------------------
c          ! rotate tminvl to global frame if it's in local
c          if (localmode) then
c            call tripmt(dmatl(1,1,li),tminvl(1,1,li),dmatlp(1,1,li),
c     >                  kmymax,kmymax,kmymaxp)
c          end if
c         -----------------------------------------------------------
        else
          if(.not.bulk.and..not.cpamatinl) then
            write(6,'(/'' Check cpamatinl !!!'')')
            stop
          end if
        endif
      enddo
c
      do li=1,ninprcr
        cpalay=(1.0d0-dabs(concr(li))).gt.tiny
        if(.not.cpalay) then
c         -----------------------------------------------------------
          if (.not.localmode) then
             rbtmp = rbr(:,li)
          else 
             rbtmp = rz
          end if
          call tmat(
     >    ce,c,sxr(li),lmax,idpotra(li),v0,E_Fermi,deltara(1,li),
     >    vrra(1,li),brra(1,li),
     >    dxr(li),rsr(li),nsr(li),tx,tminvr(1,1,li),alpharkkr(0,li),
     >    alpharkkrh(0,li),rbtmp)
c         -----------------------------------------------------------
c          ! rotate tminvr to global frame if it's in local
c          if (localmode) then
c            call tripmt(dmatr(1,1,li),tminvr(1,1,li),dmatrp(1,1,li),
c     >                  kmymax,kmymax,kmymaxp)
c          end if
c         -----------------------------------------------------------
        else
          if(.not.bulk.and..not.cpamatinr) then
            write(6,'(/'' Check cpamatinr !!!'')')
            stop
          end if
        endif
      enddo
c
      do li=1,nintfc
        cpalay=(1.d0-conc(li)).gt.tiny
c       ---------------------------------------------------------------
        if (.not.localmode) then
           rbtmp = rba(:,li)
        else 
           rbtmp = rz
        end if
        call tmat(
     >  ce,c,sxa(li),lmax,idpota(li),v0,E_Fermi,deltaa(1,li),
     >  vra(1,li),bra(1,li),
     >  dx(li),rs(li),ns(li),ptminva(1,1,li),tminva(1,1,li),
     >  alphaintkkr(0,li),alphaintkkrh(0,li),
     >  rbtmp)
        ! local version: tmat works in local frame;
        !                ptminv left there, tminv rotated to global
        ! global version: tmat works in global frame;
        !                 tminv left there, ptminv inverse rotated!
        ! TODO: make this consistent in global case, watch out for singlesite case!
c       ---------------------------------------------------------------
c        if (.not.localmode) then
          !write(6,*) '  screened tm1 matrix in global frame'
          !call outmat1(tminva(1,1,li),kmymax,kmymax,kmymaxp,tol,6)
          !write(6,*) '  unscreened tm1 matrix in global frame'
          !call outmat1(ptminva(1,1,li),kmymax,kmymax,kmymaxp,tol,6)
c          call tripmt(dmatpa(1,1,li),ptminva(1,1,li),dmata(1,1,li),
c     >                kmymax,kmymax,kmymaxp)
          !write(6,*) '  unscreened tm1 matrix rotated to local frame'
          !call outmat1(ptminva(1,1,li),kmymax,kmymax,kmymaxp,tol,6)
c        else
          !write(6,*) '  screened tm1 matrix in local frame'
          !call outmat1(tminva(1,1,li),kmymax,kmymax,kmymaxp,tol,6)
c          call tripmt(dmata(1,1,li),tminva(1,1,li),dmatpa(1,1,li),
c     >                kmymax,kmymax,kmymaxp)
          !write(6,*) '  screened tm1 matrix rotated to global frame'
          !call outmat1(tminva(1,1,li),kmymax,kmymax,kmymaxp,tol,6)
c        end if
c       ---------------------------------------------------------------
c=============
c== CPA if ===
c=============
        if (cpalay) then
          if (.not.localmode) then
             rbtmp = rbb(:,li)
          else 
             rbtmp = rz
          end if
          call tmat(
     >    ce,c,sxb(li),lmax,idpotb(li),v0,E_Fermi,deltab(1,li),
     >    vrb(1,li),brb(1,li),
     >    dx(li),rs(li),ns(li),ptminvb(1,1,li),tminvb(1,1,li),
     >    alphaintkkr(0,li),alphaintkkrh(0,li),
     >    rbtmp)
          else
             tminvb(:,:,li) = tminva(:,:,li)
        end if
c=====================
c=== End of CPA if ===
c=====================
c       ---------------------------------------------------------------
c        if (.not.localmode) then
c          call tripmt(dmatpb(1,1,li),ptminvb(1,1,li),dmatb(1,1,li),
c     >                kmymax,kmymax,kmymaxp)
c        else
c          call tripmt(dmatb(1,1,li),tminvb(1,1,li),dmatpb(1,1,li),
c     >                kmymax,kmymax,kmymaxp)
c        end if
c       ---------------------------------------------------------------
c
c        tminv(:,:,li)=tminva(:,:,li)
c
        if (.not.(cpamatin.and.cpalay)) then
        write(6,*) 'ATA approximation used for t-matrix in tmatini'
c       -------------------------------------------------------
        call tata(conc(li),2*kmymax,tminv(1,1,li),tminva(1,1,li),
     >            tminvb(1,1,li))
c       -------------------------------------------------------
        end if 
      end do
c
      return
      end
