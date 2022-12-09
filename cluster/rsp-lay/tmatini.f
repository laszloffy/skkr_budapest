c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine tmatini(
c     ==================
     > ce,lmax,nintfc,iscreen,vscreen,v0,
     > bulk,wrel,sxl,sxr,sxa,sxb,
     > concl,concr,conc,cpamatin,cpamatinl,cpamatinr,
c
     > idpotla,vrla,brla,boprla,
     > idpotlb,vrlb,brlb,boprlb,rsl,dxl,nsl,
     > idpotra,vrra,brra,boprra,
     > idpotrb,vrrb,brrb,boprrb,rsr,dxr,nsr,
     > idpota,vra,bra,bopra,idpotb,vrb,brb,boprb,dx,ns,rs,
c
     > dmatl,dmatlp,dmatr,dmatrp,dmata,dmatb,dmatpa,dmatpb,
c
     > tminvl,tminvr,tminv,tminva,tminvb,ptminva,ptminvb)
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      logical wrel,bulk,cpalay
      logical unsctest
      logical cpamatin,cpamatinl,cpamatinr
c
      character*10 idpota(mintfc),idpotla(minprc),idpotra(minprc)
      character*10 idpotb(mintfc),idpotlb(minprc),idpotrb(minprc)
      character*30 for006
c
      dimension vra(nrad,mintfc),bra(nrad,mintfc)
      dimension bopra(nrad,2,mintfc)
      dimension vrb(nrad,mintfc),brb(nrad,mintfc)
      dimension boprb(nrad,2,mintfc)
      dimension rs(mintfc),dx(mintfc),ns(mintfc)
      dimension vrla(nrad,minprc),brla(nrad,minprc)
      dimension boprla(nrad,2,minprc)
      dimension vrlb(nrad,minprc),brlb(nrad,minprc)
      dimension boprlb(nrad,2,minprc)
      dimension rsl(minprc),dxl(minprc),nsl(minprc)
      dimension vrra(nrad,minprc),brra(nrad,minprc)
      dimension boprra(nrad,2,minprc)
      dimension vrrb(nrad,minprc),brrb(nrad,minprc)
      dimension boprrb(nrad,2,minprc)
      dimension rsr(minprc),dxr(minprc),nsr(minprc)
c
      dimension conc(mintfc),sxa(mintfc),sxb(mintfc)
      dimension concl(minprc),sxl(minprc)
      dimension concr(minprc),sxr(minprc)
c
      complex*16 ce,we,psq
c
      complex*16 tx(kmymaxp,kmymaxp)
      complex*16 tminv(kmymaxp,kmymaxp,mintfc)
      complex*16 tminva(kmymaxp,kmymaxp,mintfc)
      complex*16 tminvb(kmymaxp,kmymaxp,mintfc)
      complex*16 ptminva(kmymaxp,kmymaxp,mintfc)
      complex*16 ptminvb(kmymaxp,kmymaxp,mintfc)
      complex*16 tminvl(kmymaxp,kmymaxp,minprc)
      complex*16 tminvr(kmymaxp,kmymaxp,minprc)
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
      common/scrpar/alphalkkr,alpharkkr,alphaintkkr
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
c
      ninprcl = ninprc(0)
      ninprcr = ninprc(nprc+1)
c
c ********************
c screening parameters
c ********************
c
c zero.f is real*8 => double arrays size!
      call czero(alphalkkr,minprc*(lmaxp+1))
      call czero(alpharkkr,minprc*(lmaxp+1))
      call czero(alphaintkkr,mintfc*(lmaxp+1))
      if(iscreen.ge.1) then
c Use square well t-matrix
c        ------------------------------------------------------
         call alphamat(psq,vscreen,lmax,nintfc,ninprcl,ninprcr)
c        ------------------------------------------------------
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
          call tmat(ce,psq,lmax,idpotla(li),v0,vrla(1,li),brla(1,li),
     >              boprla(1,1,li),dxl(li),nsl(li),rsl(li),
     >              tx,tminvl(1,1,li),sxl(li),alphalkkr(0,li))
c         -----------------------------------------------------------
          call tripmt(dmatl(1,1,li),tminvl(1,1,li),dmatlp(1,1,li),
     >                kmymax,kmymax,kmymaxp)
c         -----------------------------------------------------------
        else
          if(.not.bulk.and..not.cpamatinl) then
            write(6,'(/'' <tmatini>: STOP: Check cpamatinl !!!'')')
            stop
          end if
        endif
      enddo
c
      do li=1,ninprcr
        cpalay=(1.0d0-dabs(concr(li))).gt.tiny
        if(.not.cpalay) then
c         -----------------------------------------------------------
          call tmat(ce,psq,lmax,idpotra(li),v0,vrra(1,li),brra(1,li),
     >              boprra(1,1,li),dxr(li),nsr(li),rsr(li),
     >              tx,tminvr(1,1,li),sxr(li),alpharkkr(0,li))
c         -----------------------------------------------------------
          call tripmt(dmatr(1,1,li),tminvr(1,1,li),dmatrp(1,1,li),
     >                kmymax,kmymax,kmymaxp)
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
        call tmat(ce,psq,lmax,idpota(li),v0,vra(1,li),bra(1,li),
     >       bopra(1,1,li),dx(li),ns(li),rs(li),
     >       ptminva(1,1,li),tminva(1,1,li),sxa(li),alphaintkkr(0,li))
c       ---------------------------------------------------------------
        call tripmt(dmata(1,1,li),tminva(1,1,li),dmatpa(1,1,li),
     >              kmymax,kmymax,kmymaxp)
        call tripmt(dmata(1,1,li),ptminva(1,1,li),dmatpa(1,1,li),
     >              kmymax,kmymax,kmymaxp)
c DEBUG
      unsctest=.false.
      if(unsctest) then
         write(6,*) '<tmatini>: unscreened tm1 matrix',li
         call outmat1(ptminva(1,1,li),kmymax,kmymax,kmymaxp,tol,6)
      end if
c DEBUG
c       ---------------------------------------------------------------
        call tmat(ce,psq,lmax,idpotb(li),v0,vrb(1,li),brb(1,li),
     >       boprb(1,1,li),dx(li),ns(li),rs(li),
     >       ptminvb(1,1,li),tminvb(1,1,li),sxb(li),alphaintkkr(0,li))
c       ---------------------------------------------------------------
        call tripmt(dmatb(1,1,li),tminvb(1,1,li),dmatpb(1,1,li),
     >              kmymax,kmymax,kmymaxp)
        call tripmt(dmatb(1,1,li),ptminvb(1,1,li),dmatpb(1,1,li),
     >              kmymax,kmymax,kmymaxp)
c       ---------------------------------------------------------------
        if (.not.(cpamatin.and.cpalay)) 
c       -------------------------------------------------------
     >  call tata(conc(li),kmymax,tminv(1,1,li),tminva(1,1,li),
     >            tminvb(1,1,li))
c       -------------------------------------------------------
      end do
c
      return
      end
