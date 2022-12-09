c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine tmatiniimp(
c     ==================
     > ce,lmax,nintfc,nimp,nposimp,iscreen,vscreen,v0,
     > wrel,sxcl,idpota,vra,bra,bopra,dx,ns,rs,
     > dmata,dmatpa,tminva,ptminva,ptminv0)
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      integer mmax
      parameter (mmax=mimp)
c
      logical wrel
c
      character*10 idpota(mmax)
      character*30 for006
c
      integer nposimp(3,mimp)
      integer ns(mmax)
      integer nintfc
      integer nimp
c
      integer imp
      integer li
c
      real*8 vra(nrad,mmax)
      real*8 bra(nrad,mmax)
      real*8 bopra(nrad,2,mmax)
c
      real*8 rs(mmax)
      real*8 dx(mmax)
      real*8 sxcl(mmax)
c
      complex*16 ce,we,psq
c
      complex*16 tx(kmymaxp,kmymaxp)
      complex*16 tminv(kmymaxp,kmymaxp,mmax)
      complex*16 tminva(kmymaxp,kmymaxp,mmax)
      complex*16 ptminva(kmymaxp,kmymaxp,mmax)
      complex*16 ptminv0(kmymaxp,kmymaxp,mmax)
c
      complex*16 dmata(kmymaxp,kmymaxp,mmax)
      complex*16 dmatpa(kmymaxp,kmymaxp,mmax)
c
c
      complex*16 alphaintkkr(0:lmaxp,mintfc)
      common/scrpar/alphaintkkr
c
      integer itest
      common/test/itest
c
      real*8 tol
      data tol/1.0d-8/ 
c
      real*8 tiny
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
c
c ********************
c screening parameters
c ********************
c
      call czero(alphaintkkr,mintfc*(lmaxp+1))
c
!     if(itest.ge.2) then
!       write(6,*) '<tmatiniimp>: czero ready'
!       write(6,*) '<tmatiniimp>: nintfc=',nintfc
!       write(6,*) '<tmatiniimp>: mintfc=',mintfc
!     end if
c
      if(iscreen.ge.1) then
c Use square well t-matrix
c        ------------------------------------------------------
         call alphamat1(psq,vscreen,lmax,nintfc)
c        ------------------------------------------------------
      end if
c
!     if(itest.ge.2) then
!       write(6,*) '<tmatiniimp>: alphamat1 ready'
!     end if
c
c **********************************************************
c t-matrix for bulk and interface: screened and rotated to
c                                  global frame of reference
c **********************************************************
c
      do imp=1,nimp
!       write(6,*) '<tmatiniimp>: imp=',imp
        li=nposimp(3,imp)
!     if(itest.ge.2) then
!       write(6,*) '<tmatiniimp>: li=',li
!       write(6,*) '<tmatiniimp>: ce=',ce
!       write(6,*) '<tmatiniimp>: psq=',psq
!       write(6,*) '<tmatiniimp>: idpota=',idpota(imp)
!       write(6,*) '<tmatiniimp>: v0=',v0
!       write(6,*) '<tmatiniimp>: dx=',dx(imp)
!       write(6,*) '<tmatiniimp>: ns=',ns(imp)
!       write(6,*) '<tmatiniimp>: rs=',rs(imp)
!       write(6,*) '<tmatiniimp>: sxcl=',sxcl(imp)
!     end if
c       ---------------------------------------------------------------
        call tmat(ce,psq,lmax,idpota(imp),v0,vra(1,imp),bra(1,imp),
     >       bopra(1,1,imp),dx(imp),ns(imp),rs(imp),
     >       ptminva(1,1,imp),tminva(1,1,imp),sxcl(imp),
     >       alphaintkkr(0,li))
c       ---------------------------------------------------------------
c
!     if(itest.ge.2) then
!       write(6,*) '<tmatiniimp>: tmat ready',imp
!     end if
c
c save ptminv, inverso of t-matrix in the local frame, for derivatives
        call repl(ptminv0(1,1,imp),ptminva(1,1,imp),kmymax,kmymaxp)
        call tripmt(dmata(1,1,imp),tminva(1,1,imp),dmatpa(1,1,imp),
     >              kmymax,kmymax,kmymaxp)
        call tripmt(dmata(1,1,imp),ptminva(1,1,imp),dmatpa(1,1,imp),
     >              kmymax,kmymax,kmymaxp)
c
!     if(itest.ge.2) then
!       write(6,*) '<tmatiniimp>: tripmt ready',imp
!     end if
!
c       ---------------------------------------------------------------
      end do
c
!     do li=1,nintfc
!       cpalay=(1.d0-conc(li)).gt.tiny
c       ---------------------------------------------------------------
!       call tmat(ce,psq,lmax,idpotb(li),v0,vrb(1,li),brb(1,li),
!    >       boprb(1,1,li),dx(li),ns(li),rs(li),
!    >       ptminvb(1,1,li),tminvb(1,1,li),sxb(li),alphaintkkr(0,li))
c       ---------------------------------------------------------------
!       call tripmt(dmatb(1,1,li),tminvb(1,1,li),dmatpb(1,1,li),
!    >              kmymax,kmymax,kmymaxp)
!       call tripmt(dmatb(1,1,li),ptminvb(1,1,li),dmatpb(1,1,li),
!    >              kmymax,kmymax,kmymaxp)
c       ---------------------------------------------------------------
!       if (.not.(cpamatin.and.cpalay)) 
c       -------------------------------------------------------
!    >  call tata(conc(li),kmymax,tminv(1,1,li),tminva(1,1,li),
!    >            tminvb(1,1,li))
c       -------------------------------------------------------
!     end do
c
      return
      end
