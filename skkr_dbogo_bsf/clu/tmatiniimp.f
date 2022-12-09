c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine tmatiniimp(
c     ==================
     > ce,lmax,tdim,nintfc,nimp,nposimp,iscreen,vscreen,v0,
     > wrel,sxcl,idpotimp,vrimp,brimp,boprimp,dx,ns,rs,
c    > dmata,dmatpa,
     > tminvcl,ptminvcl,!ptminv0,
     > rbcl,deltaimp,E_Fermi,c_light,singratimp,uratimp,dratimp)
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      logical wrel
      logical tautest
c
      integer nimp
      character*10 idpotimp(nimp)
      character*30 for006
c
      integer nposimp(3,nimp)
      integer ns(nimp)
      integer nintfc
      integer tdim
c
      integer imp
      integer li
c
      real*8 vrimp(nrad,nimp)
      real*8 brimp(nrad,nimp)
      real*8 boprimp(nrad,2,nimp)
      complex*16 deltaimp(nrad,nimp)
      dimension singratimp(nimp),uratimp(nimp),dratimp(nimp)
c
      real*8 rs(nimp)
      real*8 dx(nimp)
      real*8 sxcl(nimp)
      real*8 rbcl(3,nimp)
c
      complex*16 ce,we,psq
c
c      complex*16 tx(kmymaxp,kmymaxp)
c     complex*16 tminv(dbogomaxp,dbogomaxp,nimp)
      complex*16 tminvhelp(dbogomaxp,dbogomaxp)
      complex*16 ptminvhelp(dbogomaxp,dbogomaxp)
      complex*16 tminvcl(tdim,tdim,nimp)
      complex*16 ptminvcl(tdim,tdim,nimp)
c     complex*16 ptminv0(tdim,tdim,nimp)
c
c     complex*16 dmata(kmymaxp,kmymaxp,nimp)
c     complex*16 dmatpa(kmymaxp,kmymaxp,nimp)
c
c
      complex*16 alphaintkkr(0:lmaxp,mintfc)
      complex*16 alphaintkkrh(0:lmaxp,mintfc)
      common/scrpar/alphaintkkr
      common/scrparh/aphaintkkrh
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
c c initialize constants
c initialize constants
c ********************
c
c---> c in rydberg units:
c      c=274.072d0
      c=c_light
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
      call czero(alphaintkkrh,mintfc*(lmaxp+1))
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
         call alphamat1(psq,vscreen,lmax,nintfc) ! TODO check this
c         call alphamat(psq,iscreen,vscreen,lmax,E_Fermi,
c     >                 nintfc,ninprcl,ninprcr,c)
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
        call tmat(
     >  ce,c_light,sxcl(imp),lmax,idpotimp(imp),v0,E_Fermi,
     >  deltaimp(1,imp),
     >  singratimp(imp),uratimp(imp),dratimp(imp),
     >  vrimp(1,imp),brimp(1,imp),
     >  dx(imp),rs(imp),ns(imp),ptminvhelp,tminvhelp,
     >  alphaintkkr(0,li),alphaintkkrh(0,li),rbcl(1,imp))
c       ---------------------------------------------------------------
        call repldim(tminvcl(1,1,imp),tminvhelp,tdim,tdim,dbogomaxp)
        call repldim(ptminvcl(1,1,imp),ptminvhelp,tdim,tdim,dbogomaxp)
c
c
      tautest=.false.
      if(tautest) then
       write(6,*) ' <tmatiniimp> : tautest, after tmat for a component'
       write(6,*) "ce=",ce
       write(6,*) "psq=",psq
       write(6,*) "sxcl=",sxcl(imp)
       write(6,*) "idpotimp(imp)=",idpotimp(imp)
       write(6,*) "v0=",v0
       write(6,*) "E_Fermi=",E_Fermi
       write(6,*) "deltaimp(1,imp)=",deltaimp(1:nrad,imp)
       write(6,*) "vrimp(1,imp)=",vrimp(1:nrad,imp)
       write(6,*) "brimp(1,imp)=",brimp(1:nrad,imp)
       write(6,*) "dx(imp)=",dx(imp)
       write(6,*) "rs(imp)=",rs(imp)
       write(6,*) "ns(imp)=",ns(imp)
       write(6,*) ' <tmatiniimp> : tautest, tminvhelp imp=',imp
       call outmat1(tminvhelp(1,1),tdim,tdim,dbogomaxp,
     >              tol,6)
       write(6,*) ' <tmatiniimp> : tautest, ptminvhelp imp=',imp
       call outmat1(ptminvhelp(1,1),tdim,tdim,dbogomaxp,
     >              tol,6)
       write(6,*) ' <tmatiniimp> : tautest, tminvcl imp=',imp
       call outmat1(tminvcl(1,1,imp),tdim,tdim,tdim,
     >              tol,6)
       write(6,*) ' <tmatiniimp> : tautest, ptminvcl imp=',imp
       call outmat1(ptminvcl(1,1,imp),tdim,tdim,tdim,
     >              tol,6)
       write(6,*) "alphaintkkr(0,li)=",alphaintkkr(0,li)
       write(6,*) "alphaintkkrh(0,li)=",alphaintkkrh(0,li)
       write(6,*) "rbcl(1,imp)=",rbcl(1,imp)
      end if
c
!     if(itest.ge.2) then
!       write(6,*) '<tmatiniimp>: tmat ready',imp
!     end if
c
c save ptminv, inverso of t-matrix in the local frame, for c        call repl(ptminv0(1,1,imp),ptminva(1,1,imp),kmymax,kmymaxp)
c        call repl(ptminv0(1,1,imp),ptminva(1,1,imp),kmymax,kmymaxp)
c        call tripmt(dmata(1,1,imp),tminva(1,1,imp),dmatpa(1,1,imp),
c     >              kmymax,kmymax,kmymaxp)
c        call tripmt(dmata(1,1,imp),ptminva(1,1,imp),dmatpa(1,1,imp),
c     >              kmymax,kmymax,kmymaxp)
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
