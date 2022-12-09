c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine localrot2imp(
c     ===================
     > lmax,nintfc,rba,vecna,phia,dmata,dmatpa,
     > ddpha,ddphpa,ddtha,ddthpa,rmata,rmatpa,
     > theta0,phi0,
     > d2dph, d2dphp, d2dth, d2dthp, d2dthph, d2dthphp)
c    > rbl,rbr, rbb, vecnb, phib,
c    > dmatl,dmatlp,dmatr,dmatrp,
c    > dmatb, dmatpb, ddphb, ddphpb, ddthb, ddthpb,rmatb,rmatpb)
c
c-----------------------------------------------------------------------------
c Routine to calculate the rotation matrices and their derivatives with
c respect of the polar and asimuthal angles phi and theta.
c First the R0 rotation is generated given by theta0 and phi0, then the
c R rotation is formed where the z-dir will be the (theta0, phi0)
c direction. The whole rotation is then RxR0.
c-----------------------------------------------------------------------------
      implicit none
c
      include '../param.h'
c
      integer mmax
      parameter (mmax=mimp)
c
c     character*30 for006
c
      complex*16 dz(3)
      real*8     tvec(3)
      real*8     x,y,rr
c
      real*8 zdir(3)
      real*8 rba(3,mmax)
      real*8 vecna(3,mmax)
      real*8 phia(mmax)
      real*8 theta0,phi0
c
      complex*16 dmata(kmymaxp,kmymaxp,mmax)
      complex*16 dmatpa(kmymaxp,kmymaxp,mmax)
      complex*16 rmata(lmsup,lmsup,mmax)
      complex*16 rmatpa(lmsup,lmsup,mmax)
c
      complex*16 d2dph(kmymaxp,kmymaxp,mmax)
      complex*16 d2dphp(kmymaxp,kmymaxp,mmax)
      complex*16 d2dth(kmymaxp,kmymaxp,mmax)
      complex*16 d2dthp(kmymaxp,kmymaxp,mmax)
      complex*16 d2dthph(kmymaxp,kmymaxp,mmax)
      complex*16 d2dthphp(kmymaxp,kmymaxp,mmax)
c
      complex*16 ddpha(kmymaxp,kmymaxp,mmax)
      complex*16 ddphpa(kmymaxp,kmymaxp,mmax)
      complex*16 ddtha(kmymaxp,kmymaxp,mmax)
      complex*16 ddthpa(kmymaxp,kmymaxp,mmax)
c
      complex*16 dmat0(kmymaxp,kmymaxp)
      complex*16 dmat0p(kmymaxp,kmymaxp)
      complex*16 rmat0(lmsup,lmsup)
      complex*16 rmat0p(lmsup,lmsup)
c      
      complex*16 dummy1(kmymaxp,kmymaxp)
      complex*16 dummy2(kmymaxp,kmymaxp)
      complex*16 wrk(kmymaxp,kmymaxp),wrk1(kmymaxp,kmymaxp)
      complex*16 wrk2(lmsup,lmsup),wrk3(lmsup,lmsup)
c
      real*8     tol, tiny      
      real*8     cost,sint,cosp,sinp
      real*8     r0(3)
      integer    j,li,nintfc,lmax,itest,kmymax
      integer    lmaxs, lmmaxs
c    
      common/test/itest
c
      data tol/1.0d-8/ 
      data tiny/1.0d-6/
c
      kmymax = 2*(lmax+1)**2
      lmaxs=2*lmax
      lmmaxs=(lmaxs+1)*(lmaxs+1)
      write(6,*) 'lmmaxs:',lmmaxs

c
c **************************************
c initialize rotation matrices:  z --> B
c **************************************
c
      zdir(1)=0.d0
      zdir(2)=0.d0
      zdir(3)=1.d0
c forms unit vectors parallel to the axes of the new rference system
      cost = dcos(theta0)
      sint = dsin(theta0)
      cosp = dcos(phi0)
      sinp = dsin(phi0)
c
      r0(1) = sint*cosp
      r0(2) = sint*sinp
      r0(3) = cost
      write(6,'('' orientation of the new frame:'',2x,3f8.4)'),
     >                                            (r0(j),j=1,3)
c R0 rotation matrix
      call matrot(zdir,r0,lmax,dmat0,dmat0p,
     >            rmat0,rmat0p,vecna,phia)

c       -------------------------------------------------------------
c
      do li=1,nintfc
c       -------------------------------------------------------------
        call matrot(r0,rba(1,li),lmax,dmata(1,1,li),dmatpa(1,1,li),
     >              rmata(1,1,li),rmatpa(1,1,li),vecna(1,li),phia(li))
        call doubmt(dmata(1,1,li),dmat0,kmymax,kmymaxp)
        call doubmt(rmata(1,1,li),rmat0,lmmaxs,lmsup)
        call doubmtm(dmat0p,dmatpa(1,1,li),kmymax,kmymaxp)
        call doubmtm(rmat0p,rmatpa(1,1,li),lmmaxs,lmsup)
c --Test--
c       write(6,*) '-------- Test dmat --------'
c       call doubmt1(dmatpa(1,1,li),dmata(1,1,li),wrk,kmymax,kmymaxp)        
c       call outmatc(wrk,kmymax,kmymax,kmymaxp,6,wrk1)
c       write(6,*) '-------- Test rmat --------'
c       call doubmt1(rmatpa(1,1,li),rmata(1,1,li),wrk2,lmmaxs,lmsup)        
c       call outmatc(wrk2,lmmaxs,lmmaxs,lmsup,6,wrk3)
c -- Test vege --       
c
        x=rba(1,li)
        y=rba(2,li)
        rr=x*x+y*y
c       -------------------------------------------------------------
c       if(rr.gt.tiny) 
        call matrot2(r0,rba(1,li),lmax,dummy1,dummy2,
     >              dz,tvec,phia(li),
     >              ddpha(1,1,li),ddphpa(1,1,li),
     >              d2dph(1,1,li),d2dphp(1,1,li),
     >              ddtha(1,1,li),ddthpa(1,1,li),
     >              d2dth(1,1,li),d2dthp(1,1,li),
     >              d2dthph(1,1,li),d2dthphp(1,1,li))
        call doubmt(ddpha(1,1,li),dmat0,kmymax,kmymaxp)
        call doubmt(d2dph(1,1,li),dmat0,kmymax,kmymaxp)
        call doubmt(ddtha(1,1,li),dmat0,kmymax,kmymaxp)
        call doubmt(d2dth(1,1,li),dmat0,kmymax,kmymaxp)
        call doubmt(d2dthph(1,1,li),dmat0,kmymax,kmymaxp)
c
        call doubmtm(dmat0p,ddphpa(1,1,li),kmymax,kmymaxp)
        call doubmtm(dmat0p,d2dphp(1,1,li),kmymax,kmymaxp)
        call doubmtm(dmat0p,ddthpa(1,1,li),kmymax,kmymaxp)
        call doubmtm(dmat0p,d2dthp(1,1,li),kmymax,kmymaxp)
        call doubmtm(dmat0p,d2dthphp(1,1,li),kmymax,kmymaxp)
c       -------------------------------------------------------------
      enddo
c Test eleje
c       write(6,*) 'Test rotation matrices'
c       call testrot(dmata,dmatpa,ddpha,ddphpa,nintfc,kmymax,1)
c       call testrot(dmata,dmatpa,ddtha,ddthpa,nintfc,kmymax,0)
c Test vege
c
      return
      end
      subroutine testrot(dmata,dmatpa,dda,ddpa,nimp,kmymax,itest)
      include '../param.h'
      integer mmax
      parameter (mmax=mimp)
      complex*16 dda(kmymaxp,kmymaxp,mmax)
      complex*16 ddpa(kmymaxp,kmymaxp,mmax)
      complex*16 dmata(kmymaxp,kmymaxp,mmax)
      complex*16 dmatpa(kmymaxp,kmymaxp,mmax)
      complex*16 amat(kmymaxp,kmymaxp)
      complex*16 bmat(kmymaxp,kmymaxp)
c
      do ii = 1,nimp
c
c      write(6,*) ' Impurity:',ii
       if(itest.eq.1)then
c       write(6,*) '  dmat*dmatp'       
        call doubmt1(dmata(1,1,ii),dmatpa(1,1,ii),amat,kmymax,kmymaxp)
c       call outmatc(amat,kmymax,kmymax,kmymaxp,6,bmat)
        call outmat1(amat,kmymax,kmymax,kmymaxp,1.d-6,6)
       endif
c
c      write(6,*) '  test rot:',ii
       call doubmt1(dda(1,1,ii),dmatpa(1,1,ii),amat,kmymax,kmymaxp)
       call doubmt1(dmata(1,1,ii),ddpa(1,1,ii),bmat,kmymax,kmymaxp)
       call addmat(amat,bmat,kmymax,kmymaxp)
c      call outmatc(amat,kmymax,kmymax,kmymaxp,6,bmat)
       call outmat1(amat,kmymax,kmymax,kmymaxp,1.d-6,6)
      end do
      return
      end
