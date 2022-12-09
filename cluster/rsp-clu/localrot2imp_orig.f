c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine localrot2imp(
c     ===================
     > lmax,nintfc,rba,vecna,phia,dmata,dmatpa,
     > ddpha,ddphpa,ddtha,ddthpa,rmata,rmatpa)
c    > rbl,rbr, rbb, vecnb, phib,
c    > dmatl,dmatlp,dmatr,dmatrp,
c    > dmatb, dmatpb, ddphb, ddphpb, ddthb, ddthpb,rmatb,rmatpb)
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      integer mmax
      parameter (mmax=mimp)
c
      character*30 for006
c
      complex*16 dz(3)
      real*8     tvec(3)
      real*8     x,y,rr
c
!     real*8 rbl(3,minprc)
!     real*8 rbr(3,minprc),
!     real*8 vecnb(3,mmax)
!     real*8 phib(mmax)
!     real*8 rbb(3,mmax)
      real*8 zdir(3)
      real*8 rba(3,mmax)
      real*8 vecna(3,mmax)
      real*8 phia(mmax)
c
c     complex*16 dmatl(kmymaxp,kmymaxp,minprc)
c     complex*16 dmatlp(kmymaxp,kmymaxp,minprc)
c     complex*16 dmatr(kmymaxp,kmymaxp,minprc)
c     complex*16 dmatrp(kmymaxp,kmymaxp,minprc)
c
c     complex*16 dmatb(kmymaxp,kmymaxp,mmax)
c     complex*16 dmatpb(kmymaxp,kmymaxp,mmax)
c     complex*16 rmatb(lmsup,lmsup,mmax)
c     complex*16 rmatpb(lmsup,lmsup,mmax)
c
c     complex*16 ddphb(kmymaxp,kmymaxp,mmax)
c     complex*16 ddphpb(kmymaxp,kmymaxp,mmax)
c     complex*16 ddthb(kmymaxp,kmymaxp,mmax)
c     complex*16 ddthpb(kmymaxp,kmymaxp,mmax)
c
      complex*16 dmata(kmymaxp,kmymaxp,mmax)
      complex*16 dmatpa(kmymaxp,kmymaxp,mmax)
      complex*16 rmata(lmsup,lmsup,mmax)
      complex*16 rmatpa(lmsup,lmsup,mmax)
c
      complex*16 ddph(kmymaxp,kmymaxp)
      complex*16 ddphp(kmymaxp,kmymaxp)
      complex*16 d2dph(kmymaxp,kmymaxp)
      complex*16 d2dphp(kmymaxp,kmymaxp)
      complex*16 ddth(kmymaxp,kmymaxp)
      complex*16 ddthp(kmymaxp,kmymaxp)
      complex*16 d2dth(kmymaxp,kmymaxp)
      complex*16 d2dthp(kmymaxp,kmymaxp)
      complex*16 d2dthph(kmymaxp,kmymaxp)
      complex*16 d2dthphp(kmymaxp,kmymaxp)
c
      complex*16 ddpha(kmymaxp,kmymaxp,mmax)
      complex*16 ddphpa(kmymaxp,kmymaxp,mmax)
      complex*16 ddtha(kmymaxp,kmymaxp,mmax)
      complex*16 ddthpa(kmymaxp,kmymaxp,mmax)
c    
      common/test/itest
c     common/lay2d/cvec(mtotal,3),nextra,nbulkl,nbulkr,
c    &             nprc,ninprc(0:mprc+1)
c
      data tol/1.0d-8/ 
      data tiny/1.0d-6/
c
c **************************************
c initialize rotation matrices:  z --> B
c **************************************
c
c     ninprcl = ninprc(0)
c     ninprcr = ninprc(nprc+1)
c
      zdir(1)=0.d0
      zdir(2)=0.d0
      zdir(3)=1.d0
c
c     do li=1,ninprcl
c         -------------------------------------------------------------
c         call matrot(zdir,rbl(1,li),lmax,dmatl(1,1,li),dmatlp(1,1,li),
c    >                rmata(1,1,li),rmatpa(1,1,li),vecna(1,1),phia(1))
c         -------------------------------------------------------------
c     enddo
c
c     do li=1,ninprcr
c         -------------------------------------------------------------
c         call matrot(zdir,rbr(1,li),lmax,dmatr(1,1,li),dmatrp(1,1,li),
c    >                rmata(1,1,li),rmatpa(1,1,li),vecna(1,1),phia(1))
c         -------------------------------------------------------------
c     enddo
c
      do li=1,nintfc
c       -------------------------------------------------------------
        call matrot(zdir,rba(1,li),lmax,dmata(1,1,li),dmatpa(1,1,li),
     >              rmata(1,1,li),rmatpa(1,1,li),vecna(1,li),phia(li))
        x=rba(1,li)
        y=rba(2,li)
        rr=x*x+y*y
        if(rr.gt.tiny) 
     >  call matrot2(zdir,rba(1,li),lmax,dmata(1,1,li),dmatpa(1,1,li),
     >              dz,tvec,phia(li),
     >              ddpha(1,1,li),ddphpa(1,1,li),d2dph,d2dphp,
     >              ddtha(1,1,li),ddthpa(1,1,li),d2dth,d2dthp,
     >              d2dthph,d2dthphp)
c       -------------------------------------------------------------
c       call matrot(zdir,rbb(1,li),lmax,dmatb(1,1,li),dmatpb(1,1,li),
c    >              rmatb(1,1,li),rmatpb(1,1,li),vecnb(1,li),phib(li))
c       x=rbb(1,li)
c       y=rbb(2,li)
c       rr=x*x+y*y
c       if(rr.gt.tiny) 
c    >  call matrot2(zdir,rbb(1,li),lmax,dmatb(1,1,li),dmatpb(1,1,li),
c    >              dz,tvec,phib(li),
c    >              ddphb(1,1,li),ddphpb(1,1,li),d2dph,d2dphp,
c    >              ddthb(1,1,li),ddthpb(1,1,li),d2dth,d2dthp,
c    >              d2dthph,d2dthphp)
c       -------------------------------------------------------------
      enddo
c
      return
      end
