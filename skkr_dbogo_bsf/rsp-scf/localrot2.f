c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine localrot2(
c     ===================
     > lmax,nintfc,rbl,rbr,rba,rbb,vecna,vecnb,phia,phib,
     > dmatl,dmatlp,dmatr,dmatrp,dmata,dmatb,dmatpa,dmatpb,
     > ddpha,ddtha,d2dpha,d2dtha,d2dthpha,
     > ddphpa,ddthpa,d2dphpa,d2dthpa,d2dthphpa,
     > ddphb,ddthb,d2dphb,d2dthb,d2dthphb,
     > ddphpb,ddthpb,d2dphpb,d2dthpb,d2dthphpb,
     > rmata,rmatpa,rmatb,rmatpb)
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      dimension rbl(3,minprc),rbr(3,minprc),zdir(3)
      dimension rba(3,mintfc),rbb(3,mintfc)
      dimension vecna(3,mintfc),phia(mintfc)
      dimension vecnb(3,mintfc),phib(mintfc)
c
      complex*16 dmatl(kmymaxp,kmymaxp,minprc)
      complex*16 dmatlp(kmymaxp,kmymaxp,minprc)
      complex*16 dmatr(kmymaxp,kmymaxp,minprc)
      complex*16 dmatrp(kmymaxp,kmymaxp,minprc)
      complex*16 dmata(kmymaxp,kmymaxp,mintfc)
      complex*16 dmatpa(kmymaxp,kmymaxp,mintfc)
      complex*16 dmatb(kmymaxp,kmymaxp,mintfc)
      complex*16 dmatpb(kmymaxp,kmymaxp,mintfc)
      complex*16 rmata(lmsup,lmsup,mintfc),rmatpa(lmsup,lmsup,mintfc)
      complex*16 rmatb(lmsup,lmsup,mintfc),rmatpb(lmsup,lmsup,mintfc)
c
      complex*16 ddpha(kmymaxp,kmymaxp,mintfc)
      complex*16 ddtha(kmymaxp,kmymaxp,mintfc)
      complex*16 d2dpha(kmymaxp,kmymaxp,mintfc)
      complex*16 d2dtha(kmymaxp,kmymaxp,mintfc)
      complex*16 d2dthpha(kmymaxp,kmymaxp,mintfc)
      complex*16 ddphpa(kmymaxp,kmymaxp,mintfc)
      complex*16 ddthpa(kmymaxp,kmymaxp,mintfc)
      complex*16 d2dphpa(kmymaxp,kmymaxp,mintfc)
      complex*16 d2dthpa(kmymaxp,kmymaxp,mintfc)
      complex*16 d2dthphpa(kmymaxp,kmymaxp,mintfc)
c
      complex*16 ddphb(kmymaxp,kmymaxp,mintfc)
      complex*16 ddthb(kmymaxp,kmymaxp,mintfc)
      complex*16 d2dphb(kmymaxp,kmymaxp,mintfc)
      complex*16 d2dthb(kmymaxp,kmymaxp,mintfc)
      complex*16 d2dthphb(kmymaxp,kmymaxp,mintfc)
      complex*16 ddphpb(kmymaxp,kmymaxp,mintfc)
      complex*16 ddthpb(kmymaxp,kmymaxp,mintfc)
      complex*16 d2dphpb(kmymaxp,kmymaxp,mintfc)
      complex*16 d2dthpb(kmymaxp,kmymaxp,mintfc)
      complex*16 d2dthphpb(kmymaxp,kmymaxp,mintfc)
c    
      common/test/itest
      common/lay2d/cvec(mtotal,3),nextra,nbulkl,nbulkr,
     &             nprc,ninprc(0:mprc+1)
c
      data tol/1.0d-8/ 
      data tiny/1.0d-6/
c
c **************************************
c initialize rotation matrices:  z --> B
c **************************************
c
      ninprcl = ninprc(0)
      ninprcr = ninprc(nprc+1)
c
      zdir(1)=0.d0
      zdir(2)=0.d0
      zdir(3)=1.d0
c
      do li=1,ninprcl
c         -------------------------------------------------------------
          call matrot(zdir,rbl(1,li),lmax,dmatl(1,1,li),dmatlp(1,1,li),
     >                rmata(1,1,li),rmatpa(1,1,li),vecna(1,1),phia(1))
c         -------------------------------------------------------------
      enddo
c
      do li=1,ninprcr
c         -------------------------------------------------------------
          call matrot(zdir,rbr(1,li),lmax,dmatr(1,1,li),dmatrp(1,1,li),
     >                rmata(1,1,li),rmatpa(1,1,li),vecna(1,1),phia(1))
c         -------------------------------------------------------------
      enddo
c
      do li=1,nintfc
c       -------------------------------------------------------------
c       write(6,*) rba(1,li),rba(2,li),rba(3,li)
        call matrot(zdir,rba(1,li),lmax,dmata(1,1,li),dmatpa(1,1,li),
     >              rmata(1,1,li),rmatpa(1,1,li),vecna(1,li),phia(li))
        x=rba(1,li)
        y=rba(2,li)
        if(x*x+y*y.gt.tiny)
     >  call matrot2(
     >       zdir,rba(1,li),lmax,dmata(1,1,li),dmatpa(1,1,li),phia(li),
     >       ddpha(1,1,li),ddphpa(1,1,li),
     >       d2dpha(1,1,li),d2dphpa(1,1,li),
     >       ddtha(1,1,li),ddthpa(1,1,li),
     >       d2dtha(1,1,li),d2dthpa(1,1,li),
     >       d2dthpha(1,1,li),d2dthphpa(1,1,li))
c      write(6,*) x*x+y*y
c      write(6,*) 'ddtha'
c      call outmat1(ddtha(1,1,li),kmymaxp,kmymaxp,kmymaxp,1.0d-10,6)
c      write(6,*) 'ddthpa'
c      call outmat1(ddthpa(1,1,li),kmymaxp,kmymaxp,kmymaxp,1.0d-10,6)
c       -------------------------------------------------------------
        call matrot(zdir,rbb(1,li),lmax,dmatb(1,1,li),dmatpb(1,1,li),
     >              rmatb(1,1,li),rmatpb(1,1,li),vecnb(1,li),phib(li))
        x=rbb(1,li)
        y=rbb(2,li)
        if(x*x+y*y.gt.tiny)
     >  call matrot2(
     >       zdir,rbb(1,li),lmax,dmatb(1,1,li),dmatpb(1,1,li),phib(li),
     >       ddphb(1,1,li),ddphpb(1,1,li),
     >       d2dphb(1,1,li),d2dphpb(1,1,li),
     >       ddthb(1,1,li),ddthpb(1,1,li),
     >       d2dthb(1,1,li),d2dthpb(1,1,li),
     >       d2dthphb(1,1,li),d2dthphpb(1,1,li))
c       -------------------------------------------------------------
      enddo
c
      return
      end
