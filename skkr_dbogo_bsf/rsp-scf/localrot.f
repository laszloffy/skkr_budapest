c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine localrot(
c     ===================
     > lmax,nintfc,rbl,rbr,rba,rbb,vecna,vecnb,phia,phib,
     > dmatl,dmatlp,dmatr,dmatrp,dmata,dmatb,dmatpa,dmatpb,
     > rmata,rmatpa,rmatb,rmatpb)
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      character*30 for006
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
c rmata,rmatpa,vecna,phia dummy here !!!
c         -------------------------------------------------------------
          call matrot(zdir,rbl(1,li),lmax,dmatl(1,1,li),dmatlp(1,1,li),
     >                rmata(1,1,li),rmatpa(1,1,li),vecna(1,1),phia(1))
c         -------------------------------------------------------------
      enddo
c
      do li=1,ninprcr
c rmata,rmatpa,vecna,phia dummy here !!!
c         -------------------------------------------------------------
          call matrot(zdir,rbr(1,li),lmax,dmatr(1,1,li),dmatrp(1,1,li),
     >                rmata(1,1,li),rmatpa(1,1,li),vecna(1,1),phia(1))
c         -------------------------------------------------------------
      enddo
c
      do li=1,nintfc
c       -------------------------------------------------------------
        call matrot(zdir,rba(1,li),lmax,dmata(1,1,li),dmatpa(1,1,li),
     >              rmata(1,1,li),rmatpa(1,1,li),vecna(1,li),phia(li))
c       -------------------------------------------------------------
        call matrot(zdir,rbb(1,li),lmax,dmatb(1,1,li),dmatpb(1,1,li),
     >              rmatb(1,1,li),rmatpb(1,1,li),vecnb(1,li),phib(li))
c       -------------------------------------------------------------
      enddo
c
      return
      end
