c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine vgen1(itscf,itscfcont,nimp,nposimp,
     > orbpol,opot,madmax,
     > za,qca,qva,
     > rhoca,rhova,rhospa,rhodspa,dx,ns,rs,
     > vra,bra,bopra,lza,
     > ferr1,
     > enpota,enela,enxca,enmaga,
     > vmadih,vmadich,vmadid,vmadiq)
c
c =======================================================================
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      integer mmax
      parameter (mmax=mimp)
c
      parameter (mbig=4*nrad*mmax+2)
c
      logical orbpol
      logical opot
c
      real*8 lza(2,mmax)
      real*8 dx(mmax)
      real*8 rs(mmax)
      real*8 vra(nrad,mmax)
      real*8 bra(nrad,mmax)
      real*8 bopra(nrad,2,mmax)
      real*8 gvra(nrad,mmax)
      real*8 gbra(nrad,mmax)
      real*8 gbopra(nrad,2,mmax)
      real*8 rhoca(nrad,mmax)
      real*8 rhova(nrad,mmax)
      real*8 rhospa(nrad,2,mmax)
      real*8 rhodspa(nrad,2,mmax)
      real*8 za(mmax)
      real*8 qca(mmax)
      real*8 qva(mmax)
      real*8 enpota(mmax)
      real*8 enmaga(mmax)
      real*8 enela(mmax)
      real*8 enxca(mmax)
c
      real*8 xbr(mbig)
      real*8 fbr(mbig)
      real*8 dfbr(mbig)
      real*8 pfbr(mbig)
      real*8 pxbr(mbig)
      real*8 denbr(mbig)
      real*8 defbr(mbig)
c
      real*8 x
      real*8 r
c
      integer nposimp(3,mmax)
      integer ns(mmax)
      integer nimp
      integer lquad
      integer madmax
      integer itscf
      integer itscfcont
c
      integer nbig
      integer nbig1
      integer nbig2
      integer nbig3
      integer iimp
      integer ibig
      integer ii
      integer i
c
      character*2  tempr
      character*50 temp1
      character*50 temp2
      character*50 temp3
      integer root, myrank, nprocs
      common/mpi/root,myrank,nprocs
c
      real*8 vmadih(mmax)
      real*8 vmadich(mmax)
      real*8 vmadid(mmax)
      real*8 vmadiq(mmax)
c
      common/test/itest
      common/broypot/cmix,wbr,nbrmax(100)
c
      data anul/1.0d-15/,tiny/1.0d-6/
c
      save ibr,ibrmax,ibrcase
c =====================================================================
c
      if(madmax.ge.2) lquad=1
c     ----------------------------------------------------------------
      call newpot1(itscfcont,orbpol,nimp,nposimp,
     >            za,qca,qva,
     >            rhoca,rhova,rhospa,dx,ns,rs,
     >            gvra,gbra,gbopra,
     >            enpota,enela,enxca,enmaga,
     >            lquad,vmadih,vmadich,vmadid,vmadiq)
c     ----------------------------------------------------------------
c
      if(orbpol) then
      do iimp=1,nimp
        if(itscf.eq.1.and.(.not.opot)) then
          do ii=1,ns(iimp) 
             bopra(ii,1,iimp)=gbopra(ii,1,iimp)
             bopra(ii,2,iimp)=gbopra(ii,2,iimp)
          end do
        end if
      end do
      end if
c
c - Potential mixing
c
      ibig=0
      do iimp=1,nimp
        x=dlog(rs(iimp))-(ns(iimp)-1)*dx(iimp)
        do i=1,ns(iimp)
          r=dexp(x)
          ibig=ibig+1
          xbr(ibig)=vra(i,iimp)
          fbr(ibig)=gvra(i,iimp)-vra(i,iimp)
          x=x+dx(iimp)
        end do
      end do
c
      nbig1=ibig
c
      do iimp=1,nimp
        x=dlog(rs(iimp))-(ns(iimp)-1)*dx(iimp)
        do i=1,ns(iimp)
          r=dexp(x)
          ibig=ibig+1
          xbr(ibig)=bra(i,iimp)
          fbr(ibig)=gbra(i,iimp)-bra(i,iimp)
          x=x+dx(iimp)
        end do
      end do
      nbig2=ibig
c
      do iimp=1,nimp
        x=dlog(rs(iimp))-(ns(iimp)-1)*dx(iimp)
        do i=1,ns(iimp)
          r=dexp(x)
          ibig=ibig+1
          xbr(ibig)=bopra(i,1,iimp)
          fbr(ibig)=gbopra(i,1,iimp)-bopra(i,1,iimp)
          ibig=ibig+1
          xbr(ibig)=bopra(i,2,iimp)
          fbr(ibig)=gbopra(i,2,iimp)-bopra(i,2,iimp)
          x=x+dx(iimp)
        end do
      end do
      nbig3=ibig
c
      nbig=ibig
c
      call dott(xbr,xbr,err1,nbig1)
      call dott(fbr,fbr,t1,nbig1)
      ferr1=dsqrt(t1/err1)
      write(6,'(/'' Error for V(r) :'',3d17.6)')  t1,err1,ferr1
c
      call dott(xbr(nbig1+1),xbr(nbig1+1),err2,nbig2-nbig1)
      call dott(fbr(nbig1+1),fbr(nbig1+1),t1,nbig2-nbig1)
      if(err2.lt.anul) then
        ferr2=0.d0
      else
        ferr2=dsqrt(t1/err2)
      end if
      write(6,'('' Error for Bsp(r) :'',3d17.6)')  t1,err2,ferr2
c
      call dott(xbr(nbig2+1),xbr(nbig2+1),err2,nbig3-nbig2)
      call dott(fbr(nbig2+1),fbr(nbig2+1),t1,nbig3-nbig2)
      if(err2.lt.anul) then
        ferr2=0.d0
      else
        ferr2=dsqrt(t1/err2)
      end if
      write(6,'('' Error for Bop(r) :'',3d17.6)')  t1,err2,ferr2
c
c     ----------------------------------------------------------------
c
      if(myrank.le.9) then
         write(tempr,'("0",i1)') myrank
      else
         write(tempr,'(i2)') myrank
      end if

      temp1 = 'ftn21.'//tempr
      temp2 = 'ftn22.'//tempr
      temp3 = 'ftn23.'//tempr
c
      if(itscf.eq.1) then
        ibrcase=1
        ibrmax=nbrmax(1)
        ibr=1
      end if
c
      write(6,'(/'' itscf,ibrcase,ibrmax,ibr:'',4i4/)')
     >itscf,ibrcase,ibrmax,ibr
c
      call broyd(xbr,fbr,0.d0,ibr,nbig,dfbr,pfbr,pxbr,denbr,defbr,
     >           cmix,wbr,21,22,23,'ftn21','ftn22','ftn23')
c     ----------------------------------------------------------------
c
      if(ibr.ge.ibrmax) then
        ibr=1
        ibrcase=ibrcase+1
        ibrmax=nbrmax(ibrcase)
      else
        ibr=ibr+1
      end if
c
      ibig=0
      do iimp=1,nimp
        x=dlog(rs(iimp))-(ns(iimp)-1)*dx(iimp)
        do i=1,ns(iimp)
          r=dexp(x)
          ibig=ibig+1
          vra(i,iimp)=xbr(ibig)
        end do
      end do
c
      do iimp=1,nimp
        x=dlog(rs(iimp))-(ns(iimp)-1)*dx(iimp)
        do i=1,ns(iimp)
          r=dexp(x)
          ibig=ibig+1
          bra(i,iimp)=xbr(ibig)
          x=x+dx(iimp)
        end do
      end do
c
      do iimp=1,nimp
        x=dlog(rs(iimp))-(ns(iimp)-1)*dx(iimp)
        do i=1,ns(iimp)
          r=dexp(x)
          ibig=ibig+1
          bopra(i,1,iimp)=xbr(ibig)
          ibig=ibig+1
          bopra(i,2,iimp)=xbr(ibig)
          x=x+dx(iimp)
        end do
      end do
c
      return
      end
