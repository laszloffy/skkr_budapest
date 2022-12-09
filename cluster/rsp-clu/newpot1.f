c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine newpot1(itscf,orbpol,nimp,nposimp,
     >            za,qca,qva,
     >            rhoca,rhova,rhospa,dx,ns,rs,
     >            vra,bra,bopra,
     >            enpota,enela,enxca,enmaga,
     >            lquad,vmadih,vmadich,vmadid,vmadiq)
c
c ======================================================================
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      logical orbpol
c
      real*8 dx(mimp)
      real*8 rs(mimp)
      real*8 za(mimp)
c
      real*8 vra(nrad,mimp)
      real*8 bra(nrad,mimp)
      real*8 bopra(nrad,2,mimp)
c
      real*8 rhoca(nrad,mimp)
      real*8 rhova(nrad,mimp)
      real*8 rhospa(nrad,2,mimp)
      real*8 rhodspa(nrad,2,mimp)
c
      real*8 lza(2,mimp)
      real*8 sxa(mimp)
      real*8 qca(mimp)
      real*8 qva(mimp)
      real*8 qimp(mimp)
      real*8 dimp(mimp)
c
      real*8 enpota(mimp)
      real*8 enela(mimp)
      real*8 enxca(mimp)
      real*8 enmaga(mimp)
c
      real*8 rgrid(nrad)
      real*8 rhota(nrad)
      real*8 pota(nrad)
      real*8 xpot1(nrad)
      real*8 xpot2(nrad)
      real*8 rxa(nrad)
      real*8 rxb(nrad)
      real*8 rtmp(nrad)
c
      real*8 rdiff(3)
      real*8 rdiff2
      real*8 qa 
c
      real*8 vmadih(mimp)
      real*8 vmadich(mimp)
      real*8 vmadid(mimp)
      real*8 vmadiq(mimp)
c
      integer ns(mimp)
      integer nimp
      integer nposimp(3,mimp)
      integer lquad
c
      integer liimp
      integer ljimp
      integer iimp
      integer jimp
      integer n0
      integer i
c
c --- Common blocks ---
      real*8 a1
      real*8 a2
      common/brav2d/a1(2),a2(2)
c
      real*8 cvec
      integer nextra
      integer nbulkl
      integer nbulkr
      integer nprc
      integer ninprc
      common/lay2d/cvec(mtotal,3),nextra,nbulkl,nbulkr,
     &             nprc,ninprc(0:mprc+1)
c
      common/test/itest
!     common/latt2d/rlay0(3),i2dlat,ninprc,dperp,a2d,b2d
      common/xch/ixch
c ======================================================================
c
      write(6,'(/2x,''<newpot1>: routine NEWPOT>''/)')
c
      pi=4.d0*datan(1.d0)
      fpi=4.d0*pi
c
c charges for Madelung potential
c
      qif=0.d0
      do iimp=1,nimp
        qa=qca(iimp)+qva(iimp)-za(iimp)
        qimp(iimp)=qa
        if(itest.ge.2)
     >  write(6,'('' I'',i3,''  Layer'',i2,''  q ='',f15.10)')
     >            itscf,iimp,qimp(iimp)
        qif=qif+qimp(iimp)
      end do
      write(6,'('' I'',i3,''  Cluster  q='',f15.10)') itscf,qif
c
c =====================================================================
c          ***************************************************
c          * Madelung potential for impurities indexed by 00 *
c          ***************************************************
c
      n0 = ninprc(0)*(nextra+1)
c
      do iimp=1,nimp
        liimp=nposimp(3,iimp)
        vmadich(iimp)=0.d0
          do jimp=1,nimp
          ljimp=nposimp(3,jimp)
            if(iimp.ne.jimp) then 
              rdiff(1)=(nposimp(1,iimp)-nposimp(1,jimp))*a1(1)+
     >                 (nposimp(2,iimp)-nposimp(2,jimp))*a2(1)+
     >                 (cvec(liimp+n0,1)-cvec(ljimp+n0,1))
c
              rdiff(2)=(nposimp(1,iimp)-nposimp(1,jimp))*a1(2)+
     >                 (nposimp(2,iimp)-nposimp(2,jimp))*a2(2)+
     >                 (cvec(liimp+n0,2)-cvec(ljimp+n0,2))
c
              rdiff(3)=cvec(liimp+n0,3)-cvec(ljimp+n0,3)
c
              rdiff2=dsqrt(rdiff(1)*rdiff(1)+rdiff(2)*rdiff(2)+
     >               rdiff(3)*rdiff(3))
c
!         vmadi(iimp)=vmadi(iimp)+(2.d0*qimp(jimp))/rdiff2
          vmadich(iimp)=vmadich(iimp)+(2.d0*qimp(jimp))/rdiff2
        end if
        end do
      end do
c
c =====================================================================
c
c solve Poisson equation and
c put together muffin tin potential
c
      if(itest.ge.1) write(6,*)
      do iimp=1,nimp
        den1a=0.d0
        den2a=0.d0
        r1=rs(iimp)
        n1=ns(iimp)
        dx1=dx(iimp)
        xa=dlog(r1)
        x0=dlog(r1)-(n1-1)*dx1
c
c smooth charge density near to origo (for vacuum layers useful)
        x=xa+dx1
        do i=n1,1,-1
          x=x-dx1
          r=dexp(x)
          fpirr=fpi*r*r
c -A
          den=rhospa(i,1,iimp)/fpirr
          if(den.gt.0.d0) then
            den1a=den
          else
            den=den1a
          end if
          rhospa(i,1,iimp)=fpirr*den
          den=rhospa(i,2,iimp)/fpirr
          if(den.gt.0.d0) then
            den2a=den
          else
            den=den2a
          end if
          rhospa(i,2,iimp)=fpirr*den
          rhova(i,iimp)=rhospa(i,1,iimp)+rhospa(i,2,iimp)
        end do
c
        x=x0
        do i=1,n1
          rgrid(i)=dexp(x)
          rhota(i)=rhoca(i,iimp)+rhova(i,iimp)
          x=x+dx1
        end do
c       --------------------
        call rzero(pota,nrad)
c       --------------------
        qa=qca(iimp)+qva(iimp)
c       -----------------------------------------
        call poissond(rhota,qa,pota,dx1,x0,r1,n1)
c       -----------------------------------------
c
c constant to be added to the electrostatic potential
c
      if(lquad.eq.0) then
        vcoul=vmadih(iimp)+vmadich(iimp)+vmadid(iimp)
      else
        vcoul=vmadih(iimp)+vmadich(iimp)+vmadid(iimp)+vmadiq(iimp)
      end if
c
        if(itest.ge.1) then
        write(6,'('' I'', i3,''  Imp.'',i3,''  Mad. pot. '',
     >            ''  Vm(h)='',f15.10)') itscf,iimp,vmadih(iimp)
        write(6,'('' I'', i3,''  Imp.'',i3,''  Mad. pot. '',
     >            ''  Vm(0)='',f15.10)') itscf,iimp,vmadich(iimp)
        write(6,'('' I'', i3,''  Imp.'',i3,''  Mad. pot. '',
     >            ''  Vm(1)='',f15.10)') itscf,iimp,vmadid(iimp)
        write(6,'('' I'', i3,''  Imp.'',i3,''  Mad. pot. '',
     >            ''  Vm(h+0+1)='',f15.10)') itscf,iimp,
     >            vmadih(iimp)+vmadich(iimp)+vmadid(iimp)
        write(6,'('' I'', i3,''  Imp.'',i3,''  Mad. pot. '',
     >            ''  Vm(2)='',f15.10)') itscf,iimp,vmadiq(iimp)
        write(6,'('' I'', i3,''  Imp.'',i3,''  Mad. pot. '',
     >            ''  Vm(t)='',f15.10,
     >            ''  Vc='',f15.10)') itscf,iimp,vcoul,vcoul
        write(6,*)
        end if
c
c electrostatic potential for impurity iimp
c
        do i=1,n1
          pota(i)=pota(i)+vcoul
          rxa(i)=rhota(i)*pota(i)
          vnucla=-2.d0*za(iimp)/rgrid(i)
          pota(i)=pota(i)+vnucla
        end do
c
c electrostatic energy for layer iimp
        enela(iimp)=rsimp(rxa,rgrid,n1,dx1)+za(iimp)*vcoul
        enela(iimp)=-0.5d0*enela(iimp)
c
        if(itest.ge.3) then
          write(6,'(/'' Layer'',i2)') iimp
          write(6,'('' Electrostatic potential'')')
          write(6,'(4d20.12)') (pota(i),i=1,n1)
          write(6,*)
        end if
c
c exchange-correlation potential and energy for layer iimp
        do i=1,n1
          r=rgrid(i)
          fpir2=fpi*r*r
c
          rho1=(0.5d0*rhoca(i,iimp)+rhospa(i,1,iimp))/fpir2
          rho2=(0.5d0*rhoca(i,iimp)+rhospa(i,2,iimp))/fpir2
          rho=rho1+rho2
c         -------------------------------------------------------
          call exchange(rho1,rho2,rho,ixch,xpot1(i),xpot2(i),exc)
c         -------------------------------------------------------
          vxc=0.5d0*(xpot1(i)+xpot2(i))
          bxc=0.5d0*(xpot2(i)-xpot1(i))
          if(dabs(bxc).lt.1.d-15) bxc=0.d0
          pota(i)=pota(i)+vxc
          vra(i,iimp)=r*pota(i)
          bra(i,iimp)=r*bxc
          rxa(i)=rhota(i)*exc
          rtmp(i)=rhota(i)*vxc
          rxb(i)=(rhospa(i,2,iimp)-rhospa(i,1,iimp))*bxc
        end do
        enxca(iimp)=rsimp(rxa,rgrid,n1,dx1)
        enpota(iimp)=-rsimp(rtmp,rgrid,n1,dx1)
        enmaga(iimp)=-rsimp(rxb,rgrid,n1,dx1)
c
c  rxb has been used for temporary storage!
c
        if(itest.ge.3) then
          write(6,'('' A-type spin-up x-c potential'')')
          write(6,'(4d20.12)') (xpot1(i),i=1,n1)
          write(6,'('' A-type spin-down x-c potential'')')
          write(6,'(4d20.12)') (xpot2(i),i=1,n1)
        end if
c
        if(itest.ge.3) then
          write(6,'('' Total potential'')')
          write(6,'(4d20.12)') (pota(i),i=1,n1)
          write(6,'('' magnetic field'')')
          write(6,'(4d20.12)') (bra(i,iimp)/rgrid(i),i=1,n1)
        end if
c
      if(orbpol) then
c
        call bop(rs(iimp),ns(iimp),dx(iimp),rhodspa(1,1,iimp),
     >           lza(1,iimp),bopra(1,1,iimp))
c
        if(itest.ge.3) then
          write(6,'(/ ''LAYER '',i3)') iimp
          write(6,'(/ ''BOPR A SPIN 1'')')
          write(6,'(4d20.10)') (bopra(j,1,iimp),j=1,ns(iimp))
          write(6,'(/ ''BOPR A SPIN 2'')')
          write(6,'(4d20.10)') (bopra(j,2,iimp),j=1,ns(iimp))
        end if
c
      end if
c
      end do
c
      return
      end
