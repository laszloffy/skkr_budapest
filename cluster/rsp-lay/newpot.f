c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine newpot(
     >            itscf,nintfc,lmax,sigma,rightm,bulk,
     >            orbpol,conc,concl,concr,
     >            za,qca,qva,qmoma,qmomla,qmomra,lza,
     >            zb,qcb,qvb,qmomb,qmomlb,qmomrb,lzb,
     >            rhoca,rhova,rhospa,rhodspa,rhocb,rhovb,rhospb,rhodspb,
     >            dx,ns,rs,vvac,vra,vrb,bra,brb,bopra,boprb,
     >            enpota,enela,enxca,enpotb,enelb,enxcb,enmaga,enmagb,
c ---cluster---
     >            nimp,nposimp,vmadih)
c23456789012345678901234567890123456789012345678901234567890123456789012
c        1         2         3         4         5         6         7
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      logical bulk,orbpol
      character*1 rightm
c
      real*8 lza(2,mintfc),lzb(2,mintfc)
      dimension vra(nrad,mintfc),vrb(nrad,mintfc)
      dimension bra(nrad,mintfc),brb(nrad,mintfc)
      dimension bopra(nrad,2,mintfc),boprb(nrad,2,mintfc)
      dimension dx(mintfc),ns(mintfc),rs(mintfc)
      dimension rhoca(nrad,mintfc),rhocb(nrad,mintfc)
      dimension rhova(nrad,mintfc),rhovb(nrad,mintfc)
      dimension rhospa(nrad,2,mintfc),rhospb(nrad,2,mintfc)
      dimension rhodspa(nrad,2,mintfc),rhodspb(nrad,2,mintfc)
      dimension za(mintfc),qca(mintfc),qva(mintfc)
      dimension zb(mintfc),qcb(mintfc),qvb(mintfc)
      dimension conc(mintfc),concr(minprc),concl(minprc)
      dimension vmad(mintfc)
      dimension enpota(mintfc),enela(mintfc),enxca(mintfc)
      dimension enpotb(mintfc),enelb(mintfc),enxcb(mintfc)
      dimension enmaga(mintfc),enmagb(mintfc)
c
      dimension qlay(mintfc)
      dimension rgrid(nrad),rhota(nrad),rhotb(nrad),pota(nrad),
     >potb(nrad),xpot1(nrad),xpot2(nrad)
      dimension rxa(nrad),rxb(nrad),rtmp(nrad)
c
      complex*16 qmoma(lmsup,mintfc)
      complex*16 qmomb(lmsup,mintfc)
      complex*16 qmomla(lmsup,minprc)
      complex*16 qmomlb(lmsup,minprc)
      complex*16 qmomra(lmsup,minprc)
      complex*16 qmomrb(lmsup,minprc)
      complex*16 qmoml(lmsup,minprc)
      complex*16 qmomr(lmsup,minprc)
      complex*16 qmomlay(lmsup,mintfc)
c -------------------------------------------
c variables for cluster calculation
      logical fullmad
      integer nimp
      integer nposimp(3,mimp)
c
      real*8 vmadih(mimp)
      real*8 vmadi(mimp)
c -------------------------------------------
c
      common/lay2d/clay(mtotal,3),nextra,nbulkl,nbulkr,
     &             nprc,ninprc(0:mprc+1)
      common/test/itest
      common/xch/ixch
c
      write(6,'(/2x,''routine NEWPOT>''/)')
c
      pi=4.d0*datan(1.d0)
      fpi=4.d0*pi
      lmaxs=2*lmax
      lmmaxs=(lmaxs+1)*(lmaxs+1)
c
c charges for Madelung potential
c
      qif=0.d0
      do li=1,nintfc
        qa=qca(li)+qva(li)-za(li)
        qb=qcb(li)+qvb(li)-zb(li)
        qlay(li)=conc(li)*qa+(1.d0-conc(li))*qb
        qmomlay(1,li)=qlay(li)
        do i=2,lmmaxs
          qmomlay(i,li)=conc(li)*qmoma(i,li)+(1.d0-conc(li))*qmomb(i,li)
        end do
        if(itest.ge.2)
     >  write(6,'('' I'',i3,''  Layer'',i2,''  q_A ='',f15.10,
     >            ''  q_B ='',f15.10,''   q='',f15.10)') 
     >            itscf,li,qa,qb,qlay(li)
        qif=qif+qlay(li)
      end do
      write(6,'('' I'',i3,''  Interf.  q='',f15.10)') itscf,qif
c       
      if(.not.bulk) then
        do li=1,ninprc(0)
          do i=1,lmmaxs
            qmoml(i,li)=concl(li)*qmomla(i,li)+
     >                  (1.d0-concl(li))*qmomlb(i,li)
          end do 
        end do 
        if(rightm.eq.'B') then 
          do li=1,ninprc(nprc+1) 
            do i=1,lmmaxs 
              qmomr(i,li)=concr(li)*qmomra(i,li)+
     >                    (1.d0-concr(li))*qmomrb(i,li) 
            end do 
          end do 
        end if 
      else 
        do li=1,ninprc(0) 
          do i=1,lmmaxs 
            qmoml(i,li)=qmomlay(i,li) 
          end do 
        end do 
      end if
c       
c compute layer dependent Madelung potentials
c
c     ------------------------------------------------------------------
      call madelungb(qmoml,nintfc,lmax,sigma,bulk,'L',vmad,vleft,vright)
c     ------------------------------------------------------------------
      if(bulk) goto 10
c
      if(rightm.eq.'B')
c     ------------------------------------------------------------------
     >call madelungb(qmomr,nintfc,lmax,sigma,bulk,'R',vmad,vleft,vright)
c     ------------------------------------------------------------------
c
      if(rightm.eq.'V') then 
c
c       --------------------------------------------------------------
        call madelungv(qmoml,qmomlay,nintfc,lmax,sigma,vleft,vmad,vvac)
c       --------------------------------------------------------------
c
c print out new vacuum potential level 
c
c        if(itest.ge.1) write(6,'('' V0 for vacuum:   '',f15.10)') vvac
c
      else if(rightm.eq.'B') then
c
c       -----------------------------------------------------
        call madelungi(qmoml,qmomr,qmomlay,nintfc,lmax,sigma,
     >                 vleft,vright,vmad)
c       -----------------------------------------------------
c
      end if
  10  continue
c
c ----------------------------------------------------------------------
c Generate the real-space Madelung contants for cluster calculation
c
        fullmad=.true.  !DEBUG OPTION!
        write(6,*) '<newpot>: fullmad=',fullmad
        if(.not.fullmad) write(6,*) '<newpot>: ************************'
        if(.not.fullmad) write(6,*) '<newpot>: ******  CAUTION!  ******'
        if(.not.fullmad) write(6,*) '<newpot>: fullmad=F is a TEST opt.'
        if(.not.fullmad) write(6,*) '<newpot>: ************************'
c
        call madcolay(1,fullmad,nimp,nposimp,qmomlay,vmadi)  
c
        do iimp=1,nimp
          li=nposimp(3,iimp)
          vmadih(iimp)=vmad(li)-vmadi(iimp)
c 
          write(6,*) '<newpot>: vmadih',iimp,vmadih(iimp)
c
          if(itest.ge.2) then
            write(6,*) '<newpot>: vmad',li,vmad(li)
            write(6,*) '<newpot>: vmadi',iimp,vmadi(iimp)
          end if
c 
        end do
c ----------------------------------------------------------------------
c
c solve Poisson equation and
c put together muffin tin potential
c
      den=0.d0 !DEBUG
      den1a=0.d0 !DEBUG
      den2a=0.d0 !DEBUG
      den1b=0.d0 !DEBUG
      den2b=0.d0 !DEBUG
      if(itest.ge.1) write(6,*)
      do li=1,nintfc
        r1=rs(li)
        n1=ns(li) 
        dx1=dx(li)
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
          den=rhospa(i,1,li)/fpirr
          if(den.gt.0.d0) then
            den1a=den
          else
            den=den1a
          end if
          rhospa(i,1,li)=fpirr*den
          den=rhospa(i,2,li)/fpirr
          if(den.gt.0.d0) then
            den2a=den
          else
            den=den2a
          end if
          rhospa(i,2,li)=fpirr*den
          rhova(i,li)=rhospa(i,1,li)+rhospa(i,2,li)
c -B
          den=rhospb(i,1,li)/fpirr
          if(den.gt.0.d0) then
            den1b=den
          else
            den=den1b
          end if
          rhospb(i,1,li)=fpirr*den
          den=rhospb(i,2,li)/fpirr
          if(den.gt.0.d0) then
            den2b=den
          else
            den=den2b
          end if
          rhospb(i,2,li)=fpirr*den
          rhovb(i,li)=rhospb(i,1,li)+rhospb(i,2,li)
        end do
c
        x=x0
        do i=1,n1
          rgrid(i)=dexp(x)
          rhota(i)=rhoca(i,li)+rhova(i,li)
          rhotb(i)=rhocb(i,li)+rhovb(i,li)
          x=x+dx1
        end do
c       --------------------
        call rzero(pota,nrad)
        call rzero(potb,nrad)
c       --------------------
        qa=qca(li)+qva(li)
        qb=qcb(li)+qvb(li)
c       -----------------------------------------
        call poissond(rhota,qa,pota,dx1,x0,r1,n1)
        call poissond(rhotb,qb,potb,dx1,x0,r1,n1)
c       -----------------------------------------
c
c constant to be added to the electrostatic potential
c
        if(itest.ge.1) 
     >  write(6,'('' I'', i3,''  Mad. pot.'',i2,''  Vm='',f15.10)')
     >  itscf,li,vmad(li)
c
c electrostatic potential for layer li
c
        do i=1,n1
          vnucla=-2.d0*za(li)/rgrid(i)
          vnuclb=-2.d0*zb(li)/rgrid(i)
          rxa(i)=rhota(i)*(pota(i)+2.0d0*vnucla)
          rxb(i)=rhotb(i)*(potb(i)+2.0d0*vnuclb)
          pota(i)=pota(i)+vnucla+vmad(li)
          potb(i)=potb(i)+vnuclb+vmad(li)
        end do
c electrostatic energy for layer li
        enela(li)=rsimp(rxa,rgrid,n1,dx1)
        enelb(li)=rsimp(rxb,rgrid,n1,dx1)
        enela(li)=0.5d0*enela(li)
        enelb(li)=0.5d0*enelb(li)
c
        if(itest.ge.3) then
          write(6,'(/'' Layer'',i2)') li
          write(6,'('' Electrostatic potential'')')
          write(6,'(4d20.12)') (pota(i),i=1,n1)
          write(6,*)
          write(6,'(4d20.12)') (potb(i),i=1,n1)
        end if
c
c exchange-correlation potential and energy for layer li
        do i=1,n1
          r=rgrid(i)
          fpir2=fpi*r*r
c
          rho1=(0.5d0*rhoca(i,li)+rhospa(i,1,li))/fpir2
          rho2=(0.5d0*rhoca(i,li)+rhospa(i,2,li))/fpir2
          rho=rho1+rho2
c         -------------------------------------------------------
          call exchange(rho1,rho2,rho,ixch,xpot1(i),xpot2(i),exc)
c         -------------------------------------------------------
          vxc=0.5d0*(xpot1(i)+xpot2(i))
          bxc=0.5d0*(xpot2(i)-xpot1(i))
          if(dabs(bxc).lt.1.d-15) bxc=0.d0
          pota(i)=pota(i)+vxc
          vra(i,li)=r*pota(i)
          bra(i,li)=r*bxc
          rxa(i)=rhota(i)*exc
          rtmp(i)=rhota(i)*pota(i)
          rxb(i)=(rhospa(i,2,li)-rhospa(i,1,li))*bxc
        end do
        enxca(li)=rsimp(rxa,rgrid,n1,dx1)
        enpota(li)=-rsimp(rtmp,rgrid,n1,dx1)
        enmaga(li)=-rsimp(rxb,rgrid,n1,dx1)
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
        do i=1,n1
          r=rgrid(i)
          fpir2=fpi*r*r
c
          rho1=(0.5d0*rhocb(i,li)+rhospb(i,1,li))/fpir2
          rho2=(0.5d0*rhocb(i,li)+rhospb(i,2,li))/fpir2
          rho=rho1+rho2
c         -------------------------------------------------------
          call exchange(rho1,rho2,rho,ixch,xpot1(i),xpot2(i),exc)
c         -------------------------------------------------------
          vxc=0.5d0*(xpot1(i)+xpot2(i))
          bxc=0.5d0*(xpot2(i)-xpot1(i))
          if(dabs(bxc).lt.1.d-15) bxc=0.d0
          potb(i)=potb(i)+vxc
          vrb(i,li)=r*potb(i)
          brb(i,li)=r*bxc
          rxb(i)=rhotb(i)*exc
          rtmp(i)=rhotb(i)*potb(i)
          rxa(i)=(rhospb(i,2,li)-rhospb(i,1,li))*bxc
        end do
        enxcb(li)=rsimp(rxb,rgrid,n1,dx1)
        enpotb(li)=-rsimp(rtmp,rgrid,n1,dx1)
        enmagb(li)=-rsimp(rxa,rgrid,n1,dx1)
c
c  rxa has been used for temporary storage!
c
        if(itest.ge.3) then
          write(6,'('' B-type spin-up x-c potential'')')
          write(6,'(4d20.12)') (xpot1(i),i=1,n1)
          write(6,'('' B-type spin-down x-c potential'')')
          write(6,'(4d20.12)') (xpot2(i),i=1,n1)
        end if
c
        if(itest.ge.3) then
          write(6,'('' Total potential'')')
          write(6,'(4d20.12)') (pota(i),i=1,n1)
          write(6,'('' magnetic field'')')
          write(6,'(4d20.12)') (bra(i,li)/rgrid(i),i=1,n1)
          write(6,*)
          write(6,'(4d20.12)') (potb(i),i=1,n1)
          write(6,'('' magnetic field'')')
          write(6,'(4d20.12)') (brb(i,li)/rgrid(i),i=1,n1)
        end if
c
        if(orbpol) then
c
          call bop(rs(li),ns(li),dx(li),rhodspa(1,1,li),lza(1,li),
     >             bopra(1,1,li))
          call bop(rs(li),ns(li),dx(li),rhodspb(1,1,li),lzb(1,li),
     >             boprb(1,1,li))
c
          if(itest.ge.3) then
            write(6,'(/ ''LAYER '',i3)') li
            write(6,'(/ ''BOPR A SPIN 1'')')
            write(6,'(4d20.10)') (bopra(j,1,li),j=1,ns(li))
            write(6,'(/ ''BOPR A SPIN 2'')')
            write(6,'(4d20.10)') (bopra(j,2,li),j=1,ns(li))
            write(6,'(/ ''BOPR B SPIN 1'')')
            write(6,'(4d20.10)') (boprb(j,1,li),j=1,ns(li))
            write(6,'(/ ''BOPR B SPIN 2'')')
            write(6,'(4d20.10)') (boprb(j,2,li),j=1,ns(li))
          end if
c
        end if
c
      end do
c
      return
      end
