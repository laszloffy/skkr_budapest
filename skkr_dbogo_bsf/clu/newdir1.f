c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine newdir1(
c     =================
     > itscf,itscfcont,nintfc,lfix1,lfix2,
     > rba,vecna,phia,lmax,
     > spin_magvpa,spin_magva,
     > orb_magvpa,orb_magva,lza,
     > spinmoma,orbmoma,th0a,th1a,ph0a,ph1a)
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      integer mmax
      parameter(mmax=mimp)
c
      parameter(mbrdir=4*mmax)
c
      real*8 lza(2,mmax)
      real*8 rhospa(nrad,2,mmax)
      real*8 rhodspa(nrad,2,mmax)
      real*8 bopra(nrad,2,mmax)
c
      real*8 spin_magvpa(kmymaxp,mmax,3)
      real*8 spin_magva(mmax,3),spin_magvb(mmax,3)
      real*8 orb_magvpa(kmymaxp,mmax,3)
      real*8 orb_magva(mmax,3)
c
      real*8 rba(3,mmax)
      real*8 vecna(3,mmax)
      real*8 phia(mmax)
      real*8 vec1(3)
      real*8 vec2(3)
c
      real*8 spinmoma(mmax)
      real*8 orbmoma(mmax)
      real*8 th0a(mmax)
      real*8 th1a(mmax)
      real*8 th2a(mmax)
      real*8 ph0a(mmax)
      real*8 ph1a(mmax)
      real*8 ph2a(mmax)
c
      real*8 xbr(mbrdir)
      real*8 fbr(mbrdir)
      real*8 dfbr(mbrdir)
      real*8 pfbr(mbrdir)
      real*8 pxbr(mbrdir)
      real*8 denbr(mbrdir)
      real*8 defbr(mbrdir)
c
      character*3  tempr
      character*50 temp1
      character*50 temp2
      character*50 temp3
      integer root, myrank, nprocs
      common/mpi/root,myrank,nprocs

c
      common/broydir/dmix,wbrdir,nbrdir
c
      data tol/1.0d-10/
c
      save ibr,ibrmax,ibrcase
c
      lmmax=(lmax+1)*(lmax+1)
c
c calculate first spin-resolution of z-component of
c the d-like orbital moments in the (old) local frame of reference
c
      do li=1,nintfc
         lza(1,li)=0.d0
         lza(2,li)=0.d0
!        lzb(1,li)=0.d0
!        lzb(2,li)=0.d0
         do lm=5,9
           lza(1,li)=lza(1,li)+orb_magvpa(lm,li,3)
           lza(2,li)=lza(2,li)+orb_magvpa(lm+lmmax,li,3)
!          lzb(1,li)=lzb(1,li)+orb_magvpb(lm,li,3)
!          lzb(2,li)=lzb(2,li)+orb_magvpb(lm+lmmax,li,3)
         end do
      end do
c
c calculate spin- and orbital moments in the global frame
c of reference
c
      do li=1,nintfc
         do i=1,3
           vec1(i)=spin_magva(li,i)
         end do
         call rotvec(vec1,vec2,vecna(1,li),phia(li),1)
         do i=1,3
           spin_magva(li,i)=vec2(i)
         end do
c
!        do i=1,3
!          vec1(i)=spin_magvb(li,i)
!        end do
!        call rotvec(vec1,vec2,vecnb(1,li),phib(li),1)
!        do i=1,3
!          spin_magvb(li,i)=vec2(i)
!        end do
c
         do i=1,3
           vec1(i)=orb_magva(li,i)
         end do
         call rotvec(vec1,vec2,vecna(1,li),phia(li),1)
         do i=1,3
           orb_magva(li,i)=vec2(i)
         end do
c
!        do i=1,3
!          vec1(i)=orb_magvb(li,i)
!        end do
!        call rotvec(vec1,vec2,vecnb(1,li),phib(li),1)
!        do i=1,3
!          orb_magvb(li,i)=vec2(i)
!        end do
      end do
c
c calculate theta and phi angles of old and new spin-directions
c
      do li=1,nintfc
         x0=rba(1,li)
         y0=rba(2,li)
         z0=rba(3,li)
         call findangles(x0,y0,z0,theta0,phi0)
         th0a(li)=theta0
         ph0a(li)=phi0
         x1=spin_magva(li,1)
         y1=spin_magva(li,2)
         z1=spin_magva(li,3)
         rmod=dsqrt(x1*x1+y1*y1+z1*z1)
         spinmoma(li)=rmod
         orbmoma(li)=dsqrt(orb_magva(li,1)*orb_magva(li,1)
     >                    +orb_magva(li,2)*orb_magva(li,2)
     >                    +orb_magva(li,3)*orb_magva(li,3))
         if(rmod.gt.tol) then
           x1=x1/rmod
           y1=y1/rmod
           z1=z1/rmod
           call findangles(x1,y1,z1,theta1,phi1)
           th1a(li)=theta1
           ph1a(li)=phi1
         else  
           th1a(li)=th0a(li)
           ph1a(li)=ph0a(li)
         end if
c
!        x0=rbb(1,li)
!        y0=rbb(2,li)
!        z0=rbb(3,li)
!        call findangles(x0,y0,z0,theta0,phi0)
!        th0b(li)=theta0
!        ph0b(li)=phi0 
!        x1=spin_magvb(li,1)
!        y1=spin_magvb(li,2)
!        z1=spin_magvb(li,3)
!        rmod=dsqrt(x1*x1+y1*y1+z1*z1)
!        spinmomb(li)=rmod
!        orbmomb(li)=dsqrt(orb_magvb(li,1)*orb_magvb(li,1)
!    >                    +orb_magvb(li,2)*orb_magvb(li,2)
!    >                    +orb_magvb(li,3)*orb_magvb(li,3)) 
!        if(rmod.gt.tol) then
!          x1=x1/rmod
!          y1=y1/rmod
!          z1=z1/rmod
!          call findangles(x1,y1,z1,theta1,phi1)
!          th1b(li)=theta1
!          ph1b(li)=phi1
!        else
!          th1b(li)=th0b(li)
!          ph1b(li)=ph0b(li)
!        end if
      end do
c
c and mix them to get the next guess for the spin-direction
c
      i=0
      do li=1,nintfc
        if(lfix1.eq.li.or.lfix2.eq.li) goto 10
          i=i+1
          xbr(i)=th0a(li)
          fbr(i)=th1a(li)-th0a(li)
!         i=i+1
!         xbr(i)=th0b(li)
!         fbr(i)=th1b(li)-th0b(li)
          i=i+1
          xbr(i)=ph0a(li)
          fbr(i)=ph1a(li)-ph0a(li)
!         i=i+1
!         xbr(i)=ph0b(li)
!         fbr(i)=ph1b(li)-ph0b(li)
  10    continue
      end do
      nbrdim=i
c
      if(dmix.lt.tol) goto 100
c
      call dott(xbr,xbr,xnorm2,nbrdim)
      call dott(fbr,fbr,err2,nbrdim)
      if(xnorm2.lt.tol) then
        ferr=0.d0
      else
        ferr=dsqrt(err2/xnorm2)
      end if
      write(6,'(/'' Error for spin directions :'',3d17.6)')  
     >err2,xnorm2,ferr
c
      ibrmax=nbrdir
      if(itscf.eq.1) then
        ibrcase=1
        ibr=1
      end if
      write(6,'('' itscf,ibrcase,ibrmax,ibr:'',4i4)')
     >itscf,ibrcase,ibrmax,ibr
c
      write(tempr,'(i3.3)') myrank

      temp1 = 'ftn31.'//tempr
      temp2 = 'ftn32.'//tempr
      temp3 = 'ftn33.'//tempr
c     ----------------------------------------------------------------
      call broyd(xbr,fbr,0.d0,ibr,nbrdim,dfbr,pfbr,pxbr,denbr,defbr,
     >           dmix,wbrdir,31,32,33,temp1,temp2,temp3)
c     ----------------------------------------------------------------
      if(ibr.ge.ibrmax) then
        ibr=1
        ibrcase=ibrcase+1
      else
        ibr=ibr+1
      end if
c
  100 continue
      i=0
      do li=1,nintfc
        if(li.ge.lfix1.and.li.le.lfix2) then    
          th1a(li)=th0a(li)
          ph1a(li)=ph0a(li)
!         th1b(li)=th0b(li)
!         ph1b(li)=ph0b(li)
        else
          i=i+1
          th1a(li)=xbr(i)
!         i=i+1
!         th1b(li)=xbr(i)
          i=i+1
          ph1a(li)=xbr(i)
!         i=i+1
!         ph1b(li)=xbr(i)
        end if        
        rba(1,li)=dsin(th1a(li))*dcos(ph1a(li))
        rba(2,li)=dsin(th1a(li))*dsin(ph1a(li))
        rba(3,li)=dcos(th1a(li))
!       rbb(1,li)=dsin(th1b(li))*dcos(ph1b(li))
!       rbb(2,li)=dsin(th1b(li))*dsin(ph1b(li))
!       rbb(3,li)=dcos(th1b(li))
      end do
c
      return
      end
