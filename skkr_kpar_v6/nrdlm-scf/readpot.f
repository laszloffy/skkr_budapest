c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine readpot(idpot,conc,vr,br,dxnew,nsnew,rsnew,z,potzero,
     >                   imag,rspkkrin)
c=======================
c
c read in spherical potential and interpolate to W-S radius
c
c input   dxnew,nsnew,rsnew - logarithmic increment
c                             number of radial points
c                             wigner-seitz radius
c         imag  - read (or not) input magnetic field
c
c output  idpot - identifier
c         conc - concentration
c         vr - potential*r
c         br - effective field*r
c         dx - logarithmic increment
c         ns - number of radial points    
c         z  - atomic number
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      logical opot,rspkkrin
      character*10 idpot
      dimension vr(nrad),br(nrad)
      dimension rs(nrad),xs(nrad),vr0(0:nrad),br0(0:nrad),rs0(0:nrad)
      common/test/itest
c
      if(itest.ge.3) write(6,'('' < READPOT'')')
c
      ctol=1.0d8
c
c     write(6,*) 'imag=',imag 
      read(20,'(a10)') idpot
c     write(6,'(a10)') idpot
      read(20,*) conc
c     write(6,*) conc
      read(20,*) z,potzero,dx0,ns0,rws
c     write(6,*) z,potzero,dx0,ns0,rws
      if(ns0.gt.nrad) stop
     &     'Increase parameter nrad according to I/pot.mesh'
      read(20,*) (vr0(j),j=1,ns0)
c     write(6,*) (vr0(j),j=1,ns0)
      if(imag.eq.1) read(20,*) (br0(j),j=1,ns0)
c
      rs0(0)=0.d0
      vr0(0)=-2.d0*z
      br0(0)=0.0d0
      xws=dlog(rws)
      do j=1,ns0
        x=xws-(ns0-j)*dx0
        rs0(j)=dexp(x)
      end do
      if(rspkkrin) vr0(1:ns0)=vr0(1:ns0)*rs0(1:ns0)
      if(rspkkrin) br0(1:ns0)=br0(1:ns0)*rs0(1:ns0)
c
      do j=1,nsnew
         xs(j)=dlog(rsnew)-(nsnew-j)*dxnew
         rs(j)=dexp(xs(j))
         vr(j)=ylag(rs(j),rs0(0),vr0(0),0,3,ns0+1,iex)
         if(imag.eq.1) br(j)=ylag(rs(j),rs0(0),br0(0),0,3,ns0+1,iex)
      end do
      do j=nsnew+1,nrad
         xs(j)=dlog(rsnew)-(nsnew-j)*dxnew
         rs(j)=dexp(xs(j))
         vr(j)=0.0d0
         br(j)=0.0d0 
      end do
c
      if(itest.le.2) return 
      write(6,'(a10)') idpot
      write(6,'(f10.5)') conc
      write(6,'(f6.2,f10.5,f7.3,i5,f15.10)') 
     >z,potzero,dxnew,nsnew,rsnew
      write(6,'(4d20.10)') (vr(j),j=1,nsnew)
      if(imag.eq.1) write(6,'(4d20.10)') (br(j),j=1,nsnew)
c
      return
      end
