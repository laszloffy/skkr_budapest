c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine readpot(idpot,conc,vr,br,bopr,nomag,opot,
     >                   dxnew,nsnew,rsnew,z,potzero,rspkkrin)
c=======================
c
c read in spherical potential
c
c input   dxnew,nsnew,rsnew - logarithmic increment 
c                             number of radial points
c                             wigner-seitz radius
c output  idpot - identifier
c         conc - concentration
c         vr - potential*r
c         br - effective field*r
c         bopr - vector potential*r
c         opot - if 'true' read in bopr
c         z  - atomic number
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      logical opot,rspkkrin
      character*10 idpot
      dimension vr(nrad),br(nrad),bopr(nrad,2)
      dimension rs(nrad),xs(nrad),vr0(0:nrad),br0(0:nrad),rs0(0:nrad)
      dimension bopr1(0:nrad),bopr2(0:nrad)
      common/test/itest
      data tiny/1.0d-6/
c
      if(itest.ge.2) write(6,'('' < READPOT'')')
c
      read(20,'(a10)') idpot
      read(20,*) conc
      read(20,*) z,potzero,dx0,ns0,rws
      read(20,*) (vr0(j),j=1,ns0)
      if(nomag.eq.0) then
        read(20,*) (br0(j),j=1,ns0)
        if(opot) then
          read(20,*) (bopr1(j),j=1,ns0)
          read(20,*) (bopr2(j),j=1,ns0)
        end if
      end if
c
      rs0(0)=0.0d0
      vr0(0)=-2.0d0*z
      br0(0)=0.0d0
      bopr1(0)=0.0d0
      bopr2(0)=0.0d0
      xws=dlog(rws)
      do j=1,ns0
        x=xws-(ns0-j)*dx0
        rs0(j)=dexp(x)
      end do
      if(rspkkrin) vr0(1:ns0)=vr0(1:ns0)*rs0(1:ns0)
c
c -- interpolation --
c
      do j=1,nsnew
         xs(j)=dlog(rsnew)-(nsnew-j)*dxnew
         rs(j)=dexp(xs(j))
         vr(j)=ylag(rs(j),rs0(0),vr0(0),0,3,ns0+1,iex)
         if(nomag.eq.1) then
           br(j)=0.0d0
           bopr(j,1)=0.0d0
           bopr(j,2)=0.0d0
         else   
           br(j)=ylag(rs(j),rs0(0),br0(0),0,3,ns0+1,iex)
           if(opot) then
             bopr(j,1)=ylag(rs(j),rs0(0),bopr1(0),0,3,ns0+1,iex)
             bopr(j,2)=ylag(rs(j),rs0(0),bopr2(0),0,3,ns0+1,iex)
           else
             bopr(j,1)=0.0d0
             bopr(j,2)=0.0d0
           endif
         endif
      end do
c
      do j=nsnew+1,nrad
         xs(j)=dlog(rsnew)-(nsnew-j)*dxnew
         rs(j)=dexp(xs(j))
         vr(j)=0.0d0
         br(j)=0.0d0
         bopr(j,1)=0.0d0
         bopr(j,2)=0.0d0
      end do
c
      if(itest.le.2) return 
      write(6,'(a10)') idpot
      write(6,'(f10.5)') conc
      write(6,'(f6.2,f10.5,f7.3,i5,f15.10)') 
     >z,potzero,dxnew,nsnew,rsnew
      write(6,'(4d20.10)') (vr(j),j=1,nsnew)
      write(6,'(4d20.10)') (br(j),j=1,nsnew)
      write(6,'(4d20.10)') (bopr(j,1),j=1,nsnew)
      write(6,'(4d20.10)') (bopr(j,2),j=1,nsnew)
c
      return
      end
