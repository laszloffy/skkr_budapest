c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine readpot(idpot,conc,lmax,qmom,vr,br,bopr,opot,
     >                   dx,ns,rsnew,z,for006,imom)
c=======================
c
c read in spherical potential
c
c input   rsnew - the actual (wigner-seitz) radius;
c                 the number of points ns and the log. increment dx
c                 stay the same as in the read potential (ns0,dx0) !!
c output  idpot - identifier
c         conc - concentration
c         vr - potential*r
c         br - effective field*r
c         bopr - vector potential*r
c         opot - if 'true' read in bopr
c         dx - logarithmic increment
c         ns - number of radial points    
c         z  - atomic number
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      logical opot
      character*10 idpot
      character*30 for006
      dimension vr(nrad),br(nrad),bopr(nrad,2)
      dimension rs(nrad),xs(nrad),vr0(0:nrad),br0(0:nrad),rs0(0:nrad)
      dimension bopr1(0:nrad),bopr2(0:nrad)
      dimension qq(2*lmsup)
      complex*16 qmom(lmsup)
      common/test/itest
      data tiny/1.0d-6/
c
      if(itest.ge.2) write(6,'('' < READPOT'')')
c
      lmaxs=2*lmax
      lmmaxs=(lmaxs+1)*(lmaxs+1)
c
      read(20,'(a10)') idpot
      read(20,*) conc
      if(imom.eq.1) read(20,*) (qq(i),i=1,2*lmmaxs)
      do i=1,lmmaxs
        qmom(i)=dcmplx(qq(2*i-1),qq(2*i))
      end do
      read(20,*) z,potzero,dx0,ns0,rws
      read(20,*) (vr0(j),j=1,ns0)
      read(20,*) (br0(j),j=1,ns0)
      if(opot) then
        read(20,*) (bopr1(j),j=1,ns0)
        read(20,*) (bopr2(j),j=1,ns0)
      else
        call rzero(bopr,2*nrad)
      end if
c
      rs0(0)=0.d0
      vr0(0)=-2.d0*z
      br0(0)=0.0d0
      bopr1(0)=0.0d0
      bopr2(0)=0.0d0
      xws=dlog(rws)
      do j=1,ns0
        x=xws-(ns0-j)*dx0
        rs0(j)=dexp(x)
        vr0(j)=vr0(j)-potzero*rs0(j)
      end do
      potzero=0.d0
c
c -- initialize the interpolation --
c
      ns=ns0
      dx=dx0
c
      do j=1,ns
         xs(j)=dlog(rsnew)-(ns-j)*dx
         rs(j)=dexp(xs(j))
         vr(j)=ylag(rs(j),rs0(0),vr0(0),0,3,ns0+1,iex)
         br(j)=ylag(rs(j),rs0(0),br0(0),0,3,ns0+1,iex)
         if(opot) then
            bopr(j,1)=ylag(rs(j),rs0(0),bopr1(0),0,3,ns0+1,iex)
            bopr(j,2)=ylag(rs(j),rs0(0),bopr2(0),0,3,ns0+1,iex)
         endif
      end do
c
      if(itest.le.2) return 
      write(6,'(a10)') idpot
      write(6,'(f10.5)') conc
      write(6,'(f6.2,f10.5,f7.3,i5,f15.10)') z,potzero,dx,ns,rsnew
      write(6,'(4d20.10)') (vr(j),j=1,ns)
      write(6,'(4d20.10)') (br(j),j=1,ns)
      write(6,'(4d20.10)') (bopr(j,1),j=1,ns)
      write(6,'(4d20.10)') (bopr(j,2),j=1,ns)
c
      return
      end
