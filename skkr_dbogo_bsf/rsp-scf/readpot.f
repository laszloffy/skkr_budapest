c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine readpot(idpot,conc,vr,br,bopr,potin,nomag,opot,
     >    dx,ns,rsnew,z,potzero,delta,rlambda,singrat,urat,drat)
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
c         delta(r) = lambda*chi(r)
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      logical opot,potin
      character*10 idpot,head
      dimension vr(nrad),br(nrad),bopr(nrad,2)
      dimension rs(nrad),xs(nrad),vr0(0:nrad),br0(0:nrad),rs0(0:nrad)
      dimension bopr1(0:nrad),bopr2(0:nrad)
      real*8 delta_re,delta_im,concd
      dimension chi_re(nrad)
      dimension chi_im(nrad)
      complex*16 delta(nrad)
      integer idelta
      common/test/itest
      data tiny/1.0d-6/
c
      if(itest.ge.2) write(6,'('' < READPOT'')')
c 
c      read(22,*)  head,concd !!! Header added to delta file
c      if(itest.ge.2) write(6,*) head
c      if(itest.ge.2) write(6,'(f10.5)') concd
c      read(22,*) delta_re,delta_im,rlambda
c      if(itest.ge.2) write(6,*) delta_re,delta_im,rlambda
      vr0 = (0.d0)
      br0 = (0.d0)
      read(20,'(a10)') idpot
      if(itest.ge.2) write(6,*) idpot
c      write(6,*) idpot
      read(20,*) conc
      read(20,*) z,potzero,dx0,ns0,rws,idelta
c      write(6,*) delta
      if (idelta.eq.0) then 
         write(6,*) "Normal mode! > delta=0"
         do j=1,nrad
           delta(j) = CMPLX(0.d0,0.d0)
         end do
      end if
      singrat=1.0d0
      urat=0.0d0
      drat=0.0d0
      if (idelta.eq.1) then
         write(6,*) "Constant delta potential"
         read(20,*) chi_re(1),chi_im(1),rlambda   
         do j=1,nrad
           delta(j)=rlambda*CMPLX(chi_re(1),chi_im(1))
         end do
      end if
      if (idelta.eq.3) then
         write(6,*) "Constant delta potential with triplet"
         read(20,*) chi_re(1),chi_im(1),rlambda   
         read(20,*) singrat,urat,drat
         do j=1,nrad
           delta(j)=rlambda*CMPLX(chi_re(1),chi_im(1))
         end do
         write(6,*) "Singlet Delta ratio = ",singrat
         write(6,*) "Up-up Triplet Delta ratio = ",urat
         write(6,*) "Down-down Triplet Delta ratio = ",drat
      end if
      read(20,*) (vr0(j),j=1,ns0)
      if(potin) then
        do j=1,ns0
          vr0(j)=vr0(j)-2.0d0*z
        end do
      end if
      if(nomag.eq.0) then
        read(20,*) (br0(j),j=1,ns0)
        if(opot) then
          read(20,*) (bopr1(j),j=1,ns0)
          read(20,*) (bopr2(j),j=1,ns0)
        end if
      end if
      if (idelta.eq.2) then
         write(6,*) "Radial delta potential"
         read(20,*) rlambda
         read(20,*) (chi_re(j),j=1,ns0)
         read(20,*) (chi_im(j),j=1,ns0)
         do j=1,nrad
           delta(j)=rlambda*CMPLX(chi_re(j),chi_im(j))
         end do
      end if
c
      if(itest.ge.2) write(6,*) 'vr,br read succesful'
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
c
c -- interpolation --
c
      ns=ns0
      dx=dx0
c
      do j=1,ns
         xs(j)=dlog(rsnew)-(ns-j)*dx
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
      do j=ns+1,nrad
         xs(j)=dlog(rsnew)-(ns-j)*dx
         rs(j)=dexp(xs(j))
         vr(j)=0.0d0
         br(j)=0.0d0
         bopr(j,1)=0.0d0
         bopr(j,2)=0.0d0
      end do
c
c
c  ----- generate the delta(r) function
c
c      do j=1,nrad
c        write(6,*) j
c        delta(j)=rlambda*CMPLX(delta_re,delta_im)
c      end do
c
      if(itest.le.2) return 
      write(6,'(a10)') idpot
      write(6,'(f10.5)') conc
      write(6,'(f6.2,f10.5,f7.3,i5,f15.10)') z,potzero,dx,ns,rsnew
      write(6,*) 'Delta'
      write(6,'(4d20.10)') (delta(j),j=1,ns)
      write(6,*) 'vr'
      write(6,'(4d20.10)') (vr(j),j=1,ns)
      write(6,*) 'br'
      write(6,'(4d20.10)') (br(j),j=1,ns)
      write(6,*) 'bopr1'
      write(6,'(4d20.10)') (bopr(j,1),j=1,ns)
      write(6,*) 'bopr2'
      write(6,'(4d20.10)') (bopr(j,2),j=1,ns)
c
      return
      end
