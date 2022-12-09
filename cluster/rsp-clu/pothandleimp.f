c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
c***********************************************************************
c Now, the radial mesh of the potential within a layer is chosen as the 
c finest grid between both CPA species.
c
c For bulk calculation, assumed: ninprcl=ninprcr=ninprc(1)  &  nprc=1 
c***********************************************************************
      subroutine pothandleimp(
     > lmax,nintfc,v00,ivacpot,for006,opot,laypot,dx,ns,rs,idpota,
     > vra,bra,bopra,za,ib0,b0,layb0,layb0f,isigba)
c
c -read Left,Right and Layer potentials
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
      integer mmax
      parameter (mmax=mimp)
c
      logical opot
      character*10 idpot0
      character*10 idpota(mmax)
      character*30 laypot
      character*30 for006
c
      real*8 vra(nrad,mmax)
      real*8 bra(nrad,mmax)
      real*8 bopra(nrad,2,mmax)
      real*8 dx(mmax)
      real*8 rs(mmax)
      real*8 za(mmax)
c
      integer isigba(mmax)
      integer ns(mmax)
c
      complex*16 qmom(lmsup)
c
      real*8 potshift
      real*8 b0
c
      integer layb0
      integer layb0f
      integer ib0
      integer imom
c
      data tiny/1.0d-6/
c
      lmaxs=2*lmax
      lmmaxs=(lmaxs+1)*(lmaxs+1)
c
      imom=0
c -- rsl/~r is the Wigner-Seitz radius for the left and right bulk and
c    rs(i),i=1,nintfc is the WS-radius for the interface !!
c
      open(20,file=laypot,status='old')
      read(20,*)
      write(6,*) 'imppot=',laypot
      do li=1,nintfc
        call readpot(idpot0,conc0,lmax,qmom,vra(1,li),bra(1,li),
     >       bopra(1,1,li),opot,dx(li),ns(li),rs(li),za(li),
     >       for006,imom)
c
        idpota(li)=idpot0
c
c       if((dabs(dx1-dx(li)).gt.tiny).or.(ns1.ne.ns(li))) then
c         write(*,*) 
c    &    ' WARNING - POTHANDLE: radial grid mismatch in layer ',li
c         ns0 = max(ns1,ns(li))
c         dx0 = min(dx1,dx(li))
c         call interpot(ns(li),dx(li),ns0,dx0,za(li),rs(li),
c    >                  vra(1,li),bra(1,li),bopra(1,1,li),opot)
c         ns(li) = ns0
c         dx(li) = dx0
c       endif
c
        if((layb0.eq.0).or.((li.ge.layb0).and.(li.le.layb0f))) then
          if(ib0.eq.0) then
            x=dlog(rs(li))-(ns(li)-1)*dx(li)
            do i=1,ns(li)
              r=dexp(x)
              bra(i,li)=bra(i,li)+r*b0
              x=x+dx(li)
            end do
          else
            x=dlog(rs(li))-(ns(li)-1)*dx(li)
            do i=1,ns(li)
              r=dexp(x)
              bra(i,li)=r*b0
              x=x+dx(li)
            end do
          end if
        end if
c
        if(isigba(li).ne.1) then
          do i=1,ns(li)
            bra(i,li)=-bra(i,li)
          end do
          write(6,'(''<pothandleimp> : Sign of effective 
     >          field changed'')') li
        end if
c
        write(6,'(2x,a3)')  idpota(li)
c
      end do
      close(20)
c
      write(6,'(''<pothandleimp> : pothandleimp processed!'')')
c
c
c     if(.not.bulk) return
c
c limitation to 3D periodicity:
c
c     do li=1,nintfc
c       idpotla(li)=idpota(li)
c       idpotra(li)=idpota(li)
c       idpotlb(li)=idpotb(li)
c       idpotrb(li)=idpotb(li)
c       dxl(li)=dx(li)
c       dxr(li)=dx(li)
c       rsl(li)=rs(li)
c       rsr(li)=rs(li)
c       nsl(li)=ns(li)
c       nsr(li)=ns(li)
c       zla(li)=za(li)
c       zra(li)=za(li)
c       zlb(li)=zb(li)
c       zrb(li)=zb(li)
c       concl(li)=conc(li)
c       concr(li)=conc(li)
c     enddo
c
      return
      end
c===============
      subroutine interpot(ns0,dx0,ns,dx,z,rs,vr,br,bopr,opot)
c
c Interpolate input potential (in grid defined by ns0,dx0) to ns,dx.
c
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
      logical opot
      dimension vr(nrad),br(nrad),bopr(nrad,2)
      dimension vr0(0:nrad),br0(0:nrad),bopr1(0:nrad),bopr2(0:nrad),
     &          rs0(0:nrad)
c
      rs0(0) = 0.d0
      vr0(0) = -2.d0*z
      br0(0) = 0.d0
      bopr1(0) = 0.d0
      bopr2(0) = 0.d0
      do j=1,ns0
        xws=dlog(rs)-(ns0-j)*dx0
        rs0(j)=dexp(xws)
        vr0(j) = vr(j)
        br0(j) = br(j)
        if(opot) then
          bopr1(j) = bopr(j,1)
          bopr2(j) = bopr(j,2)
        endif
      enddo
c
      do j=1,ns
        xws=dlog(rs)-(ns-j)*dx
        rws=dexp(xws)
        vr(j)=ylag(rws,rs0(0),vr0(0),0,3,ns0+1,iex)
        br(j)=ylag(rws,rs0(0),br0(0),0,3,ns0+1,iex)
        if(opot) then
          bopr(j,1)=ylag(rws,rs0(0),bopr1(0),0,3,ns0+1,iex)
          bopr(j,2)=ylag(rws,rs0(0),bopr2(0),0,3,ns0+1,iex)
        endif
      enddo
c
      return
      end
