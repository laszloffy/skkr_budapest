c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
c***********************************************************************
c REVISE OBLIQUE LATTICE !!!!!!!!!!!!
c***********************************************************************
      subroutine  speck (spkx,spky,w,nspk,lat)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Routine to determine sine sets of special points in an irreducible
c  sector of a 2-d Brillouin Zone.
c  Determines the smallest number of special points greater than or
c  equal to the requested nspk.
c  Sum of (equal) weights normalized to unity.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit real*8 (a-h,o-z)
      parameter (nspkm=136)
c
      character*4 iogroup
      dimension spkx(nspk),spky(nspk),w(nspk)
      dimension ospkx(nspkm),ospky(nspkm),x(2),y(2),hy(6),hx(6),ir(2)
c
      common/cunnbz/ rlvax,rlvay,rlvbx,rlvby,phi,iogroup
c
      data (hx(i),i=1,6) /.1e1,.50,-.50,-.1e1,-.50,.50/
      data (ir(i),i=1,2) /1,0/
      data tol /1.0d-6/
c
      sr3=sqrt(3.0d+00)
      sr3d2=0.5*sr3
      hy(1)=0.0
      hy(2)=sr3d2
      hy(3)=sr3d2
      hy(4)=0.0
      hy(5)=-sr3d2
      hy(6)=-sr3d2
      stol=tol
      if ((lat.lt.1).or.(lat.gt.5)) then
      write (6,*) ' error: speck invalid lattice type.'
      stop
      endif
      if (nspk.gt.nspkm) then
      call error (nspk,nspkm,'nspkm in speck(0)')
      stop
      endif
      go to (10,60,60,120,160), lat
c
c----- special points for oblique lattice
c
   10 delta=-rlvay/rlvby
      beta=rlvax/rlvby
      q=1.0/sqrt(delta*delta+beta*beta)
      x(2)=rlvax/4.0
      y(2)=0.0
      x(1)=0.0
      y(1)=rlvby/4.0
      spkx(1)=x(1)
      spky(1)=y(1)
      w(1)=1.0
      if (nspk.eq.1) return
      nt=int(log(dfloat(nspk))/log(2.0)+1.0-stol)
      nspk=2**(nt+1)/2
      if (nspk.gt.nspkm) then
      call error (nspk,nspkm,'nspkm in speck(1)')
      stop
      endif
      nsp=1
   20 do 30 i=1,nsp
      ospkx(i)=spkx(i)
   30 ospky(i)=spky(i)
      ic=int(dfloat(ir(2)+1)*q)
      if (ir(1).lt.ic) then
      nc=1
      else
      nc=2
      endif
      ir0=ir(nc)+1
      ir(nc)=ir0
      fac=2.0/(2.0**ir0)
      xkg=x(nc)*fac
      ykg=y(nc)*fac
      do 40 i=1,nsp
      j=2*i
      spkx(j-1)=ospkx(i)+xkg
      spky(j-1)=ospky(i)+ykg
      spkx(j)=ospkx(i)-xkg
   40 spky(j)=ospky(i)-ykg
      nsp=nsp*2
      if (nsp.lt.nspk) go to 20
      alpha=1.0/dfloat(nsp)
      do 50 i=1,nsp
   50 w(i)=alpha
      return
c
c----- special points for centred and primitive rectangular lattices
c
   60 xkg=rlvax/2.0
      ykg=rlvby/8.0
      if (lat.eq.3) then
      xkg=xkg/2.0
      ykg=ykg*2.0
      endif
      nsp=1
      spkx(1)=xkg
      spky(1)=ykg
      w(1)=1.0
      if (nspk.eq.1) return
      nt=int(log(dfloat(nspk))/log(4.0)+1.0-stol)
      nspk=4**(nt+1)/4
      if (nspk.gt.nspkm) then
      call error (nspk,nspkm,'nspkm in speck(2,3)')
      stop
      endif
   70 do 80 i=1,nsp
      ospkx(i)=spkx(i)
   80 ospky(i)=spky(i)
      xkg=xkg/2.0
      ykg=ykg/2.0
      do 90 i=1,nsp
      j=i*4
      x(1)=ospkx(i)+xkg
      y(1)=ospky(i)+ykg
      x(2)=ospkx(i)-xkg
      y(2)=ospky(i)-ykg
      spkx(j-3)=x(1)
      spky(j-3)=y(1)
      spkx(j-2)=x(1)
      spky(j-2)=y(2)
      spkx(j-1)=x(2)
      spky(j-1)=y(1)
      spkx(j)=x(2)
   90 spky(j)=y(2)
      nsp=nsp*4
      if (nsp.lt.nspk) go to 70
      alpha=1.0/dfloat(nsp)
      do 100 i=1,nsp
  100 w(i)=alpha
      if (lat.eq.3) return
      beta=-rlvax/rlvay
      do 110 i=1,nsp
      t=spky(i)+0.5*spkx(i)/beta-0.75*rlvax/beta
      if (t.le.0.0) go to 110
      spkx(i)=-spkx(i)+rlvax
      spky(i)=-spky(i)-rlvay
  110 continue
      return
c
c----- special points for square lattice
c
  120 xkg=rlvax/4.0
      ykg=rlvax/4.0
      spkx(1)=xkg
      spky(1)=ykg
      w(1)=1.0
      if (nspk.le.1) return
      nt=int(log((-1.0+sqrt(1.0+8.0*dfloat(nspk)))/2.0)/log(2.0)+1.0-
     1 stol)
      nu=2**(nt+1)/2
      nspk=(nu*nu+nu)/2
      if (nspk.gt.nspkm) then
      call error (nspk,nspkm,'nspkm in speck(4)')
      stop
      endif
      insp=0
  130 nt=(2**(insp+1))/2
      nsp=(nt*nt+nt)/2
      if (nsp.ge.nspk) return
      alpha=0.5/(2.0*nt*nt)
      do 140 i=1,nsp
      ospkx(i)=spkx(i)
  140 ospky(i)=spky(i)
      xkg=xkg/2.0
      ykg=ykg/2.0
      j=1
      do 150 i=1,nsp
      x(1)=ospkx(i)+xkg
      y(1)=ospky(i)+ykg
      x(2)=ospkx(i)-xkg
      y(2)=ospky(i)-ykg
      do 150 k1=1,2
      do 150 k2=1,2
      if (x(k1).lt.y(k2)) go to 150
      spkx(j)=x(k1)
      spky(j)=y(k2)
      w(j)=alpha
      if (x(k1).gt.y(k2)) w(j)=2.0*alpha
      j=j+1
  150 continue
      insp=insp+1
      go to 130
c
c----- special points for hexagonal lattice
c
  160 sa=0.5*rlvax
      x(2)=2.0*sa/3.0
      y(2)=x(2)/sr3
      x(1)=4.0*sa/3.0
      y(1)=0.0
      spkx(1)=x(2)
      spky(1)=y(2)
      w(1)=1.0
      if (nspk.le.1) return
      if (nspk.gt.nspkm) then
      call error (nspk,nspkm,'nspkm in speck(5)')
      stop
      endif
      nsp=1
      alpha=1.0
      in=2
  170 in=2-mod(in+1,2)
      x(in)=x(in)/3.0
      y(in)=y(in)/3.0
      xkg=x(in)
      ykg=y(in)
      do 180 i=1,nsp
      ospkx(i)=spkx(i)
  180 ospky(i)=spky(i)
      k=1
      do 210 i=1,nsp
      xa=ospkx(i)
      ya=ospky(i)
      do 210 j=1,6
      jm=mod(j,6)+1
      xap=hx(jm)*xkg+hy(jm)*ykg+xa
      yap=-hy(jm)*xkg+hx(jm)*ykg+ya
      if (yap.lt.-tol) then
       goto 210
      elseif ((yap-xap/sr3).gt.tol) then
       goto 210
      elseif ((yap+xap*sr3-2.0d+00/sr3*rlvax).gt.tol) then
       goto 210
      endif
      if (k.le.1) go to 200
      do 190 kc=1,k-1
      if ((abs(xap-spkx(kc)).lt.tol).and.(abs(yap-spky(kc)).lt.tol)) go
     1 to 210
  190 continue
  200 spkx(k)=xap
      spky(k)=yap
      k=k+1
  210 continue
      nsp=k-1
      alpha=alpha*3.0
      if (nsp.lt.nspk) go to 170
      nspk=nsp
      do 220 i=1,nsp
      w(i)=2.0/alpha
      xa=spkx(i)
      ya=spky(i)
      if ((abs(ya).lt.tol).or.(abs(ya-xa/sr3).lt.tol).or.(abs(ya+sr3*xa-
     1 2.00/sr3*rlvax).lt.tol)) w(i)=w(i)/2.0
  220 continue
      return
      end
