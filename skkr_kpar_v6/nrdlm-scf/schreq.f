c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine outwrd (rel,e,l,vr,dx,xnot,rn,n,a,ap,logder)
c ====================
c     ******************************************************************
c     ***   performs the outward integration.  the runge-kutta       ***
c     ***   procedure is used to provide the starting values for     ***
c     ***   the adams integration.                                   ***
c     ***   spin-orbitless relativistic equations used.              ***
c     ***   see harmon & koelling, j.phys.c 10,2107(1977)            ***
c     ******************************************************************
c
c  input:  rel  - logical variable for non-relativistic or
c                 scalar-relativistic calculation
c          e    - complex*16 energy variable
c          l    - angular momentum quantumnumber
c          vr   - array for r*potential
c          dx   - logarithmic increment for radial scale
c          xnot - first point on logarithmic mesh
c          rn   - radius at which the logarithmic derivative has to be
c                 computed
c          n    - index corresponding to rn
c
c  output: a      - array for r*wavefuction
c          ap     - array for derivative of a
c          derlog - logarithmic derivative at rn
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
      parameter (nr=nrad)
c
      logical rel
      dimension vr(nr)
      complex*16 a(nr),b(nr),ap(nr),bp(nr),e,logder
      complex*16 pi,qi,remv,br,k1,m1,k2,m2,k3,m3,k4,m4
      complex*16 p0,q0,p1,q1,apip1,bpip1,asum,bsum
      data tol/1.0d-10/,itmax/100/,astart/1.0d-20/,imode/0/
c
      dkoef1=475.d0/502.d0
      dkoef2= 27.d0/502.d0
c
c---> c in rydberg units:
      c=274.072d0
      csq=c*c
      csqin=1./csq
      if(.not.rel) csqin=0.d0
      fl=dfloat(l)
      xk=fl*(fl+1.d0)
      rnot=dexp(xnot)
      dx2=0.5d0*dx
c
c---> set up starting values for regular solutions
c
      a(1)=astart
      if(.not.rel) then
         b(1)=a(1)*(-vr(1)/(fl+1)+fl/rnot)
      else
         tzoc=-vr(1)/c
         fact=xk-tzoc*tzoc
         fact=dsqrt(fact+1.d0)-1.d0
         fact=c*fact/tzoc
         b(1)=fact*a(1)
      endif
c
      v=vr(1)
      r=rnot
      remv=e*r-v
      br=csqin*remv+r
      ap(1)=br*b(1)+a(1)
      bp(1)=(xk/br-remv)*a(1)-b(1)
c
c---> start runge-kutta procedure (points 2, ... , 6)
c
      do i=1,5
c
         x=xnot+(i-1)*dx
         qi=b(i)
         pi=a(i)
c
         v=vr(i)
         r=dexp(x)
         remv=e*r-v
         br=csqin*remv+r
         k1=br*qi+pi
         m1=(xk/br-remv)*pi-qi
c
         v=0.5d0*(vr(i)+vr(i+1))
         r=dexp(x+dx2)
         remv=e*r-v
         br=csqin*remv+r
         k2=br*(qi+dx2*m1)+(pi+dx2*k1)
         m2=(xk/br-remv)*(pi+dx2*k1)-(qi+dx2*m1)
c
         k3=br*(qi+dx2*m2)+(pi+dx2*k2)
         m3=(xk/br-remv)*(pi+dx2*k2)-(qi+dx2*m2)
c
         v=vr(i+1)
         r=dexp(x+dx)
         remv=e*r-v
         br=csqin*remv+r
         k4=br*(qi+dx*m3)+(pi+dx*k3)
         m4=(xk/br-remv)*(pi+dx*k3)-(qi+dx*m3)
c
         a(i+1)=pi+dx*(k1+2.d0*k2+2.d0*k3+k4)/6.d0
         b(i+1)=qi+dx*(m1+2.d0*m2+2.d0*m3+m4)/6.d0
c
         ap(i+1)=br*b(i+1)+a(i+1)
         bp(i+1)=(xk/br-remv)*a(i+1)-b(i+1)
c
      end do
c
c---> begin adams procedure (points 7, 8, ... , n)
c
      do i=6,n-1
c
         x=xnot+(i-1)*dx
         v=vr(i+1)
         r=dexp(x+dx)
         remv=e*r-v
         br=csqin*remv+r
c
         asum=646.d0*ap(i)-264.d0*ap(i-1)+106.d0*ap(i-2)-19.d0*ap(i-3)
         bsum=646.d0*bp(i)-264.d0*bp(i-1)+106.d0*bp(i-2)-19.d0*bp(i-3)
c
c  predict for point i+1
c
         p0=a(i)+dx*(251.d0*ap(i-4)-1274.d0*ap(i-3)+2616.d0*ap(i-2)-
     >               2774.d0*ap(i-1)+1901.d0*ap(i))/720.d0
         q0=b(i)+dx*(251.d0*bp(i-4)-1274.d0*bp(i-3)+2616.d0*bp(i-2)-
     >               2774.d0*bp(i-1)+1901.d0*bp(i))/720.d0
c
c  correct
c
         if(imode.eq.1) then
c
           do it=1,itmax
              apip1=br*q0+p0
              bpip1=(xk/br-remv)*p0-q0
              p1=a(i)+dx*(251.d0*apip1+asum)/720.d0
              q1=b(i)+dx*(251.d0*bpip1+bsum)/720.d0
c
c  compare predictor with corrector
c
              if(cdabs(p1-p0).lt.tol*cdabs(p0).and.
     >           cdabs(q1-q0).lt.tol*cdabs(q0))  goto 10
c
              p0=p1
              q0=q1
           end do
c
           write(6,*) i+1,r,' not converged'
c
   10      continue
c
         else
c
           apip1=br*q0+p0
           bpip1=(xk/br-remv)*p0-q0
           p1=a(i)+dx*(251.d0*apip1+asum)/720.d0
           q1=b(i)+dx*(251.d0*bpip1+bsum)/720.d0
           p1=dkoef1*p1+dkoef2*p0
           q1=dkoef1*q1+dkoef2*q0
c
         endif
c
         a(i+1)=p1
         b(i+1)=q1
         ap(i+1)=br*q1+p1
         bp(i+1)=(xk/br-remv)*p1-q1
c
      end do
c
c---> logarithmic derivarive at rn
c
      logder = ( ap(n)/a(n) - 1.d0 ) / rn
c
      return
      end
      subroutine inwrd (rel,e,l,vr,dx,rn,n,a,ap)
c ====================
c     ******************************************************************
c     ***   performs the inward integration.                         ***
c     ***   the adams integration has been employed for solving the  ***
c     ***   spin-orbitless relativistic equations.                   ***
c     ***   see harmon & koelling, j.phys.c 10,2107(1977)            ***
c     ******************************************************************
c
c  input:  rel - logical variable for non-relativistic or
c                scalar-relativistic calculation
c          e   - complex*16 energy variable
c          l   - angular momentum quantumnumber
c          vr  - array for r*potential
c          dx  - logarithmic increment for radial scale
c          rn  - radius at which the integration starts inward
c          n   - index corresponding to rn
c
c  output: a   - array for r*wavefuction
c          ap  - array for derivative of a
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
      parameter (nr=nrad)
c
      logical rel
      dimension vr(nr),xlag(nr),vrlag(nr)
      complex*16 a(nr),b(nr),ap(nr),bp(nr)
      complex*16 e,esqrt
      complex*16 pi,qi,remv,br,k1,m1,k2,m2,k3,m3,k4,m4
      complex*16 p0,q0,p1,q1,apim1,bpim1,asum,bsum
      complex*16 fb(0:l1maxp),fn(0:l1maxp),fh(0:l1maxp)
      data tol/1.0d-10/,itmax/100/,imode/0/,nout/5/
c
      dkoef1=475.d0/502.d0
      dkoef2= 27.d0/502.d0
c
c---> c in rydberg units:
      c=274.072d0
      csq=c*c
      csqin=1./csq
      if(.not.rel) csqin=0.d0
      fl=dfloat(l)
      xk=fl*(fl+1.d0)
      dx2=0.5d0*dx
      xn=dlog(rn)
      do i=1,n
        xlag(i)=xn+(i-n)*dx
        vrlag(i)=vr(i)      
      end do
      vrlag(n+1)=0.d0
      xlag(n+1)=xn+nout*dx
      nlag=n+1
      ilag=3
c
c---> set up the starting values outside the muffin tin
c     corresponding to the boundary condition
c
      re=dreal(e)
      esqrt=cdsqrt(e)
      if(dimag(esqrt).lt.0.d0) esqrt=-esqrt
      x=xn+nout*dx 
      r=dexp(x)
      call csbf(l+1,l1maxp,esqrt,r,fb,fn,fh)
      v=0.d0
      remv=e*r-v
      br=csqin*remv+r
      a(n+nout)= r*fb(l)
      b(n+nout)= fl*fb(l)-esqrt*r*fb(l+1)
      ap(n+nout)=br*b(n+nout)+a(n+nout)
      bp(n+nout)=(xk/br-remv)*a(n+nout)-b(n+nout)
c
c---> start runge-kutta procedure (points n+nout-1, ... , n+nout-5)
c
      do i=n+nout,n+nout-4,-1
c
         x=xn+(i-n)*dx
         qi=b(i)
         pi=a(i)
c
         r=dexp(x)
         rmid=dexp(x-dx2)
         rm1=dexp(x-dx)
         v=ylag(x,xlag,vrlag,0,ilag,nlag,iex)
         vmid=ylag(x-dx2,xlag,vrlag,0,ilag,nlag,iex)
         vm1=ylag(x-dx,xlag,vrlag,0,ilag,nlag,iex)
c
         remv=e*r-v
         br=csqin*remv+r
         k1=br*qi+pi
         m1=(xk/br-remv)*pi-qi
c
         remv=e*rmid-vmid
         br=csqin*remv+rmid
         k2=br*(qi-dx2*m1)+(pi-dx2*k1)
         m2=(xk/br-remv)*(pi-dx2*k1)-(qi-dx2*m1)
c
         k3=br*(qi-dx2*m2)+(pi-dx2*k2)
         m3=(xk/br-remv)*(pi-dx2*k2)-(qi-dx2*m2)
c
         remv=e*rm1-vm1
         br=csqin*remv+rm1
         k4=br*(qi-dx*m3)+(pi-dx*k3)
         m4=(xk/br-remv)*(pi-dx*k3)-(qi-dx*m3)
c
         a(i-1)=pi-dx*(k1+2.d0*k2+2.d0*k3+k4)/6.d0
         b(i-1)=qi-dx*(m1+2.d0*m2+2.d0*m3+m4)/6.d0
         ap(i-1)=br*b(i-1)+a(i-1)
         bp(i-1)=(xk/br-remv)*a(i-1)-b(i-1)
c
      end do
c
c---> begin adams procedure (points n+nout-6, ... , 1)
c
      do i=n+nout-5,2,-1
c
         x=xn+(i-n)*dx
         r=dexp(x-dx)
         if(i.gt.n+1) then
           v=ylag(x-dx,xlag,vrlag,0,ilag,nlag,iex)
         else
           v=vr(i-1)
         endif
         remv=e*r-v
         br=csqin*remv+r
c
         asum=646.d0*ap(i)-264.d0*ap(i+1)+106.d0*ap(i+2)-19.d0*ap(i+3)
         bsum=646.d0*bp(i)-264.d0*bp(i+1)+106.d0*bp(i+2)-19.d0*bp(i+3)
c
c  predict for point i-1
c
         p0=a(i)-dx*(251.d0*ap(i+4)-1274.d0*ap(i+3)+2616.d0*ap(i+2)-
     >               2774.d0*ap(i+1)+1901.d0*ap(i))/720.d0
         q0=b(i)-dx*(251.d0*bp(i+4)-1274.d0*bp(i+3)+2616.d0*bp(i+2)-
     >               2774.d0*bp(i+1)+1901.d0*bp(i))/720.d0
c
c  correct
c
         if(imode.eq.1) then
c
           do it=1,itmax
              apim1=br*q0+p0
              bpim1=(xk/br-remv)*p0-q0
              p1=a(i)-dx*(251.d0*apim1+asum)/720.d0
              q1=b(i)-dx*(251.d0*bpim1+bsum)/720.d0
c
c  compare predictor with corrector
c
              if(cdabs(p1-p0).lt.tol*cdabs(p0).and.
     >           cdabs(q1-q0).lt.tol*cdabs(q0))  goto 10
c
              p0=p1
              q0=q1
           end do
           write(6,*) i-1,r,' not converged'
c
         else
c
           apim1=br*q0+p0
           bpim1=(xk/br-remv)*p0-q0
           p1=a(i)-dx*(251.d0*apim1+asum)/720.d0
           q1=b(i)-dx*(251.d0*bpim1+bsum)/720.d0
           p1=dkoef1*p1+dkoef2*p0
           q1=dkoef1*q1+dkoef2*q0
c
         endif
c
   10    continue
         a(i-1)=p1
         b(i-1)=q1
         ap(i-1)=br*q1+p1
         bp(i-1)=(xk/br-remv)*p1-q1    
c
      end do
c
      return
      end
