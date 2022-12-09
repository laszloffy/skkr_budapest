c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine poissonl(den,z,pot,dx,x0,j)
c========================
c
c  NUMERICAL SOLUTION OF POISSON'S EQUATION
c      (TERRY LOUCKS: APW METHOD, PP.98-104)
c
c den : array containing charge density ( times 4*pi*r**2 ) on 
c       equidistant logarithmic mesh
c z   : total charge
c pot : result potential 
c dx  : step on logarithmic scale
c x0  : first point on log. scale
c j   : index referring to a point where den is negligible
c 
      implicit real*8(a-h,o-z)
c
      include '../param.h'
c
      dimension den(nrad),pot(nrad),e(nrad),f(nrad),w(nrad)
c
      a=1.d0-dx*dx/48.d0
      b=-2.d0*(1.d0+5.d0*dx*dx/48.d0)
      edl=dexp(0.5d0*dx)
      c=dx*dx/6.d0
      c2=-b/a
      e(1)=0.d0
      f(1)=edl
c
      x=x0+dx
      itop=j-1
      do i=2,itop
        d=c*dexp(0.5d0*x)*(edl*den(i+1)+10.d0*den(i)+den(i-1)/edl)
        f(i)=c2-1.d0/f(i-1)
        e(i)=(d/a+e(i-1))/f(i)
        x=x+dx
      end do
      w(j)=2.d0*z*dexp(-0.5d0*x)
      do i=1,itop
        jv=j-i
        w(jv)=e(jv)+w(jv+1)/f(jv)
      end do
c
      x=x0
      do jv=1,j
        rh=dexp(0.5d0*x)
        r=rh*rh
        pot(jv)=w(jv)*rh/r
        x=x+dx
      end do
c
      return
      end
      subroutine poissond(den,z,pot,dx,x0,sws,j)
c========================
c
c  NUMERICAL SOLUTION OF POISSON'S EQUATION BY MEANS
c  OF DIRECT INTEGRATION
c
c den : array containing charge density ( times 4*pi*r**2 ) on 
c       equidistant logarithmic mesh
c z   : total charge
c pot : result potential 
c dx  : step on logarithmic scale
c x0  : first point on log. scale
c sws : radius to match
c j   : index just before sws
c 
      implicit real*8(a-h,o-z)
c
      include '../param.h'
      dimension den(nrad),pot(nrad),work(nrad),rsc(nrad)
c
      x=x0
      do i=1,j+1
        rsc(i)=dexp(x)
        work(i)=den(i)/rsc(i)
        x=x+dx
      end do
c
      q1=rsimp(work,rsc,j-1,dx)
      q2=rsimp(work,rsc,j,dx)
      q3=rsimp(work,rsc,j+1,dx)
      w1=rsc(j-1)
      w2=rsc(j)
      w3=rsc(j+1)
      d1=(sws-w2)*(sws-w3)/((w1-w2)*(w1-w3))
      d2=(sws-w1)*(sws-w3)/((w2-w1)*(w2-w3))
      d3=(sws-w1)*(sws-w2)/((w3-w1)*(w3-w2))
      qq=q1*d1+q2*d2+q3*d3
c
      pot(1)=2.d0*qq
      do i=2,j+1
        qq1=rsimp(den,rsc,i,dx)
        qq2=rsimp(work,rsc,i,dx)
        pot(i)=2.d0*(qq1+(qq-qq2)*rsc(i))/rsc(i)
      end do
c
      return
      end
