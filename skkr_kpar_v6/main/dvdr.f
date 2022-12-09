c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine dvdr(n,r,v,dv)
c     ===============
c
      implicit real*8 (a-h,o-z)
      include '../param.h'       
c
c     calculate derivative of the potential
c     dv/dr=(d(r*v)/dr-v)/r
c
c     j.redinger april 1985
c
c
      dimension r(nrad),v(nrad),rv(nrad),dv(nrad)
c
      do i=1,n
        rv(i)=v(i)*r(i)
      end do
c
      call derspl(n,r,rv,dv)
c
      do i=1,n
        dv(i)=(dv(i)-v(i))/r(i)
      end do       
c
      do i=n+1,nrad
        dv(i)=0.d0
      end do
c
      return
      end
      subroutine derspl(n,x,f,d)
c     =================
c
      implicit real*8 (a-h,o-z)
      include '../param.h'       
      parameter (na=nrad*3)
c
      dimension x(n),f(n),d(n),a(na)
c
c     f(i) are the function values at the points x(i) for i=1,n
c     and the spline derivatives d(i) are found.
c     the dimension of a must not be less than 3*n.
c
      do 5 i=2,n
      if(x(i)-x(i-1))1,1,5
1     write(6,3)i
3     format(' return from derspl  ',i3,' out of order')
      a(1)=1.d0
      return
5     continue
      do 30 i=1,n
      j=2
      if(i-1)6,10,6
6     j=n-1
      if(i.eq.n)go to 10
      h1=1.d0/(x(i)-x(i-1))
      h2=1.d0/(x(i+1)-x(i))
      a(3*i-2)=h1
      a(3*i-1)=2.d0*(h1+h2)
      a(3*i)=h2
      d(i)=3*(f(i+1)*h2*h2+f(i)*(h1*h1-h2*h2)-f(i-1)*h1*h1)
      go to 30
10    h1=1.d0/(x(j)-x(j-1))
      h2=1.d0/(x(j+1)-x(j))
      a(3*i-2)=h1*h1
      a(3*i-1)=h1*h1-h2*h2
      a(3*i)=-h2*h2
      d(i)=2.d0*(f(j)*(h2*h2*h2+h1*h1*h1)-f(j+1)*h2*h2*h2-
     -     f(j-1)*h1*h1*h1)
30    continue
      p=a(4)/a(1)
      a(5)=a(5)-p*a(2)
      a(6)=a(6)-p*a(3)
      d(2)=d(2)-p*d(1)
      do 50 i=3,n
      k=3*i-4
      p=a(k+2)/a(k)
      a(k+3)=a(k+3)-p*a(k+1)
      d(i)=d(i)-p*d(i-1)
      if(i.ne.n-1)go to 50
      p=a(k+5)/a(k)
      a(k+5)=a(k+6)-p*a(k+1)
      a(k+6)=a(k+7)
      d(n)=d(n)-p*d(n-2)
50    continue
      d(n)=d(n)/a(3*n-1)
      do 60 i=3,n
      j=n+2-i
60    d(j)=(d(j)-a(3*j)*d(j+1))/a(3*j-1)
      d(1)=(d(1)-d(2)*a(2)-d(3)*a(3))/a(1)
      a(1)=0.d0
      return
      end
