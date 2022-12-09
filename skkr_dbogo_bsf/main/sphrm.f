c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine sphrm (ylm,nn,ct,st,cf,lmax)
c     ================
c
c----- calculate spherical harmonics ylm(theta,fi) for
c      complex arguments.  see pendry appendix a.  the ylm are ordered
c      (lm) = (00),(1-1),(10),(11),(2-2),(2-1),.......
c
c      ct=cos(theta), st=sin(theta), cf=exp(i*fi)
c      lmax  maximum order of l
c ---------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex *16 ylm(nn),ct,st,cf,sf,sa,fac
c
c----- set ylm(00)
c
      pi=4.0d0*atan(1.0d0)
      pii4=0.25d0/pi
      ylm(1)=sqrt(pii4)
      if (lmax.eq.0) return
c
c----- set ylm (m=l,m=l-1) using explicit expressions (a.16) and (a.17)
c
      a=1.0
      b=1.0
      asg=1.0
      sf=1.0
      sa=1.0
      lp=1
      do 10 l=1,lmax
      fl=dfloat(l)
      a=0.5*a*fl*(2.0*fl-1.0)
      b=fl*b
      asg=-asg
      lm=lp+1
      lp=lp+l+l+1
      sf=sf*cf
      fac=sqrt((2.0*fl+1.0)*a/(4.0*pi*b*b))*sa
      ylm(lm)=fac*st/sf
      ylm(lp)=asg*fac*st*sf
      fac=sqrt(2.0*fl)*fac*ct
      ylm(lm+1)=fac*cf/sf
      ylm(lp-1)=-asg*fac*sf/cf
   10 sa=sa*st
c
      if (lmax.eq.1) return
c
c----- set remaining ylm using recurrence relation in l (a.14)
c
      do 20 m=2,lmax
      mm=m+m-4
      fm=dfloat(m-2)
      a=sqrt(1.0/(fm+fm+3.0))
      ln=m*m-1
      lm=ln-m-m+2
      do 20 l=m,lmax
      fl=dfloat(l)
      b=sqrt((fl+fm)*(fl-fm)/((fl+fl+1.0)*(fl+fl-1.0)))
      lp=ln+l+l
      ylm(lp)=(ct*ylm(ln)-a*ylm(lm))/b
      ylm(lp-mm)=(ct*ylm(ln-mm)-a*ylm(lm-mm))/b
      a=b
      lm=ln
   20 ln=lp
c
      return
      end

