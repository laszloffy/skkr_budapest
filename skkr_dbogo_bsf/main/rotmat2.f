c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine rotmat2 (d,tr,ddthe,d2dthe,ij0,lm,ig)
c rotation matrix for spherical harmonics with angular momentum
c number j
c Phi= angle of rotation (theta)
c input: ij0 = 2j+1 (j may be integer or half-integer)
c        lm ... lm(1)=cos(Phi/2)
c               lm(2..4)=sin(Phi/2)*n(x..z), where n=normalized
c               direction vector
c        ig ... >0 proper rotations, <0 improper rotations
c output: d ... rotation matrix in condon-shortley convention
c         tr ... trace of d
c         ddthe ... derivative of the rotation matrix via Phi
c         d2dthe ... second derivative of the rot. matrix via Phi
c
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
      parameter (jdim=2*lmaxp+2)
c
      complex*16 d(jdim,jdim),f1,f2,f3,f4,fk,c
      complex*16 ddthe(jdim,jdim),d2dthe(jdim,jdim)
      real*8 lm(4),fac(jdim+2),sp,ft,ft2,ffac
c
c     sp=sin(Phi/2)
      sp=dsin(dacos(lm(1)))
      fac(1)=1.d0
      do 10 i=1,jdim+1
   10    fac(i+1)=fac(i)*i
      f1=dcmplx(lm(1),-lm(4))
      f2=dconjg(f1)
      f3=dcmplx(lm(3),-lm(2))
      f4=dconjg(f3)
      ij=iabs(ij0)
      do 100 in=1,ij
         do 100 im=1,ij
         km=max0(1,in-im+1)
         kx=min0(ij-im+1,in)
         d(in,im)=(0.d0,0.d0)
         ddthe(in,im)=(0.d0,0.d0)
         d2dthe(in,im)=(0.d0,0.d0) 
         do 90 k=km,kx
            k0=k-1
            fk=(1.d0,0.d0)
            ffac=fac(ij-im-k0+1)*fac(in-k0)*fac(k)*fac(im-in+k)
c          ft=(-1)**(im-in+k0)*((im-in+2*k0)*(f1**2)-(ij-1-
c     1      (im-in+2*k0))*(sp**2))/ffac
          ft2=(-1)**(im-in+k0)*((im-in+2*k0)*(f1**2)*((im-in+2*k0-1)*(f1
     1      **2)-(ij-(im-in+2*k0))*(sp**2))-(ij-1-(im-in+2*k0))*(sp**2)
     1      *((im-in+2*k0+1)*(f1**2)-(sp**2)*(ij-2-(im-in+2*k0))))/ffac
            if (ij-im-k0.eq.0) go to 20
            fk=f1**(ij-im-k0)
   20       if (in-k.eq.0) go to 30
            fk=fk*f2**(in-k)
   30       if (k0.eq.0) go to 40
            if (lm(2).eq.0.d0) fk=fk*lm(3)**k0
            if (lm(2).ne.0.d0) fk=fk*f3**k0
   40       if (im-in+k0.eq.0) go to 50
            if (lm(2).eq.0.d0) fk=fk*(-lm(3))**(im-in+k0)
            if (lm(2).ne.0.d0) fk=fk*(-f4)**(im-in+k0)
   50       fk=fk/ffac
            d(in,im)=d(in,im)+fk
c            if (im-in+2*k0.eq.1) go to 60
c            ft=ft*(sp**(im-in+2*k0-1))
c   60       if (im-in+2*k0-ij+2.eq.0) go to 70
c            ft=ft*(f1**(ij-2-(im-in+2*k0))) 
c   70       ddthe(in,im)=ddthe(in,im)+ft/2
            ddthe(in,im)=ddthe(in,im)+fk*((im-in+2*k0)*lm(1)/sp-
     1      (ij-1+in-im-2*k0)*sp/lm(1))/2
            if (im-in+2*k0.eq.2) go to 80
            ft2=ft2*(sp**(im-in+2*k0-2))
   80       if (im-in+2*k0-ij+3.eq.0) go to 90
            ft2=ft2*(f1**(ij-3-(im-in+2*k0))) 
   90       d2dthe(in,im)=d2dthe(in,im)+ft2/4
          d(in,im)=d(in,im)*dsqrt(fac(ij-in+1)*fac(in)*fac(ij-im+1)*
     1    fac(im))
          ddthe(in,im)=ddthe(in,im)*dsqrt(fac(ij-in+1)*fac(ij-im+1)*
     1    fac(in)*fac(im))
          d2dthe(in,im)=d2dthe(in,im)*dsqrt(fac(ij-in+1)*fac(ij-im+1)*
     1    fac(in)*fac(im))*f4**(im-in)*sp**(in-im)
         if (mod(ij,4).eq.3.and.isign(1,ig).ne.1) d(in,im)=-d(in,im)
         if (ij0.lt.0.and.isign(1,ig).ne.1) d(in,im)=-d(in,im)
         if (dabs(dreal(d(in,im))).lt.1.d-14) d(in,im)=dcmplx(0.d0,
     1    dimag(d(in,im)))
  100    if (dabs(dimag(d(in,im))).lt.1.d-14) d(in,im)=dreal(d(in,im))
      c=(0.d0,0.d0)
      do 110 i=1,ij
  110    c=c+d(i,i)
      if (dabs(dimag(c)).gt.1.d-10) stop 10
      tr=c
      return
      end
