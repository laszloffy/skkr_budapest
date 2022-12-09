c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      function blm(il1,im1,il2,im2,il3,im3)
      implicit real*8(a-h,o-z)
      parameter (pi=3.1415926535897932384626d0) 
      blm=0.0d0
      onm=-1.0d0
      l1=il1
      l2=il2
      l3=il3
      m1=iabs(im1)
      m2=iabs(im2)
      m3=iabs(im3)
2     m=iabs(im1+im2+im3)+mod(l1+l2+l3,2)
      if(m.ne.0.or.l1.gt.(l2+l3).or.l1.lt.iabs(l2-l3))return
      l=l1
      m=m1
      m1=max0(m1,m2,m3)
      if(m1.eq.m)goto 5
      if(m1.eq.m2)goto 3
      if(m1.eq.m3)goto 4
3     l1=l2
      l2=l
      m2=m
      goto 5
4     l1=l3
      l3=l
      m3=m
5     if(l2.ge.l3)goto 6
      l=l2
      m=m2
      l2=l3
      m2=m3
      l3=l
      m3=m
6     is=(l1+l2+l3)/2
      blm=(onm)**(is-l2-m3+m1)*sqrt(dfloat((2*l1+1)*(2*l2+1)*(2*l3+1)) 
     & /pi)/dfloat(2*(2*is+1))*sqrt(fac(l1+m1,m1+m1))/fac(is-l1,is-l1)
     & *sqrt(fac(l2+m2,m2+m2))/fac(is-l2,is-l2)*sqrt(fac(l3+m3,m3+m3)) 
     & *fac(max0(l2+l3-m1,l2-l3+m1),iabs(2*(l3-m1)))**isign(1,l3-m1)
      if(l3.eq.0)return
      do 7 i=1,l3
7     blm=blm/dfloat(2*(2*(is-i)+1))
      a=1.0d0
      aa=1.0d0
      i1=l1+m1
      i2=max0(1,l2+l3-m1)
      i3=0
      i4=l1-m1
      i5=l2-l3+m1
      i6=l3-m3
      it=min0(i4,i6)
      if(it.eq.0)return
      do 8 i=1,it
      i1=i1+1
      i3=i3+1
      i5=i5+1
      aa=-aa*dfloat(i1*i4*i6)/dfloat(i2*i3*i5)
      a=a+aa
      i2=i2-1
      i4=i4-1
8     i6=i6-1
      blm=a*blm
      return
      end
