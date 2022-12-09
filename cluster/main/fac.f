c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      function fac(n,m)
      implicit real*8(a-h,o-z)
      fac=1.0d+00
      if(m.eq.0)return
      m1=n-m+1
      do 1 i=m1,n
1     fac=fac*float(i)
      return
      end
