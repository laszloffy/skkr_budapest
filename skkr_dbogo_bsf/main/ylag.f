c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
c
      function ylag(xi,x,y,ind1,n1,imax,iex)
      implicit real*8 (a-h,o-z)
c
c     program authors a.a.brooks and e.c.long,
c     computing technology center, union carbide corp.,
c     nuclear div., oak ridge, tenn.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     lagrangian interpolation
c     xi is interpolated entry into x-array
c     n is the order of lagrangian interpolation
c     y is array from which ylag is obtained by interpolation
c     ind is the min-i for x(i).gt.xi ==> if ind=0, x-array is searched
c     imax is max index of x-and y-arrays
c     extrapolation can occur ==> iex=-1 or +1
c
c
      dimension x(imax),y(imax)
      ind=ind1
      n=n1
      iex=0
      if (n.le.imax) goto 10
      n=imax
      iex=n
   10 if (ind.gt.0) goto 40
      do 20 j=1,imax
      if (xi-x(j)) 30,130,20
   20 continue
      iex=1
      goto 70
   30 ind=j
   40 if(ind.gt.1) goto 50
      iex=-1
   50 inl=ind-(n+1)/2
      if (inl.gt.0) goto 60
      inl=1
   60 inu=inl+n-1
      if(inu.le.imax) goto 80
   70 inl=imax-n+1
      inu=imax
   80 s=0.d0
      p=1.d0
      do 110 j=inl,inu
      p=p*(xi-x(j))
      d=1.d0
      do 100 i=inl,inu
      if (i.ne.j) goto 90
      xd=xi
      goto 100
   90 xd=x(j)
  100 d=d*(xd-x(i))
  110 s=s+y(j)/d
      ylag=s*p
  120 return
  130 ylag=y(j)
      goto 120
      end
