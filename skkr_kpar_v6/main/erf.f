c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      function erf(a)
c=================
c
      implicit real*8 (a-h,o-z)
      parameter(tol=1.d-14)
c
      pi=4.d0*datan(1.d0)
      if(a.lt.1.0e-3) then
        erf=1.d0-a**2/3.d0+a**4/10.d0-a**6/42.d0
        erf=2.d0*a*erf/dsqrt(pi)
        return
      end if
      ee=dexp(a*a)
      erf=1.d0
      term=1.d0
        do n=1,1000
          term=2.d0*term*a*a/(2*n+1)
          erf=erf+term
          if(dabs(term).lt.tol) goto 2
        end do
      stop'erf: non converging'
    2 erf=2.d0*a*erf/ee/dsqrt(pi)
c
      return
      end
