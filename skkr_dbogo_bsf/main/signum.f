c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      real*8 function signum(x)
      real*8 x,tol
      data tol/1.0d-14/
      if(dabs(x).lt.tol) signum=0.0d0
      if(x.gt.tol) signum=1.0d0
      if(x.lt.-tol) signum=-1.0d0
      return
      end
       
