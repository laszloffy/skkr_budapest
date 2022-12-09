c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine denblloyd(denbll,we,lmax,
     >                     tau,tinv,dmat,dmatp,ddmat,ddmatp)
c=======================
c
c Derivative of the band energy of an impurity via Lloyd formula 
c
c input:  we           ... weighting factor for energy integration 
c         lmax         ... actual max. value of angular momentum index l
c         tau          ... unscreened tau matrix 
c         tinv         ... inverse of unscreened impurity t-matrix 
c                          (tau and tinv are calculated in the global frame!!!)
c         dmat         ... rotation matrix
c         dmatp        ... its adjungate 
c         ddmat        ... derivative of the rotation matrix
c         ddmatp       ... its adjungate
c output: denbll       ... integrand for the derivative of the band energy
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      complex*16 we
      complex*16 tau(kmymaxp,kmymaxp)
      complex*16 tinv(kmymaxp,kmymaxp)
      complex*16 dtinv(kmymaxp,kmymaxp)
      complex*16 dmat(kmymaxp,kmymaxp)
      complex*16 dmatp(kmymaxp,kmymaxp)
      complex*16 ddmat(kmymaxp,kmymaxp)
      complex*16 ddmatp(kmymaxp,kmymaxp)
      complex*16 xmat(kmymaxp,kmymaxp)
      complex*16 ymat(kmymaxp,kmymaxp)
c
      data tol/1.0d-10/
c
      pi=4.0d0*datan(1.0d0)
      lmmax=(lmax+1)**2
      kmymax=2*lmmax
c
      call tripmt1(ddmat,dmatp,tinv,xmat,kmymax,kmymax,kmymaxp)
      call tripmt1(tinv,dmat,ddmatp,ymat,kmymax,kmymax,kmymaxp)
      call addmat1(xmat,ymat,dtinv,kmymax,kmymaxp)
      call doubmt1(dtinv,tau,xmat,kmymax,kmymaxp)
      denbll=0.0d0
      do kmy=1,kmymax
          denbll=denbll+dimag(we*xmat(kmy,kmy))/pi
      end do
c     write(6,'(d15.6)') denbll
c
      return
      end
