c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine denblloyd(we,lmax,tau,tinv,dmat,dmatp,
     >ddph,ddphp,ddth,ddthp,d2dph,d2dphp,d2dth,d2dthp,d2dthph,d2dthphp,
     >ddphenb,ddthenb,d2dphenb,d2dthenb,d2dthphenb)
c=======================
c
c Derivative of the band energy of an impurity via Lloyd formula 
c
c input:  we           ... weighting factor for energy integration 
c         lmax         ... actual max. value of angular momentum index l
c         tau          ... unscreened tau matrix in global frame
c         tinv         ... inverse of unscreened t-matrix in local frame
c         dmat         ... rotation matrix
c         dmatp        ... its adjungate 
c         ddph,...     ... derivatives of the rotation matrix
c         ddphp,...    ... adjungate matrices
c output: ddphenb,...  ... integrands for the derivatives of the band energy
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      complex*16 we
      complex*16 tau(kmymaxp,kmymaxp)
      complex*16 tinv(kmymaxp,kmymaxp)
      complex*16 dmat(kmymaxp,kmymaxp)
      complex*16 dmatp(kmymaxp,kmymaxp)
      complex*16 ddph(kmymaxp,kmymaxp)
      complex*16 ddth(kmymaxp,kmymaxp)
      complex*16 d2dph(kmymaxp,kmymaxp)
      complex*16 d2dth(kmymaxp,kmymaxp)
      complex*16 d2dthph(kmymaxp,kmymaxp)
      complex*16 ddphp(kmymaxp,kmymaxp)
      complex*16 ddthp(kmymaxp,kmymaxp)
      complex*16 d2dphp(kmymaxp,kmymaxp)
      complex*16 d2dthp(kmymaxp,kmymaxp)
      complex*16 d2dthphp(kmymaxp,kmymaxp)
      complex*16 ddthtinv(kmymaxp,kmymaxp)
      complex*16 ddphtinv(kmymaxp,kmymaxp)
      complex*16 d2dthtinv(kmymaxp,kmymaxp)
      complex*16 d2dphtinv(kmymaxp,kmymaxp)
      complex*16 d2dthphtinv(kmymaxp,kmymaxp)
      complex*16 mat1(kmymaxp,kmymaxp)
      complex*16 mat2(kmymaxp,kmymaxp)
      complex*16 mat3(kmymaxp,kmymaxp)
      complex*16 mat4(kmymaxp,kmymaxp)
c
      data tol/1.0d-10/
c
      pi=4.0d0*datan(1.0d0)
      lmmax=(lmax+1)**2
      kmymax=2*lmmax
c
      call tripmt1(ddth,tinv,dmatp,mat1,kmymax,kmymax,kmymaxp)
      call tripmt1(dmat,tinv,ddthp,mat2,kmymax,kmymax,kmymaxp)
      call addmat1(mat1,mat2,ddthtinv,kmymax,kmymaxp)
c
      call tripmt1(ddph,tinv,dmatp,mat1,kmymax,kmymax,kmymaxp)
      call tripmt1(dmat,tinv,ddphp,mat2,kmymax,kmymax,kmymaxp)
      call addmat1(mat1,mat2,ddphtinv,kmymax,kmymaxp)
c
      call tripmt1(d2dth,tinv,dmatp,mat1,kmymax,kmymax,kmymaxp)
      call tripmt1(dmat,tinv,d2dthp,mat2,kmymax,kmymax,kmymaxp)
      call tripmt1(ddth,tinv,ddthp,mat3,kmymax,kmymax,kmymaxp)
      call addmat1(mat1,mat2,d2dthtinv,kmymax,kmymaxp)
      call addmat(d2dthtinv,mat3,kmymax,kmymaxp)
      call addmat(d2dthtinv,mat3,kmymax,kmymaxp)
c
      call tripmt1(d2dph,tinv,dmatp,mat1,kmymax,kmymax,kmymaxp)
      call tripmt1(dmat,tinv,d2dphp,mat2,kmymax,kmymax,kmymaxp)
      call tripmt1(ddph,tinv,ddphp,mat3,kmymax,kmymax,kmymaxp)
      call addmat1(mat1,mat2,d2dphtinv,kmymax,kmymaxp)
      call addmat(d2dphtinv,mat3,kmymax,kmymaxp)
      call addmat(d2dphtinv,mat3,kmymax,kmymaxp)
c
      call tripmt1(d2dthph,tinv,dmatp,mat1,kmymax,kmymax,kmymaxp)
      call tripmt1(dmat,tinv,d2dthphp,mat2,kmymax,kmymax,kmymaxp)
      call tripmt1(ddth,tinv,ddphp,mat3,kmymax,kmymax,kmymaxp)
      call tripmt1(ddph,tinv,ddthp,mat4,kmymax,kmymax,kmymaxp)
      call addmat1(mat1,mat2,d2dthphtinv,kmymax,kmymaxp)
      call addmat(d2dthphtinv,mat3,kmymax,kmymaxp)
      call addmat(d2dthphtinv,mat4,kmymax,kmymaxp)
c
      call doubmt1(tau,ddthtinv,mat1,kmymax,kmymaxp)
      ddthenb=0.0d0
      do kmy=1,kmymax
          ddthenb=ddthenb+dimag(we*mat1(kmy,kmy))/pi
      end do
c
      call doubmt1(tau,ddphtinv,mat1,kmymax,kmymaxp)
      ddphenb=0
      do kmy=1,kmymax
          ddphenb=ddphenb+dimag(we*mat1(kmy,kmy))/pi
      end do
c
      call doubmt1(tau,d2dthtinv,mat1,kmymax,kmymaxp)
      call doubmt1(tau,ddthtinv,mat2,kmymax,kmymaxp)
      call doubmt1(mat2,mat2,mat3,kmymax,kmymaxp)
      d2dthenb=0.0d0
      do kmy=1,kmymax
          d2dthenb=d2dthenb+dimag(we*(mat1(kmy,kmy)-mat3(kmy,kmy)))/pi
      end do
c
      call doubmt1(tau,d2dphtinv,mat1,kmymax,kmymaxp)
      call doubmt1(tau,ddphtinv,mat2,kmymax,kmymaxp)
      call doubmt1(mat2,mat2,mat3,kmymax,kmymaxp)
      d2dphenb=0.0d0
      do kmy=1,kmymax
          d2dphenb=d2dphenb+dimag(we*(mat1(kmy,kmy)-mat3(kmy,kmy)))/pi
      end do
c
      call doubmt1(tau,d2dthphtinv,mat1,kmymax,kmymaxp)
      call doubmt1(tau,ddthtinv,mat2,kmymax,kmymaxp)
      call doubmt1(tau,ddphtinv,mat3,kmymax,kmymaxp)
      call doubmt1(mat2,mat3,mat4,kmymax,kmymaxp)
      d2dthphenb=0.0d0
      do kmy=1,kmymax
          d2dthphenb=d2dthphenb+
     >               dimag(we*(mat1(kmy,kmy)-mat4(kmy,kmy)))/pi
      end do
c
      return
      end
