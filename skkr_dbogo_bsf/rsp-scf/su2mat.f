c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine SU2mat(rb,Smat,Smatinv)
!--------------------------------------------------------------------
! compute SU(2) matrix and its inverse
! S transforms rb to the z direction,
! Smat transforms the spinor eigenvectors of (sigma*rb) to eigenvectors of sigma_z
! Smatinv == Smat^(-1) == Smat^+
!
!

      implicit none
      include '../param.h'
! -------------------------------------------------------------------
      real*8,     intent(in)  :: rb(3)
      complex*16, intent(out) :: Smat(2,2), Smatinv(2,2)
! -------------------------------------------------------------------
      real*8  :: theta, phi, pi, nvec(3), cos_thhalf, sin_thhalf
      real*8,     parameter :: tol=1.0d-12
      complex*16, parameter :: im1=(0.d0,1.d0)

      pi = 4.0d0*atan(1.0d0)

      Smat = (0.d0,0.d0)
      Smatinv = (0.d0,0.d0)
! test for local mode or z orientation, return unit matrices if no rotation needed
      if (localmode .or. abs(rb(3)-1.d0)<tol) then
         ! no need to rotate
         Smat(1,1) = (1.d0,0.d0)
         Smat(2,2) = (1.d0,0.d0)
         Smatinv(1,1) = (1.d0,0.d0)
         Smatinv(2,2) = (1.d0,0.d0)

         ! DEBUG, TODO: REMOVE!!!!!
         theta = 0.d0
         phi = 0.d0
         nvec = (/1.d0,1.d0,1.d0/)
      else
         ! determine axis of rotation and angle
         ! this is not really hard: we need to rotate to (theta=0) from (th,ph)
         ! so the angle is simply -theta of rb, and the axis is in plane
         ! so the axis has (theta=pi/2,phi=ph+pi/2)
         call findangles(rb(1),rb(2),rb(3),theta,phi)
         nvec(1) = cos(phi+pi/2.d0)
         nvec(2) = sin(phi+pi/2.d0)
         nvec(3) = 0.d0  ! not actually used, left for clarity
         ! theta is used for constructing the rotation matrix elements
         cos_thhalf = cos(-theta/2.d0) ! good sign
         sin_thhalf = sin(-theta/2.d0) ! good sign
         ! DEBUG TRY: change sign of theta, inverse transformation...
         !cos_thhalf = cos(theta/2.d0) ! changed sign
         !sin_thhalf = sin(theta/2.d0) ! changed sign
         
         ! rotation matrix from rb to z
         ! Smat = exp(-i sigma nvec alpha/2) = cos(alpha/2) - i sigma nvec sin(alpha/2)
         ! also make use of the fact that nvec(3)==0
         ! note that S is the transformation that maps rb to z --> rotate with -theta actually, see above
         ! and note that we're working in the (down,up) basis, so the 1<->2 swapped of the Pauli matrices is needed
         Smat(2,2) = cos_thhalf                           ! with Pauli (up,up)
         Smat(1,2) = (-im1*nvec(1) + nvec(2))*sin_thhalf  ! with Pauli (down.up)
         Smat(2,1) = (-im1*nvec(1) - nvec(2))*sin_thhalf  ! with Pauli (up,down)
         Smat(1,1) = cos_thhalf                           ! with Pauli (down,down)

         ! the adjoint/inverse of Smat, making use of that nvec(3)==0 --> diagonal of Smat is real
         ! and off-diagonal is proportional to sin(-theta/2)
         Smatinv(1,1) = Smat(1,1)
         Smatinv(2,1) = -Smat(2,1)
         Smatinv(1,2) = -Smat(1,2)
         Smatinv(2,2) = Smat(2,2)
      end if

c      !DEBUG, TODO: REMOVE!!!!
c      write(401,*) 'exchange field direction:'
c      write(401,*) rb
c      write(401,*) 'theta,phi in degrees:'
c      write(401,*) theta/pi*180,phi/pi*180
c      write(401,*) 'rotation axis:'
c      write(401,*) nvec
c      write(401,*) 'rotation matrix Smat:'
c      write(401,*) Smat(1,:)
c      write(401,*) Smat(2,:)
c      write(401,*) 'inverse rotation matrix Smatinv:'
c      write(401,*) Smatinv(1,:)
c      write(401,*) Smatinv(2,:)
c
      end subroutine SU2mat
