c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine gjinv(a,n,nmax,detl)
!=====================
!
      implicit none
!-------------------------------------------------
      integer,    intent(in)    :: n, nmax
      complex*16, intent(inout) :: a(nmax,nmax)
      complex*16, intent(out)   :: detl
!-------------------------------------------------
      complex*16, parameter :: ci = (0.d0,1.d0)
      complex*16, parameter :: ipi=(0.d0,3.1415926535897932d0)
      complex*16          :: atemp(n,n)
      integer             :: ipiv(n), info, lwork, i, rowswap
      complex*16              :: dlwork   !dummy double complex for lwork in query call to zgetri
      complex*16, allocatable :: work(:)

      atemp = a(1:n,1:n)

!     calculate pivoted LU decomposition, overwrite matrix
      call zgetrf(n,n,atemp,n,ipiv,info) !actually A=P*L*U, atemp=L\U
      if (info /= 0) then
!      if (info < 0) then !TODO check and remove
        write(6,'("Error in gjinv->dgetrf: info=",i4)') info
        stop
      end if

      detl = (0.d0,0.d0)  !=ln det A
      !ln det U
      do i=1,n
        !det A = det U = prod(diag(atemp))
        detl = detl + log(abs(atemp(i,i))) +
     >             ci*atan2(aimag(atemp(i,i)),real(atemp(i,i)))
      end do

      rowswap = 0
      !ln det P, TODO: check if this is OK
      do i=1,n
        if (ipiv(i) /= i) then  !rows where changed
          rowswap = rowswap + 1
        end if
      end do
      if (mod(rowswap,2) == 1) then !det P = -1
        detl = detl - ipi  !ln det P = ln (-1) = i pi, sign according to gjinv.f
      end if

!     get optimal lwork value
      call zgetri(n,atemp,n,ipiv,dlwork,-1,info)  !lwork has optimal lwork value
      lwork = dble(dlwork)  !implicit double -> int

!     allocate work array
      allocate(work(lwork))

!     calculate inverse from LU, overwrite matrix
      call zgetri(n,atemp,n,ipiv,work,lwork,info)
      if (info /= 0) then
        write(6,'("Error in gjinv->dgetri: info=",i4)') info
        stop
      end if

!     deallocate work array
      deallocate(work)

      a(1:n,1:n) = atemp

      end subroutine gjinv
