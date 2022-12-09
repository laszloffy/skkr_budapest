c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
          subroutine delta_bcast(ferr1,rlambda,delta,qveha,qvhea,
     >                           entotifc)
c
      implicit real*8 (a-h,o-z)
      include '../param.h'
#ifdef MPIP
      include 'mpif.h'
#endif
c
      integer root, myrank, nprocs, ierror
      common/mpi/root,myrank,nprocs
c
      real*8    ferr1, entotifc
      dimension rlambda(mintfc)
      complex*16 delta(nrad,mintfc)
      complex*16 qveha(mintfc),qvhea(mintfc)
c
#ifdef MPIP
      call mpi_bcast(ferr1,1,mpi_real8,root,mpi_comm_world,ierror)
      call mpi_bcast(entotifc,1,mpi_real8,root,mpi_comm_world,ierror)
      call mpi_bcast(rlambda,mintfc,mpi_real8,root,
     >               mpi_comm_world,ierror)
      call mpi_bcast(delta,nrad*mintfc,mpi_double_complex,
     >    root,mpi_comm_world,ierror)
      call mpi_bcast(qveha,mintfc,mpi_double_complex,
     >    root,mpi_comm_world,ierror)
      call mpi_bcast(qvhea,mintfc,mpi_double_complex,
     >    root,mpi_comm_world,ierror)
#endif
c
      end subroutine delta_bcast
 
