c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine jij_bcast(vij,gradph,gradth)
c          
        implicit real*8 (a-h,o-z)
        include '../param.h'
        include 'mpif.h'
          !
        integer root, myrank, nprocs, ierror
        common/mpi0/root,myrank,nprocs
        real*8 vij(4,mimp,mimp)
        real*8 gradth(mimp)
        real*8 gradph(mimp)      
c
        call mpi_bcast(vij,4*mimp*mimp,mpi_real8,
     >          root,mpi_comm_world,ierror)
        call mpi_bcast(gradph,mimp,mpi_real8,root,mpi_comm_world,ierror)
        call mpi_bcast(gradth,mimp,mpi_real8,root,mpi_comm_world,ierror)

       end subroutine jij_bcast
