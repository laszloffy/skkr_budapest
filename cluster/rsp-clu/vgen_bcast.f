c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine vgen_bcast(vra,bra,
     > v0,
c     > dv0,
     > ferr1,entotifc,
     > vmadih,vmadich,vmadid,vmadiq) 
c          
        implicit real*8 (a-h,o-z)
        include '../param.h'
        include 'mpif.h'
          !
        integer root, myrank, nprocs, ierror
        common/mpi0/root,myrank,nprocs
        real*8     vra(nrad,mimp)
        real*8     bra(nrad,mimp)
        real*8     v0
c        real*8     dv0,
        real*8     ferr1, entotifc
        real*8     vmadih(mimp)
        real*8     vmadich(mimp)
        real*8     vmadid(mimp)
        real*8     vmadiq(mimp)
c        call mpi_bcast(efermi,1,mpi_real8,root,mpi_comm_world,ierror)
c
        call mpi_bcast(v0,1,mpi_real8,root,mpi_comm_world,ierror)
c        call mpi_bcast(dv0,1,mpi_real8,root,mpi_comm_world,ierror)
        call mpi_bcast(ferr1,1,mpi_real8,root,mpi_comm_world,ierror)
        call mpi_bcast(entotifc,1,mpi_real8,root,mpi_comm_world,ierror)
c
c
        call mpi_bcast(vra,nrad*mimp,mpi_real8,
     >    root,mpi_comm_world,ierror)
        call mpi_barrier(mpi_comm_world,ierror)
        call mpi_bcast(bra,nrad*mimp,mpi_real8,
     >    root,mpi_comm_world,ierror)
        call mpi_barrier(mpi_comm_world,ierror)
        call mpi_bcast(vmadih,mimp,mpi_real8,root,mpi_comm_world,ierror)
        call mpi_barrier(mpi_comm_world,ierror)
        call mpi_bcast(vmadich,mimp,mpi_real8,
     >    root,mpi_comm_world,ierror)
        call mpi_barrier(mpi_comm_world,ierror)
        call mpi_bcast(vmadid,mimp,mpi_real8,root,mpi_comm_world,ierror)
        call mpi_barrier(mpi_comm_world,ierror)
        call mpi_bcast(vmadiq,mimp,mpi_real8,root,mpi_comm_world,ierror)
        call mpi_barrier(mpi_comm_world,ierror)
c
       end subroutine vgen_bcast
