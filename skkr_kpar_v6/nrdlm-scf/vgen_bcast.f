c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
        subroutine vgen_bcast(efermi,defermi,vra,vrb,bra,brb,rs,rsl,rsr,
     >                        v0,dv0,ferr1,ferr2,entotifc) 
c          
        implicit real*8 (a-h,o-z)
        include '../param.h'
        include 'mpif.h'
          !
        integer root, myrank, nprocs, ierror
        common/mpi/root,myrank,nprocs
        real*8     efermi
        real*8     defermi
        real*8     vra(nrad,mintfc),vrb(nrad,mintfc)
        real*8     bra(nrad,mintfc), brb(nrad,mintfc)
        real*8     rs(mintfc),rsl(minprc),rsr(minprc)
        real*8     v0
        real*8     dv0, ferr1, ferr2, entotifc
c
        call mpi_bcast(efermi,1,mpi_real8,root,mpi_comm_world,ierror)
        call mpi_bcast(defermi,1,mpi_real8,root,mpi_comm_world,ierror)
        call mpi_bcast(v0,1,mpi_real8,root,mpi_comm_world,ierror)
        call mpi_bcast(dv0,1,mpi_real8,root,mpi_comm_world,ierror)
        call mpi_bcast(ferr1,1,mpi_real8,root,mpi_comm_world,ierror)
        call mpi_bcast(ferr2,1,mpi_real8,root,mpi_comm_world,ierror)
        call mpi_bcast(entotifc,1,mpi_real8,root,mpi_comm_world,ierror)
        call mpi_bcast(rs,mintfc,mpi_real8,root,mpi_comm_world,ierror)
        call mpi_bcast(rsl,minprc,mpi_real8,root,mpi_comm_world,ierror)
        call mpi_bcast(rsr,minprc,mpi_real8,root,mpi_comm_world,ierror)
        call mpi_bcast(vra,nrad*mintfc,mpi_real8,
     >    root,mpi_comm_world,ierror)
        call mpi_bcast(vrb,nrad*mintfc,mpi_real8,
     >    root,mpi_comm_world,ierror)
        call mpi_bcast(bra,nrad*mintfc,mpi_real8,
     >    root,mpi_comm_world,ierror)
        call mpi_bcast(brb,nrad*mintfc,mpi_real8,
     >    root,mpi_comm_world,ierror)

       end subroutine vgen_bcast
