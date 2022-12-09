c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
          subroutine vgen_bcast(efermi,defermi,vra,vrb,bra,brb,
     > bopra,boprb,v0,dv0,rerr1,
     > enpota,enela,enxca,enpotb,enelb,enxcb,enmaga,enmagb,vmadih) 
c          
        implicit none
        include '../param.h'
        include 'mpif.h'
          !
        integer root, myrank, nprocs, ierror
        common/mpi/root,myrank,nprocs
        real*8   efermi
        real*8   defermi
        real*8   vra(nrad,mintfc),vrb(nrad,mintfc)
        real*8   bra(nrad,mintfc), brb(nrad,mintfc)
        real*8   bopra(nrad,2,mintfc)
        real*8   boprb(nrad,2,mintfc)
        real*8   v0, dv0, rerr1
        real*8   vmadih(mimp)
        real*8   enpota(mintfc),enpotb(mintfc)
        real*8   enela(mintfc),enelb(mintfc)
        real*8   enxca(mintfc),enxcb(mintfc)
        real*8   enmaga(mintfc),enmagb(mintfc)
c
        call mpi_bcast(efermi,1,mpi_real8,root,mpi_comm_world,ierror)
        call mpi_bcast(defermi,1,mpi_real8,root,mpi_comm_world,ierror)
        call mpi_bcast(v0,1,mpi_real8,root,mpi_comm_world,ierror)
        call mpi_bcast(dv0,1,mpi_real8,root,mpi_comm_world,ierror)
        call mpi_bcast(rerr1,1,mpi_real8,root,mpi_comm_world,ierror)
        call mpi_bcast(vmadih,mimp,mpi_real8,root,mpi_comm_world,ierror)
        call mpi_bcast(vra,nrad*mintfc,mpi_real8,
     >    root,mpi_comm_world,ierror)
        call mpi_bcast(vrb,nrad*mintfc,mpi_real8,
     >    root,mpi_comm_world,ierror)
        call mpi_bcast(bra,nrad*mintfc,mpi_real8,
     >    root,mpi_comm_world,ierror)
        call mpi_bcast(brb,nrad*mintfc,mpi_real8,
     >    root,mpi_comm_world,ierror)
      call mpi_bcast(enpota,mintfc,mpi_real8,
     >    root,mpi_comm_world,ierror)
      call mpi_bcast(enpotb,mintfc,mpi_real8,
     >    root,mpi_comm_world,ierror)
      call mpi_bcast(enela,mintfc,mpi_real8,
     >    root,mpi_comm_world,ierror)
      call mpi_bcast(enelb,mintfc,mpi_real8,
     >    root,mpi_comm_world,ierror)
      call mpi_bcast(enxca,mintfc,mpi_real8,
     >    root,mpi_comm_world,ierror)
      call mpi_bcast(enxcb,mintfc,mpi_real8,
     >    root,mpi_comm_world,ierror)
      call mpi_bcast(enmaga,mintfc,mpi_real8,
     >    root,mpi_comm_world,ierror)
      call mpi_bcast(enmagb,mintfc,mpi_real8,
     >    root,mpi_comm_world,ierror)
        call mpi_bcast(bopra,2*mintfc*nrad, mpi_real8,
     >    root,mpi_comm_world,ierror)
        call mpi_bcast(boprb,2*mintfc*nrad, mpi_real8,
     >    root,mpi_comm_world,ierror)

       end subroutine vgen_bcast
