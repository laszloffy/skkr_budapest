c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine collect_data(kmymax,nimp,ne,
     > dosa,qvpa,qva,qvdiffa,nqva,
     > vmadid,vmadiq,enba,enbdiffa,denba,nenba,enorba,qmoma,
     > rhova,rhospa,rhodspa,rhomaga,
     > dosmaga,spin_magvpa,spin_magva,orb_magvpa,orb_magva,
     > vij,gradph,gradth,
     > ebl,detl2,enpba)
c
      implicit none
      include '../param.h'
      include 'mpif.h'
c
c MPI parameters (esimon)
      integer root, myrank, nprocs, ierror
      common/mpi0/root,myrank,nprocs
c
      integer kmymax,nimp,ne
c
      real*8 dosa(kmymaxp,mimp,me)
      real*8 qvpa(kmymaxp,mimp)
      real*8 qva(mimp)
      complex*16 nqva(mimp,me)
      real*8 qvdiffa(mimp,me)
      real*8 vmadid(mimp)
      real*8 vmadiq(mimp)
      real*8 enba(mimp)
      real*8 enpba(kmymaxp,mimp)
      real*8 enbdiffa(mimp,me)
      real*8 denba(mimp)
      real*8 enorba(mimp)
      complex*16 qmoma(lmsup,mimp)
      complex*16 detl2(me)
      complex*16 nenba(mimp,me)
      real*8 rhova(nrad,mimp)
      real*8 rhospa(nrad,2,mimp)
      real*8 rhodspa(nrad,2,mimp)
      real*8 rhomaga(nrad,mimp)
      real*8 dosmaga(kmymaxp,mimp,me)
      real*8 spin_magvpa(kmymaxp,mimp,3)
      real*8 spin_magva(mimp,3)
      real*8 orb_magvpa(kmymaxp,mimp,3)
      real*8 orb_magva(mimp,3)
      real*8 vij(4,mimp,mimp)
      real*8 gradth(mimp)
      real*8 gradph(mimp)
      real*8 ebl(1)
c
      integer*4 idat,idatc
      integer*4 cidat,cidatc
      integer*4 i1,i2,i3
c
c     parameter (idat=
c    > kmymaxp*mimp*me+kmymaxp*mimp+mimp+mimp*me+mimp+mimp+mimp
c    > +mimp*me+mimp+mimp+nrad*mimp+nrad*2*mimp+nrad*2*mimp
c    > +nrad*mimp+kmymaxp*mimp*me+kmymaxp*mimp*3+mimp*3
c    > +kmymaxp*mimp*3+mimp*3+4*mimp*mimp+mimp+mimp+1)
c
c     real*8 datain(idat),dataout(idat)
c
      real*8, allocatable :: datain(:)
      real*8, allocatable :: dataout(:)
      complex*16, allocatable :: cdatain(:)
      complex*16, allocatable :: cdataout(:)
c
      integer AllocateStatus
c
      idat=kmymax*nimp*ne+kmymax*nimp+nimp+nimp*ne+nimp+nimp+nimp
     > +nimp*ne+nimp+nimp+nrad*nimp+nrad*2*nimp+nrad*2*nimp
     > +nrad*nimp+kmymax*nimp*ne+kmymax*nimp*3+nimp*3
     > +kmymax*nimp*3+nimp*3+4*nimp*nimp+nimp+nimp+1
     > +kmymax*nimp
c
      cidat=lmsup*nimp+ne+nimp*ne+nimp*ne
c
      ALLOCATE ( datain(idat),STAT = AllocateStatus)
      call alloccheck(AllocateStatus,'datain in skkr')
      ALLOCATE ( dataout(idat),STAT = AllocateStatus)
      call alloccheck(AllocateStatus,'dataout in skkr')
      ALLOCATE ( cdatain(cidat),STAT = AllocateStatus)
      call alloccheck(AllocateStatus,'cdatain in skkr')
      ALLOCATE ( cdataout(cidat),STAT = AllocateStatus)
      call alloccheck(AllocateStatus,'cdataout in skkr')
c
c ----------------------------------------------
c
      idatc=0
c
      do i1=1,kmymax
      do i2=1,nimp
      do i3=1,ne
       idatc=idatc+1
       datain(idatc)=dosa(i1,i2,i3)
      end do
      end do
      end do
c
      do i1=1,kmymax
      do i2=1,nimp
       idatc=idatc+1
       datain(idatc)=qvpa(i1,i2)
      end do
      end do
c
      do i1=1,nimp
       idatc=idatc+1
       datain(idatc)=qva(i1)
      end do
c
      do i1=1,nimp
      do i2=1,ne
       idatc=idatc+1
       datain(idatc)=qvdiffa(i1,i2)
      end do
      end do
c
      do i1=1,nimp
       idatc=idatc+1
       datain(idatc)=vmadid(i1)
      end do
c
      do i1=1,nimp
       idatc=idatc+1
       datain(idatc)=vmadiq(i1)
      end do
c
      do i1=1,nimp
       idatc=idatc+1
       datain(idatc)=enba(i1)
      end do
c
      do i1=1,nimp
      do i2=1,ne
       idatc=idatc+1
       datain(idatc)=enbdiffa(i1,i2)
      end do
      end do
c
      do i1=1,nimp
       idatc=idatc+1
       datain(idatc)=denba(i1)
      end do
c
      do i1=1,nimp
       idatc=idatc+1
       datain(idatc)=enorba(i1)
      end do
c
      do i1=1,nrad
      do i2=1,nimp
       idatc=idatc+1
       datain(idatc)=rhova(i1,i2)
      end do
      end do
c
      do i1=1,nrad
      do i2=1,2
      do i3=1,nimp
       idatc=idatc+1
       datain(idatc)=rhospa(i1,i2,i3)
      end do
      end do
      end do
c
      do i1=1,nrad
      do i2=1,2
      do i3=1,nimp
       idatc=idatc+1
       datain(idatc)=rhodspa(i1,i2,i3)
      end do
      end do
      end do
c
      do i1=1,nrad
      do i2=1,nimp
       idatc=idatc+1
       datain(idatc)=rhomaga(i1,i2)
      end do
      end do
c
      do i1=1,kmymax
      do i2=1,nimp
      do i3=1,ne
       idatc=idatc+1
       datain(idatc)=dosmaga(i1,i2,i3)
      end do
      end do
      end do
c
      do i1=1,kmymax
      do i2=1,nimp
      do i3=1,3
       idatc=idatc+1
       datain(idatc)=spin_magvpa(i1,i2,i3)
      end do
      end do
      end do
c
      do i1=1,nimp
      do i2=1,3
       idatc=idatc+1
       datain(idatc)=spin_magva(i1,i2)
      end do
      end do
c
      do i1=1,kmymax
      do i2=1,nimp
      do i3=1,3
       idatc=idatc+1
       datain(idatc)=orb_magvpa(i1,i2,i3)
      end do
      end do
      end do
c
      do i1=1,nimp
      do i2=1,3
       idatc=idatc+1
       datain(idatc)=orb_magva(i1,i2)
      end do
      end do
c
      do i1=1,4
      do i2=1,nimp
      do i3=1,nimp
       idatc=idatc+1
       datain(idatc)=vij(i1,i2,i3)
      end do
      end do
      end do
c
      do i1=1,nimp
       idatc=idatc+1
       datain(idatc)=gradth(i1)
      end do
c
      do i1=1,nimp
       idatc=idatc+1
       datain(idatc)=gradph(i1)
      end do
c
      do i1=1,1
       idatc=idatc+1
       datain(idatc)=ebl(i1)
      end do
c
      do i1=1,kmymax
      do i2=1,nimp
        idatc=idatc+1
        datain(idatc)=enpba(i1,i2)
      end do
      end do
c
c
c     complex*16 qmoma(lmsup,mimp)
c     complex*16 detl2(me)
c     complex*16 nenba(mimp,me)
c
      cidatc=0
      do i1=1,lmsup
      do i2=1,nimp
        cidatc=cidatc+1
        cdatain(cidatc)=qmoma(i1,i2)
      end do
      end do
c
      do i1=1,ne
        cidatc=cidatc+1
        cdatain(cidatc)=detl2(i1)
      end do
c
      do i1=1,nimp
      do i2=1,ne
        cidatc=cidatc+1
        cdatain(cidatc)=nenba(i1,i2)
      end do
      end do
c
      do i1=1,nimp
      do i2=1,ne
        cidatc=cidatc+1
        cdatain(cidatc)=nqva(i1,i2)
      end do
      end do
c
c --------------------------------------
      do i1=1,idat
        dataout(i1) = 0.0
      end do
      do i1=1,cidat
        cdataout(i1) = (0.0,0.0)
      end do
c
       if(myrank.eq.root) write(6,*) 'idat:',idat, idatc
       if(myrank.eq.root) write(6,*) 'cidat:',cidat, cidatc
c
c     call mpi_barrier(mpi_comm_world,ierror)
c
      call mpi_allreduce(datain,dataout,idat,
     > mpi_real8,mpi_sum,mpi_comm_world,ierror)
      call mpi_allreduce(cdatain,cdataout,cidat,
     > mpi_complex16,mpi_sum,mpi_comm_world,ierror)
c
      idatc=0
c
      do i1=1,kmymax
      do i2=1,nimp
      do i3=1,ne
       idatc=idatc+1
       dosa(i1,i2,i3)=dataout(idatc)
      end do
      end do
      end do
c
      do i1=1,kmymax
      do i2=1,nimp
       idatc=idatc+1
       qvpa(i1,i2)=dataout(idatc)
      end do
      end do
c
      do i1=1,nimp
       idatc=idatc+1
       qva(i1)=dataout(idatc)
      end do
c
      do i1=1,nimp
      do i2=1,ne
       idatc=idatc+1
       qvdiffa(i1,i2)=dataout(idatc)
      end do
      end do
c
      do i1=1,nimp
       idatc=idatc+1
       vmadid(i1)=dataout(idatc)
      end do
c
      do i1=1,nimp
       idatc=idatc+1
       vmadiq(i1)=dataout(idatc)
      end do
c
      do i1=1,nimp
       idatc=idatc+1
       enba(i1)=dataout(idatc)
      end do
c
      do i1=1,nimp
      do i2=1,ne
       idatc=idatc+1
       enbdiffa(i1,i2)=dataout(idatc)
      end do
      end do
c
      do i1=1,nimp
       idatc=idatc+1
       denba(i1)=dataout(idatc)
      end do
c
      do i1=1,nimp
       idatc=idatc+1
       enorba(i1)=dataout(idatc)
      end do
c
      do i1=1,nrad
      do i2=1,nimp
       idatc=idatc+1
       rhova(i1,i2)=dataout(idatc)
      end do
      end do
c
      do i1=1,nrad
      do i2=1,2
      do i3=1,nimp
       idatc=idatc+1
       rhospa(i1,i2,i3)=dataout(idatc)
      end do
      end do
      end do
c
      do i1=1,nrad
      do i2=1,2
      do i3=1,nimp
       idatc=idatc+1
       rhodspa(i1,i2,i3)=dataout(idatc)
      end do
      end do
      end do
c
      do i1=1,nrad
      do i2=1,nimp
       idatc=idatc+1
       rhomaga(i1,i2)=dataout(idatc)
      end do
      end do
c
      do i1=1,kmymax
      do i2=1,nimp
      do i3=1,ne
       idatc=idatc+1
       dosmaga(i1,i2,i3)=dataout(idatc)
      end do
      end do
      end do
c
      do i1=1,kmymax
      do i2=1,nimp
      do i3=1,3
       idatc=idatc+1
       spin_magvpa(i1,i2,i3)=dataout(idatc)
      end do
      end do
      end do
c
      do i1=1,nimp
      do i2=1,3
       idatc=idatc+1
       spin_magva(i1,i2)=dataout(idatc)
      end do
      end do
c
      do i1=1,kmymax
      do i2=1,nimp
      do i3=1,3
       idatc=idatc+1
       orb_magvpa(i1,i2,i3)=dataout(idatc)
      end do
      end do
      end do
c
      do i1=1,nimp
      do i2=1,3
       idatc=idatc+1
       orb_magva(i1,i2)=dataout(idatc)
      end do
      end do
c
      do i1=1,4
      do i2=1,nimp
      do i3=1,nimp
       idatc=idatc+1
       vij(i1,i2,i3)=dataout(idatc)
      end do
      end do
      end do
c
      do i1=1,nimp
       idatc=idatc+1
       gradth(i1)=dataout(idatc)
      end do
c
      do i1=1,nimp
       idatc=idatc+1
       gradph(i1)=dataout(idatc)
      end do
c
      do i1=1,1
       idatc=idatc+1
       ebl(i1)=dataout(idatc)
      end do
c
      do i1=1,kmymax
      do i2=1,nimp
        idatc=idatc+1
        enpba(i1,i2)=dataout(idatc)
      end do
      end do
c
c
      cidatc=0
      do i1=1,lmsup
      do i2=1,nimp
        cidatc=cidatc+1
        qmoma(i1,i2)=cdataout(cidatc)
      end do
      end do
c
      do i1=1,ne
        cidatc=cidatc+1
        detl2(i1)=cdataout(cidatc)
      end do
c
      do i1=1,nimp
      do i2=1,ne
        cidatc=cidatc+1
        nenba(i1,i2)=cdataout(cidatc)
      end do
      end do
c
      do i1=1,nimp
      do i2=1,ne
        cidatc=cidatc+1
        nqva(i1,i2)=cdataout(cidatc)
      end do
      end do
c
      deallocate(datain,stat=AllocateStatus)
      call alloccheck( AllocateStatus,'datain dealloc in collectdat' )
      deallocate(dataout,stat=AllocateStatus)
      call alloccheck( AllocateStatus,'dataout dealloc in collectdat' )
      deallocate(cdatain,stat=AllocateStatus)
      call alloccheck( AllocateStatus,'cdatain dealloc in collectdat' )
      deallocate(cdataout,stat=AllocateStatus)
      call alloccheck( AllocateStatus,'cdataout dealloc in collectdat' )
c
      return
      end
c
c-------------------------------------------
