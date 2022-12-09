c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine collect_data(kmymax,nimp,ne,
     > dosimp,doshimp,dosehimp,
     > dosheimp,dostehimp,dostheimp,
     > qvpimp,qvhpimp,qvimp,qvhimp,qvehimp,qvheimp,
     > qvtehimp,qvtheimp,
c     qva,
c     > vmadid,vmadiq,enba,enbdiffa,denba,nenba,enorba,qmoma,
c     > rhova,rhospa,rhodspa,rhomaga,
c     > dosmaga,spin_magvpa,spin_magva,orb_magvpa,orb_magva,
c     > vij,gradph,gradth,
c     > ebl,detl2)
     > tripletmat)
c
      implicit none
      include '../param.h'
#ifdef MPIP
      include 'mpif.h'
#endif
c
c MPI parameters (esimon)
      integer root, myrank, nprocs, ierror
      common/mpi0/root,myrank,nprocs
c
      integer kmymax,nimp,ne
c
      real*8 dosimp(kmymaxp,nimp,me)
      real*8 doshimp(kmymaxp,nimp,me)
      real*8 dosehimp(kmymaxp,nimp,me)
      real*8 dosheimp(kmymaxp,nimp,me)
      real*8 dostehimp(kmymaxp,nimp,me)
      real*8 dostheimp(kmymaxp,nimp,me)
      real*8 qvpimp(kmymaxp,nimp)
      real*8 qvhpimp(kmymaxp,nimp)
      real*8 qvimp(nimp)
      real*8 qvhimp(nimp)
      real*8 qvehimp(nimp)
      real*8 qvheimp(nimp)
      real*8 qvtehimp(nimp)
      real*8 qvtheimp(nimp)
      real*8 tripletmat(lmmaxp,lmmaxp,3,nimp,me)
c      complex*16 nqva(mimp,me)
c      real*8 qvdiffa(mimp,me)
c      real*8 vmadid(mimp)
c      real*8 vmadiq(mimp)
c      real*8 enba(mimp)
c      real*8 enbdiffa(mimp,me)
c      real*8 denba(mimp)
c      real*8 enorba(mimp)
c      complex*16 qmoma(lmsup,mimp)
c      complex*16 detl2(me)
c      complex*16 nenba(mimp,me)
c      real*8 rhova(nrad,mimp)
c      real*8 rhospa(nrad,2,mimp)
c      real*8 rhodspa(nrad,2,mimp)
c      real*8 rhomaga(nrad,mimp)
c      real*8 dosmaga(kmymaxp,mimp,me)
c      real*8 spin_magvpa(kmymaxp,mimp,3)
c      real*8 spin_magva(mimp,3)
c      real*8 orb_magvpa(kmymaxp,mimp,3)
c      real*8 orb_magva(mimp,3)
c      real*8 vij(4,mimp,mimp)
c      real*8 gradth(mimp)
c      real*8 gradph(mimp)
c      real*8 ebl(1)
c
      integer*4 idat,idatc
      integer*4 cidat,cidatc
      integer*4 i1,i2,i3,i4,i5
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
c      complex*16, allocatable :: cdatain(:)
c      complex*16, allocatable :: cdataout(:)
c
      integer AllocateStatus
c
#ifdef MPIP
      idat = kmymaxp*nimp*ne*6+ !dos dosh doseh doshe dosteh dosthe
     > kmymaxp*nimp*2+ !qvp,qvph
     > nimp*6+ ! qv qvh qveh qvhe qvteh qvthe
     > lmmaxp*lmmaxp*3*nimp*ne
c      idat=kmymax*nimp*ne+kmymax*nimp+nimp+nimp*ne+nimp+nimp+nimp
c     > +nimp*ne+nimp+nimp+nrad*nimp+nrad*2*nimp+nrad*2*nimp
c     > +nrad*nimp+kmymax*nimp*ne+kmymax*nimp*3+nimp*3
c     > +kmymax*nimp*3+nimp*3+4*nimp*nimp+nimp+nimp+1
c
c      cidat=lmsup*nimp+ne+nimp*ne+nimp*ne
c
      ALLOCATE ( datain(idat),STAT = AllocateStatus)
      ALLOCATE ( dataout(idat),STAT = AllocateStatus)
c      ALLOCATE ( cdatain(cidat),STAT = AllocateStatus)
c      ALLOCATE ( cdataout(cidat),STAT = AllocateStatus)
c
c ----------------------------------------------
c
      idatc=0
c
      do i1=1,kmymaxp
      do i2=1,nimp
      do i3=1,ne
         idatc=idatc+1
         datain(idatc)=dosimp(i1,i2,i3)
      end do
      end do
      end do
c
      do i1=1,kmymaxp
      do i2=1,nimp
      do i3=1,ne
         idatc=idatc+1
         datain(idatc)=doshimp(i1,i2,i3)
      end do
      end do
      end do
c
      do i1=1,kmymaxp
      do i2=1,nimp
      do i3=1,ne
         idatc=idatc+1
         datain(idatc)=dosehimp(i1,i2,i3)
      end do
      end do
      end do
c
      do i1=1,kmymaxp
      do i2=1,nimp
      do i3=1,ne
         idatc=idatc+1
         datain(idatc)=dosheimp(i1,i2,i3)
      end do
      end do
      end do
c
      do i1=1,kmymaxp
      do i2=1,nimp
      do i3=1,ne
         idatc=idatc+1
         datain(idatc)=dostehimp(i1,i2,i3)
      end do
      end do
      end do
c
      do i1=1,kmymaxp
      do i2=1,nimp
      do i3=1,ne
         idatc=idatc+1
         datain(idatc)=dosehimp(i1,i2,i3)
      end do
      end do
      end do
c
      do i1=1,kmymaxp
      do i2=1,nimp
         idatc=idatc+1
         datain(idatc)=qvpimp(i1,i2)
      end do
      end do
c
      do i1=1,kmymaxp
      do i2=1,nimp
         idatc=idatc+1
         datain(idatc)=qvhpimp(i1,i2)
      end do
      end do
c
      do i1=1,nimp
         idatc=idatc+1
         datain(idatc)=qvimp(i1)
      end do
c
      do i1=1,nimp
         idatc=idatc+1
         datain(idatc)=qvhimp(i1)
      end do
c
      do i1=1,nimp
         idatc=idatc+1
         datain(idatc)=qvehimp(i1)
      end do
c
      do i1=1,nimp
         idatc=idatc+1
         datain(idatc)=qvheimp(i1)
      end do
c
      do i1=1,nimp
         idatc=idatc+1
         datain(idatc)=qvtehimp(i1)
      end do
c
      do i1=1,nimp
         idatc=idatc+1
         datain(idatc)=qvtheimp(i1)
      end do
      do i1=1,lmmaxp
      do i2=1,lmmaxp
      do i3=1,3
      do i4=1,nimp
      do i5=1,ne
         idatc=idatc+1
         datain(idatc)=tripletmat(i1,i2,i3,i4,i5)
      end do
      end do
      end do
      end do
      end do
c --------------------------------------
      do i1=1,idat
        dataout(i1) = 0.0
      end do
c
       if(myrank.eq.root) write(6,*) 'idat:',idat, idatc
c
      call mpi_allreduce(datain,dataout,idat,
     > mpi_real8,mpi_sum,mpi_comm_world,ierror)
c
      idatc=0
c
      do i1=1,kmymax
      do i2=1,nimp
      do i3=1,ne
       idatc=idatc+1
       dosimp(i1,i2,i3)=dataout(idatc)
      end do
      end do
      end do
c
      do i1=1,kmymax
      do i2=1,nimp
      do i3=1,ne
       idatc=idatc+1
       doshimp(i1,i2,i3)=dataout(idatc)
      end do
      end do
      end do
c
      do i1=1,kmymax
      do i2=1,nimp
      do i3=1,ne
       idatc=idatc+1
       dosehimp(i1,i2,i3)=dataout(idatc)
      end do
      end do
      end do
c
      do i1=1,kmymax
      do i2=1,nimp
      do i3=1,ne
       idatc=idatc+1
       dosheimp(i1,i2,i3)=dataout(idatc)
      end do
      end do
      end do
c
      do i1=1,kmymax
      do i2=1,nimp
      do i3=1,ne
       idatc=idatc+1
       dostehimp(i1,i2,i3)=dataout(idatc)
      end do
      end do
      end do
c
      do i1=1,kmymax
      do i2=1,nimp
      do i3=1,ne
       idatc=idatc+1
       dostheimp(i1,i2,i3)=dataout(idatc)
      end do
      end do
      end do
c
      do i1=1,kmymax
      do i2=1,nimp
       idatc=idatc+1
       qvpimp(i1,i2)=dataout(idatc)
      end do
      end do
c
      do i1=1,kmymax
      do i2=1,nimp
       idatc=idatc+1
       qvhpimp(i1,i2)=dataout(idatc)
      end do
      end do
c
      do i1=1,nimp
       idatc=idatc+1
       qvimp(i1)=dataout(idatc)
      end do
c
      do i1=1,nimp
       idatc=idatc+1
       qvhimp(i1)=dataout(idatc)
      end do
c
      do i1=1,nimp
       idatc=idatc+1
       qvehimp(i1)=dataout(idatc)
      end do
c
      do i1=1,nimp
       idatc=idatc+1
       qvheimp(i1)=dataout(idatc)
      end do
c
      do i1=1,nimp
       idatc=idatc+1
       qvtehimp(i1)=dataout(idatc)
      end do
c
      do i1=1,nimp
       idatc=idatc+1
       qvtheimp(i1)=dataout(idatc)
      end do
      do i1=1,lmmaxp
      do i2=1,lmmaxp
      do i3=1,3
      do i4=1,nimp
      do i5=1,ne
         idatc=idatc+1
         tripletmat(i1,i2,i3,i4,i5)=dataout(idatc)
      end do
      end do
      end do
      end do
      end do
#endif
c
      return
      end
c
c-------------------------------------------
