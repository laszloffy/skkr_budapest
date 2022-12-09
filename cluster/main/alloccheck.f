c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine alloccheck(AllocateStat,arrayname)
c
      integer*4 AllocateStat
      character*50 arrayname
c
      if(AllocateStat.ne.0) then
        write(6,*) 'Insufficient space for array: ',arrayname
        stop
      endif
      end
