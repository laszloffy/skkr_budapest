c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine initia
c
      implicit real*8(a-h,o-z)
c
      call lmfill
      call init
      call clebsch
      call gafill
c
      return
      end
