c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine lmfill
c======================
c
c fill the lml vector by 'l' values, the lmm by 'm'
c values belonging to  lm=(l,m) .
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      dimension lml(lmmaxp),lmm(lmmaxp)
      common/mome/lml,lmm
c
      lm = 0
      do 1 l=0,lmaxp
      do 1 m=-l,l
      lm=lm+1
      lml(lm) = l
    1 lmm(lm) = m
c
      return
      end
