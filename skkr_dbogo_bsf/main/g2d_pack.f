c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine g2d_pack(g2d,g2de,g2dh,lmmax,mmprc1)
c
c     Set up the screened BdG structure constant matrix from the
c     electron-like and hole-like screened structure constants.
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      complex*16 g2d(bogomaxp,bogomaxp,-minprc:minprc,minprc,0:mmprc1+1)
      complex*16 g2de(lmmaxp,lmmaxp,-minprc:minprc,minprc,0:mmprc1+1)
      complex*16 g2dh(lmmaxp,lmmaxp,-minprc:minprc,minprc,0:mmprc1+1)
c
        g2d(:,:,:,:,:)=0.0d0
        g2d(1:lmmax,1:lmmax,:,:,:)=g2de(:,:,:,:,:)
        g2d((lmmax+1):(2*lmmax),(lmmax+1):(2*lmmax),:,:,:)=
     >                                   g2dh(:,:,:,:,:)
c

c
c
      return
      end
