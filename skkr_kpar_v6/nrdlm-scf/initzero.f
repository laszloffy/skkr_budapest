c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine initzero(qvpa,qvpb,qva,qvb,enba,enbb,qmoma,qmomb,
     & rhova,rhovb)
c
      implicit real*8(a-h,o-z)
      include '../param.h'
c
      dimension qva(mintfc),qvb(mintfc)
      dimension qvpa(0:lmaxp,mintfc),qvpb(0:lmaxp,mintfc)   
      dimension enba(mintfc),enbb(mintfc)
      dimension rhova(nrad,mintfc),rhovb(nrad,mintfc)
c
      complex*16 qmoma(lmsup,mintfc),qmomb(lmsup,mintfc)
c
c     -----------------------------
      call rzero(qvpa,l1maxp*mintfc)
      call rzero(qvpb,l1maxp*mintfc)
      call rzero(qva,mintfc)
      call rzero(qvb,mintfc)
      call rzero(enba,mintfc)
      call rzero(enbb,mintfc)
      call czero(qmoma,lmsup*mintfc)
      call czero(qmomb,lmsup*mintfc)
      call rzero(rhova,nrad*mintfc)
      call rzero(rhovb,nrad*mintfc)
c     -----------------------------
c
      return
      end
