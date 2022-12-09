c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine phystau(tau,tinv1,tinv2,alphakkr,lmax,idiag)
c=======================
c
c convert diagonal tau matrix from screened representation
c into unscreened one
c
c input:  tau          ... screened layer diagonal tau matrix
c                         (will be replaced by unscreened one on output)
c         tinv1,tinv2  ... inverse of screened scattering matrix
c         alphakkr     ... screening parameters 
c         lmax         ... actual max. value of angular momentum index l
c (output: tau)
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      complex*16 tau,tinv1,tinv2,alphakkr,lmat,rmat,cmat
      complex*16 alphal,alphalp
      dimension tau(lmmaxp,lmmaxp),alphakkr(0:lmaxp)
      dimension lmat(lmmaxp,lmmaxp),rmat(lmmaxp,lmmaxp)
      dimension cmat(lmmaxp,lmmaxp)
      dimension tinv1(lmmaxp,lmmaxp),tinv2(lmmaxp,lmmaxp)
      dimension lml(lmmaxp),lmm(lmmaxp)
      common/mome/lml,lmm
c
      lmmax=(lmax+1)**2
c
      call czero(lmat,lmmaxp*lmmaxp)
      call czero(cmat,lmmaxp*lmmaxp)
      call czero(rmat,lmmaxp*lmmaxp)
c
      do l=1,lmmax
        alphal=alphakkr(lml(l))
        do lp=1,lmmax
          alphalp=alphakkr(lml(lp))
          lmat(l,lp)=alphal*tinv1(l,lp)
          rmat(l,lp)=tinv2(l,lp)*alphalp
          cmat(l,lp)=alphal*tinv1(l,lp)*alphalp
        end do
        lmat(l,l)=1.d0+lmat(l,l)
        rmat(l,l)=1.d0+rmat(l,l)
        cmat(l,l)=alphal+cmat(l,l)
      end do
      call tripmt(lmat,tau,rmat,lmmax,lmmax,lmmaxp)
c
      if(idiag.ne.1) return
      do l=1,lmmax
      do lp=1,lmmax
        tau(l,lp)=tau(l,lp)-cmat(l,lp)
      end do
      end do
c
      return
      end
