c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine phystau2(tau,tinv1,tinv2,alphakkr1,alphakkr2,
     >                    lmax,idiag)
c=======================
c
c convert principal layer diagonal tau matrix from screened 
c representation into unscreened one
c
c input:  tau          ... screened tau matrix
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
      complex*16 tau,tinv1,tinv2,alphakkr1,alphakkr2,lmat,rmat,cmat
      complex*16 alphal,alphalp
      dimension tau(kmymaxp,kmymaxp),alphakkr1(0:lmaxp)
      dimension alphakkr2(0:lmaxp)
      dimension lmat(kmymaxp,kmymaxp),rmat(kmymaxp,kmymaxp)
      dimension cmat(kmymaxp,kmymaxp)
      dimension tinv1(kmymaxp,kmymaxp),tinv2(kmymaxp,kmymaxp)
      dimension ldex(50)
      data ldex/0,0,
     *          1,1,1,1,1,1,
     *          2,2,2,2,2,2,2,2,2,2,
     *          3,3,3,3,3,3,3,3,3,3,3,3,3,3,
     *          4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4/
c
      lmmax=(lmax+1)**2
      kmymax=2*lmmax
c
      call czero(lmat,kmymaxp*kmymaxp)
      call czero(cmat,kmymaxp*kmymaxp)
      call czero(rmat,kmymaxp*kmymaxp)
c
      do kmy=1,kmymax
        alphal=alphakkr1(ldex(kmy))
        do kmyp=1,kmymax
          alphalp=alphakkr2(ldex(kmyp))
          lmat(kmy,kmyp)=alphal*tinv1(kmy,kmyp)
          rmat(kmy,kmyp)=tinv2(kmy,kmyp)*alphalp
          cmat(kmy,kmyp)=alphal*tinv1(kmy,kmyp)*alphalp
        end do
        lmat(kmy,kmy)=1.d0+lmat(kmy,kmy)
        rmat(kmy,kmy)=1.d0+rmat(kmy,kmy)
        cmat(kmy,kmy)=alphal+cmat(kmy,kmy)
      end do
      call tripmt(lmat,tau,rmat,kmymax,kmymax,kmymaxp)
c
      if(idiag.ne.1) return
      do kmy=1,kmymax
      do kmyp=1,kmymax
        tau(kmy,kmyp)=tau(kmy,kmyp)-cmat(kmy,kmyp)
      end do
      end do
c
      return
      end
