c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine phystau2(tau,ntaudim,tinv1,tinv2,alphakkr1,alphakkrh1,
     >                    alphakkr2,alphakkrh2,
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
      complex*16 tau,tinv1,tinv2,alphakkr1,alphakkr2,lmat,rmat,cmat,help
      complex*16 alphakkrh1,alphakkrh2
      complex*16 alphal,alphalp
      dimension tau(ntaudim,ntaudim)
      dimension alphakkr1(0:lmaxp),alphakkrh1(0:lmaxp)
      dimension alphakkr2(0:lmaxp),alphakkrh2(0:lmaxp)
      dimension lmat(dbogomaxp,dbogomaxp),rmat(dbogomaxp,dbogomaxp)
      dimension cmat(dbogomaxp,dbogomaxp),help(dbogomaxp,dbogomaxp)
      dimension tinv1(dbogomaxp,dbogomaxp),tinv2(dbogomaxp,dbogomaxp)
c     dimension ldex(50)
c      data ldex/0,0,
c     *          1,1,1,1,1,1,
c     *          2,2,2,2,2,2,2,2,2,2,
c     *          3,3,3,3,3,3,3,3,3,3,3,3,3,3,
c     *          4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4/
c
      lmmax=(lmax+1)**2
      kmymax=2*lmmax
c
      call czero(lmat,dbogomaxp*dbogomaxp)
      call czero(cmat,dbogomaxp*dbogomaxp)
      call czero(rmat,dbogomaxp*dbogomaxp)
c
      do kmy=1,(2*kmymax)
        if (kmy>kmymax) then
          alphal=alphakkrh1(ldex(kmy-kmymax))
        else
          alphal=alphakkr1(ldex(kmy))
        end if
        do kmyp=1,(2*kmymax)
          help(kmy,kmyp)=tau(kmy,kmyp)
          if (kmyp>kmymax) then
            alphalp=alphakkrh2(ldex(kmyp-kmymax))
          else
            alphalp=alphakkr2(ldex(kmyp))
          end if
          lmat(kmy,kmyp)=alphal*tinv1(kmy,kmyp)
          rmat(kmy,kmyp)=tinv2(kmy,kmyp)*alphalp
          cmat(kmy,kmyp)=alphal*tinv1(kmy,kmyp)*alphalp
        end do
        lmat(kmy,kmy)=1.d0+lmat(kmy,kmy)
        rmat(kmy,kmy)=1.d0+rmat(kmy,kmy)
        cmat(kmy,kmy)=alphal+cmat(kmy,kmy)
      end do
      call tripmt(lmat,help,rmat,2*kmymax,2*kmymax,dbogomaxp)
      call repldim(tau,help,ntaudim,ntaudim,dbogomaxp)
c
      if(idiag.ne.1) return
      do kmy=1,(2*kmymax)
      do kmyp=1,(2*kmymax)
        tau(kmy,kmyp)=help(kmy,kmyp)-cmat(kmy,kmyp)
      end do
      end do
c
      return
      end
