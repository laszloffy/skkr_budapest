c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine phystau(tau,ntaudim,tinv1,tinv2,alphakkr,alphakkrh,
     >                   lmax,idiag)
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
      logical tautest
      complex*16 tau,tinv1,tinv2,alphakkr,lmat,rmat,cmat,help
      complex*16 alphal,alphalp,alphakkrh
      dimension tau(ntaudim,ntaudim)
      dimension alphakkr(0:lmaxp),alphakkrh(0:lmaxp)
      dimension lmat(dbogomaxp,dbogomaxp),rmat(dbogomaxp,dbogomaxp)
      dimension cmat(dbogomaxp,dbogomaxp),help(dbogomaxp,dbogomaxp)
      dimension tinv1(dbogomaxp,dbogomaxp),tinv2(dbogomaxp,dbogomaxp)
c
      data tol/1.0d-10/
c
      lmmax=(lmax+1)**2
      kmymax=2*lmmax
c
      tautest=.false.
      if(tautest) then
      write(6,*) '<phystau> : ntaudim=',ntaudim
      write(6,*) '<phystau> : dbogomaxp=',dbogomaxp
      write(6,*) '<phystau> : 2*kmymax=',2*kmymax
      end if
c
      call czero(lmat,dbogomaxp*dbogomaxp)
      call czero(cmat,dbogomaxp*dbogomaxp)
      call czero(rmat,dbogomaxp*dbogomaxp)
c
      do kmy=1,(2*kmymax)
        if (kmy>kmymax) then
          alphal=alphakkrh(ldex(kmy-kmymax))
        else
          alphal=alphakkr(ldex(kmy))
        end if
        do kmyp=1,(2*kmymax)
          help(kmy,kmyp)=tau(kmy,kmyp)
          if (kmyp>kmymax) then
            alphalp=alphakkrh(ldex(kmyp-kmymax))
          else
             alphalp=alphakkr(ldex(kmyp))
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
c
      if(idiag.ne.1) then
        tau(1:2*kmymax,1:2*kmymax)=help(1:2*kmymax,1:2*kmymax)
        return
      end if
      do kmy=1,(2*kmymax)
      do kmyp=1,(2*kmymax)
        tau(kmy,kmyp)=help(kmy,kmyp)-cmat(kmy,kmyp)
      end do
      end do
c
c      write(6,*) 'phystau: '
c      call outmat1(tau,2*kmymax,2*kmymax,dbogomaxp,tol,6)
c
      return
      end
