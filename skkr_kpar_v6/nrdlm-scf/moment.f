c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine moment(
     &lmax,rs,dx,ns,yz,yj,coeff,lambda,tau,zmom,irreg,para)
c======================
c
c input:  lmax  - maximum of angular momentum index
c         rs,dx,ns - as in 'readpot'
c         yz,yj   - regular and irregular radial solution*r
c         coeff - matrix of gaunt coefficients
c         lambda - corresponding angular momentum quantumnumber
c         tau - site--diagonal tau matrix
c output: zmom - (l,m)-like resolution of the moment's density 
c                corresponding to rws
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      logical para
c
      dimension rrr(nrad),rir(nrad)
      dimension rri(nrad),rii(nrad)
      dimension rad(nrad)
c
      complex*16 yz(nrad,0:lmaxp),yj(nrad,0:lmaxp)
      complex*16 tau(lmmaxp,lmmaxp)
      complex*16 utr(lmmaxp,lmmaxp),utrp(lmmaxp,lmmaxp)
      complex*16 coeff(lmmaxp,lmmaxp),cxi,cxr,cox
      complex*16 sum(lmmaxp,lmmaxp,2)
      complex*16 gf(lmmaxp,lmmaxp)
      complex*16 zmom(lmmaxp)
c
      data small/1.0d-12/
c
      pi=4.d0*datan(1.d0)
      fac=-1.d0/pi
      if(para) fac=2.0d0*fac
      lmmax=(lmax+1)*(lmax+1)
      x=dlog(rs)-(ns-1)*dx
      do ii=1,ns
         rad(ii)=dexp(x)
         x=x+dx
      end do
c
c calculate matrices of radial integrals
c
      do l1=0,lmax
      do l2=0,lmax
c
        call rzero(rrr,ns)
        call rzero(rir,ns)
        call rzero(rri,ns)
        call rzero(rii,ns)
        do i=1,ns
            rlam=rad(i)**lambda
            cxr=yz(i,l1)*yz(i,l2)*rlam
            cxi=yz(i,l1)*yj(i,l2)*rlam
            rrr(i)=rrr(i)+dreal(cxr)
            rir(i)=rir(i)+dimag(cxr)
            rri(i)=rri(i)+dreal(cxi)
            rii(i)=rii(i)+dimag(cxi)
        end do
        r1rr=rsimp(rrr,rad,ns,dx)
        r1ir=rsimp(rir,rad,ns,dx)
        r1ri=rsimp(rri,rad,ns,dx)
        r1ii=rsimp(rii,rad,ns,dx)
        cxr=dcmplx(r1rr,r1ir)
        cxi=dcmplx(r1ri,r1ii)
c
        do m1=-l1,l1
             lm1=l1*(l1+1)+m1+1
        do m2=-l2,l2
             lm2=l2*(l2+1)+m2+1
             sum(lm1,lm2,1)=fac*cxr*coeff(lm1,lm2)
             sum(lm1,lm2,2)=fac*cxi*coeff(lm1,lm2)
        end do
        end do
c
      end do
      end do
c
c multiply regular matrix by tau-matrix and substract irregular matrix
c (note that the off-diagonal elements of the irregular matrix are dummy!!!)
c
      call repl(gf,sum(1,1,1),lmmax,lmmaxp)
      call doubmt(gf,tau,lmmax,lmmaxp)
      if(irreg.eq.1) call submat(gf,sum(1,1,2),lmmax,lmmaxp)
c
c transform to real spherical harmonics for lambda=0
c
      if(lambda.ne.0) goto 100
      call ytrafo(utr,utrp,lmax)
      call tripmt(utr,gf,utrp,lmmax,lmmax,lmmaxp)
c
  100 continue
c find diagonals 
c
      do lm=1,lmmax
        zmom(lm)=gf(lm,lm)
      end do
c
      return
      end
