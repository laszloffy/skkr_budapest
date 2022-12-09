c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine moment(lmax,rs,dx,ns,gz,fz,gj,fj,
     >nuz,indz,coeff,lambda,tau,zmom,lms,irreg)    
c======================
c
c input:  lmax  - maximum of angular momentum index
c         rs,dx,ns - as in 'readpot'
c         gz,fz   - big and small component of regular radial solution*r
c         gj,fj   - big and small component of irreg. radial solution*r
c         nuz - no. of (kap',my') components for (kap,my)
c         indz - selects (kap',my') for (kap,my)
c         coeff - matrix of generalized gaunt coefficients
c         lambda - corresponding angular momentum quantumnumber
c         tau - site--diagonal tau matrix
c         lms - if .true. output in (l,m,s) representation
c output: zmommt, zmomws - 
c         (kappa,my)-like or (l,m,s)-like resolution of the moment's
c         density corresponding to rws.
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      logical lms
      dimension nuz(kmymaxp),indz(nuzp,kmymaxp)
      dimension rrr(nrad),rir(nrad)
      dimension rri(nrad),rii(nrad)
      dimension rad(nrad)
c
      complex*16 gz(nrad,nuzp,kmymaxp),fz(nrad,nuzp,kmymaxp)
      complex*16 gj(nrad,nuzp,kmymaxp),fj(nrad,nuzp,kmymaxp)
      complex*16 tau(kmymaxp,kmymaxp)
      complex*16 coeff(kmymaxp,kmymaxp),cxi,cxr,cox
      complex*16 sum(kmymaxp,kmymaxp,2)
      complex*16 gf(kmymaxp,kmymaxp)
      complex*16 zmom(kmymaxp),zmomfnr(kmymaxp,kmymaxp)
      complex*16 zmomf(kmymaxp,kmymaxp)
      complex*16 utr(lmmaxp,lmmaxp)
      complex*16 utrp(lmmaxp,lmmaxp)
      complex*16 zdostmp1(lmmaxp,lmmaxp)
      complex*16 zdostmp2(lmmaxp,lmmaxp)
c
      data small/1.0d-12/
c
      nl=lmax+1
      nl2=nl*nl
      kmax=2*lmax+1
      kmymax=2*nl2
c
      pi=4.d0*datan(1.d0)
      fac=-1.d0/pi
c      kmymax=2*(lmax+1)*(lmax+1)
      x=dlog(rs)-(ns-1)*dx
      do ii=1,ns
         rad(ii)=dexp(x)
         x=x+dx
      end do
c
c calculate matrices of radial integrals
c
      do kmy=1,kmymax
      do kmyp=1,kmymax
c
        call rzero(rrr,ns)
        call rzero(rir,ns)
        call rzero(rri,ns)
        call rzero(rii,ns)
        do nu=1,nuz(kmy)
          kmy1=indz(nu,kmy)
        do nup=1,nuz(kmyp)
          kmyp1=indz(nup,kmyp)
          cox=fac*coeff(kmy1,kmyp1)
          if(cdabs(cox).gt.small) then
          do i=1,ns
            rlam=rad(i)**lambda
            cxr=(gz(i,nu,kmy)*gz(i,nup,kmyp)+
     >          fz(i,nu,kmy)*fz(i,nup,kmyp))*cox*rlam
            cxi=(gz(i,nu,kmy)*gj(i,nup,kmyp)+
     >          fz(i,nu,kmy)*fj(i,nup,kmyp))*cox*rlam
            rrr(i)=rrr(i)+dreal(cxr)
            rir(i)=rir(i)+dimag(cxr)
            rri(i)=rri(i)+dreal(cxi)
            rii(i)=rii(i)+dimag(cxi)
          end do
          end if
c
        end do
        end do
c
        r1rr=rsimp(rrr,rad,ns,dx)
        r1ir=rsimp(rir,rad,ns,dx)
        r1ri=rsimp(rri,rad,ns,dx)
        r1ii=rsimp(rii,rad,ns,dx)
c
        sum(kmy,kmyp,1)=dcmplx(r1rr,r1ir)
        sum(kmy,kmyp,2)=dcmplx(r1ri,r1ii)
c         
c next kmyp!
      end do
c next kmy!
      end do
c
c multiply regular matrix by tau-matrix and substract irregular matrix
c
      call repl(gf,sum(1,1,1),kmymax,kmymaxp)
      call doubmt(gf,tau,kmymax,kmymaxp)
      if(irreg.eq.1) call submat(gf,sum(1,1,2),kmymax,kmymaxp)
c
c find diagonals in the required representation
c
      if(lms) then
        do k1=1,kmymax
        do k2=1,kmymax
              zmomf(k1,k2)=gf(k1,k2)
        end do
        end do
c convert zdosf to lms-representation with respect to the complex spherical harmonics
          call replmsf(zmomfnr,zmomf,lmax)
c initialize transformation matricces from the complex Y_lm to the real Y_lm
          call ytrafo(utr,utrp,lmax)   
c Set temporary matrices for the spin resolved zdos
           do iii=1,nl2
           do jjj=1,nl2
            zdostmp1(iii,jjj)=zmomfnr(iii,jjj)  !spin 1
            zdostmp2(iii,jjj)=zmomfnr(lmmaxp+iii,lmmaxp+jjj) !spin 2
           end do
           end do
c transformation of the two spin channels to the real Y_lm
          call tripmt(utr,zdostmp1,utrp,nl2,nl2,lmmaxp)  ! spin 1
          call tripmt(utr,zdostmp2,utrp,nl2,nl2,lmmaxp) ! spin 2
c write back into the original matrices
            do iii=1,nl2
            do jjj=1,nl2
             zmomfnr(iii,jjj)=zdostmp1(iii,jjj)
             zmomfnr(lmmaxp+iii,lmmaxp+jjj)=zdostmp2(iii,jjj)
            end do
            end do
        do k12=1,kmymax
          zmom(k12)=zmomfnr(k12,k12)
        end do
      else
        do kmy=1,kmymax
          zmom(kmy)=gf(kmy,kmy)
        end do
      end if
c
      return
      end
