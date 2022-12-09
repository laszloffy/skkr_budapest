c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine magnet(lmax,rs,dx,ns,gz,fz,gj,fj,
     >       nuz,indz,tau,gcoeff,gbcoeff,zmom,lms,irreg)    
c======================
c
c input:  lmax  - maximum of angular momentum index
c         rs,dx,ns - as in 'readpot'
c         gz,fz - big and small component of regular radial solution * r
c         gj,fj - big and small component of irreg. radial solution * r
c         nuz - no. of (kap',my') components for (kap,my)
c         indz - selects (kap',my') for (kap,my)
c         gcoeff, gbcoeff - angular momentum matrixelements
c         tau - site--diagonal tau matrix
c         lms - if .true. output in (l,m,s) representation
c         irreg - if 1, subtract irregular rad.int. part (do, if tau is
c                 a single tau matrix and not a difference!)
c output: zmom
c         (kappa,my)-like or (l,m,s)-like resolution of the 
c         magnetization density corresponding to rs.
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
      complex*16 cxr,cxi,gcox,gbcox
      complex*16 sum(kmymaxp,kmymaxp,2)
      complex*16 gf(kmymaxp,kmymaxp)
      complex*16 zmom(kmymaxp)
      complex*16 gcoeff(kmymaxp,kmymaxp),gbcoeff(kmymaxp,kmymaxp)
c
      data small/1.0d-12/
c
      pi=4.d0*datan(1.d0)
      fac=-1.d0/pi
      kmymax=2*(lmax+1)*(lmax+1)
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
        call rzero(rrr,nrad)
        call rzero(rir,nrad)
        call rzero(rri,nrad)
        call rzero(rii,nrad)
        do nu=1,nuz(kmy)
          kmy1=indz(nu,kmy)
        do nup=1,nuz(kmyp)
          kmyp1=indz(nup,kmyp)
          gcox=fac*gcoeff(kmy1,kmyp1)
          gbcox=fac*gbcoeff(kmy1,kmyp1)
c
          if(cdabs(gcox).gt.small) then
          do i=1,ns
            cxr=gz(i,nu,kmy)*gz(i,nup,kmyp)*gcox
            cxi=gz(i,nu,kmy)*gj(i,nup,kmyp)*gcox
            rrr(i)=rrr(i)+dreal(cxr)
            rir(i)=rir(i)+dimag(cxr)
            rri(i)=rri(i)+dreal(cxi)
            rii(i)=rii(i)+dimag(cxi)
          end do
          end if
          if(cdabs(gbcox).gt.small) then
          do i=1,ns
            cxr=fz(i,nu,kmy)*fz(i,nup,kmyp)*gbcox
            cxi=fz(i,nu,kmy)*fj(i,nup,kmyp)*gbcox
            rrr(i)=rrr(i)-dreal(cxr)
            rir(i)=rir(i)-dimag(cxr)
            rri(i)=rri(i)-dreal(cxi)
            rii(i)=rii(i)-dimag(cxi)
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
c sum(.,.,1) contains (big - small) component of the regular solution
c of the Dirac equation
c sum(.,.,2) contains the same for the irregular solution of the Dirac
c equation
c
        sum(kmy,kmyp,1)=dcmplx(r1rr,r1ir)
        sum(kmy,kmyp,2)=dcmplx(r1ri,r1ii)
c         
c next kmyp!
      end do
c next kmy!
      end do
c
c multiply regular matrix by tau-matrix and substract
c irregular matrix
c
      call repl(gf,sum(1,1,1),kmymax,kmymaxp)
      call doubmt(gf,tau,kmymax,kmymaxp)
      if(irreg.eq.1) call submat(gf,sum(1,1,2),kmymax,kmymaxp)
c
c find diagonals in the required representation
c
      if(lms) then 
        call replms(zmom,gf,lmax)
      else
        do kmy=1,kmymax
          zmom(kmy)=gf(kmy,kmy)
        end do
      end if
c
      return
      end
