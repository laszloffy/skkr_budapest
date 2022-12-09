c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine dens(lmax,rs,dx,ns,gz,fz,gj,fj,nuz,indz,tau,bopr,
     >                zrho,zrhosp,zrhodsp,zmagrho,zenorb)
c===================
c
c input:  lmax  - maximum of angular momentum index
c         rs,dx,ns - as in 'readpot'
c         gz,fz - big and small component of regular radial solution * r
c         gj,fj - big and small component of irreg. radial solution * r
c         nuz - no. of (kap',my') components for (kap,my)
c         indz - selects (kap',my') for (kap,my)
c         tau - site-diagonal tau matrix
c         bopr- vector potential for orbital polarization
c output: zrho - radial density distribution
c         zrhosp - spin resolution of radial density distribution
c         zrhodsp - only l=2 subspace is considered
c         zmagrho - radial magnetic density distribution
c         zenorb - integrand of orbital polarization energy
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      dimension nuz(kmymaxp),indz(nuzp,kmymaxp)
      dimension bopr(nrad,2),f1(nrad),f2(nrad),r(nrad)
c
      complex*16 gz(nrad,nuzp,kmymaxp),fz(nrad,nuzp,kmymaxp)
      complex*16 gj(nrad,nuzp,kmymaxp),fj(nrad,nuzp,kmymaxp)
      complex*16 zrho(nrad),zrhosp(nrad,2),zrhodsp(nrad,2),zmagrho(nrad)
      complex*16 tau(kmymaxp,kmymaxp)
      complex*16 sum1,sum2,zenorb
      complex*16 gmat(kmymaxp,kmymaxp),fmat(kmymaxp,kmymaxp)
      complex*16 gmat1(kmymaxp,kmymaxp),fmat1(kmymaxp,kmymaxp)
      complex*16 rmat(kmymaxp,kmymaxp),rlms(kmymaxp),cf(nrad)
c
      complex*16 sxcoeff(kmymaxp,kmymaxp),sxbcoeff(kmymaxp,kmymaxp)
      complex*16 sycoeff(kmymaxp,kmymaxp),sybcoeff(kmymaxp,kmymaxp)
      complex*16 szcoeff(kmymaxp,kmymaxp),szbcoeff(kmymaxp,kmymaxp)
      common/sigmat/sxcoeff,sxbcoeff,sycoeff,sybcoeff,szcoeff,szbcoeff
c
      complex*16 lxcoeff(kmymaxp,kmymaxp),lxbcoeff(kmymaxp,kmymaxp)
      complex*16 lycoeff(kmymaxp,kmymaxp),lybcoeff(kmymaxp,kmymaxp)
      complex*16 lzcoeff(kmymaxp,kmymaxp),lzbcoeff(kmymaxp,kmymaxp)
      common/lmat/lxcoeff,lxbcoeff,lycoeff,lybcoeff,lzcoeff,lzbcoeff
c
      data small/1.0d-12/
c
      pi=4.d0*datan(1.d0)
      fac=-1.d0/pi
      lmmax=(lmax+1)*(lmax+1)
      kmymax=2*lmmax
      call czero(zrho,nrad)
      call czero(zrhosp,2*nrad)
      call czero(zrhodsp,2*nrad)
      call czero(zmagrho,nrad)
c
      do i=1,ns
c
        do kmy1=1,kmymax
        do kmyp1=1,kmymax
          sum1=(0.d0,0.d0)
          sum2=(0.d0,0.d0)
          do kmy=1,kmymax
            do nu=1,nuz(kmy)
              if(kmy1.ne.indz(nu,kmy)) goto 101
          do kmyp=1,kmymax
            do nup=1,nuz(kmyp)
              if(kmyp1.ne.indz(nup,kmyp)) goto 100
c
              sum1=sum1+gz(i,nu,kmy)*tau(kmy,kmyp)*gz(i,nup,kmyp)
              sum2=sum2+fz(i,nu,kmy)*tau(kmy,kmyp)*fz(i,nup,kmyp)
              if(kmy.eq.kmyp) then
                sum1=sum1-gz(i,nu,kmy)*gj(i,nup,kmy)
                sum2=sum2-fz(i,nu,kmy)*fj(i,nup,kmy)
              end if
c
  100       continue
            end do
            end do
  101       continue
          end do
          end do
          gmat(kmy1,kmyp1)=sum1
          fmat(kmy1,kmyp1)=sum2
        end do
        end do
c
        call addmat1(gmat,fmat,rmat,kmymax,kmymaxp)
        call replms(rlms,rmat,lmax)
c
c charge density
        sum1=(0.d0,0.d0)
        do kmy=1,kmymax
          sum1=sum1+rmat(kmy,kmy)
        end do
        zrho(i)=fac*sum1
c
c spin-projected charge density
        sum1=(0.d0,0.d0)
        sum2=(0.d0,0.d0)
        do k=1,lmmax
          sum1=sum1+rlms(k)
          sum2=sum2+rlms(k+lmmax)
        end do
        zrhosp(i,1)=fac*sum1
        zrhosp(i,2)=fac*sum2
c
c spin-projected d-like (l=2) charge density
        sum1=(0.d0,0.d0)
        sum2=(0.d0,0.d0)
        do k=5,9
          sum1=sum1+rlms(k)
          sum2=sum2+rlms(k+lmmax)
        end do
        zrhodsp(i,1)=fac*sum1
        zrhodsp(i,2)=fac*sum2
c
c spin-only magnetization density
        call doubmt1(gmat,szcoeff,gmat1,kmymax,kmymaxp)
        call doubmt1(fmat,szbcoeff,fmat1,kmymax,kmymaxp)
        call submat1(gmat1,fmat1,rmat,kmymax,kmymaxp)
        sum1=(0.d0,0.d0)
        do kmy=1,kmymax
          sum1=sum1+rmat(kmy,kmy)
        end do
        zmagrho(i)=fac*sum1
c
c integrand to orbital polarization energy
        call doubmt1(gmat,lzcoeff,gmat1,kmymax,kmymaxp)
        call doubmt1(fmat,lzbcoeff,fmat1,kmymax,kmymaxp)
        call submat1(gmat1,fmat1,rmat,kmymax,kmymaxp)
        call replms(rlms,rmat,lmax)
        sum1=(0.d0,0.d0)
        sum2=(0.d0,0.d0)
        do lm=5,9
          sum1=sum1+rlms(lm)
          sum2=sum2+rlms(lm+lmmax)
        end do
        cf(i)=fac*(sum1*bopr(i,1)+sum2*bopr(i,2))
c
c next radial point
      end do
c
      x=dlog(rs)-(ns-1)*dx
      do i=1,ns
        r(i)=dexp(x)
        f1(i)=dreal(cf(i))/r(i)
        f2(i)=dimag(cf(i))/r(i)
        x=x+dx
      end do
      zenorb=0.5d0*dcmplx(rsimp(f1,r,ns,dx),rsimp(f2,r,ns,dx))
c
      return
      end
