c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine magnet(lmax,rs,dx,ns,gz,fz,gj,fj,
     >                  tau,gcoeff,gbcoeff,zmom,lms)    
c======================
c
c input:  lmax  - maximum of angular momentum index
c         rs,dx,ns - as in 'readpot'
c         gz,fz - large and small component of regular radial solution * r
c         gj,fj - large and small component of irreg. radial solution * r
c         gcoeff, gbcoeff - angular momentum matrixelements
c         tau - site--diagonal tau matrix
c         lms - if .true. output in (l,m,s) representation
c output: zmom
c         (kappa,my)-like or (l,m,s)-like resolution of the 
c         magnetization density corresponding to rs.
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      logical lms
      dimension rrr(nrad),rir(nrad)
      dimension rri(nrad),rii(nrad)
      dimension rad(nrad)
c
      complex*16 gz(kmymaxp,kmymaxp,nrad)
      complex*16 fz(kmymaxp,kmymaxp,nrad)
      complex*16 gj(kmymaxp,kmymaxp,nrad)
      complex*16 fj(kmymaxp,kmymaxp,nrad)
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
        rrr(1:ns)=0.d0
        rir(1:ns)=0.d0
        rri(1:ns)=0.d0
        rii(1:ns)=0.d0
        do kmy1=1,kmymax
        do kmyp1=1,kmymax
          gcox=fac*gcoeff(kmy1,kmyp1)
          gbcox=fac*gbcoeff(kmy1,kmyp1)
c
          if(cdabs(gcox).gt.small) then
          do i=1,ns
            cxr=gz(kmy1,kmy,i)*gz(kmyp1,kmyp,i)*gcox
            cxi=gz(kmy1,kmy,i)*gj(kmyp1,kmyp,i)*gcox
            rrr(i)=rrr(i)+dreal(cxr)
            rir(i)=rir(i)+dimag(cxr)
            rri(i)=rri(i)+dreal(cxi)
            rii(i)=rii(i)+dimag(cxi)
          end do
          end if
          if(cdabs(gbcox).gt.small) then
          do i=1,ns
            cxr=fz(kmy1,kmy,i)*fz(kmyp1,kmyp,i)*gbcox
            cxi=fz(kmy1,kmy,i)*fj(kmyp1,kmyp,i)*gbcox
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
      call submat(gf,sum(1,1,2),kmymax,kmymaxp)
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
