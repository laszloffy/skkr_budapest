c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine wafu(e,lmax,idpot,v0,vr,dx,nmt,rmt,iflag,rel,
     >                tminv,fmt,gmt,fpint,yl1,yl2)
c====================
c
c input:  e - any complex energy
c         lmax  - maximum of angular momentum index
c         v0 - vacuum potential
c         idpot,vr,nmt,rmt - as in 'readpot'
c         iflag - if 0 compute tminv only, if 1 fint,gint as well
c         rel   - if 'false' non-relativistic, if 'true' scalar-rel.
c output: tminv - inverse of single-site scattering matrix
c         fmt  - radial integral of the regular solution upto rmt
c         gmt  - radial integral of the irregular solution upto rmt
c         fpint - radial integrals for dipole moments
c         yl1   - regular radial solution
c         yl2   - irregular radial solution 
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      logical rel
      character*10 idpot
c
      dimension vr(nrad),rs(nrad)
      dimension re1(nrad),ri1(nrad),re2(nrad),ri2(nrad)
c
      complex*16 tminv(0:lmaxp)
      complex*16 fmt(0:lmaxp),gmt(0:lmaxp),fpint(lmaxp)
      complex*16 yl1(nrad,0:lmaxp),yl2(nrad,0:lmaxp)
      complex*16 yl1p(nrad),yl2p(nrad)
      complex*16 fb(0:l1maxp),fn(0:l1maxp),fh(0:l1maxp)
      complex*16 fmb(0:l1maxp),fmn(0:l1maxp),fmh(0:l1maxp)
      complex*16 fmodb(0:l1maxp),fmodn(0:l1maxp),fmodh(0:l1maxp)
      complex*16 e,p,tand,react,xldl,sqrm1
      complex*16 ynorm,y1rm,y1prm,y2rm,y2prm,y1,y2
      complex*16 alpha,beta,delta,acoeff,bcoeff
      complex*16 wronsky,br,bmt,evac,pvac
c
      data sqrm1/(0.d0,1.d0)/
c
c---> c in rydberg units:
      clight=274.072d0
      csqinv=1.d0/(clight*clight)
      p=cdsqrt(e)
      if(dimag(p).lt.0.d0) p=-p
      if(.not.rel) csqinv=0.d0
c
      if(idpot.eq.'Vacuum    ') then
        evac=dcmplx(e-v0)
        pvac=cdsqrt(evac)
        if(dimag(pvac).lt.0.d0) pvac=-pvac
        call csbf(lmax+1,l1maxp,p,rmt,fmb,fmn,fmh)
        call csbf(lmax+1,l1maxp,pvac,rmt,fmodb,fmodn,fmodh)
      else
        call csbf(lmax+1,l1maxp,p,rmt,fmb,fmn,fmh)
      end if
c
      xmt=dlog(rmt)
      x0=xmt-(nmt-1)*dx
      x=x0
      do i=1,nmt
        rs(i)=dexp(x)
        x=x+dx
      end do
c
      do 1 l=0,lmax
        fl=dfloat(l)
        if(idpot.eq.'Vacuum    ') then
          tand=(p*fmb(l+1)*fmodb(l)-pvac*fmb(l)*fmodb(l+1))/
     >         (p*fmn(l+1)*fmodb(l)-pvac*fmn(l)*fmodb(l+1))
        else
          call outwrd(rel,e,l,vr,dx,x0,rmt,nmt,yl1(1,l),yl1p,xldl)
          tand=(p*fmb(l+1)+(xldl-fl/rmt)*fmb(l))/
     >         (p*fmn(l+1)+(xldl-fl/rmt)*fmn(l))
        endif
        react=-tand/p
        tminv(l)=1.d0/react+sqrm1*p
        if(iflag.eq.0) goto 1
c
        call inwrd(rel,e,l,vr,dx,rmt,nmt,yl2(1,l),yl2p)
c
        y1rm=yl1(nmt,l)
        y1prm=yl1p(nmt)
        y2rm=yl2(nmt,l)
        y2prm=yl2p(nmt)
        delta=y1rm*y2prm-y2rm*y1prm
        alpha=rmt*fmb(l)
        beta=rmt*((fl+1.d0)*fmb(l)-p*rmt*fmb(l+1))
        acoeff= (alpha*y2prm-beta*y2rm)/delta
        bcoeff=-(alpha*y1prm-beta*y1rm)/delta
c
        y1rm=(p*fmn(l)+fmb(l)/react)*rmt
        ynorm=y1rm/yl1(nmt,l)
c
        do i=1,nmt
          yl2(i,l)=acoeff*yl1(i,l)+bcoeff*yl2(i,l)
          yl2p(i) =acoeff*yl1p(i)+bcoeff*yl2p(i)
          yl1(i,l)=yl1(i,l)*ynorm
          yl1p(i)=yl1p(i)*ynorm
          y1=yl1(i,l)*yl1(i,l)
          y2=yl1(i,l)*yl2(i,l)
          re1(i)=dreal(y1)
          ri1(i)=dimag(y1)
          re2(i)=dreal(y2)
          ri2(i)=dimag(y2)
        end do
c
        r1r=rsimp(re1,rs,nmt,dx)
        r1i=rsimp(ri1,rs,nmt,dx)
        r2r=rsimp(re2,rs,nmt,dx)
        r2i=rsimp(ri2,rs,nmt,dx)
        fmt(l)=dcmplx(r1r,r1i)
        gmt(l)=dcmplx(r2r,r2i)
c
    1 continue
c
      do l1=1,lmax
        do i=1,nmt
          y1=yl1(i,l1-1)*yl1(i,l1)*rs(i)
          re1(i)=dreal(y1)
          ri1(i)=dimag(y1)
        end do
        r1r=rsimp(re1,rs,nmt,dx)
        r1i=rsimp(ri1,rs,nmt,dx)
        fpint(l1)=dcmplx(r1r,r1i)
      end do 
c
  100 format(i3,2x,f10.8,2x,6(d15.8))
      return
      end
