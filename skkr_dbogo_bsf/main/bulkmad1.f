c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
c     ---------------------------------------------------------
      subroutine bulkmad1(natom,maxnatom,lmax,sigma,vol,rhoatom,
     >                    maxk,veck,distk,maxr,rvec,distr,
     >                    conmad)
c ======================
c
c subroutine to calculate sum of madelung constants
c
c input:
c   natom    - number of atoms per unit cell
c   sigma    - ewald parameter (arbitrary positive constant)
c   lmax     - self-explaining
c
c   rvec,distr,maxr in common/dirvec/ - explanation in routine vector
c   veck,distk,maxk in common/dirvec/ - explanation in routine vector
c
c output:
c   conmad    madelungs constant A^L,L'_mu,nu for L=1,lmmax
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      dimension ffact(0:lshape+1),powm1(0:lshape),pow2(0:lshape)
      dimension rhoatom(3,maxnatom)
      dimension veck(3,mrec),rvec(3,mdir)
      dimension ak(3),jndk(mrec),jshk(mk)
      dimension ar(3),jndr(mdir),jshr(mr)
c
      complex*16 ci,qfacg,qfacr,eikrho,zarg
      complex*16 sum,sumdir,sumrec,gsum(lmshape)
      complex*16 powi(0:lshape),gam(0:lshape+1)
      complex*16 ylm(-lshape:lshape,0:lshape)
c
      complex*16 conmad(lmsup,lmsup,maxnatom,maxnatom)
c
      common /test/itest
      data ci/(0.0d0,1.0d0)/
      data tol/1.0d-12/
c
      lmaxs=2*lmax
      lmaxf=4*lmax
c 
c initialize constants
c
      pi=4.d0*datan(1.d0)
      facg=16.0d0*pi**2/vol
      facr=4.0d0*dsqrt(pi)
      subg=(4.d0*pi)**(3.0d0/2.0d0)*sigma*sigma/vol
c
      ffact(0)=1.0d0
      pow2(0)=1.0d0
      powm1(0)=1.0d0
      powi(0)=(1.0d0,0.0d0)
      do l=1,lmaxf
        ffact(l)=dfloat(2*l-1)*ffact(l-1)
        pow2(l)=2.0d0*pow2(l-1)
        powm1(l)=-powm1(l-1)
        powi(l)=ci*powi(l-1)
      end do
      ffact(lmaxf+1)=dfloat(2*lmaxf+1)*ffact(lmaxf)
c
c sort 3D reciprocal space vectors around ak=0
c
      ak(1)=0.d0
      ak(2)=0.d0
      ak(3)=0.d0
      call vecsrt3d(maxk,mk,veck,ak,distk,jndk,jshk,jmxshk)
c
c sum over atoms in unit cell
      do mu=1,natom
      do nu=1,natom
c
         if(itest.ge.2)
     >   write(6,'(''Atoms:'',2i4)') mu,nu
c
         ar(1)=rhoatom(1,mu)-rhoatom(1,nu)
         ar(2)=rhoatom(2,mu)-rhoatom(2,nu)
         ar(3)=rhoatom(3,mu)-rhoatom(3,nu)
c
c sort 3D direct space vectors around ar
         call vecsrt3d(maxr,mr,rvec,ar,distr,jndr,jshr,jmxshr)
c
c initialize real-space lattice summation:
c leave out first term for rho = zero
c
         rr=dsqrt(ar(1)*ar(1)+ar(2)*ar(2)+ar(3)*ar(3))
         if(rr.lt.1.0d-06) then
           imin=2
           n1=1
         else
           imin=1
           n1=0
         endif
c
c loop over L's
c
         do l=0,lmaxf
c
         qfacg=facg*powi(l)/ffact(l+1)
         qfacr=facr*pow2(l)*powm1(l)/ffact(l+1)
c
         do m=-l,l
            lm=l*(l+1)+m+1
c
            if(itest.ge.5)
     >      write(6,'(''L,M,(L,M):'',3i4)') l,m,lm
c
            sumrec=(0.d0,0.d0)
            sumdir=(0.d0,0.d0)
c
c 1. Summation in reciprocal space
c
c sum over shells (leave out G=0)
c
            n=1
            do i=2,jmxshk
              ns=jshk(i)
              iv=jndk(n+1)
              rsq=veck(1,iv)**2.d0
     &           +veck(2,iv)**2.d0
     &           +veck(3,iv)**2.d0
              rmod=dsqrt(rsq)
              exk=dexp(-rsq*sigma*sigma)
              fac=exk*rmod**dfloat(l-2)
              if (dabs(fac) .lt. tol) goto 1
c
c sum over reciprocal vectors in the shell i
c
              sum=(0.0d0,0.0d0)
              do j=1,ns
                n=n+1
                iv=jndk(n)
                ak(1)=veck(1,iv)
                ak(2)=veck(2,iv)
                ak(3)=veck(3,iv)
                arg=ar(1)*ak(1) + ar(2)*ak(2) + ar(3)*ak(3)
                eikrho=dcmplx(dcos(arg),dsin(arg))
                call spherh(l,ak(1),ak(2),ak(3),ylm)
                sum=sum+eikrho*dconjg(ylm(m,l))
              end do
c
              sumrec=sumrec+fac*sum
              if(itest.ge.5)
     >        write(6,'('' sh,sumg ::'',t30,i6,2d20.12)') i,sumrec
            end do
c
    1       continue
            if(itest.ge.5)
     >      write(6,'('' sumrec ::'',t36,2d20.12)') sumrec
c
            sumrec=sumrec*qfacg
            if(lm.eq.1) sumrec=sumrec-subg
c
c
c 2. Summation in direct space
c
c sum over shells 
c
            n=n1
            imax=jmxshr
            do i=imin,imax
              ns=jshr(i)
              iv=jndr(n+1)
              rsq=(ar(1)+rvec(1,iv))**2.d0
     &           +(ar(2)+rvec(2,iv))**2.d0
     &           +(ar(3)+rvec(3,iv))**2.d0
              rmod=dsqrt(rsq)
              zarg=rsq/(4.0d0*sigma*sigma)
              call incgamma(gam,-l,zarg)
              fac=dreal(gam(l))/rmod**dfloat(l+1)
              if (dabs(fac) .lt. tol) goto 2
c
c sum over real-space vectors in the shell i
c
              sum=(0.d0,0.d0)
              do j=1,ns
                n=n+1
                iv=jndr(n)
                x0=ar(1)+rvec(1,iv)
                y0=ar(2)+rvec(2,iv)
                z0=ar(3)+rvec(3,iv)
                call spherh(l,x0,y0,z0,ylm)
                sum=sum+dconjg(ylm(m,l))
              end do
c
              sumdir=sumdir+fac*sum
              if(itest.ge.5)
     >        write(6,'('' sh,sumr ::'',t30,i6,2d20.12)') i,sumdir
            end do
c
    2       continue
            if(itest.ge.5)
     >      write(6,'('' sumdir ::'',t36,2d20.12)') sumdir
c
            sumdir=sumdir*qfacr
c
            gsum(lm)=sumrec+sumdir
            if(mu.eq.nu.and.lm.eq.1) gsum(lm)=gsum(lm)-2.0d0/sigma
c
            if(itest.ge.2.and.cdabs(gsum(lm)).gt.tol)
     >      write(6,'(3i4,2x,2d15.7,2x,2d15.7,2x,2d15.7)') 
     >      l,m,lm,sumrec,sumdir,gsum(lm)
c
         end do
         end do
c end of loop over L's
c
c Put together Madelung constants
c
         do l=0,lmaxs
         do m=-l,l
           lm=l*l+l+m+1
         do lp=0,lmaxs
         do mp=-lp,lp
           lmp=lp*lp+lp+mp+1
c Laszlo's convention
           ilm=(l+lp)*(l+lp+1)+(mp-m)+1
c Jan's convention
c          ilm=(l+lp)*(l+lp+1)+m+mp+1
c
c Laszlo's convention
           conmad(lm,lmp,mu,nu)=2.0d0*powm1(lp)*ffact(l+lp+1)*
     >     gaunt(lp,l,l+lp,mp,m,mp-m)*gsum(ilm)/(ffact(l)*ffact(lp))
c Jan's convention
c          conmad(lm,lmp,mu,nu)=powm1(lp)*ffact(l+lp+1)*
c    >     gaunt(l+lp,l,lp,m+mp,m,mp)*gsum(ilm)/(ffact(l)*ffact(lp))
c
         end do
         end do
         end do
         end do
c
      end do
      end do
c end of loop over atoms
c
      return
      end
