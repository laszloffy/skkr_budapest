c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
c Build up kambe structure constant
c=======================================================================
c Input:
c      e    : energy
c     kx,ky : k-point
c     g2    : reciprocal 2D-lattice vectors
c     lmax  : maximal angular momentum
c     lmmax : (lmax+1)^2
c     ng2   : number of reciprocal lattice vectors
c     nr2   : number of lattice vectors
c     r2    : 2D-lattice vectors
c     v2    : volume of 2D-cell   (common/crystl2d/)
c     ts    : separation of the 2 layers
c
c Output: 
c    alm    : Kambe's structure constant
c=======================================================================
c
      subroutine kamstr(alm,e,kx,ky,lmax,lmmax,ngsh2,gsh2,indg,g2,nrsh2, 
     &     rsh2,indr,r2,ts,gau,factor)
c
      implicit none              
      integer lmax,lmmax,indg(*),indr(*),nrsh2,rsh2(*),ngsh2,gsh2(*),
     &        itest
      double precision kx,ky,g2(2,*),r2(2,*),v2,volbz,ts(0:3),gau(*),
     &                 factor,ctol
      complex*16 alm(lmmax,lmmax),e
c
c local
      integer i,ig,il,ir,is,iis,j,jl,js,k,l,m,mm,lm,n,lmx,lmmx,lmax2,   
     &     lm1,lm2,l1,l2,m1,m2,ipar,o,ll,no,nf,nl,lmmax2
      parameter (lmx=12, lmmx=(lmx+1)*(lmx+1))
      double precision d,dc,dl,pi,phi,gkx,gky,sx,sy,sz,r,rk,rkx,rky,rz,
     &     rx,ry,km,rs,rsx,rsy,rsz,fac(0:lmx*2),flms(lmmx),s,ds,gx,gy,
     &     dmt,flm(lmmx),tin,tst
      parameter (tin=1.d-10)
!
! Modified by Laszlo Szunyogh on 4 July 2002
      parameter (ctol=1.d-4)
!
      parameter (pi=3.1415926535897932384626d0) 
      complex*16 a,b,c0,c1,ci,c,xc,zc,kp,del(0:lmx),ex,ex1,ex2,csum, 
     &     kappa,cs1,cs2,cm(0:lmx),cp(0:lmx),cs(0:lmx),alpha,       
     &     gam(0:lmx),ylm(169),intil(0:lmx),cx,cy,cz,d003,dlm1(lmmx),
c !!!&     gam(0:lshape+1),ylm(lmmx),intil(0:lmx),cx,cy,cz,d003,dlm1(lmmx),
     &     dlm2(lmmx),w
      parameter (c0=(0.d0,0.d0),c1=(1.d0,0.d0),ci=(0.d0,1.d0))
      logical f35,f40,f47,f50,f52
c
c external
      complex*16 cerf
      external cerf,incgamma,ilkam
c
      common/crystl2d/v2,volbz
      common/test/itest
c
c=======================================================================
      lmax2=lmax*2
      lmmax2=(lmax2+1)**2
      kappa=sqrt(e)
      cs1=-c1/v2/kappa
      cs2=-kappa/2.d0/sqrt(pi)
      alpha=factor*e*v2/4.d0/pi
      w=alpha*2.d0/e
!     write(6,'(a,2d18.9)')'kappa',kappa
!     write(6,'(a,8d18.9)')'kappa,cs1,cs2,alpha',kappa,cs1,cs2,alpha
!
!=========================== factorial of lmax*2 =======================
      fac(0)=1.d0
      do i=1,lmax2+lmax2
         fac(i)=fac(i-1)*float(i)
      enddo
!====================================== F(lm) ======= appr. (5.29)-Ch ==
      d=2.d0
      do l=0,lmax2
         d=d/2.d0
         dl=2.d0*l+1.d0
         do m=0,l
            a=d*dsqrt(dl*fac(l+m)*fac(l-m))
            lm=l*(l+1)+m+1
            flm(lm)=a
            lm=l*(l+1)-m+1
            flm(lm)=a
         enddo
      enddo
!                                       _
!====================================== F(lm) ============= (5.47)-Ch ==
      d=4.d0
      do l=0,lmax2
         d=d/4.d0
         dl=2.d0*l+1.d0
         do m=0,l
            if(mod(l-m,2).eq.0) then
               dc=(-1.d0)**l*(-1.d0)**((l-m)/2)*d
               k=(l-m)/2
               n=(l+m)/2
               i=l+m
               j=l-m
               r=dc*dsqrt(dl*fac(i)*fac(j))/fac(k)/fac(n)
               lm=l*(l+1)+m+1
               flms(lm)=r
               lm=l*(l+1)-m+1
               flms(lm)=r
            else
               lm=l*(l+1)+m+1
               flms(lm)=0.d0
               lm=l*(l+1)-m+1
               flms(lm)=0.d0
            endif
         enddo
      enddo
!
!================================= Initialize ==========================
      do i=1,lmmax
         do j=1,lmmax
            alm(i,j)=c0
         enddo
      enddo
!===================================== formula 55 ========= (5.48)-Ch ==
      a=sqrt(alpha)
      c=ci*(exp(-alpha)-cerf(a))-c1/sqrt(pi)/a
      d003=-kappa*(exp(alpha)*c-ci)/2.d0/dsqrt(pi)
!=============================== end of formula 55 =====================
      sx= ts(1)
      sy= ts(2)
      sz= ts(3)
      s = ts(0)
c
!================== Initialize dlm's  ==================================
      do lm=1,lmmx
         dlm1(lm)=c0
         dlm2(lm)=c0
      enddo
      if(s.lt.ctol) dlm1(1)=dlm1(1)+d003
!
! Modified by Laszlo Szunyogh on 4 July 2002
      if(s.gt.ctol.and.abs(sz).lt.ctol) then
        sz=1.1d0*ctol
        s=sqrt(sx*sx+sy*sy+sz*sz)
      end if
!
!================== Select conditions according to cvec ================
      f35=.false.
      if(s.gt.ctol.and.abs(sz).gt.ctol) f35=.true.
      f40=.false.
      if(s.gt.ctol.and.abs(sz).lt.ctol) f40=.true.
      f47=.false.
      if(s.gt.ctol) f47=.true.
      f50=.false.
      f52=.false.
      if(s.lt.ctol) then 
         f50=.true.
         f52=.true.
      endif
!
!======================== (-kappa*sz)**(2n-s) ==========================
      c=kappa*sz
!     a=c1/c
      a=c1
      cs(0)=a
      do l=1,lmax2
         a=a*c
         cs(l)=a
c        write(6,*)'l,cs(l)',l,cs(l)
      enddo
!=======================================================================
      no=0
      tst=1.d+30
      do is=1,ngsh2 
         do iis=1,gsh2(is)
            no=no+1
            ig=indg(no)
            gx=g2(1,ig)
            gy=g2(2,ig)
!----------------------
            gkx=gx+kx
            gky=gy+ky
            km=gkx*gkx+gky*gky
            kp=sqrt(e-km)
            km=sqrt(km)
            ex2=1.0d0
            if (km.gt.1.d-10) ex2=dcmplx(gkx,gky)/km
            xc=exp(-ci*pi)*(e-km*km)*w/2.d0
            call incgamma(gam,lmax2,xc)
!============================== formula 35,40 ==========================
!======================== (kp/kappa)**(2n-1) ===========================
            c=kp/kappa
            c=c*c
            a=(kp/kappa)**(-3)
            do l=0,lmax2
               a=a*c
               cp(l)=a
            enddo
!======================== (km/kappa)**(l-s) ============================
            c=km/kappa
            a=c1/c
            do l=0,lmax2
               a=a*c
               cm(l)=c**l       !a
            enddo
!=======================================================================
!=================================== delta =============================
            if(f40.or.f35) then
               zc=kp*sz
               zc=zc*zc
               if(f35) call delkam(del,lmax2,xc,zc)
!=======================================================================
               ex1=exp(ci*(gkx*sx+gky*sy))*cs1   
               do l=0,lmax2
                  ex=ex2**(l+1)*ex1
                  do m=-l,l
                     ex=ex/ex2
                     lm=l*(l+1)+m+1
                     mm=abs(m)
                     if(f40) goto 11
!============================ formula 35 ================= (5.28)-Ch ===
                     csum=c0
                     do n=0,l-mm
                        a=cp(n)*del(n)
                        b=c0
                        ipar=mod(l-mm,2)
                        do k=n+mod(n+ipar,2),min(n+n,l-mm),2
                           i=n+n-k
                           j=l-k
                           c=cs(i)*cm(j)/fac(i)/fac(k-n)
                           i=(l-mm-k)/2
                           j=(l+mm-k)/2
                           b=b+c/fac(i)/fac(j)
                        enddo
                        csum=csum+a*b
                     enddo
                     a=ex*flm(lm)*csum*ci**(1-m)
                     goto 20 
!================================ or formula 40 ========== (5.33)-Ch ===
 11                  continue
                     csum=c0
                     do n=0,(l-mm)/2
                        i=l-2*n
                        a=cp(n)*cm(i)*gam(n)/fac(n)
                        i=(l-mm-2*n)/2
                        j=(l+mm-2*n)/2
                        csum=csum+a/fac(i)/fac(j)
                     enddo
                     a=ex*csum*flm(lm)
!=======================================================================
 20                  dlm1(lm)=dlm1(lm)+a
                  enddo
               enddo
            elseif(f50) then
!===================================== formula 50 ========= (5.43)-Ch ==
               do l=0,lmax2
                  ex=ex2**(l+1)*cs1
                  do m=-l,l
                     ex=ex/ex2
                     lm=l*(l+1)+m+1
                     mm=abs(m)
                     csum=c0
                     if(mod(l-mm,2).eq.0) then
                        do n=0,(l-mm)/2
                           i=l-2*n
                           a=cp(n)*cm(i)/fac(n)*gam(n)
                           i=(l-mm-2*n)/2
                           j=(l+mm-2*n)/2
                           csum=csum+a/fac(i)/fac(j)
                        enddo
                     endif
                     a=ex*csum*flm(lm)*ci**(1-m)
                     dlm1(lm)=dlm1(lm)+a
                  enddo
               enddo
            endif
         enddo   !iis
!=========================== test of convergence =======================
         d=0.d0
         do lm=1,lmmax2
            d=d+abs(dlm1(lm))
         enddo

         if(dabs(tst-d).lt.tin.and.is.ge.3) goto 800
         tst=d
      enddo   !is
      write(6,*)'g-loop not converged!'
      STOP
 800  continue
      if(itest.gt.5) write(6,*)'numbers of g-vectors:',no
!=============================== end of formula 50 =====================
!=======================================================================
!================================== formula 47,52 ======================
      no=0
      tst=1.d+30
      do is=1,nrsh2 
         do iis=1,rsh2(is)
            no=no+1
            ir=indr(no)
            rx=r2(1,ir)
            ry=r2(2,ir)
            r=sqrt(rx*rx+ry*ry)
            rkx=rx*kx
            rky=ry*ky
            rk=rkx+rky
            ex1=exp(-ci*rk)
            ex2=1.0d0
            if (r.gt.1.d-10) ex2=dcmplx(rx,ry)/r
            f47=.true.
            if(f47) then
!=======================================================================
!================================== formula 47 =========== (5.40)-Ch ===
               rsx=rx+sx
               rsy=ry+sy
               rsz=sz
               rs=sqrt(rsx*rsx+rsy*rsy+rsz*rsz)
               if(rs.lt.1.d-5) goto 70
               call ylmc(rsx,rsy,rsz,lmax2,ylm)
               call ilkam(intil,lmax2,alpha,rs,kappa)
               a=-kappa*rs/2.d0
               b=ex1*cs2/a
               do l=0,lmax2
                  b=b*a
                  c=intil(l)*b
                  do m=-l,l
                     lm=l*(l+1)+m+1
                     dlm2(lm)=dlm2(lm)+c*conjg(ylm(lm))
                  enddo
               enddo
 70            continue
!=============================== end of formula 47 =====================
!================================== formula 52 ============ (5.45)-Ch ==
            f52=.false.
            elseif(f52) then
               if(r.lt.1.d-4) goto 30
               call ilkam(intil,lmax2,alpha,r,kappa)
               a=kappa*r
               b=ex1*cs2/a/2.d0/sqrt(pi)
               do l=0,lmax2
                  b=b*a                     
                  ex=ex2**(l+1)*b*intil(l)
                  do m=-l,l
                     lm=l*(l+1)+m+1
                     ex=ex/ex2
                     c=ex*flms(lm)
                     dlm2(lm)=dlm2(lm)+c
                  enddo
               enddo
 30            continue
            endif
         enddo    !iis
!========================== test of convergence ========================
         d=0.d0
         do lm=1,lmmax2
            d=d+abs(dlm2(lm))
         enddo
         if(dabs(tst-d).lt.tin.and.is.ge.3) goto 801
         tst=d
      enddo    !is
      write(6,*)'r-loop not converged!'
      STOP
 801  continue
      if(itest.gt.5) write(6,*)'numbers of r-vectors:',no
!=======================================================================
!     do lm=1,lmmax2
!        l=ll(lm)
!        dlm1(lm)=dlm1(lm)*(-1)**l
!        dlm2(lm)=dlm2(lm)*(-1)**l
!     enddo
      o=0
!     write(6,'(3f12.6)')ts(1),ts(2),ts(3)
      do l1=0,lmax
         do m1=-l1,l1
            lm1=l1*(l1+1)+m1+1
            do l2=0,lmax
               do m2=-l2,l2
                  lm2=l2*(l2+1)+m2+1
                  csum=c0
                  do l=0,lmax2
                     do m=-l,l
                        o=o+1
                        lm=l*(l+1)+m+1
                        csum=csum+gau(o)*(dlm1(lm)+dlm2(lm))
                     enddo
                  enddo
                  alm(lm1,lm2)=csum
!                 write(6,'(2d18.9,2i5)')alm(lm1,lm2),lm1,lm2
               enddo
            enddo
         enddo
      enddo
!     write(6,*)
!     write(6,'(''sep:'',3f10.5,5x,l1,1x,l1,1x,l1,2x,l1,1x,l1))') 
!    <(ts(i),i=1,3),f35,f40,f50,f47,f52
      return
      end
