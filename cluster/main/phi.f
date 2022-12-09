c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine phi(rho,sigma,fi,psi,xi)
c ==================
c input:  rho - vector argument
c         sigma -  Ewald parameter
c output: fi,psi,xi
c
      implicit real*8 (a-h,o-z)
c
      character*4 idgroup
      dimension rho(3),rhopar(2)
      dimension fi1(5),psi1(5),xi1(5),x(5)
c
      common/latt2d/idgroup,ilat,alat,beta,delta
c
      data tiny/1.0d-8/
c
      small=0.05
      if(ilat.le.3) small=0.1 
c
      rhopar(1)=rho(1)
      rhopar(2)=rho(2)
      rhoperp=rho(3)
      dpar=rhopar(1)*rhopar(1)+rhopar(2)*rhopar(2)
      dpar=dsqrt(dpar)
c
      if(dabs(rhoperp).lt.small) then
c
        call phipar(rhopar,sigma,fi)
c
        if(dabs(rhoperp).gt.tiny) then
c
          x(1)=0.d0
          fi1(1)=fi
          do i=2,5
            x(i)=(1.0d0+0.1d0*(i-1))*small
            if(rhoperp.lt.0.d0) x(i)=-x(i)
            call phiperp(rhopar,x(i),fi1(i),psi1(i),xi1(i))
            d=dpar*dpar+x(i)*x(i)
            d=dsqrt(d)
            if(dpar.lt.tiny) fi1(i)=fi1(i)-2.0d0/d
          end do
          fi=ylag(rhoperp,x,fi1,0,4,5,iex)
          if(dpar.lt.tiny) then
            d=dpar*dpar+rhoperp*rhoperp
            d=dsqrt(d)
            fi=fi+2.0d0/d
          end if
          psi=ylag(rhoperp,x(2),psi1(2),0,3,4,iex)
          xi=ylag(rhoperp,x(2),xi1(2),0,3,4,iex)
c
        else
          psi=0.0d0
          xi=0.0d0
        end if
c
      else
c
        call phiperp(rhopar,rhoperp,fi,psi,xi)
c
      end if
c
      return
      end
      subroutine phiperp(rho,rhoperp,fi,psi,xi)
c ======================
c
c "out-of-layer" terms to Madelung constant as defined in
c J. Kudrnovsky et al., PRB ... (1993), (A1) and (A2)
c
c input:  rho - paralell components of vector argument
c         rhoperp -  it's perpendicular component
c output: fi, psi
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      dimension rho(2)
      common/recvec2d/veck(2,mrec),distk,maxk
      common/vindk2d/indk(mrec),nshk(mk),maxshk
      common/crystl2d/vol2d,volbz2d
      data tol/1.d-10/
c
      pi=4.d0*datan(1.d0)
      fac=4.d0*pi/vol2d
      rhop=dabs(rhoperp)
c
c summation in reciprocal space
c leave out first term
c
      sum1=0.d0
      sum2=0.d0
      sum3=0.d0
      i=1
      do ish=2,maxshk-1
        ind=indk(i+1)
        d=dsqrt(veck(1,ind)*veck(1,ind)+veck(2,ind)*veck(2,ind))
        arg=d*rhop
        if(arg.gt.25.) goto 2
        term=dexp(-arg)
        sumin=0.d0
        nsh=nshk(ish)
        do insh=1,nsh
          i=i+1
          ind=indk(i)
          rg=rho(1)*veck(1,ind)+rho(2)*veck(2,ind)
          sumin=sumin+dcos(rg)
        end do
        sum1=sum1+term*sumin/d
        sum2=sum2+term*sumin
        sum3=sum3+term*sumin*d
c       write(6,'(i6,3e20.10)') ish,d,arg,term
        if(dabs(term).lt.tol) goto 2
      end do
      write(6,'(/i6,3e20.10)') ish,d,fac*term*sumin/d,fac*sum1
      write(6,'('' WARNING!!! PHIPERP: ''/'' reciprocal lattice '',
     >''sum did not converge up to shell no.'',i4/)') ish
      stop
    2 continue
c
      fi=fac*sum1
      psi=fac*sum2
      xi=fac*sum3
      if(rhoperp.lt.0.d0) psi=-psi
c
      return 
      end
      subroutine phipar(rho,sigma,fi)
c =====================
c
c "in-layer" component of Madelung constant
c generalization of (A3) of J. Kudrnovsky et al., PRB ... (1993)
c
c input:  rho - vector argument, paralell to the layer
c         sigma - Ewald parameter
c output  fi
c
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
      dimension rho(2)
      common/dirvec2d/rvec(2,mdir),distr,maxr
      common/vindr2d/indr(mdir),nshr(mr),maxshr
      common/recvec2d/veck(2,mrec),distk,maxk
      common/vindk2d/indk(mrec),nshk(mk),maxshk
      common/crystl2d/vol2d,volbz2d
      data tol/1.d-10/,tiny/1.0d-6/
c
      pi=4.d0*datan(1.d0)
      fac=4.0d0*pi/vol2d
      piroot=dsqrt(pi)
c
c summation in direct space
c
c sort real space vectors around rho
      call vecsrt2d(maxr,mr,rvec,rho,distr,indr,nshr,maxshr)
c
      j=indr(1)
      d1=dsqrt((rho(1)+rvec(1,j))**2+(rho(2)+rvec(2,j))**2)
      if(d1.lt.tiny) then
        i0=1
        ish0=2
        rcorr=-2.d0/sigma/piroot
      else
        i0=0
        ish0=1
        rcorr=0.0d0
      end if
      sumr=rcorr
      i=i0
      do ish=ish0,maxshr-1
        nsh=nshr(ish)
        ind=indr(i+1)
        d=dsqrt((rho(1)+rvec(1,ind))**2+(rho(2)+rvec(2,ind))**2)
        arg=0.5d0*d/sigma
        if(arg.gt.25.d0) goto 1
        termr=2.0d0*nsh*(1.d0-erf(arg))/d
        sumr=sumr+termr
        if(dabs(termr).lt.tol) goto 1
        i=i+nsh
      end do
      write(6,'(i6,3e20.10)') ish,d,termr,sumr
      write(6,'(/'' WARNING!!! PHIPAR: ''/'' direct lattice '',
     >''sum did not converge up to shell no.'',i4/)') ish
      stop
    1 continue
c
c summation in reciprocal space
c leave out first term
c
      sumk=0.d0
      i=1
      do ish=2,maxshk-1
        ind=indk(i+1)
        d=dsqrt(veck(1,ind)*veck(1,ind)+veck(2,ind)*veck(2,ind))
        arg=d*sigma
        if(arg.gt.10.d0) goto 2
        termk=(1.d0-erf(arg))/d
        sumink=0.d0
        nsh=nshk(ish)
        do insh=1,nsh
          i=i+1
          ind=indk(i)
          rg=rho(1)*veck(1,ind)+rho(2)*veck(2,ind)
          sumink=sumink+dcos(rg)
        end do
        sumk=sumk+termk*sumink
        if(dabs(termk).lt.tol) goto 2
      end do
      write(6,'(i6,3e20.10)') ish,d,fac*termk*sumink,fac*sumk
      write(6,'(/'' WARNING!!! PHIPAR: ''/'' reciprocal lattice '',
     >''sum did not converge up to shell no.'',i4/)') ish
      stop
    2 continue
c
      fi=sumr-8.d0*sigma*piroot/vol2d+fac*sumk
c
      return 
      end
