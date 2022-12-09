c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine madelungi(qmoml,qmomr,qmomlay,nintfc,lmax,sigma0,
     >                     vleft,vright,vmad)
c=========================
c
c Calculate layer dependent Madelung potentials in monopole
c and dipole approximation for an interlayer 
c
c input:  
c         qmoml - moments of charge densities in the left bulk
c         qmomr - moments of charge densities in the right bulk
c         qmomlay - moments of charge densities in the central region
c         nintfc - number of interfacial layers
c         lmax - l-cutoff                        
c         sigma - Ewald parameter
c         vleft - Madelung potential for the leftmost layer as
c                 calculated by 'madelungb' previously
c         vright - Madelung potential for the rightmost layer as
c                 calculated by 'madelungb' previously
c     in common/lay2d/
c         clay - layer generating vector
c         nbulkl - number of layers in a bulk unit to the left
c         nbulkr - number of layers in a bulk unit to the right
c     in common/crystl2d/
c         vol2d - volume of the 2D unit cell 
c output: vmad - Madelung potentials for the interfacial layers
c
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
      dimension cvec(3),clayb(3)
      dimension qlay(mintfc),dlay(mintfc),vmad(mintfc)
      dimension qbulkl(mbulk),dbulkl(mbulk),cbulkl(3,mbulk),cbulkl0(3)
      dimension qbulkr(mbulk),dbulkr(mbulk),cbulkr(3,mbulk),cbulkr0(3)
c
      complex*16 qmoml(lmsup,minprc)
      complex*16 qmomr(lmsup,minprc)
      complex*16 qmomlay(lmsup,mintfc)
c
      common/lay2d/clay(mtotal,3),nextra,nbulkl,nbulkr,
     &             nprc,ninprc(0:mprc+1) 
      common/crystl2d/vol2d,volbz2d
      common/test/itest
      common/madelung_energy/emad
c
      data tol/1.0d-8/,ctol/1.0d-6/,kmax/20/
c
      emad=0.0d0
      pi=4.d0*datan(1.d0)
      fac=4.d0*pi/vol2d
      t1=dsqrt(3.d0)
      scale=dsqrt(vol2d)
      sigma=sigma0*scale
      lmaxs=2*lmax
      lmmaxs=(lmaxs+1)*(lmaxs+1)
c Index of the uppermost left bulk layer (not part of the I-region)
      nshiftl=(nextra+1)*ninprc(0)
c Index of the uppermost layer in the intermediate region
      nshiftr=(nextra+1)*ninprc(0)+nintfc
c
      nunitsl=ninprc(0)/nbulkl
      nunitsr=ninprc(nprc+1)/nbulkr
c
      do li=1,nintfc
        qlay(li)=qmomlay(1,li)
        dlay(li)=qmomlay(3,li)*t1
      end do
c
c Assign moments and non-primitive translation vectors for the bulk
c
      do ibulk=1,nbulkl
        qbulkl(ibulk)=0.d0
        dbulkl(ibulk)=0.d0
        do iunit=1,nunitsl
          ilay=(iunit-1)*nbulkl+ibulk
          qbulkl(ibulk)=qbulkl(ibulk)+qmoml(1,ilay)
          dbulkl(ibulk)=dbulkl(ibulk)+qmoml(3,ilay)*t1
        end do
        qbulkl(ibulk)=qbulkl(ibulk)/nunitsl
        dbulkl(ibulk)=dbulkl(ibulk)/nunitsl
        do i=1,3
            cbulkl(i,ibulk)=clay(nshiftl-nbulkl+ibulk,i)-
     >                      clay(nshiftl-nbulkl+1,i)
        end do
      end do
      do ibulk=1,nbulkr
        qbulkr(ibulk)=0.d0
        dbulkr(ibulk)=0.d0
        do iunit=1,nunitsr
          ilay=(iunit-1)*nbulkr+ibulk
          qbulkr(ibulk)=qbulkr(ibulk)+qmomr(1,ilay)
          dbulkr(ibulk)=dbulkr(ibulk)+qmomr(3,ilay)*t1
        end do
        qbulkr(ibulk)=qbulkr(ibulk)/nunitsr
        dbulkr(ibulk)=dbulkr(ibulk)/nunitsr
        do i=1,3
            cbulkr(i,ibulk)=clay(nshiftr+ibulk,i)-
     >                      clay(nshiftr+1,i)
        end do
      end do
c
c Generating vectors of the bulk units
c
      do i=1,3
        cbulkl0(i)=clay(nshiftl,i)-clay(nshiftl-nbulkl,i)
        cbulkr0(i)=clay(nshiftr+nbulkr+1,i)-clay(nshiftr+1,i)
      end do
c
c --------------------------------------------------------
c Calculate normalization terms for the Madelung constants
c
      v0=0.0d0
      vnp1=0.0d0
c
      do jbulk=1,nbulkl
c       write(6,'(/'' jbulk='',i3)') jbulk
c Loop over inequivalent sublattices in the bulk to the left
c
        sum1=0.d0
        sum2=0.d0
        do k=0,-kmax,-1
c Loop over bulk units downwards
c
c Separation vector between layers no. 0 and (jbulk,k)
c Note, that the (artificial) left bulk position 
c vector corresponding to index (nbulkl,0) is equivalent
c to clay(nshiftl,*)!
c
          do i=1,3
            clayb(i)=clay(nshiftl,i)-cbulkl(i,nbulkl)+
     >               cbulkl(i,jbulk)+k*cbulkl0(i)
            cvec(i)=clay(nshiftl,i)-clayb(i)
          end do
c
c         ------------------------------
          call phi(cvec,sigma,fi,psi,xi)
c         ------------------------------
c         write(6,'(i3,3f10.5,2d15.7)') k,(cvec(i),i=1,3),fi,psi    
c
          sum1=sum1+fi
          sum2=sum2+psi
          if(dabs(fi).lt.tol) goto 10
c
        end do
        write(6,'(/'' WARNING!!! Madelungi: V(0) ''/'' layer '',
     >    ''left sum did not converge up to bulk unit no.'',i4/)') k
   10   v0=v0+sum1*qbulkl(jbulk)+sum2*dbulkl(jbulk)
c
        sum1=0.d0
        sum2=0.d0
        do k=0,-kmax,-1
c Loop over bulk units downwards
c
c Separation vector between layers no. nintfc+1 and (jbulk,k)
c
          do i=1,3
            clayb(i)=clay(nshiftl,i)-cbulkl(i,nbulkl)+
     >               cbulkl(i,jbulk)+k*cbulkl0(i)
            cvec(i)=clay(nshiftr+1,i)-clayb(i)
          end do
c
c         ------------------------------
          call phi(cvec,sigma,fi,psi,xi)
c         ------------------------------
c         write(6,'(i3,3f10.5,2d15.7)') k,(cvec(i),i=1,3),fi,psi    
c
          sum1=sum1+fi  
          sum2=sum2+psi 
          if(dabs(fi).lt.tol) goto 20
c
        end do
        write(6,'(/'' WARNING!!! Madelungi: V(n+1)''/'' layer '',
     >    ''left sum did not converge up to bulk unit no.'',i4/)') k
   20   vnp1=vnp1+sum1*qbulkl(jbulk)+sum2*dbulkl(jbulk)
c
      end do
c
      do jbulk=1,nbulkr
c       write(6,'(/'' jbulk='',i3)') jbulk
c Loop over inequivalent sublattices in the bulk to the right
c
        sum1=0.d0
        sum2=0.d0
        do k=0,kmax
c Loop over bulk units upwards
c
c Separation vector between layers no. 0 and (jbulk,k)
c Note, that the (artificial) right bulk position 
c vector corresponding to index (1,0) is equivalent
c to clay(nshiftr+1)!
c
          do i=1,3
            clayb(i)=clay(nshiftr+1,i)-cbulkr(i,1)+
     >               cbulkr(i,jbulk)+k*cbulkr0(i)
            cvec(i)=clay(nshiftl,i)-clayb(i)
          end do
c
c         ------------------------------
          call phi(cvec,sigma,fi,psi,xi)
c         ------------------------------
c         write(6,'(i3,3f10.5,2d15.7)') k,(cvec(i),i=1,3),fi,psi    
c
          sum1=sum1+fi   
          sum2=sum2+psi  
          if(dabs(fi).lt.tol) goto 30
c
        end do
        write(6,'(/'' WARNING!!! Madelungi: V(0)''/'' layer '',
     >    ''right sum did not converge up to bulk unit no.'',i4/)') k
   30   v0=v0+sum1*qbulkr(jbulk)+sum2*dbulkr(jbulk)
c
        sum1=0.d0
        sum2=0.d0
        do k=0,kmax
c Loop over bulk units upwards
c
c Separation vector between layers no. nintfc+1 and (jbulk,k)
c
          do i=1,3
            clayb(i)=clay(nshiftr+1,i)-cbulkr(i,1)+
     >               cbulkr(i,jbulk)+k*cbulkr0(i)
            cvec(i)=clay(nshiftr+1,i)-clayb(i)
          end do
c
c         ------------------------------
          call phi(cvec,sigma,fi,psi,xi)
c         ------------------------------
c         write(6,'(i3,3f10.5,2d15.7)') k,(cvec(i),i=1,3),fi,psi    
c
          sum1=sum1+fi
          sum2=sum2+psi
          if(dabs(fi).lt.tol) goto 40
c
        end do
        write(6,'(/'' WARNING!!! Madelungi: V(n+1)''/'' layer '',
     >    ''right sum did not converge up to bulk unit no.'',i4/)') k
   40   vnp1=vnp1+sum1*qbulkr(jbulk)+sum2*dbulkr(jbulk)
c
      end do
c
c Loop over interfacial layers
      do iq=1,nintfc
c
c Separation vector between layers no. 0 and iq
c
        cvec(1)=clay(nshiftl,1)-clay(nshiftl+iq,1)
        cvec(2)=clay(nshiftl,2)-clay(nshiftl+iq,2)
        cvec(3)=clay(nshiftl,3)-clay(nshiftl+iq,3)
c
c       ------------------------------
        call phi(cvec,sigma,fi,psi,xi)
c       ------------------------------
c
        v0=v0+qlay(iq)*fi+dlay(iq)*psi
c
c Separation vector between layers no. nintfc+1 and iq
c
        cvec(1)=clay(nshiftr+1,1)-clay(nshiftl+iq,1)
        cvec(2)=clay(nshiftr+1,2)-clay(nshiftl+iq,2)
        cvec(3)=clay(nshiftr+1,3)-clay(nshiftl+iq,3)
c
c       ------------------------------
        call phi(cvec,sigma,fi,psi,xi)
c       ------------------------------
c
        vnp1=vnp1+qlay(iq)*fi+dlay(iq)*psi
c
      end do
c
      bcoeff=vleft-v0
      do iq=1,nintfc
        bcoeff=bcoeff+
     >  fac*qlay(iq)*(clay(nshiftl+iq,3)-clay(nshiftl,3))
      end do
      do ibulk=1,nbulkl
        if(dabs(cbulkl(3,ibulk)-cbulkl(3,nbulkl)).lt.ctol)
     >  bcoeff=bcoeff+fac*dbulkl(ibulk)
      end do
c
      acoeff=vright-bcoeff-vnp1
      do iq=1,nintfc
         acoeff=acoeff+
     >   fac*qlay(iq)*(clay(nshiftr+1,3)-clay(nshiftl+iq,3))
     >   -2.0d0*fac*dlay(iq)
      end do
      do ibulk=1,nbulkr
        if(dabs(cbulkr(3,ibulk)-cbulkr(3,1)).lt.ctol)
     >  acoeff=acoeff-fac*dbulkr(ibulk)
      end do
      acoeff=acoeff/(clay(nshiftr+1,3)-clay(nshiftl,3))
c
c -------------------------------------------------------
c
c Now calculate Madelung potential for all layers
c
      do ip=1,nintfc
c       write(6,'(/'' Layer'',i3)') ip
c
        vmad(ip)=bcoeff+acoeff*(clay(nshiftl+ip,3)-clay(nshiftl,3))
c
c
c 1. G.ne.0 terms
c
        do jbulk=1,nbulkl
c         write(6,'(/'' jbulk='',i3)') jbulk
c Loop over inequivalent sublattices in the bulk to the left
c
          sum1=0.d0
          sum2=0.d0
          do k=0,-kmax,-1
c Loop over bulk units downwards
c
c Separation vector between layers no. ip and (jbulk,k)
c
            do i=1,3
              clayb(i)=clay(nshiftl,i)-cbulkl(i,nbulkl)+
     >                 cbulkl(i,jbulk)+k*cbulkl0(i)
              cvec(i)=clay(nshiftl+ip,i)-clayb(i)
            end do
c
c           ------------------------------
            call phi(cvec,sigma,fi,psi,xi)
c           ------------------------------
c           write(6,'(i3,3f10.5,2d15.7)') k,(cvec(i),i=1,3),fi,psi    
c
            sum1=sum1+fi
            sum2=sum2+psi
            if(dabs(fi).lt.tol) goto 50
c
          end do
          write(6,'(/'' WARNING!!! Madelungi: V('',i2,'')''/
     >    ''sum did not converge up to bulk unit no.'',i4/)') ip,k
   50     vmad(ip)=vmad(ip)+sum1*qbulkl(jbulk)+sum2*dbulkl(jbulk)
c
        end do
c
        do jbulk=1,nbulkr
c         write(6,'(/'' jbulk='',i3)') jbulk
c Loop over inequivalent sublattices in the bulk to the right
c
          sum1=0.d0
          sum2=0.d0
          do k=0,kmax
c Loop over bulk units upwards
c
c Separation vector between layers no. ip and (jbulk,k)
c
            do i=1,3
              clayb(i)=clay(nshiftr+1,i)-cbulkr(i,1)+
     >                 cbulkr(i,jbulk)+k*cbulkr0(i)
              cvec(i)=clay(nshiftl+ip,i)-clayb(i)
            end do
c
c           ------------------------------
            call phi(cvec,sigma,fi,psi,xi)
c           ------------------------------
c           write(6,'(i3,3f10.5,2d15.7)') k,(cvec(i),i=1,3),fi,psi    
c
            sum1=sum1+fi
            sum2=sum2+psi
            if(dabs(fi).lt.tol) goto 60
c
          end do
          write(6,'(/'' WARNING!!! Madelungi: V('',i2,'')''/
     >    ''sum did not converge up to bulk unit no.'',i4/)') ip,k
   60     vmad(ip)=vmad(ip)+sum1*qbulkr(jbulk)+sum2*dbulkr(jbulk)
c
        end do
c
        do iq=1,nintfc
c Loop over interfacial layers
c
c Separation vector between layers no. ip and iq
c
          cvec(1)=clay(nshiftl+ip,1)-clay(nshiftl+iq,1)
          cvec(2)=clay(nshiftl+ip,2)-clay(nshiftl+iq,2)
          cvec(3)=clay(nshiftl+ip,3)-clay(nshiftl+iq,3)
c
c         ------------------------------
          call phi(cvec,sigma,fi,psi,xi)
c         ------------------------------
c
          vmad(ip)=vmad(ip)+qlay(iq)*fi+dlay(iq)*psi
c
          emad=emad+qlay(ip)*fi*qlay(iq)+
     >              qlay(ip)*psi*dlay(iq)
     >             -dlay(ip)*psi*qlay(iq)
c
        end do
c
c 2. G=0 terms
c
        do iq=1,nintfc
c
          cperp=clay(nshiftl+ip,3)-clay(nshiftl+iq,3)
c
c  monopole + dipole terms
          vmad(ip)=vmad(ip)-fac*qlay(iq)*dabs(cperp)
          if(cperp.gt.ctol) vmad(ip)=vmad(ip)+2.0d0*fac*dlay(iq)
          if(dabs(cperp).lt.ctol) vmad(ip)=vmad(ip)+fac*dlay(iq)
c
          emad=emad-fac*qlay(ip)*dabs(cperp)*qlay(iq)
          if(cperp.gt.ctol) emad=emad+2.0d0*fac*qlay(ip)*dlay(iq)
     >                               -2.0d0*fac*dlay(ip)*qlay(iq)
          if(cperp.lt.ctol) emad=emad-2.0d0*fac*qlay(ip)*dlay(iq)
     >                               +2.0d0*fac*dlay(ip)*qlay(iq)
c
        end do
c
        if(itest.ge.2)
     >  write(6,'('' MAD total:'',t30,f20.10)') vmad(ip)
c
      end do
c -------------------------------------------------------
c        
      emad=emad/2.0d0
c
      return
      end
