c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine madelungv(qmoml,qmomlay,nintfc,lmax,sigma0,
     >                     vleft,vmad,vvac)
c=========================
c
c Calculate layer dependent Madelung potentials in monopole
c and dipole approximation for a free surface
c
c input:  
c         qmoml - moments of charge in the lhs bulk 
c         qmomlay - moments of charge densities in the interface
c         nintfc - number of layers in the interface
c         lmax - l-cutoff
c         sigma0 - Ewald parameter
c         vleft - Madelung potential for the first layer as
c                 calculated by 'madelungb' previously
c     in common/lay2d/
c         clay - layer generating vector
c         nbulk - number of layers in a bulk unit to the left
c     in common/extf/
c         qc - charge per 2D cell of external capacitor
c     in common/crystl2d/
c         vol2d - volume of the 2D unit cell 
c output: vmad - Madelung potentials for the interfacial layers
c         vvac - Madelung potential corresponding in the vacuum
c
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
      dimension cvec(3),clayb(3)
      dimension cbulk(3,mbulk),cbulk0(3)
      dimension qlay(mintfc),dlay(mintfc),vmad(mintfc)
      dimension qbulk(mbulk),dbulk(mbulk)
      complex*16 qmomlay(lmsup,mintfc),qmoml(lmsup,minprc)
c
      common/lay2d/clay(mtotal,3),nextra,nbulkl,nbulkr,
     &             nprc,ninprc(0:mprc+1)
      common/crystl2d/vol2d,volbz2d
      common/test/itest
      common/extf/qc
      common/madelung_energy/emad
c
      data tol/1.0d-8/,ctol/1.0d-6/,kmax/100/
c
      if(itest.ge.2) then
         write(6,'('' MADELUNGV> '')')
         call flush(6)
      end if
c
      emad=0.0d0
      pi=4.d0*datan(1.d0)
      fac=4.d0*pi/vol2d
      t1=dsqrt(3.d0)
      scale=dsqrt(vol2d)
      sigma=sigma0*scale
c Index of the uppermost left bulk layer (not part of the I-region)
      ninprcb=ninprc(0)
      nshift=(nextra+1)*ninprcb
      nbulk=nbulkl
      nbunit=ninprcb/nbulk
c
      do li=1,nintfc
        qlay(li)=qmomlay(1,li)
        dlay(li)=qmomlay(3,li)*t1
      end do
c
c Assign moments and non-primitive translation vectors for a bulk unit
c
      do ibulk=1,nbulk
        qbulk(ibulk)=0.d0
        dbulk(ibulk)=0.d0
        do iunit=1,nbunit
          ilay=(iunit-1)*nbulk+ibulk
          qbulk(ibulk)=qbulk(ibulk)+qmoml(1,ilay)
          dbulk(ibulk)=dbulk(ibulk)+qmoml(3,ilay)*t1
        end do
        qbulk(ibulk)=qbulk(ibulk)/nbunit    
        dbulk(ibulk)=dbulk(ibulk)/nbunit      
        do i=1,3
            cbulk(i,ibulk)=clay(nshift-nbulk+ibulk,i)-
     >                     clay(nshift-nbulk+1,i)
        end do
      end do
c
c Generating vector of the bulk units
c
      do i=1,3
        cbulk0(i)=clay(nshift,i)-clay(nshift-nbulk,i)
      end do
c
c -------------------------------------------------------
      if(itest.ge.2) then
       write(6,'(''   cbulk0'',3f10.7)') (cbulk0(j),j=1,3)
       do i=1,nbulk
         write(6,'(3f10.7)') (cbulk(j,i),j=1,3)
       enddo
       write(6,*) ' Charges and dipoles for L-region'
       do i=1,nbulk  
         write(6,'(2f15.10)') qbulk(i),dbulk(i)
       enddo
       write(6,*) ' Charges and dipoles for I-region'
       do i=1,nintfc
         write(6,'(2f15.10)') qlay(i),dlay(i)
       enddo
      endif
c -------------------------------------------------------
c
c Calculate normalization term for the Madelung constants
c
      if(itest.ge.2) 
     >  write(6,'(''Normalization term for Madelung potentials'')')
c
      vnorm=vleft
      do ibulk=1,nbulk
        if(dabs(cbulk(3,ibulk)-cbulk(3,nbulk)).lt.ctol) 
     >  vnorm=vnorm+fac*dbulk(ibulk)
      end do
c
      do jbulk=1,nbulk
c Loop over inequivalent sublattices in bulk
c
        sum1=0.d0
        sum2=0.d0
        do k=0,-kmax,-1
c Loop over bulk units downwards
c
c Separation vector between layers no. 0 and (jbulk,k)
c Note, that the (artificial) left bulk position 
c vector corresponding to index (nbulk,0) is equivalent
c to clay(nshift,*)!
c
          do i=1,3
            clayb(i)=clay(nshift,i)-cbulk(i,nbulk)
     >                             +cbulk(i,jbulk)+k*cbulk0(i)
            cvec(i)=clay(nshift,i)-clayb(i)
          end do
c
c         ------------------------------
          call phi(cvec,sigma,fi,psi,xi)
c         ------------------------------
c
          sum1=sum1+fi
          sum2=sum2+psi
          if(dabs(fi).lt.tol) goto 10
c
        end do   
        write(6,'(/'' WARNING!!! Madelungv: ''/'' layer '',
     >    ''sum did not converge up to bulk unit no.'',i4/)') k
   10   vnorm=vnorm-sum1*qbulk(jbulk)-sum2*dbulk(jbulk)
      if(itest.ge.2) write(6,'(''   left-bulk converges at layer'',
     &   i5,f20.10)') k,vnorm
c
      end do
c
      do iq=1,nintfc
c Loop over interfacial layers
c
c Separation vector between layers no. 1 and iq
c
        cvec(1)=clay(nshift,1)-clay(nshift+iq,1)
        cvec(2)=clay(nshift,2)-clay(nshift+iq,2)
        cvec(3)=clay(nshift,3)-clay(nshift+iq,3)
c
c       ------------------------------
        call phi(cvec,sigma,fi,psi,xi)
c       ------------------------------
c
        vnorm=vnorm-qlay(iq)*fi-dlay(iq)*psi
c
      end do
      if(itest.ge.2) write(6,'(''   interface contribution '',
     &   f20.10)') vnorm
c
      if(itest.ge.2)  
     >  write(6,'(''  In total:'',f20.10)') vnorm
c
c -------------------------------------------------------
c
c Now calculate Madelung potential for all layers
c
      do ip=1,nintfc
c
        vmad(ip)=vnorm 
        if(itest.ge.2) 
     &   write(6,'(/,'' Madelung constant for layer'',i3)') ip
c
c 1. G.ne.0 terms
c
        do jbulk=1,nbulk
c Loop over inequivalent sublattices in bulk
c
          sum1=0.d0
          sum2=0.d0
          do k=0,-kmax,-1
c Loop over bulk units downwards
c
c Separation vector between layers no. ip and (jbulk,k)
c
            do i=1,3
              clayb(i)=clay(nshift,i)-cbulk(i,nbulk)
     >                               +cbulk(i,jbulk)+k*cbulk0(i)
              cvec(i)=clay(nshift+ip,i)-clayb(i)
            end do
c
c           ------------------------------
            call phi(cvec,sigma,fi,psi,xi)
c           ------------------------------
c
            sum1=sum1+fi
            sum2=sum2+psi
            if(dabs(fi).lt.tol) goto 20
c
          end do
          write(6,'(/'' WARNING!!! Madelungv: ''/'' layer '',
     >    ''sum did not converge up to bulk unit no.'',i4/)') k
   20     vmad(ip)=vmad(ip)+sum1*qbulk(jbulk)+sum2*dbulk(jbulk)
      if(itest.ge.2) write(6,'(''   left-bulk converges at layer'',
     &   i5,2f20.10)') k,sum1,sum2
c
        end do
c
        ssum=0.d0
        do iq=1,nintfc
c Loop over interfacial layers
c
c Separation vector between layers no. ip and iq
c
          cvec(1)=clay(nshift+ip,1)-clay(nshift+iq,1)
          cvec(2)=clay(nshift+ip,2)-clay(nshift+iq,2)
          cvec(3)=clay(nshift+ip,3)-clay(nshift+iq,3)
c
c         ------------------------------
          call phi(cvec,sigma,fi,psi,xi)
c         ------------------------------
c
          vmad(ip)=vmad(ip)+qlay(iq)*fi+dlay(iq)*psi
          ssum = ssum +qlay(iq)*fi+dlay(iq)*psi 
c
          emad=emad+qlay(ip)*fi*qlay(iq)+
     >              qlay(ip)*psi*dlay(iq)
     >             -dlay(ip)*psi*qlay(iq)
c
        end do
        if(itest.ge.2) write(6,'(''   contribution from interface '',
     &   f20.10)') ssum
c
c 2. G=0 terms
c
        ssum = 0.d0
        do iq=1,nintfc 
c
c  monopole term
          vmad(ip)=vmad(ip)+2.0d0*fac*qlay(iq)*
     >    (dmin1(clay(nshift+ip,3),clay(nshift+iq,3))-clay(nshift,3))
          emad=emad+2.d0*fac*qlay(ip)*qlay(iq)*
     >    dmin1(clay(nshift+ip,3),clay(nshift+iq,3))
c
c  dipole term
          diffc=clay(nshift+ip,3)-clay(nshift+iq,3)
          if(diffc.gt.ctol) then
            vmad(ip)=vmad(ip)+2.0d0*fac*dlay(iq)
            ssum = ssum + 2.0d0*fac*dlay(iq)
            emad=emad+2.d0*fac*qlay(ip)*dlay(iq)
     >               -2.d0*fac*dlay(ip)*qlay(iq)
          else if(dabs(diffc).lt.ctol) then
            vmad(ip)=vmad(ip)+fac*dlay(iq)
            ssum = ssum + fac*dlay(iq)
          end if
c
        end do
        if(itest.ge.2) write(6,'(''   G=0 terms                   '',
     &   f20.10)') ssum
c
c  add potential from external capacitor
c
        cperpc=clay(nshift+nintfc,3) ! external capacitor is in the last vac layer (Szunyogh 18.05.2018)
        vmad(ip)=vmad(ip)+2.0d0*fac*qc*
     >  (dmin1(clay(nshift+ip,3),cperpc)-clay(nshift,3))
c
        if(itest.ge.2)
     >  write(6,'('' MAD total:'',t30,f20.10)') vmad(ip)
c
      end do
c
c  vacuum potential level
c
      vvac=vmad(nintfc)
      do iq=1,nintfc
        if(dabs(clay(nshift+iq,3)-clay(nshift+nintfc,3)).lt.ctol)
     >  vvac=vvac+fac*dlay(iq)
      end do
      if(itest.ge.2)
     >write(6,'('' Vvac      '',t30,f20.10)') vvac
c
c -------------------------------------------------------
c        
      emad=emad/2.d0
c
      return
      end
