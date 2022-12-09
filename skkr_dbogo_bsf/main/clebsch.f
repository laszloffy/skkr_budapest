c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine clebsch
c     ==================
c
      implicit real*8 (a-h,o-z)
c
      dimension u1(72),ind1(72)
      dimension u2(72),ind2(72)
      common/cgc/u1,u2,ind1,ind2
c
c clebsh-gordan coefficients between (lms) and (kappa,my) bases
c
      u1=0.d0
      u2=0.d0
      ind1=0
      ind2=0
c
      inr=0
      do l=0,5
         twolp1=dfloat(2*l+1)
         do m=-l,l
            inr=inr+1
c
c j=l-1/2
c
            kap=l
            if(kap.eq.0) goto 2
c
c ms=-1/2
            if(m+l.eq.0) goto 1
c
            ir=2*kap*kap+kap+m
c           write(6,'(5i4)') l,m,inr,kap,ir
            u1(ir)=dsqrt((l+m)/twolp1)
            ind1(ir)=inr
c
    1       continue
c ms=+1/2
c
            ir=2*kap*kap+kap+m+1
c           write(6,'(5i4)') l,m,inr,kap,ir
            u2(ir)=-dsqrt((l-m)/twolp1)
            ind2(ir)=inr
    2       continue
c
c j=l+1/2
c
            kap=-l-1
c
c ms=-1/2
c
            ir=2*kap*kap+kap+m
c           write(6,'(5i4)') l,m,inr,kap,ir
            u1(ir)=dsqrt((l-m+1)/twolp1)
            ind1(ir)=inr
c
c ms=+1/2
c
            ir=2*kap*kap+kap+m+1
c           write(6,'(5i4)') l,m,inr,kap,ir
            u2(ir)=dsqrt((l+m+1)/twolp1)
            ind2(ir)=inr
c
         end do         
      end do         
c
c      write(6,*)
       do ir=1,72
         if(ind1(ir).eq.0) ind1(ir)=1
         if(ind2(ir).eq.0) ind2(ir)=1
c        write(6,'(i3,6x,i3,f18.10,6x,i3,f18.10)') 
c    >   ir,ind1(ir),u1(ir),ind2(ir),u2(ir)
       end do
c
      return
      end
