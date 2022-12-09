c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine clebsch
c=======================
      implicit real*8 (a-h,o-z)
c
      dimension u1(72),ind1(72)
      dimension u2(72),ind2(72)
      common/cgc/u1,u2,ind1,ind2
      data tiny/1.0d-12/
c
c clebsch-gordan rectangular matrices to transform from (lm) to
c (kappa,my) basis
c
      do 1 i=1,72
      u1(i)=0.d0
      u2(i)=0.d0
      ind1(i)=0
      ind2(i)=0
    1 continue
c
      inr=0
      do 3 l=0,5
         twolp1=dfloat(2*l+1)
         do 3 m=-l,l
            inr=inr+1
c
c j=l-1/2
c
            kap=l
            if(kap.eq.0) goto 2
c
c ms=-1/2
c
            ir=2*kap*kap+kap+m
c           write(6,'(5i4)') l,m,inr,kap,ir
            u1(ir)=dsqrt((l+m)/twolp1)
            ind1(ir)=inr
c
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
    3 continue
c
c      write(6,*)
       do ir=1,72
         if(ind1(ir).eq.0) ind1(ir)=1
         if(ind2(ir).eq.0) ind2(ir)=1
c        write(6,'(i3,2x,i3,2x,i3,2x,2f20.10)') 
c    >   ir,ind1(ir),ind2(ir),u1(ir),u2(ir)
       end do
c
      return
      end
      subroutine ruthi(kap,my,ii,norder)
c=====================
c     find index for (kapa-my)-representation
c
      dimension kapdex(72),mydex(72)
      data kapdex/-1,-1,1,1,-2,-2,-2,-2,2,2,2,2,-3,-3,-3,-3,-3,-3
     * ,3,3,3,3,3,3,-4,-4,-4,-4,-4,-4,-4,-4,
     * 4,4,4,4,4,4,4,4,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,
     * 5,5,5,5,5,5,5,5,5,5,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6/
      data mydex/-1,1,-1,1,-3,-1,1,3,-3,-1,1,3,-5,-3,-1,1,3,5,
     * -5,-3,-1,1,3,5,-7,-5,-3,-1,1,3,5,7,
     * -7,-5,-3,-1,1,3,5,7,-9,-7,-5,-3,-1,1,3,5,7,9,
     * -9,-7,-5,-3,-1,1,3,5,7,9,-11,-9,-7,-5,-3,-1,1,3,5,7,9,11/
c
      do 100 i=1,norder
      if(kapdex(i).eq.kap) go to 150
  100 continue
      go to 500
  150 do 170 j=i,norder
      if(mydex(j).eq.my) go to 200
  170 continue
      go to 500
c
  200 ii=j
      return
  500 stop
      end
