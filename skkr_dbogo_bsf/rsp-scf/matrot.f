c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine matrot(r0,r1,lmax,dmat,dmatp,rmat,rmatp,tvec,phi)
c =====================
c
c if r1 = R * r0, dmat = D(R) and dmatp = D(R)+ 
c
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
      dimension r0(3),r1(3),rr(3),drot(3,3),tvec(3),tlm(4)
      complex*16 dmat(kmymaxp,kmymaxp),dmatp(kmymaxp,kmymaxp)
      complex*16 cmat(kmymaxp,kmymaxp)
      complex*16 rmat(lmsup,lmsup),rmatp(lmsup,lmsup)
c
      common/test/itest
c
      pi=4.d0*datan(1.d0)
      sth=dsqrt(3.d0)
c
      r0mod=0.d0
      r1mod=0.d0
      do i=1,3
        r0mod=r0mod+r0(i)*r0(i)
        r1mod=r1mod+r1(i)*r1(i)
      end do
      r0mod=dsqrt(r0mod)
      r1mod=dsqrt(r1mod)
      do i=1,3
        r0(i)=r0(i)/r0mod
        r1(i)=r1(i)/r1mod
      end do
c
c Find normal vector and angle of rotation
c 
      cosphi=r0(1)*r1(1)+r0(2)*r1(2)+r0(3)*r1(3)
      if(dabs(cosphi-1.d0).lt.1.d-10) then
        phi=0.d0
        tvec(1)=r0(1) 
        tvec(2)=r0(2) 
        tvec(3)=r0(3) 
        goto 10
      else if(dabs(cosphi+1.d0).lt.1.d-10) then
        phi=pi
        r0xy=dsqrt(r0(1)**2+r0(2)**2)
        if(r0xy.gt.1.d-10) then
          tvec(1)=-r0(2)/r0xy
          tvec(2)=r0(1)/r0xy
          tvec(3)=0.d0
        else
          tvec(1)=1.d0
          tvec(2)=0.d0
          tvec(3)=0.d0
        end if
        goto 10
      end if
      phi=dacos(cosphi)
c
c T=(R0xR1)/|R0xR1|
c
      tvec(1)=r0(2)*r1(3)-r0(3)*r1(2)
      tvec(2)=r0(3)*r1(1)-r0(1)*r1(3)
      tvec(3)=r0(1)*r1(2)-r0(2)*r1(1)
      tmod=dsqrt(tvec(1)**2+tvec(2)**2+tvec(3)**2)
      tvec(1)=tvec(1)/tmod
      tvec(2)=tvec(2)/tmod
      tvec(3)=tvec(3)/tmod
c
   10 continue
c     
      if(itest.gt.3) then
      write(6,'(''   r0= '',3f10.4)') r0
      write(6,'(''   r1= '',3f10.4)') r1
      write(6,'('' Normal vector: '',3f10.4)') tvec
      write(6,'('' Angle in degree:'',f10.4)') phi*180.d0/pi
      end if
c
      cosp2=dcos(phi/2.d0)
      sinp2=dsin(phi/2.d0)
      tlm(1)=cosp2
      tlm(2)=sinp2*tvec(1)
      tlm(3)=sinp2*tvec(2)
      tlm(4)=sinp2*tvec(3)
c
c Set up matrix of rotation for spinor spherical harmonics 
      call matr(tlm,lmax,dmat,dmatp)
c Set up matrix of rotation for spherical harmonics 
      call matnr(tlm,lmax,rmat,rmatp)
c
      return
      end
