c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine inouh(dp,dq,dr,dq1,dfl,dv,z,test,nuc,jc,
     >                 dep,deq,db,dvc,dsal,dk,dm,dpno,dqno)
c
c     ****************************************
c     * start values for outward integration *
c     ****************************************
c
      implicit real*8(a-h,o-z)
c
      include '../param.h'
c
      dimension dp(nrad),dq(nrad),dr(nrad)
      dimension dep(10),deq(10)
      dimension dpno(4,iorb),dqno(4,iorb)
c
      do i=1,10
        dp(i)=0.0d0
        dq(i)=0.0d0
      end do
c
      if(nuc.le.0) then
c
        dval=z/dvc
        deva1=-dval
        deva2=dv/dvc+dval/dr(1)-db
        deva3=0.0d0
        if(dk.le.0.) then
          dbe=(dk-dfl)/dval
        else
          dbe=dval/(dk+dfl)
        end if
        dq(10)=dq1
        dp(10)=dbe*dq1
c
      else
c
        dval=dv+z*(3.0d0-dr(1)*dr(1)/(dr(nuc)*dr(nuc)))
     *             /(dr(nuc)+dr(nuc))
        deva1=0.0d0
        deva2=(dval-3.0d0*z/(dr(nuc)+dr(nuc)))/dvc-db
        deva3=z/(dr(nuc)*dr(nuc)*dr(nuc)*dsal)
        if(dk.le.0.) then
          dp(10)=dq1
        else
          dq(10)=dq1
        end if
c        
      end if
c
      do i=1,5
        dp(i)=dp(10)
        dq(i)=dq(10)
        dep(i)=dp(i)*dfl
        deq(i)=dq(i)*dfl
      end do
c
      m=1
 10   dm=m+dfl
      dsum=dm*dm-dk*dk+deva1*deva1
      dqr=(dsal-deva2)*dq(m+9)-deva3*dq(m+7)
      dpr=deva2*dp(m+9)+deva3*dp(m+7)
      dval=((dm-dk)*dqr-deva1*dpr)/dsum
      dsum=((dm+dk)*dpr+deva1*dqr)/dsum
c
      j=-1
      do i=1,5
        dpr=dr(i)**m
        dqr=dsum*dpr
        dpr=dval*dpr
        if(m.gt.1) then
          if(dabs(dpr/dp(i)).le.test.and.abs(dqr/dq(i)).le.test) j=1
        end if
        dp(i)=dp(i)+dpr
        dq(i)=dq(i)+dqr
        dep(i)=dep(i)+dpr*dm
        deq(i)=deq(i)+dqr*dm
      end do
      if (j.eq.1) go to 99
      dp(m+10)=dval
      dq(m+10)=dsum
      m=m+1
c
      if (m.le.20) go to 10
c
      write(6,'('' core level='',i3)') jc
      write(6,'('' expansion at the origin does not converge'')')
      stop ' in inouh'
c
 99   do i=1,4
        dpno(i,jc)=dp(i+9)
        dqno(i,jc)=dq(i+9)
      end do
c
      return
      end
