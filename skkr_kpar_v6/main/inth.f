c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine inth(dp,dq,dv,dr,dep,deq,db,dvc,dsal,dk,dm)
c
c     *****************************************************
c     * integration of the dirac equation using a 5 point *
c     * adams method                                      *
c     *****************************************************
c
      implicit real*8(a-h,o-z)
c
      dimension dep(10),deq(10)
c
c     dpas    -    delx, logarithmic increment for radius
c     dm      -    dpas/720
c     dkoef1  -    475/502
c     dkoef2  -    27/502
c
      dkoef1 = 475.0d0/502.0d0
      dkoef2 = 27.0d0/502.0d0
c
      dpr=dp+dm*((251.0d0*dep(1)+2616.0d0*dep(3)+1901.0d0*dep(5))
     *    -(1274.0d0*dep(2)+2774.0d0*dep(4)))
      dqr=dq+dm*((251.0d0*deq(1)+2616.0d0*deq(3)+1901.0d0*deq(5))
     *    -(1274.0d0*deq(2)+2774.0d0*deq(4)))
c
      do 13 i=2,5
      dep(i-1)=dep(i)
   13 deq(i-1)=deq(i)
c
      dsum=(db-dv/dvc)*dr
      dep(5)=-dk*dpr+(dsal*dr+dsum)*dqr
      deq(5)=dk*dqr-dsum*dpr
c
      dp=dp+dm*((106.0d0*dep(2)+646.0d0*dep(4)+251.0d0*dep(5))
     *   -(19.0d0*dep(1)+264.0d0*dep(3)))
      dq=dq+dm*((106.0d0*deq(2)+646.0d0*deq(4)+251.0d0*deq(5))
     *   -(19.0d0*deq(1)+264.0d0*deq(3)))
c
      dp=dkoef1*dp+dkoef2*dpr
      dq=dkoef1*dq+dkoef2*dqr
      dep(5)=-dk*dp+(dsal*dr+dsum)*dq
      deq(5)=dk*dq-dsum*dp
c
      return
      end
