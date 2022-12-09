c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
c***vector.spg  processed by SPAG 5.11R  at 12:38 on 17 Feb 2000
      subroutine vector3d(a1,a2,a3,nbulk,r1,r2,r3,
     >                    rvec,vol,distr,maxr,
     >                    veck,volbz,distk,maxk)
c --------------------------------------------------------------------**
c determine basis vectors of the reciprocal lattice
c
c input : a1,a2,a3 - basis translations in direct space
c         nbulk - number of inequivalent positions in bulk unit cell
c output: r1,r2,r3 - basis translations in rec. space
c         rvec - translation vectors in direct space
c         vol  - volumen of unit cell in direct space
c         distr - diameter of sphere where rvec's take place
c         maxr - number of rvecs's
c         veck,volbz,distk,maxk - same as above for rec. space
c                                                                     **
c ------------------------- common variables -------------------------**
c                                                                     **
c  modifies    distk   distr   sws     vol                            **
c  passes arg  maxk    maxr    rvec    veck                           **
c  uses value  distk   distr   itest   maxk    maxr    rvec    sws    **
c              veck    vol                                            **
c                                                                     **
c ----------------------- external subprograms -----------------------**
c                                                                     **
c  calls       latgen3d                                               **
c  called by   genlatt3d                                              **
c                                                                     **
c --------------------------------------------------------------------**
c
      implicit real*8 (a-h,o-z)
c     implicit none
C*** Start of declarations inserted by SPAG
c     real*8 a1,a2,a3,det,diagk,diagkm,diagr,diagrm,distk1,distk2,
c    &       distk3,distk,distr1,distr2,distr3,distr,dk1,dk2,dk3,dk4,
c    &       dr1,dr2,dr3,dr4,pi
c     real*8 r1,r2,r3,rvec,sws,totk,totr,twopi,veck,vol,volbz
c     integer i,itest,maxk,maxr,md,mdr1,mdk1,nk1,nk2,nk3,nr1,nr2,nr3
      logical overflow
C*** End of declarations inserted by SPAG
c
      include '../param.h'
c
      parameter (pi=3.1415926535897932D0,twopi=pi+pi)
c
      dimension a1(3),a2(3),a3(3),r1(3),r2(3),r3(3)
      dimension rvec(3,mdir)
      dimension veck(3,mrec)
c
      common /test  / itest
      common /swscheck/ sws
c
      r1(1)=a2(2)*a3(3)-a2(3)*a3(2)
      r2(1)=a3(2)*a1(3)-a3(3)*a1(2)
      r3(1)=a1(2)*a2(3)-a1(3)*a2(2)
      r1(2)=a2(3)*a3(1)-a2(1)*a3(3)
      r2(2)=a3(3)*a1(1)-a3(1)*a1(3)
      r3(2)=a1(3)*a2(1)-a1(1)*a2(3)
      r1(3)=a2(1)*a3(2)-a2(2)*a3(1)
      r2(3)=a3(1)*a1(2)-a3(2)*a1(1)
      r3(3)=a1(1)*a2(2)-a1(2)*a2(1)
c
      det=a1(1)*r1(1)+a1(2)*r1(2)+a1(3)*r1(3)
      vol=dabs(det)
      sws=(3.D0/(4.D0*pi)*vol/nbulk)**(1.D0/3.D0)
      volbz=twopi**3/vol
c
      do i=1,3
        r1(i)=r1(i)*twopi/det
        r2(i)=r2(i)*twopi/det
        r3(i)=r3(i)*twopi/det
      enddo
c
c determine the longest diagonal of the real and reciprocal
c unit cells. add the diagonal to the maximum length to cover
c all possible nonprimitive vectors.
c
      dr1=0.D0
      dr2=0.D0
      dr3=0.D0
      dr4=0.D0
c
      dk1=0.D0
      dk2=0.D0
      dk3=0.D0
      dk4=0.D0
c
      do i=1,3
        dr1=dr1+(a1(i)+a2(i)+a3(i))**2
        dr2=dr2+(a1(i)-a2(i)+a3(i))**2
        dr3=dr3+(a1(i)+a2(i)-a3(i))**2
        dr4=dr4+(a1(i)-a2(i)-a3(i))**2
c
        dk1=dk1+(r1(i)+r2(i)+r3(i))**2
        dk2=dk2+(r1(i)-r2(i)+r3(i))**2
        dk3=dk3+(r1(i)+r2(i)-r3(i))**2
        dk4=dk4+(r1(i)-r2(i)-r3(i))**2
      enddo
c
      diagr=dsqrt(dmax1(dr1,dr2,dr3,dr4))
      diagk=dsqrt(dmax1(dk1,dk2,dk3,dk4))
      diagrm=dsqrt(dmin1(dr1,dr2,dr3,dr4))
      diagkm=dsqrt(dmin1(dk1,dk2,dk3,dk4))
c
      distr1=twopi/dsqrt(r1(1)**2+r1(2)**2+r1(3)**2)
      distr2=twopi/dsqrt(r2(1)**2+r2(2)**2+r2(3)**2)
      distr3=twopi/dsqrt(r3(1)**2+r3(2)**2+r3(3)**2)
c
      distk1=twopi/dsqrt(a1(1)**2+a1(2)**2+a1(3)**2)
      distk2=twopi/dsqrt(a2(1)**2+a2(2)**2+a2(3)**2)
      distk3=twopi/dsqrt(a3(1)**2+a3(2)**2+a3(3)**2)
c
c estimate the minimum of shells in both spaces   
c     mdr1=int(0.5D0*(mdir*vol)**(1.D0/3.D0)/diagr)
      mdr1=1
c     mdk1=int(0.5D0*(mrec*volbz)**(1.D0/3.D0)/diagk)
      mdk1=1
c
      do md=mdr1,mdir
        distr=md*diagrm
        totr =distr+diagr
c determine number of cells required to cover sphere with radius totr
        nr1=idint(totr/distr1)+1
        nr2=idint(totr/distr2)+1
        nr3=idint(totr/distr3)+1
c generate real space lattice vectors
        call latgen3d(nr1,nr2,nr3,mdir,rvec,maxr,totr,a1,a2,a3,overflow)
        if (overflow) goto 100
      enddo
 100  continue
      distr=(md-1)*diagrm
      totr =distr+diagr
c determine number of cells required to cover sphere with radius totr
      nr1=idint(totr/distr1)+1
      nr2=idint(totr/distr2)+1
      nr3=idint(totr/distr3)+1
c generate real space lattice vectors
      call latgen3d(nr1,nr2,nr3,mdir,rvec,maxr,totr,a1,a2,a3,overflow)
      if (overflow) 
     &stop 'ERROR in <vector3D>: unexpected overflow in direct space !'

      do md=mdk1,mrec
        distk=md*diagkm
        totk=distk+diagk
c determine number of cells required to cover sphere with radius totk
        nk1=idint(totk/distk1)+1
        nk2=idint(totk/distk2)+1
        nk3=idint(totk/distk3)+1
c generate reciprocal space lattice vectors
        call latgen3d(nk1,nk2,nk3,mrec,veck,maxk,totk,r1,r2,r3,overflow)
        if (overflow) goto 200
      enddo
 200  continue
      distk=(md-1)*diagkm
      totk=distk+diagk
c determine number of cells required to cover sphere with radius totk
      nk1=idint(totk/distk1)+1
      nk2=idint(totk/distk2)+1
      nk3=idint(totk/distk3)+1
c generate reciprocal space lattice vectors
      call latgen3d(nk1,nk2,nk3,mrec,veck,maxk,totk,r1,r2,r3,overflow)
      if (overflow) 
     &stop 
     &'ERROR in <vector3D>: unexpected overflow in reciprocal space !'

c
      if (itest.ge.1) then
c       write (6,'(/2x,''routine VECTOR3D>''/)')
        write (6,'(/a,1x,2f10.5)') 
     &                         ' diagonals of real and rec. unit cells:'
     &                         ,diagr,diagk
        write (6,'(a,1x,2f10.5/)') 
     &                         ' max. length of real and rec. vectors :'
     &                         ,distr,distk
        write (6,'(a,3i5)') 'nr: ',nr1,nr2,nr3
        write (6,'(a,3i5)') 'nk: ',nk1,nk2,nk3
        write (6,910) vol,sws,a1,a2,a3
        write (6,920) r1,r2,r3
        write (6,*) ' max-r: ',maxr
        write (6,*) ' max-k: ',maxk
      endif
c
      return
  910 format (/' unit cell volume : ',f20.12,/,' average ws-radius: ',
     &        f20.12,/,' real lattice unit cell: '/3(/3(1x,f15.10)))
  920 format (/' reciprocal lattice unit cell: '/3(/3(1x,f15.10))/)
      end
