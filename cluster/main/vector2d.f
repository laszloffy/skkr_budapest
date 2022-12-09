c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine vector2d(a1,a2,b1,b2)
c======================
c
c determine basis vectors of the reciprocal lattice
c
c input : a1,a2 - basis translations in direct space
c         b1,b2 - basis translations in reciprocal space
c output:     rvec,distr,maxr in common/dirvec/ -
c             rvec : translation vectors in direct space
c             distr: diameter of sphere where rvec's take place
c             maxr : number of rvecs's
c         veck,distk,maxk in common/recvec/ - same as in dirvec for rec. space 
c         vol,volbz in common/crystl2d/ - volumen of unit cell in direct
c                                         and reciprocal space
c 
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      parameter(pi=3.1415926535897932d0,twopi=pi+pi)
c
      dimension a1(2),a2(2),b1(2),b2(2)
      logical overflow
c
      common/dirvec2d/rvec(2,mdir),distr,maxr
      common/recvec2d/veck(2,mrec),distk,maxk
      common/crystl2d/vol,volbz
      common/test/itest
c
      vol  =dabs(a1(1)*a2(2)-a1(2)*a2(1))
      volbz=dabs(b1(1)*b2(2)-b1(2)*b2(1))
c
c
c determine the longest diagonal of the real and reciprocal
c unit cells. add the diagonal to the maximum length to cover
c all possible nonprimitive vectors.
c
      dr1=0.d0
      dr2=0.d0
c
      dk1=0.d0
      dk2=0.d0
c
      do 2 i=1,2
      dr1=dr1+(a1(i)+a2(i))**2
      dr2=dr2+(a1(i)-a2(i))**2
c
      dk1=dk1+(b1(i)+b2(i))**2
      dk2=dk2+(b1(i)-b2(i))**2
    2 continue
c
      diagr=dsqrt(dmax1(dr1,dr2))
      diagk=dsqrt(dmax1(dk1,dk2))
      diagrm=dsqrt(dmin1(dr1,dr2))
      diagkm=dsqrt(dmin1(dk1,dk2))
c
      distr1=twopi/dsqrt(b1(1)**2+b1(2)**2)
      distr2=twopi/dsqrt(b2(1)**2+b2(2)**2)
c
      distk1=twopi/dsqrt(a1(1)**2+a1(2)**2)
      distk2=twopi/dsqrt(a2(1)**2+a2(2)**2)
c
      mdr1=1
      mdk1=1
c
      do md=mdr1,mdir
        distr=md*diagrm
        totr =distr+diagr
c determine number of cells required to cover sphere with radius totr
        nr1=idint(totr/distr1)+1
        nr2=idint(totr/distr2)+1
c generate real space lattice vectors
        call latgen2d(nr1,nr2,mdir,rvec,maxr,totr,a1,a2,overflow)
        if(overflow) goto 100
      end do
  100 continue
      distr=(md-1)*diagrm
      totr =distr+diagr
c determine number of cells required to cover sphere with radius totr
      nr1=idint(totr/distr1)+1
      nr2=idint(totr/distr2)+1
c generate real space lattice vectors
      call latgen2d(nr1,nr2,mdir,rvec,maxr,totr,a1,a2,overflow)
      if (overflow)
     &stop 'ERROR in <vector2D>: unexpected overflow in direct space !'
c
      do md=mdk1,mrec
        distk=md*diagkm
        totk =distk+diagk
c determine number of cells required to cover sphere with radius totk
        nk1=idint(totk/distk1)+1
        nk2=idint(totk/distk2)+1
c generate reciprocal space lattice vectors
        call latgen2d(nk1,nk2,mrec,veck,maxk,totk,b1,b2,overflow)
        if(overflow) goto 200
      end do
  200 continue
      distk=(md-1)*diagkm
      totk =distk+diagk
c determine number of cells required to cover sphere with radius totk
      nk1=idint(totk/distk1)+1
      nk2=idint(totk/distk2)+1
c generate reciprocal space lattice vectors
      call latgen2d(nk1,nk2,mrec,veck,maxk,totk,b1,b2,overflow)
      if (overflow)
     &stop 
     &'ERROR in <vector2D>: unexpected overflow in reciprocal space !'
c
      if(itest.ge.1) then
      write(6,'(/2x,''routine VECTOR2D>''/)')
      write(6,'(a,1x,2f10.5)')
     >' diagonals of 2D real and rec. unit cells:',diagr,diagk
      write(6,'(a,1x,2f10.5/)')
     >' max. length of 2D real and rec. vectors :',distr,distk
      write(6,'(a,2i5)') 'nr in 2d: ',nr1,nr2
      write(6,'(a,2i5)') 'nk in 2d: ',nk1,nk2
      write(6,902) vol,a1,a2
      write(6,901) b1,b2
      write(6,*) ' max-r in 2d: ',maxr
      write(6,*) ' max-k in 2d: ',maxk
      end if
c
      return
  902 format(/' unit cell area: ',f15.8//
     * ' real lattice unit cell in 2D: '/2(/2(1x,f15.10)))
  901 format(/' reciprocal lattice unit cell in 2D: '/2(/2(1x,
     *      f15.10))/)
      end
