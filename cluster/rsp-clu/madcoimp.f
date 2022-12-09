c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine madcoimp(madmax,iimp,nimp,nposimp,mlmij)
c
c **********************************************************
c * input:                                                 *
c *        madmax - max momentum of the Madelung expansion *  
c *        nimp - no. of the impurities                    *
c *        nposimp - positions of impurities               *
c * output:                                                *
c *        mlmij - the Madelung const. corresponding to    * 
c *                the qn. lm between the ith anh jth imp. *
c *--------------------------------------------------------*
c *                                                        *
c * calls: rotvec - rotate the lattice vectors, due to the *
c *                 rotated basis                          *
c *        spherh - calculate spherical harmonics          *
c **********************************************************  
c ===================================================
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
c -----------------------------
c for impurity calculation --> param.h
      parameter(maxmadi=lmsup)
c for off-diagonal calculation
c     parameter(maxmadi=4)
c -----------------------------
c
      real*8 rdiff(3)
      real*8 srdiff(3)
      real*8 rdiff2
c
      integer nimp
      integer nposimp(3,mimp)
      integer iimp
      integer jimp
      integer kimp
      integer i
c
      complex*16 ylm(-lshape:lshape,0:lshape)
c
      complex*16 mlmij(maxmadi,mimp)
c
c - common blocks -
      real*8 cvec
      integer nextra
      integer nbulkl
      integer nbulkr
      integer nprc
      integer ninprc
      common/lay2d/cvec(mtotal,3),nextra,nbulkl,nbulkr,
     &             nprc,ninprc(0:mprc+1)
c
      real*8 a1
      real*8 a2
      common/brav2d/a1(2),a2(2)
c
      common/test/itest
c ===================================================
c
      if(itest.ge.2) write(6,*) '<madcoimp>: Beginn'
      pi=4.d0*datan(1.d0)
      sqrpi=dsqrt(pi)
      fac=4.d0*sqrpi
c
      n0 = ninprc(0)*(nextra+1)
c
      jimp=iimp
      do kimp=1,nimp
        lj=nposimp(3,jimp)
        lk=nposimp(3,kimp)
c
        if(kimp.ne.jimp) then
            rdiff(1)=(nposimp(1,kimp)-nposimp(1,jimp))*a1(1)+
     >               (nposimp(2,kimp)-nposimp(2,jimp))*a2(1)+
     >               (cvec(lk+n0,1)-cvec(lj+n0,1))
c
            rdiff(2)=(nposimp(1,kimp)-nposimp(1,jimp))*a1(2)+
     >               (nposimp(2,kimp)-nposimp(2,jimp))*a2(2)+
     >               (cvec(lk+n0,2)-cvec(lj+n0,2))
c
            rdiff(3)=(cvec(lk+n0,3)-cvec(lj+n0,3))
c
          rdiff2=dsqrt(rdiff(1)*rdiff(1)+rdiff(2)*rdiff(2)+
     >          rdiff(3)*rdiff(3))
c
          call spherh(madmax,rdiff(1),rdiff(2),rdiff(3),ylm) 
c
          i=0
          do l=0,madmax
             do m=-l,l
               i=i+1 
               mlmij(i,kimp)=fac*dconjg(ylm(m,l))/(rdiff2**(l+1))
             end do
          end do
c
        end if
      end do
c
      return
      end
