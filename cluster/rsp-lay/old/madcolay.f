c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine madcolay(madmax,nimp,nposimp,mlmij)
c
c **********************************************************
c * input:                                                 *
c *        madmax - max momentum of the Madelung expansion *  
c *        maxmad - momentum dimension of mlmij            *
c *        nimp - no. of the impurities                    *
c *        nposimp - positions of impurities               *
c *        vec - the normal vector in the rotation z --> B *
c *        phi - the angle of rotation                     *
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
c for impurity calculation
c     parameter(maxmad=9)
c for off-diagonal calculation
      parameter(maxmadl=4)
c -----------------------------
c
c     real*8 vecn(3,mintfc)
c     real*8 phi(mintfc)
c
      real*8 rdiff(3)
c     real*8 srdiff(3)
c     real*8 vecna(3)
      real*8 rdiff2
      real*8 sqrpi
      real*8 fac
      real*8 pi
c
      integer nposimp(3,mimp)
      integer nimp
      integer n0
c
      integer iimp
      integer jimp
      integer lj
      integer li
      integer ii
      integer ll
      integer mm
c
      complex*16 ylm(-lsup:lsup,0:lsup)
c
      complex*16 mlmij(maxmadl,mimp,mimp)
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
c ===================================================
c TEST
c     write(6,*) 'maxmad=',maxmad
c     call flush(6)
c TEST
c
c
      pi=4.d0*datan(1.d0)
      sqrpi=dsqrt(pi)
      fac=4.d0*sqrpi
c
      n0 = ninprc(0)*(nextra+1)
c
c
      do iimp=1,nimp
      do jimp=1,nimp
        lj=nposimp(3,jimp)
        li=nposimp(3,iimp)
c
          if(iimp.ne.jimp) then
            rdiff(1)=(nposimp(1,iimp)-nposimp(1,jimp))*a1(1)+
     >               (nposimp(2,iimp)-nposimp(2,jimp))*a2(1)+
     >               (cvec(li+n0,1)-cvec(lj+n0,1))
c
            rdiff(2)=(nposimp(1,iimp)-nposimp(1,jimp))*a1(2)+
     >               (nposimp(2,iimp)-nposimp(2,jimp))*a2(2)+
     >               (cvec(li+n0,2)-cvec(lj+n0,2))
c
            rdiff(3)=(cvec(li+n0,3)-cvec(lj+n0,3))
c
            rdiff2=dsqrt(rdiff(1)*rdiff(1)+rdiff(2)*rdiff(2)+
     >            rdiff(3)*rdiff(3))
c
c         srdiff(1)=0.d0
c         srdiff(2)=0.d0
c         srdiff(3)=0.d0
c
c         vecna(1)=vecn(1,lj)
c         vecna(2)=vecn(2,lj)
c         vecna(3)=vecn(3,lj)
c
c -----------------------------------------------
c inverse of the rotation into the magnetic field
c
c       call rotvec(rdiff,srdiff,vecna,phi(ljj),-1)
c -----------------------------------------------
c
c       call spherh(madmax,srdiff(1),srdiff(2),srdiff(3),ylm) !original version
        call spherh(madmax,rdiff(1),rdiff(2),rdiff(3),ylm) 
c
        ii=0
        do ll=0,madmax
           do mm=-ll,ll
           ii=ii+1 
              mlmij(ii,iimp,jimp)=
     >                 1/(rdiff2**(ll+1))*fac*dconjg(ylm(mm,ll))
           end do
        end do
c
        end if
      end do
      end do
c
      return
      end
