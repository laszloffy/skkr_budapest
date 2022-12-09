c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine dtat(tclinv,tiinv,nposimp,lmax,tdim,deltat,nimp,nintfc)
c
c ***********************************************
c * Calculate the delta-t matrices              *
c *-------------------------------------------- *
c *                                             *
c *  input:  tclinv - inverse of the t-matrices *
c *                   of the impuritiies        *
c *          tiinv  - inverse of the t-matrices *
c *                   of the host (layers ind)  *
c *          lmax   - lmax                      *
c *                                             *
c *  output: deltat - the difference of         *
c *                   tclinv and tiinv for each *
c *                   impurities                *
c *                                             *
c ***********************************************
c
c ================================================
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
      integer iimp,lay
      integer nposimp(3,nimp)
c
      integer k1,k2
      integer lmax
      integer tdim
c
      complex*16 tclinv(tdim,tdim,nimp)
      complex*16 tiinv(tdim,tdim,nintfc)
      complex*16 deltat(tdim,tdim,nimp)
      data tol/1.0d-25/
c
      nl=lmax+1
      nl2=nl*nl
      kmax=2*lmax+1
      kmymax=2*nl2
c
      do iimp=1,nimp
        lay=nposimp(3,iimp)
        do k1=1,tdim
        do k2=1,tdim
c          write(6,*) "k1,k2,iimp,lay,tdim: ",k1,k2,iimp,lay,tdim
           deltat(k1,k2,iimp)=tiinv(k1,k2,lay)-tclinv(k1,k2,iimp)
        end do
        end do
      end do
c
      return
      end
