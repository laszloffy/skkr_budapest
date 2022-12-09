c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine dtat(tclinv,tiinv,nposimp,lmax,deltat,nimp)
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
      integer nposimp(3,mimp)
c
      integer k1,k2
      integer lmax
c
      complex*16 tclinv(kmymaxp,kmymaxp,mimp)
      complex*16 tiinv(kmymaxp,kmymaxp,mintfc)
      complex*16 deltat(kmymaxp,kmymaxp,mimp)
      complex*16 detl
c
      logical local
c
c     common/test/itest
c
      data tol/1.0d-25/
c ================================================
c
      local=.false. 
c
      nl=lmax+1
      nl2=nl*nl
      kmax=2*lmax+1
      kmymax=2*nl2
c
      do iimp=1,nimp
        lay=nposimp(3,iimp)
c
        if(local) then
          write(6,*) '  unscreened m - matrix of host'
          write(6,*) 'IMP=',iimp
          call outmat1(tiinv(1,1,lay),kmymax,kmymax,kmymaxp,tol,6)
          write(6,*) 'IMP=',iimp
          write(6,*) '  unscreened m - matrix of imp'
          call outmat1(tclinv(1,1,iimp),kmymax,kmymax,kmymaxp,tol,6)
        end if
c
        do k1=1,kmymax
        do k2=1,kmymax
           deltat(k1,k2,iimp)=tiinv(k1,k2,lay)-tclinv(k1,k2,iimp)
        end do
        end do
c
        if(local) then
          write(6,*) 'dtat'
          write(6,*) 'IMP=',iimp
c         do k1=1,kmymax
c         do k2=1,kmymax
c            write(6,'(2d35.30)') deltat(k1,k2,iimp)
c         end do
c         end do
c
          call outmat1(deltat(1,1,iimp),kmymax,kmymax,kmymaxp,tol,6)
        end if
      end do
c
      return
      end
