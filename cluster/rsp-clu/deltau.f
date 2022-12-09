c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine deltau(deltat,tauij,dttau,lmax,nimp)
c
c *************************************************
c * Multiplies the deltat and the tauhost         *
c *-----------------------------------------------*
c * input: deltat - the difference of tclinv and  *
c *                 tiinv for each impurities     *
c *        tauij  - tauhost                       *
c *        lmax   - lmax                          *
c *        nimp   - true. mpurities               *
c *                                               *
c * output: dttau - deltat*tauij                  *
c *                                               *
c *************************************************
c
c =================================================
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
      complex*16 tauij(kmymaxp,kmymaxp,mimp,mimp)
      complex*16 tautmp(kmymaxp,kmymaxp)
      complex*16 dtautmp(kmymaxp,kmymaxp)
      complex*16 dttau(kmymaxp,kmymaxp,mimp,mimp)
      complex*16 deltat(kmymaxp,kmymaxp,mimp)
c
      complex*16 sum
c
      integer nimp,lmax
c
      integer k1,k2
      integer iimp,jimp 
c
      logical tautest
      logical local
c
      data tol/1.0d-20/
c
c =================================================
c
      tautest=.false. 
      local=.false.
      nl=lmax+1
      nl2=nl*nl
      kmax=2*lmax+1
      kmymax=2*nl2
      l2=2*lmax
c
c
c ***************
c *  ZERO dttau *
c ***************
c
      do iimp=1,nimp
      do jimp=1,nimp
      do k1=1,kmymax
      do k2=1,kmymax
         dttau(k1,k2,iimp,jimp)=(0.d0,0.d0)
      end do
      end do
      end do
      end do
c
c *************** 
c ------------------------------Old----------------------------------
c       do iimp=1,nimp
c       do jimp=1,nimp
c          do k1=1,kmymax
c          do k2=1,kmymax
c             sum=(0.d0,0.d0)
c             do k3=1,kmymax
c               sum=sum+deltat(k1,k3,iimp)*tauij(k3,k2,iimp,jimp)
c             end do
c             dttau(k1,k2,iimp,jimp)=-sum
c          end do
c          end do
c       end do
c       end do
c       do iimp=1,nimp
c         do k1=1,kmymax
c           dttau(k1,k1,iimp,iimp)=(1.0d0,0.0d0)+dttau(k1,k2,iimp,jimp)
c         end do
c        if(tautest) then
c         write(6,*) 'deltau'
c         write(6,*) 'IMP=',iimp
c         call outmat1(dttau(1,1,iimp,iimp),kmymax,kmymax,kmymaxp,tol,6)
c        end if
c       end do
c -----------------------------------------------------------------------
c
c -----------------------New-----------------------------
      do iimp=1,nimp
      do jimp=1,nimp
         do k1=1,kmymax
         do k2=1,kmymax
           tautmp(k1,k2)=tauij(k1,k2,iimp,jimp)
         end do
         end do
c
        call doubmt1(deltat(1,1,iimp),tautmp,
     >  dtautmp,kmymax,kmymaxp)
c
        do k1=1,kmymax 
        do k2=1,kmymax 
           dttau(k1,k2,iimp,jimp)=-dtautmp(k1,k2)
        end do
        end do
      end do
      end do
c
      do iimp=1,nimp
         do k1=1,kmymax
         dttau(k1,k1,iimp,iimp)=(1.0d0,0.0d0)+dttau(k1,k1,iimp,iimp) 
         end do
      end do
c ----------------------------------------------------------------------
c
c ------------------------------Test-------------------------------------
      if(tautest.OR.local) then
        do iimp=1,nimp
        do jimp=1,nimp
          write(6,*) 'deltau'
          write(6,*) 'IIMP=',iimp
          write(6,*) 'JIMP=',jimp
          call outmat1(dttau(1,1,iimp,jimp),kmymax,kmymax,kmymaxp,tol,6)
        end do
        end do
      end if
c ----------------------------------------------------------------------
c
c ----------------Test-----------------------------------
c     diff2=0.d0
c     do k1=1,kmymax
c     do k2=1,kmymax
c        diff=cdabs(dttau(k1,k2,1,1)-dttau(k1,k2,2,2))
c        if(diff.gt.diff2) diff2=diff
c     end do
c     end do
c     write(6,*) 'DTauproblem:',diff2
c --------------------------------------------------------
c
       return
       end

