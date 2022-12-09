c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine invdeltau(dttau,nimp,lmax)
c
c ***************************************************
c * Invert the [1-dt(E)**(-1)*tauhost(E)]           *
c * matirx, therefore makes it reindexed            *
c * kappa1,kappa2,iimp,jimp -> i,j                  *
c * at the end:                                     *
c * i,j ->  kappa1,kappa2,iimp,jimp                 *
c *-------------------------------------------------*
c * input:  dttau - [1-dt(E)**(-1)*tauhost(E)]      *
c *         nimp  - # of impurities                 *
c *         lmax  - lmax                            *
c *                                                 *
c * output: dttau: [1-dt(E)**(-1)*tauhost(E)]**(-1) *
c *-------------------------------------------------*
c * calls: gjinv                                    *
c *                                                 *
c ***************************************************
c
c ===================================================
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
c
      integer mdimimp
      parameter (mdimimp=mimp*kmymaxp)
      complex*16 dttau(kmymaxp,kmymaxp,mimp,mimp)
c ----- Allocatable arrays ----------------------------------------
c     complex*16 help(mdimimp,mdimimp)
      integer    AllocateStat
      complex*16, allocatable :: help(:,:)
c
      complex*16 detl
c
      integer nimp
      integer lmax
      integer big
c
      integer iimp,jimp
      integer k1,k2
c
      logical tautest,local
c
      data tol/1.0d-15/
c
c ====================================================
c Allocata memories ---------------------------------------------------------
      allocate(help(mdimimp,mdimimp), stat = AllocateStat)
      call alloccheck(AllocateStat,'help in invdeltau')
c
      tautest=.false.
      local=.false.
c
      nl=lmax+1
      nl2=nl*nl
      kmax=2*lmax+1
      kmymax=2*nl2
      mdim=mimp*kmymaxp
c
      big=nimp*kmymax
c              
      do iimp=1,nimp
      do jimp=1,nimp
         do k1=1,kmymax
         do k2=1,kmymax      
           li=(iimp-1)*kmymaxp
           lj=(jimp-1)*kmymaxp
           help(li+k1,lj+k2)=dttau(k1,k2,iimp,jimp)
         end do
         end do
      end do
      end do
c
c -----------------------Test------------------------------
c     diff2=0.d0
c       do k1=1,kmymax
c       do k2=1,kmymax
c         diff=cdabs(dttau(k1,k2,1,1)-dttau(k1,k2,2,2))
c         if(diff.ge.diff2) diff2=diff
c       end do
c       end do
c      write(6,'(a,d35.30)') 'Tauproblem:',diff2
c ----------------------------------------------------------
         
c     ----------------------------
      call gjinv(help,big,mdim,detl)
c     ----------------------------
      do iimp=1,nimp
      do jimp=1,nimp
         do k1=1,kmymax
         do k2=1,kmymax
           li=(iimp-1)*kmymaxp
           lj=(jimp-1)*kmymaxp
           dttau(k1,k2,iimp,jimp)=help(li+k1,lj+k2) 
         end do
         end do
c
        if(tautest.OR.local) then
          write(6,*) 'invdeltau'
          write(6,*) 'IIMP=',iimp,'  JIMP=',jimp
          call outmat1(dttau(1,1,iimp,jimp),kmymax,kmymax,kmymaxp,tol,6)
        end if
c 
      end do
      end do 
c
      deallocate(help,stat=AllocateStat)
      call alloccheck(AllocateStat,'help dealloc in invdeltau')
c
c -----------------------Test------------------------------
c     diff2=0.d0
c       do k1=1,kmymax
c       do k2=1,kmymax
c         diff=cdabs(dttau(k1,k2,1,1)-dttau(k1,k2,2,2))
c         if(diff.ge.diff2) diff2=diff
c       end do
c       end do
c      write(6,'(a,d35.30)') 'InvdelTauproblem:',diff2
c ----------------------------------------------------------
c
      return
      end
