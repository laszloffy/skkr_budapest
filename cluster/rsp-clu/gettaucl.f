c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine gettaucl1(tauij,deltat,nimp,lmax,tauclij,taucl,detl)
c
c **************************************************************************     
c * construct tau for the embedded cluster
C *  taucl = (tauhost(E)**(-1) - dt(E)**(-1))**(-1)
c  input:
c         nimp     - number of impurities
c         tauij    - tauhost
c         delatat  - dt**(-1)
c  output:
c         tauclij    - tau cluster
c **************************************************************************     
      include '../param.h'
c
      logical tautest
      integer mdimimp
      parameter (mdimimp=mimp*kmymaxp)
c
      complex*16 tauij(kmymaxp,kmymaxp,mimp,mimp)
      complex*16 deltat(kmymaxp,kmymaxp,mimp)
      complex*16 tauclij(mdimimp,mdimimp)
      complex*16 taucl(kmymaxp,kmymaxp,mimp)
c     complex*16 help(mdimimp,mdimimp)
      complex*16 detl
      complex*16 sum, zpi
      integer    big, ipiv(mdimimp)
      integer    AllocateStat
c ----- Allocatable arrays ----------------------------------------
c     complex*16 help(mdimimp,mdimimp)
      complex*16, allocatable :: help(:,:)
c
      parameter (zpi = (0.d0,3.14159265359d0))
c
      real*8 tol
      data tol/1.0d-8/
c
      nl=lmax+1
      nl2=nl*nl
      kmax=2*lmax+1
      kmymax=2*nl2
      big = kmymax*nimp
c Allocata memories ---------------------------------------------------------
      allocate(help(mdimimp,mdimimp), stat = AllocateStat)
      call alloccheck(AllocateStat,'help in gettaucl')
c
      tautest=.false.
      if(tautest) then
       do iimp=1,nimp
       do jimp=1,nimp
        write(6,*) ' <gettaucl> : tautest, tauij iimp, jimp=',iimp,jimp
        call outmat1(tauij(1:kmymax,1:kmymax,iimp,jimp),
     >               kmymax,kmymax,kmymax,
     >               tol,6)
       end do
       end do
       do iimp=1,nimp
        write(6,*) ' <gettaucl> : tautest, deltat iimp=',iimp
        call outmat1(deltat(1:kmymax,1:kmymax,iimp),
     >               kmymax,kmymax,kmymax,
     >               tol,6)
       end do
      end if
c------------------------------------------------------------------------------------
c   tau = (tau_h*delta t^{-1} + 1)^{-1} tau_h 
c------------------------------------------------------------------------------------
      do iimp=1,nimp
      do jimp=1,nimp
         do k1=1,kmymax
         do k2=1,kmymax
           li=(iimp-1)*kmymaxp
           lj=(jimp-1)*kmymaxp
           sum = (0.d0,0.d0)
           do k3 = 1,kmymax
            sum = sum - tauij(k1,k3,iimp,jimp)*deltat(k3,k2,jimp)
           end do
           help(li+k1,lj+k2) = sum
           tauclij(li+k1,lj+k2) = tauij(k1,k2,iimp,jimp)
         end do
         end do
      end do
      end do
c
      do i = 1,big
       help(i,i) = help(i,i) + (1.d0,0.d0)
      end do
c
      if(tautest) then
       do iimp=1,nimp
       do jimp=1,nimp
        write(6,*) ' <gettaucl> : tautest, tauclij iimp, jimp=',
     >             iimp,jimp
        li=(iimp-1)*kmymaxp
        lj=(jimp-1)*kmymaxp
        call outmat1(tauclij(li+1:li+kmymax,lj+1:lj+kmymax),
     >               kmymax,kmymax,kmymax,
     >               tol,6)
       end do
       end do
      end if
c --------------------------------------------------------------
      call ZGETRF(big,big,help,mdimimp,IPIV,INFO)
      detl = (0d0,0.d0)
      isign = 1
      do i=1,big
        if (ipiv(i).ne.i) then
         isign = -isign
        endif
        detl = detl + zlog(help(i,i))
      end do
      detl = detl - zpi
      CALL ZGETRS('N',big,big,help,mdimimp,IPIV,tauclij,mdimimp,INFO)  
c --------------------------------------------------------------
c
      if(tautest) then
       do iimp=1,nimp
       do jimp=1,nimp
        write(6,*) ' <gettaucl> : tautest, help iimp, jimp=',
     >             iimp,jimp
        li=(iimp-1)*kmymaxp
        lj=(jimp-1)*kmymaxp
        call outmat1(help(li+1:li+kmymax,lj+1:lj+kmymax),
     >               kmymax,kmymax,kmymax,
     >               tol,6)
       end do
       end do
      end if
c
      do i = 1,big
       do j = 1,big
        help(i,j) = tauclij(i,j)
       end do
      end do
c
      call packtaucl(tauclij,help,mdimimp,nimp,kmymax)
c
      do iimp=1,nimp
         do k1=1,kmymax
         do k2=1,kmymax
           li=(iimp-1)*kmymaxp
           lj=(iimp-1)*kmymaxp
c          taucl(k1,k2,iimp) = tauclij(li+k1,lj+k2)
           taucl(k1,k2,iimp) = help(li+k1,lj+k2)
         end do
         end do
      end do
c
      deallocate(help,stat=AllocateStat)
      call alloccheck( AllocateStat,'help dealloc in gettaucl')
c
      return
      end
c
      subroutine packtaucl(tauclij,help,md,nimp,kmymax)
      include '../param.h'
c
      complex*16 tauclij(kmymaxp,kmymaxp,mimp,mimp)
      complex*16 help(md,md)
c
      do iimp=1,nimp
      do jimp=1,nimp
         do k1=1,kmymax
         do k2=1,kmymax
           li=(iimp-1)*kmymaxp
           lj=(jimp-1)*kmymaxp
           tauclij(k1,k2,iimp,jimp) = help(li+k1,lj+k2)
         end do
         end do
      end do
      end do
      return
      end
