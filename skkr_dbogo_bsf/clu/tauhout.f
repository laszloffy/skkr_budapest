c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine tauhout(ie,ce,npair,npair1,kpair,nintfc,nk,tdim,
     >                   tau_ij,tau_ji,tau_ii,taupath)
c
c ***********************************************
c * write out (inverse of) t-matrices of host   *
c *---------------------------------------------*
c * input: ce - energy                          *
c *        npair - # of pairs                   *
c *        kmymax                               *
c *        tau_ij - off-diag tau-matrix         *
c *        nk - # of k-points                   *
c *                                             *
c *                                             *
c ***********************************************
c
c ===============================================
c      implicit none
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
      character*4  temp
      character*50 tauout_file
      character*300 taupath
c
      integer ktauout
      integer tdim
      integer nintfc
      integer npair,npair1
      integer ipair
      integer ij
      integer nk
      integer li
      integer ie
      integer k1
      integer k2
      integer kpair(4,npair1)
c      integer n
c
c      integer ndim
c      parameter(ndim=kmymaxp)
c
      complex*16 ce
      complex*16 tau_ij(tdim,tdim,npair1)
      complex*16 tau_ji(tdim,tdim,npair1)
      complex*16 tau_ii(tdim,tdim,nintfc)
c
      logical tautest
      parameter(ktauout=32)
c
      real*8 tol
      data tol/1.0d-15/
c ===============================================
c
      write(6,*) ' ie=', ie, 'in tauhout'
c     call flush(6)
      write(temp,'(i4.4)') ie
      ij=index(taupath,' ')-1
      tauout_file=taupath(1:ij)//temp//'.tau'
      write(6,*) ' Unscreened tau is written', tauout_file
c
c     Open binary file for the tau-matrices of the host
c
      open(UNIT=ktauout,FILE=tauout_file,
     >     FORM='unformatted',STATUS='unknown')
c
      write(ktauout) ie
      write(ktauout) ce
      write(ktauout) nk
      write(ktauout) npair
c
      n=tdim
c
      if(npair.gt.0) then
        do ipair=1,npair
          write(ktauout) ipair
          write(ktauout) (kpair(k1,ipair),k1=1,4)
          write(ktauout) ((tau_ij(k1,k2,ipair),k1=1,n),k2=1,n)
          write(ktauout) ((tau_ji(k1,k2,ipair),k1=1,n),k2=1,n)
c -TEST BEGINN
          tautest=.false.
          if(tautest) then
            write(6,*) 'Off-Diagonal elements: tau_ij ',ipair
            write(6,*) 'ipair=',ipair
            call outmat1(tau_ij(1,1,ipair),n,n,n,tol,6)
            write(6,*) 'Off-Diagonal elements: tau_ji ',ipair
            write(6,*) 'ipair=',ipair
            call outmat1(tau_ji(1,1,ipair),n,n,n,tol,6)
          end if
c -TEST END
        end do
      end if
c
      do li=1,nintfc
        write(ktauout) li
        write(ktauout) ((tau_ii(k1,k2,li),k1=1,n),k2=1,n)
      end do
c
      close(ktauout)
c
      return
      end
