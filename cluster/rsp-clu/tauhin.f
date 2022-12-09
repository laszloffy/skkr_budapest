c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine tauhin(ie,ze,npair,nintfc,kmymax,
     >                  tau_ij,tau_ji,tau_ii,taupath,imesh)
c
c ***********************************************
c * read in the tau matrix of the cluster       *
c *---------------------------------------------*
c * inpur:  ze - energy array                   *
c *         ne - # of energy points             *
c *         nimp - # of impurities              *
c *         n - kmymax                          *
c *                                             *
c * output:                                     *
c *         tauij - off-diag tau-matrix         *
c *         nk - number of k-points             *
c *                                             *
c ***********************************************
c
c ======================================================
      implicit none
      include '../param.h'
c
      character*4  temp
      character*250 tauhin_file
      character*200 taupath
c
      integer ie,iein
      integer ipair
      integer jpair
      integer ij
      integer li
      integer lj
c     integer jimp
      integer k1
      integer k2
      integer ne
      integer nk
      integer kmymax
      integer n
c
      integer nintfc
      integer npair
      integer nimp
      integer ndim
      integer ktauin
c
c     parameter(ndim=kmymaxp)
      parameter(ktauin=32)
c
      real*8 tol
c
      complex*16 ce
      complex*16 ze
      complex*16 tau_ij(kmymaxp,kmymaxp,mpair)
      complex*16 tau_ji(kmymaxp,kmymaxp,mpair)
      complex*16 tau_ii(kmymaxp,kmymaxp,mintfc)
c     complex*16 tauij(ndim,ndim,mimp,mimp)
c
      logical tautest
      integer imesh
c
c
      data tol/1.0d-15/
c
c ======================================================
c
      tautest=.false.
c
      write(temp,'(i4.4)') ie
      ij=index(taupath,' ')-1
      tauhin_file=taupath(1:ij)//temp//'.tau'
      write(6,*) '<tauhin>: the tau read from ', tauhin_file
c
c     Open binary file for the tau-matrices of the host
c
      open(UNIT=ktauin,FILE=tauhin_file,
     >     FORM='unformatted',STATUS='old')
c
      read(ktauin) iein
      read(ktauin) ce
      read(ktauin) nk
      write(6,*) '<tauhin> imesh:',imesh
c     if(ze.ne.ce) write(6,*) '<tauhin> : STOP: energy mismatch',ie
      if(zabs(ze-ce).gt.tol) then
       if(imesh.ne.0) then
         write(6,*) '<tauhin> : STOP: energy mismatch',ie
         stop '<tauhin> : energy mismatch'
       else
        write(6,'(''Energy mismatch: ie='',i3,'' ce='',2e14.6,
     >            '' ze='',2e14.6)'),ie,ce,ze
        ze = ce
       endif
      endif
c
      n=kmymax
c
      if(npair.gt.0) then
        do ipair=1,npair
          read(ktauin) jpair
          read(ktauin) ((tau_ij(k1,k2,ipair),k1=1,n),k2=1,n)
          read(ktauin) ((tau_ji(k1,k2,ipair),k1=1,n),k2=1,n)
c -TEST BEGINN
          if(tautest) then
            write(6,*) 'Tau elements: tau_ij ',ipair
            write(6,*) 'ipair=',ipair
            call outmat1(tau_ij(1,1,ipair),n,n,kmymaxp,tol,6)
            write(6,*) 'Tau elements: tau_ji ',ipair
            write(6,*) 'ipair=',ipair
            call outmat1(tau_ji(1,1,ipair),n,n,kmymaxp,tol,6)
          end if
c -TEST END
        end do
      end if
c
      do li=1,nintfc
        read(ktauin) lj
        read(ktauin) ((tau_ii(k1,k2,li),k1=1,n),k2=1,n)
c
c -TEST BEGINN
        if(tautest) then
          write(6,*) 'Diagonal elements: tau_ii',li
          write(6,*) 'li=',li,'  ie=',ie
          call outmat1(tau_ii(1,1,li),n,n,kmymaxp,tol,6)
        end if
c -TEST END
c
      end do
c
      close(ktauin)
c
      return
      end
