c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine tauhin(ie,ze,npair,npair1,kpair,nimp,tdim,nintfc,
     >                  tau_ij,tau_ji,tau_ii,taupath)
c
c ***********************************************
c * read in the tau matrix of the cluster       *
c *---------------------------------------------*
c * inpur:  ze - energy array                   *
c *         ne - # of energy points             *
c *         nimp - # of impurities              *
c *         n - tdim - nbogomax                 *
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
      character*50 tauhin_file
      character*300 taupath
c
      integer ie
      integer ieread
      integer ipair
      integer jpair
      integer kpair(4,npair1)
      integer kpairin(4)
      integer ij
      integer li
      integer lj
c     integer jimp
      integer k1
      integer k2
      integer ne
      integer nk
      integer n
c
      integer nintfc
      integer npair,npair1,npairin
      integer nimp,nimp1
      integer ktauin
      integer tdim
c
      integer itest
c
      parameter(ktauin=32)
c
      real*8 tol
c
      complex*16 ce
      complex*16 ze
      complex*16 tau_ij(tdim,tdim,npair1)
      complex*16 tau_ji(tdim,tdim,npair1)
      complex*16 tau_ijin(tdim,tdim)
      complex*16 tau_jiin(tdim,tdim)
      complex*16 tau_ii(tdim,tdim,nintfc)
c
      logical tautest
c
c
      data tol/1.0d-15/
c
      common/test/itest
c ======================================================
c
      tautest=.false. 
c
      write(temp,'(i4.4)') ie
      ij=index(taupath,' ')-1
      tauhin_file=taupath(1:ij)//temp//'.tau'
      if(itest.ge.2) then
        write(6,*) '<tauhin>: the tau read from ', tauhin_file
      end if
c
c     Open binary file for the tau-matrices of the host
c
      open(UNIT=ktauin,FILE=tauhin_file,
     >     FORM='unformatted',STATUS='old')
c
      read(ktauin) ieread
      read(ktauin) ce
      read(ktauin) nk
      read(ktauin) npairin
      if(zabs(ze-ce).gt.tol) then
          write(6,*) '<tauhin> : STOP: energy mismatch',ie,ze,ce
           stop '<tauhin> : energy mismatch'
      end if
c
      n=tdim
c
      if(npair.gt.0) then
        ipair=0
        do while (ipair.lt.npair.or.jpair.lt.npairin)
          read(ktauin) jpair
c          write(6,*) "jpair=",jpair
          read(ktauin) (kpairin(k1),k1=1,4)
          read(ktauin) ((tau_ijin(k1,k2),k1=1,n),k2=1,n)
          read(ktauin) ((tau_jiin(k1,k2),k1=1,n),k2=1,n)
          do li=1,npair
         if(kpair(1,li).eq.kpairin(1).and.kpair(2,li).eq.kpairin(2).and.
     >   kpair(3,li).eq.kpairin(3).and.kpair(4,li).eq.kpairin(4)) then
c              write(6,*) '<tauhin>',li,kpair(:,li)
c             write(6,*) '<tauhin> : kpairin=',kpairin
              tau_ij(:,:,li)=tau_ijin(:,:)
              tau_ji(:,:,li)=tau_jiin(:,:)
              ipair=ipair+1
              exit
           elseif(kpair(1,li).eq.-kpairin(1).and.
     >             kpair(2,li).eq.-kpairin(2).and.
     >             kpair(3,li).eq.kpairin(4).and.
     >             kpair(4,li).eq.kpairin(3)) then
c              write(6,*) '<tauhin>',li,kpair(:,li)
c             write(6,*) '<tauhin> : kpairin=',kpairin
              tau_ij(:,:,li)=tau_jiin(:,:)
              tau_ji(:,:,li)=tau_ijin(:,:)
              ipair=ipair+1
              exit
           end if
          end do
c -TEST BEGINN
          if(tautest) then
            write(6,*) 'Diagonal elements: tau_ij ',ipair
            write(6,*) 'ipair=',ipair
            call outmat1(tau_ij(1,1,ipair),n,n,n,tol,6)
            write(6,*) 'Diagonal elements: tau_ji ',ipair
            write(6,*) 'ipair=',ipair
            call outmat1(tau_ji(1,1,ipair),n,n,n,tol,6)
          end if
c -TEST END
c        write(6,*) "ipair,npair",ipair,npair
c        write(6,*) "ipair<npair",ipair.lt.npair
        end do
        if(ipair.ne.npair.and.jpair.ne.npairin) then
          write(6,*) '<tauhin> : STOP: ',ipair,npair,jpair,npairin
          stop '<tauhin> : pair mismatch'
        end if
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
          call outmat1(tau_ii(1,1,li),n,n,n,tol,6)
        end if
c -TEST END
c
      end do
c
      close(ktauin)
c
      return
      end
