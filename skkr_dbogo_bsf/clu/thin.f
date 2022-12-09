c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine thin(ie,ze,nimp,tdim,nintfc,tiinv,taupath)
c
c ***********************************************
c * read in  (inverse of) t-matrices of host    *
c *---------------------------------------------*
c * input: ze - energy points                   *
c *        ne - # of energy points              *
c *        nintfc - # of layers                 *
c *        n - tdim - nbogomax                  *
c *                                             *
c * output:                                     *
c *        tiinv - inverse t-matrices of layers *
c *                                             *
c ***********************************************
c
c =============================================================
      implicit none
c
      include '../param.h'
c
      character*4  temp
      character*50 thin_file
      character*300 taupath 
c
      complex*16 ce
      complex*16 ze
c
      real*8 tol
c
      integer tdim
      integer n
      integer nintfc
      integer ktin
      integer ie,ieread
      integer k1
      integer k2
      integer li
      integer li0
      integer ij
      integer nimp
c
      logical ttest
c 
      complex*16 tiinv(tdim,tdim,nintfc)
c
      parameter(ktin=32)
c
      data tol/1.0d-15/ 
c =============================================================
c
      ttest=.false.
c
      write(temp,'(i4.4)') ie
      ij=index(taupath,' ')-1
      thin_file=taupath(1:ij)//temp//'.tm'
c
c     Open binary file for the t-matrices of the host
c
      open(ktin,FILE=thin_file,FORM='unformatted',STATUS='old')
c
      read(ktin) ieread
      read(ktin) ce
c
      if(zabs(ze-ce).gt.tol)  then 
          write(6,*) '<thin>: STOP: energy mismatch'
          stop '<thin>: STOP: energy mismatch'
      end if
c
      n=tdim
c
      do li=1,nintfc
           read(ktin) li0
           read(ktin) ((tiinv(k1,k2,li),k1=1,n),k2=1,n)
c
           if(ttest) then
             write(6,'(/'' tminv host'')')
             write(6,*) 'ie=',ie
             write(6,*) 'li=',li0
             call outmat1(tiinv(1,1,li),n,n,n,tol,6)
           end if
c
      end do
c
      close(ktin)
c
      return
      end
