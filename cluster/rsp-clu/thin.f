c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine thin(ie,ze,nintfc,kmymax,tiinv,taupath)
c
c ***********************************************
c * read in  (inverse of) t-matrices of host    *
c *---------------------------------------------*
c * input: ze - energy points                   *
c *        ne - # of energy points              *
c *        nintfc - # of layers                 *
c *        n - kmymax                           *
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
      character*250 thin_file
      character*200 taupath 
c
      complex*16 ce
      complex*16 ze
c
      real*8 tol
c
      integer kmymax
      integer n
      integer nintfc
      integer ndim
      integer ktin
      integer ie,iein
      integer k1
      integer k2
      integer li
      integer li0
      integer ij
c
      logical ttest
c 
      parameter(ndim=kmymaxp)
      complex*16 tiinv(ndim,ndim,mintfc)
c
      parameter(ktin=32)
c
      data tol/1.0d-10/ 
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
      read(ktin) iein
      read(ktin) ce
c
      if(zabs(ze-ce).gt.tol) then
        write(6,*) '<thin>: STOP: energy mismatch',ze,ce
        stop '<thin>: STOP: energy mismatch'
      end if
c
      n=kmymax
c
      do li=1,nintfc
           read(ktin) li0
           read(ktin) ((tiinv(k1,k2,li),k1=1,n),k2=1,n)
c
           if(ttest) then
             write(6,'(/'' tminv host'')')
             write(6,*) 'ie=',ie
             write(6,*) 'li=',li0
             call outmat1(tiinv(1,1,li),n,n,kmymaxp,tol,6)
           end if
c
      end do
c
      close(ktin)
c
      return
      end
