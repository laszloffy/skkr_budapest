c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
       subroutine thout(ie,ce,nintfc,nbogomax,tmatz,taupath)
c
c ***********************************************
c * write out (inverse of) effective t-matrices *
c * NON-SCREENED!!!!                            *
c *---------------------------------------------*
c * input: ce - energy                          *
c *        nintfc - # of layers                 *
c *        tminv - inverse t-matrices of layers *
c *                                             *
c *                                             *
c *                                             *
c ***********************************************
c
c ======================================================
      implicit none
c      implicit real*8 (a-h,o-z)
      include '../param.h'
c
c
      character*4  temp
      character*50 tout_file
      character*300 taupath
c
c     integer ndim
c     parameter(ndim=kmymaxp)
c      integer tdim
      integer nbogomax,nintfc
c
      complex*16 tmatz(nbogomax,nbogomax,nintfc)
      complex*16 ce
c
      integer*4 ie
      integer ij
      integer li
      integer k1
      integer k2
      integer n
c
      integer ktout
      parameter(ktout=32)
c
c ======================================================
c
      write(6,*) 'In thout'
      call flush(6)
      write(6,*) 'ce: ', ce
      call flush(6)
      write(6,*) 'ie: ', ie
      call flush(6)
      write(temp,'(i4.4)') ie
      ij=index(taupath,' ')-1
      tout_file=taupath(1:ij)//temp//'.tm'
      write(6,*) ' Unscreened cpa tm is written',ie,tout_file
      call flush(6)
c
c     Open binary file for the t-matrices of the host
c
      open(ktout,FILE=tout_file,
     >     FORM='unformatted',STATUS='unknown')
c
      write(ktout) ie
      write(ktout) ce
c
      n=nbogomax
c
      call flush(ktout)
c
      do li=1,nintfc
        write(ktout) li
        write(ktout) ((tmatz(k1,k2,li),k1=1,n),k2=1,n)
      end do 
c
      close(ktout)
c
      return
      end 
