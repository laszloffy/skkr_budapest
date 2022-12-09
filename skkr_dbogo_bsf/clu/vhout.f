c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine vhout(nimp,vmadih,v0,iscreen,vscreen,kpairind,taupath)
c
c =======================================
      implicit none
      include '../param.h'
c
      real*8 vmadih(mimp)
      real*8 v0
      real*8 vscreen
c
      character*50 vh_file
      character*300 taupath
c
      integer kpairind(mimp,mimp)
      integer nimp
      integer iimp
      integer jimp
      integer ij
      integer iscreen
c
      integer kvhout
      parameter(kvhout=32)
c =======================================
c
      ij=index(taupath,' ')-1
      vh_file=taupath(1:ij)//'mad.vh'
c
c     Open binary file for the madelung pot. from the host
c
      open(UNIT=kvhout,FILE=vh_file,FORM='unformatted',STATUS='unknown')
c
      write(kvhout) (vmadih(iimp),iimp=1,nimp)
      write(kvhout) ((kpairind(iimp,jimp),iimp=1,nimp),jimp=1,nimp)
      write(kvhout) v0
      write(kvhout) vscreen
      write(kvhout) iscreen
c
      close(kvhout)
c
      return
      end
