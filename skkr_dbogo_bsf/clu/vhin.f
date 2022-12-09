c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine vhin(nimp,vmadih,v0,iscreen,vscreen,kpairind,taupath)
c
c ==============================================================
      implicit none
c
      include '../param.h'
c
      real*8 vmadih(nimp)
      real*8 v0
      real*8 vscreen0
      real*8 vscreen
c
      character*50 vh_file
      character*300 taupath
c
      integer kpairind(nimp,nimp)
      integer ij
      integer iimp
      integer jimp
      integer nimp
      integer kvhin
c
      integer iscreen0
      integer iscreen
c
      integer itest
c
      parameter(kvhin=32)
      common/test/itest
c ==============================================================
c
      ij=index(taupath,' ')-1
      vh_file=taupath(1:ij)//'mad.vh'
c
c     Open binary file for the madelung pot. from the host
c
      open(UNIT=kvhin,FILE=vh_file,FORM='unformatted',STATUS='old')
c 
      read(kvhin) (vmadih(iimp),iimp=1,nimp)
      read(kvhin) ((kpairind(iimp,jimp),iimp=1,nimp),jimp=1,nimp)
      read(kvhin) v0
      read(kvhin) vscreen0
      read(kvhin) iscreen0
c
      if(vscreen0.ne.vscreen) then
        write(6,*) '<vhin>: WARNING: vscreen0.ne.vscreen'
        write(6,*) '<vhin>: vscreen modified to vscreen=vscreen0'
        vscreen=vscreen0
      end if
c
      iscreen=iscreen0
      if(iscreen0.ne.0) then
        write(6,*) '<vhin>: CLUSTER CALCULATION IN SCREENED REPR.'
        write(6,*) '<vhin>: vscreen=',vscreen 
      else
       if(itest.ge.2) then
        write(6,*) '<vhin>: CLUSTER CALCULATION IN PHYSICAL REPR.'
       end if
      end if
c
      if(iscreen0.ne.0) then
        write(6,*) '<vhin>: STOP:: LUSTER CALCULATION IN SCREENED REPR.'
        write(6,*) '<vhin>:        NOT YET IMPLEMENTED'
        write(6,*) '<vhin>:  iscreen=',iscreen
        stop
      end if
c
      close(kvhin)
c
      return
      end
