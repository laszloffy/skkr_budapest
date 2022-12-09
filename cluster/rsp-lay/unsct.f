c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine unsct(tmin,tmout,kmymax,alphakkr,li)
c ======================================================================
c
c     implicit real*8 (a-h,o-z)
      implicit none
c
      include '../param.h'
c
      logical unsctest
c
      complex*16 tmin(kmymaxp,kmymaxp)  !screened tm1 matrix
      complex*16 tmout(kmymaxp,kmymaxp) !unscreened tm1 matrix
      complex*16 ttmp(kmymaxp,kmymaxp)
      complex*16 detl
c
      complex*16 alphakkr(0:lmaxp)
c     complex*16 alphalkkr(0:lmaxp,minprc)
c     complex*16 alpharkkr(0:lmaxp,minprc)
c     complex*16 alphaintkkr(0:lmaxp,mintfc)
c     common/scrpar/alphalkkr,alpharkkr,alphaintkkr
c
      integer kmymax
      integer kmy
      integer li
      integer ll

      integer itest
      common/test/itest

      real*8 tol
      data tol/1.0d-15/
c
      integer ldex(50)
      data ldex/0,0,
     >          1,1,1,1,1,1,
     >          2,2,2,2,2,2,2,2,2,2,
     >          3,3,3,3,3,3,3,3,3,3,3,3,3,3,
     >          4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4/
c
c ======================================================================
c
c DEBUG
      unsctest=.false.
c DEBUG
c
c inverse screened t-matrix
      if(itest.gt.3) then
         write(6,*) '  screened tm1 - matrix',li
         call outmat1(tmin,kmymax,kmymax,kmymaxp,tol,6)
      end if
c
c screened t-matrix
      call repl(ttmp,tmin,kmymax,kmymaxp)
      call gjinv(ttmp,kmymax,kmymaxp,detl)
c
      if(itest.gt.3) then
         write(6,*) '  screened t - matrix',li
         call outmat1(ttmp,kmymax,kmymax,kmymaxp,tol,6)
      end if
c
c unscreening-transformation on t-matrix
      call repl(tmout,ttmp,kmymax,kmymaxp)
      do kmy=1,kmymax
         ll=ldex(kmy)
         tmout(kmy,kmy)=tmout(kmy,kmy)+alphakkr(ll)
      end do
c
      if(itest.gt.3) then
         write(6,*) '  physical t - matrix',li
         call outmat1(tmout,kmymax,kmymax,kmymaxp,tol,6)
      end if
c
c  unscreened tm1 matrix'
      call gjinv(tmout,kmymax,kmymaxp,detl)
c
      if((itest.gt.3).OR.unsctest) then
         write(6,*) '<unsct>:  unscreened tm1 matrix',li
         call outmat1(tmout,kmymax,kmymax,kmymaxp,tol,6)
      end if
c
      return
      end
