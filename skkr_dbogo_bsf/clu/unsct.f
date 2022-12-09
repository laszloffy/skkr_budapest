c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine unsct(tmin,tmout,tdim,alphakkr,alphakkrh,li)
c ======================================================================
c
c     implicit real*8 (a-h,o-z)
      implicit none
c
      include '../param.h'
c
      logical unsctest
c
      complex*16 tmin(dbogomaxp,dbogomaxp)  !screened tm1 matrix
      integer tdim
      complex*16 tmout(tdim,tdim) !unscreened tm1 matrix
      complex*16 ttmp(dbogomaxp,dbogomaxp)
      complex*16 detl
c
      complex*16 alphakkr(0:lmaxp)
      complex*16 alphakkrh(0:lmaxp)
      complex*16 alphal
c     complex*16 alphalkkr(0:lmaxp,minprc)
c     complex*16 alpharkkr(0:lmaxp,minprc)
c     complex*16 alphaintkkr(0:lmaxp,mintfc)
c     common/scrpar/alphalkkr,alpharkkr,alphaintkkr
c
      integer kmy
      integer li
      integer ll

      integer itest
      common/test/itest

      real*8 tol
      data tol/1.0d-15/
c
c      integer ldex(50)
c      data ldex/0,0,
c     >          1,1,1,1,1,1,
c     >          2,2,2,2,2,2,2,2,2,2,
c     >          3,3,3,3,3,3,3,3,3,3,3,3,3,3,
c     >          4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4/
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
         call outmat1(tmin,tdim,tdim,dbogomaxp,tol,6)
      end if
c
c screened t-matrix
      call repl(ttmp,tmin,tdim,dbogomaxp)
      call gjinv(ttmp,tdim,dbogomaxp,detl)
c
      if(itest.gt.3) then
         write(6,*) '  screened t - matrix',li
         call outmat1(ttmp,tdim,tdim,dbogomaxp,tol,6)
      end if
c
c unscreening-transformation on t-matrix
      do kmy=1,tdim
        if (kmy>tdim/2) then
          alphal=alphakkrh(ldex(kmy-tdim/2))
        else
          alphal=alphakkr(ldex(kmy))
        end if
        tmout(kmy,kmy)=ttmp(kmy,kmy)+alphal
      end do
c
      if(itest.gt.3) then
         write(6,*) '  physical t - matrix',li
         call outmat1(tmout,tdim,tdim,dbogomaxp,tol,6)
      end if
c
c  unscreened tm1 matrix'
      call gjinv(tmout,tdim,dbogomaxp,detl)
c
      if((itest.gt.3).OR.unsctest) then
         write(6,*) '<unsct>:  unscreened tm1 matrix',li
         call outmat1(tmout,tdim,tdim,dbogomaxp,tol,6)
      end if
c
      return
      end
