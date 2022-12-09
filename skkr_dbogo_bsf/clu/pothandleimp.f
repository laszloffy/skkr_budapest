c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
c***********************************************************************
c Now, the radial mesh of the potential within a layer is chosen as the 
c finest grid between both CPA species.
c
c For bulk calculation, assumed: ninprcl=ninprcr=ninprc(1)  &  nprc=1 
c***********************************************************************
      subroutine pothandleimp(
     > lmax,nimp,v00,ivacpot,opot,imppot,dx,ns,rs,idpotimp,
     > vrimp,brimp,boprimp,zimp,ib0,b0,impb0,impb0f,isigb,
     > deltaimp,rlambdaimp,singratimp,uratimp,dratimp)
c
c -read Left,Right and Layer potentials
c bnyari: TODO chech these dummy variables 
c                 v00,ivacpot 
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c     integer mmax
c     parameter (mmax=nimp)
c
      logical opot,potin
      character*10 idpot0
      character*10 idpotimp(nimp)
      character*30 imppot
c
      real*8 vrimp(nrad,nimp)
      real*8 brimp(nrad,nimp)
      real*8 boprimp(nrad,2,nimp)
      complex*16 deltaimp(nrad,nimp)
      dimension rlambdaimp(nimp)
      dimension singratimp(nimp),uratimp(nimp),dratimp(nimp)
      real*8 dx(nimp)
      real*8 rs(nimp)
      real*8 zimp(nimp)
c
      integer isigb(nimp)
      integer ns(nimp)
c
c      complex*16 qmom(lmsup)
c
      real*8 potshift
      real*8 b0
c
      integer impb0
      integer impb0f
      integer ib0
      integer imom
c
      data tiny/1.0d-6/
      common/test/itest
c
      lmaxs=2*lmax
      lmmaxs=(lmaxs+1)*(lmaxs+1)
c
      imom=0
c -- rsl/~r is the Wigner-Seitz radius for the left and right bulk and
c    rs(i),i=1,nintfc is the WS-radius for the interface !!
c
      open(20,file=imppot,status='old')
      read(20,*) potin
      write(6,*) 'imppot=',imppot
      write(6,*) 'potin=',potin
      do li=1,nimp
c        write(6,*) ''
        call readpot(idpot0,conc0,vrimp(1,li),brimp(1,li),
     >       boprimp(1,1,li),potin,0,opot,dx(li),ns(li),rs(li),zimp(li),
     >       potzero,deltaimp(1,li),rlambdaimp(li),
     >       singratimp(li),uratimp(li),dratimp(li))
c        write(6,*)  '<pothandleimp> vrimp(:,li)'
c        write(6,*) vrimp(:,li)
c        write(6,*)  '<pothandleimp> brimp(:,li)'
c        write(6,*) brimp(:,li)
c
        idpotimp(li)=idpot0
c
        if((impb0.eq.0).or.((li.ge.impb0).and.(li.le.impb0f))) then
          if(ib0.eq.0) then
            x=dlog(rs(li))-(ns(li)-1)*dx(li)
            do i=1,ns(li)
              r=dexp(x)
              brimp(i,li)=brimp(i,li)+r*b0
              x=x+dx(li)
            end do
          else
            x=dlog(rs(li))-(ns(li)-1)*dx(li)
            do i=1,ns(li)
              r=dexp(x)
              brimp(i,li)=r*b0
              x=x+dx(li)
            end do
          end if
        end if
c
        if(isigb(li).ne.1) then
          do i=1,ns(li)
            brimp(i,li)=-brimp(i,li)
          end do
          write(6,'(''<pothandleimp> : Sign of effective 
     >          field changed'')') li
        end if
c
        write(6,'(i5,2x,a12)') li,idpotimp(li)
c
      end do
      close(20)
c
      if(itest.ge.2) then
      write(6,'(''<pothandleimp> : pothandleimp processed!'')')
      end if
c
c
      return
      end
