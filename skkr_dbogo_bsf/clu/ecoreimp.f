c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine ecoreimp(nimp,npair,npair1,nintfc,
     >           iscreen,tau_ij,tau_ji,tau_ii,tminvcl,
     >           tminvh,lmax,tdim,nposimp,kpairind,!dmata,dmatpa,
     >           taucl,detl)
c
      implicit none
c
      include '../param.h' 
c
      integer tdim
c
      logical tautest
c
      complex*16 alphaintkkr(0:lmaxp,mintfc)
      complex*16 tminvh(tdim,tdim,nintfc)
      complex*16 tminvcl(tdim,tdim,nimp)
c     complex*16 tminvclscr(dbogomaxp,dbogomaxp,nimp)
c     complex*16 deltat(dbogomaxp,dbogomaxp,nimp)
      complex*16, allocatable :: deltat(:,:,:)
      complex*16 tau_ij(tdim,tdim,npair1)
      complex*16 tau_ji(tdim,tdim,npair1)
      complex*16 tau_ii(tdim,tdim,nintfc)
c     complex*16 tauh(dbogomaxp,dbogomaxp,nimp,nimp)
      complex*16, allocatable :: tauh(:,:,:,:)
c     complex*16 dttau(dbogomaxp,dbogomaxp,nimp,nimp)
      complex*16 taucl(tdim,tdim,nimp)
c     complex*16 tauclij(tdim,tdim,nimp,nimp)
      complex*16, allocatable :: tauclij(:,:,:,:)
c     complex*16 gtaucl(dbogomaxp,dbogomaxp,nimp)
c     complex*16 dmata(dbogomaxp,dbogomaxp,nimp)
c     complex*16 dmatpa(dbogomaxp,dbogomaxp,nimp)
      complex*16 detl
c
      integer itest
      integer nimp
      integer npair,npair1
      integer iscreen
      integer nposimp(3,nimp)
      integer kpairind(nimp,nimp)
      integer li
      integer iimp
      integer jimp
      integer mipair
      integer ipair
      integer k1
      integer k2
      integer kmax
      integer kmymax
      integer lmax
      integer nintfc
      integer l2
      integer nl
      integer nl2
      integer AllocateStatus
c
      common/scrpar/alphaintkkr
      common/test/itest
c
      real*8 tol
      data tol/1.0d-15/
c
      if(itest.ge.2) write(6,*) '<ecoreimp>: BEGINN'
c
      nl=lmax+1
      nl2=nl*nl
      kmax=2*lmax+1
      kmymax=2*nl2
      l2=2*lmax
c
      ALLOCATE ( deltat(tdim,tdim,nimp),
     >                  STAT = AllocateStatus)
      call alloccheck( AllocateStatus,
     >     'deltat in ecoreimp                                ' )
      ALLOCATE ( tauh(tdim,tdim,nimp,nimp),
     >                  STAT = AllocateStatus)
      call alloccheck( AllocateStatus,
     >     'tauh in ecoreimp                                  ' )
      ALLOCATE ( tauclij(tdim,tdim,nimp,nimp),
     >                  STAT = AllocateStatus)
      call alloccheck( AllocateStatus,
     >     'tauclij in ecoreimp                               ' )
c
      do iimp=1,nimp-1
        do jimp=(iimp+1),nimp
          ipair=kpairind(iimp,jimp)
          if(ipair.gt.0) then
            do k1=1,2*kmymax
              do k2=1,2*kmymax
                tauh(k1,k2,iimp,jimp)=tau_ij(k1,k2,ipair)
                tauh(k1,k2,jimp,iimp)=tau_ji(k1,k2,ipair)
              end do
            end do
          else
            mipair=-1*ipair
            do k1=1,2*kmymax
              do k2=1,2*kmymax
                tauh(k1,k2,iimp,jimp)=tau_ji(k1,k2,mipair)
                tauh(k1,k2,jimp,iimp)=tau_ij(k1,k2,mipair)
              end do
            end do
          end if
        end do
      end do 
c
      do iimp=1,nimp
        li=nposimp(3,iimp)
        do k1=1,2*kmymax
          do k2=1,2*kmymax
            tauh(k1,k2,iimp,iimp)=tau_ii(k1,k2,li)
          end do
        end do
      end do
c
      tautest=.false.
      if(tautest) then
       do iimp=1,nimp
        write(6,*) ' <gettaucl> : tautest, tminvcl iimp=',iimp
        call outmat1(tminvcl(1:kmymax,1:kmymax,iimp),
     >               kmymax,kmymax,kmymax,
     >               tol,6)
       end do
       do iimp=1,nintfc
        write(6,*) ' <gettaucl> : tautest, tminvh il=',iimp
        call outmat1(tminvh(1:kmymax,1:kmymax,iimp),
     >               kmymax,kmymax,kmymax,
     >               tol,6)
       end do
      end if
c ------------------------------------------------------------------
c --- Calculate delta-t matrices, dt = thost**-1 - tclu**-1  -------
c ------------------------------------------------------------------
         call dtat(tminvcl,tminvh,nposimp,lmax,tdim,deltat,
     >             nimp,nintfc)
c ------------------------------------------------------------------
c --- Construct taucl ----------------------------------------------
c --- taucl = (tauhost**(-1) - dt**(-1))**(-1) --------------
c ------------------------------------------------------------------
         call gettaucl1(tauh,deltat,nimp,lmax,tdim,tdim*nimp,
     >                  tauclij,taucl,detl)
        tautest=.false.
        if(tautest) then
          write(6,*) 'ecoreimp tminvcl(:,:,1)'
          write(6,*) tminvcl(:,:,1)
          write(6,*) 'ecoreimp tminvh(:,:,1)'
          write(6,*) tminvh(:,:,1)
          write(6,*) 'ecoreimp deltat(:,:,1)'
          write(6,*) deltat(:,:,1)
          write(6,*) 'ecoreimp tauh(:,:,1,1)'
          write(6,*) tauh(:,:,1,1)
          write(6,*) 'ecoreimp taucl(:,:,1)'
          write(6,*) taucl(:,:,1)
        end if
c -------------------------------------------------------------------
c --- loop over layers to transform impurity diagonal tau matrices --
c --- into physical representation and rotate to the local ----------
c --- frame of reference --------------------------------------------
c -------------------------------------------------------------------
c      do iimp=1,nimp
cc    !!! no screening for impurity calculation !!!
c        if(iscreen.ne.0) then
c          li=nposimp(3,iimp)
cc         --------------------------------------------------------
c          call phystau(taucl(1,1,iimp),tminvcl(1,1,iimp),
c     >                 tminvcl(1,1,iimp),alphaintkkr(0,li),lmax,1)
cc         --------------------------------------------------------
c        end if
cc rotate tau-matrix to local frame of reference
cc        call repl(gtaucl(1,1,iimp),taucl(1,1,iimp),kmymax,kmymaxp)
cc       ------------------------------------------------------
cc        call tripmt(dmatpa(1,1,iimp),taucl(1,1,iimp),dmata(1,1,iimp),
cc     >              kmymax,kmymax,kmymaxp)
cc       ------------------------------------------------------
c        if(itest.gt.2) then
c          write(6,'(/'' tau -impurity '',i3)') iimp
c          call outmat1(taucl(1,1,iimp),kmymax,kmymax,kmymaxp,tol,6)
c        end if
c      end do
      DEALLOCATE ( deltat,
     >                  STAT = AllocateStatus)
      call alloccheck( AllocateStatus,
     >     'deltat dealloc in ecoreimp                        ' )
      DEALLOCATE ( tauh,
     >                  STAT = AllocateStatus)
      call alloccheck( AllocateStatus,
     >     'tauh dealloc in ecoreimp                          ' )
      DEALLOCATE ( tauclij,
     >                  STAT = AllocateStatus)
      call alloccheck( AllocateStatus,
     >     'tauclij dealloc in ecoreimp                       ' )
c
       return
       end
