c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
c Read geometry of the cluster and cluster-dependent input
c
      subroutine readclu(nintfc,nimp,nposimp,npair,npair0,npair1,kpair,
     >                   kpairind,
c    >                   iscreencl,
     >                   bthcl,bphcl,rbcl,sxcl,taupath,
     >                   imppot,madmax,lfix1,lfix2,
     >                   ib0,b0,impb0,impb0f,isigb)
c
c
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
c -----------immpurities----------
      character*30 imppot
      character*300 taupath
c
      logical uniform
      logical dpair
      logical dpairgl
      logical dummylogic
c
      integer iscreen0
c     integer iscreencl
c
      integer nposimp(3,nimp)
      integer kpairind(nimp,nimp)
      integer madmax
      integer nimp
      integer iimp
      integer limp
c
      integer kpair(4,npair0)
      integer npair,npair0,npair1
      integer kpairtmp(4)
      integer ipair
c
      real*8 rbcl(3,nimp)
      real*8 bthcl(nimp)
      real*8 bphcl(nimp)
      real*8 sxcl(nimp)
      real*8 rimp1
      real*8 rimp2
      real*8 rimp3
      real*8 fact
      real*8 b0
c
      integer ib0
      integer impb0
      integer impb0i
      integer impb0f
      integer isigb(nimp)
c
      real*8 a1
      real*8 a2
      real*8 cvec
      integer nextra
      integer nbulkl
      integer nbulkr
      integer nprc
      integer ninprc
c --- Common blocks ---
c
      common/brav2d/a1(2),a2(2)
      common/test/itest
      common/lay2d/cvec(mtotal,3),nextra,nbulkl,nbulkr,
     &             nprc,ninprc(0:mprc+1)
c
      real*8 pi
      data pi/3.1415926535897932d0/
c --------------------Impurities--------------------
c
      write(6,'(/,/'' ======================================= '')')
      write(6,'(''          CLUSTER CONFIGURATION '')')
      write(6,'('' ======================================= '',/)')
c
      open(unit=5,file='cluster_geo.in',status='old')
c
c      read (5,*) nimp
c    bnyari: we need a dummy row, nimp is in input_rsp.in
c    laszloffy: I do not support it, that is why we have cluster_geo.in.
c               An other solution is implemented in readini, where we read
c               the first line of cluster_geo.in
      read(5,*) 
c
c      if(nimp.gt.mimp) then
c        write(6,*) '<readclu>: STOP: Too many impurities (nimp > mimp)'
c       call flush(6)
c        stop '<readclu>: STOP: too many impurities'
c      end if
c
      do iimp=1,nimp
      read(5,*) nposimp(1,iimp),nposimp(2,iimp),nposimp(3,iimp)
        if(nposimp(3,iimp).le.0.OR.nposimp(3,iimp).gt.nintfc) then
         write (6,*) '<readclu>: STOP: Bad layer position'
c        call flush(6)
         stop '<readclu>: STOP: Bad layer position'
        end if
      end do
c
      if(nimp.gt.1) then
      kpair(1,1)=nposimp(1,1)-nposimp(1,2)
      kpair(2,1)=nposimp(2,1)-nposimp(2,2)
      kpair(3,1)=nposimp(3,1)
      kpair(4,1)=nposimp(3,2)
      end if
c
      do iimp=1,nimp
        do jimp=1,nimp
          kpairind(iimp,jimp)=0
        end do
      end do
c
      if(nimp.gt.1) then
        npair=1
        do iimp=1,nimp-1
          do jimp=iimp+1,nimp
            dpairgl=.false.
            kpairtmp(1)=nposimp(1,iimp)-nposimp(1,jimp)
            kpairtmp(2)=nposimp(2,iimp)-nposimp(2,jimp)
            kpairtmp(3)=nposimp(3,iimp)
            kpairtmp(4)=nposimp(3,jimp)
            do ipair=1,npair
              dpair=(kpairtmp(1).eq.kpair(1,ipair)).AND. 
     >              (kpairtmp(2).eq.kpair(2,ipair)).AND. 
     >              (kpairtmp(3).eq.kpair(3,ipair)).AND. 
     >              (kpairtmp(4).eq.kpair(4,ipair))
              if(dpair) then
                kpairind(iimp,jimp)=ipair
                kpairind(jimp,iimp)=-ipair
                dpairgl=.true.
              end if
              dpair=(kpairtmp(1).eq.-kpair(1,ipair)).AND. 
     >              (kpairtmp(2).eq.-kpair(2,ipair)).AND. 
     >              (kpairtmp(3).eq.kpair(4,ipair)).AND. 
     >              (kpairtmp(4).eq.kpair(3,ipair))
              if(dpair) then
                kpairind(iimp,jimp)=-ipair
                kpairind(jimp,iimp)=ipair
                dpairgl=.true.
              end if
              dpairgl=dpairgl.OR.dpair
            end do
            if(.NOT.dpairgl) then
              npair=npair+1
              kpair(1,npair)=kpairtmp(1)
              kpair(2,npair)=kpairtmp(2)
              kpair(3,npair)=kpairtmp(3)
              kpair(4,npair)=kpairtmp(4)
              kpairind(iimp,jimp)=ipair
              kpairind(jimp,iimp)=-ipair
c             if(npair.gt.mpair) then
c               write (6,*) '<readclu>: STOP: npair > mpair'
c               write (6,*) '<readclu>: STOP: increase mpair in param.h'
c               call flush(6)
c               stop '<readclu>: STOP: npair > mpair'
c             end if
            end if
          end do
        end do
        npair1=npair
      else
        npair=0
        npair1=1
      end if
c
      write(6,'(/'' THE POSITIONS OF THE CLUSTER SITES'')')
      write(6,'(/'' In relative coordinates:'')')
      write(6,*)
      do iimp=1,nimp
         write(6,'(''  Imp.'',i3,''  : '',3i3)')
     >           iimp,nposimp(1,iimp),nposimp(2,iimp),nposimp(3,iimp)
      end do
c
      n0 = ninprc(0)*(nextra+1)
c
      write(6,*)
      write(6,'(/'' In absolute coordinates:'')')
      write(6,*)
      do iimp=1,nimp
c
         limp=nposimp(3,iimp)
c
         rimp1=nposimp(1,iimp)*a1(1)+nposimp(2,iimp)*a2(1)+
     >         cvec(limp+n0,1)
         rimp2=nposimp(1,iimp)*a1(2)+nposimp(2,iimp)*a2(2)+
     >         cvec(limp+n0,2)
         rimp3=cvec(limp+n0,3)
c
         write(6,'(''  Imp.'',i3,'': '',(3f15.8))')
     >           iimp,rimp1,rimp2,rimp3
      end do
c
      write(6,*)
      write(6,'(/'' The different pairs'')')
      write(6,*)
      do ipair=1,npair
         write(6,'(''  ipair.'',i3,''  : '',4i3)')
     >           ipair,kpair(1,ipair),kpair(2,ipair),kpair(3,ipair),
     >           kpair(4,ipair)
      end do
c
c DEBUG BEGINN
c     do iimp=1,nimp
c       do jimp=1,nimp
c         ipair=abs(kpairind(iimp,jimp))
c         write(6,'(''  iimp.'',i3,''  jimp.'',i3,''  : '',11i3)')
c    >      iimp,jimp,kpairind(iimp,jimp),
c    >      nposimp(1,iimp),nposimp(2,iimp),nposimp(3,iimp),
c    >      nposimp(1,jimp),nposimp(2,jimp),nposimp(3,jimp),
c    >      kpair(1,ipair),kpair(2,ipair),kpair(3,ipair),kpair(4,ipair)
c       end do
c     end do
c DEBUG END
c
c read in orientation of magnetic field in the cluster
c
      read(5,*) uniform
      read(5,*) bthcl(1),bphcl(1)
      if(uniform) then
        do iimp=1,nimp
          bthcl(iimp)=bthcl(1)
          bphcl(iimp)=bphcl(1)
        end do
      else
        do iimp=2,nimp
          read(5,*) bthcl(iimp),bphcl(iimp)
        end do
      end if
c
      fact=pi/180.d0
      do iimp=1,nimp
        bthcl(iimp)=fact*bthcl(iimp)
        bphcl(iimp)=fact*bphcl(iimp)
      end do
c
      do iimp=1,nimp
        rbcl(1,iimp)=dsin(bthcl(iimp))*dcos(bphcl(iimp))
        rbcl(2,iimp)=dsin(bthcl(iimp))*dsin(bphcl(iimp))
        rbcl(3,iimp)=dcos(bthcl(iimp))
      end do
c
c  read in scaling factors for SOC
c
      read(5,*) sxcl(1)
      do iimp=1,nimp
        sxcl(iimp)=sxcl(1)
      end do
c
      write(6,'(/'' Orientation of magnetic fields'',
     >           '' in the cluster:'')')
      write(6,*)
      do iimp=1,nimp
        write(6,'(i3,3f15.8)') iimp,(rbcl(i,iimp),i=1,3)
      end do
c
      read(5,*) iscreencl
      if(iscreencl.le.0) then
        write(6,'(/'' No screening in cluster calculation !'')')
      else
        write(6,*) ' Cluster calculation in screened reprezentation'
        write(6,*) ' Not yet implementes'
        write(6,*) ' iscreencl=0'
        iscreencl=0
      end if
c
c  read in taupath
c     read(5,'(a15)') taupath
      read(5,*) dummylogic
      read(5,*) taupath
        write(6,'('' Name of the output tau files: '',a)') taupath
      read(5,*) imppot
        if(itest.ge.2) then
          write(6,'('' IMPPOT= '',a)') imppot
        end if
      read(5,*) madmax
       write(6, '('' Maximal momentum of the
     > Madelung expansion: '',i5)') madmax
      read(5,*) lfix1,lfix2
      if(lfix1.ge.1) then
       if(lfix2.ge.1) then
         write(6,'('' Keep orientation of magnetization fixed '',
     >    ''at site:'',i3,i5)') lfix1,lfix2
        else
          write(6,'('' Keep orientation of magnetization fixed '',
     >    ''at site:'',i3)') lfix1
        end if
      end if
      read(5,*) ib0,b0
      read(5,*) impb0,impb0f,isigb(1)
      if(impb0.gt.impb0f) then
        impb0i = impb0f
        impb0f=impb0
        impb0=impb0i
      endif
      do iimp=1,nimp
        isigb(iimp)=1
      end do
      if(isigb(1).ne.1) then
        do iimp=impb0,impb0f
          isigb(iimp)=-1
        end do
      end if
c
c
c DEBUG BEGINN
      if (itest.ge.2) then
        write(6,'(/''<readclu>: readclu processed!'')')
        call flush(6)
      end if
c DEBUG END
      return
      end
