c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine tcpa_main(
     > leftmat,rightmat,laymat,bulk,dos,kset,nintfc,ninprcl,ninprcr,
     > conc,concl,concr,cpain,cpamatin,cpamatinl,cpamatinr)
c=======================
c
c -decide about opening binary files containing (inverse of) 
c  effective t-matrices 
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      logical dos,bulk
      logical cpamatin,cpamatinl,cpamatinr,cpain
      logical cpai,cpal,cpar,cpa
c
      character*30 leftmat,rightmat,laymat
c
      dimension conc(mintfc),concl(minprc),concr(minprc)
c
      data tiny/1.0d-6/
c
      write(6,'(/''  routine TCPA_MAIN> '')')
      cpamatin=.false.
      cpamatinl=.false.
      cpamatinr=.false.
c
c confirm whether there is any disorder
c
      cpal=.false.
      do li=1,ninprcl
        cpal=cpal.or.(1.d0-concl(li).gt.tiny)
      end do
      cpar=.false.
      do li=1,ninprcr
        cpar=cpar.or.(1.d0-concr(li).gt.tiny)
      end do
      cpai=.false.
      do li=1,nintfc
        cpai=cpai.or.(1.d0-conc(li).gt.tiny)
      end do
      cpa=cpal.or.cpar.or.cpai
      write(6,'('' cpal '',l1,5x,''cpar '',l1,5x,''cpai '',l1)') 
     >cpal,cpar,cpai
      if(.not.cpa) goto 1000
c
      if(.not.dos) then
c
c SCF calculation
c
        if(.not.bulk) then
c
c    SURFACE (INTERFACE)
c
          write(6,'(/'' Surface or interface SCF'')')
          if(cpal) then
            cpamatinl=.true.
      write(6,'('' Effective t-matrices read in from '',a)') leftmat
          end if
          if(cpar) then
            cpamatinr=.true.
      write(6,'('' Effective t-matrices read in from '',a)') rightmat
          end if
        end if
c
      else
c
c DOS calculation
c
        if(kset.eq.0) then
c
c    SPECTRAL DOS
c        
          if(bulk) then
c       
c        BULK
c
            write(6,'(/'' Bulk spectral DOS'')')
            if(cpa) then
              cpamatin=.true.
              write(6,'('' Effective t-matrices read in from '',a)') 
     >        laymat
              write(6,*) cpamatin
            end if
c
          else
c
c        SURFACE (INTERFACE) 
c
            write(6,'(/'' Surface or interface spectral DOS'')')
            if(cpal) then
              cpamatinl=.true.
              write(6,'('' Effective t-matrices read in from '',a)') 
     >        leftmat
            end if
            if(cpar) then
              cpamatinr=.true.
              write(6,'('' Effective t-matrices read in from '',a)') 
     >        rightmat
            end if
            if(cpai) then
              cpamatin=.true.
              write(6,'('' Effective t-matrices read in from '',a)') 
     >        laymat
            end if
c
          end if
c
        else
c
c    TOTAL DOS
c        
          if(.not.bulk) then
c
c        SURFACE (INTERFACE) 
c
            write(6,'(/'' Surface or interface total DOS'')')
            if(cpal) then
              cpamatinl=.true.
              write(6,'('' Effective t-matrices read in from '',a)') 
     >        leftmat
            end if
            if(cpar) then
              cpamatinr=.true.
              write(6,'('' Effective t-matrices read in from '',a)') 
     >        rightmat
            end if
c
          end if
c
          if(cpai.and.cpain) then
              cpamatin=.true.
              write(6,'('' Effective t-matrices read in from '',a)') 
     >        laymat
          end if
c
        end if
c
      end if
c
 1000 continue
      write(6,'('' cpamatinl '',l1)') cpamatinl
      write(6,'('' cpamatinr '',l1)') cpamatinr
      write(6,'('' cpamatin  '',l1)') cpamatin 
      write(6,*)
      call flush(6)
c
      return
      end
      subroutine tcpa_open(
     > for010,leftmat,rightmat,laymat,cpamatin,cpamatinl,cpamatinr)
c=========================
c
c -open binary files containing (inverse of) effective t-matrices 
c
      implicit real*8 (a-h,o-z)
c
      parameter(nle=15,nri=16,nlay=17)
c
      character*30 for010,leftmat,rightmat,laymat
      logical cpamatin,cpamatinl,cpamatinr
c
      open(10,file=for010,form='unformatted',status='unknown')
      if(cpamatin) 
     >open(nlay,file=laymat,form='unformatted',status='old')
      if(cpamatinl) 
     >open(nle,file=leftmat,form='unformatted',status='old')
      if(cpamatinr) 
     >open(nri,file=rightmat,form='unformatted',status='old')
c
      return
      end
      subroutine tcpa_close(
     > leftmat,rightmat,laymat,cpamatin,cpamatinl,cpamatinr)
c==========================
c
c -close binary files containing (inverse of) effective t-matrices 
c
      implicit real*8 (a-h,o-z)
c
      parameter(nle=15,nri=16,nlay=17)
c
      character*30 leftmat,rightmat,laymat
      logical cpamatin,cpamatinl,cpamatinr
c
      close(9)
      if(cpamatin) close(nlay)
      if(cpamatinl) close(nle)
      if(cpamatinr) close(nri)
c
      return
      end
      subroutine tcpa_in(leftmat,rightmat,laymat,ce,nintfc,
     > ninprcl,ninprcr,conc,concl,concr,cpamatin,cpamatinl,
     > cpamatinr,n,tminv,tminvl,tminvr,
     > tau,taua,taub,ndim)
c=======================
c
c -read in (inverse of) effective t-matrices if necessary
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
      parameter(nle=15,nri=16,nlay=17)
c
      logical cpamatin,cpamatinl,cpamatinr,cpalay
c
      character*30 leftmat,rightmat,laymat
c
      dimension conc(mintfc),concl(minprc),concr(minprc)
c
      complex*16 ce,cein
      complex*16 tminvl(ndim,ndim,minprc)
      complex*16 tminvr(ndim,ndim,minprc)
      complex*16 tminv(ndim,ndim,mintfc)
      complex*16 tau(ndim,ndim,mintfc)
      complex*16 taua(ndim,ndim,mintfc)
      complex*16 taub(ndim,ndim,mintfc)
c
      data tiny/1.0d-6/
c
      if(cpamatinl) then
        write(6,*) 'leftmats read from ', leftmat
        read(nle) cein
        if(cdabs(ce-cein).gt.1.0d-12) 
     >  stop ' check energy array leftmat'
        do li=1,ninprcl
          cpalay=1.d0-concl(li).gt.tiny
          if(cpalay) then
            read(nle) ((tminvl(k1,k2,li),k1=1,n),k2=1,n)
          end if
        enddo
      end if
      if(cpamatinr) then
        write(6,*) 'rightmats read from ', rightmat
        read(nri) cein
        if(cdabs(ce-cein).gt.1.0d-08) 
     >  stop ' check energy array rightmat'
        do li=1,ninprcr
          cpalay=1.d0-concr(li).gt.tiny
          if(cpalay) 
     >    read(nri) ((tminvr(k1,k2,li),k1=1,n),k2=1,n)
        enddo
      end if
      if(cpamatin) then
        write(6,*) 'laymats read from ', laymat
        read(nlay) cein
        if(cdabs(ce-cein).gt.1.0d-12)
     >  stop ' check energy array laymat'
        do li=1,nintfc
          cpalay=1.d0-conc(li).gt.tiny
          if(cpalay) then
            read(nlay) ((tminv(k1,k2,li),k1=1,n),k2=1,n)
            read(nlay) ((tau(k1,k2,li),k1=1,n),k2=1,n)
            read(nlay) ((taua(k1,k2,li),k1=1,n),k2=1,n)
            read(nlay) ((taub(k1,k2,li),k1=1,n),k2=1,n)
c          tminvl(1:n,1:n,li) = tminv(1:n,1:n,li)
c          tminvr(1:n,1:n,li) = tminv(1:n,1:n,li)
          end if
        end do
      end if 
c
      return
      end
      subroutine tcpa_out(ce,nintfc,conc,n,tminv,tau,taua,taub,ndim)
c=======================
c
c -write out (inverse of) effective t-matrices
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      logical cpalay
c
      dimension conc(mintfc)
c
      complex*16 ce
      complex*16 tminv(ndim,ndim,mintfc)
      complex*16 tau(ndim,ndim,mintfc)
      complex*16 taua(ndim,ndim,mintfc)
      complex*16 taub(ndim,ndim,mintfc)
c
      data tiny/1.0d-6/
c
      write(10) ce
      do li=1,nintfc
        cpalay=1.d0-conc(li).gt.tiny
        if(cpalay) then 
           write(10) ((tminv(k1,k2,li),k1=1,n),k2=1,n)
           write(10) ((tau(k1,k2,li),k1=1,n),k2=1,n)
           write(10) ((taua(k1,k2,li),k1=1,n),k2=1,n)
           write(10) ((taub(k1,k2,li),k1=1,n),k2=1,n)
        end if
      end do
c
      return
      end
