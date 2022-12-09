c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine restart_in1(itscf0,emax,v0,opot,
     >                      vra,bra,bopra,rba,
     >                      nintfc,idpota,za,laycore,for006)
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      integer mmax
      parameter (mmax=mimp)
c
      logical rest,opot
      character*30 for012,for006,laycore
      character*10 idpota(mmax)
      character*3  idel
      character*4  titrea(iorb,mmax)
!     character*4  titreb(iorb,mmax)
      character*4  oplus(4),ominus(3)
c
      dimension rba(3,mmax)
      dimension vra(nrad,mmax)
      dimension bra(nrad,mmax)
      dimension bopra(nrad,2,mmax)
      dimension norba(mmax)
      dimension dena(iorb,mmax)
      dimension nela(iorb,mmax)
      dimension nqna(iorb,mmax)
      dimension nka(iorb,mmax)
      dimension nqla(iorb,mmax)
      dimension dfla(iorb,mmax)
      dimension dq1a(iorb,mmax)
      dimension za(mmax)
c
      common/core1a/dena,dfla,dq1a,norba,nela,nqna,nka,nqla
!     common/core1b/denb,dflb,dq1b,norbb,nelb,nqnb,nkb,nqlb
      common/core2/test,nuc,nes,iskip
      common/core3/titrea!,titreb
c
      data oplus/'s1/2','p3/2','d5/2','f7/2'/
      data ominus/'p1/2','d3/2','f5/2'/
      data tiny/1.0d-6/
c
      dvc=137.036d0
      ij=index(for006,'.')-1
      for012=for006(1:ij)//'.update_imp'
c
      inquire(file=for012,EXIST=rest)
      if(rest) then
         open(12,file=for012,status='OLD',form='UNFORMATTED')
         read(12) relnum
         read(12) itscf0
         write(6,'(/'' RESTARTING from iteration '',i3, ''!!!''/)') 
     >   itscf0
!        if(bulk) then
!           emax=relnum
!        else
            v0=relnum
!        end if
         do li=1,nintfc
!          cpalay=1.d0-conc(li).gt.tiny
           read(12) (vra(i,li),i=1,nrad)
           read(12) (bra(i,li),i=1,nrad)
           if(opot) read(12) ((bopra(i,is,li),i=1,nrad),is=1,2)
           read(12) (rba(i,li),i=1,3)
!          if(cpalay) then
!            read(12) (vrb(i,li),i=1,nrad)
!            read(12) (brb(i,li),i=1,nrad)
!            if(opot) read(12) ((boprb(i,is,li),i=1,nrad),is=1,2)
!            read(12) (rbb(i,li),i=1,3)
!          else
!            do i=1,nrad
!              vrb(i,li)=vra(i,li)
!              brb(i,li)=bra(i,li)
!              boprb(i,1,li)=bopra(i,1,li)
!              boprb(i,2,li)=bopra(i,2,li)
!            end do
!            do i=1,3
!              rbb(i,li)=rba(i,li)
!            end do
!          end if
         end do
c read core info from update file
         read(12) nuc,nes,test,iskip
         do li=1,nintfc
!          cpalay=1.d0-conc(li).gt.tiny
           read(12) norba(li)
           if(norba(li).gt.iorb) stop 'Increase parameter iorb'
           read(12) (dena(i,li),nqna(i,li),nka(i,li),nela(i,li),
     >               i=1,norba(li))
!          if(cpalay) then
!            read(12) norbb(li)
!            if(norbb(li).gt.iorb) stop 'Increase parameter iorb'
!            read(12) (denb(i,li),nqnb(i,li),nkb(i,li),nelb(i,li),
!    >                 i=1,norbb(li))
!          else
!            norbb(li)=norba(li)
!            do i=1,norbb(li)
!              denb(i,li)=dena(i,li)
!              nqnb(i,li)=nqna(i,li)
!              nkb(i,li)=nka(i,li)
!              nelb(i,li)=nela(i,li)
!            end do
!          end if
         end do
         close(12)
      else
        itscf0=0
c read core info from general input deck
        open(20,file=laycore,status='old')
        read(20,*) nuc,nes,test,iskip
        do 1 li=1,nintfc
c
           if(idpota(li)(1:3).eq.'Vac') goto 2
c
           rewind 20
           read(20,*)
   10      read(20,'(a3)') idel
           read(20,*) norbit
           if(idel.ne.idpota(li)(1:3)) then
             do k=1,norbit
               read(20,*)
             end do
             goto 10
           end if
           if(iskip.ge.1) write(6,'(/a10,a3)') idpota(li),idel
           norba(li)=norbit
           if(norba(li).gt.iorb) stop 'Increase parameter iorb'
c
c     read orbital information
c
           do i=1,norba(li)
             read(20,*) dena(i,li),nqna(i,li),nka(i,li),nela(i,li)
           end do
c
    2      continue
!          if(idpotb(li)(1:3).eq.'Vac') goto 1
c
!          rewind 20
!          read(20,*)
!  15      read(20,'(a3)') idel
!          read(20,*) norbit
!          if(idel.ne.idpotb(li)(1:3)) then
!            do k=1,norbit
!              read(20,*)
!            end do
!            goto 15
!          end if
!          if(iskip.ge.1) write(6,'(/a10,a3)') idpotb(li),idel
!          norbb(li)=norbit
!          if(norbb(li).gt.iorb) stop 'Increase parameter iorb'
c
c     read orbital information
c
!         do i=1,norbb(li)
!           read(20,*) denb(i,li),nqnb(i,li),nkb(i,li),nelb(i,li)
!         end do
c
    1   continue 
      close(20)
      end if
c
c initialize calculation of core energies
c
      do 101 li=1,nintfc
        if(idpota(li)(1:3).eq.'Vac') goto 102
c
        do i=1,norba(li)
c
          dval=za(li)/dvc
          nqla(i,li)=iabs(nka(i,li))
          if (nka(i,li).lt.0) nqla(i,li)=nqla(i,li)-1
          dfla(i,li)=nka(i,li)*nka(i,li)
          dfla(i,li)=dsqrt(dfla(i,li)-dval*dval)
          l=1
          j=nqna(i,li)-nqla(i,li)
          if((j-2*(j/2)).eq.0) l=-l
          dq1a(i,li)=l*nka(i,li)/iabs(nka(i,li))
          if(nka(i,li).lt.0) go to 25
c
c       j  = l - 1/2
c
          do 20 ll=1,3
            if(ll.ne.nka(i,li)) go to 20
            titrea(i,li)=ominus(ll)
            go to 40
   20     continue
c
c       j = l + 1/2
c
   25     do 30 ll=1,4
            lp1=-ll
            if(lp1.ne.nka(i,li)) go to 30
            titrea(i,li)=oplus(ll)
            go to 40
   30     continue
c
   40     continue
          if(iskip.ge.2) write(6,'(2i6,5x,a4,i5,d20.10)')
     >    nqna(i,li),nka(i,li),titrea(i,li),nela(i,li),dena(i,li)
c
        end do
c
        if(iskip.ge.2) write(6,*)
c
  102   continue
!       if(idpotb(li)(1:3).eq.'Vac') goto 101
!       do i=1,norbb(li)
c
!         dval=zb(li)/dvc
!         nqlb(i,li)=iabs(nkb(i,li))
!         if (nkb(i,li).lt.0) nqlb(i,li)=nqlb(i,li)-1
!         dflb(i,li)=nkb(i,li)*nkb(i,li)
!         dflb(i,li)=dsqrt(dflb(i,li)-dval*dval)
!         l=1
!         j=nqnb(i,li)-nqlb(i,li)
!         if((j-2*(j/2)).eq.0) l=-l
!         dq1b(i,li)=l*nkb(i,li)/iabs(nkb(i,li))
!         if(nkb(i,li).lt.0) go to 125
!
c       j  = l - 1/2
c
!         do 120 ll=1,3
!           if(ll.ne.nkb(i,li)) go to 120
!           titreb(i,li)=ominus(ll)
!           go to 140
! 120     continue
c
c       j = l + 1/2
c
! 125     do 130 ll=1,4
!           lp1=-ll
!           if(lp1.ne.nkb(i,li)) go to 130
!           titreb(i,li)=oplus(ll)
!           go to 140
! 130     continue
c
! 140     continue
!         if(iskip.ge.2) write(6,'(2i6,5x,a4,i5,d20.10)')
!    >    nqnb(i,li),nkb(i,li),titreb(i,li),nelb(i,li),denb(i,li)
!
!       end do
c
        if(iskip.ge.2) write(6,*)
c
 101  continue
      return
      end
c
c==========================
      subroutine restart_out1(itscf,emax,v0,opot,
     >                       vra,bra,bopra,rba,
     >                       nintfc,for006)
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      integer mmax
      parameter (mmax=mimp)
c
      logical opot
      character*30 for012,for006
c
      dimension rba(3,mmax)
      dimension vra(nrad,mmax)
      dimension bra(nrad,mmax)
      dimension bopra(nrad,2,mmax)
      dimension norba(mmax)
      dimension dena(iorb,mmax)
      dimension nela(iorb,mmax)
      dimension nqna(iorb,mmax)
      dimension nka(iorb,mmax)
      dimension nqla(iorb,mmax)
      dimension dfla(iorb,mmax)
      dimension dq1a(iorb,mmax)
c
      common/core1a/dena,dfla,dq1a,norba,nela,nqna,nka,nqla
!     common/core1b/denb,dflb,dq1b,norbb,nelb,nqnb,nkb,nqlb
      common/core2/test,nuc,nes,iskip
c
      data tiny/1.0d-6/
c
      ij=index(for006,' ')-5
      for012=for006(1:ij)//'.update'
c
      open(12,file=for012,status='UNKNOWN',form='UNFORMATTED')
      rewind(12)
c
!     if(bulk) then
!        relnum=emax
!     else
         relnum=v0
!     end if
      write(12) relnum
      write(12) itscf
c
      do li=1,nintfc
!       cpalay=1.d0-conc(li).gt.tiny
        write(12) (vra(i,li),i=1,nrad)
        write(12) (bra(i,li),i=1,nrad)
        if(opot) write(12) ((bopra(i,is,li),i=1,nrad),is=1,2)
        write(12) (rba(i,li),i=1,3)
!       if(cpalay) then
!          write(12) (vrb(i,li),i=1,nrad)
!          write(12) (brb(i,li),i=1,nrad)
!          if(opot) write(12) ((boprb(i,is,li),i=1,nrad),is=1,2)
!          write(12) (rbb(i,li),i=1,3)
!       end if
      end do
c
      write(12) nuc,nes,test,iskip
      do li=1,nintfc
!       cpalay=1.d0-conc(li).gt.tiny
        write(12) norba(li)
        write(12) (dena(i,li),nqna(i,li),nka(i,li),nela(i,li),
     >             i=1,norba(li))
!       if(cpalay) then
!         write(12) norbb(li)
!         write(12) (denb(i,li),nqnb(i,li),nkb(i,li),nelb(i,li),
!    >               i=1,norbb(li))
!       end if
      end do
c
      close(12)
      return
      end
