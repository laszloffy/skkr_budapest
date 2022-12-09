c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine printres1(
     >           for007,for008,for009,it,lmax,nintfc,ne,cear,efermi,
     >           dosa,qvpa,qva,dosmaga,
     >           spin_magvpa,spin_magva,
     >           orb_magvpa,orb_magva,
     >           arot,qmoma,za,qca,
     >           qvdiffa,enbdiffa,
     >           entota,entot,entotifc,
     >           idpota,vra,bra,bopra,
     >           rs,dx,ns,lms,orbpol,
     >           enpba)
c ========================
c
c Print out results 
C
C All the momentum- and spin-resolved results concern the local frame 
C of reference, while the unresolved vector quantities refer to the
C global frame of reference rotated by 'arot'
C
c
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
      integer mmax
      parameter (mmax=mimp)
c
      logical lms
      logical cpasum
      logical orbpol
c
      character*30 for007
      character*30 for008
      character*30 for009
      character*10 idpota(mmax)
      character*1 ab
c
      real*8 vra(nrad,mmax)
      real*8 bra(nrad,mmax)
      real*8 bopra(nrad,2,mmax)
      real*8 rs(mmax)
      real*8 dx(mmax)
      real*8 dosa(kmymaxp,mmax,me)
      real*8 dosmaga(kmymaxp,mmax,me)
      real*8 qvdiffa(mmax,me)
      real*8 enbdiffa(mmax,me)
      real*8 qva(mmax)
      real*8 qca(mmax)
      real*8 qvpa(kmymaxp,mmax)
      real*8 za(mmax)
      real*8 spin_magvpa(kmymaxp,mmax,3)
      real*8 spin_magva(mmax,3)
      real*8 orb_magvpa(kmymaxp,mmax,3)
      real*8 orb_magva(mmax,3)
      real*8 entota(mmax)
      real*8 entot(mmax)
      real*8 conc(mmax)
c
      integer ns(mmax)
      integer ii
c
      complex*16 cear(me)
      complex*16 qmoma(lmsup,mmax)
c
      real*8 enpba(kmymaxp,mmax)
c
      common/test/itest
      common/madelung_energy/emad
c
      data tiny/1.0d-6/
      data tol/1.0d-10/
c
      do ii=1,mmax
       conc(ii)=1.0d0
      end do
      kmax=2*lmax+1
      lmmax=(lmax+1)*(lmax+1)
      kmymax=2*lmmax
      lmaxs=2*lmax
      lmmaxs=(lmaxs+1)*(lmaxs+1)
c
c - Print potential in file for007
c
      open(7,file=for007,status='unknown')
c     rewind(7)
      write(7,'(''SCI: '',i4,2x)') it
c
      zeromt=0.d0
      do li=1,nintfc
         write(7,'(a)') idpota(li)
         write(7,'(f8.5)') 1.00000
         write(7,'(3f8.3,i5,f18.12)') 
     >                  za(li),zeromt,dx(li),ns(li),rs(li)
         write(7,'(4d20.12)') (vra(irad,li),irad=1,ns(li))
         write(7,'(4d20.12)') (bra(irad,li),irad=1,ns(li))
         if(orbpol) then
           write(7,'(4d20.12)') (bopra(irad,1,li),irad=1,ns(li))
           write(7,'(4d20.12)') (bopra(irad,2,li),irad=1,ns(li))
         end if
c        if(cpalay) then
c          write(7,'(a)') idpotb(li)
c          write(7,'(f8.5)') 1.0-conc(li)
c          if(bulk) write(7,'(4d20.12)') qmomb(1,li)+qcb(li)-zb(li),
c    >     (qmomb(i,li),i=2,lmmaxs)
c          write(7,'(3f8.3,i5,f18.12)') 
c    >       zb(li),zeromt,dx(li),ns(li),rs(li)
c          write(7,'(4d20.12)') (vrb(irad,li),irad=1,ns(li))
c          write(7,'(4d20.12)') (brb(irad,li),irad=1,ns(li))
c          if(orbpol) then
c             write(7,'(4d20.12)') (boprb(irad,1,li),irad=1,ns(li))
c             write(7,'(4d20.12)') (boprb(irad,2,li),irad=1,ns(li))
c          end if
c        end if
      end do
      close(7)
c
c - Print DOS
c
      open(8,file=for008,status='unknown')
      open(9,file=for009,status='unknown')
c     rewind(8)
      write(8,'(''SCI: '',i4,2x)') it
      write(9,'(''SCI: '',i4,2x)') it
c
      do li=1,nintfc
        write(8,'(''##'',i3,t10,''DOS   A'')') li
        if (lms)
     >  write(9,'(''##'',i3,t10,''PDOS   A'')') li
c       -------------------------------------------
        call dosprint(lms,ich,lmax,li,ne,cear,dosa)
c       -------------------------------------------
c       if(cpalay) then
c         write(8,'(''##'',i3,t10,''DOS   B'')') li
c         -------------------------------------------
c         call dosprint(lms,ich,lmax,li,ne,cear,dosb)
c         -------------------------------------------
c       end if
      end do
c
c -Print MDOS
c
      write(8,*)
      do li=1,nintfc
c       cpalay=1.d0-conc(li).gt.tiny
        write(8,'(''##'',i3,t10,''MDOS   A'')') li
        if (lms)
     >  write(9,'(''##'',i3,t10,''MPDOS   A'')') li
c       ----------------------------------------------
        call dosprint(lms,ich,lmax,li,ne,cear,dosmaga)
c       ----------------------------------------------
c       if(cpalay) then
c         write(8,'(''##'',i3,t10,''MDOS   B'')') li
c         ----------------------------------------------
c         call dosprint(lms,ich,lmax,li,ne,cear,dosmagb)
c         ----------------------------------------------
c       end if
      end do
c
c - Print No. of electrons in the layers
c
      cpasum=.false.
c     do li=1,nintfc
c        cpalay=1.d0-conc(li).gt.tiny
c        cpasum=cpasum.or.cpalay
c     end do
c
c     if(cpasum) then
c       write(8,'(/2x,''No. of electrons  - A type - :''/)') 
c     else
        write(8,'(/2x,''No. of electrons :''/)')
c     end if
c
      if(.not.lms) then
         do li=1,nintfc
c          --------------------------------------
           call nelprint1(li,kmax,qvpa(1,li))
c          --------------------------------------
         end do
c        -------------------- do B-type if there was cpalay ------
c        if(cpasum) then
c          write(8,'(/2x,''No. of electrons  - B type - :''/)')
c          do li=1,nintfc
c             cpalay=1.d0-conc(li).gt.tiny
c             if(cpalay) then
c                --------------------------------------
c                call nelprint1(li,kmax,qvpb(1,li))
c                --------------------------------------
c             end if
c          end do
c        end if
c        write(8,150)
         do li=1,nintfc
c           cpalay=1.d0-conc(li).gt.tiny
c           if(cpalay) then
c              ab='A'
c           else
               ab=' '
c           end if
c          ------------------------------------------
            call nelprint2(li,lmax,qvpa,qva,lms,ab,1)
c          ------------------------------------------
!           if(cpalay) then
!              ab='B'
!             -----------------------------------------
!             call nelprint2(li,lmax,qvpb,qvb,lms,ab,1)
!             -----------------------------------------
!           end if
         end do
      else
         write(8,160)
         do is=1,2
            write(8,'(/'' Spin:'',i2)') is
            do li=1,nintfc
c              cpalay=1.d0-conc(li).gt.tiny
c              if(cpalay) then
c                 ab='A'
c              else
                  ab=' '
c              end if
c              ------------------------------------------
               call nelprint2(li,lmax,qvpa,qva,lms,ab,is)
c              ------------------------------------------
c              if(cpalay) then
c                 ab='B'
c                 ------------------------------------------
c                 call nelprint2(li,lmax,qvpb,qvb,lms,ab,is)
c                 ------------------------------------------
c              end if
            end do
         end do
      end if
c
c - Print magnetization in the layers
c
      if(dabs(arot-0.d0).gt.tol) write(8,'(/''WARNING!!'',/,
     > ''Magnetic moments reference frame rotated by'',f12.8)') arot
c
      write(8,'(/29x,''Spin magnetic moments :'',
     >          /29x,''----------------------''/ 
     >          /29x,''         X            '')') 
c     ------------------------------------------------
      call printmag(nintfc,lmax,conc,
     >             spin_magva(1,1),spin_magvpa(1,1,1),
     >             spin_magva(1,1),spin_magvpa(1,1,1),.true.,cpasum)
c     ------------------------------------------------
      write(8,'(/29x,''         Y            '')') 
c     ------------------------------------------------
      call printmag(nintfc,lmax,conc,
     >             spin_magva(1,2),spin_magvpa(1,1,2),
     >             spin_magva(1,2),spin_magvpa(1,1,2),.true.,cpasum)
c     ------------------------------------------------
      write(8,'(/29x,''         Z            '')') 
c     ------------------------------------------------
      call printmag(nintfc,lmax,conc,
     >             spin_magva(1,3),spin_magvpa(1,1,3),
     >             spin_magva(1,3),spin_magvpa(1,1,3),.true.,cpasum)
c     ------------------------------------------------
c
      write(8,'(/29x,''Orbital magnetic moments :'',
     >          /29x,''----------------------''/ 
     >          /29x,''         X            '')') 
c     ------------------------------------------------
      call printmag(nintfc,lmax,conc,
     >             orb_magva(1,1),orb_magvpa(1,1,1),
     >             orb_magva(1,1),orb_magvpa(1,1,1),.true.,cpasum)
c     ------------------------------------------------
      write(8,'(/29x,''         Y            '')') 
c     ------------------------------------------------
      call printmag(nintfc,lmax,conc,
     >             orb_magva(1,2),orb_magvpa(1,1,2),
     >             orb_magva(1,2),orb_magvpa(1,1,2),.true.,cpasum)
c     ------------------------------------------------
      write(8,'(/29x,''         Z            '')') 
c     ------------------------------------------------
      call printmag(nintfc,lmax,conc,
     >             orb_magva(1,3),orb_magvpa(1,1,3),
     >             orb_magva(1,3),orb_magvpa(1,1,3),.true.,cpasum)
c     ------------------------------------------------
c
c - Print moments of charge densities (in global frame of reference)
c
      write(8,'(//''Moments of charge densities:''
     >           /''---------------------------''/)')
      i=0
      do l=0,2*lmax 
      do m=-l,l     
         write(8,'(/10x,''(l,m)'',2i3/)') l,m
         i=i+1
!        if(cpasum) then
!           write(8,'(t26,''A'',t66,''B'')')
!        end if
         do li=1,nintfc
!          cpalay=1.d0-conc(li).gt.tiny
!          if(cpalay) then
!            write(8,'(2x,'' L'',i2,t10,2d15.7,5x,2d15.7)') li,
!    >       qmoma(i,li),qmomb(i,li)
!          else
             write(8,'(2x,'' L'',i4,t10,2d15.7)') li,qmoma(i,li)
!          end if
         end do
      end do
      end do
c
c - Print total energies
c
      write(8,'(//29x,''Total energy:''/,
     >           29x,''--------------'',//,''Layer'')')
      cpasum=.false.
      etka=0.0
      etkb=0.0
      do li=1,nintfc
!       cpalay=1.d0-conc(li).gt.tiny
        ab=' '
!       if(cpalay) ab='A'
        etka=etka+entota(li)
        write(8,'(1x,i4,1x,a1,1x,f15.7)') li,ab,entota(li)
!       if(cpalay) then
!          ab='B'
!          cpasum=cpasum.or.cpalay
!          etkb=etkb+entotb(li)
!          write(8,'(4x,a1,1x,f15.7)') ab,entotb(li)
!       end if
      end do
!     if(cpasum) then
!     write(8,'(/'' A-type Ifc '',f15.7)') etka
!     write(8,'(/'' B-type Ifc '',f15.7)') etkb
!     write(8,'(/'' Madelung   '',f15.7)') emad
!     write(8,'(/'' Sum    Ifc '',f15.7)') entotifc+emad
!     else
      write(8,'(/'' Ifc  '',e16.7)') entotifc+emad
!     end if
c
c - Print energy-resolved contributions to band energy  
c
      write(8,'(//''Band-energy:''
     >           /''------------''/)')
      omint=0.d0
      do ie=ne,1,-1 
         om=0.d0
         do li=1,nintfc
           om=om+(enbdiffa(li,ie)-efermi*qvdiffa(li,ie))
!          om=om+conc(li)*(enbdiffa(li,ie)-efermi*qvdiffa(li,ie))
!    >   +(1.d0-conc(li))*(enbdiffb(li,ie)-efermi*qvdiffb(li,ie))
         end do
         omint=omint+om
         write(8,'(i4,2x,2d15.7,5x,2d17.7)') ie,cear(ie),om,omint
      end do
c
      write(8,'(//''Partial band-energy per site:''
     >           /''-----------------------------''/)')
      do li=1,nintfc
         do k=1,kmymax
            enpba(k,li)=enpba(k,li)-qvpa(k,li)*efermi
         end do
         write(8,'(i4,32e16.8)') li,(enpba(k,li),k=1,kmymax)
      end do
c
      close(8)
      if (lms) close(9)
c
      return
 150  format(/t15,
     >          '      s    ',
     >          '    p 1/2  ',
     >          '    p 3/2  ',
     >          '    d 3/2  ',
     >          '    d 5/2  ',
     >          '    total  ')
 160  format(/t20,
     >          '         s     ',
     >          '         p     ',
     >          '         d     ',
     >          '       total   ')
      end
c     ================================================
      subroutine dosprint(lms,ich,lmax,li,ne,cear,dos)
c
c     ================================================
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
      integer mmax
      parameter (mmax=mimp)
c
      logical lms
c
      dimension dos(kmymaxp,mmax,me)
      dimension ddd(kmaxp),ppp(0:lmaxp,2)
c
      complex*16 cear(me)
c
      kmax=2*lmax+1
      lmmax=(lmax+1)*(lmax+1)
      kmymax=2*lmmax
c
      do ie=1,ne
        er=dreal(cear(ie))
        if(.not.lms) then
          sum=0.d0
          kmy=0
          do k=1,kmax
             l=k/2
             if(2*l.eq.k) then
               kap=l
             else
               kap=-l-1
             end if
             nk=2*iabs(kap)
             ddd(k)=0.d0
             do j=1,nk
                kmy=kmy+1
                sum=sum+dos(kmy,li,ie)
                ddd(k)=ddd(k)+dos(kmy,li,ie)
             end do
          end do
          write(8,'(f10.6,6e15.8)') er,(ddd(k),k=1,kmax),sum
        else
          ii=0
          do is=1,2
             sum=0.d0
             do l=0,lmax
                ppp(l,is)=0.d0
                do m=-l,l
                   ii=ii+1
                   sum=sum+dos(ii,li,ie)
                   ppp(l,is)=ppp(l,is)+dos(ii,li,ie)
                end do
             end do
             write(8,'(f10.6,6e15.8)') er,(ppp(l,is),l=0,lmax),sum
          end do
          write(9,'(f10.6,50e15.8)') er,(dos(l,li,ie),l=1,2*lmmax)
        end if
      end do
c
      return
      end
c
c     =========================================================
      subroutine printmag(nintfc,lmax,
     >             conc,qmaga,qmagpa,qmagb,qmagpb,lms,cpasum)
c
c     =========================================================
c
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
      integer mmax
      parameter (mmax=mimp)
c
      logical lms,cpasum,cpalay
      character*1 ab
      dimension qmagpa(kmymaxp,mmax),qmaga(mmax)
      dimension qmagpb(kmymaxp,mmax),qmagb(mmax)
      dimension conc(mmax)
c
      data tiny/1.0d-6/
c
      kmax=2*lmax+1
      lmmax=(lmax+1)*(lmax+1)
      kmymax=2*lmmax
c
      if(.not.lms) then
        if(cpasum) then
          write(8,'(2x,'' -- A type --'')')
        end if
        do li=1,nintfc
c          ------------------------------------
           call nelprint1(li,kmax,qmagpa(1,li))
c          ------------------------------------
        end do
        if(cpasum) then
          write(8,'(2x,'' -- B type --'')')
          do li=1,nintfc
             cpalay=1.d0-conc(li).gt.tiny
             if(cpalay) then
c              --------------------------------
               call nelprint1(li,kmax,qmagpb(1,li))
c              --------------------------------
             end if
          end do
        end if
        write(8,150)
        do li=1,nintfc
           cpalay=1.d0-conc(li).gt.tiny
           if(cpalay) then
              ab='A'
           else
              ab=' '
           end if
c          ---------------------------------------------
           call nelprint2(li,lmax,qmagpa,qmaga,lms,ab,1)
c          ---------------------------------------------
           if(cpalay) then
              ab='B'
c             ---------------------------------------------
              call nelprint2(li,lmax,qmagpb,qmagb,lms,ab,1)
c             ---------------------------------------------
           end if
        end do
      else
        write(8,160)
        do is=1,2
           write(8,'('' Spin:'',i2)') is
           do li=1,nintfc
              cpalay=1.d0-conc(li).gt.tiny
               if(cpalay) then
                  ab='A'
               else
                  ab=' '
               end if
c              ---------------------------------------------
               call nelprint2(li,lmax,qmagpa,qmaga,lms,ab,is)
c              ---------------------------------------------
               if(cpalay) then
                  ab='B'
c                 ----------------------------------------------
                  call nelprint2(li,lmax,qmagpb,qmagb,lms,ab,is)
c                 ----------------------------------------------
               end if
           end do
        end do
           write(8,'('' Sum :'')')
           do li=1,nintfc
              cpalay=1.d0-conc(li).gt.tiny
               if(cpalay) then
                  ab='A'
               else
                  ab=' '
               end if
c              ---------------------------------------------
               call nelprint2(li,lmax,qmagpa,qmaga,lms,ab,0)
c              ---------------------------------------------
               if(cpalay) then
                  ab='B'
c                 ----------------------------------------------
                  call nelprint2(li,lmax,qmagpb,qmagb,lms,ab,0)
c                 ----------------------------------------------
               end if
           end do
      end if
c
      return
 150  format(/t15,
     >          '      s    ',
     >          '    p 1/2  ',
     >          '    p 3/2  ',
     >          '    d 3/2  ',
     >          '    d 5/2  ',
     >          '    total  ')
 160  format(t20,
     >          '         s     ',
     >          '         p     ',
     >          '         d     ',
     >          '       total   ')
      end
c     ========================================
      subroutine nelprint1(li,kmax,qvp)
c
c     ========================================
c
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
      dimension qvp(kmymaxp)
      dimension ppp(0:lmaxp,2)
c
      write(8,'('' Layer'',i4)') li
      kmy0=1
      do k=1,kmax
         l=k/2
         if(2*l.eq.k) then
             kap=l
         else
             kap=-l-1
         end if
         nk=2*iabs(kap)
         kmy1=kmy0+nk-1
         write(8,'(t10,6f12.8)') (qvp(kmy),kmy=kmy0,kmy1)
         kmy0=kmy1+1
      end do
c
      return
      end
c
c     ==============================================
      subroutine nelprint2(li,lmax,qvp,qv,lms,ab,is)
c
c     ==============================================
c
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
      integer mmax
      parameter (mmax=mimp)
c
      logical lms
      character*1 ab
      dimension qv(mmax)
      dimension qvp(kmymaxp,mmax)
      dimension ddd(kmaxp),ppp(0:lmaxp)
c
      kmax=2*lmax+1
      lmmax=(lmax+1)*(lmax+1)
      kmymax=2*lmmax
c
      if(.not.lms) then
        kmy=0
        do k=1,kmax
          l=k/2
          if(2*l.eq.k) then
              kap=l
          else
              kap=-l-1
          end if
          nk=2*iabs(kap)
          ddd(k)=0.d0
          do j=1,nk
              kmy=kmy+1
              ddd(k)=ddd(k)+qvp(kmy,li)
          end do
        end do
        if(ab.eq.'A'.or.ab.eq.' ') then
        write(8,'('' Layer'',i4,2x,a1,2x,t15,6f12.8)') li,
     >              ab,(ddd(k),k=1,kmax),qv(li)
        else
        write(8,'(''      '',4x,a1,2x,t15,6f12.8)')
     >              ab,(ddd(k),k=1,kmax),qv(li)
        end if
      else
        sum=0.d0
        if(is.le.1) ii=0
        if(is.eq.2) ii=lmmax
        do l=0,lmax
           ppp(l)=0.d0
           do m=-l,l
              ii=ii+1
              sum=sum+qvp(ii,li)
              ppp(l)=ppp(l)+qvp(ii,li)
              if(is.eq.0) sum=sum+qvp(ii+lmmax,li)
              if(is.eq.0) ppp(l)=ppp(l)+qvp(ii+lmmax,li)
           end do
        end do
        if(ab.eq.'A'.or.ab.eq.' ') then
        write(8,'('' Layer'',i4,2x,a1,2x,t20,4f15.10)') li,
     >              ab,(ppp(l),l=0,lmax),sum
        else
        write(8,'(''      '',4x,a1,2x,t20,4f15.10)')
     >              ab,(ppp(l),l=0,lmax),sum
        end if
      end if
c
      return
      end
