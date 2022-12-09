c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine printres1(
     >     dosout,pdosout,it,lmax,nimp,ne,cear,efermi,
     >     dosimp,doshimp,dosehimp,dosheimp,dostehimp,dostheimp,
     >     qvpimp,qvhpimp,
     >     qvimp,qvhimp,qvehimp,qvheimp,qvtehimp,qvtheimp,
     >     dos,lms)
!     >           for007,for008,it,lmax,nintfc,ne,cear,efermi,
!     >           dosa,qvpa,qva,dosmaga,
!     >           spin_magvpa,spin_magva,
!     >           orb_magvpa,orb_magva,
!     >           arot,qmoma,za,qca,
!     >           qvdiffa,enbdiffa,
!     >           entota,entot,entotifc,
!     >           idpota,vra,bra,bopra,
!     >           rs,dx,ns,lms,orbpol)
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
c      integer nimp
c      parameter (mmax=mimp)
c
      logical dos,lms
      logical cpasum
      logical orbpol
c
c      character*30 for007
      character*30 dosout
      character*34 pdosout
      character*10 idpota(nimp)
      character*1 ab
c
c      real*8 vra(nrad,mmax)
c      real*8 bra(nrad,mmax)
c      real*8 bopra(nrad,2,mmax)
c      real*8 rs(mmax)
c      real*8 dx(mmax)
      real*8 dosimp(kmymaxp,nimp,me)
      real*8 doshimp(kmymaxp,nimp,me)
      real*8 dosehimp(kmymaxp,nimp,me)
      real*8 dosheimp(kmymaxp,nimp,me)
      real*8 dostehimp(kmymaxp,nimp,me)
      real*8 dostheimp(kmymaxp,nimp,me)
      real*8 qvimp(nimp),qvhimp(nimp),qvehimp(nimp)
      real*8 qvheimp(nimp),qvtehimp(nimp),qvtheimp(nimp)
      real*8 qvpimp(kmymaxp,nimp),qvhpimp(kmymaxp,nimp)
c
      integer ii
c
      complex*16 cear(me)
c
      common/test/itest
      common/madelung_energy/emad
c
      data tiny/1.0d-6/
      data tol/1.0d-10/
c
      kmax=2*lmax+1
      lmmax=(lmax+1)*(lmax+1)
      kmymax=2*lmmax
      lmaxs=2*lmax
      lmmaxs=(lmaxs+1)*(lmaxs+1)
c
c - Print potential in file for007
c
c      open(7,file=for007,status='unknown')
c     rewind(7)
c      write(7,'(''SCI: '',i4,2x)') it
c
c      zeromt=0.d0
c      do li=1,nintfc
c         write(7,'(a)') idpota(li)
c         write(7,'(f8.5)') 1.00000
c         write(7,'(3f8.3,i5,f18.12)') 
c     >                  za(li),zeromt,dx(li),ns(li),rs(li)
c         write(7,'(4d20.12)') (vra(irad,li),irad=1,ns(li))
c         write(7,'(4d20.12)') (bra(irad,li),irad=1,ns(li))
c         if(orbpol) then
c           write(7,'(4d20.12)') (bopra(irad,1,li),irad=1,ns(li))
c           write(7,'(4d20.12)') (bopra(irad,2,li),irad=1,ns(li))
c         end if
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
c      end do
c      close(7)
c
c - Print DOS
c
      open(8,file=dosout,status='unknown')
      if(dos) open(9,file=pdosout,status='unknown')
c     rewind(8)
      write(8,'(''SCI: '',i4,2x)') it
      if(dos) write(9,'(''SCI: '',i4,2x)') it
c
      do li=1,nimp
        write(8,'(''##'',i3,t10,''DOS   e'')') li
        if(dos) write(9,'(''##'',i3,t10,''PDOS   e'')') li
c       -------------------------------------------
        call dosprint1(nimp,lms,lmax,li,ne,cear,dosimp)
        if(dos) call pdosprint1(nimp,lms,lmax,li,ne,cear,dosimp)
c       -------------------------------------------
        write(8,'(''##'',i3,t10,''DOS   h'')') li
        if(dos) write(9,'(''##'',i3,t10,''PDOS   h'')') li
c       -------------------------------------------
        call dosprint1(nimp,lms,lmax,li,ne,cear,doshimp)
        if(dos) call pdosprint1(nimp,lms,lmax,li,ne,cear,doshimp)
c       -------------------------------------------
        write(8,'(''##'',i3,t10,''DOS   s'')') li
        if(dos) write(9,'(''##'',i3,t10,''PDOS   s'')') li
c       -------------------------------------------
        call dosprint1(nimp,lms,lmax,li,ne,cear,dosehimp)
        if(dos) call pdosprint1(nimp,lms,lmax,li,ne,cear,dosehimp)
c       -------------------------------------------
        write(8,'(''##'',i3,t10,''DOS   t0(z)'')') li
        if(dos) write(9,'(''##'',i3,t10,''PDOS   t0(z)'')') li
c       -------------------------------------------
        call dosprint1(nimp,lms,lmax,li,ne,cear,dosheimp)
        if(dos) call pdosprint1(nimp,lms,lmax,li,ne,cear,dosheimp)
c       -------------------------------------------
        write(8,'(''##'',i3,t10,''DOS   tup(x)'')') li
        if(dos) write(9,'(''##'',i3,t10,''PDOS   tup(x)'')') li
c       -------------------------------------------
        call dosprint1(nimp,lms,lmax,li,ne,cear,dostehimp)
        if(dos) call pdosprint1(nimp,lms,lmax,li,ne,cear,dostehimp)
c       -------------------------------------------
        write(8,'(''##'',i3,t10,''DOS   tdown(y)'')') li
        if(dos) write(9,'(''##'',i3,t10,''PDOS   tdown(y)'')') li
c       -------------------------------------------
        call dosprint1(nimp,lms,lmax,li,ne,cear,dostheimp)
        if(dos) call pdosprint1(nimp,lms,lmax,li,ne,cear,dostheimp)
c       -------------------------------------------
c       if(cpalay) then
c         write(8,'(''##'',i3,t10,''DOS   B'')') li
c         -------------------------------------------
c         call dosprint(lms,lmax,li,ne,cear,dosb)
c         -------------------------------------------
c       end if
      end do
c
c -Print MDOS
c
c      write(8,*)
c      do li=1,nintfc
c       cpalay=1.d0-conc(li).gt.tiny
c        write(8,'(''##'',i3,t10,''MDOS   A'')') li
c       ----------------------------------------------
c        call dosprint(lms,lmax,li,ne,cear,dosmaga)
c       ----------------------------------------------
c       if(cpalay) then
c         write(8,'(''##'',i3,t10,''MDOS   B'')') li
c         ----------------------------------------------
c         call dosprint(lms,lmax,li,ne,cear,dosmagb)
c         ----------------------------------------------
c       end if
c      end do
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
         do li=1,nimp
c          --------------------------------------
           call nelprint1(li,kmax,qvpimp(1,li),8)
c          --------------------------------------
           call nelprint1(li,kmax,qvhpimp(1,li),8)
c          --------------------------------------
         end do
c         do li=1,nimp
c               ab=' '
c          ------------------------------------------
c            call nelprint2(li,lmax,qvpa,qva,lms,ab,1)
c          ------------------------------------------
c         end do
      else
         write(8,160)
         do is=1,2
            write(8,'(/'' Spin:'',i2)') is
            do li=1,nimp
                  ab=' '
c              ------------------------------------------
               call nelprintimp2(nimp,li,lmax,qvpimp,qvimp,lms,ab,is)
c              ------------------------------------------
               call nelprintimp2(nimp,li,lmax,qvhpimp,qvimp,lms,ab,is)
c              ------------------------------------------
            end do
         end do
      end if
c
c - Print magnetization in the layers
c
c      if(dabs(arot-0.d0).gt.tol) write(8,'(/''WARNING!!'',/,
c     > ''Magnetic moments reference frame rotated by'',f12.8)') arot
c
c      write(8,'(/29x,''Spin magnetic moments :'',
c     >          /29x,''----------------------''/ 
c     >          /29x,''         X            '')') 
c     ------------------------------------------------
c      call printmag(nintfc,lmax,conc,
c     >             spin_magva(1,1),spin_magvpa(1,1,1),
c     >             spin_magva(1,1),spin_magvpa(1,1,1),.true.,cpasum)
c     ------------------------------------------------
c      write(8,'(/29x,''         Y            '')') 
c     ------------------------------------------------
c      call printmag(nintfc,lmax,conc,
c     >             spin_magva(1,2),spin_magvpa(1,1,2),
c     >             spin_magva(1,2),spin_magvpa(1,1,2),.true.,cpasum)
c     ------------------------------------------------
c      write(8,'(/29x,''         Z            '')') 
c     ------------------------------------------------
c      call printmag(nintfc,lmax,conc,
c     >             spin_magva(1,3),spin_magvpa(1,1,3),
c     >             spin_magva(1,3),spin_magvpa(1,1,3),.true.,cpasum)
c     ------------------------------------------------
c
c      write(8,'(/29x,''Orbital magnetic moments :'',
c     >          /29x,''----------------------''/ 
c     >          /29x,''         X            '')') 
c     ------------------------------------------------
c      call printmag(nintfc,lmax,conc,
c     >             orb_magva(1,1),orb_magvpa(1,1,1),
c     >             orb_magva(1,1),orb_magvpa(1,1,1),.true.,cpasum)
c     ------------------------------------------------
c      write(8,'(/29x,''         Y            '')') 
c     ------------------------------------------------
c      call printmag(nintfc,lmax,conc,
c     >             orb_magva(1,2),orb_magvpa(1,1,2),
c     >             orb_magva(1,2),orb_magvpa(1,1,2),.true.,cpasum)
c     ------------------------------------------------
c      write(8,'(/29x,''         Z            '')') 
c     ------------------------------------------------
c      call printmag(nintfc,lmax,conc,
c     >             orb_magva(1,3),orb_magvpa(1,1,3),
c     >             orb_magva(1,3),orb_magvpa(1,1,3),.true.,cpasum)
c     ------------------------------------------------
c
c - Print moments of charge densities (in global frame of reference)
c
c      write(8,'(//''Moments of charge densities:''
c     >           /''---------------------------''/)')
c      i=0
c      do l=0,2*lmax 
c      do m=-l,l     
c         write(8,'(/10x,''(l,m)'',2i3/)') l,m
c         i=i+1
!        if(cpasum) then
!           write(8,'(t26,''A'',t66,''B'')')
!        end if
c         do li=1,nintfc
!          cpalay=1.d0-conc(li).gt.tiny
!          if(cpalay) then
!            write(8,'(2x,'' L'',i2,t10,2d15.7,5x,2d15.7)') li,
!    >       qmoma(i,li),qmomb(i,li)
!          else
c             write(8,'(2x,'' L'',i4,t10,2d15.7)') li,qmoma(i,li)
!          end if
c         end do
c      end do
c      end do
c
c - Print total energies
c
c      write(8,'(//29x,''Total energy:''/,
c     >           29x,''--------------'',//,''Layer'')')
c      cpasum=.false.
c      etka=0.0
c      etkb=0.0
c      do li=1,nintfc
!       cpalay=1.d0-conc(li).gt.tiny
c        ab=' '
!       if(cpalay) ab='A'
c        etka=etka+entota(li)
c        write(8,'(1x,i4,1x,a1,1x,f15.7)') li,ab,entota(li)
!       if(cpalay) then
!          ab='B'
!          cpasum=cpasum.or.cpalay
!          etkb=etkb+entotb(li)
!          write(8,'(4x,a1,1x,f15.7)') ab,entotb(li)
!       end if
c      end do
!     if(cpasum) then
!     write(8,'(/'' A-type Ifc '',f15.7)') etka
!     write(8,'(/'' B-type Ifc '',f15.7)') etkb
!     write(8,'(/'' Madelung   '',f15.7)') emad
!     write(8,'(/'' Sum    Ifc '',f15.7)') entotifc+emad
!     else
c      write(8,'(/'' Ifc  '',e16.7)') entotifc+emad
!     end if
c
c - Print energy-resolved contributions to band energy  
c
c      write(8,'(//''Band-energy:''
c     >           /''------------''/)')
c      omint=0.d0
c      do ie=ne,1,-1 
c         om=0.d0
c         do li=1,nintfc
c           om=om+(enbdiffa(li,ie)-efermi*qvdiffa(li,ie))
!          om=om+conc(li)*(enbdiffa(li,ie)-efermi*qvdiffa(li,ie))
!    >   +(1.d0-conc(li))*(enbdiffb(li,ie)-efermi*qvdiffb(li,ie))
c         end do
c         omint=omint+om
c         write(8,'(i4,2x,2d15.7,5x,2d17.7)') ie,cear(ie),om,omint
c      end do
c
      close(8)
      if(dos) close(9)
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
      subroutine dosprint1(nimp,lms,lmax,li,ne,cear,dos)
c
c     ================================================
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
      logical lms
c
      real*8 dos(kmymaxp,nimp,me)
      real*8 ddd(kmaxp),ppp(0:lmaxp,2)
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
          write(8,'(f16.10,6e18.8)') er,(ddd(k),k=1,kmax),sum
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
             write(8,'(f16.10,6e18.8)') er,(ppp(l,is),l=0,lmax),sum
          end do
c         write(8,'(f14.10,32e15.8)') er,(dos(l,li,ie),l=1,2*lmmax)
        end if
      end do
c
      return
      end
c
c
c     ============================================
      subroutine pdosprint1(nimp,lms,lmax,li,ne,cear,dos)
c     ============================================
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
      logical lms
c
      dimension dos(kmymaxp,nimp,me)
      dimension ddd(kmaxp),ppp(0:lmmaxp,2)
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
          write(9,'(f14.10,2x,19f15.8)') er,(dos(k,li,ie),k=1,kmymax),sum
        else
          ii=0
          do is=1,2
             iip=0
             sum=0.d0
             do l=0,lmax
                do m=-l,l
                   ii=ii+1
                   iip=iip+1
                   sum=sum+dos(ii,li,ie)
                   ppp(iip,is)=dos(ii,li,ie)
                end do
             end do
             write(9,'(f14.10,2x,19f15.8)') er,(ppp(l,is),l=0,lmmax),sum
          end do
        end if
      end do
c
      return
      end
c     =========================================================
c      subroutine printmag(nintfc,lmax,
c     >             conc,qmaga,qmagpa,qmagb,qmagpb,lms,cpasum)
c
c     =========================================================
c
c      implicit real*8 (a-h,o-z)
c      include '../param.h'
c
c      integer mmax
c      parameter (mmax=mimp)
c
c      logical lms,cpasum,cpalay
c      character*1 ab
c      dimension qmagpa(kmymaxp,mmax),qmaga(mmax)
c      dimension qmagpb(kmymaxp,mmax),qmagb(mmax)
c      dimension conc(mmax)
c
c      data tiny/1.0d-6/
c
c      kmax=2*lmax+1
c      lmmax=(lmax+1)*(lmax+1)
c      kmymax=2*lmmax
c
c      if(.not.lms) then
c        if(cpasum) then
c          write(8,'(2x,'' -- A type --'')')
c        end if
c        do li=1,nintfc
c          ------------------------------------
c           call nelprint1(li,kmax,qmagpa(1,li))
c          ------------------------------------
c        end do
c        if(cpasum) then
c          write(8,'(2x,'' -- B type --'')')
c          do li=1,nintfc
c             cpalay=1.d0-conc(li).gt.tiny
c             if(cpalay) then
c              --------------------------------
c               call nelprint1(li,kmax,qmagpb(1,li))
c              --------------------------------
c             end if
c          end do
c        end if
c        write(8,150)
c        do li=1,nintfc
c           cpalay=1.d0-conc(li).gt.tiny
c           if(cpalay) then
c              ab='A'
c           else
c              ab=' '
c           end if
c          ---------------------------------------------
c           call nelprint2(li,lmax,qmagpa,qmaga,lms,ab,1)
c          ---------------------------------------------
c           if(cpalay) then
c              ab='B'
c             ---------------------------------------------
c              call nelprint2(li,lmax,qmagpb,qmagb,lms,ab,1)
c             ---------------------------------------------
c           end if
c        end do
c      else
c        write(8,160)
c        do is=1,2
c           write(8,'('' Spin:'',i2)') is
c           do li=1,nintfc
c              cpalay=1.d0-conc(li).gt.tiny
c               if(cpalay) then
c                  ab='A'
c               else
c                  ab=' '
c               end if
c              ---------------------------------------------
c               call nelprint2(li,lmax,qmagpa,qmaga,lms,ab,is)
c              ---------------------------------------------
c               if(cpalay) then
c                  ab='B'
c                 ----------------------------------------------
c                  call nelprint2(li,lmax,qmagpb,qmagb,lms,ab,is)
c                 ----------------------------------------------
c               end if
c           end do
c        end do
c           write(8,'('' Sum :'')')
c           do li=1,nintfc
c              cpalay=1.d0-conc(li).gt.tiny
c               if(cpalay) then
c                  ab='A'
c               else
c                  ab=' '
c               end if
c              ---------------------------------------------
c               call nelprint2(li,lmax,qmagpa,qmaga,lms,ab,0)
c              ---------------------------------------------
c               if(cpalay) then
c                  ab='B'
c                 ----------------------------------------------
c                  call nelprint2(li,lmax,qmagpb,qmagb,lms,ab,0)
c                 ----------------------------------------------
c               end if
c           end do
c      end if
c
c      return
c 150  format(/t15,
c     >          '      s    ',
c     >          '    p 1/2  ',
c     >          '    p 3/2  ',
c     >          '    d 3/2  ',
c     >          '    d 5/2  ',
c     >          '    total  ')
c 160  format(t20,
c     >          '         s     ',
c     >          '         p     ',
c     >          '         d     ',
c     >          '       total   ')
c      end
c     ========================================
c      subroutine nelprint1(li,kmax,qvp)
c
c     ========================================
c
c      implicit real*8 (a-h,o-z)
c      include '../param.h'
c
c      dimension qvp(kmymaxp)
c      dimension ppp(0:lmaxp,2)
c
c      write(8,'('' Layer'',i4)') li
c      kmy0=1
c      do k=1,kmax
c         l=k/2
c         if(2*l.eq.k) then
c             kap=l
c         else
c             kap=-l-1
c         end if
c         nk=2*iabs(kap)
c         kmy1=kmy0+nk-1
c         write(8,'(t10,6f11.7)') (qvp(kmy),kmy=kmy0,kmy1)
c         kmy0=kmy1+1
c      end do
c
c      return
c      end
c
c     ==============================================
      subroutine nelprintimp2(nimp,li,lmax,qvp,qv,lms,ab,is)
c
c     ==============================================
c
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
      integer nimp
c      parameter (mmax=mimp)
c
      logical lms
      character*1 ab
      dimension qv(nimp)
      dimension qvp(kmymaxp,nimp)
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
        write(8,'('' Layer'',i4,2x,a1,2x,t15,6f11.7)') li,
     >              ab,(ddd(k),k=1,kmax),qv(li)
        else
        write(8,'(''      '',4x,a1,2x,t15,6f11.7)')
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
