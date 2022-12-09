c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine printres(
     >           dosout,pdosout,potout,momout,
     >           it,lmax,nintfc,ne,cear,efermi,
     >           conc,dosa,dosb,dosha,doshb,
     >           doseha,doshea,dosteha,dosthea,
     >           dosehb,dosheb,dostehb,dostheb,
     >           qvpa,qvpb,qva,qvb,dosmaga,dosmagb,dosmagha,
     >           spin_magvpa,spin_magvpb,spin_magva,spin_magvb,
     >           orb_magvpa,orb_magvpb,orb_magva,orb_magvb,
     >           arot,qmoma,qmomb,za,zb,qca,qcb,
     >           qvdiffa,qvdiffb,enbdiffa,enbdiffb,
     >           entota,entotb,entot,entotifc,
     >           idpota,idpotb,vra,bra,bopra,vrb,brb,boprb,
     >           rs,dx,ns,vrsbulk,lms,orbpol,bulk)
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
      logical lms,cpalay,cpasum,orbpol,bulk
c
      character*30 dosout,pdosout,potout,momout
      character*10 idpota(mintfc),idpotb(mintfc)
      character*1 ab
c
      dimension vra(nrad,mintfc),bra(nrad,mintfc)
      dimension vrb(nrad,mintfc),brb(nrad,mintfc)
      dimension bopra(nrad,2,mintfc),boprb(nrad,2,mintfc)
      dimension rs(mintfc),dx(mintfc),ns(mintfc)
      dimension dosa(kmymaxp,mintfc,me),dosb(kmymaxp,mintfc,me)
      dimension dosha(kmymaxp,mintfc,me),doshb(kmymaxp,mintfc,me)
      dimension doseha(kmymaxp,mintfc,me),dosehb(kmymaxp,mintfc,me)
      dimension doshea(kmymaxp,mintfc,me),dosheb(kmymaxp,mintfc,me)
      dimension dosteha(kmymaxp,mintfc,me),dostehb(kmymaxp,mintfc,me)
      dimension dosthea(kmymaxp,mintfc,me),dostheb(kmymaxp,mintfc,me)
      dimension dosmaga(kmymaxp,mintfc,me),dosmagb(kmymaxp,mintfc,me)
      dimension dosmagha(kmymaxp,mintfc,me),dosmaghb(kmymaxp,mintfc,me)
      dimension qvdiffa(mintfc,me),qvdiffb(mintfc,me)
      dimension enbdiffa(mintfc,me),enbdiffb(mintfc,me)
      dimension qva(mintfc),qvb(mintfc)
      dimension qca(mintfc),qcb(mintfc)
      dimension qvpa(kmymaxp,mintfc),qvpb(kmymaxp,mintfc)
      dimension za(mintfc),zb(mintfc)
      dimension spin_magvpa(kmymaxp,mintfc,3)
      dimension spin_magvpb(kmymaxp,mintfc,3)
      dimension spin_magva(mintfc,3),spin_magvb(mintfc,3)
      dimension orb_magvpa(kmymaxp,mintfc,3)
      dimension orb_magvpb(kmymaxp,mintfc,3)
      dimension orb_magva(mintfc,3),orb_magvb(mintfc,3)
      dimension entota(mintfc),entotb(mintfc),entot(mintfc)
      dimension conc(mintfc)
c
      complex*16 cear(me)
      complex*16 qmoma(lmsup,mintfc),qmomb(lmsup,mintfc)
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
c - Print potential in file potout
c
      open(7,file=potout,status='unknown')
      write(7,'(''SCI: '',i4)') it
      do li=1,nintfc
         cpalay=1.d0-conc(li).gt.tiny
         write(7,'(a)') idpota(li)
         write(7,'(f8.5)') conc(li)
         write(7,'(f8.3,f15.10,d20.10,i5,f18.12)')
     >   za(li),vrsbulk,dx(li),ns(li),rs(li)
         write(7,'(4d20.12)') (vra(irad,li),irad=1,ns(li))
         write(7,'(4d20.12)') (bra(irad,li),irad=1,ns(li))
         if(orbpol) then
           write(7,'(4d20.12)') (bopra(irad,1,li),irad=1,ns(li))
           write(7,'(4d20.12)') (bopra(irad,2,li),irad=1,ns(li))
         end if
         if(cpalay) then
            write(7,'(a)') idpotb(li)
            write(7,'(f8.5)') 1.0-conc(li)
            write(7,'(f8.3,f15.10,d20.10,i5,f18.12)')
     >      zb(li),vrsbulk,dx(li),ns(li),rs(li)
            write(7,'(4d20.12)') (vrb(irad,li),irad=1,ns(li))
            write(7,'(4d20.12)') (brb(irad,li),irad=1,ns(li))
            if(orbpol) then
              write(7,'(4d20.12)') (boprb(irad,1,li),irad=1,ns(li))
              write(7,'(4d20.12)') (boprb(irad,2,li),irad=1,ns(li))
            end if
         end if
      end do
      close(7)
c
c - Store moments in file momout
c
      open(7,file=momout,status='unknown')
      write(7,'(''SCI: '',i4)') it
      write(7,'(i3)') lmax
      do li=1,nintfc
         cpalay=1.d0-conc(li).gt.tiny
         write(7,'(a)') idpota(li)
         write(7,'(f8.5)') conc(li)
         write(7,'(4d20.12)') qmoma(1,li)+qca(li)-za(li),
     >   (qmoma(i,li),i=2,lmmaxs)
         if(cpalay) then
         write(7,'(a)') idpotb(li)
         write(7,'(f8.5)') 1.d0-conc(li)
         write(7,'(4d20.12)') qmomb(1,li)+qcb(li)-zb(li),
     >   (qmomb(i,li),i=2,lmmaxs)
         end if
      end do
      close(7)
c
c - Print DOS
c
      open(7,file=dosout,status='unknown')
      write(7,'(''SCI: '',i4,2x)') it
      open(8,file=pdosout,status='unknown')
      write(8,'(''SCI: '',i4,2x)') it
c
      do li=1,nintfc
        cpalay=1.d0-conc(li).gt.tiny
        write(7,'(''##'',i3,t10,''DOS electron   A'')') li
        write(8,'(''##'',i3,t10,a,''PDOS electron   A'')') li,idpotb(li)
c       ---------------------------------------
        call dosprint(lms,lmax,li,ne,cear,dosa)
        call pdosprint(lms,lmax,li,ne,cear,dosa) 
c       ---------------------------------------
        if(cpalay) then
          write(7,'(''##'',i3,t10,''DOS electron  B'')') li
          write(8,'(''##'',i3,t10,a,''PDOS electron B'')') li,idpotb(li)
c         ---------------------------------------
          call dosprint(lms,lmax,li,ne,cear,dosb)
          call pdosprint(lms,lmax,li,ne,cear,dosb) 
c         ---------------------------------------
        end if
      end do
c
c -- DOS hole
c
      do li=1,nintfc
        cpalay=1.d0-conc(li).gt.tiny
        write(7,'(''##'',i3,t10,''DOS hole   A'')') li
        write(8,'(''##'',i3,t10,a,''PDOS hole   A'')') li,idpotb(li)
c       ---------------------------------------
        call dosprint(lms,lmax,li,ne,cear,dosha)
        call pdosprint(lms,lmax,li,ne,cear,dosha)
c       ---------------------------------------
        if(cpalay) then
          write(7,'(''##'',i3,t10,''DOS hole  B'')') li
          write(8,'(''##'',i3,t10,a,''PDOS hole B'')') li,idpotb(li)
c         ---------------------------------------
          call dosprint(lms,lmax,li,ne,cear,doshb)
          call pdosprint(lms,lmax,li,ne,cear,doshb) 
c         ---------------------------------------
        end if

      end do
c
c -- DOS singlet
c
      do li=1,nintfc
       cpalay=1.d0-conc(li).gt.tiny
       write(7,'(''##'',i3,t10,''DOS singlet A'')') li
       write(8,'(''##'',i3,t10,a,''PDOS singlet A'')') li,idpotb(li)
c      ---------------------------------------
       call dosprint(lms,lmax,li,ne,cear,doseha)
       call pdosprint(lms,lmax,li,ne,cear,doseha)
c      ---------------------------------------
       if (cpalay) then
         write(7,'(''##'',i3,t10,''DOS singlet B'')') li
         write(8,'(''##'',i3,t10,a,''PDOS singlet B'')') li,idpotb(li)
c        ---------------------------------------
         call dosprint(lms,lmax,li,ne,cear,dosehb)
         call pdosprint(lms,lmax,li,ne,cear,dosehb)
c        ---------------------------------------
      end if
      end do
c
      if(itest.ge.2) then
c -- DOS T0
c
      do li=1,nintfc
        cpalay=1.d0-conc(li).gt.tiny
        write(7,'(''##'',i3,t10,''DOS triplet S0'')') li
        write(8,'(''##'',i3,t10,a,''PDOS triplet S0'')') li,idpotb(li)
c       ---------------------------------------
        call dosprint(lms,lmax,li,ne,cear,doshea)
        call pdosprint(lms,lmax,li,ne,cear,doshea)
c       ---------------------------------------
        if (cpalay) then
          write(7,'(''##'',i3,t10,''DOS triplet S0'')') li
          write(8,'(''##'',i3,t10,a,''PDOS triplet S0'')') li,idpotb(li)
c         ---------------------------------------
          call dosprint(lms,lmax,li,ne,cear,dosheb)
          call pdosprint(lms,lmax,li,ne,cear,dosheb)
c         ---------------------------------------
        end if
      end do
c
c -- DOS T_up
c
      do li=1,nintfc
        cpalay=1.d0-conc(li).gt.tiny
        write(7,'(''##'',i3,t10,''DOS triplet up'')') li
        write(8,'(''##'',i3,t10,a,''PDOS triplet up'')') li, idpotb(li)
c       ---------------------------------------
        call dosprint(lms,lmax,li,ne,cear,dosteha)
        call pdosprint(lms,lmax,li,ne,cear,dosteha)
c       ---------------------------------------
        if (cpalay) then
          write(7,'(''##'',i3,t10,''DOS triplet up'')') li
          write(8,'(''##'',i3,t10,a,''PDOS triplet up'')') li,idpotb(li)
c         ---------------------------------------
          call dosprint(lms,lmax,li,ne,cear,dostehb)
          call pdosprint(lms,lmax,li,ne,cear,dostehb)
c         ---------------------------------------
        end if
      end do
c
c -- DOS T_down
c
      do li=1,nintfc
        cpalay=1.d0-conc(li).gt.tiny
        write(7,'(''##'',i3,t10,''DOS triplet down'')') li
        write(8,'(''##'',i3,t10,a,''PDOS down'')') li, idpotb(li)
c       ---------------------------------------
        call dosprint(lms,lmax,li,ne,cear,dosthea)
        call pdosprint(lms,lmax,li,ne,cear,dosthea)
c       ---------------------------------------
        if (cpalay) then
          write(7,'(''##'',i3,t10,''DOS triplet down'')') li
          write(8,'(''##'',i3,t10,a,''PDOS triplet down'')') 
     >          li,idpotb(li)
c         ---------------------------------------
          call dosprint(lms,lmax,li,ne,cear,dostheb)
          call pdosprint(lms,lmax,li,ne,cear,dostheb)
c         ---------------------------------------
        end if
      end do
      end if ! for itest .gt. 2
      close(8)
c
c
c -Print MDOS for electrons
c
      write(7,*)
      do li=1,nintfc
        cpalay=1.d0-conc(li).gt.tiny
        write(7,'(''##'',i3,t10,''MDOS electron-electron part  A'')') li
c       ------------------------------------------
        call dosprint(lms,lmax,li,ne,cear,dosmaga)
c       ------------------------------------------
        if(cpalay) then
          write(7,'(''##'',i3,t10,''MDOS electron-electron part  B'')') li
c         ------------------------------------------
          call dosprint(lms,lmax,li,ne,cear,dosmagb)
c         ------------------------------------------
        end if
      end do
c
c -Print MDOS for holes
c
      write(7,*)
      do li=1,nintfc
        cpalay=1.d0-conc(li).gt.tiny
        write(7,'(''##'',i3,t10,''MDOS hole-hole part  A'')') li
c       ------------------------------------------
        call dosprint(lms,lmax,li,ne,cear,dosmagha)
c       ------------------------------------------
        if(cpalay) then
          write(7,'(''##'',i3,t10,''MDOS hole-hole part  B'')') li
c         ------------------------------------------
          call dosprint(lms,lmax,li,ne,cear,dosmaghb)
c         ------------------------------------------
        end if
      end do
c
c - Print No. of electrons in the layers
c
      cpasum=.false.
      do li=1,nintfc
         cpalay=1.d0-conc(li).gt.tiny
         cpasum=cpasum.or.cpalay
      end do
c
      if(cpasum) then
        write(7,'(/2x,''No. of electrons  - A type - :''/)') 
      else
        write(7,'(/2x,''No. of electrons :''/)')
      end if
c
      if(.not.lms) then
         do li=1,nintfc
c          --------------------------------------
           call nelprint1(li,kmax,qvpa(1,li),7)
c          --------------------------------------
         end do
c        -------------------- do B-type if there was cpalay ------
         if(cpasum) then
           write(7,'(/2x,''No. of electrons  - B type - :''/)')
           do li=1,nintfc
              cpalay=1.d0-conc(li).gt.tiny
              if(cpalay) then
c                --------------------------------------
                 call nelprint1(li,kmax,qvpb(1,li),7)
c                --------------------------------------
              end if
           end do
         end if
         if(lmax.eq.0) write(7,150)
         if(lmax.eq.1) write(7,151)
         if(lmax.eq.2) write(7,152)
         if(lmax.eq.3) write(7,153)
         if(lmax.eq.4) write(7,154)
         do li=1,nintfc
            cpalay=1.d0-conc(li).gt.tiny
            if(cpalay) then
               ab='A'
            else
               ab=' '
            end if
c          ------------------------------------------
            call nelprint2(li,lmax,qvpa,qva,lms,ab,1)
c          ------------------------------------------
            if(cpalay) then
               ab='B'
c             -----------------------------------------
              call nelprint2(li,lmax,qvpb,qvb,lms,ab,1)
c             -----------------------------------------
            end if
         end do
      else
         if(lmax.eq.0) write(7,160)
         if(lmax.eq.1) write(7,161)
         if(lmax.eq.2) write(7,162)
         if(lmax.eq.3) write(7,163)
         if(lmax.eq.4) write(7,164)
         do is=1,2
            write(7,'(/'' Spin:'',i2)') is
            do li=1,nintfc
               cpalay=1.d0-conc(li).gt.tiny
               if(cpalay) then
                  ab='A'
               else
                  ab=' '
               end if
c              ------------------------------------------
               call nelprint2(li,lmax,qvpa,qva,lms,ab,is)
c              ------------------------------------------
               if(cpalay) then
                  ab='B'
c                 ------------------------------------------
                  call nelprint2(li,lmax,qvpb,qvb,lms,ab,is)
c                 ------------------------------------------
               end if
            end do
         end do
      end if
c
c - Print magnetization in the layers
c
      if(dabs(arot-0.d0).gt.tol) write(7,'(/''WARNING!!'',/,
     > ''Magnetic moments reference frame rotated by'',f12.8)') arot
c
      write(7,'(/29x,''Spin magnetic moments (just electron) :'',
     >          /29x,''----------------------''/ 
     >          /29x,''         X            '')') 
c     ------------------------------------------------
      call printmag(nintfc,lmax,conc,
     >             spin_magva(1,1),spin_magvpa(1,1,1),
     >             spin_magvb(1,1),spin_magvpb(1,1,1),.true.,cpasum)
c     ------------------------------------------------
      write(7,'(/29x,''         Y            '')') 
c     ------------------------------------------------
      call printmag(nintfc,lmax,conc,
     >             spin_magva(1,2),spin_magvpa(1,1,2),
     >             spin_magvb(1,2),spin_magvpb(1,1,2),.true.,cpasum)
c     ------------------------------------------------
      write(7,'(/29x,''         Z            '')') 
c     ------------------------------------------------
      call printmag(nintfc,lmax,conc,
     >             spin_magva(1,3),spin_magvpa(1,1,3),
     >             spin_magvb(1,3),spin_magvpb(1,1,3),.true.,cpasum)
c     ------------------------------------------------
c
      write(7,'(/29x,''Orbital magnetic moments (just electron) :'',
     >          /29x,''----------------------''/ 
     >          /29x,''         X            '')') 
c     ------------------------------------------------
      call printmag(nintfc,lmax,conc,
     >             orb_magva(1,1),orb_magvpa(1,1,1),
     >             orb_magvb(1,1),orb_magvpb(1,1,1),.true.,cpasum)
c     ------------------------------------------------
      write(7,'(/29x,''         Y            '')') 
c     ------------------------------------------------
      call printmag(nintfc,lmax,conc,
     >             orb_magva(1,2),orb_magvpa(1,1,2),
     >             orb_magvb(1,2),orb_magvpb(1,1,2),.true.,cpasum)
c     ------------------------------------------------
      write(7,'(/29x,''         Z            '')') 
c     ------------------------------------------------
      call printmag(nintfc,lmax,conc,
     >             orb_magva(1,3),orb_magvpa(1,1,3),
     >             orb_magvb(1,3),orb_magvpb(1,1,3),.true.,cpasum)
c     ------------------------------------------------
c
c - Print moments of charge densities (in global frame of reference)
c
      write(7,'(//''Moments of charge densities (just electron) :''
     >           /''---------------------------''/)')
      i=0
      do l=0,2*lmax 
      do m=-l,l     
         write(7,'(/10x,''(l,m)'',2i3/)') l,m
         i=i+1
         if(cpasum) then
            write(7,'(t26,''A'',t66,''B'')')
         end if
         do li=1,nintfc
           cpalay=1.d0-conc(li).gt.tiny
           if(cpalay) then
             write(7,'(2x,'' L'',i2,t10,2d15.7,5x,2d15.7)') li,
     >       qmoma(i,li),qmomb(i,li)
           else
             write(7,'(2x,'' L'',i2,t10,2d15.7)') li,qmoma(i,li)
           end if
         end do
      end do
      end do
c
c - Print total energies
c
      write(7,'(//29x,''Total energy (just electron):''/,
     >           29x,''--------------'',//,''Layer'')')
      cpasum=.false.
      etka=0.0
      etkb=0.0
      do li=1,nintfc
        cpalay=1.d0-conc(li).gt.tiny
        ab=' '
        if(cpalay) ab='A'
        etka=etka+entota(li)
        write(7,'(1x,i2,1x,a1,1x,f15.7)') li,ab,entota(li)
        if(cpalay) then
           ab='B'
           cpasum=cpasum.or.cpalay
           etkb=etkb+entotb(li)
           write(7,'(4x,a1,1x,f15.7)') ab,entotb(li)
        end if
      end do
      if(cpasum) then
      write(7,'(/'' A-type Ifc '',f15.7)') etka
      write(7,'(/'' B-type Ifc '',f15.7)') etkb
      write(7,'(/'' Madelung   '',f15.7)') emad
      write(7,'(/'' Sum    Ifc '',f15.7)') entotifc+emad
      else
      write(7,'(/'' Ifc  '',f15.7)') entotifc+emad
      end if
c
c - Print energy-resolved contributions to band energy  
c
      write(7,'(//''Band-energy:''
     >           /''------------''/)')
      omint=0.d0
      do ie=ne,1,-1 
         om=0.d0
         do li=1,nintfc
           om=om+conc(li)*(enbdiffa(li,ie)-efermi*qvdiffa(li,ie))
     >   +(1.d0-conc(li))*(enbdiffb(li,ie)-efermi*qvdiffb(li,ie))
         end do
         omint=omint+om
         write(7,'(i4,2x,2d15.7,5x,2d22.12)') ne-ie+1,cear(ie),om,omint
      end do
c
      close(7)
c
      return
 150  format(/t15,
     >          '    s 1/2  ',
     >          '    total  ')
 151  format(/t15,
     >          '    s 1/2  ',
     >          '    p 1/2  ',
     >          '    p 3/2  ',
     >          '    total  ')
 152  format(/t15,
     >          '    s 1/2  ',
     >          '    p 1/2  ',
     >          '    p 3/2  ',
     >          '    d 3/2  ',
     >          '    d 5/2  ',
     >          '    total  ')
 153  format(/t15,
     >          '    s 1/2  ',
     >          '    p 1/2  ',
     >          '    p 3/2  ',
     >          '    d 3/2  ',
     >          '    d 5/2  ',
     >          '    f 5/2  ',
     >          '    f 7/2  ',
     >          '    total  ')
 154  format(/t15,
     >          '    s 1/2  ',
     >          '    p 1/2  ',
     >          '    p 3/2  ',
     >          '    d 3/2  ',
     >          '    d 5/2  ',
     >          '    f 5/2  ',
     >          '    f 7/2  ',
     >          '    g 7/2  ',
     >          '    g 9/2  ',
     >          '    total  ')
 160  format(t15,
     >          '     s     ',
     >          '   total   ')
 161  format(t15,
     >          '     s     ',
     >          '     p     ',
     >          '   total   ')
 162  format(t15,
     >          '     s     ',
     >          '     p     ',
     >          '     d     ',
     >          '   total   ')
 163  format(t15,
     >          '     s     ',
     >          '     p     ',
     >          '     d     ',
     >          '     f     ',
     >          '   total   ')
 164  format(t15,
     >          '     s     ',
     >          '     p     ',
     >          '     d     ',
     >          '     f     ',
     >          '     g     ',
     >          '   total   ')
      end
c     ============================================
      subroutine dosprint(lms,lmax,li,ne,cear,dos)
c     ============================================
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
      logical lms
c
      dimension dos(kmymaxp,mintfc,me)
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
          write(7,'(f16.10,4x,18e18.8)') er,(ddd(k),k=1,kmax),sum
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
             write(7,'(f16.10,4x,18e18.8)') er,(ppp(l,is),l=0,lmax),sum
          end do
        end if
      end do
c
      return
      end
c
c     ============================================
      subroutine pdosprint(lms,lmax,li,ne,cear,dos)
c     ============================================
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
      logical lms
c
      dimension dos(kmymaxp,mintfc,me)
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
          write(8,'(f14.10,2x,15f15.5)') er,(dos(k,li,ie),k=1,kmymax),sum
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
             write(8,'(f14.10,2x,15f15.5)') er,(ppp(l,is),l=0,lmmax),sum
          end do
        end if
      end do
c
      return
      end
c
c     =========================================================
      subroutine printmag(nintfc,lmax,
     >             conc,qmaga,qmagpa,qmagb,qmagpb,lms,cpasum)
c     =========================================================
c
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
      logical lms,cpasum,cpalay
      character*1 ab
      dimension qmagpa(kmymaxp,mintfc),qmaga(mintfc)
      dimension qmagpb(kmymaxp,mintfc),qmagb(mintfc)
      dimension conc(mintfc)
c
      data tiny/1.0d-6/
c
      kmax=2*lmax+1
      lmmax=(lmax+1)*(lmax+1)
      kmymax=2*lmmax
c
      if(.not.lms) then
        if(cpasum) then
          write(7,'(2x,'' -- A type --'')')
        end if
        do li=1,nintfc
c          ------------------------------------
           call nelprint1(li,kmax,qmagpa(1,li),7)
c          ------------------------------------
        end do
        if(cpasum) then
          write(7,'(2x,'' -- B type --'')')
          do li=1,nintfc
             cpalay=1.d0-conc(li).gt.tiny
             if(cpalay) then
c              --------------------------------
               call nelprint1(li,kmax,qmagpb(1,li),7)
c              --------------------------------
             end if
          end do
        end if
        if(lmax.eq.0) write(7,150)
        if(lmax.eq.1) write(7,151)
        if(lmax.eq.2) write(7,152)
        if(lmax.eq.3) write(7,153)
        if(lmax.eq.4) write(7,154)
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
        if(lmax.eq.0) write(7,160)
        if(lmax.eq.1) write(7,161)
        if(lmax.eq.2) write(7,162)
        if(lmax.eq.3) write(7,163)
        if(lmax.eq.4) write(7,164)
        do is=1,2
           write(7,'('' Spin:'',i2)') is
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
           write(7,'('' Sum :'')')
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
     >          '    s 1/2  ',
     >          '    total  ')
 151  format(/t15,
     >          '    s 1/2  ',
     >          '    p 1/2  ',
     >          '    p 3/2  ',
     >          '    total  ')
 152  format(/t15,
     >          '    s 1/2  ',
     >          '    p 1/2  ',
     >          '    p 3/2  ',
     >          '    d 3/2  ',
     >          '    d 5/2  ',
     >          '    total  ')
 153  format(/t15,
     >          '    s 1/2  ',
     >          '    p 1/2  ',
     >          '    p 3/2  ',
     >          '    d 3/2  ',
     >          '    d 5/2  ',
     >          '    f 5/2  ',
     >          '    f 7/2  ',
     >          '    total  ')
 154  format(/t15,
     >          '    s 1/2  ',
     >          '    p 1/2  ',
     >          '    p 3/2  ',
     >          '    d 3/2  ',
     >          '    d 5/2  ',
     >          '    f 5/2  ',
     >          '    f 7/2  ',
     >          '    g 7/2  ',
     >          '    g 9/2  ',
     >          '    total  ')
 160  format(t20,
     >          '     s     ',
     >          '   total   ')
 161  format(t20,
     >          '     s     ',
     >          '     p     ',
     >          '   total   ')
 162  format(t20,
     >          '     s     ',
     >          '     p     ',
     >          '     d     ',
     >          '   total   ')
 163  format(t20,
     >          '     s     ',
     >          '     p     ',
     >          '     d     ',
     >          '     f     ',
     >          '   total   ')
 164  format(t20,
     >          '     s     ',
     >          '     p     ',
     >          '     d     ',
     >          '     f     ',
     >          '     g     ',
     >          '   total   ')
      end
c     ========================================
      subroutine nelprint1(li,kmax,qvp,nf)
c     ========================================
c
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
      dimension qvp(kmymaxp)
      dimension ppp(0:lmaxp,2)
c
      write(nf,'('' Layer'',i2)') li
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
         write(nf,'(t10,10f12.8)') (qvp(kmy),kmy=kmy0,kmy1)
         kmy0=kmy1+1
      end do
c
      return
      end
c
c     ==============================================
      subroutine nelprint2(li,lmax,qvp,qv,lms,ab,is)
c     ==============================================
c
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
      logical lms
      character*1 ab
      dimension qv(mintfc)
      dimension qvp(kmymaxp,mintfc)
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
        write(7,'('' Layer'',i2,2x,a1,2x,t15,10f12.8)') li,
     >              ab,(ddd(k),k=1,kmax),qv(li)
        else
        write(7,'(''      '',4x,a1,2x,t15,10f12.8)')
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
        write(7,'('' Layer'',i2,2x,a1,2x,t15,10f12.8)') li,
     >              ab,(ppp(l),l=0,lmax),sum
        else
        write(7,'(''      '',4x,a1,2x,t15,10f12.8)')
     >              ab,(ppp(l),l=0,lmax),sum
        end if
      end if
c
      return
      end
