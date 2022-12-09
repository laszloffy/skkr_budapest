c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine printres(
     >potout,dosout,pdosout,momout,
     >it,lmax,nintfc,ne,cear,efermi,conc,dlm,bulk,
     >dosa,dosb,doslma,doslmb,qvpa,qvpb,qva,qvb,
     >qmagvpa,qmagvpb,qmagva,qmagvb,qmoma,qmomb,za,zb,qca,qcb,
     >enba,enbb,enkina,enkinb,enela,enelb,enxca,enxcb,
     >idpota,idpotb,vra,vrb,bra,brb,dx,ns,rs,vrsbulk)
c23456789012345678901234567890123456789012345678901234567890123456789012
c ========================
c
c Print out results 
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      logical cpalay,dlm,bulk
      character*30 potout,dosout,pdosout,momout
      character*10 idpota(mintfc),idpotb(mintfc)
c
      dimension vra(nrad,mintfc),vrb(nrad,mintfc)
      dimension bra(nrad,mintfc),brb(nrad,mintfc)
      dimension dx(mintfc),ns(mintfc),rs(mintfc)
      dimension dosa(0:lmaxp,mintfc,me,2),dosb(0:lmaxp,mintfc,me,2)
      dimension doslma(lmmaxp,mintfc,me,2),doslmb(lmmaxp,mintfc,me,2)
      dimension za(mintfc),zb(mintfc)
      dimension qva(mintfc),qvb(mintfc)
      dimension qca(mintfc),qcb(mintfc)
      dimension qmagva(mintfc),qmagvb(mintfc)
      dimension qvpa(0:lmaxp,mintfc),qvpb(0:lmaxp,mintfc)
      dimension qmagvpa(0:lmaxp,mintfc),qmagvpb(0:lmaxp,mintfc)
      dimension enba(mintfc),enbb(mintfc)
      dimension enkina(mintfc),enkinb(mintfc)
      dimension enela(mintfc),enelb(mintfc)
      dimension enxca(mintfc),enxcb(mintfc)
      dimension conc(mintfc)
c
      complex*16 cear(me)
      complex*16 qmoma(lmsup,mintfc),qmomb(lmsup,mintfc)
c
      common/madelung_energy/emad
c
      data tiny/1.0d-6/
c
      lmaxs=2*lmax
      lmmaxs=(lmaxs+1)*(lmaxs+1)
      lmmax=(lmax+1)*(lmax+1)
c
c - Print DOS
c
      open(8,file=dosout,status='unknown')
      open(9,file=pdosout,status='unknown')
      write(8,'(''SCI: '',i4)') it
      write(9,'(''SCI: '',i4)') it
c
      write(8,'(''## Spin 1'')')
      write(9,'(''## Spin 1'')')
      do li=1,nintfc
        cpalay=1.d0-conc(li).gt.tiny
        write(8,'(''##'',i3,t10,''A'')') li
        write(9,'(''##'',i3,t10,''A'')') li
        do ie=1,ne
          er=dreal(cear(ie))
          sum=0.d0
          do l=0,lmax
            sum=sum+dosa(l,li,ie,1)
          end do
          write(8,'(f10.6,5f10.4)') er,(dosa(l,li,ie,1),
     >                                 l=0,lmax),sum
        end do
        do ie=1,ne
          er=dreal(cear(ie))
          sum=0.d0
          do i=1,lmmax
            sum=sum+doslma(i,li,ie,1)
          end do
          write(9,'(f10.6,17f10.4)') er,(doslma(i,li,ie,1),
     >                                 i=1,lmmax),sum
        end do
        if(cpalay) then
        write(8,'(''##'',t10,''B'')')
        write(9,'(''##'',t10,''B'')')
        do ie=1,ne
          er=dreal(cear(ie))
          sum=0.d0
          do l=0,lmax
            sum=sum+dosb(l,li,ie,1)
          end do
          write(8,'(f10.6,5f10.4)') er,(dosb(l,li,ie,1),
     >                                 l=0,lmax),sum
        end do
        do ie=1,ne
          er=dreal(cear(ie))
          sum=0.d0
          do i=1,lmmax
            sum=sum+doslmb(i,li,ie,1)
          end do
          write(9,'(f10.6,17f10.4)') er,(doslmb(i,li,ie,1),
     >                                 i=1,lmmax),sum
        end do
        end if
      end do
      write(8,'(''## Spin 2'')')
      write(9,'(''## Spin 2'')')
      do li=1,nintfc
        cpalay=1.d0-conc(li).gt.tiny
        write(8,'(''##'',i3,t10,''A'')') li
        write(9,'(''##'',i3,t10,''A'')') li
        do ie=1,ne
          er=dreal(cear(ie))
          sum=0.d0
          do l=0,lmax
            sum=sum+dosa(l,li,ie,2)
          end do
          write(8,'(f10.6,5f10.4)') er,(dosa(l,li,ie,2),
     >                                 l=0,lmax),sum
        end do
        do ie=1,ne
          er=dreal(cear(ie))
          sum=0.d0
          do i=1,lmmax
            sum=sum+doslma(i,li,ie,2)
          end do
          write(9,'(f10.6,17f10.4)') er,(doslma(i,li,ie,2),
     >                                 i=1,lmmax),sum
        end do
        if(cpalay) then
        write(8,'(''##'',t10,''B'')')
        write(9,'(''##'',t10,''B'')')
        do ie=1,ne
          er=dreal(cear(ie))
          sum=0.d0
          do l=0,lmax
            sum=sum+dosb(l,li,ie,2)
          end do
          write(8,'(f10.6,5f10.4)') er,(dosb(l,li,ie,2),
     >                                 l=0,lmax),sum
        end do
        do ie=1,ne
          er=dreal(cear(ie))
          sum=0.d0
          do i=1,lmmax
            sum=sum+doslmb(i,li,ie,2)
          end do
          write(9,'(f10.6,17f10.4)') er,(doslmb(i,li,ie,2),
     >                                 i=1,lmmax),sum
        end do
        end if
      end do
c
c
c - Print No. of electrons in the layers
c
      if(lmax.eq.1) then
      write(8,'(/2x,''No. of electrons:''/t15,
     >          ''         s     '',
     >          ''         p     '',
     >          ''       total   '')')
      end if
      if(lmax.eq.2) then
      write(8,'(/2x,''No. of electrons:''/t15,
     >          ''         s     '',
     >          ''         p     '',
     >          ''         d     '',
     >          ''       total   '')')
      end if
      if(lmax.eq.3) then
      write(8,'(/2x,''No. of electrons:''/t15,
     >          ''         s     '',
     >          ''         p     '',
     >          ''         d     '',
     >          ''         f     '',
     >          ''       total   '')')
      end if
      do li=1,nintfc
        cpalay=1.d0-conc(li).gt.tiny
        write(8,'(2x,''L'',i3,''  A'',t15,5f15.10)') 
     >  li,(qvpa(l,li),l=0,lmax),qva(li)
        if(cpalay) then
        write(8,'(6x,''  B'',t15,5f15.10)') 
     >  (qvpb(l,li),l=0,lmax),qvb(li)
        end if
      end do
c
c - Print magnetic moments in the layers
c
      if(lmax.eq.1) then
      write(8,'(/2x,''Magnetic moments:''/t15,
     >          ''         s     '',
     >          ''         p     '',
     >          ''       total   '')')
      end if
      if(lmax.eq.2) then
      write(8,'(/2x,''Magnetic moments:''/t15,
     >          ''         s     '',
     >          ''         p     '',
     >          ''         d     '',
     >          ''       total   '')')
      end if
      if(lmax.eq.3) then
      write(8,'(/2x,''Magnetic moments:''/t15,
     >          ''         s     '',
     >          ''         p     '',
     >          ''         d     '',
     >          ''         f     '',
     >          ''       total   '')')
      end if
      do li=1,nintfc
        cpalay=1.d0-conc(li).gt.tiny
        write(8,'(2x,''L'',i3,''  A'',t15,5f15.10)') 
     >  li,(qmagvpa(l,li),l=0,lmax),qmagva(li)
        if(cpalay) then
        write(8,'(6x,''  B'',t15,5f15.10)') 
     >  (qmagvpb(l,li),l=0,lmax),qmagvb(li)
        end if
      end do
c
c - Print band energies
c
      write(8,'(/2x,''Band energies:''/)')
      etb=0.
      do li=1,nintfc
        cpalay=1.d0-conc(li).gt.tiny
        write(8,'(2x,''L'',i3,''  A'',t15,f15.10)') li,enba(li)
        if(cpalay) then
        write(8,'(6x,''  B'',t15,f15.10)') enbb(li)
        end if
        etb=etb+conc(li)*enba(li)+(1.0-conc(li))*enbb(li)
      end do
      write(8,'(/2x,''Total'',t15,f15.10)') etb
c
c - Print moments
c
      write(8,'(/2x,''Moments of charge densities:''/)')
      i=0
      do l=0,2*lmax 
      do m=-l,l     
         write(8,'(/10x,''(l,m)'',2i3/)') l,m
         i=i+1
         do li=1,nintfc
           cpalay=1.d0-conc(li).gt.tiny
           if(cpalay) then
             write(8,'(2x,'' L'',i2,t10,2d15.7,5x,2d15.7)') li,
     >       qmoma(i,li),qmomb(i,li)
           else
             write(8,'(2x,'' L'',i2,t10,2d15.7)') li,qmoma(i,li)
           end if
         end do
      end do
      end do
c
c - Print total energy
c
      write(8,'(/2x,''Total energy:''/
     >          ''Layer'',
     >          ''     Kinetic'',
     >          ''          El.stat.'',
     >          ''        Exc.-corr.'',
     >          ''       Total'')')
      etk=0.d0
      etu=0.d0
      etx=0.d0
      do li=1,nintfc
        cpalay=1.d0-conc(li).gt.tiny
        etot=enkina(li)+enela(li)+enxca(li)
        write(8,'(i2,'' A '',f15.7,3f17.7)') li,
     >  enkina(li),enela(li),enxca(li),etot
        etk=etk+conc(li)*enkina(li)
        etu=etu+conc(li)*enela(li)
        etx=etx+conc(li)*enxca(li)
        if(cpalay) then
        etot=enkinb(li)+enelb(li)+enxcb(li)
        write(8,'(2x,'' B '',f15.7,3f17.7)') 
     >  enkinb(li),enelb(li),enxcb(li),etot
        etk=etk+(1.d0-conc(li))*enkinb(li)
        etu=etu+(1.d0-conc(li))*enelb(li)
        etx=etx+(1.d0-conc(li))*enxcb(li)
        end if
      end do
      write(8,'(/'' Ifc '',f15.7,3f17.7)') etk,etu,etx,etk+etu+etx
      write(8,'(/'' Ifc '',f15.7,3f17.7)')
     >etk,etu+emad,etx,etk+etu+etx+emad
c
      close(8)
      close(9)
c
c - Print potential in file potout
c
      open(7,file=potout,status='unknown')
      write(7,'(i6,d20.8,4x,''.false.'',d20.8,4x, ''SCI: '',i4)')
     >ns(1),dx(1),efermi,it
c
      do li=1,nintfc
         cpalay=1.d0-conc(li).gt.tiny
         cc=conc(li)
         if(dlm) cc=1.0d0
         write(7,'(a)') idpota(li)
         write(7,'(f8.5)') cc
         write(7,'(f8.3,f15.10,d20.10,i5,f18.12)')
     >   za(li),vrsbulk,dx(li),ns(li),rs(li)
         write(7,'(4d20.12)') (vra(irad,li),irad=1,ns(li))
         write(7,'(4d20.12)') (bra(irad,li),irad=1,ns(li))
         if(cpalay.and.(.not.dlm)) then
         write(7,'(a)') idpotb(li)
         write(7,'(f8.5)') 1.d0-conc(li)
         write(7,'(f8.3,f15.10,d20.10,i5,f18.12)')
     >   zb(li),vrsbulk,dx(li),ns(li),rs(li)
         write(7,'(4d20.12)') (vrb(irad,li),irad=1,ns(li))
         write(7,'(4d20.12)') (brb(irad,li),irad=1,ns(li))
         end if
      end do
c
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
c
      close(7)
c
      return
      end
