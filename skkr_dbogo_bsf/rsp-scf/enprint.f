c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine enprint(itscf,nintfc,conc,enca,encb,enba,enbb,
     >    enela,enelb,enxca,enxcb,enpota,enpotb,enmaga,enmagb,
     >    enorba,enorbb,entota,entotb,entot,entotifc)
c=======================
c
c -calculate and print out energies
c
      implicit real*8(a-h,o-z)
c
      include '../param.h'
c
      logical cpalay
c
      dimension enca(mintfc),encb(mintfc),enc(mintfc)
      dimension enba(mintfc),enbb(mintfc),enb(mintfc)
      dimension enela(mintfc),enelb(mintfc),enel(mintfc)
      dimension enxca(mintfc),enxcb(mintfc),enxc(mintfc)
      dimension enpota(mintfc),enpotb(mintfc),enpot(mintfc)
      dimension enmaga(mintfc),enmagb(mintfc),enmag(mintfc)
      dimension enorba(mintfc),enorbb(mintfc),enorb(mintfc)
      dimension entota(mintfc),entotb(mintfc),entot(mintfc)
      dimension conc(mintfc)
c
      common/test/itest
      common/madelung_energy/emad
c
      data tiny/1.0d-6/
c
      entotifc=0.d0
      enbifc=0.d0
      do li=1,nintfc
c
        enc(li)=conc(li)*enca(li)+(1.d0-conc(li))*encb(li)
        enb(li)=conc(li)*enba(li)+(1.d0-conc(li))*enbb(li)
        enel(li)=conc(li)*enela(li)+(1.d0-conc(li))*enelb(li)
        enxc(li)=conc(li)*enxca(li)+(1.d0-conc(li))*enxcb(li)
        enpot(li)=conc(li)*enpota(li)+(1.d0-conc(li))*enpotb(li)
        enmag(li)=conc(li)*enmaga(li)+(1.d0-conc(li))*enmagb(li)
        enorb(li)=conc(li)*enorba(li)+(1.d0-conc(li))*enorbb(li)
c
        entota(li)=enca(li)+enba(li)+enela(li)+enxca(li)+
     >             enpota(li)+enmaga(li)+enorba(li)
        entotb(li)=encb(li)+enbb(li)+enelb(li)+enxcb(li)+
     >             enpotb(li)+enmagb(li)+enorbb(li)
        entot(li)=conc(li)*entota(li)+(1.d0-conc(li))*entotb(li)
        enbifc=enbifc+enb(li)
        entotifc=entotifc+entot(li)
      end do
c
      if(itest.ge.2) then
      write(6,*)
c
      do li=1,nintfc
        cpalay=1.d0-conc(li).gt.tiny
        if(cpalay) then
        write(6,'('' I'',i4,'' L'',i4,''  Ecore'',t20,
     >  ''A '',d23.15,''  B '',d23.15,3x,d23.15)')
     >  itscf,li,enca(li),encb(li),enc(li)
        else
        write(6,'('' I'',i4,'' L'',i4,''  Ecore'',t22,d23.15)')
     >  itscf,li,enca(li)
        end if
      end do
      write(6,*)
c
      do li=1,nintfc
        cpalay=1.d0-conc(li).gt.tiny
        if(cpalay) then
        write(6,'('' I'',i4,'' L'',i4,''  Eband'',t20,
     >  ''A '',d23.15,''  B '',d23.15,3x,d23.15)')
     >  itscf,li,enba(li),enbb(li),enb(li)
        else
        write(6,'('' I'',i4,'' L'',i4,''  Eband'',t22,d23.15)')
     >  itscf,li,enba(li)
        end if
      end do
      write(6,*)
c
      do li=1,nintfc
        cpalay=1.d0-conc(li).gt.tiny
        if(cpalay) then
        write(6,'('' I'',i4,'' L'',i4,''  Eelst'',t20,
     >  ''A '',d23.15,''  B '',d23.15,3x,d23.15)')
     >  itscf,li,enela(li),enelb(li),enel(li)
        else
        write(6,'('' I'',i4,'' L'',i4,''  Eelst'',t22,d23.15)')
     >  itscf,li,enela(li)
        end if
      end do
      write(6,*)
c
      do li=1,nintfc
        cpalay=1.d0-conc(li).gt.tiny
        if(cpalay) then
        write(6,'('' I'',i4,'' L'',i4,''  Exc  '',t20,
     >  ''A '',d23.15,''  B '',d23.15,3x,d23.15)')
     >  itscf,li,enxca(li),enxcb(li),enxc(li)
        else
        write(6,'('' I'',i4,'' L'',i4,''  Exc  '',t22,d23.15)')
     >  itscf,li,enxca(li)
        end if
      end do
      write(6,*)
c
      do li=1,nintfc
        cpalay=1.d0-conc(li).gt.tiny
        if(cpalay) then
        write(6,'('' I'',i4,'' L'',i4,''  Epot '',t20,
     >  ''A '',d23.15,''  B '',d23.15,3x,d23.15)')
     >  itscf,li,enpota(li),enpotb(li),enpot(li)
        else
        write(6,'('' I'',i4,'' L'',i4,''  Epot '',t22,d23.15)')
     >  itscf,li,enpota(li)
        end if
      end do
      write(6,*)
c
      do li=1,nintfc
        cpalay=1.d0-conc(li).gt.tiny
        if(cpalay) then
        write(6,'('' I'',i4,'' L'',i4,''  Emag '',t20,
     >  ''A '',d23.15,''  B '',d23.15,3x,d23.15)')
     >  itscf,li,enmaga(li),enmagb(li),enmag(li)
        else
        write(6,'('' I'',i4,'' L'',i4,''  Emag '',t22,d23.15)')
     >  itscf,li,enmaga(li)
        end if
      end do
      write(6,*)
c
      do li=1,nintfc
        cpalay=1.d0-conc(li).gt.tiny
        if(cpalay) then
        write(6,'('' I'',i4,'' L'',i4,''  Eorb '',t20,
     >  ''A '',d23.15,''  B '',d23.15,3x,d23.15)')
     >  itscf,li,enorba(li),enorbb(li),enorb(li)
        else
        write(6,'('' I'',i4,'' L'',i4,''  Eorb '',t22,d23.15)')
     >  itscf,li,enorba(li)
        end if
      end do
      write(6,*)
c
      endif
c
      do li=1,nintfc
        cpalay=1.d0-conc(li).gt.tiny
        if(cpalay) then
        write(6,'('' I'',i4,'' L'',i4,''  Etot (just electron)'',t20,
     >  ''A '',d23.15,''  B '',d23.15,3x,d23.15)')
     >  itscf,li,entota(li),entotb(li),entot(li)
        else
        write(6,'('' I'',i4,'' L'',i4,
     >        ''  Etot (just electron)'',t22,d23.15)')
     >  itscf,li,entota(li)
        end if
      end do
c
      write(6,'(/'' I'',i4,''  Emad_ifc (just electron)'',t22,d23.15)')
     >itscf,emad
c
      entotifc=entotifc+emad
      write(6,'(/'' I'',i4,''  Etot_ifc (just electron)'',t22,d23.15)')
     >itscf,entotifc
c
      return
      end
