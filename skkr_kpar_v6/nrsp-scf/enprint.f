c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine enprint(itscf,nintfc,conc,enca,encb,enba,enbb,
     >    enela,enelb,enxca,enxcb,enpota,enpotb,enkina,enkinb,
     >    enmaga,enmagb)
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
      dimension enba(mintfc),enbb(mintfc)
      dimension enca(mintfc),encb(mintfc)
      dimension enela(mintfc),enelb(mintfc)
      dimension enxca(mintfc),enxcb(mintfc)
      dimension enpota(mintfc),enpotb(mintfc)
      dimension enkina(mintfc),enkinb(mintfc)
      dimension enmaga(mintfc),enmagb(mintfc)
      dimension conc(mintfc)
c
      common/test/itest
      common/madelung_energy/emad
c
      data tiny/1.0d-6/
c
      entotifc=0.d0
      do li=1,nintfc
        cpalay=1.d0-conc(li).gt.tiny
c
        enc=conc(li)*enca(li)+(1.d0-conc(li))*encb(li)
        enb=conc(li)*enba(li)+(1.d0-conc(li))*enbb(li)
        enpot=conc(li)*enpota(li)+(1.d0-conc(li))*enpotb(li)
        enmag=conc(li)*enmaga(li)+(1.d0-conc(li))*enmagb(li)
        enel=conc(li)*enela(li)+(1.d0-conc(li))*enelb(li)
        enxc=conc(li)*enxca(li)+(1.d0-conc(li))*enxcb(li)

        enkina(li)=enba(li)+enca(li)+enpota(li)+enmaga(li)
        enkinb(li)=enbb(li)+encb(li)+enpotb(li)+enmagb(li)

        enkin=conc(li)*enkina(li)+(1.d0-conc(li))*enkinb(li)

        entota=enkina(li)+enela(li)+enxca(li)
        entotb=enkinb(li)+enelb(li)+enxcb(li)
        entot=conc(li)*entota+(1.d0-conc(li))*entotb
        entotifc=entotifc+entot
c
        if(itest.ge.2) then
c
        if(cpalay) then
        write(6,'(/'' I'',i4,'' L'',i3,''  Ecore'',t20,
     >  ''A '',d20.12,''  B '',d20.12,3x,d20.10)')
     >  itscf,li,enca(li),encb(li),enc
        write(6,'('' I'',i4,'' L'',i3,''  Eband'',t20,
     >  ''A '',d20.12,''  B '',d20.12,3x,d20.10)')
     >  itscf,li,enba(li),enbb(li),enb
        write(6,'('' I'',i4,'' L'',i3,''  Epot '',t20,
     >  ''A '',d20.12,''  B '',d20.12,3x,d20.10)')
     >  itscf,li,enpota(li),enpotb(li),enpot
        write(6,'('' I'',i4,'' L'',i3,''  Emag '',t20,
     >  ''A '',d20.12,''  B '',d20.12,3x,d20.10)')
     >  itscf,li,enmaga(li),enmagb(li),enmag
        write(6,'('' I'',i4,'' L'',i3,''  Ekin '',t20,
     >  ''A '',d20.12,''  B '',d20.12,3x,d20.10)')
     >  itscf,li,enkina(li),enkinb(li),enkin
        write(6,'('' I'',i4,'' L'',i3,''  Eelst'',t20,
     >  ''A '',d20.12,''  B '',d20.12,3x,d20.10)')
     >  itscf,li,enela(li),enelb(li),enel
        write(6,'('' I'',i4,'' L'',i3,''  Exc  '',t20,
     >  ''A '',d20.12,''  B '',d20.12,3x,d20.10)')
     >  itscf,li,enxca(li),enxcb(li),enxc
        else
        write(6,'(/'' I'',i4,'' L'',i3,''  Ecore'',t22,d20.12)')
     >  itscf,li,enca(li)
        write(6,'('' I'',i4,'' L'',i3,''  Eband'',t22,d20.12)')
     >  itscf,li,enba(li)
        write(6,'('' I'',i4,'' L'',i3,''  Epot '',t22,d20.12)')
     >  itscf,li,enpota(li)
        write(6,'('' I'',i4,'' L'',i3,''  Emag '',t22,d20.12)')
     >  itscf,li,enmaga(li)
        write(6,'('' I'',i4,'' L'',i3,''  Ekin '',t22,d20.12)')
     >  itscf,li,enkina(li)
        write(6,'('' I'',i4,'' L'',i3,''  Eelst'',t22,d20.12)')
     >  itscf,li,enela(li)
        write(6,'('' I'',i4,'' L'',i3,''  Exc  '',t22,d20.12)')
     >  itscf,li,enxca(li)
        end if
c
        end if
c
        if(cpalay) then
        write(6,'('' I'',i4,'' L'',i3,''  Etot '',t20,
     >  ''A '',d20.12,''  B '',d20.12,3x,d20.10)')
     >  itscf,li,entota,entotb,entot
        else
        write(6,'('' I'',i4,'' L'',i3,''  Etot '',t22,d20.12)')
     >  itscf,li,entota
        end if
c
      end do
c
      write(6,'(/'' I'',i4,''  Emad_ifc '',t22,d20.12)')
     >itscf,emad
c
      entotifc=entotifc+emad
      write(6,'(/'' I'',i4,''  Etot_ifc '',t22,d20.12)')
     >itscf,entotifc
c
      return
      end
