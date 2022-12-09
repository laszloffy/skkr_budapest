c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine enprint1(itscf,nintfc,enca,enba,
     >    enela,enxca,enpota,enmaga,
     >    enorba,entota,entot,entotifc)
c=======================
c
c -calculate and print out energies
c
      implicit real*8(a-h,o-z)
c
      include '../param.h'
c
      integer mmax
      parameter (mmax=mimp)
c
      real*8 enca(mmax)
      real*8 enc(mmax)
      real*8 enba(mmax)
      real*8 enb(mmax)
      real*8 enela(mmax)
      real*8 enel(mmax)
      real*8 enxca(mmax)
      real*8 enxc(mmax)
      real*8 enpota(mmax)
      real*8 enpot(mmax)
      real*8 enmaga(mmax)
      real*8 enmag(mmax)
      real*8 enorba(mmax)
      real*8 enorb(mmax)
      real*8 entota(mmax)
      real*8 entot(mmax)
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
        enc(li)=enca(li)
        enb(li)=enba(li)
        enel(li)=enela(li)
        enxc(li)=enxca(li)
        enpot(li)=enpota(li)
        enmag(li)=enmaga(li)
        enorb(li)=enorba(li)
c
        entota(li)=enca(li)+enba(li)+enela(li)+enxca(li)+
     >             enpota(li)+enmaga(li)+enorba(li)
        entot(li)=entota(li)
        enbifc=enbifc+enb(li)
        entotifc=entotifc+entot(li)
      end do
c
      if(itest.ge.2) then
      write(6,*)
c
      do li=1,nintfc
        write(6,'('' I'',i4,'' L'',i4,''  Ecore'',t22,d23.15)')
     >  itscf,li,enca(li)
      end do
      write(6,*)
c
      do li=1,nintfc
        write(6,'('' I'',i4,'' L'',i4,''  Eband'',t22,d23.15)')
     >  itscf,li,enba(li)
      end do
      write(6,*)
c
      do li=1,nintfc
        write(6,'('' I'',i4,'' L'',i4,''  Eelst'',t22,d23.15)')
     >  itscf,li,enela(li)
      end do
      write(6,*)
c
      do li=1,nintfc
        write(6,'('' I'',i4,'' L'',i4,''  Exc  '',t22,d23.15)')
     >  itscf,li,enxca(li)
      end do
      write(6,*)
c
      do li=1,nintfc
        write(6,'('' I'',i4,'' L'',i4,''  Epot '',t22,d23.15)')
     >  itscf,li,enpota(li)
      end do
      write(6,*)
c
      do li=1,nintfc
        write(6,'('' I'',i4,'' L'',i4,''  Emag '',t22,d23.15)')
     >  itscf,li,enmaga(li)
      end do
      write(6,*)
c
      do li=1,nintfc
        write(6,'('' I'',i4,'' L'',i4,''  Eorb '',t22,d23.15)')
     >  itscf,li,enorba(li)
      end do
      write(6,*)
c
      endif
c
      do li=1,nintfc
        write(6,'('' I'',i4,'' L'',i4,''  Etot '',t22,d23.15)')
     >  itscf,li,entota(li)
      end do
c
      write(6,'(/'' I'',i4,''  Emad_ifc '',t22,d23.15)')
     >itscf,emad
c
      entotifc=entotifc+emad
      write(6,'(/'' I'',i4,''  Etot_ifc '',t22,d23.15)')
     >itscf,entotifc
c
      return
      end
