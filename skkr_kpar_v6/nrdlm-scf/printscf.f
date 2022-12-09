c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine printscf(
     & nintfc,itscfcont,conc,qva,qvb,qmagva,qmagvb)
c
      implicit real*8 (a-h,o-z)
      include '../param.h' 
c
      logical cpalay
      dimension conc(mintfc),qva(mintfc),qvb(mintfc)
      dimension qmagva(mintfc),qmagvb(mintfc)
      data tiny/1.0d-6/ 
c
c -----------------------------------------------------------------------
c -- Write layer charges ------------------------------------------------
c -----------------------------------------------------------------------
      write(6,'(/"    Layer charges")')
      do li=1,nintfc
         if(1.d0-conc(li).lt.tiny) then
           write(6,'('' I'',i4,''  L'',i3,''   Q   '',f18.10)')
     >     itscfcont,li,qva(li)
         else
           write(6,'('' I'',i4,''  L'',i3,''   Q_A '',f18.10,
     >     ''   Q_B '',d18.10)') itscfcont,li,qva(li),qvb(li)
         endif
      enddo
      call flush(6)
c
c -----------------------------------------------------------------------
c -- Write spin moments -------------------------------------------------
c -----------------------------------------------------------------------
      write(6,'(/"    Magnetic moments")')
      do li=1,nintfc
        if(1.d0-conc(li).lt.tiny) then
          write(6,'('' I'',i4,''  L'',i3,''   M   '',f18.10)')
     >    itscfcont,li,qmagva(li)
        else
          write(6,'('' I'',i4,''  L'',i3,''   M_A '',f18.10,
     >    ''   M_B '',d18.10)') itscfcont,li,qmagva(li),qmagvb(li)
        endif
      enddo
      call flush(6)
c
      return
      end
