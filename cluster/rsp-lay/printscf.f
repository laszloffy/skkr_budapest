c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine printscf(
     > linbw,nintfc,itscfcont,conc,enbifc,qvifc,omifc,
     > qva,qca,za,qvb,qcb,zb,enba,enbb,denba,denbb,efermi,
     > arot,rba,rbb,spin_magva,spin_magvb,orb_magva,orb_magvb,
     > spinmoma,orbmoma,th0a,th1a,ph0a,ph1a,
     > spinmomb,orbmomb,th0b,th1b,ph0b,ph1b)
c
      implicit real*8 (a-h,o-z)
      include '../param.h' 
c
      logical cpalay,linbw
      dimension rba(3,mintfc),rbb(3,mintfc)
      dimension conc(mintfc),qva(mintfc),qvb(mintfc)
      dimension qca(mintfc),qcb(mintfc)
      dimension za(mintfc),zb(mintfc)
      dimension enba(mintfc),enbb(mintfc)
      dimension denba(mintfc),denbb(mintfc)
      dimension spin_magva(mintfc,3),spin_magvb(mintfc,3)
      dimension orb_magva(mintfc,3),orb_magvb(mintfc,3)    
      dimension veca(mintfc,3),vecb(mintfc,3),vec(2)    
c
      dimension spinmoma(mintfc),orbmoma(mintfc)
      dimension th0a(mintfc),th1a(mintfc)
      dimension ph0a(mintfc),ph1a(mintfc)
      dimension spinmomb(mintfc),orbmomb(mintfc)
      dimension th0b(mintfc),th1b(mintfc)
      dimension ph0b(mintfc),ph1b(mintfc)
c
      data tiny/1.0d-6/ 
      pi=dacos(-1.0d0)
      fac=180.d0/pi
c
      if(linbw) then
        write(6,'(/''   '',t13,''Eband'',t33,''  Q'',t53,'' Omega'')')
        write(6,'(3x,3d20.12)') enbifc,qvifc,omifc
      end if
c
        write(6,'(/''   '',t13,''Eb A'',t33,''  Eb B  '',t53,''Eb'')')
c
      qvifc=0.d0
      enbifc=0.d0
      omifc=0.d0
      do li=1,nintfc
          q=conc(li)*qva(li)+(1.0d0-conc(li))*qvb(li)
          enb=conc(li)*enba(li)+(1.0d0-conc(li))*enbb(li)
          om=enb-efermi*q
c         om=enb
          qvifc=qvifc+q
          enbifc=enbifc+enb
          omifc=omifc+om
c         if(li.eq.4) then
c            cth=rba(3,li)
c            sth=dsqrt(1.0d0-cth*cth)
c            th=dacos(cth)
c            cph=rba(1,li)/sth
c            sph=rba(2,li)/sth
c            ph=dacos(cph)
c            if(sph.lt.0.d0) ph=2*pi-ph
c            write(6,'(i3,2d15.6,d20.12)') li,th,ph,
c    >       (enba(li)-efermi*qva(li))*rba(3,li)
c    >                               enbb(li)-efermi*qvb(li),om
c         end if
          write(6,'(i3,3d20.12)') li,
     >       enba(li)-efermi*qva(li),
     >       enbb(li)-efermi*qvb(li),om
c         write(6,'(i3,3d20.12)') li,enba(li),enbb(li),om
      end do
c
      write(6,'(''-----------------------------------------------'',
     >''----------------'')')
      write(6,'(t44,d20.12)') omifc
      call flush(6)
c
      write(6,'(/''   '',t13,''DEb A'',t33,'' DEb B  '')')
      do li=1,nintfc
          write(6,'(i3,2d20.12)') li,denba(li),denbb(li)
      end do
c
      if(linbw) return
c
      write(6,*)
      do li=1,nintfc
         if(1.d0-conc(li).lt.tiny) then
           write(6,'('' I'',i4,''  L'',i4,''   Q   '',d18.10)')
     >     itscfcont,li,qva(li)
         else
           write(6,'('' I'',i4,''  L'',i4,''   Q_A '',d18.10,
     >     ''   Q_B '',d18.10)') itscfcont,li,qva(li),qvb(li)
         end if
      end do
      call flush(6)
c
      write(6,*)
      do li=1,nintfc
         cpalay=1.d0-conc(li).gt.tiny
         if(.not.cpalay) then
           write(6,'('' I'',i4,''  L'',i4,''  Bth   '',2f15.5)')
     >     itscfcont,li,th0a(li)*fac,th1a(li)*fac
         else
           write(6,'('' I'',i4,''  L'',i4,''  Bth A '',2f15.5,
     >     '' B '',2f15.5)') itscfcont,li,
     >     th0a(li)*fac,th1a(li)*fac,
     >     th0b(li)*fac,th1b(li)*fac
         end if
      end do
      call flush(6)
c
      write(6,*)
      do li=1,nintfc
         cpalay=1.d0-conc(li).gt.tiny
         if(.not.cpalay) then
           write(6,'('' I'',i4,''  L'',i4,''  Bph   '',2f15.5)')
     >     itscfcont,li,ph0a(li)*fac,ph1a(li)*fac
         else
           write(6,'('' I'',i4,''  L'',i4,''  Bph A '',2f15.5,
     >     '' B '',2f15.5)') itscfcont,li,
     >     ph0a(li)*fac,ph1a(li)*fac,
     >     ph0b(li)*fac,ph1b(li)*fac
         end if
      end do
      call flush(6)
c
      write(6,*)
      do li=1,nintfc
         cpalay=1.d0-conc(li).gt.tiny
         if(.not.cpalay) then
           write(6,'('' I'',i4,''  L'',i4,''  MS  '',d18.10)')
     >     itscfcont,li,spinmoma(li)
         else
           write(6,'('' I'',i4,''  L'',i4,''  MS_A'',d18.10,
     >     ''  MS_B'',d18.10)') itscfcont,li,spinmoma(li),spinmomb(li)
         end if
      end do
      call flush(6)
c
      write(6,*)
      do li=1,nintfc
         cpalay=1.d0-conc(li).gt.tiny
         if(.not.cpalay) then
           write(6,'('' I'',i4,''  L'',i4,''  MO  '',d18.10)')
     >     itscfcont,li,orbmoma(li)
         else
           write(6,'('' I'',i4,''  L'',i4,''  MO_A'',d18.10,
     >     ''  MO_B'',d18.10)') itscfcont,li,orbmoma(li),orbmomb(li)
         end if
      end do
      call flush(6)
c
      return
      end
