c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine printscf1(
     > nimp,itscfcont,
     > qvimp,qvhimp,qvehimp,qvheimp,qvtehimp,qvtheimp)
c
      implicit real*8 (a-h,o-z)
      include '../param.h' 
c
c      logical cpalay,linbw
c      dimension rba(3,mintfc),rbb(3,mintfc)
c      dimension conc(mintfc),qva(mintfc),qvb(mintfc)
c      dimension qvha(mintfc),qvhb(mintfc)
c      dimension qca(mintfc),qcb(mintfc)
c      dimension za(mintfc),zb(mintfc)
c      dimension enba(mintfc),enbb(mintfc)
c      dimension ddphenba(mintfc),ddphenbb(mintfc)
c      dimension ddthenba(mintfc),ddthenbb(mintfc)
c      dimension d2dphenba(mintfc),d2dphenbb(mintfc)
c      dimension d2dthenba(mintfc),d2dthenbb(mintfc)
c      dimension d2dthphenba(mintfc),d2dthphenbb(mintfc)
c      dimension spin_magva(mintfc,3),spin_magvb(mintfc,3)
c      dimension orb_magva(mintfc,3),orb_magvb(mintfc,3)    
c      dimension veca(mintfc,3),vecb(mintfc,3),vec(2)    
c
c      dimension spinmoma(mintfc),orbmoma(mintfc)
c      dimension th0a(mintfc),th1a(mintfc)
c      dimension ph0a(mintfc),ph1a(mintfc)
c      dimension spinmomb(mintfc),orbmomb(mintfc)
c      dimension th0b(mintfc),th1b(mintfc)
c      dimension ph0b(mintfc),ph1b(mintfc)
c
      real*8 qvimp(nimp),qvhimp(nimp)
      real*8 qvehimp(nimp),qvheimp(nimp)
      real*8 qvtehimp(nimp),qvtheimp(nimp)
c
      data tiny/1.0d-6/ 
      pi=dacos(-1.0d0)
      fac=180.d0/pi
c
c      if(linbw) then
c        write(6,'(/''   '',t13,''Eband'',t33,''  Q'',t53,'' Omega'')')
c        write(6,'(3x,3d20.12)') enbifc,qvifc,omifc
c      end if
c
c        write(6,'(/''   '',t13,''Eb A'',t33,''  Eb B  '',t53,''Eb'')')
c
c      qvifc=0.d0
c      enbifc=0.d0
c      omifc=0.d0
c      do li=1,nintfc
c          q=conc(li)*qva(li)+(1.0d0-conc(li))*qvb(li)
c          enb=conc(li)*enba(li)+(1.0d0-conc(li))*enbb(li)
c          om=enb-efermi*q
c         om=enb
c          qvifc=qvifc+q
c          enbifc=enbifc+enb
c          omifc=omifc+om
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
c         write(6,'(i3,3d20.12)') li,enba(li),enbb(li),enb
c          write(6,'(i3,3d20.12)') li,
c     >       enba(li)-efermi*qva(li),
c     >       enbb(li)-efermi*qvb(li),om
c      end do
c
c      write(6,'(''-----------------------------------------------'',
c     >''----------------'')')
c      write(6,'(t44,d20.12)') omifc
c      call flush(6)
c
c      write(6,'(/t13,''DEb A'',t33,'' DEb B  '')')
c      write(6,'(''-------------------------------------------'')')
c      write(6,'(t13,''th  A'',t33,'' th  B  '')')
c      dda=0.0d0
c      ddb=0.0d0
c      do li=1,nintfc
c          write(6,'(i3,3d20.12)') li,ddthenba(li),ddthenbb(li)
c          dda=dda+ddthenba(li)
c          ddb=ddb+ddthenbb(li)
c      end do
c      write(6,'(''-------------------------------------------'')')
c      write(6,'(3x,3d20.12)') dda,ddb
c      write(6,'(''-------------------------------------------'')')
c      write(6,'(t13,''ph  A'',t33,'' ph  B  '')')
c      do li=1,nintfc
c          write(6,'(i3,3d20.12)') li,ddphenba(li),ddphenbb(li)
c      end do
c      write(6,'(''-------------------------------------------'')')
c      write(6,'(t13,''ththA'',t33,'' ththB  '')')
c      do li=1,nintfc
c          write(6,'(i3,3d20.12)') li,d2dthenba(li),d2dthenbb(li)
c      end do
c      write(6,'(''-------------------------------------------'')')
c      write(6,'(t13,''phphA'',t33,'' phphB  '')')
c      do li=1,nintfc
c          write(6,'(i3,3d20.12)') li,d2dphenba(li),d2dphenbb(li)
c      end do
c      write(6,'(''-------------------------------------------'')')
c      write(6,'(t13,''thphA'',t33,'' thphB  '')')
c      do li=1,nintfc
c          write(6,'(i3,3d20.12)') li,d2dthphenba(li),d2dthphenbb(li)
c      end do
c      write(6,'(''-------------------------------------------'')')
c
c      if(linbw) return
c
c------------------------------
c++++ Write charges to prn ++++
c------------------------------
c
      write(6,*)
      do li=1,nimp
        write(6,'('' I'',i4,''  imp'',i4,''  electron,  Q  '',d18.10)')
     >     itscfcont,li,qvimp(li)
      end do
      call flush(6)
c
      write(6,*)
      do li=1,nimp
          write(6,'('' I'',i4,''  imp'',i4,''  hole,   Q  '',d18.10)')
     >     itscfcont,li,qvhimp(li)
      end do
c
      write(6,*)
      do li=1,nimp
           write(6,'('' I'',i4,''  imp'',i4,'' S,  Q  '',d18.10)')
     >      itscfcont,li,qvehimp(li)
      end do
c
      write(6,*)
      do li=1,nimp
           write(6,'('' I'',i4,''  imp'',i4,'' T 0, Q  '',d18.10)')
     >      itscfcont,li,qvheimp(li)
      end do
c
      write(6,*)
      do li=1,nimp
           write(6,'('' I'',i4,''  imp'',i4,'' T 1,  Q  '',d18.10)')
     >     itscfcont,li,qvtehimp(li)
      end do
c
      write(6,*)
      do li=1,nimp
           write(6,'('' I'',i4,''  imp'',i4,'' T -1,  Q '',d18.10)')
     >     itscfcont,li,qvtheimp(li)
      end do
      call flush(6)
c
c      write(6,*)
c      do li=1,nintfc
c         cpalay=1.d0-conc(li).gt.tiny
c         if(.not.cpalay) then
c           write(6,'('' I'',i4,''  L'',i4,''  Bth   '',2f15.5)')
c     >     itscfcont,li,th0a(li)*fac,th1a(li)*fac
c         else
c           write(6,'('' I'',i4,''  L'',i4,''  Bth A '',2f15.5,
c     >     '' B '',2f15.5)') itscfcont,li,
c     >     th0a(li)*fac,th1a(li)*fac,
c     >     th0b(li)*fac,th1b(li)*fac
c         end if
c      end do
c      call flush(6)
c
c      write(6,*)
c      do li=1,nintfc
c         cpalay=1.d0-conc(li).gt.tiny
c         if(.not.cpalay) then
c           write(6,'('' I'',i4,''  L'',i4,''  Bph   '',2f15.5)')
c     >     itscfcont,li,ph0a(li)*fac,ph1a(li)*fac
c         else
c           write(6,'('' I'',i4,''  L'',i4,''  Bph A '',2f15.5,
c     >     '' B '',2f15.5)') itscfcont,li,
c     >     ph0a(li)*fac,ph1a(li)*fac,
c     >     ph0b(li)*fac,ph1b(li)*fac
c         end if
c      end do
c      call flush(6)
c
c     write(6,*)
c     do li=1,nintfc
c        cpalay=1.d0-conc(li).gt.tiny
c        if(.not.cpalay) then
c         write(6,'('' I'',i4,''  L'',i4,''  MS (electron)  '',d18.10)')
c    >     itscfcont,li,spinmoma(li)
c        else
c          write(6,'('' I'',i4,''  L'',i4,''  MS_A'',d18.10,
c    >     ''  MS_B'',d18.10)') itscfcont,li,spinmoma(li),spinmomb(li)
c        end if
c     end do
c     call flush(6)
c
c     write(6,*)
c     do li=1,nintfc
c        cpalay=1.d0-conc(li).gt.tiny
c        if(.not.cpalay) then
c         write(6,'('' I'',i4,''  L'',i4,''  MO (electron)  '',d18.10)')
c    >     itscfcont,li,orbmoma(li)
c        else
c          write(6,'('' I'',i4,''  L'',i4,''  MO_A'',d18.10,
c    >     ''  MO_B'',d18.10)') itscfcont,li,orbmoma(li),orbmomb(li)
c        end if
c     end do
c     call flush(6)
c
      return
      end
