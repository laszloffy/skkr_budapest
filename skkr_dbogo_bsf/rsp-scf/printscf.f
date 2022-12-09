c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine printscf(
     > linbw,nintfc,itscfcont,conc,enbifc,qvifc,omifc,
     > qva,qvha,qveha,qvhea,qvteha,qvthea,
     > qvb,qvhb,qvehb,qvheb,qvtehb,qvtheb,
     > qca,za,qcb,zb,enba,enbb,efermi,
     > ddphenba,ddthenba,d2dphenba,d2dthenba,d2dthphenba,
     > ddphenbb,ddthenbb,d2dphenbb,d2dthenbb,d2dthphenbb,
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
      dimension qvha(mintfc),qvhb(mintfc)
      dimension qca(mintfc),qcb(mintfc)
      dimension za(mintfc),zb(mintfc)
      dimension enba(mintfc),enbb(mintfc)
      dimension ddphenba(mintfc),ddphenbb(mintfc)
      dimension ddthenba(mintfc),ddthenbb(mintfc)
      dimension d2dphenba(mintfc),d2dphenbb(mintfc)
      dimension d2dthenba(mintfc),d2dthenbb(mintfc)
      dimension d2dthphenba(mintfc),d2dthphenbb(mintfc)
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
      complex*16 qveha(mintfc),qvhea(mintfc)
      complex*16 qvehb(mintfc),qvheb(mintfc)
      complex*16 qvteha(mintfc),qvthea(mintfc)
      complex*16 qvtehb(mintfc),qvtheb(mintfc)
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
      qvifc=0.d0
      enbifc=0.d0
      omifc=0.d0
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
      do li=1,nintfc
         if(1.d0-conc(li).lt.tiny) then
          write(6,'('' I'',i4,''  L'',i4,''  electron,  Q  '',d18.10)')
     >     itscfcont,li,qva(li)
         else
           write(6,'('' I'',i4,''  L'',i4,'' electron,  Q_A '',d18.10,
     >     ''   Q_B '',d18.10)') itscfcont,li,qva(li),qvb(li)
         end if
      end do
      call flush(6)
c
      write(6,*)
      do li=1,nintfc
         if(1.d0-conc(li).lt.tiny) then
          write(6,'('' I'',i4,''  L'',i4,''  hole,   Q  '',d18.10)')
     >     itscfcont,li,qvha(li)
         else
           write(6,'('' I'',i4,''  L'',i4,'' hole,  Q_A '',d18.10,
     >     ''   Q_B '',d18.10)') itscfcont,li,qvha(li),qvhb(li)
         end if

      end do
      call flush(6)
c
      write(6,*)
      do li=1,nintfc
         if(1.d0-conc(li).lt.tiny) then
           write(6,'('' I'',i4,''  L'',i4,'' singlet,  Q  '',d18.10)')
     >      itscfcont,li,dreal(qveha(li))
         else
           write(6,'('' I'',i4,''  L'',i4,'' singlet,  Q_A '',d18.10,
     >     ''   Q_B '',d18.10)') itscfcont,li,
     >       dreal(qveha(li)),dreal(qvehb(li))
c           write(6,*) 'I ',itscfcont, 'L singlet',li,qveha(li),qvehb(li)
         end if

      end do
      call flush(6)
c
      write(6,*)
      do li=1,nintfc
         if(1.d0-conc(li).lt.tiny) then
           write(6,'('' I'',i4,''  L'',i4,'' T 0, Q  '',d18.10)')
     >      itscfcont,li,dreal(qvhea(li))
         else
           write(6,'('' I'',i4,''  L'',i4,'' T 0,  Q_A '',d18.10,
     >     ''   Q_B '',d18.10)') itscfcont,li,
     >       dreal(qvhea(li)),dreal(qvheb(li))
c           write(6,*) 'I ',itscfcont, 'L triplet 0 ',
c     >      li,qvhea(li),qvheb(li)
         end if

      end do
      call flush(6)
c
      write(6,*)
      do li=1,nintfc
         if(1.d0-conc(li).lt.tiny) then
           write(6,'('' I'',i4,''  L'',i4,'' T 1,  Q  '',d18.10)')
     >     itscfcont,li,dreal(qvteha(li))
         else
           write(6,'('' I'',i4,''  L'',i4,'' T 1,  Q_A '',d18.10,
     >     ''   Q_B '',d18.10)') itscfcont,li,
     >       dreal(qvteha(li)),dreal(qvtehb(li))
c           write(6,*) 'I ',itscfcont, 'L  triplet up',
c     >       li,qvteha(li),qvtehb(li)
         end if

      end do
      call flush(6)
c
      write(6,*)
      do li=1,nintfc
         if(1.d0-conc(li).lt.tiny) then
           write(6,'('' I'',i4,''  L'',i4,'' T -1,  Q '',d18.10)')
     >     itscfcont,li,dreal(qvthea(li))
         else
           write(6,'('' I'',i4,''  L'',i4,'' T -1,  Q_A '',d18.10,
     >     ''   Q_B '',d18.10)') itscfcont,li, 
     >        dreal(qvthea(li)),dreal(qvtheb(li))
c           write(6,*) 'I ',itscfcont, 'L  triplet up',
c     >     li,qvthea(li),qvtheb(li)
         end if


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
      write(6,*)
      do li=1,nintfc
         cpalay=1.d0-conc(li).gt.tiny
         if(.not.cpalay) then
          write(6,'('' I'',i4,''  L'',i4,''  MS (electron)  '',d18.10)')
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
          write(6,'('' I'',i4,''  L'',i4,''  MO (electron)  '',d18.10)')
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
