c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine printscf(
     > linbw,nintfc,itscfcont,conc,enbifc,qvifc,omifc,
     > qva,qca,za,qvb,qcb,zb,enba,enbb,efermi,
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
      data tiny/1.0d-6/ 
      pi=dacos(-1.0d0)
      fac=180.d0/pi
c
      if(linbw) then
        write(6,'(/''   '',t13,''Eband'',t33,''  Q'',t53,'' Omega'')')
        write(6,'(3x,3d20.12)') enbifc,qvifc,omifc
      end if
c
c -----------------------------------------------------------------------
c -- Write grand canonical potentials -----------------------------------
c -----------------------------------------------------------------------
      write(6,'(/''   '',t13,''Eband A'',t30,''  
     >   Eband B  '',t53,''Eband'')')
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
          write(6,'(i3,3f20.12)') li,
     >       enba(li)-efermi*qva(li),
     >       enbb(li)-efermi*qvb(li),om
      end do
c
      write(6,'(''-----------------------------------------------'',
     >''----------------'')')
      write(6,'(t30,"OmegaSum:",f20.12)') omifc
      call flush(6)
c
      write(6,'(/t13,''DEb A'',t33,'' DEb B  '')')
      write(6,'(''-------------------------------------------'')')
      write(6,'(t13,''th  A'',t33,'' th  B  '')')
      dda=0.0d0
      ddb=0.0d0
      do li=1,nintfc
          write(6,'(i3,3d20.12)') li,ddthenba(li),ddthenbb(li)
          dda=dda+ddthenba(li)
          ddb=ddb+ddthenbb(li)
      end do
      write(6,'(''-------------------------------------------'')')
      write(6,'(3x,3f20.12)') dda,ddb
      write(6,'(''-------------------------------------------'')')
      write(6,'(t13,''ph  A'',t33,'' ph  B  '')')
      do li=1,nintfc
          write(6,'(i3,3f20.12)') li,ddphenba(li),ddphenbb(li)
      end do
      write(6,'(''-------------------------------------------'')')
      write(6,'(t13,''ththA'',t33,'' ththB  '')')
      do li=1,nintfc
          write(6,'(i3,3f20.12)') li,d2dthenba(li),d2dthenbb(li)
      end do
      write(6,'(''-------------------------------------------'')')
      write(6,'(t13,''phphA'',t33,'' phphB  '')')
      do li=1,nintfc
          write(6,'(i3,3f20.12)') li,d2dphenba(li),d2dphenbb(li)
      end do
      write(6,'(''-------------------------------------------'')')
      write(6,'(t13,''thphA'',t33,'' thphB  '')')
      do li=1,nintfc
          write(6,'(i3,3f20.12)') li,d2dthphenba(li),d2dthphenbb(li)
      end do
      write(6,'(''-------------------------------------------'')')
c
      if(linbw) return
c
c -----------------------------------------------------------------------
c -- Write charges ------------------------------------------------------
c -----------------------------------------------------------------------
      write(6,'(/"    Layer charges")')
      do li=1,nintfc
         if(1.d0-conc(li).lt.tiny) then
           write(6,'('' I'',i4,''  L'',i4,''   Q   '',f18.10)')
     >     itscfcont,li,qva(li)
         else
           write(6,'('' I'',i4,''  L'',i4,''   Q_A '',f18.10,
     >     ''   Q_B '',d18.10)') itscfcont,li,qva(li),qvb(li)
         end if
      end do
      call flush(6)
c
c -----------------------------------------------------------------------
c -- Write spin axes ----------------------------------------------------
c -----------------------------------------------------------------------
      write(6,'(/"    Spin axes")')
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
c -----------------------------------------------------------------------
c -- Write spin moments -------------------------------------------------
c -----------------------------------------------------------------------
      write(6,'(/"    Spin moments")')
      do li=1,nintfc
         cpalay=1.d0-conc(li).gt.tiny
         if(.not.cpalay) then
           write(6,'('' I'',i4,''  L'',i4,''  MS  '',f18.10)')
     >     itscfcont,li,spinmoma(li)
         else
           write(6,'('' I'',i4,''  L'',i4,''  MS_A'',f18.10,
     >     ''  MS_B'',d18.10)') itscfcont,li,spinmoma(li),spinmomb(li)
         end if
      end do
      call flush(6)
c
c -----------------------------------------------------------------------
c -- Write orbital moments -------------------------------------------------
c -----------------------------------------------------------------------
      write(6,'(/"    Orbital moments")')
      do li=1,nintfc
         cpalay=1.d0-conc(li).gt.tiny
         if(.not.cpalay) then
           write(6,'('' I'',i4,''  L'',i4,''  MO  '',f18.10)')
     >     itscfcont,li,orbmoma(li)
         else
           write(6,'('' I'',i4,''  L'',i4,''  MO_A'',f18.10,
     >     ''  MO_B'',d18.10)') itscfcont,li,orbmoma(li),orbmomb(li)
         end if
      end do
      call flush(6)
c
      return
      end
