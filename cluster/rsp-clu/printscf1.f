c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine printscf1(
     > nintfc,itscfcont,enbifc,qvifc,omifc,
     > qva,qca,za,enba,denba,efermi,
     > rba,spin_magva,orb_magva,
     > spinmoma,orbmoma,th0a,th1a,ph0a,ph1a)
c
!     implicit real*8 (a-h,o-z)
      implicit none
      include '../param.h' 
c
      integer mmax
      parameter (mmax=mimp)
c
      real*8 rba(3,mmax)
      real*8 qva(mmax)
      real*8 qca(mmax)
      real*8 za(mmax)
      real*8 enba(mmax)
      real*8 denba(mmax)
      real*8 spin_magva(mmax,3)
      real*8 orb_magva(mmax,3)
      real*8 veca(mmax,3)
c
      real*8 spinmoma(mmax)
      real*8 orbmoma(mmax)
      real*8 th0a(mmax)
      real*8 th1a(mmax)
      real*8 ph0a(mmax)
      real*8 ph1a(mmax)
c
      real*8 enbifc
      real*8 efermi
      real*8 qvifc
      real*8 omifc
      real*8 enb
      real*8 fac
      real*8 pi
      real*8 om
      real*8 q
c
      integer itscfcont
      integer nintfc
      integer li
c
      real*8 tiny
      data tiny/1.0d-6/ 
c ---------------------------------------------------------------------
      pi=dacos(-1.0d0)
      fac=180.d0/pi
c
!     if(linbw) then
!       write(6,'(/''   '',t13,''Eband'',t33,''  Q'',t53,'' Omega'')')
!       write(6,'(3x,3d20.12)') enbifc,qvifc,omifc
!     end if
c
        write(6,'(/''   '',t13,''Eb'')')
c
      qvifc=0.d0
      enbifc=0.d0
      omifc=0.d0
      do li=1,nintfc
        q=qva(li)
        enb=enba(li)
        om=enb-efermi*q
c       om=enb
        qvifc=qvifc+q
        enbifc=enbifc+enb
        omifc=omifc+om
c       if(li.eq.4) then
c         cth=rba(3,li)
c         sth=dsqrt(1.0d0-cth*cth)
c         th=dacos(cth)
c         cph=rba(1,li)/sth
c         sph=rba(2,li)/sth
c         ph=dacos(cph)
c         if(sph.lt.0.d0) ph=2*pi-ph
c         write(6,'(i3,2d15.6,d20.12)') li,th,ph,
c    >    (enba(li)-efermi*qva(li))*rba(3,li)
c    >                            enbb(li)-efermi*qvb(li),om
c         end if
        write(6,'(i3,3d20.12)') li,enba(li)-efermi*qva(li),om
c       write(6,'(i3,3d20.12)') li,enba(li),enbb(li),om
      end do
c
      write(6,'(''-----------------------------------------------'',
     >''----------------'')')
      write(6,'(t44,d20.12)') omifc
      call flush(6)
c
      write(6,'(/''   '',t13,''DEb'')')
      do li=1,nintfc
        write(6,'(i3,2d20.12)') li,denba(li)
      end do
c
      write(6,*)
      do li=1,nintfc
        write(6,'('' I'',i4,''  L'',i4,''   Q   '',d18.10)')
     >  itscfcont,li,qva(li)
      end do
      call flush(6)
c
      write(6,*)
      do li=1,nintfc
        write(6,'('' I'',i4,''  L'',i4,''  Bth   '',2f15.5)')
     >  itscfcont,li,th0a(li)*fac,th1a(li)*fac
      end do
      call flush(6)
c
      write(6,*)
      do li=1,nintfc
        write(6,'('' I'',i4,''  L'',i4,''  Bph   '',2f15.5)')
     >  itscfcont,li,ph0a(li)*fac,ph1a(li)*fac
      end do
      call flush(6)
c
      write(6,*)
      do li=1,nintfc
        write(6,'('' I'',i4,''  L'',i4,''  MS  '',4d18.10)')
     >  itscfcont,li,spinmoma(li),spin_magva(li,1),
     >  spin_magva(li,2),spin_magva(li,3)
      end do
      call flush(6)
c
      write(6,*)
      do li=1,nintfc
        write(6,'('' I'',i4,''  L'',i4,''  MO  '',4d18.10)')
     >  itscfcont,li,orbmoma(li),orb_magva(li,1),
     >  orb_magva(li,2),orb_magva(li,3)
      end do
      call flush(6)
c
      return
      end
