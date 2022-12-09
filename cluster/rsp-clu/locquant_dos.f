c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine locquant_dos(
c     ===================
     > ce,lmax,nintfc,wrel,lms,sxa,v0,
     > idpota,vra,bra,bopra,dx,ns,rs,
     > taua, dosa)
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      integer mmax
      parameter (mmax=mimp)
c
      logical wrel,lms
c
      character*10 idpota(mmax)
c
      integer nuz(kmymaxp)
      integer indz(nuzp,kmymaxp)
      integer mpar(-lsup:lsup)
      integer ns(mmax)
c
      real*8 vra(nrad,mmax)
      real*8 bra(nrad,mmax)
      real*8 bopra(nrad,2,mmax)
      real*8 rs(mmax)
      real*8 dx(mmax)
c
c
      real*8 dosa(kmymaxp,mmax)
c
      real*8 sxa(mmax)
c
      complex*16 ce
      complex*16 psq
c
c ======================================================================
c
      complex*16 taua(kmymaxp,kmymaxp,mmax)
      complex*16 tm(kmymaxp,kmymaxp)
c
      complex*16 gz(nrad,nuzp,kmymaxp),fz(nrad,nuzp,kmymaxp)
      complex*16 gj(nrad,nuzp,kmymaxp),fj(nrad,nuzp,kmymaxp)
c
      complex*16 zdos(kmymaxp)
      complex*16 zqmom(kmymaxp,lmsup)
c
      complex*16 rgacoeff(kmymaxp,kmymaxp,lmsup)
      common/rggaunt/rgacoeff
c
      data tol/1.0d-8/,cih/(0.0d0,-0.5d0)/
      data tiny/1.0d-6/
c
c ********************
c initialize constants
c ********************
c
c---> c in rydberg units:
      c=274.072d0
      if(.not.wrel) then
        psq=ce+ce*ce/(c*c)
      else
        psq=ce
      end if
c
      nl=lmax+1
      nl2=nl*nl
      kmax=2*lmax+1
      kmymax=2*nl2
      lmaxs=2*lmax
      lmmaxs=(lmaxs+1)*(lmaxs+1)
c
      mpar(0)=1
      do m=1,lmaxs
        mpar(m)=-mpar(m-1)
        mpar(-m)=mpar(m)
      end do
c
c ***************************************************
c * loop over layers to compute physical quantities *
c ***************************************************
c 
      do li=1,nintfc
!        cpalay=(1.d0-conc(li)).gt.tiny
         write(8,*) 'Impurity:',li
c
c Compute scattering solutions 
c        --------------------------------------------------------
         call wafu(ce,psq,lmax,idpota(li),v0,vra(1,li),bra(1,li),
     >             bopra(1,1,li),dx(li),ns(li),rs(li),
     >             tm,gz,fz,gj,fj,nuz,indz,1,sxa(li))
c        --------------------------------------------------------
c
c Density of multipole moments
         i=0
         do lam=0,lmaxs
         do mu=-lam,lam
           i=i+1
c        --------------------------------------------
           call moment(lmax,rs(li),dx(li),ns(li),
     >                 gz,fz,gj,fj,nuz,indz,
     >                 rgacoeff(1,1,i),lam,taua(1,1,li),
     >                 zqmom(1,i),lms,1)
c        --------------------------------------------
         end do
         end do
         do kmy=1,kmymax
            zdos(kmy)=zqmom(kmy,1)
         enddo
c
c        -------------------------------------------------
c
         do k=1,kmymax
           dosa(k,li)=dimag(zdos(k))
         end do
         write(8,'(2e14.6,2x,18e14.6)')
     >           dreal(ce),dimag(ce),(dosa(k,li),k=1,kmymax)
      end do
c ************************
c * end loop over layers *
c ************************
c
      return
      end
