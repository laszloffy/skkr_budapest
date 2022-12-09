c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine tmat(ce,cl,soc,lmax,idpot,v0,E_Fermi,delta,
     >                singrat,urat,drat,vr,br,
     >                dx,rs,ns,ptminv,tminv,alphakkr,alphakkrh,rb)
c======================
c
      implicit none
c
      include '../param.h'
c
      character*10 idpot
      integer lmax,itest,lmmax,kmax,kmymax,l,lm,m,ns,kmy
      real*8 E_Fermi
      real*8 soc,cl,v0,pi,tol,tiny,dx,rs
      real*8 vr(nrad),br(nrad),rb(3)
      real*8 singrat,urat,drat
      complex*16 ce,detl,alphakkr(0:lmaxp),sum1,sum2
      complex*16 alphakkrh(0:lmaxp)
      complex*16 tminv(dbogomaxp,dbogomaxp)
      complex*16 tmatx(dbogomaxp,dbogomaxp),ptminv(dbogomaxp,dbogomaxp)
      complex*16 gz(dbogomaxp,dbogomaxp,nrad)
      complex*16 fz(dbogomaxp,dbogomaxp,nrad)
      complex*16 gj(dbogomaxp,dbogomaxp,nrad)
      complex*16 fj(dbogomaxp,dbogomaxp,nrad)
      complex*16 glz(dbogomaxp,dbogomaxp,nrad)
      complex*16 flz(dbogomaxp,dbogomaxp,nrad)
      complex*16 glj(dbogomaxp,dbogomaxp,nrad)
      complex*16 flj(dbogomaxp,dbogomaxp,nrad)
      complex*16 xlms(lmmaxp,lmmaxp,2,2)
      complex*16 delta(nrad)
c
      logical tautest
c
c     xlms dimenzioja !!!!
c
      common/test/itest
c
      data tol/1.0d-15/,tiny/0.001/
c
      pi=4.d0*datan(1.d0)
      lmmax=(lmax+1)*(lmax+1)
      kmax=2*lmax+1
      kmymax=2*lmmax
c
      ptminv=(0.d0,0.d0)
      tminv=(0.d0,0.d0)
      xlms=(0.d0,0.d0)
c
      tautest=.false.
c
      if(tautest) then
       write(6,*) ' <tmat> : input parameters'
       write(6,*) 'ce=',ce
       write(6,*) 'cl=',cl
       write(6,*) 'soc=',soc
       write(6,*) 'lmax=',lmax
       write(6,*) 'idpot=',idpot
       write(6,*) 'v0=',v0
       write(6,*) 'E_Fermi=',E_Fermi
       write(6,*) 'delta=',delta
       write(6,*) 'vr=',vr
       write(6,*) 'br=',br
       write(6,*) 'dx=',dx
       write(6,*) 'rs=',rs
       write(6,*) 'ns=',ns
       write(6,*) 'alphakkr=',alphakkr
       write(6,*) 'alphakkrh=',alphakkrh
       write(6,*) 'rb=',rb
      end if
c     write(6,'(2x,a,''      e='',2f10.6)') idpot,ce
c
c calculate inverse of physical t-matrix
c     ----------------------------------------------------
      call wafu(ce,cl,soc,lmax,idpot,v0,E_Fermi,
     >          singrat,urat,drat,delta,
     >          vr,br,rb,dx,rs,ns,
     >          ptminv,gz,fz,gj,fj,glz,flz,glj,flj,0)
c     ----------------------------------------------------
c
      !let's not call symmat, ptminv is not supposed to be symmetric unless B in xz plane
      !call symmat(ptminv,kmymax,kmymaxp)
c
c      if(soc.lt.tiny) then
c        call matlms(xlms,ptminv,lmax)
c        lm=0
c        do l=0,lmax
c          sum1=(0.d0,0.d0)
c          sum2=(0.d0,0.d0)
c          do m=-l,l
c            lm=lm+1
c            sum1=sum1+xlms(lm,lm,1,1)
c            sum2=sum2+xlms(lm,lm,2,2)
c          end do
c          lm=lm-2*l-1
c          do m=-l,l
c            lm=lm+1
c            xlms(lm,lm,1,1)=sum1/(2.0d0*float(l)+1)
c            xlms(lm,lm,2,2)=sum2/(2.0d0*float(l)+1)
c          end do
c        end do
c        call matkmy(ptminv,xlms,lmax)
c      end if
c
      if(itest.ge.2.or.tautest) then
         write(6,'(2x,a,''      e='',3f10.6)') idpot,ce,v0
         write(6,*) ' SOC:',soc
         write(6,*) '  unscreened tm1 matrix'
         call outmat1(ptminv,2*kmymax,2*kmymax,dbogomaxp,tol,6)
      end if
c
c physical t-matrix
      call repl(tmatx,ptminv,2*kmymax,dbogomaxp)
      call gjinv(tmatx,2*kmymax,dbogomaxp,detl)
c
c screening-transformation on t-matrix
      do kmy=1,kmymax
         l=ldex(kmy)
         tmatx(kmy,kmy)=tmatx(kmy,kmy)-alphakkr(l)
         tmatx(kmymax+kmy,kmymax+kmy)=
     >           tmatx(kmymax+kmy,kmymax+kmy)-alphakkrh(l)
      end do
c
      if(itest.ge.2) then
         write(6,*) '  screened t - matrix'
         call outmat1(tmatx,2*kmymax,2*kmymax,dbogomaxp,tol,6)
      end if
c
c inverse screened t-matrix
      call repl(tminv,tmatx,2*kmymax,dbogomaxp)
      call gjinv(tminv,2*kmymax,dbogomaxp,detl)
c
      if(itest.ge.2.or.tautest) then
         write(6,*) '  screened tm1 matrix'
         call outmat1(tminv,2*kmymax,2*kmymax,dbogomaxp,tol,6)
      end if
c
      return
      end
