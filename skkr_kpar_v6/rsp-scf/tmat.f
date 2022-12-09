c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine tmat(ce,psq,lmax,idpot,v0,vr,br,bopr,
     >                dx,ns,rs,ptminv,tminv,socsc,alphakkr)
c======================
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      character*10 idpot
      dimension vr(nrad),br(nrad),bopr(nrad,2)
      dimension nuz(kmymaxp),indz(nuzp,kmymaxp)
      dimension ldex(50)
      complex*16 ce,psq,detl,alphakkr(0:lmax),sum1,sum2
      complex*16 tminv(kmymaxp,kmymaxp),tmk(kmaxp),tmpk(kmaxp)
      complex*16 tmatx(kmymaxp,kmymaxp),ptminv(kmymaxp,kmymaxp)
      complex*16 gz(nrad,nuzp,kmymaxp),fz(nrad,nuzp,kmymaxp)
      complex*16 gj(nrad,nuzp,kmymaxp),fj(nrad,nuzp,kmymaxp)
      complex*16 xlms(lmmaxp,lmmaxp,2,2)
c
      common/test/itest
      data ldex/0,0,
     *          1,1,1,1,1,1,
     *          2,2,2,2,2,2,2,2,2,2,           
     *          3,3,3,3,3,3,3,3,3,3,3,3,3,3,        
     *          4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4/      
      data tol/1.0d-15/,tiny/0.001/
c
      pi=4.d0*datan(1.d0)
      lmmax=(lmax+1)*(lmax+1)
      kmax=2*lmax+1
      kmymax=2*lmmax
c
      call czero(ptminv,kmymaxp*kmymaxp)
      call czero(tminv,kmymaxp*kmymaxp)
c
c inverse of physical t-matrix
c     ------------------------------------------------------
      call wafu(ce,psq,lmax,idpot,v0,vr,br,bopr,dx,ns,rs,
     >          ptminv,gz,fz,gj,fj,nuz,indz,0,socsc)
c     ------------------------------------------------------
c
      call symmat(ptminv,kmymax,kmymaxp)
c
      if(socsc.lt.tiny) then
        call matlms(xlms,ptminv,lmax)
        call czero(xlms(1,1,1,2),lmmaxp*lmmaxp)
        call czero(xlms(1,1,2,1),lmmaxp*lmmaxp)
        lm=0
        do l=0,lmax
          sum1=(0.d0,0.d0)
          sum2=(0.d0,0.d0)
          do m=-l,l
            lm=lm+1
            sum1=sum1+xlms(lm,lm,1,1)
            sum2=sum2+xlms(lm,lm,2,2)
          end do
          lm=lm-2*l-1
          do m=-l,l
            lm=lm+1
            xlms(lm,lm,1,1)=sum1/(2.0d0*float(l)+1)
            xlms(lm,lm,2,2)=sum2/(2.0d0*float(l)+1)
          end do
        end do
        call matkmy(ptminv,xlms,lmax)
      end if
c
      if(itest.ge.2) then
         write(6,'(2x,a,''      e='',3f10.6)') idpot,ce,v0
         write(6,*) ' SOC:',socsc
         write(6,*) '  unscreened tm1 matrix'
         call outmat1(ptminv,kmymax,kmymax,kmymaxp,tol,6)
      end if
c
c physical t-matrix
      call repl(tmatx,ptminv,kmymax,kmymaxp)
      call gjinv(tmatx,kmymax,kmymaxp,detl)
c
c screening-transformation on t-matrix
      do kmy=1,kmymax
         l=ldex(kmy)
         tmatx(kmy,kmy)=tmatx(kmy,kmy)-alphakkr(l)
      end do
c
      if(itest.ge.2) then
         write(6,*) '  screened t - matrix'
         call outmat1(tmatx,kmymax,kmymax,kmymaxp,tol,6)
      end if
c
c inverse screened t-matrix
      call repl(tminv,tmatx,kmymax,kmymaxp)
      call gjinv(tminv,kmymax,kmymaxp,detl)
c
      if(itest.ge.2) then
         write(6,*) '  screened tm1 matrix'
         call outmat1(tminv,kmymax,kmymax,kmymaxp,tol,6)
      end if
c
      return
      end
