c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine rcore(for006,idpot,z,za,nws,rws,dx,
     >                 nuc,iskip,norb,nes,test,den,nqn,nk,nel,
     >                 zcore,titre,nql,dfl,dq1,rhoc,encore)
c
      implicit real*8(a-h,o-z)
c
      include '../param.h'
c
      character*4 titre(iorb)
      character*10 idpot
      character*30 for006
      dimension za(nrad),rhoc(nrad)
      dimension den(iorb),dq1(iorb),dfl(iorb),nqn(iorb),nql(iorb),
     >          nk(iorb),nel(iorb)
      dimension dv(nrad),dr(nrad),dp(nrad),dq(nrad)
      dimension dgc(nrad,iorb),dpc(nrad,iorb)
      dimension xx(nrad),yy(nrad)
      dimension cwf(iorb)
c
c
c     ***************************************************************
c     * calculate fully relativistically one-electron orbitals      *
c     * energies for core one-electron orbitals                     *
c     * core charge density                                         *
c     * kinetic energy contribution from the core                   *
c     ***************************************************************
c
      imax0=0
c     np=nrad
      np=nws+50
      if(np.ge.nrad) np=nrad
c
      dpas=dx
      x=dlog(rws)-(nws-1)*dx
      do i=1,np
        dr(i)=dexp(x)
        dv(i)=0.5d0*za(i)/dr(i)
        rhoc(i)=0.0d0
        x=x+dx
      end do
c
      if(iskip.ge.2) then
        write(6,'( '' Potential for:'',a10)') idpot
        write(6,'(4d20.10)') (dr(k),dv(k),k=1,np)
      end if
c
c     now do core part
c
      do i=1,np
        xx(i)=0.0d0
        do j=1,norb
          dgc(i,j)=0.0d0
          dpc(i,j)=0.0d0
        end do
      end do
c
c
c loop over core orbitals
c
      zcore=0.d0
      do j=1,norb
        if(iskip.ge.2) then
          write(6,'('' orbital:'',i3)') j
        end if
c
        zcore=zcore+nel(j)
        do kk=1,np
          dp(kk)=0.0d0
          dq(kk)=0.0d0
        end do
c
        call resld(nqn(j),nql(j),nk(j),imax,den(j),dfl(j),dq1(j),j,
     >             dv,dr,dp,dq,dpas,z,nes,test,np,nuc,iskip)
c
        if(imax.gt.imax0) imax0=imax
c
c   calculate core charge density
c
        do i=1,imax
          dgc(i,j)=dp(i)
          dpc(i,j)=dq(i)
          xx(i)=xx(i)+nel(j)*(dp(i)*dp(i)+dq(i)*dq(i))
          rhoc(i)=xx(i)
        end do
c
      end do
c
      np=imax0
c
      do i=imax+1,np
        rhoc(i)=0.0d0
      end do
c
      if(iskip.ge.1) then
c
c     write out the precious results
c
      write(6,'(/10x,''one electron energies'')')
      do i=1,norb
        write(6,'(a10,'' e-core'',i3,d20.10,3i4,3x,i1,a4)') 
     >  idpot,i,den(i),nqn(i),nk(i),nel(i),nqn(i),titre(i)
      end do
      write(6,*)
      end if
c
c     orthonormality relations
c
      if(iskip.ge.2) then
        write(6,'(/10x,''orthogonality relations'')')
        do 50 i=1,norb
        do 50 j=i,norb
          if(nql(i).ne.nql(j)) go to 50
          if(nk(i).ne.nk(j)) go to 50
          do k=1,np
            yy(k)=dpc(k,i)*dpc(k,j)+dgc(k,i)*dgc(k,j)
          end do
          anorm=rsimp(yy,dr,nws,dpas)
          bnorm=rsimp(yy,dr,np-1,dpas)
          if(i.eq.j) cwf(i)=anorm
          if(iskip.eq.0) go to 50
          write(6,'(''('',i2,a4,'','',i2,a4,'')'',3x,2e15.8)') 
     >    nqn(i),titre(i),nqn(j),titre(j),bnorm,anorm
   50   continue
      end if
c
c     calculate number of core electrons
c
      do k=1,np
        xx(k)=rhoc(k)
      end do
      anorm=rsimp(xx,dr,nws,dpas)
      if(iskip.ge.1) 
     >write(6,'(/10x,''number of core-electrons in atomic sphere ''
     >          /10x,''before normalization: '', e15.8)') anorm
c
c     normalize charge density to zcore within WS sphere
c
      anorm=zcore/anorm
      do k=1,np
        rhoc(k)=rhoc(k)*anorm
        xx(k)=rhoc(k)
      end do
      anorm=rsimp(xx,dr,nws,dpas)
      if(iskip.ge.1) 
     >write(6,'(10x,''after normalization:  '',e15.8/)') anorm
c
c     calculate core part of total energy
c
      encore=0.d0
      do i=1,norb
        encore = encore + nel(i)*den(i)
      end do
c
c  write out core charge densities
c
c     write(6,'(4d20.12)') (rhoc(j),j=1,nws)
c
c  write out core wavefunctions
c
c      do i=1,norb
c        write(6,'(i1,a4)') nqn(i),titre(i)
c        write(6,'(d20.10)') den(i)
c        write(6,'(i4)') nk(i)
c        write(6,'(i4)') nws+1
c        do j=1,nws+1
c          write(6,'(3e20.10)') dr(j),dgc(j,i),dpc(j,i)
c        end do
c      end do
c
      return
      end
