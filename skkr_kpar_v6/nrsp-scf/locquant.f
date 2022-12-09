c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine locquant(
     > para,ce,we,lmax,nintfc,rel,conc,v0,
     > idpota,vra,idpotb,vrb,dx,ns,rs,
     > taua,taub,taukdiag,nk,wk,
     > dosa,doslma,qvpa,qva,enba,qmoma,rhova,
     > dosb,doslmb,qvpb,qvb,enbb,qmomb,rhovb,
     > spdos,sdos)
c
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
      parameter(mdim=mdimnr)
c
      logical para,rel,cpalay
c 
      character*10 idpota(mintfc),idpotb(mintfc)
c
      dimension conc(mintfc)
c
      dimension mpar(-lsup:lsup)
      dimension vra(nrad,mintfc),vrb(nrad,mintfc)
      dimension dx(mintfc),ns(mintfc),rs(mintfc)
      dimension dosa(0:lmaxp,mintfc),dosb(0:lmaxp,mintfc)
      dimension doslma(lmmaxp,mintfc),doslmb(lmmaxp,mintfc)
      dimension spdos(mkpar,0:lmaxp,mintfc),sdos(0:lmaxp,mintfc)
      dimension qvpa(0:lmaxp,mintfc),qvpb(0:lmaxp,mintfc)
      dimension qva(mintfc),qvb(mintfc)
      dimension enba(mintfc),enbb(mintfc)
      dimension rhova(nrad,mintfc),rhovb(nrad,mintfc)
      dimension wk(mkpar)
c
      complex*16 taua(lmmaxp,lmmaxp,mintfc)
      complex*16 taub(lmmaxp,lmmaxp,mintfc)
      complex*16 taukdiag(mkpar,0:lmaxp,mintfc)
c
      complex*16 ce,we,cih,ttmp(0:lmaxp)
      complex*16 fint(0:lmaxp),gint(0:lmaxp),fpint(lmaxp)
      complex*16 yl1(nrad,0:lmaxp),yl2(nrad,0:lmaxp) 
      complex*16 zqmom(lmmaxp,lmsup),zdos(0:lmaxp),zrho(nrad)
      complex*16 qmoma(lmsup,mintfc),qmomb(lmsup,mintfc)
      complex*16 u(lmmaxp,lmmaxp),up(lmmaxp,lmmaxp)
c
      common/test/itest
      common/lay2d/cvec(mtotal,3),nextra,nbulkl,nbulkr,
     &             nprc,ninprc(0:mprc+1)
c
      complex*16 gacoeff(lmmaxp,lmmaxp,lmsup)
      common/ggaunt/gacoeff
c
      data outtol/1.0d-12/
      data tiny/1.0d-6/
      data cih/(0.0d0,-0.5d0)/
c
c ********************
c initialize constants
c ********************
c
      pi=4.0d0*datan(1.0d0)
      fact=-1.0d0/pi
      nl=lmax+1
      lmmax=nl*nl
      lmaxs=2*lmax
      lmmaxs=(lmaxs+1)*(lmaxs+1)
c
      mpar(0)=1
      do m=1,lmaxs
        mpar(m)=-mpar(m-1)
        mpar(-m)=mpar(m)
      end do
c
      call ytrafo(u,up,lmax)
c
c ***************************************************
c * loop over layers to compute physical quantities *
c ***************************************************
c
      do li=1,nintfc
         cpalay=(1.d0-conc(li)).gt.tiny
c
         call tripmt(u,taua(1,1,li),up,lmmax,lmmax,lmmaxp)
c
c Compute scattering solutions
c        --------------------------------------------------------
         call wafu(ce,lmax,idpota(li),v0,vra(1,li),dx(li),ns(li),
     >   rs(li),1,rel,ttmp,fint,gint,fpint,yl1,yl2)
c        --------------------------------------------------------
c
c Radial distribution of charge
c        -----------------------------------------------------
         call dens(taua(1,1,li),yl1,yl2,zrho,ns(li),lmax,para)
c        -----------------------------------------------------
c
c Spectral density for non-CPA 
         do l=0,lmax
           sdos(l,li)=0.0d0 
           do ik=1,nk
             spdos(ik,l,li)=
     >       fact*dimag(taukdiag(ik,l,li)*fint(l)-gint(l))
             sdos(l,li)=sdos(l,li)+wk(ik)*spdos(ik,l,li)
           end do
         end do
c
c Density of multipole moments
         i=0
         do lam=0,lmaxs
         do mu=-lam,lam
            i=i+1
c           --------------------------------------------------
            call moment(lmax,rs(li),dx(li),ns(li),yl1,yl2,
     >      gacoeff(1,1,i),lam,taua(1,1,li),zqmom(1,i),1,para)
c           --------------------------------------------------
         end do
         end do
c
         i=0
         do l=0,lmax
            zdos(l)=(0.0d0,0.d0)
            do m=-l,l
               i=i+1
               zdos(l)=zdos(l)+zqmom(i,1)
               doslma(i,li)=dimag(zqmom(i,1))
            end do
            dosa(l,li)=dimag(zdos(l))
         enddo
c      
c add contribution to contour integral with respect to the energy
c
         do l=0,lmax
            qvpa(l,li)=qvpa(l,li)+dimag(we*zdos(l))
            qva(li)=qva(li)+dimag(we*zdos(l))
            enba(li)=enba(li)+dimag(we*ce*zdos(l))
         end do
         do irad=1,ns(li)
            rhova(irad,li)=rhova(irad,li)+dimag(we*zrho(irad))
         end do
c
         do lam=0,lmaxs
         do mu=-lam,lam
            i=lam*(lam+1)+mu+1
            ii=lam*(lam+1)-mu+1
            do k=1,lmmax
               qmoma(i,li)=qmoma(i,li)+cih*(we*zqmom(k,i)-
     >                     mpar(mu)*dconjg(we*zqmom(k,ii)))
            end do
         end do
         end do
c
         if(itest.ge.2) then
           write(6,'('' DOS: Layer'',i2)') li
           write(6,'(4d15.6)') (dosa(l,li),l=0,lmax)
         end if
c
         if(cpalay) then
c+------------+
c+ BIG CPA IF +
c+------------+
c
         call tripmt(u,taub(1,1,li),up,lmmax,lmmax,lmmaxp)
c
c Compute scattering solutions
c        --------------------------------------------------------
         call wafu(ce,lmax,idpotb(li),v0,vrb(1,li),dx(li),ns(li),
     >   rs(li),1,rel,ttmp,fint,gint,fpint,yl1,yl2)
c        --------------------------------------------------------
c
c Radial distribution of charge
c        -----------------------------------------------------
         call dens(taub(1,1,li),yl1,yl2,zrho,ns(li),lmax,para)
c        -----------------------------------------------------
c
c Density of multipole moments
         i=0
         do lam=0,lmaxs
         do mu=-lam,lam
            i=i+1
c           --------------------------------------------------
            call moment(lmax,rs(li),dx(li),ns(li),yl1,yl2,
     >      gacoeff(1,1,i),lam,taub(1,1,li),zqmom(1,i),1,para)
c           --------------------------------------------------
         end do
         end do
c
         i=0
         do l=0,lmax
            zdos(l)=(0.0d0,0.d0)
            do m=-l,l
               i=i+1
               zdos(l)=zdos(l)+zqmom(i,1)
               doslmb(i,li)=dimag(zqmom(i,1))
            end do
            dosb(l,li)=dimag(zdos(l))
         enddo
c
c add contribution to contour integral with respect to the energy
c
         do l=0,lmax
            qvpb(l,li)=qvpb(l,li)+dimag(we*zdos(l))
            qvb(li)=qvb(li)+dimag(we*zdos(l))
            enbb(li)=enbb(li)+dimag(we*ce*zdos(l))
         end do
         do irad=1,ns(li)
            rhovb(irad,li)=rhovb(irad,li)+dimag(we*zrho(irad))
         end do
c
         do lam=0,lmaxs
         do mu=-lam,lam
            i=lam*(lam+1)+mu+1
            ii=lam*(lam+1)-mu+1
            do k=1,lmmax
               qmomb(i,li)=qmomb(i,li)+cih*(we*zqmom(k,i)-
     >                     mpar(mu)*dconjg(we*zqmom(k,ii)))
            end do
         end do
         end do
c
         else
c
         do i=1,lmmax
           doslmb(i,li)=doslma(i,li)
         end do
         do l=0,lmax
           dosb(l,li)=dosa(l,li)
           qvpb(l,li)=qvpa(l,li)
         end do
         qvb(li)=qva(li)
         enbb(li)=enba(li)
         do i=1,lmmaxs
            qmomb(i,li)=qmoma(i,li)
         end do
         do irad=1,ns(li)
            rhovb(irad,li)=rhova(irad,li)
         end do

         end if
c+----------------+
c+ END BIG CPA IF +
c+----------------+ 
c
      end do
c ************************
c * end loop over layers *
c ************************  
c
      return
      end
