c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine locquant(
     > para,ce,we,lmax,nintfc,rel,conc,v0,
     > idpota,vra,bra,idpotb,vrb,brb,dx,ns,rs,taua,taub,
     > dosa,doslma,qvpa,qva,enba,qmoma,rhova,
     > dosb,doslmb,qvpb,qvb,enbb,qmomb,rhovb)
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
      dimension wr(nrad)
      dimension vra(nrad,mintfc),vrb(nrad,mintfc)
      dimension bra(nrad,mintfc),brb(nrad,mintfc)
      dimension dx(mintfc),ns(mintfc),rs(mintfc)
      dimension dosa(0:lmaxp,mintfc,2),dosb(0:lmaxp,mintfc,2)
      dimension doslma(lmmaxp,mintfc,2),doslmb(lmmaxp,mintfc,2)
      dimension qvpa(0:lmaxp,mintfc,2),qvpb(0:lmaxp,mintfc,2)
      dimension qva(mintfc,2),qvb(mintfc,2)
      dimension enba(mintfc,2),enbb(mintfc,2)
      dimension rhova(nrad,mintfc,2),rhovb(nrad,mintfc,2)
c
      complex*16 taua(lmmaxp,lmmaxp,2,mintfc)
      complex*16 taub(lmmaxp,lmmaxp,2,mintfc)
c
      complex*16 ce,we,cih,ttmp(0:lmaxp)
      complex*16 fint(0:lmaxp),gint(0:lmaxp),fpint(lmaxp)
      complex*16 yl1(nrad,0:lmaxp),yl2(nrad,0:lmaxp) 
      complex*16 zqmom(lmmaxp,lmsup),zdos(0:lmaxp),zrho(nrad)
      complex*16 qmoma(lmsup,mintfc,2),qmomb(lmsup,mintfc,2)
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
c **************************
c * loop over spin indices *
c **************************
      do is=1,2
c
c ********************
c * loop over layers * 
c ********************
      do li=1,nintfc
c
         cpalay=1.0d0-conc(li).gt.tiny
         do i=1,ns(li)
           wr(i)=vra(i,li)+(3-2*is)*bra(i,li)
         end do
c
c Compute scattering solutions
c        --------------------------------------------------------
         call wafu(ce,lmax,idpota(li),v0,wr,dx(li),ns(li),
     >   rs(li),1,rel,ttmp,fint,gint,fpint,yl1,yl2)
c        --------------------------------------------------------
c
c Radial distribution of charge
c        -----------------------------------------------------
         call dens(taua(1,1,is,li),yl1,yl2,zrho,ns(li),lmax,para)
c        -----------------------------------------------------
c
c Density of multipole moments
         i=0
         do lam=0,lmaxs
         do mu=-lam,lam
            i=i+1
c           --------------------------------------------------
            call moment(lmax,rs(li),dx(li),ns(li),yl1,yl2,
     >      gacoeff(1,1,i),lam,taua(1,1,is,li),zqmom(1,i),1,para)
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
               doslma(i,li,is)=dimag(zqmom(i,1))
            end do
            dosa(l,li,is)=dimag(zdos(l))
         enddo
c      
c add contribution to contour integral with respect to the energy
c
         do l=0,lmax
            qvpa(l,li,is)=qvpa(l,li,is)+dimag(we*zdos(l))
            qva(li,is)=qva(li,is)+dimag(we*zdos(l))
            enba(li,is)=enba(li,is)+dimag(we*ce*zdos(l))
         end do
         do irad=1,ns(li)
            rhova(irad,li,is)=rhova(irad,li,is)+dimag(we*zrho(irad))
         end do
c
         do lam=0,lmaxs
         do mu=-lam,lam
            i=lam*(lam+1)+mu+1
            ii=lam*(lam+1)-mu+1
            do k=1,lmmax
               qmoma(i,li,is)=qmoma(i,li,is)+cih*(we*zqmom(k,i)-
     >                        mpar(mu)*dconjg(we*zqmom(k,ii)))
            end do
         end do
         end do
c
         if(itest.ge.2) then
           write(6,'('' DOS: Layer'',i2)') li
           write(6,'(4d15.6)') (dosa(l,li,1),l=0,lmax)
           write(6,'(4d15.6)') (dosa(l,li,2),l=0,lmax)
         end if
c
         if(cpalay) then
c+------------+
c+ BIG CPA IF +
c+------------+
c
         do i=1,ns(li)
           wr(i)=vrb(i,li)+(3-2*is)*brb(i,li)
         end do
c
c Compute scattering solutions
c        --------------------------------------------------------
         call wafu(ce,lmax,idpotb(li),v0,wr,dx(li),ns(li),
     >   rs(li),1,rel,ttmp,fint,gint,fpint,yl1,yl2)
c        --------------------------------------------------------
c
c Radial distribution of charge
c        -----------------------------------------------------
         call dens(taub(1,1,is,li),yl1,yl2,zrho,ns(li),lmax,para)
c        -----------------------------------------------------
c
c Density of multipole moments
         i=0
         do lam=0,lmaxs
         do mu=-lam,lam
            i=i+1
c           --------------------------------------------------
            call moment(lmax,rs(li),dx(li),ns(li),yl1,yl2,
     >      gacoeff(1,1,i),lam,taub(1,1,is,li),zqmom(1,i),1,para)
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
               doslmb(i,li,is)=dimag(zqmom(i,1))
            end do
            dosb(l,li,is)=dimag(zdos(l))
         enddo
c
c add contribution to contour integral with respect to the energy
c
         do l=0,lmax
            qvpb(l,li,is)=qvpb(l,li,is)+dimag(we*zdos(l))
            qvb(li,is)=qvb(li,is)+dimag(we*zdos(l))
            enbb(li,is)=enbb(li,is)+dimag(we*ce*zdos(l))
         end do
         do irad=1,ns(li)
            rhovb(irad,li,is)=rhovb(irad,li,is)+dimag(we*zrho(irad))
         end do
c
         do lam=0,lmaxs
         do mu=-lam,lam
            i=lam*(lam+1)+mu+1
            ii=lam*(lam+1)-mu+1
            do k=1,lmmax
               qmomb(i,li,is)=qmomb(i,li,is)+cih*(we*zqmom(k,i)-
     >                        mpar(mu)*dconjg(we*zqmom(k,ii)))
            end do
         end do
         end do
c
         else
c
         do i=1,lmmax
           doslmb(i,li,is)=doslma(i,li,is)
         end do
         do l=0,lmax
           dosb(l,li,is)=dosa(l,li,is)
           qvpb(l,li,is)=qvpa(l,li,is)
         end do
         qvb(li,is)=qva(li,is)
         enbb(li,is)=enba(li,is)
         do i=1,lmmaxs
            qmomb(i,li,is)=qmoma(i,li,is)
         end do
         do irad=1,ns(li)
            rhovb(irad,li,is)=rhova(irad,li,is)
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
      end do
c ******************************
c * end loop over spin indices *
c ******************************

      return
      end
