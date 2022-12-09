c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine vgen(
     > itscf,itscfcont,lmax,nintfc,sigma,iesublatt,
     > rightm,bulk,vrsbulk,vrsh,nspheres,volint,
     > conc,concl,concr,rsl,rsr,
     > za,qca,qva,qmoma,qmomla,qmomra,
     > zb,qcb,qvb,qmomb,qmomlb,qmomrb,
     > rhoca,rhova,rhocb,rhovb,rhomaga,rhomagb,dx,ns,rs,
     > efermi,defermi,vra,vrb,bra,brb,v0,dv0,ferr1,ferr2,newvac,
     > enpota,enela,enxca,enpotb,enelb,enxcb,enmaga,enmagb,dlm)
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      parameter (mbig=4*nrad*mintfc+2)
c
      integer root, myrank, nprocs
      common/mpi/root,myrank,nprocs
c
      logical bulk,dlm,nspheres
      character*1 rightm
      character*2 tempr
      character*50 temp1,temp2,temp3
c
      dimension iesublatt(mintfc)
c
      dimension vra0(0:nrad),vrb0(0:nrad)
      dimension bra0(0:nrad),brb0(0:nrad)
      dimension r0(0:nrad)
c
      dimension dx(mintfc),rs(mintfc),ns(mintfc)
      dimension rsnew(mintfc)
      dimension rsr(minprc),rsl(minprc)
      dimension vra(nrad,mintfc),bra(nrad,mintfc)
      dimension vrb(nrad,mintfc),brb(nrad,mintfc)
      dimension gvra(nrad,mintfc),gvrb(nrad,mintfc)
      dimension gbra(nrad,mintfc),gbrb(nrad,mintfc)
      dimension rhoca(nrad,mintfc),rhova(nrad,mintfc)
      dimension rhocb(nrad,mintfc),rhovb(nrad,mintfc)
      dimension rhomaga(nrad,mintfc),rhomagb(nrad,mintfc)
      dimension za(mintfc),zb(mintfc)
      dimension qca(mintfc),qcb(mintfc)
      dimension qva(mintfc),qvb(mintfc)
      dimension conc(mintfc),concl(minprc),concr(minprc)
      dimension enpota(mintfc),enmaga(mintfc),enela(mintfc),
     >          enxca(mintfc)
      dimension enpotb(mintfc),enmagb(mintfc),enelb(mintfc),
     >          enxcb(mintfc)
      dimension xbr(mbig),fbr(mbig)
      dimension dfbr(mbig),pfbr(mbig),pxbr(mbig),denbr(mbig),defbr(mbig)
c
      complex*16 qmoma(lmsup,mintfc),qmomb(lmsup,mintfc)
      complex*16 qmomla(lmsup,minprc),qmomlb(lmsup,minprc)
      complex*16 qmomra(lmsup,minprc),qmomrb(lmsup,minprc)
c
      common/test/itest
      common/broy/cmix,wbr,nbrmax(100),nbrcase
c
      data tolv0/0.1d0/,anul/1.0d-15/
c
      save ibr,ibrmax,ibrcase
c
      if(bulk.and.cmix.lt.1.0d-8) return
c
c     ----------------------------------------------------------------
      call newpot(itscfcont,lmax,nintfc,sigma,
     >            rightm,bulk,vrsbulk,vrsh,nspheres,volint,
     >            conc,concl,concr,dx,ns,rs,rsnew,
     >            za,qca,qva,qmoma,qmomla,qmomra,
     >            zb,qcb,qvb,qmomb,qmomlb,qmomrb,
     >            rhoca,rhova,rhocb,rhovb,rhomaga,rhomagb,
     >            gv0,gvra,gvrb,gbra,gbrb,
     >            enpota,enmaga,enela,enxca,
     >            enpotb,enmagb,enelb,enxcb)
c     ----------------------------------------------------------------
c     write(6,'('' Newpot ready'')')
c     call flush(6)
c
      if(dlm) then
        do li=1,nintfc
        do i=1,ns(li)
          gvra(i,li)=0.5d0*(gvra(i,li)+gvrb(i,li))
          gvrb(i,li)=gvra(i,li)
          gbra(i,li)=0.5d0*(gbra(i,li)-gbrb(i,li))
          gbrb(i,li)=-gbra(i,li)
        end do
        end do
      end if
c
c     -----------------------------------------
      call sublattpot(nintfc,ns,gvra,iesublatt)
      call sublattpot(nintfc,ns,gvrb,iesublatt)
      call sublattpot(nintfc,ns,gbra,iesublatt)
      call sublattpot(nintfc,ns,gbrb,iesublatt)
c     -----------------------------------------
c
c - Potential mixing
c
      ibig=0
      do li=1,nintfc
        do i=1,ns(li)
          ibig=ibig+1
          xbr(ibig)=vra(i,li)
          fbr(ibig)=gvra(i,li)-vra(i,li)
          ibig=ibig+1
          xbr(ibig)=vrb(i,li)
          fbr(ibig)=gvrb(i,li)-vrb(i,li)
        end do
      end do
      if(bulk) then
        ibig=ibig+1
        xbr(ibig)=efermi-defermi
        fbr(ibig)=defermi
      else if(rightm.eq.'V'.and.newvac.eq.1) then
        ibig=ibig+1
        xbr(ibig)=v0
        fbr(ibig)=gv0-v0
      end if
      nbig1=ibig
c
      do li=1,nintfc
        x=dlog(rs(li))-(ns(li)-1)*dx(li)
        do i=1,ns(li)
          r=dexp(x)
          ibig=ibig+1
          xbr(ibig)=bra(i,li)/r
          fbr(ibig)=(gbra(i,li)-bra(i,li))/r
          ibig=ibig+1
          xbr(ibig)=brb(i,li)/r
          fbr(ibig)=(gbrb(i,li)-brb(i,li))/r
          x=x+dx(li)
        end do
      end do
      nbig=ibig
c
      call dott(xbr,xbr,t1,nbig1)
      call dott(fbr,fbr,err1,nbig1)
      ferr1=dsqrt(err1/t1)
      write(6,'(/'' Error for V(r) :'',3d17.6)')  t1,err1,ferr1
      call dott(xbr(nbig1+1),xbr(nbig1+1),t1,nbig-nbig1)
      call dott(fbr(nbig1+1),fbr(nbig1+1),err2,nbig-nbig1)
      if(t1.lt.anul) then
        ferr2=0.d0
      else
        ferr2=dsqrt(err2/t1)
      end if
      write(6,'('' Error for B(r) :'',3d17.6)')  t1,err2,ferr2
c
      if(myrank.le.9) then
         write(tempr,'("0",i1)') myrank
      else
         write(tempr,'(i2)') myrank
      end if

      temp1 = 'ftn21.'//tempr
      temp2 = 'ftn22.'//tempr
      temp3 = 'ftn23.'//tempr

      if(itscf.eq.1) then
        ibrcase=1
        ibrmax=nbrmax(1)
        ibr=1
      end if
      write(6,'(/'' itscf,ibrcase,ibrmax,ibr:'',4i4)')
     >itscf,ibrcase,ibrmax,ibr
c     ----------------------------------------------------------------
      call broyd(xbr,fbr,0.d0,ibr,nbig,dfbr,pfbr,pxbr,denbr,defbr,
     >           cmix,wbr,21,22,23,'ftn21','ftn22','ftn23')
c     ----------------------------------------------------------------
      if(ibr.ge.ibrmax) then
        ibr=1
        ibrcase=ibrcase+1
        ibrmax=nbrmax(ibrcase)
      else
        ibr=ibr+1
      end if
c
      ibig=0
      do li=1,nintfc
        do i=1,ns(li)
          ibig=ibig+1
          vra(i,li)=xbr(ibig)
          ibig=ibig+1
          vrb(i,li)=xbr(ibig)
        end do
      end do
      if(bulk) then
        ibig=ibig+1
        efermi1=xbr(ibig)
      else if(rightm.eq.'V'.and.newvac.eq.1) then
        ibig=ibig+1
        v01=xbr(ibig)
      end if
c
      do li=1,nintfc
        x=dlog(rs(li))-(ns(li)-1)*dx(li)
        do i=1,ns(li)
          r=dexp(x)
          ibig=ibig+1
          bra(i,li)=xbr(ibig)*r
          ibig=ibig+1
          brb(i,li)=xbr(ibig)*r
          x=x+dx(li)
        end do
      end do
c
      if(nspheres.and.ibr.eq.1) then
        write(6,'(/'' WS radii:'')')
        do li=1,nintfc
          write(6,'('' I'', i3,''  Sws_old, Sws_new'',2f17.8)')
     >    li,rs(li),rsnew(li)
          r0(0)=0.d0
          do j=1,ns(li)
            x=dlog(rs(li))-(ns(li)-j)*dx(li)
            r0(j)=dexp(x)
          end do
          vra0(0)=-2.d0*za(li)
          vrb0(0)=-2.d0*zb(li)
          bra0(0)=0.0d0
          brb0(0)=0.0d0
c interpolate potentials
          do j=1,ns(li)
            vra0(j)=vra(j,li)
            vrb0(j)=vrb(j,li)
            bra0(j)=bra(j,li)
            brb0(j)=brb(j,li)
          end do
          do j=1,ns(li)
            x=dlog(rsnew(li))-(ns(li)-j)*dx(li)
            r=dexp(x)
            vra(j,li)=ylag(r,r0(0),vra0(0),0,3,ns(li)+1,iex)
            vrb(j,li)=ylag(r,r0(0),vrb0(0),0,3,ns(li)+1,iex)
            bra(j,li)=ylag(r,r0(0),bra0(0),0,3,ns(li)+1,iex)
            brb(j,li)=ylag(r,r0(0),brb0(0),0,3,ns(li)+1,iex)
          end do
          rs(li)=rsnew(li)
          if(bulk.and.(li.lt.minprc)) then
            rsl(li)=rs(li)
            rsr(li)=rs(li)
          end if
        end do
      end if
        if(dabs(vrsbulk).gt.tol) then    
c muffin-tin zero according to SCA of L. Vitos
          vrsb=0.d0
          rs2=0.0d0
          do li=1,nintfc
            vrsb=vrsb+(conc(li)*vra(ns(li),li)+
     >      (1.0d0-conc(li))*vrb(ns(li),li))*rs(li)
            rs2=rs2+rs(li)*rs(li)
          end do
          vrsb=vrsb/rs2
          vrsb=vrsb-vrsh
c shift potentials by vrsbulk
          do li=1,nintfc
            r1=rs(li)
            n1=ns(li)
            dx1=dx(li)
            x0=dlog(r1)-(n1-1)*dx1
            x=x0
            do i=1,n1
              r=dexp(x)
              vra(i,li)=vra(i,li)-r*vrsb
              vrb(i,li)=vrb(i,li)-r*vrsb
              x=x+dx1
            end do
          end do
c         vrsbulk=vrsbulk+vrsb
          vrsbulk=vrsb
        end if
c
c     do li=1,nintfc
c        write(6,'('' Layer'',i3)') li
c        write(6,'(4d20.10)') (vra(j,li),j=1,ns(li))
c        write(6,'(4d20.10)') (bra(j,li),j=1,ns(li))
c     end do
c      
c -Overwrite vacuum potential if required    
c
      if(bulk) then
        write(6,'(/'' I'',i4,''  VRS'',f15.10,''  EF0'',f15.10,
     >  3x,'' DEF'',f17.13,3x,'' EF'',f17.13)')
     >  itscfcont,vrsbulk,efermi,defermi,efermi1
        efermi=efermi1
      else if(rightm.eq.'V') then
        write(6,'(/'' Vacuum potential level:'',3d20.10/
     >             '' Work function:'',d20.10,'' eV''/)')
     >  gv0,v0,v01,(v01-efermi)*13.606d0
c       dv0=dabs((v0-v01)/v0)
c       if(newvac.eq.1.and.dv0.lt.tolv0) v0=v01
        if(newvac.eq.1) v0=v01
      end if
      call flush(6)
c
      return
      end
