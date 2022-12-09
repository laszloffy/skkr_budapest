c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine vgen(
     > itscf,itscfcont,nintfc,lmax,sigma,rightm,bulk,
     > orbpol,opot,conc,concl,concr,
     > za,qca,qva,qmoma,qmomla,qmomra,
     > zb,qcb,qvb,qmomb,qmomlb,qmomrb,
     > rhoca,rhova,rhospa,rhodspa,rhocb,rhovb,rhospb,rhodspb,dx,ns,rs,
     > efermi,defermi,vra,vrb,bra,brb,bopra,boprb,lza,lzb,
     > v0,dv0,ferr1,newvac,
     > enpota,enela,enxca,enpotb,enelb,enxcb,enmaga,enmagb)
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      parameter (mbig=8*nrad*mintfc+2)
c
      logical bulk,cpalay,orbpol,opot
      character*1 rightm
c
      real*8 lza(2,mintfc),lzb(2,mintfc)
      dimension dx(mintfc),rs(mintfc),ns(mintfc)
      dimension vra(nrad,mintfc),bra(nrad,mintfc),bopra(nrad,2,mintfc)
      dimension vrb(nrad,mintfc),brb(nrad,mintfc),boprb(nrad,2,mintfc)
      dimension gvra(nrad,mintfc),gvrb(nrad,mintfc)
      dimension gbra(nrad,mintfc),gbrb(nrad,mintfc)
      dimension gbopra(nrad,2,mintfc),gboprb(nrad,2,mintfc)
      dimension rhoca(nrad,mintfc),rhova(nrad,mintfc)
      dimension rhocb(nrad,mintfc),rhovb(nrad,mintfc)
      dimension rhospa(nrad,2,mintfc),rhospb(nrad,2,mintfc)
      dimension rhodspa(nrad,2,mintfc),rhodspb(nrad,2,mintfc)
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
      common/broypot/cmix,wbr,nbrmax(100)
c
      data tolv0/0.1d0/
      data anul/1.0d-15/,tiny/1.0d-6/
c
      save ibr,ibrmax,ibrcase
c
c     ----------------------------------------------------------------
      call newpot(itscfcont,nintfc,lmax,sigma,rightm,bulk,
     >            orbpol,conc,concl,concr,
     >            za,qca,qva,qmoma,qmomla,qmomra,lza,
     >            zb,qcb,qvb,qmomb,qmomlb,qmomrb,lzb,
     >            rhoca,rhova,rhospa,rhodspa,rhocb,rhovb,rhospb,rhodspb,
     >            dx,ns,rs,gv0,gvra,gvrb,gbra,gbrb,gbopra,gboprb,
     >            enpota,enela,enxca,enpotb,enelb,enxcb,enmaga,enmagb)
c     ----------------------------------------------------------------
c
      if((itscf.eq.1).and.(orbpol)) then
        do li=1,nintfc
          do i=1,ns(li) 
             bopra(i,1,li)=gbopra(i,1,li)
             bopra(i,2,li)=gbopra(i,2,li)
             boprb(i,1,li)=gboprb(i,1,li)
             boprb(i,2,li)=gboprb(i,2,li)
          end do
        end do
      end if
c
c - Potential mixing
c
      ibig=0
      do li=1,nintfc
        x=dlog(rs(li))-(ns(li)-1)*dx(li)
        do i=1,ns(li)
          r=dexp(x)
          ibig=ibig+1
          xbr(ibig)=vra(i,li)
          fbr(ibig)=gvra(i,li)-vra(i,li)
c         xbr(ibig)=(vra(i,li)+2.d0*za(li))/r
c         fbr(ibig)=(gvra(i,li)-vra(i,li))/r
          ibig=ibig+1
          xbr(ibig)=vrb(i,li)
          fbr(ibig)=gvrb(i,li)-vrb(i,li)
c         xbr(ibig)=(vrb(i,li)+2.d0*zb(li))/r
c         fbr(ibig)=(gvrb(i,li)-vrb(i,li))/r
          x=x+dx(li)
        end do
      end do
c
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
          xbr(ibig)=bra(i,li)
          fbr(ibig)=gbra(i,li)-bra(i,li)
c         xbr(ibig)=bra(i,li)/r
c         fbr(ibig)=(gbra(i,li)-bra(i,li))/r
          ibig=ibig+1
          xbr(ibig)=brb(i,li)
          fbr(ibig)=gbrb(i,li)-brb(i,li)
c         xbr(ibig)=brb(i,li)/r
c         fbr(ibig)=(gbrb(i,li)-brb(i,li))/r
          x=x+dx(li)
        end do
      end do
      nbig2=ibig
c
      do li=1,nintfc
        x=dlog(rs(li))-(ns(li)-1)*dx(li)
        do i=1,ns(li)
          r=dexp(x)
          ibig=ibig+1
          xbr(ibig)=bopra(i,1,li)
          fbr(ibig)=gbopra(i,1,li)-bopra(i,1,li)
          ibig=ibig+1
          xbr(ibig)=bopra(i,2,li)
          fbr(ibig)=gbopra(i,2,li)-bopra(i,2,li)
          ibig=ibig+1
          xbr(ibig)=boprb(i,1,li)
          fbr(ibig)=gboprb(i,1,li)-boprb(i,1,li)
          ibig=ibig+1
          xbr(ibig)=boprb(i,2,li)
          fbr(ibig)=gboprb(i,2,li)-boprb(i,2,li)
          x=x+dx(li)
        end do
      end do
      nbig3=ibig
c
      nbig=ibig
c
      call dott(xbr,xbr,err1,nbig1)
      call dott(fbr,fbr,t1,nbig1)
      ferr1=dsqrt(t1/err1)
      write(6,'(/'' Error for V(r) :'',3d17.6)')  t1,err1,ferr1
c
      call dott(xbr(nbig1+1),xbr(nbig1+1),err2,nbig2-nbig1)
      call dott(fbr(nbig1+1),fbr(nbig1+1),t1,nbig2-nbig1)
      if(err2.lt.anul) then
        ferr2=0.d0
      else
        ferr2=dsqrt(t1/err2)
      end if
      write(6,'('' Error for Bsp(r) :'',3d17.6)')  t1,err2,ferr2
c
      call dott(xbr(nbig2+1),xbr(nbig2+1),err2,nbig3-nbig2)
      call dott(fbr(nbig2+1),fbr(nbig2+1),t1,nbig3-nbig2)
      if(err2.lt.anul) then
        ferr2=0.d0
      else
        ferr2=dsqrt(t1/err2)
      end if
      write(6,'('' Error for Bop(r) :'',3d17.6)')  t1,err2,ferr2
c
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
        x=dlog(rs(li))-(ns(li)-1)*dx(li)
        do i=1,ns(li)
          r=dexp(x)
          ibig=ibig+1
          vra(i,li)=xbr(ibig)
c         vra(i,li)=xbr(ibig)*r-2.d0*za(li)
          ibig=ibig+1
          vrb(i,li)=xbr(ibig)
c         vrb(i,li)=xbr(ibig)*r-2.d0*zb(li)
          x=x+dx(li)
        end do
      end do
c
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
          bra(i,li)=xbr(ibig)
c         bra(i,li)=xbr(ibig)*r
          ibig=ibig+1
          brb(i,li)=xbr(ibig)
c         brb(i,li)=xbr(ibig)*r
          x=x+dx(li)
        end do
      end do
c
      do li=1,nintfc
        x=dlog(rs(li))-(ns(li)-1)*dx(li)
        do i=1,ns(li)
          r=dexp(x)
          ibig=ibig+1
          bopra(i,1,li)=xbr(ibig)
          ibig=ibig+1
          bopra(i,2,li)=xbr(ibig)
          ibig=ibig+1
          boprb(i,1,li)=xbr(ibig)
          ibig=ibig+1
          boprb(i,2,li)=xbr(ibig)
          x=x+dx(li)
        end do
      end do
c
c -Overwrite vacuum potential if required    
c
      if(bulk) then
        write(6,'(/'' I'',i4,''  EF0'',f17.13,
     >  3x,'' DEF'',f17.13,3x,'' EF'',f17.13)')
     >  itscfcont,efermi,defermi,efermi1
        efermi=efermi1
      else if(rightm.eq.'V') then
        write(6,'(/'' Vacuum potential level:'',3d20.10/
     >             '' Work function:'',d20.10,'' eV''/)')
     >  gv0,v0,v01,(v01-efermi)*13.606d0
c       dv0=dabs((v0-v01)/v0)
c       if(newvac.eq.1.and.dv0.lt.tolv0) v0=v01
        if(newvac.eq.1) v0=v01
      end if
c
      return
      end
