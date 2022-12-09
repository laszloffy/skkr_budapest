c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
c
c  WARNING !!!!! Revise usage of cpamatinl,r
c
      subroutine tmatini(
     > nintfc,ce,lmax,rel,v0,bulk,iscreen,vscreen,
     > concl,concr,conc,cpamatin,cpamatinl,cpamatinr,
     > idpotla,vrla,idpotlb,vrlb,dxl,nsl,rsl,
     > idpotra,vrra,idpotrb,vrrb,dxr,nsr,rsr,
     > idpota,vra,idpotb,vrb,dx,ns,rs,
     > tminvl,tminvr,tminv,tminva,tminvb)
c
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
      parameter(mdim=mdimnr)
c
      logical bulk,rel
      logical cpalay,cpamatin,cpamatinl,cpamatinr
c 
      character*10 idpota(mintfc),idpotla(minprc),idpotra(minprc)
      character*10 idpotb(mintfc),idpotlb(minprc),idpotrb(minprc)
c
      dimension conc(mintfc),concl(minprc),concr(minprc)
c
      dimension vra(nrad,mintfc),vrb(nrad,mintfc)
      dimension dx(mintfc),ns(mintfc),rs(mintfc)
      dimension vrla(nrad,minprc),vrlb(nrad,minprc)
      dimension dxl(minprc),nsl(minprc),rsl(minprc)
      dimension vrra(nrad,minprc),vrrb(nrad,minprc)
      dimension dxr(minprc),nsr(minprc),rsr(minprc)
c
      complex*16 tminv(lmmaxp,lmmaxp,mintfc)
      complex*16 tminva(lmmaxp,lmmaxp,mintfc)
      complex*16 tminvb(lmmaxp,lmmaxp,mintfc)
      complex*16 tminvl(lmmaxp,lmmaxp,minprc),
     &           tminvr(lmmaxp,lmmaxp,minprc)
c
      complex*16 ce,alphalkkr,alpharkkr,alphaintkkr
c
      common/lay2d/cvec(mtotal,3),nextra,nbulkl,nbulkr,
     &             nprc,ninprc(0:mprc+1)
      common/scrpar/alphalkkr(0:lmaxp,minprc),alpharkkr(0:lmaxp,minprc),
     &              alphaintkkr(0:lmaxp,mintfc)
c
      data tiny/1.0d-6/
c
c ********************
c initialize constants
c ********************
c 
      nl=lmax+1
      nl2=nl*nl
      l2=2*lmax
      ninprcl=ninprc(0)
      ninprcr=ninprc(nprc+1)
c
c ********************
c screening parameters
c ********************
c
      call czero(alphalkkr,minprc*(lmaxp+1))      
      call czero(alpharkkr,minprc*(lmaxp+1))      
      call czero(alphaintkkr,mintfc*(lmaxp+1))
      if (iscreen.ge.1) then
c Use square well t-matrix
c        -------------------------------------------------------------
         call alphamat(ce,iscreen,vscreen,lmax,nintfc,ninprcl,ninprcr)
c        -------------------------------------------------------------
      end if
c
c *******************************
c t-matrix for bulk and interface
c *******************************
c
      do li=1,ninprcl
        cpalay=(1.0d0-dabs(concl(li))).gt.tiny
        if(.not.cpalay) then
c         ----------------------------------------------
          call tmat(ce,lmax,idpotla(li),v0,vrla(1,li),
     >              dxl(li),nsl(li),rsl(li),rel,
     >              tminvl(1,1,li),alphalkkr(0,li))
c         ----------------------------------------------
        else
          if(.not.bulk.and..not.cpamatinl) then
            write(6,'(/'' Check cpamatinl !!!'')')
            stop
          end if
        end if
      enddo
c
      do li=1,ninprcr
        cpalay=(1.0d0-dabs(concr(li))).gt.tiny
        if(.not.cpalay) then
c         ----------------------------------------------
          call tmat(ce,lmax,idpotra(li),v0,vrra(1,li),
     >              dxr(li),nsr(li),rsr(li),rel,
     >              tminvr(1,1,li),alpharkkr(0,li))
c         ----------------------------------------------
        else
          if(.not.bulk.and..not.cpamatinr) then
            write(6,'(/'' Check cpamatinr !!!'')')
            stop
          end if
        endif
      enddo              
c
      do li=1,nintfc
        cpalay=(1.d0-conc(li)).gt.tiny
        call tmat(ce,lmax,idpota(li),v0,vra(1,li),dx(li),
     >            ns(li),rs(li),rel,tminva(1,1,li),alphaintkkr(0,li))
        call tmat(ce,lmax,idpotb(li),v0,vrb(1,li),dx(li),
     >            ns(li),rs(li),rel,tminvb(1,1,li),alphaintkkr(0,li))
        if(.not.(cpamatin.and.cpalay))
     >  call tata(conc(li),nl2,tminv(1,1,li),tminva(1,1,li),
     >            tminvb(1,1,li))
      end do
c
      return
      end
