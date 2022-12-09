c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
c***********************************************************************
c CHECK DECAY OF STRUCTURE CONSTANTS FOR FIRST ITERATION (and maybe, for
c some iterations further) ==> REVISE PRINCIPAL LAYER SIZE
c***********************************************************************
      subroutine gstore(park,psq,factor,bulkgeo,side,
     >                  lmax,nintfc_true,g2d,ihole,c_light,mmprc1)
c
c -After generating the 2d lattice, calculate 2d screened 
c  structure constants in k-space using Ewalds method (subroutine kamstr
c  which has been programmed by Arthur Ernst using the papers of Kambe).
c
c  The layer indices (j,i) have to be reversed to (i,j) to correspond
c  with Laszlo's 2d structure constants!!.
c
c  mnpnp: short for M(n+1,n+1)
c  mnnp:  short for M(n,n+1)
c  mnpn:  short for M(n+1,n)
c
c  Extended for holes --> If ihole=true then transform the structure constant
c
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
      integer lgaunt
      parameter (lgaunt=(2*lmaxp+1)**2)
c
      logical bulkgeo
      logical ihole
c
      character*1 side
c
      integer ipl(mnum)
      integer ninprc(0:mprc+1)
c
      real*8 park(2),sep(0:3)
      real*8 gau(lmmaxp*lmmaxp*lgaunt),rho(2)
      real*8 cvec(mtotal,3),cbulk(3)
c
      complex*16 psq
      complex*16 g2d(lmmaxp,lmmaxp,-minprc:minprc,minprc,0:mmprc1+1)
      complex*16 gtotal(lmmaxp,lmmaxp,mtotal1,mtotal1)
      complex*16 gtotalp(lmmaxp,lmmaxp,mtotal1,mtotal1)
c
      common/test/itest
      common/lay2d/cvec_true(mtotal,3),nextra,nbulkl,nbulkr,
     &             nprc_true,ninprc_true(0:mprc+1)
c
      common/dirvec2d/rvec(2,mdir),distr,maxr
      common/vindr2d/indr(mdir),nshr(mr),maxshr
c
      common/recvec2d/veck(2,mrec),distk,maxk
      common/vindk2d/indk(mrec),nshk(mk),maxshk
      common/brill2d/bx(2),by(2)
      data tol/1.0D-10/
c
c---> c in rydberg units:
c      c=274.072d0 -- from input
      c=c_light
c
      nl2=(lmax+1)*(lmax+1)
      kmymax=2*nl2
c       
      rlvax=bx(1)
      rlvay=bx(2)
      rlvbx=by(1)
      rlvby=by(2)
      rlmin=min(rlvax*rlvax+rlvay*rlvay,rlvbx*rlvbx+rlvby*rlvby)
      rlmin=sqrt(rlmin)
c
c - calculate gaunt coefficients for routine kamstr
c
c     ------------------------
      call gauntcalc(lmax,gau)
c     ------------------------
c
c sort reciprocal space vectors around (park(1),park(2))
c     call vecsrt2d(maxk,mk,veck,park,distk,indk,nshk,maxshk)
c
c --The outermost bulk principal layers will also be involved in 
c   the screening procedure (-> nprc+2+2*nextra; nextra determined in
c   routine readini.f)
c
      if(.not.bulkgeo) then
          nprc = nprc_true
          do iprc=0,nprc+1
            ninprc(iprc)=ninprc_true(iprc)
          end do
          nintfc = nintfc_true
          ninprcl = ninprc(0)
          ninprcr = ninprc(nprc+1)
          num  = nintfc + ninprcl + ninprcr
          ntotal  = num + nextra*(ninprcl+ninprcr)
          do i=1,ntotal
            cvec(i,1)=cvec_true(i,1)
            cvec(i,2)=cvec_true(i,2)
            cvec(i,3)=cvec_true(i,3)
          end do
      end if
c
      if(bulkgeo.and.(side.eq.'L')) then
          ninprcl = ninprc_true(0)
          ninprcr = ninprcl
          nprc=1
          do iprc=0,nprc+1
            ninprc(iprc)=ninprcl
          end do
          nintfc = ninprcl         
          ntotal  = (2*nextra+3)*ninprcl
          num  = 3*ninprcl
          cbulk(1)=cvec_true(ninprcl+1,1)-cvec_true(1,1)
          cbulk(2)=cvec_true(ninprcl+1,2)-cvec_true(1,2)
          cbulk(3)=cvec_true(ninprcl+1,3)-cvec_true(1,3)
          do i=1,ninprcl
            cvec(i,1)=cvec_true(i,1)-cvec_true(1,1)
            cvec(i,2)=cvec_true(i,2)-cvec_true(1,2)
            cvec(i,3)=cvec_true(i,3)-cvec_true(1,3)
          end do
          do iprc=2,2*nextra+3
            do inprc=1,ninprcl
              i=(iprc-1)*ninprcl+inprc
              cvec(i,1)=cvec(inprc,1)+(iprc-1)*cbulk(1)
              cvec(i,2)=cvec(inprc,2)+(iprc-1)*cbulk(2)
              cvec(i,3)=cvec(inprc,3)+(iprc-1)*cbulk(3)
            end do
          end do
c         do i=1,ntotal
c           write(6,'(3d20.10)') cvec(i,1),cvec(i,2),cvec(i,3)
c         end do
      end if
c
      if(bulkgeo.and.(side.eq.'R')) then
          ninprcl = ninprc_true(nprc_true+1)
          ninprcr = ninprcl
          nprc=1
          do iprc=0,nprc+1
            ninprc(iprc)=ninprcl
          end do
          nintfc = ninprcl         
          ntotal  = (2*nextra+3)*ninprcl
          num  = 3*ninprcl
          minind=(nextra+1)*ninprc_true(0)+nintfc_true
          cbulk(1)=cvec_true(minind+ninprcl+1,1)-cvec_true(minind+1,1)
          cbulk(2)=cvec_true(minind+ninprcl+1,2)-cvec_true(minind+1,2)
          cbulk(3)=cvec_true(minind+ninprcl+1,3)-cvec_true(minind+1,3)
          do i=1,ninprcl
            cvec(i,1)=cvec_true(minind+i,1)-cvec_true(minind+1,1)
            cvec(i,2)=cvec_true(minind+i,2)-cvec_true(minind+1,2)
            cvec(i,3)=cvec_true(minind+i,3)-cvec_true(minind+1,3)
          end do
          do iprc=2,2*nextra+3
            do inprc=1,ninprcl
              i=(iprc-1)*ninprcl+inprc
              cvec(i,1)=cvec(inprc,1)+(iprc-1)*cbulk(1)
              cvec(i,2)=cvec(inprc,2)+(iprc-1)*cbulk(2)
              cvec(i,3)=cvec(inprc,3)+(iprc-1)*cbulk(3)
            end do
          end do
c         do i=1,ntotal
c           write(6,'(3d20.10)') cvec(i,1),cvec(i,2),cvec(i,3)
c         end do
      end if
c
c -calculate 2D screened, layer indexed structure constants
c 
c     write(6,*) ' istep=',istep
      do i=1,ntotal
         do j=1,ntotal
            sep(1)=cvec(j,1)-cvec(i,1) 
            sep(2)=cvec(j,2)-cvec(i,2) 
            sep(3)=cvec(j,3)-cvec(i,3) 
            sep(0)=dsqrt(sep(1)**2+sep(2)**2+sep(3)**2)
c           write(6,*) 'Layers:',i,j
c           write(6,*)'sep:',(sep(k),k=1,3)
c
c - determine whether to use the series expansion or calculate the
c   structure constants directly depending if the spheres overlap or
c   not.
c
            if (abs(sep(3))*rlmin.ge.5.0d0) then
c                   write(6,*) ' Direct method'
              if (j.gt.i) then
c               -------------------------------------------------
                call directg(park,sep(1),psq,lmax,
     >                       gtotal(1,1,j,i),gtotal(1,1,i,j))
c               -------------------------------------------------
              end if
            else
c                   write(6,*) ' Kamstr'
c
c sort real space vectors around (sep(1),sep(2))
c
              rho(1)=sep(1)
              rho(2)=sep(2)
              call vecsrt2d(maxr,mr,rvec,rho,distr,indr,nshr,maxshr)
c
c             --------------------------------------------------
              call kamstr(gtotal(1,1,j,i),psq,park(1),park(2),
     >                    lmax,lmmaxp,maxshk,nshk,indk,veck,
     >                    maxshr,nshr,indr,rvec,sep,gau,factor)
c             --------------------------------------------------
            end if
c
         end do
      end do
c
      if(ihole) then
         gtotalp=gtotal
         do i=1,mtotal1
         do j=1,mtotal1
            gtotal(:,:,i,j)=-conjg(transpose(gtotalp(:,:,j,i)))
            call lmtolmconj(gtotal(:,:,i,j))
         end do
         end do
      end if
c
c     write(6,'('' Unscreened structure constants'')')
c     do i=1,ntotal
c        do j=1,ntotal
c           write(6,'('' layersU:'',2i3)') i,j
c           call outmat1(gtotal(1,1,i,j),nl2,nl2,lmmaxp,tol,6)
c        end do
c     end do
c
c - do the screening transformation / ihole separates electron and hole parts
c
c     ------------------------------------------------------------------
      call scrstr(gtotal,nl2,ntotal,nextra,ninprcl,ninprcr,bulkgeo,
     >            side,ihole)
c     -----------------------------------------------------------------
c
c     write(6,'('' Screened structure constants'')')
c     do i=1,ntotal
c        do j=1,ntotal
c           write(6,'('' layersS:'',2i3)') i,j
c           call outmat1(gtotal(1,1,i,j),nl2,nl2,lmmaxp,tol,6)
c        end do
c     end do
c
c--Drop the outermost (nextra)*ninprc layers---
c
c     Assign to each layer its principal layer
c
      do il=1,num
        if(il.le.ninprcl) then
          ipl(il) = 0
          icpl = 1
          icount = 0
        elseif(il.gt.(num-ninprcr)) then
          ipl(il) = nprc + 1
        else
          icount = icount + 1
          ipl(il)= icpl
          if(icount.eq.ninprc(icpl)) then
            icpl = icpl+1
            icount = 0
          endif
        endif
      enddo
c
c     Store only the tridiagonal blocks of g2d
c
      imin = nextra*ninprcl
      icount = 0
      do il=1,num  
         i= imin + il
         icpl = ipl(il)
c        if(istep.eq.2) 
c    >   write(6,'(''i='',i3,'' il='',i3,'' pl='',i3)') i,il,icpl
         icount = icount+1
         if(il.le.ninprcl) then 
           jmin = -ninprcl
           jmax = max(ninprcl,ninprc(icpl+1))
         elseif(il.gt.(num-ninprcr)) then
           jmin = -max(ninprcr,ninprc(icpl-1))
           jmax = ninprcr
         else
           if(ninprc(icpl).le.ninprc(icpl-1)) then
             jmin = -ninprc(icpl-1)
           else
             jmin = -min(ninprc(icpl),ninprc(icpl-1)+icount-1)
           endif
           if(ninprc(icpl).le.ninprc(icpl+1)) then
             jmax = ninprc(icpl+1)
           else
             jmax = min(ninprc(icpl),ninprc(icpl+1)+ninprc(icpl)-icount)
           endif
         endif
         do 20 jj=jmin,jmax
           j = i+jj      
           do 20 n=1,nl2
           do 20 m=1,nl2
              g2d(m,n,jj,icount,icpl)=gtotal(m,n,i,j)
  20     continue
         if(icount.eq.ninprc(icpl)) icount=0
      enddo
c
c
      return
      end
