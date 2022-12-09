c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine getg2d_io(g2din,g2dout,bulkgeo,lmax,nintfc_true,
     >                     iwrite,kg2d)
c
c---------------------------------------
c     g2dout=g2din
c     iwrite =-1 -> read  g2din from file
c     iwrite = 1 -> write g2dout in file
c     else no i/o operation occurs 
c---------------------------------------
c
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
      logical bulkgeo
      dimension ipl(mnum)
      integer kg2d
      complex*16
     & g2din(lmmaxp,lmmaxp,-minprc:minprc,minprc,0:mprc1+1),
     & g2dout(lmmaxp,lmmaxp,-minprc:minprc,minprc,0:mprc1+1)
c
      common/lay2d/cvec(mtotal,3),nextra,nbulkl,nbulkr,
     &             nprc_true,ninprc(0:mprc+1)  
c
      nintfc=nintfc_true
      nprc=nprc_true
      if(bulkgeo) then
        nintfc=ninprc(0)
        nprc=1
      end if
c
      nl2=(lmax+1)*(lmax+1)
      num= nintfc + ninprc(0) + ninprc(nprc+1)
c
      if(iwrite.eq.-1) then
        read(kg2d) (((((g2din(k1,k2,k3,k4,k5),k1=1,lmmaxp),k2=1,lmmaxp),
     >             k3=-minprc,minprc),k4=1,minprc),k5=0,mprc1+1)
      end if
c
      do il=1,num
        if(il.le.ninprc(0)) then
          ipl(il) = 0
          icpl = 1
          icount = 0
        elseif(il.gt.(num-ninprc(nprc+1))) then
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
      imin = (nextra)*ninprc(0)
      icount = 0
      do il=1,num
         i= imin + il
         icpl = ipl(il)
         icount = icount+1
         if(il.le.ninprc(0)) then
           jmin = -ninprc(0)
           jmax = max(ninprc(0),ninprc(icpl+1))
         elseif(il.gt.(num-ninprc(nprc+1))) then
           jmin = -max(ninprc(nprc+1),ninprc(icpl-1))
           jmax = ninprc(nprc+1)
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
              g2dout(m,n,jj,icount,icpl)=g2din(m,n,jj,icount,icpl)
  20     continue
         if(icount.eq.ninprc(icpl)) icount=0
      enddo                                           
c
      if(iwrite.eq.1) then
       write(kg2d) 
     >     (((((g2dout(k1,k2,k3,k4,k5),k1=1,lmmaxp),k2=1,lmmaxp),
     >     k3=-minprc,minprc),k4=1,minprc),k5=0,mprc1+1)
      end if
c
      return
      end
