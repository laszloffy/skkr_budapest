c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine scrstr(gtotal,nl2,nlayer,nextra,ninprcl,ninprcr,
     >                  bulkgeo,side)
c =====================
c
c     Screening Transformation
c     --------------------------
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      logical bulkgeo
c
      character*1 side
c
      dimension lml(lmmaxp),lmm(lmmaxp)
c
      complex*16 sum,detl
c     complex*16 g1(nl2*nlayer,nl2*nlayer),c1(nl2*nlayer)
c     complex*16 g(nl2*nlayer,nl2*nlayer)
      complex*16 g1(lmmaxp*mtotal1,lmmaxp*mtotal1),c1(lmmaxp*mtotal1)
      complex*16 g(lmmaxp*mtotal1,lmmaxp*mtotal1)
      complex*16 gtotal(lmmaxp,lmmaxp,mtotal1,mtotal1)
      complex*16 alpha,alphalkkr,alpharkkr,alphaintkkr
c
      common/mome/lml,lmm
      common/scrpar/alphalkkr(0:lmaxp,minprc),alpharkkr(0:lmaxp,minprc),
     &              alphaintkkr(0:lmaxp,mintfc)
c
c--write gtotal(l,lp,i,j) to matrix form---
c
      do i=1,nl2*nlayer
        il=(i/nl2)+1
        ilm=mod(i,nl2)
        if(ilm.eq.0) then
           ilm=nl2
           il=il-1
        end if
c
        do j=1,nl2*nlayer
           jl=(j/nl2)+1
           jlm=mod(j,nl2)
           if(jlm.eq.0) then
              jlm=nl2
              jl=jl-1
           end if
c
c--indices for gtotal: see gstore.f---
c
           g(i,j)=gtotal(ilm,jlm,il,jl)
c          write(6,*) "i,j,il,jl,ilm,jlm,g",i,j,il,jl,ilm,jlm,g(i,j)
        end do
      end do
c
c--set up matrix '1-alpha*G'---
c
      nshiftl=(nextra+1)*ninprcl
      nshiftr=(nextra+1)*ninprcr
      do i=1,nl2*nlayer
          il=(i/nl2)+1
          ilm=mod(i,nl2)
          if(ilm.eq.0) then
             ilm=nl2
             il=il-1
          end if
c
          if(.not.bulkgeo) then
          if((il.gt.nshiftl).and.(il.le.(nlayer-nshiftr))) then
            alpha=alphaintkkr(lml(ilm),(il-nshiftl))
          elseif(il.le.nshiftl) then
            ill = mod(il,ninprcl)
            if(ill.eq.0) ill = ninprcl
            alpha=alphalkkr(lml(ilm),ill)
          else
            ilr = mod(il-nlayer+nshiftr,ninprcr)
            if(ilr.eq.0) ilr = ninprcr
            alpha=alpharkkr(lml(ilm),ilr)
          end if
          end if
c
          if(bulkgeo.and.(side.eq.'L')) then
            ill = mod(il,ninprcl)
            if(ill.eq.0) ill = ninprcl
            alpha=alphalkkr(lml(ilm),ill)
          end if
c
          if(bulkgeo.and.(side.eq.'R')) then
            ill = mod(il,ninprcr)
            if(ill.eq.0) ill = ninprcr
            alpha=alpharkkr(lml(ilm),ill)
          end if
c
c         write(6,*) 'alpha (each layer): ',alpha
c
          do j=1,nl2*nlayer 
             g1(i,j)=-alpha*g(i,j)
          end do
          g1(i,i)=(1.d0,0.d0)+g1(i,i)
      end do
c
c--invert matrix '1-alpha*G'---
c
      call scrinv(g1,nl2*nlayer,detl)
c
c--G' = G * (1-alpha*G)**(-1)---
c
      do i=1,nl2*nlayer
        do j=1,nl2*nlayer
          sum=(0.d0,0.d0)
          do k=1,nl2*nlayer
            sum=sum+g(i,k)*g1(k,j)
          end do
          c1(j)=sum
        end do
c
        do j=1,nl2*nlayer
           il=(i/nl2)+1
           ilm=mod(i,nl2)
           if(ilm.eq.0) then
              ilm=nl2
              il=il-1
           end if
c
           jl=(j/nl2)+1
           jlm=mod(j,nl2)
           if(jlm.eq.0) then
              jlm=nl2
              jl=jl-1
           end if
c
          gtotal(ilm,jlm,il,jl)=c1(j)
c         write(6,*) "i,j,il,jl,ilm,jlm,gtotal,c1",i,j,il,jl,ilm,jlm,
c    >                gtotal(ilm,jlm,il,jl),c1(j)
        end do
      end do
c
      return
      end
