c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
c**********************************************************************
c MAYBE, A FLAG CAN BE ADDED TO SKIP DUMMY M-MATRICES CALCULATIONS IN
c THE CALLS FROM DUGO.
c**********************************************************************
      subroutine packer(ni,nim,nip,g2d,m00,m10,m01,lmax,lkmax,mdim,irel)
c
c Build up the m00,m01,m0-1 matrices corresponding to all interactions
c of a given principal layer ("0").
c
c To be used by TAU2D and DUGO.
c
c-----------------------------------------------------------
c     (nonrel):       mdim=mdimnr,   lkmax=(lmax+1)**2
c     (rel):          mdim=mdimr,    lkmax=kmymax=2*(lmax+1)**2
c     (nonrel-bogo):  mdim=mdimbogo, lkmax=2*(lmax+1)**2
c     (rel-bogo):     mdim=mdimdbogo,lkmax=4*(lmax+1)**2
c-----------------------------------------------------------
c     icpl = current principal layer
c     nim = ninprc(icpl-1)
c     ni  = ninprc(icpl)
c     nip = ninprc(icpl+1)
c-----------------------------------------------------------
c
c
      include '../param.h'
c
      complex*16 g2d(bogomaxp,bogomaxp,-minprc:minprc,minprc)
      complex*16 rg2d(dbogomaxp,dbogomaxp)
      complex*16 m00(mdim,mdim),m01(mdim,mdim),m10(mdim,mdim)
c
      call czero(m01,mdim*mdim)
      call czero(m10,mdim*mdim)
c
      do i=1,ni
        if(ni.le.nim) then
          jmin = -nim
        else
          jmin = -min(ni,nim+i-1)
        endif
        if(ni.le.nip) then
          jmax = nip
        else
          jmax = min(ni,nip+ni-i)
        endif                  
        kki = (i-1)*lkmax
c
c Fill up M0-1
c
        do j=jmin,-i
          jp=nim+i+j
          kkj = (jp-1)*lkmax
          if(irel.eq.1) call relmtrx(g2d(1,1,j,i),rg2d,lmax)
          do k1=1,lkmax
            kk1 = k1+kki
            do k2=1,lkmax
              kk2 = k2+kkj
              if(irel.eq.1) m10(kk1,kk2)= -rg2d(k1,k2)
              if(irel.eq.0) m10(kk1,kk2)= -g2d(k1,k2,j,i)
            end do
          end do
        end do
c
c Fill up M00 
c
        do j=1-i,ni-i
          jp= j+i
          kkj= (jp-1)*lkmax
          if(irel.eq.1) call relmtrx(g2d(1,1,j,i),rg2d,lmax)
          do k1=1,lkmax
            kk1 = k1+kki
            do k2=1,lkmax
              kk2 = k2+kkj
              if(irel.eq.1) m00(kk1,kk2)= -rg2d(k1,k2)
              if(irel.eq.0) m00(kk1,kk2)= -g2d(k1,k2,j,i)
            end do
          end do
        end do
c
c Fill up M10 (M0-1)
c
        do j=ni-i+1,jmax
          jp= j+i-ni
          kkj= (jp-1)*lkmax
          if(irel.eq.1) call relmtrx(g2d(1,1,j,i),rg2d,lmax)
          do k1=1,lkmax
            kk1 = k1+kki
            do k2=1,lkmax
              kk2 = k2+kkj
              if(irel.eq.1) m01(kk1,kk2)= -rg2d(k1,k2)
              if(irel.eq.0) m01(kk1,kk2)= -g2d(k1,k2,j,i)
            end do
          end do
        end do
c
      enddo
c
      return
      end
