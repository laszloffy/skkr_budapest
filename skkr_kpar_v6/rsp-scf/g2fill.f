c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine gfill
c=====================
c
c Representations of the spin, orbital and total momentum operators
c on the basis of the spinor spherical harmonics indexed by (kappa,mu)
c
c The representations of the spin-operators are multiplied by two
c
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
      complex*16 sxcoeff(kmymaxp,kmymaxp),sxbcoeff(kmymaxp,kmymaxp)
      complex*16 sycoeff(kmymaxp,kmymaxp),sybcoeff(kmymaxp,kmymaxp)
      complex*16 szcoeff(kmymaxp,kmymaxp),szbcoeff(kmymaxp,kmymaxp)
      common/sigmat/sxcoeff,sxbcoeff,sycoeff,sybcoeff,szcoeff,szbcoeff
c
      complex*16 lxcoeff(kmymaxp,kmymaxp),lxbcoeff(kmymaxp,kmymaxp)
      complex*16 lycoeff(kmymaxp,kmymaxp),lybcoeff(kmymaxp,kmymaxp)
      complex*16 lzcoeff(kmymaxp,kmymaxp),lzbcoeff(kmymaxp,kmymaxp)
      common/lmat/lxcoeff,lxbcoeff,lycoeff,lybcoeff,lzcoeff,lzbcoeff
c
      complex*16 jxcoeff(kmymaxp,kmymaxp),jxbcoeff(kmymaxp,kmymaxp)
      complex*16 jycoeff(kmymaxp,kmymaxp),jybcoeff(kmymaxp,kmymaxp)
      complex*16 jzcoeff(kmymaxp,kmymaxp),jzbcoeff(kmymaxp,kmymaxp)
      common/jmat/jxcoeff,jxbcoeff,jycoeff,jybcoeff,jzcoeff,jzbcoeff
c
      complex*16 sx(72,72),sxb(72,72)
      complex*16 sy(72,72),syb(72,72)
      complex*16 sz(72,72),szb(72,72)
      complex*16 jx(72,72),jxb(72,72)
      complex*16 jy(72,72),jyb(72,72)
      complex*16 jz(72,72),jzb(72,72)
c
      dimension u1(72),ind1(72)
      dimension u2(72),ind2(72)
      common/cgc/u1,u2,ind1,ind2
c
      complex*16 ci
      data ci/(0.0d0,1.0d0)/,small/1.0d-12/
c
      dimension kapdex(72),ldex(72),jdex(72),mydex(72)
      data kapdex/
     * -1,-1,
     *  1, 1,
     * -2,-2,-2,-2,
     *  2, 2, 2, 2,
     * -3,-3,-3,-3,-3,-3,
     *  3, 3, 3, 3, 3, 3,
     * -4,-4,-4,-4,-4,-4,-4,-4,
     *  4, 4, 4, 4, 4, 4, 4, 4,
     * -5,-5,-5,-5,-5,-5,-5,-5,-5,-5,
     *  5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
     * -6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6/
      data ldex/
     *  0, 0,
     *  1, 1,
     *  1, 1, 1, 1,
     *  2, 2, 2, 2,
     *  2, 2, 2, 2, 2, 2,
     *  3, 3, 3, 3, 3, 3,
     *  3, 3, 3, 3, 3, 3, 3, 3,
     *  4, 4, 4, 4, 4, 4, 4, 4,
     *  4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
     *  5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
     *  5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5/
      data jdex/
     *  1, 1,
     *  1, 1,
     *  3, 3, 3, 3, 
     *  3, 3, 3, 3,
     *  5, 5, 5, 5, 5, 5,
     *  5, 5, 5, 5, 5, 5,
     *  7, 7, 7, 7, 7, 7, 7, 7,
     *  7, 7, 7, 7, 7, 7, 7, 7,
     *  9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
     *  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 
     * 11,11,11,11,11,11,11,11,11,11,11,11/
      data mydex/
     * -1, 1,
     * -1, 1,
     * -3,-1, 1, 3,
     * -3,-1, 1, 3,
     * -5,-3,-1, 1, 3, 5,
     * -5,-3,-1, 1, 3, 5,
     * -7,-5,-3,-1, 1, 3, 5, 7,
     * -7,-5,-3,-1, 1, 3, 5, 7,
     * -9,-7,-5,-3,-1, 1, 3, 5, 7, 9,
     * -9,-7,-5,-3,-1, 1, 3, 5, 7, 9,
     *-11,-9,-7,-5,-3,-1, 1, 3, 5, 7, 9,11/
c
      call czero(sx,72*72)
      call czero(sy,72*72)
      call czero(sz,72*72)
      call czero(sxb,72*72)
      call czero(syb,72*72)
      call czero(szb,72*72)
      call czero(jx,72*72)
      call czero(jy,72*72)
      call czero(jz,72*72)
      call czero(jxb,72*72)
      call czero(jyb,72*72)
      call czero(jzb,72*72)
c
      do 10 kmy1=1,72
        c1d=u1(kmy1)
        c1u=u2(kmy1)
        kap1=kapdex(kmy1)
        l1=ldex(kmy1)
        j1=jdex(kmy1)
        my1=mydex(kmy1)
      do 10 kmy2=1,72
        c2d=u1(kmy2)
        c2u=u2(kmy2)
        kap2=kapdex(kmy2)
        l2=ldex(kmy2)
        j2=jdex(kmy2)
        my2=mydex(kmy2)
        x=j2*(j2+2.0d0)-my2*(my2+2.0d0)
        y=j2*(j2+2.0d0)-my2*(my2-2.0d0)
        if(x.lt.small) then
          gamp=0.d0
        else
          gamp=0.25d0*dsqrt(x)
        end if
        if(y.lt.small) then
          gamm=0.d0
        else
          gamm=0.25*dsqrt(y)
        end if
c
        if(l1.ne.l2) goto 10
c
c 'x' and 'y' matrices
c
        if(my1.eq.my2+2) then
          sx(kmy1,kmy2)=c1u*c2d
          sy(kmy1,kmy2)=-ci*c1u*c2d
          if(j1.eq.j2) then
          jx(kmy1,kmy2)=gamp
          jy(kmy1,kmy2)=-ci*gamp
          end if
        else if(my1.eq.my2-2) then
          sx(kmy1,kmy2)=c1d*c2u
          sy(kmy1,kmy2)=ci*c1d*c2u
          if(j1.eq.j2) then
          jx(kmy1,kmy2)=gamm
          jy(kmy1,kmy2)=ci*gamm
          end if
        end if
c
c 'z' matrices
c
        if(my1.eq.my2) then
          sz(kmy1,kmy2)=c1u*c2u-c1d*c2d
          if(j1.eq.j2) jz(kmy1,kmy2)=0.5d0*my1
        end if
c
  10  continue 
c 
      do kmy1=1,32
        kap1=kapdex(kmy1)
        j1=jdex(kmy1)
        my1=mydex(kmy1)
        kapp1=-kap1
        kmyp1=2*kapp1*kapp1+kapp1+(my1+1)/2
c       write(*,'(/2x,i3,2x,3i3,5x,2i3/)') kmy1,kap1,j1,my1,kapp1,kmyp1
      do kmy2=1,32
        kap2=kapdex(kmy2)
        j2=jdex(kmy2)
        my2=mydex(kmy2)
        kapp2=-kap2
        kmyp2=2*kapp2*kapp2+kapp2+(my2+1)/2
c       write(*,'(2x,i3,2x,3i3,5x,2i3)') kmy2,kap2,j2,my2,kapp2,kmyp2
c
        sxb(kmy1,kmy2)=sx(kmyp1,kmyp2)
        syb(kmy1,kmy2)=sy(kmyp1,kmyp2)
        szb(kmy1,kmy2)=sz(kmyp1,kmyp2)
        jxb(kmy1,kmy2)=jx(kmyp1,kmyp2)
        jyb(kmy1,kmy2)=jy(kmyp1,kmyp2)
        jzb(kmy1,kmy2)=jz(kmyp1,kmyp2)
c
      end do
      end do
c
      do kmy1=1,kmymaxp
      do kmy2=1,kmymaxp
        sxcoeff(kmy1,kmy2)=sx(kmy1,kmy2)
        sycoeff(kmy1,kmy2)=sy(kmy1,kmy2)
        szcoeff(kmy1,kmy2)=sz(kmy1,kmy2)
        sxbcoeff(kmy1,kmy2)=sxb(kmy1,kmy2)
        sybcoeff(kmy1,kmy2)=syb(kmy1,kmy2)
        szbcoeff(kmy1,kmy2)=szb(kmy1,kmy2)
        jxcoeff(kmy1,kmy2)=jx(kmy1,kmy2)
        jycoeff(kmy1,kmy2)=jy(kmy1,kmy2)
        jzcoeff(kmy1,kmy2)=jz(kmy1,kmy2)
        jxbcoeff(kmy1,kmy2)=jxb(kmy1,kmy2)
        jybcoeff(kmy1,kmy2)=jyb(kmy1,kmy2)
        jzbcoeff(kmy1,kmy2)=jzb(kmy1,kmy2)
        lxcoeff(kmy1,kmy2)=jx(kmy1,kmy2)-0.5d0*sx(kmy1,kmy2)
        lycoeff(kmy1,kmy2)=jy(kmy1,kmy2)-0.5d0*sy(kmy1,kmy2)
        lzcoeff(kmy1,kmy2)=jz(kmy1,kmy2)-0.5d0*sz(kmy1,kmy2)
        lxbcoeff(kmy1,kmy2)=jxb(kmy1,kmy2)-0.5d0*sxb(kmy1,kmy2)
        lybcoeff(kmy1,kmy2)=jyb(kmy1,kmy2)-0.5d0*syb(kmy1,kmy2)
        lzbcoeff(kmy1,kmy2)=jzb(kmy1,kmy2)-0.5d0*szb(kmy1,kmy2)
      end do
      end do
c
      return
      end
