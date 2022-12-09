c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine upm(ce,E_Fermi,cl,soc,dx,rs,ns,vr,br,rb,
     >               urp,urm,couplings)
c
c     ***********************************************************
c     * Potential functions up*r and um*r/c^2 needed to solve   *
c     * the radial Dirac equation                               *
c     * Spin-orbit scaling is employed                          *
c     * Effective field points along rb; it's set to z if local *
c     ***********************************************************
c
      implicit none
c
      include '../param.h'
c
c I/O variables
      integer ns,couplings(kmymaxp,kmymaxp)
      real*8 dx,rs,vr(nrad),br(nrad),rb(3),cl
      real*8 soc,E_Fermi
      complex*16 ce,sr
      complex*16 urp(kmymaxp,kmymaxp,nrad),urm(kmymaxp,kmymaxp,nrad)
c
c internal variables
      integer i,j,k,li,ki,ilag,nlag,iex,kmymax
      real*8 x,x0,r,xki,fact,csqinv,sum,tiny
      real*8 xlag(nrad),vrlag(nrad),brlag(nrad)
      real*8 ylag
      complex*16 gc,gbc
c
      complex*16 scoeff(kmymaxp,kmymaxp),sbcoeff(kmymaxp,kmymaxp)
c
      complex*16 sxcoeff(kmymaxp,kmymaxp),sxbcoeff(kmymaxp,kmymaxp)
      complex*16 sycoeff(kmymaxp,kmymaxp),sybcoeff(kmymaxp,kmymaxp)
      complex*16 szcoeff(kmymaxp,kmymaxp),szbcoeff(kmymaxp,kmymaxp)
      common/sigmat/sxcoeff,sxbcoeff,sycoeff,sybcoeff,szcoeff,szbcoeff
c
      integer nout
      data nout/5/
c
      data tiny/1.0d-15/
c
c
      ilag=3
      vrlag(1:ns)=vr(1:ns)
      brlag(1:ns)=br(1:ns)
      x0=dlog(rs)-(ns-1)*dx
      do k=1,ns
        xlag(k)=x0+(k-1)*dx
      end do
      nlag=ns+1
      vrlag(nlag)=0.d0
      brlag(nlag)=0.d0
      xlag(nlag)=x0+(ns+nout-1)*dx
      do k=ns+1,ns+nout
        x=x0+(k-1)*dx
        vr(k)=ylag(x,xlag,vrlag,0,ilag,nlag,iex)
        br(k)=ylag(x,xlag,brlag,0,ilag,nlag,iex)
      end do
c
      urp=(0.0d0,0.0d0)
      urm=(0.0d0,0.0d0)
c

      scoeff=rb(1)*sxcoeff+rb(2)*sycoeff+rb(3)*szcoeff
      sbcoeff=rb(1)*sxbcoeff+rb(2)*sybcoeff+rb(3)*szbcoeff
c
      csqinv=1.0/(cl*cl)
      do i=1,kmymaxp
        li=ldex(i)
        ki=kapdex(i)
        xki=soc*(1.d0+ki)-1.d0
        fact=li*(li+1.d0)-xki*(xki+1.d0) 
        do j=1,kmymaxp
          gc=scoeff(j,i)
          gbc=sbcoeff(j,i)
          if(ldex(i).ne.ldex(j)) then
            gc=0.0d0
            gbc=0.0d0
          end if
          do k=1,ns+nout
            urp(j,i,k)=gc*br(k)
            urm(j,i,k)=-gbc*br(k)*csqinv
          end do
        end do
        x=x0
        do k=1,ns+nout
          r=dexp(x)
          urm(i,i,k)=vr(k)*csqinv+urm(i,i,k)
          sr=r*(1.0d0+ce*csqinv)-urm(i,i,k)
          urp(i,i,k)=vr(k)+urp(i,i,k)+fact/sr
          x=x+dx
        end do
      end do
c
c determine non-zero functions
c OLD VERSION
c      couplings=0
c      do i=1,kmymaxp
c      do j=1,kmymaxp
c        sum=0.d0
c        do k=1,ns+nout 
c          sum=sum+urp(j,i,k)*urp(j,i,k)
c     >           +urm(j,i,k)*urm(j,i,k)
c        end do
c        if(sum.gt.tiny) couplings(j,i)=1
c      end do
c      end do
c
c NEW VERSION, from skkr-ldau
! TODO: kmymax in loops above and here?? (only need lmax)
      couplings=0                           
      do i=1,kmymaxp
      do j=1,kmymaxp                                   
c                                                                                                                                                                                          
c        if (ldex(i).eq.ldex(j)) couplings(i,j)=1
       if ( i .eq. j ) couplings(i,j)=1
       if ( scoeff(i,j) .ne. 0 ) couplings(i,j)=1
       if ( sbcoeff(i,j) .ne. 0 ) couplings(i,j)=1
       if ( ldex(i) .ne. ldex(j) ) couplings(i,j)=0
c       if (couplings(i,j) .ne. 0) write(44,*) i,j
c                                                                                                                                                                                               
      end do                                 
      end do                                           
c
      return
      end
