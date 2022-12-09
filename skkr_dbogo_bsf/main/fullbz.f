c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine fullbz(i2dlat,iset,kuning,phik,xk1,yk1,nk1,
     &                                          xk2,yk2,nk2)
c =====================
c
c k points and weights for a two-dimensional lattice
c in the full brillouin zone
c
c input:  amat(0), ng(0)  - matrices and no. of point-group operations
c                           (for actual idgroup (0), and minimum IBZ)
c         iset            - selects actual symmetry for 1 k-point
c         xk1,yk1,nk1     - k points in the minimum IBZ of i2dlat
c         phik            - rotation of iset=0 for some special cases
c
c output: xk2,yk2,nk2     - k points, weights and their total number 
c                           in the full BZ  (for 1 k-point, all symmetry
c                           equivalent points in the BZ are generated)
c =====================
c
      implicit real*8 (a-h,o-z)
      include '../param.h'
      dimension xk1(mkpar),yk1(mkpar)
      dimension xk2(mkpar),yk2(mkpar)
      dimension iskip(20),veck(2)
      common/rrotmat/amat(2,2,melem),a0mat(2,2,melem),ng,ng0
      data t/1.732050807568877d0/
      data tol/1.0d-8/
c
      pi=4.d0*datan(1.d0)
      if(ng0.gt.20) stop 'FULLBZ: Increase dimension of iskip!'
c
c     Generate points symmetry equivalent to the input k-point
c     (according to idgroup).
      if(iset.eq.0) then
        if(dabs(phik-0.d0).gt.tol) then
          veck(1)=xk1(1)
          veck(2)=yk1(1)
          call vecrot(veck,phik)
          xk1(1)=veck(1)
          yk1(1)=veck(2)
        endif   
        ik2=0
        do ig=1,ng0
          ik2=ik2+1
          xk2(ik2)=a0mat(1,1,ig)*xk1(1) +a0mat(1,2,ig)*yk1(1)
          yk2(ik2)=a0mat(2,1,ig)*xk1(1) +a0mat(2,2,ig)*yk1(1)
        enddo
c       Avoid double counting of edge points.
        do ik=1,ik2
          iskip(ik)=0
        enddo
        do ik=1,ik2-1
          do jk=ik+1,ik2
            diff=dsqrt((xk2(ik)-xk2(jk))**2+(yk2(ik)-yk2(jk))**2)
            if((diff-0.d0).lt.tol) iskip(jk)=1
          enddo
        enddo
        nk2=0
        do ik=1,ik2
          if(iskip(ik).eq.0) then
            nk2=nk2+1
            xk2(nk2)=xk2(ik)
            yk2(nk2)=yk2(ik)
          endif
        enddo
c
      else
c     Generate full BZ from minimum IBZ of crystal group 
c     (determined by i2dlat).
c
        ik2=0
        do ik1=1,nk1
          ngeff=ng
c
c         Avoid double counting of edge points.
c         ------------
c         Oblique
          if(i2dlat.eq.1) then
            if(kuning.eq.0) then
              if(dabs(xk1(ik1)-0.d0).lt.tol) ngeff=1
            else
              if(dabs(yk1(ik1)-0.d0).lt.tol) ngeff=1
            endif
c         ------------
c         Rectangular
          elseif(i2dlat.eq.2.or.i2dlat.eq.3) then
            if(dabs(yk1(ik1)-0.d0).lt.tol.
     &      or.dabs(xk1(ik1)-0.d0).lt.tol) ngeff=2
c         ------------
c         Square
          elseif(i2dlat.eq.4) then
            if(dabs(xk1(ik1)-0.d0).lt.tol.
     &      or.dabs(xk1(ik1)-yk1(ik1)).lt.tol) ngeff=4
c         ------------
c         Hexagonal
          else
            r=dsqrt(xk1(ik1)**2+yk1(ik1)**2)
            phi=dasin(yk1(ik1)/r)
            if(dabs(phi-0.d0).lt.tol.
     &      or.dabs(phi-pi/6.d0).lt.tol) ngeff=6
          endif
c         ------------
c
          do ig=1,ngeff
            ik2=ik2+1
            xk2(ik2)=amat(1,1,ig)*xk1(ik1)+amat(1,2,ig)*yk1(ik1)
            yk2(ik2)=amat(2,1,ig)*xk1(ik1)+amat(2,2,ig)*yk1(ik1)
          enddo
        enddo
        nk2=ik2
      endif
c
      return
      end
