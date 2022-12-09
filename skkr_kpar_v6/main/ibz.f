c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine ibz(i2dlat,idgroup,iogroup,kuning,xk1,yk1,nk1,
     &                                             xk2,yk2,nk2)
c
c  Generate actual IBZ applying point-group symmetry operations to
c  the minimum IBZ.
c
c input:  amat, ng        - matrices and no. of point-group operations 
c         i2dlat, idgroup - translational lattice and point-group
c         iogroup         - original idgroup for some special cases
c         xk1,yk1,nk1     - k points in the minimum IBZ of ilat
c
c output: xk2,yk2,nk2     - k points, weights and their total number
c                           in the actual IBZ of idgroup (iogroup)
c
      implicit real*8 (a-h,o-z)
      include '../param.h'
      character*4 idgroup,iogroup
      dimension xk1(*),yk1(*)
      dimension xk2(*),yk2(*)
      dimension iskip(20) 
      common/rrotmat/amat(2,2,melem),a0mat(2,2,melem),ng,ng0
c
      data tol/1.0d-8/
      pi=4.d0*datan(1.d0) 
      if(ng.gt.20) stop 'IBZ: Increase dimension of iskip!'
c
c     write(6,'('' nk   ='',i3)') nk1
c     do ik=1,nk1
c       write(6,'(2f10.5)') xk1(ik),yk1(ik)
c     end do
c     call flush(6)
c
      ik2=0
      do ik1=1,nk1
        ngeff=ng
c       write(6,'(i3,2x,a)') i2dlat,idgroup
c       write(6,'(2i3,2x,a)') ng,ng0
        do ig=1,ngeff
          iskip(ig)=0
        enddo
c
c----------------------------------------------------------
c For each crystal lattice, select the symmetry operations 
c to generate the actual IBZ from the minimum one. 
c
c Avoid double counting of edge points (ngeff, iskip).
c----------------------------------------------------------
c
c       =========
c        Oblique
c       =========
c
        if(i2dlat.eq.1) then
          if(idgroup.eq.'C2') ngeff=1
          if(kuning.eq.0)then
            if(dabs(xk1(ik1)-0.d0).lt.tol) ngeff=1
          else
            if(dabs(yk1(ik1)-0.d0).lt.tol) ngeff=1
          endif
c
c       =============
c        Rectangular
c       =============
c
        elseif(i2dlat.eq.2.or.i2dlat.eq.3) then
          if(idgroup.eq.'C2v') then
            ngeff=1
          elseif(idgroup.eq.'Cs'.or.idgroup.eq.'C2') then
            if(iogroup.eq.'Csy') then
              if(dabs(yk1(ik1)-0.d0).lt.tol) then
                ngeff=1
              else
                iskip(2) = 1
                ngeff=3
              endif
            else
              if(dabs(xk1(ik1)-0.d0).lt.tol) then
                ngeff=1
              else
                iskip(2) = 1
                iskip(3) = 1
              endif
            endif
          else
            if(dabs(xk1(ik1)-0.d0).lt.tol.or.
     &         dabs(yk1(ik1)-0.d0).lt.tol) ngeff=2
          endif
c
c       ========
c        Square
c       ========
c
        elseif(i2dlat.eq.4) then
          if(idgroup.eq.'C4v') then
            ngeff=1
          elseif(idgroup.eq.'C2v'.or.idgroup.eq.'C4') then
            if(iogroup.eq.'C2vd') then
              ngeff=6
              do ig=2,5
                iskip(ig)=1
              enddo
              if(dabs(yk1(ik1)-0.d0).lt.tol) ngeff=1
            else
              if(dabs(xk1(ik1)-yk1(ik1)).lt.tol) then
                ngeff=1
              else
                do ig=2,ng-1
                  iskip(ig)=1
                enddo
              endif
            endif
          elseif(idgroup.eq.'C2'.or.idgroup.eq.'Cs') then
            iskip(3) = 1
            iskip(4) = 1
            iskip(7) = 1
            if(iogroup.eq.'Csd1'.or.iogroup.eq.'Csd2') then
              iskip(5) = 1
              if(dabs(xk1(ik1)-yk1(ik1)).lt.tol) then
                ngeff=6
              elseif(dabs(yk1(ik1)-0.d0).lt.tol) then
                ngeff=2
              endif
            else
              iskip(6) = 1
              if(dabs(xk1(ik1)-yk1(ik1)).lt.tol) then
                ngeff=2
              elseif(dabs(yk1(ik1)-0.d0).lt.tol) then
                ngeff=5
              endif
            endif
          else
            if(dabs(xk1(ik1)-yk1(ik1)).lt.tol.or.
     &         dabs(yk1(ik1)-0.d0).lt.tol) ngeff=4
          endif
c
c       ===========
c        Hexagonal
c       ===========
c
        else
          r=dsqrt(xk1(ik1)**2+yk1(ik1)**2)
          phi=dasin(yk1(ik1)/r)
          if(idgroup.eq.'C6v') then
            ngeff=1
          elseif(idgroup.eq.'C6'.or.idgroup.eq.'C3vB') then
            if(dabs(phi-0.d0).lt.tol) then
              ngeff=1
            else
              ngeff=10
              do ig=2,ngeff-1
                iskip(ig)=1
              enddo
            endif
          elseif(idgroup.eq.'C3vA') then
            if(dabs(phi-pi/6.d0).lt.tol) then
              ngeff=1
            else
              ngeff=8
              do ig=2,ngeff-1
                iskip(ig)=1
              enddo
            endif
          elseif(idgroup.eq.'C3') then
            if(dabs(phi-pi/6.d0).lt.tol) then
              ngeff=2
            else
              do ig=3,ngeff-1
                if(ig.ne.8) iskip(ig)=1
              enddo
              if(dabs(phi-0.d0).lt.tol) iskip(8)=1
            endif
          elseif(idgroup.eq.'C2v') then
            if(iogroup.eq.'C2v2'.or.iogroup.eq.'C2v3') then
              ngeff=10
              do ig=2,7
                iskip(ig)=1
              enddo
              iskip(9)=1
              if(dabs(phi-pi/6.d0).lt.tol) iskip(8)=1
              if(dabs(phi-0.d0).lt.tol) ngeff=8
            else
              if(dabs(phi-pi/6.d0).lt.tol.or.dabs(phi-0.d0).lt.tol) then
                ngeff=2
              else
                ngeff=8
                do ig=3,ngeff-1
                  iskip(ig)=1
                enddo
              endif
            endif
          elseif(idgroup.eq.'Cs'.or.idgroup.eq.'C2') then
            if(iogroup.eq.'Csy'.or.iogroup.eq.'Csy2'.or.
     &      iogroup.eq.'Csy3') then
              do ig=4,7
                iskip(ig)=1
              enddo
              ngeff=10
              if(dabs(phi-pi/6.d0).lt.tol) then
                iskip(8)=1
                iskip(10)=1
              endif
              if(dabs(phi-0.d0).lt.tol) ngeff=3
            else
              do ig=3,ngeff-1
                if(ig.ne.4.and.ig.ne.7.and.ig.ne.8) iskip(ig)=1
              enddo
              if(dabs(phi-pi/6.d0).lt.tol) ngeff=4
              if(dabs(phi-0.d0).lt.tol) ngeff=7
            endif
          else
            if(dabs(phi-pi/6.d0).lt.tol.or.dabs(phi-0.d0).lt.tol) 
     &        ngeff=6
          endif
        endif
c       write(6,'('' ngeff='',i3)') ngeff 
c       call flush(6)
c
c
c-----------------------------------------------------------------
c Now, generate the IBZ applying the selected symmetry operations.
c-----------------------------------------------------------------
c
        do ig=1,ngeff
          if(iskip(ig).eq.0) then
            ik2=ik2+1
            xk2(ik2)=amat(1,1,ig)*xk1(ik1)+amat(1,2,ig)*yk1(ik1)
            yk2(ik2)=amat(2,1,ig)*xk1(ik1)+amat(2,2,ig)*yk1(ik1)
          endif
        enddo
      enddo
      nk2=ik2
c
c     write(6,'('' nk   ='',i3)') nk2
c     do ik=1,nk2
c       write(6,'(2f10.5)') xk2(ik),yk2(ik)
c     end do
c     call flush(6)
c
      return
      end
