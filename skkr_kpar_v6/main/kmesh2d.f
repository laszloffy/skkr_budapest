c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
c*********************************************************
c  REVISE TREATMENT OF OBLIQUE LATTICES IN speck 
c*********************************************************
      subroutine kmesh2d(intbz,kunning,iset,park,xk,yk,wk,nk,mkpar)
c=======================
c
c   This routine initializes k-points and weights according 
c   to i2dlat and idgroup.
c   They are generated using Cunnigham coordinate system, and then
c   rotated into the orientation defined in primkr2d.
c
c   If kunning=0 a simple equidistant point set is used.
c   If kunning=1 the Cunningham k-points are used. 
c                    [Cunningham, PRB Vol. 10, 4988 (1974)].
c
      implicit real*8 (a-h,o-z)
      character*4 idgroup,iogroup
      dimension xk0(mkpar),yk0(mkpar),wk0(mkpar),veck(2)
      dimension xk(mkpar),yk(mkpar),wk(mkpar),park(2)
      common/latt2d/idgroup,i2dlat,a2d,beta,delta
      common/brill2d/ bx(2),by(2)
      common/cunnbz/ bxc(2),byc(2),phik,iogroup
      data anul,tol/0.0d0,1.0d-10/
c 
      pi=dacos(-1.0d0)
c 
c--------------------------------------------------------------
c   First, generate minimum irreducible wedge according to iset
c   iset=0:    1 k-point
c   iset.gt.0: oblique     (i2dlat.eq.1)   => 'C2'  => 1/2  BZ
c              rectangular (i2dlat.eq.2,3) => 'C2v' => 1/4  BZ
c              square      (i2dlat.eq.4)   => 'C4v' => 1/8  BZ
c              hexagonal   (i2dlat.eq.5)   => 'C6v' => 1/12 BZ
c--------------------------------------------------------------
      if(iset.eq.0) then
        nk0=1
        xk0(1)=park(1)
        yk0(1)=park(2)
c       xk0(1)=park(1)*pi/a2d
c       yk0(1)=park(2)*pi/a2d
        wk0(1)=1.d0
      else
c
        IF(i2dlat.EQ.1.and.kunning.eq.1) THEN
          WRITE(6,*) 'WARNING!! '
          WRITE(6,*) 
     &    'Cunningham points not used for oblique lattice !!'
          WRITE(6,*) 
     &    '==> Using equidistant point-sampling'
          kunning=0
        ENDIF
c
        if(kunning.eq.0) then
           call speck2(iset,xk0,yk0,wk0,nk0,i2dlat,a2d,beta,delta)
        else 
           call setspk(i2dlat,iset,nk0)
           call speck(xk0,yk0,wk0,nk0,i2dlat)
        end if
      end if
c
      if(nk0.gt.mkpar) then
        write(6,*) 'Increase mkpar to',nk0
        stop
      endif
      if(nk0.gt.mkpar) then
        write(6,*) 'KMESH2D: Increase parameter mkpar to',nk0
        stop
      endif
c
c     write(6,'('' intbz='',i3)') intbz
c     write(6,'('' iset ='',i3)') iset
c     write(6,'('' nk0  ='',i3)') nk0
c     do ik=1,nk0
c       write(6,'(2f10.5)') xk0(ik),yk0(ik)
c     end do
c     call flush(6)
c
c--------------------------------------------------------------
c   Now, either the actual IBZ or the full BZ is obtained.
c--------------------------------------------------------------
      if(intbz.le.1) then
        if(iset.eq.0) then
          nk=nk0
          xk(1)=xk0(1)
          yk(1)=yk0(1)
          wk(1)=wk0(1)
        else
          call ibz(i2dlat,idgroup,iogroup,kunning,xk0,yk0,nk0,
     &             xk,yk,nk)
        endif
      elseif(intbz.gt.1) then
        call fullbz(i2dlat,iset,kunning,phik,xk0,yk0,nk0,
     &              xk,yk,nk)
      end if
c
      if(nk.gt.mkpar) then
        write(6,*) 'Increase mkpar to',nk
        stop
      endif
c
c     write(6,'('' nk   ='',i3)') nk 
c     do ik=1,nk
c       write(6,'(2f10.5)') xk(ik),yk(ik)
c     end do
c     call flush(6)
c
c--------------------------------------------------------------
c   Rotate the final k-mesh to the actual lattice orientation.
c   Give new weighting factors. 
c--------------------------------------------------------------
c
      phi = phik
c
c     For some special cases skip rotation of IBZ
      if(intbz.le.1.and.iset.ne.0) then
        if(i2dlat.eq.4) then
          if(iogroup.eq.'Csy') phi=0.d0
          if(iogroup.eq.'Csd2') phi=-phik
        elseif(i2dlat.eq.5) then
          if(iogroup.eq.'Csx2'.or.iogroup.eq.'Csx3') phi=0.d0
          if(iogroup.eq.'Csy2'.or.iogroup.eq.'Csy3') phi=pi/2.d0
          if(iogroup.eq.'C2v2') phi=-phik
        endif
      endif
c
      if(intbz.le.1.or.iset.ne.0) then
        if(dabs(phi-0.d0).gt.tol) then
          do i=1,nk
            veck(1)=xk(i)
            veck(2)=yk(i)
            call vecrot(veck,phi)
            xk(i)=veck(1)
            yk(i)=veck(2)
          enddo
        endif
      endif
c
c Determine relative weights for the case of C1 (fullbz) by taking care
c of the k-point falling on the boundary of the Brillouin zone!
c
      do ik=1,nk
        ibound=0
        do 50 i=-2,2
        do 50 j=-2,2
          if(i.eq.0.and.j.eq.0) goto 50
          gx=dfloat(i)*bx(1)+dfloat(j)*by(1)
          gy=dfloat(i)*bx(2)+dfloat(j)*by(2)
          g2o2=0.5d0*(gx*gx+gy*gy)
          gk=gx*xk(ik)+gy*yk(ik)
          if(dabs(gk-g2o2).lt.tol) ibound=ibound+1
  50    continue
c inside the BZ
        if(ibound.eq.0) wk(ik)=1.0d0
c on the edge of the BZ
        if(ibound.eq.1) wk(ik)=0.5d0
c at the corner of the BZ
        if(ibound.eq.2) then
          if(i2dlat.eq.3.or.i2dlat.eq.4) then
            wk(ik)=0.25d0
          else
            wk(ik)=1.0d0/3.0d0
          end if
        end if
      end do
c
c Consider optionally symmetry groups taking care of k-points 
c at the boundaries of the actual IBZ
c
      if(intbz.le.1.and.iset.ne.0) then
        do ik=1,nk
          r=dsqrt(xk(ik)**2+yk(ik)**2)
c         phi=dasin(yk(ik)/r)          
          call findangles(xk(ik),yk(ik),anul,theta,phi)
          if(idgroup.eq.'Cs'.or.idgroup.eq.'C2') then
            if(i2dlat.eq.1.and.kunning.eq.0) then
              if(dabs(xk(ik)-0.d0).lt.tol) wk(ik)=0.5d0*wk(ik)
            else
              if(dabs(yk(ik)-0.d0).lt.tol) wk(ik)=0.5d0*wk(ik)
            end if
          elseif(idgroup.eq.'C4'.or.idgroup.eq.'C2v') then
            if(dabs(xk(ik)-0.d0).lt.tol) then
              if(dabs(yk(ik)-0.d0).lt.tol) then
                wk(ik)=0.25d0*wk(ik)
              else
                wk(ik)=0.5d0*wk(ik)
              end if
            elseif(dabs(yk(ik)-0.d0).lt.tol) then
              wk(ik)=0.5d0*wk(ik)
            end if
          elseif(idgroup.eq.'C4v') then
            if(dabs(xk(ik)-0.d0).lt.tol) then
              if(dabs(yk(ik)-0.0d0).lt.tol) then
                wk(ik)=0.125d0*wk(ik)
              else
                wk(ik)=0.5d0*wk(ik)
              end if
            elseif(dabs(xk(ik)-yk(ik)).lt.tol) then
              wk(ik)=0.5d0*wk(ik)
            end if
          elseif(idgroup.eq.'C3') then
            if(dabs(yk(ik)-0.d0).lt.tol) then
              if(dabs(xk(ik)-0.d0).lt.tol) then
                wk(ik)=wk(ik)/3.0d0
              else
                wk(ik)=0.5d0*wk(ik)
              end if
            elseif(dabs(phi-2.0d0*pi/3.0d0).lt.tol) then
              wk(ik)=0.5d0*wk(ik)
            end if
          elseif(idgroup.eq.'C3vA') then
            if(dabs(yk(ik)-0.d0).lt.tol) then
              if(dabs(xk(ik)-0.d0).lt.tol) then
                wk(ik)=wk(ik)/6.0d0
              else
                wk(ik)=0.5d0*wk(ik)
              end if
            elseif(dabs(phi-pi/3.0d0).lt.tol) then
              wk(ik)=0.5d0*wk(ik)
            end if
          elseif(idgroup.eq.'C3vB'.or.idgroup.eq.'C6') then
            if(dabs(xk(ik)-0.d0).lt.tol.and.
     &         dabs(yk(ik)-0.d0).lt.tol) then
              wk(ik)=wk(ik)/6.0d0
            elseif(dabs(phi-pi/6.0d0).lt.tol.or.
     &             dabs(phi+pi/6.0d0).lt.tol) then
              wk(ik)=0.5d0*wk(ik)
            end if
          elseif(idgroup.eq.'C6v') then
            if(dabs(xk(ik)-0.d0).lt.tol) then
              if(dabs(yk(ik)-0.d0).lt.tol) then
                wk(ik)=wk(ik)/12.0d0
              else
                wk(ik)=0.5d0*wk(ik)
              end if
            elseif(dabs(phi-pi/6.0d0).lt.tol) then
              wk(ik)=0.5d0*wk(ik)
            end if
          end if
        end do
      end if
c
c Normalize weighting factors to 1
c
      wsum=0.d0
      do ik=1,nk
        wsum=wsum+wk(ik)
      end do
      do ik=1,nk
        wk(ik)=wk(ik)/wsum
      end do
c
      return 
      end
