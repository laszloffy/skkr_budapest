c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
c ==========================================================
      subroutine primkr2d(ilat,idgroup,a2d,b2d,d2d,arot)
c ==========================================================
c
c First, generate/read the 2-D basis vectors (according to Cunningham).
c
c Then, rotate them so that idgroup matches the point-group
c operations implemented in the code (rgroup and sphbas).
c
c input:  ilat - type of 2D-lattice
c         a2d  - 2D-lattice parameter
c         b2d  - asymmetry ratio (for oblique and rectangular)
c         idgroup (changed on output) - point group symmetry
c
c output: ax,ay  - Real space lattice vectors in common brav2d
c         bx,by  - Reciprocal space lattice vectors in common brill2d
c         arot   - Rotation of input vectors to get implemented symmetry
c         bxc,byc- Reciprocal basis oriented according to Cunningham
c                  (used in kmesh2d)
c         phik   - Rotation of calculated IBZ to get implemented 
c                  symmetry (used in kmesh2d - ibz)
c
      implicit real*8 (a-h,o-z)
      logical group
      character*4 idgroup,iogroup
      character*32 brav(5)
      dimension v1(2),v2(2),axc(2),ayc(2)
      common/brav2d/ ax(2),ay(2)
      common/brill2d/ bx(2),by(2)
      common/cunnbz/ bxc(2),byc(2),phik,iogroup
      data brav /'oblique','centered rectangular',
     1           'primitive rectangular','square','hexagonal'/
      data tol/1.0d-10/
c
      pi=4.d0*datan(1.d0)
      arot=0.d0
c     if(idgroup.eq.'C2vd') stop 'Bug!!! (Symmetry group C2vd)' 
      iogroup=idgroup
c
c                     ********************
c                        Direct lattice
c                     ********************
c
      ax(1) = a2d
      ax(2) = 0.d0
      axc(1)=ax(1)
      axc(2)=ax(2)
c
c=============================================================
c                          Oblique
c=============================================================    
      if(ilat.eq.1) then          
c
        if(d2d.gt.0.5d0.or.b2d.ge.1) stop
     &  'Error in input values for beta,delta - oblique lattice'
        ay(1) = a2d*d2d
        ay(2) = a2d*b2d
        ayc(1)=ay(1)
        ayc(2)=ay(2)
c
        group=idgroup.eq.'C1'.or.idgroup.eq.'C2'
        if(.not.group) then
          write(6,*) 'Revise input idgroup:',idgroup
          write(6,*) 'For oblique lattice it should be C1 or C2!'
          stop
        end if 
c
c=============================================================
c                    Centered rectangular
c=============================================================    
      elseif(ilat.eq.2) then
c
        if(b2d.ge.1) stop
     &  'Error in input value for beta - c.rectangular lattice'
        d2d = 0.5d0
        ay(1) = a2d*d2d
        ay(2) = a2d*d2d*b2d 
        ayc(1)=ay(1)
        ayc(2)=ay(2)
c
        group=idgroup.eq.'C1'.or.idgroup.eq.'C2'.or.
     >        idgroup.eq.'Csx'.or.idgroup.eq.'Csy'.or.
     >        idgroup.eq.'C2v'.or.idgroup.eq.'C2vx'
        if(.not.group) then
          write(6,*) 'Revise input idgroup:',idgroup
          write(6,*) 'For rectangular lattice it should be:'
          write(6,*) ' C1, C2, Csx, Csy or C2v (=C2vx)!'
          stop
        end if                    
        if(idgroup.eq.'Csy') then
          arot=pi/2.d0
          call vecrot(ax,arot)
          call vecrot(ay,arot)
        endif
        if(idgroup.eq.'Csx'.or.idgroup.eq.'Csy') idgroup='Cs'
        if(idgroup.eq.'C2vx') idgroup='C2v'
c
c=============================================================
c                    Primitive rectangular
c=============================================================    
      elseif(ilat.eq.3) then
c
        if(b2d.ge.1) stop
     &  'Error in input value for beta - p.rectangular lattice'
        d2d = 0.d0
        ay(1) = a2d*d2d
        ay(2) = a2d*b2d 
        ayc(1)=ay(1)
        ayc(2)=ay(2)
c
        group=idgroup.eq.'C1'.or.idgroup.eq.'C2'.or.
     >        idgroup.eq.'Csx'.or.idgroup.eq.'Csy'.or.
     >        idgroup.eq.'C2v'.or.idgroup.eq.'C2vx'
        if(.not.group) then
          write(6,*) 'Revise input idgroup:',idgroup
          write(6,*) 'For rectangular lattice it should be:'
          write(6,*) ' C1, C2, Csx, Csy or C2v (=C2vx)!'
          stop
        end if                    
        if(idgroup.eq.'Csy') then
          arot=pi/2.d0
          call vecrot(ax,arot)
          call vecrot(ay,arot)
        endif
        if(idgroup.eq.'Csx'.or.idgroup.eq.'Csy') idgroup='Cs'
        if(idgroup.eq.'C2vx') idgroup='C2v'
c
c=============================================================
c                          Square
c=============================================================    
      elseif(ilat.eq.4) then
c
        d2d = 0.d0
        b2d = 1.d0
        ay(1) = a2d*d2d
        ay(2) = a2d*b2d 
        ayc(1)=ay(1)
        ayc(2)=ay(2)
c
        group=idgroup.eq.'C1'.or.idgroup.eq.'C2'.or.
     >        idgroup.eq.'C4'.or.idgroup.eq.'C4v'.or.
     >        idgroup.eq.'C2vx'.or.idgroup.eq.'C2vd'.or.
     >        idgroup.eq.'Csx'.or.idgroup.eq.'Csy'.or.
     >        idgroup.eq.'Csd1'.or.idgroup.eq.'Csd2'
        if(.not.group) then
          write(6,*) 'Revise input idgroup:',idgroup
          write(6,*) 'For square lattice it should be:'
       write(6,*) 'C1, C2, C4, C4v, C2vx, C2vd, Csx, Csy, Csd1 or Csd2!'
          stop
        endif
        if(idgroup.eq.'Csy') then
          arot=pi/2.d0
          call vecrot(ax,arot)
          call vecrot(ay,arot)
        elseif(idgroup.eq.'Csd1'.or.idgroup.eq.'C2vd') then
          arot=pi/4.d0
          call vecrot(ax,arot)
          call vecrot(ay,arot)
        elseif(idgroup.eq.'Csd2') then
          arot=-pi/4.d0
          call vecrot(ax,arot)
          call vecrot(ay,arot)
        endif
        if(idgroup.eq.'Csx'.or.idgroup.eq.'Csy'.or.
     >     idgroup.eq.'Csd1'.or.idgroup.eq.'Csd2') idgroup='Cs'
        if(idgroup.eq.'C2vx'.or.idgroup.eq.'C2vd') idgroup='C2v'
c
c=============================================================
c                          Hexagonal
c=============================================================    
      elseif(ilat.eq.5) then
c
        b2d = 1.d0
        ay(1) = a2d/2.d0
        ay(2) = a2d*dsqrt(3.d0)/2.d0
        ayc(1)=ay(1)
        ayc(2)=ay(2)
c
        group=idgroup.eq.'C1'.or.idgroup.eq.'C2'.or.
     >        idgroup.eq.'C3'.or.idgroup.eq.'C6'.or.
     >        idgroup.eq.'C3vA'.or.idgroup.eq.'C3vB'.or.
     >        idgroup.eq.'C6v'.or.idgroup.eq.'C2vx'.or.
     >        idgroup.eq.'C2v2'.or.idgroup.eq.'C2v3'.or.
     >        idgroup.eq.'Csx'.or.idgroup.eq.'Csx2'.or.
     >        idgroup.eq.'Csx3'.or.idgroup.eq.'Csy'.or.
     >        idgroup.eq.'Csy2'.or.idgroup.eq.'Csy3'
        if(.not.group) then
          write(6,*) 'Revise input idgroup:',idgroup
          write(6,*) 'For hexagonal lattice it should be:'
          write(6,*) 'C1, C2, C3, C6, C3vA, C3vB, C6v, '
          write(6,*)
     >    'C2vx, C2v2, C2v3, Csx, Csx2, Csx3, Csy, Csy2 or Csy3!'
          stop
        end if
        if(idgroup.eq.'Csy') then
          arot=pi/2.d0
          call vecrot(ax,arot)
          call vecrot(ay,arot)
        elseif(idgroup.eq.'Csx2') then
          arot=pi/3.d0
          call vecrot(ax,arot)
          call vecrot(ay,arot)
        elseif(idgroup.eq.'Csx3') then
          arot=-pi/3.d0
          call vecrot(ax,arot)
          call vecrot(ay,arot)
        elseif(idgroup.eq.'Csy2'.or.idgroup.eq.'C2v2') then
          arot=-pi/6.d0
          call vecrot(ax,arot)
          call vecrot(ay,arot)
        elseif(idgroup.eq.'Csy3'.or.idgroup.eq.'C2v3') then
          arot=pi/6.d0
          call vecrot(ax,arot)
          call vecrot(ay,arot)
        endif
        if(idgroup.eq.'Csx'.or.idgroup.eq.'Csy'.or.
     >     idgroup.eq.'Csx2'.or.idgroup.eq.'Csy2'.or.
     >     idgroup.eq.'Csx3'.or.idgroup.eq.'Csy3') idgroup='Cs'
        if(idgroup.eq.'C2vx'.or.idgroup.eq.'C2v2'.or.
     >     idgroup.eq.'C2v3') idgroup='C2v'
c
      else
        write(6,*) 'Error!! i2dlat should be .ge.1 and .le.5',ilat
        stop
      endif
c
      write (6,20) brav(ilat),idgroup
c
c                  **************************
c                      Reciprocal lattice
c                  **************************
c--------
c     Actual lattice
c
      v1(1) =  ay(2)
      v1(2) = -ay(1)
      v2(1) = -ax(2)
      v2(2) =  ax(1)
c
      dnt1 = ax(1)*v1(1) + ax(2)*v1(2)
      dnt2 = ay(1)*v2(1) + ay(2)*v2(2)
c
      do i=1,2
        bx(i) = 2.d0*pi*v1(i)/dnt1
        by(i) = 2.d0*pi*v2(i)/dnt2
      enddo
c
      phik=arot
c
c--------
c     Cunningham orientation to generate k-mesh
c
      v1(1) =  ayc(2)
      v1(2) = -ayc(1)
      v2(1) = -axc(2)
      v2(2) =  axc(1)
c
      dnt1 = axc(1)*v1(1) + axc(2)*v1(2)
      dnt2 = ayc(1)*v2(1) + ayc(2)*v2(2)
c
      do i=1,2
        bxc(i) = 2.d0*pi*v1(i)/dnt1
        byc(i) = 2.d0*pi*v2(i)/dnt2
      enddo
c
c     Check: brill2d should come from a rotation of arot of bxc,byc 
      do i=1,2
        v1(i)=bxc(i)
        v2(i)=byc(i)
      enddo
      call vecrot(v1,arot)
      call vecrot(v2,arot)
      do i=1,2
        dif1=dabs(bx(i)-v1(i))
        dif2=dabs(by(i)-v2(i))
        if(dif1.gt.tol.or.dif2.gt.tol) stop
     &    'PRIMKR2D: Revise rotation of reciprocal lattice'
      enddo
c
      return
   20 format(' lattice type is: ',a32,'(Point group=',a4,')')
      end
