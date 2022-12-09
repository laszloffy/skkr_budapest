c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine rgroup(i2dlat,idgroup,ng,a,melem)
c =====================
c
c Obtain real matrices of point group symmetry operations.
c Includes the 13 symmorphic 2D groups.
c
c======================
c
      implicit real*8 (a-h,o-z)
c
      logical group
      character*4 idgroup
      real*8 a(2,2,melem),a1(2,2),a2(2,2)
c
      pi=dacos(-1.d0)
      call rzero(a,4*melem) 
c
c First check compatibility of translational symmetry and point group
c
c=============================================================
c                          Oblique
c=============================================================
c
      if(i2dlat.eq.1) then
        group=idgroup.eq.'C1'.or.idgroup.eq.'C2'
        if(.not.group) then
          write(6,*) 'Revise input idgroup:',idgroup
          write(6,*) 'For oblique lattice it should be C1 or C2!'
          stop
        end if
c
c=============================================================
c                        Rectangular
c=============================================================
c
      else if(i2dlat.eq.2.or.i2dlat.eq.3) then
        group=idgroup.eq.'C1'.or.idgroup.eq.'Cs'.or.idgroup.eq.'C2'
     >        .or.idgroup.eq.'C2v'
        if(.not.group) then
          write(6,*) 'Revise input idgroup:',idgroup
          write(6,*) 'For rectangular lattice it should be C1, Cs, 
     > C2 or C2v!'
          stop
        end if
c
c=============================================================
c                         Square
c=============================================================
c
      else if(i2dlat.eq.4) then
        group=idgroup.eq.'C1'.or.idgroup.eq.'Cs'.or.idgroup.eq.'C2'
     >  .or.idgroup.eq.'C2v'.or.idgroup.eq.'C4'.or.idgroup.eq.'C4v'
        if(.not.group) then
          write(6,*) 'Revise input idgroup:',idgroup
          write(6,*) 'For square lattice it should be C1, Cs, C2, 
     > C2v, C4, or C4v!'
          stop
        end if
c
c=============================================================
c                       Hexagonal
c=============================================================
c
      else if(i2dlat.eq.5) then
        group=idgroup.eq.'C1'.or.idgroup.eq.'Cs'
     >   .or.idgroup.eq.'C2'.or.idgroup.eq.'C2v'
     >   .or.idgroup.eq.'C3'.or.idgroup.eq.'C3vA'.or.idgroup.eq.'C3vB'
     >   .or.idgroup.eq.'C6'.or.idgroup.eq.'C6v'
        if(.not.group) then
          write(6,*) 'Revise input idgroup:',idgroup
          write(6,*) 'For hexagonal lattice it should be C1, Cs, C2, 
     > C2v, C3, C3vA, C3vB, C6, or C6v!'
          stop
        end if
c
      else
c
        write(6,'('' RGROUP: unknown lattice type'',i3,
     >  '' recognized.'')') i2dlat
c
      end if
c
c
c Set up real-space rotation matrices for each of the point group
c elements
c
c     write(6,*) 'idgroup=',idgroup
c
c     ------------------------
      if(idgroup.eq.'C1') then
c     ------------------------
c
        ng=1
        if(melem.lt.ng) stop 
     &  'Increase parameter melem to 1'
c E
        a(1,1,1) = 1.d0
        a(2,2,1) = 1.d0
c
c     ----------------------------
      elseif(idgroup.eq.'Cs') then
c     ----------------------------
c
        ng=2
        if(melem.lt.ng) stop 
     &  'Increase parameter melem to 2'
c E
        a(1,1,1) = 1.d0
        a(2,2,1) = 1.d0
c Sx
        a(1,1,2) = 1.d0
        a(2,2,2) =-1.d0
c
c     ----------------------------
      elseif(idgroup.eq.'C2') then
c     ----------------------------
c
        ng=2
        if(melem.lt.ng) stop 
     &  'Increase parameter melem to 2'
c E
        a(1,1,1) = 1.d0
        a(2,2,1) = 1.d0
c C2
        a(1,1,2) = dcos(pi)
        a(1,2,2) =-dsin(pi)
        a(2,1,2) = dsin(pi)
        a(2,2,2) = dcos(pi)
c
c     ----------------------------
      elseif(idgroup.eq.'C3') then
c     ----------------------------
        ng=3
        if(melem.lt.ng) stop 
     &  'Increase parameter melem to 3'
c E
        a(1,1,1) = 1.d0
        a(2,2,1) = 1.d0
c C3+
        a(1,1,2) = dcos(2.d0*pi/3.d0)
        a(1,2,2) =-dsin(2.d0*pi/3.d0)
        a(2,1,2) = dsin(2.d0*pi/3.d0)
        a(2,2,2) = dcos(2.d0*pi/3.d0)
c C3-
        a(1,1,3) = dcos(-2.d0*pi/3.d0)
        a(1,2,3) =-dsin(-2.d0*pi/3.d0)
        a(2,1,3) = dsin(-2.d0*pi/3.d0)
        a(2,2,3) = dcos(-2.d0*pi/3.d0)
c
c     ----------------------------
      elseif(idgroup.eq.'C4') then
c     ----------------------------
c
        ng=4
        if(melem.lt.ng) stop 
     &  'Increase parameter melem to 4'
c E
        a(1,1,1) = 1.d0
        a(2,2,1) = 1.d0
c C4+
        a(1,1,2) = dcos(pi/2.d0)
        a(1,2,2) =-dsin(pi/2.d0)
        a(2,1,2) = dsin(pi/2.d0)
        a(2,2,2) = dcos(pi/2.d0)
c C4-
        a(1,1,3) = dcos(-pi/2.d0)
        a(1,2,3) =-dsin(-pi/2.d0)
        a(2,1,3) = dsin(-pi/2.d0)
        a(2,2,3) = dcos(-pi/2.d0)
c C2
        a(1,1,4) = dcos(pi)
        a(1,2,4) =-dsin(pi)
        a(2,1,4) = dsin(pi)
        a(2,2,4) = dcos(pi)
c
c     ----------------------------
      elseif(idgroup.eq.'C6') then
c     ----------------------------
        ng=6
        if(melem.lt.ng) stop 
     &  'Increase parameter melem to 6'
c E
        a(1,1,1) = 1.d0
        a(2,2,1) = 1.d0
c C6+
        a(1,1,2) = dcos(2.d0*pi/6.d0)
        a(1,2,2) =-dsin(2.d0*pi/6.d0)
        a(2,1,2) = dsin(2.d0*pi/6.d0)
        a(2,2,2) = dcos(2.d0*pi/6.d0)
c C6-
        a(1,1,3) = dcos(-2.d0*pi/6.d0)
        a(1,2,3) =-dsin(-2.d0*pi/6.d0)
        a(2,1,3) = dsin(-2.d0*pi/6.d0)
        a(2,2,3) = dcos(-2.d0*pi/6.d0)
c C3+
        a(1,1,4) = dcos(2.d0*pi/3.d0)
        a(1,2,4) =-dsin(2.d0*pi/3.d0)
        a(2,1,4) = dsin(2.d0*pi/3.d0)
        a(2,2,4) = dcos(2.d0*pi/3.d0)
c C3-
        a(1,1,5) = dcos(-2.d0*pi/3.d0)
        a(1,2,5) =-dsin(-2.d0*pi/3.d0)
        a(2,1,5) = dsin(-2.d0*pi/3.d0)
        a(2,2,5) = dcos(-2.d0*pi/3.d0)
c C2
        a(1,1,6) = dcos(pi)
        a(1,2,6) =-dsin(pi)
        a(2,1,6) = dsin(pi)
        a(2,2,6) = dcos(pi)
c
c     -----------------------------
      elseif(idgroup.eq.'C2v') then
c     -----------------------------
c
        ng=4
        if(melem.lt.ng) stop 
     &  'Increase parameter melem to 4'
c E
        a(1,1,1) = 1.d0
        a(2,2,1) = 1.d0
c C2
        a(1,1,2) = dcos(pi)
        a(1,2,2) =-dsin(pi)
        a(2,1,2) = dsin(pi)
        a(2,2,2) = dcos(pi)
c Sx
        a(1,1,3) =-1.d0
        a(2,2,3) = 1.d0
c Sy
        a(1,1,4) = 1.d0
        a(2,2,4) =-1.d0
c
c     ------------------------------
      elseif(idgroup.eq.'C3vA') then
c     ------------------------------
        ng=6
        if(melem.lt.ng) stop 
     &  'Increase parameter melem to 6'
c E
        a(1,1,1) = 1.d0
        a(2,2,1) = 1.d0
c C3+
        a(1,1,2) = dcos(2.d0*pi/3.d0)
        a(1,2,2) =-dsin(2.d0*pi/3.d0)
        a(2,1,2) = dsin(2.d0*pi/3.d0)
        a(2,2,2) = dcos(2.d0*pi/3.d0)
c C3-
        a(1,1,3) = dcos(-2.d0*pi/3.d0)
        a(1,2,3) =-dsin(-2.d0*pi/3.d0)
        a(2,1,3) = dsin(-2.d0*pi/3.d0)
        a(2,2,3) = dcos(-2.d0*pi/3.d0)
c Sv1
        a(1,1,4) = 1.d0
        a(2,2,4) =-1.d0
c Sv2
        a1(1,1) = dcos(-pi/3.d0)
        a1(1,2) =-dsin(-pi/3.d0)
        a1(2,1) = dsin(-pi/3.d0)
        a1(2,2) = dcos(-pi/3.d0)
        a(1,1,5)= 1.d0
        a(2,2,5)=-1.d0
        a2(1,1) = dcos(pi/3.d0)
        a2(1,2) =-dsin(pi/3.d0)
        a2(2,1) = dsin(pi/3.d0)
        a2(2,2) = dcos(pi/3.d0)
        call tripmtr(a1,a(1,1,5),a2,2,2,2)
c Sv3
        a1(1,1) = dcos(pi/3.d0)
        a1(1,2) =-dsin(pi/3.d0)
        a1(2,1) = dsin(pi/3.d0)
        a1(2,2) = dcos(pi/3.d0)
        a(1,1,6)= 1.d0
        a(2,2,6)=-1.d0
        a2(1,1) = dcos(-pi/3.d0)
        a2(1,2) =-dsin(-pi/3.d0)
        a2(2,1) = dsin(-pi/3.d0)
        a2(2,2) = dcos(-pi/3.d0)
        call tripmtr(a1,a(1,1,6),a2,2,2,2)
c
c     ------------------------------
      elseif(idgroup.eq.'C3vB') then
c     ------------------------------
        ng=6
        if(melem.lt.ng) stop 
     &  'Increase parameter melem to 6'
c E
        a(1,1,1) = 1.d0
        a(2,2,1) = 1.d0
c C3+
        a(1,1,2) = dcos(2.d0*pi/3.d0)
        a(1,2,2) =-dsin(2.d0*pi/3.d0)
        a(2,1,2) = dsin(2.d0*pi/3.d0)
        a(2,2,2) = dcos(2.d0*pi/3.d0)
c C3-
        a(1,1,3) = dcos(-2.d0*pi/3.d0)
        a(1,2,3) =-dsin(-2.d0*pi/3.d0)
        a(2,1,3) = dsin(-2.d0*pi/3.d0)
        a(2,2,3) = dcos(-2.d0*pi/3.d0)
c Sv1
        a(1,1,4) =-1.d0
        a(2,2,4) = 1.d0
c Sv2
        a1(1,1) = dcos(pi/6.d0)
        a1(1,2) =-dsin(pi/6.d0)
        a1(2,1) = dsin(pi/6.d0)
        a1(2,2) = dcos(pi/6.d0)
        a(1,1,5)= 1.d0
        a(2,2,5)=-1.d0
        a2(1,1) = dcos(-pi/6.d0)
        a2(1,2) =-dsin(-pi/6.d0)
        a2(2,1) = dsin(-pi/6.d0)
        a2(2,2) = dcos(-pi/6.d0)
        call tripmtr(a1,a(1,1,5),a2,2,2,2)
c Sv3
        a1(1,1) = dcos(-pi/6.d0)
        a1(1,2) =-dsin(-pi/6.d0)
        a1(2,1) = dsin(-pi/6.d0)
        a1(2,2) = dcos(-pi/6.d0)
        a(1,1,6)= 1.d0
        a(2,2,6)=-1.d0
        a2(1,1) = dcos(pi/6.d0)
        a2(1,2) =-dsin(pi/6.d0)
        a2(2,1) = dsin(pi/6.d0)
        a2(2,2) = dcos(pi/6.d0)
        call tripmtr(a1,a(1,1,6),a2,2,2,2)
c
c     -----------------------------
      elseif(idgroup.eq.'C4v') then
c     -----------------------------
        ng=8
        if(melem.lt.ng) stop 
     &  'Increase parameter melem to 8'
c E
        a(1,1,1) = 1.d0
        a(2,2,1) = 1.d0
c C4+
        a(1,1,2) = dcos(pi/2.d0)
        a(1,2,2) =-dsin(pi/2.d0)
        a(2,1,2) = dsin(pi/2.d0)
        a(2,2,2) = dcos(pi/2.d0)
c C4-
        a(1,1,3) = dcos(-pi/2.d0)
        a(1,2,3) =-dsin(-pi/2.d0)
        a(2,1,3) = dsin(-pi/2.d0)
        a(2,2,3) = dcos(-pi/2.d0)
c C2
        a(1,1,4) = dcos(pi)
        a(1,2,4) =-dsin(pi)
        a(2,1,4) = dsin(pi)
        a(2,2,4) = dcos(pi)
c Sv1(=Sy)
        a(1,1,5) =-1.d0
        a(2,2,5) = 1.d0
c Sv2(=Sx)
        a(1,1,6) = 1.d0
        a(2,2,6) =-1.d0
c Sd1
        a(1,2,7) =-1.d0
        a(2,1,7) =-1.d0
c Sd2
        a(1,2,8) = 1.d0
        a(2,1,8) = 1.d0
c
c     -----------------------------
      elseif(idgroup.eq.'C6v') then
c     -----------------------------
        ng=12
        if(melem.lt.ng) stop 
     &  'Increase parameter melem to 12'
c E
        a(1,1,1) = 1.d0
        a(2,2,1) = 1.d0
c C6+
        a(1,1,2) = dcos(2.d0*pi/6.d0)
        a(1,2,2) =-dsin(2.d0*pi/6.d0)
        a(2,1,2) = dsin(2.d0*pi/6.d0)
        a(2,2,2) = dcos(2.d0*pi/6.d0)
c C6-
        a(1,1,3) = dcos(-2.d0*pi/6.d0)
        a(1,2,3) =-dsin(-2.d0*pi/6.d0)
        a(2,1,3) = dsin(-2.d0*pi/6.d0)
        a(2,2,3) = dcos(-2.d0*pi/6.d0)
c C3+
        a(1,1,4) = dcos(2.d0*pi/3.d0)
        a(1,2,4) =-dsin(2.d0*pi/3.d0)
        a(2,1,4) = dsin(2.d0*pi/3.d0)
        a(2,2,4) = dcos(2.d0*pi/3.d0)
c C3-
        a(1,1,5) = dcos(-2.d0*pi/3.d0)
        a(1,2,5) =-dsin(-2.d0*pi/3.d0)
        a(2,1,5) = dsin(-2.d0*pi/3.d0)
        a(2,2,5) = dcos(-2.d0*pi/3.d0)
c C2
        a(1,1,6) = dcos(pi)
        a(1,2,6) =-dsin(pi)
        a(2,1,6) = dsin(pi)
        a(2,2,6) = dcos(pi)
c Sd1
        a(1,1,7) =-1.d0
        a(2,2,7) = 1.d0
c Sd2
        a1(1,1) = dcos(pi/6.d0)
        a1(1,2) =-dsin(pi/6.d0)
        a1(2,1) = dsin(pi/6.d0)
        a1(2,2) = dcos(pi/6.d0)
        a(1,1,8)= 1.d0
        a(2,2,8)=-1.d0
        a2(1,1) = dcos(-pi/6.d0)
        a2(1,2) =-dsin(-pi/6.d0)
        a2(2,1) = dsin(-pi/6.d0)
        a2(2,2) = dcos(-pi/6.d0)
        call tripmtr(a1,a(1,1,8),a2,2,2,2)
c Sd3
        a1(1,1) = dcos(-pi/6.d0)
        a1(1,2) =-dsin(-pi/6.d0)
        a1(2,1) = dsin(-pi/6.d0)
        a1(2,2) = dcos(-pi/6.d0)
        a(1,1,9)= 1.d0
        a(2,2,9)=-1.d0
        a2(1,1) = dcos(pi/6.d0)
        a2(1,2) =-dsin(pi/6.d0)
        a2(2,1) = dsin(pi/6.d0)
        a2(2,2) = dcos(pi/6.d0)
        call tripmtr(a1,a(1,1,9),a2,2,2,2)
c Sv1
        a(1,1,10) = 1.d0
        a(2,2,10) =-1.d0
c Sv2
        a1(1,1) = dcos(-pi/3.d0)
        a1(1,2) =-dsin(-pi/3.d0)
        a1(2,1) = dsin(-pi/3.d0)
        a1(2,2) = dcos(-pi/3.d0)
        a(1,1,11)= 1.d0
        a(2,2,11)=-1.d0
        a2(1,1) = dcos(pi/3.d0)
        a2(1,2) =-dsin(pi/3.d0)
        a2(2,1) = dsin(pi/3.d0)
        a2(2,2) = dcos(pi/3.d0)
        call tripmtr(a1,a(1,1,11),a2,2,2,2)
c Sv3
        a1(1,1) = dcos(pi/3.d0)
        a1(1,2) =-dsin(pi/3.d0)
        a1(2,1) = dsin(pi/3.d0)
        a1(2,2) = dcos(pi/3.d0)
        a(1,1,12)= 1.d0
        a(2,2,12)=-1.d0
        a2(1,1) = dcos(-pi/3.d0)
        a2(1,2) =-dsin(-pi/3.d0)
        a2(2,1) = dsin(-pi/3.d0)
        a2(2,2) = dcos(-pi/3.d0)
        call tripmtr(a1,a(1,1,12),a2,2,2,2)
c
      endif
c
      return
      end

