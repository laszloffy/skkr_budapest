c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
c Read geometry- and layer-dependent input
c
      subroutine readgeom(nintfc,arot,sws)
!    &          bulk,bulkgeo,rightm,intbz,
!    >          arot,
!    &          swsl,concl,rbl,sxl,
!    &          sws,conc,rba,rbb,sxa,sxb,
!    &          swsr,concr,rbr,sxr,
!    &          lfix1,lfix2,linbw,igraph,iesublatt,lsublatt) 
c
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
      logical bulk,bulkgeo,linbw
c Laszlo 15/07/2003
      logical eqsubl
c
      character*1 rightm
      character*4 idgroup,imgroup,igroup(5)
      character*5 lat3d
      character*10 idpota(mintfc),idpotla(minprc),idpotra(minprc)
      character*10 idpotb(mintfc),idpotlb(minprc),idpotrb(minprc)
      character*40 tregion
      dimension rbl(3,minprc),rbr(3,minprc)
      dimension rba(3,mintfc),rbb(3,mintfc)
      dimension swsl(minprc),sws(mintfc),swsr(minprc)
      dimension sxl(minprc),sxr(minprc)
      dimension sxa(mintfc),sxb(mintfc)
      dimension conc(mintfc),concl(minprc),concr(minprc)
!     dimension isigba(mintfc),isigbb(mintfc)
c Laszlo 07/02/2005
      dimension iesubl(mtotal),iesublatt(mintfc),lsublatt(mintfc,melem)
c
      dimension cshift(3),phibw(mintfc)
      dimension sws0(minprc),smt0(minprc),sx0(0:1)
      dimension conc0(minprc),rb0(3,minprc)
      dimension rlay0(3),rd(2),rd0(2)
c
      common/brav2d/a1(2),a2(2)
      common/brav3dl/a1l(3),a2l(3),a3l(3)
      common/brav3dr/a1r(3),a2r(3),a3r(3)
c Laszlo 16/07/2003
      common/dirvec3dl/rvecl(3,mdir),volwsl,distrl,maxrl
c
      common/rrotmat/amat(2,2,melem),a0mat(2,2,melem),ng,ng0
      common/latt2d/idgroup,i2dlat,a2d,b2d,d2d
      common/lay2d/cvec(mtotal,3),nextra,nbulkl,nbulkr,
     &             nprc,ninprc(0:mprc+1)
      common/crystl2d/vol,volbz
      common/muftin/smtl(minprc),smt(mintfc),smtr(minprc) 
c     common/pot/potshift,b0,layshift1,layshift2,layb0,layb0f,
c    &           ib0,isigba,isigbb
      common /swscheck/ sws3d
c
      data tol/1.0d-10/
      data aunit/0.5291772d0/
      data pi/3.1415926535897932d0/
      data igroup/'C2','C2v','C2v','C4v','C6v'/
c
      fact=pi/180.d0
      bulk=.false.
      intbz=1
c
c
c***************************************************************
c
c Set general conditions of input file
c
c***************************************************************  
c
c     Decide if geometry will be read in (igraph.lt.0) 
c     or build up from a 3D-lattice (igraph.ge.0).
c
c     3 options: igraph.gt.0  => generate 3D structure and start SCF
c                igraph.eq.0  => generate 3D structure and stop
c                                (output file readgeom.out)
c                igraph.lt.0  => read layer by layer structure.
c
c
c     Indicate units of lenght: Angs. (iunit=0), a.u. (iunit=1)
c
c     ----------------------
      read(5,*) igraph,iunit
c     ----------------------
c
      if(iunit.ne.1.and.iunit.ne.0) stop
     &                    'Revise iunit on input file!! '
c
c
c
c***********************************
c***********************************
c Options 1-2: build-up structure
c***********************************
c***********************************
c
      if(igraph.ge.0) then   
c
        read(5,*) nprc,nextra
        read(5,'(a5)') lat3d
        read(5,*) lface,a3d
        read(5,*) btha,bpha,sx0(0)
        if(iunit.eq.0) a3d = a3d/aunit
c
        lfix1=0
        lfix2=0
c       Optional input for a linear Bloch Wall 
        linbw=.false.
c       read(5,*) linbw,phi1,phi2,lay1,lay2,shbw
c
c       Generate layer by layer structure and 2D lattice parameters.
        call struc3d(
     &          lat3d,lface,a3d,nprc,nextra,              !input
     &          mintfc,minprc,mprc,mtotal,
     &          i2dlat,idgroup,a2d,b2d,d2d,               !output
     &          nbulkl,nbulkr,ninprc,nintfc,ntotal,
     &          concl,concr,conc,
     &          idpotla,idpotlb,idpotra,idpotrb,idpota,idpotb,
     &          cvec,ctol,swsl,swsr,sws,smtl,smtr,smt)
c Laszlo 15/07/2003
        do i=1,ntotal
          iesubl(i)=i
        end do
c
c       Set layer quantities specific to rsp-calculation. Defaults are:
c       - uniform orientation of magnetic field
c       - uniform SOC scaling 
c
        do i=1,ninprc(0)
          sxl(i) = sx0(0)
          call bvec(btha,bpha,fact,rbl(1,i))
        enddo
c
        do i=1,nintfc
          sxa(i) = sx0(0)
          sxb(i) = sx0(0)
          call bvec(btha,bpha,fact,rba(1,i))
          do j=1,3
            rbb(j,i) = rba(j,i)
          enddo
        enddo
c
        do i=1,ninprc(nprc+1)
          sxr(i) = sx0(0)
          call bvec(btha,bpha,fact,rbr(1,i))
        enddo
c
c
c***********************************
c***********************************
c Option 3: read in structure
c***********************************
c***********************************
c
      elseif(igraph.lt.0) then
c
c=============================
c      Define 2D symmetry
c=============================
c
c     -number of P.L. (nprc) and extra layers (nextra)
c
c     -input for largest 2D cell common to all layers
c      a,beta,delta defined according to Cunningham
c      aunit gives conversion factor to a.u.
c
c      Bravais lattices: a1=a*(1,0)
c
c       (i2dlat=1)    a2=a*(delta,beta)      oblique
c       (i2dlat=2)    a2=a*(1/2,beta/2)      cent.rectangular
c       (i2dlat=3)    a2=a*(0,beta)          prim.rectangular
c       (i2dlat=4)    a2=a*(0,1)             square
c       (i2dlat=5)    a2=a*(1/2,sqrt(3)/2)   hexagonal
c
c
       read(5,*) nprc,nextra
       read(5,'(a4)') idgroup
       read(5,*) i2dlat,a2d
       read(5,*) b2d,d2d
       if(iunit.eq.0) a2d = a2d/aunit
c
c -specify layers for which orientation of magnetization is kept fixed
       read(5,*) lfix1,lfix2
       if(lfix1.ge.1) then
        if(lfix2.ge.1) then
          write(6,'('' Keep orientation of magnetization fixed '',
     >    ''at layers:'',i3,i5)') lfix1,lfix2
        else
          write(6,'('' Keep orientation of magnetization fixed '',
     >    ''at layer:'',i3)') lfix1
        end if
       else
        if(lfix2.ge.1) then
          lfix1=lfix2
          lfix2=0
          write(6,'('' Keep orientation of magnetization fixed '',
     >    ''at layer:'',i3)') lfix1
        end if
       end if
c
c      Optional input for a linear Bloch Wall 
       linbw=.false.
c      read(5,*) linbw,phi1,phi2,lay1,lay2,shbw
c
c=================================================
c      Read parameters defining principal layers
c=================================================
c
c      INTERMEDIATE-region
c
       read(5,*) (ninprc(i),i=1,nprc)
       nintfc = 0
       do i=1,nprc
         nintfc = nintfc + ninprc(i)
       enddo
c
c      LEFT,RIGHT-regions
c
       read(5,*) ninprc(0),nbulkl
       read(5,*) ninprc(nprc+1),nbulkr 
c     
       ntotal=nintfc+ninprc(0)*(nextra+1)+ninprc(nprc+1)*(nextra+1)
c
c====================================================
c      Layer-dependent quantities
c====================================================
c
c
c--    Read in for ntotal layers: 
c      - PL and layer identification (dummy)
c      - CPA concentration (cross-checked in pothandle)
c      - atomic site vector (all should be referred to a common origin)
c      - Wigner-Seitz radius and scaling for muffin-tin radius
c      - orientation of magnetic field
c      - scaling factors for SOC
c      - identifier of chemical species (first 3 characters)
c
c
       read(5,*) ctol      ! tolerance for checks of cvec 
       read(5,*) runit     ! cvec and radii scale factor (iunit-units)
c
       if(iunit.eq.0) runit=runit/aunit
c
c
        ill=0
        ili=0
        ilr=0
        lmin=ninprc(0)*nextra
        lmax=ntotal-ninprc(nprc+1)*nextra
        lleft = 1
        lint  = ninprc(0)*(nextra+1)+1
        lright= lint+nintfc
c
        read(5,'(a40)') tregion       ! dummy information
        read(5,'(a40)') tregion       ! dummy information
        read(5,'(a40)') tregion       ! dummy information
        do il=1,ntotal
          if((il.eq.lleft).or.(il.eq.lint).or.(il.eq.lright)) 
     &              read(5,'(a40)') tregion
          read(5,*) plinf,linf,
     &    conc0(1),(cvec(il,i),i=1,3),sws0(1),smt0(1),
c Laszlo  8/07/2004
     &    btha,bpha,bthb,bphb,isba,isbb,sx0(0),sx0(1),iesubl(il)
c
c    &    btha,bpha,bthb,bphb,isba,isbb,sx0(0),sx0(1)
c         iesubl(il)=il
c
          do i=1,3
            cvec(il,i) = runit*cvec(il,i)
          enddo
          sws0(1)=runit*sws0(1)
          smt0(1)=runit*smt0(1)
c
c         Skip extra layers.
          if((il.gt.lmin).and.(il.le.lmax)) then
c
c======
c LEFT
c======
            if(il.le.(lmin+ninprc(0))) then
              ill=ill+1
              concl(ill)= conc0(1)
              swsl(ill) = sws0(1)
              smtl(ill) = smt0(1)
              sxl(ill)  = sx0(0)
              call bvec(btha,bpha,fact,rbl(1,ill))
c
c             Some checks.
              difs = dabs(sx0(0)-sx0(1))
              dift = dabs(btha-bthb)
              difp = dabs(bpha-bphb)
              if((difs.gt.1.0d-4).or.(dift.gt.1.0d-4).or.
     &          (difp.gt.1.0d-4)) write(6,*) 
     &          'WARNING!!! Left-side should not be CPA dependent'
c=======
c RIGHT
c=======
            elseif(il.gt.(nintfc+lmin+ninprc(0))) then
              ilr=ilr+1
              concr(ilr)= conc0(1)
              swsr(ilr) = sws0(1)
              smtr(ilr) = smt0(1)
              sxr(ilr)  = sx0(0)
              call bvec(btha,bpha,fact,rbr(1,ilr))
c
c             Some checks.
              difs = dabs(sx0(0)-sx0(1))
              dift = dabs(btha-bthb)
              difp = dabs(bpha-bphb)
              if((difs.gt.1.0d-4).or.(dift.gt.1.0d-4).or.
     &          (difp.gt.1.0d-4))  write(6,*) 
     &          'WARNING!!! Right-side should not be CPA dependent'
c==============
c INTERMEDIATE
c==============
            else
              ili=ili+1
              conc(ili)= conc0(1)
              sws(ili) = sws0(1)
              smt(ili) = smt0(1)
              sxa(ili) = sx0(0)
              iesublatt(ili)=iesubl(il)
              call bvec(btha,bpha,fact,rba(1,ili))
!             isigba(ili)=isba
c
              if(1.d0-conc0(1).lt.1.d-6) then
                call bvec(btha,bpha,fact,rbb(1,ili))
!               isigbb(ili)=isigba(ili)
                sxb(ili) = sxa(ili)
c               Some checks.
                difs = dabs(sx0(0)-sx0(1))
                dift = dabs(btha-bthb)
                difp = dabs(bpha-bphb)
                if((difs.gt.1.0d-4).or.(dift.gt.1.0d-4).or.
     &            (difp.gt.1.0d-4))  write(6,*) 
     &            'WARNING!!! I-region is not CPA dependent'
c
              else
                call bvec(bthb,bphb,fact,rbb(1,ili))
!               isigbb(ili)=isbb
                sxb(ili) = sx0(1)
              endif
            endif
          endif
        enddo
c***********************************
c***********************************
      endif
c***********************************
c***********************************
c
c
c============================================================
c  Preliminary checks
c============================================================
      if(nprc.gt.mprc) stop 'Increase mprc in param.h'
      if(nextra.gt.mextra) stop 'Increase mextra in param.h'
c
      do i=0,nprc+1
        write(6,*) 'minprc, ninprc(',i,'):',minprc,' ',ninprc(i)
        if(ninprc(i).gt.minprc) stop 'Increase minprc in param.h'
      enddo
      if(ntotal.gt.mtotal) stop 'Increase mtotal in param.h'
      if(nintfc.gt.mintfc) stop 'Increase mintfc in param.h'
      if(nintfc.lt.ninprc(0).or.nintfc.lt.ninprc(nprc+1)) then
        write(6,'(''WARNING!!''/,
     >''Revise dimension of auxiliary vectors in readgeom is big enough
     >'')')
      endif
c
      nunitl = ninprc(0)/nbulkl
      if(nunitl*nbulkl.ne.ninprc(0)) then
        write(6,'(/'' 
     >    Integer number of bulk units in L-region required!!!''/
     >    '' Present no. of layers'',i2,
     >    '' Present no. of layers in a bulk unit'',i2/)') 
     >    ninprc(0),nbulkl
        stop
      end if        
      if(bulk) then
        nunitr=nunitl
        if(nbulkl.ne.nbulkr) stop 
     >         'nbulk mismatch between Left and Right!!'
        if(nprc.ne.1) stop 
     >         'nprc should be 1!!'
        if(ninprc(1).ne.ninprc(0).or.ninprc(1).ne.ninprc(nprc+1)) stop
     >         'ninprc mismatch between Left-, Right-, I-regions!!'
      else
        nunitr = ninprc(nprc+1)/nbulkr
        if(nunitr*nbulkr.ne.ninprc(nprc+1)) then
          write(6,'(/'' 
     >      Integer number of bulk units in R-region required!!!''/
     >      '' Present no. of layers'',i2,
     >      '' Present no. of layers in a bulk unit'',i2/)') 
     >      ninprc(nprc+1),nbulkr
          stop
        end if        
      endif
c     
c     Bloch Wall 
      if(linbw) then
            phi1=phi1*pi/180.d0
            phi2=phi2*pi/180.d0
            if(lay1.lt.1) 
     >      stop 'linear Bloch Wall: lay1 less than 1'
            if(lay2.gt.nintfc) 
     >      stop 'linear Bloch Wall: lay2 larger than nintfc'
      end if
c 
c
c============================================================
c  Generate 2D direct and reciprocal lattice
c============================================================
c
      call primkr2d(i2dlat,idgroup,a2d,b2d,d2d,arot)
      call genlatt2d
c
c
c================================================================
c  Set up Bravais matrices for left and right semi-infinite bulks
c  Generate corresponding 3D lattices
c================================================================
c
      do i=1,2
        a1l(i)=a1(i)
        a2l(i)=a2(i)
        a1r(i)=a1(i)
        a2r(i)=a2(i)
      end do
      a1l(3)=0.d0
      a2l(3)=0.d0
      a1r(3)=0.d0
      a2r(3)=0.d0
      do i=1,3
        a3l(i)=cvec(nbulkl+1,i)-cvec(1,i)
        a3r(i)=cvec(ntotal,i)-cvec(ntotal-nbulkr,i)
      end do
      call genlatt3d
      if(igraph.ge.0) then
        if(dabs(sws3d-swsl(1)).ge.tol) stop
     &   'sws mismatch!! Revise construction of 3D lattice'
      endif
c
c Laszlo 16/07/2003
      volws=0.d0
      do li=1,nintfc
        volws=volws+(4.0d0*pi/3.0d0)*sws(li)**3
      end do
      volint=vol*(cvec(lright,3)-cvec(lint,3))
      dvol=dabs((volint-volws)/volws)
      write(6,'(/'' Volumen of interface region '',f15.6)') volint
      write(6,'('' Volumen of WS spheres       '',f15.6)') volws 
      write(6,'('' Volumen mismatch             '',d15.6)') dvol
      if(dvol.gt.ctol*10.d0) stop ' Real-space volumen mismatch'
c
      call flush(6)
c
c
c==============================================
c  Some manipulations on vector directions.
c==============================================
c   Move all cvec within the smallest 2-D column
c   (only for abs(igraph).gt.2)
c
      call madcheck(ctol,ntotal,igraph)
c
c   Rotate vectors to match implemented point-group operations
c 
      if(dabs(arot-0.d0).gt.tol) then
        call sphcheck(arot,rbl,rbr,rba,rbb,nintfc)
      endif
c
c   In case of a linear Bloch Wall set input directions for magnetization
      if(linbw) then
c
        do il=1,ninprc(0)
          rbl(1,il)=dcos(phi1)
          rbl(2,il)=dsin(phi1)
          rbl(3,il)=0.0d0        
        end do
        do il=1,ninprc(nprc+1)
          rbr(1,il)=dcos(phi2)
          rbr(2,il)=dsin(phi2)
          rbr(3,il)=0.0d0        
        end do
c
        shbw=shbw*a3l(3)
        laymin=ninprc(0)*(nextra+1)+lay1
        laymax=ninprc(0)*(nextra+1)+lay2
        width=cvec(laymax,3)-cvec(laymin,3)
        do il=1,nintfc
          ill=il+ninprc(0)*(nextra+1)
          d=cvec(ill,3)-cvec(laymin,3)
          phi=phi1+(phi2-phi1)*(d-shbw)/width
          if(phi.lt.phi1) phi=phi1
          if(phi.gt.phi2) phi=phi2
          phibw(il)=phi
          rba(1,il)=dcos(phi)
          rba(2,il)=dsin(phi)
          rba(3,il)=0.0d0        
          rbb(1,il)=dcos(phi)
          rbb(2,il)=dsin(phi)
          rbb(3,il)=0.0d0        
        end do
c
      end if
c
ccc
c============================================================
c  Write summary of structure
c============================================================
c
      write(6,'(/,/'' ======================================= '')')
      write(6,'(''               2D  STRUCTURE '')')
      write(6,'('' ======================================= '',/)')
      write(6,'('' Lattice constant:'',f15.8)') a2d
      write(6,'('' Asymmetry ratios'',2f15.8)') b2d,d2d             
      write(6,'('' cvec tolerance  '',e10.2)') ctol
c
      write(6,'(/'' No. of principal layers: '',i3)') nprc
      write(6,'('' No. of layers in each principal layer:'',
     >          100i2)') (ninprc(i),i=1,nprc)
      write(6,'('' Total no. of layers in the I-region: '',i3)') nintfc
      write(6,'(/'' ntotal='',i4,''  with'')') ntotal
      write(6,'(i3,'' extra layers added for screening'')') nextra
      write(6,'(i3,'' layers in a Left  principal layer'')')  ninprc(0)
      write(6,'(i3,'' layers in a Right principal layer'')')  
     >         ninprc(nprc+1)
      write(6,'(/'' No. of layers in a Left/Right bulk unit='',2i2)')  
     >         nbulkl,nbulkr
c
c
c     Write out layer by layer quantities.
      write(6,'(/,/,
     &''PL Lay  conc  cvec(x,y,z)  sws  smt  B(x,y,z)  SOC  '')')
      write(6,'(
     &''==========================================================='')')
c
      write(6,'(''   Left region              '')')
      write(6,'('' =============== '')')
      ipl=1
      do il=1,ninprc(0)
        ill=il+nextra*ninprc(0)
        write(6,'(i3,i4,f7.3,2x,3f10.5,2x,2f7.4,2x,3f7.3,f8.3,2x,i4)')
     &   ipl,il,concl(il),(cvec(ill,i),i=1,3),swsl(il),smtl(il),
     &   (rbl(j,il),j=1,3),sxl(il),iesubl(ill)
      enddo
c
      write(6,'(/''   Intermediate region              '')')
      write(6,'('' ======================= '')')
      ill=0
      ipl=1
      do il=1,nintfc
        ill = ill+1
        ilc = il+(nextra+1)*ninprc(0)
        if(linbw) then
        write(6,'(i3,i4,f7.3,2x,3f10.5,2x,2f7.4,2x,f7.3,
     &   2x,3f7.3,f8.3,2x,i4)')
     &   ipl,il,conc(il),(cvec(ilc,i),i=1,3),sws(il),smt(il),
     &   phibw(il)/fact,(rba(j,il),j=1,3),sxa(il),iesublatt(il)
        else
        write(6,'(i3,i4,f7.3,2x,3f10.5,2x,2f7.4,2x,3f7.3,f8.3,2x,i4)')
     &   ipl,il,conc(il),(cvec(ilc,i),i=1,3),sws(il),smt(il),
     &   (rba(j,il),j=1,3),sxa(il),iesublatt(il)
        end if
        if(1.d0-conc(il).gt.1.d-6) then
          write(6,'(7x,f7.3,2x,30x,2x,14x,2x,3f7.3,f8.3)') 
     &     1.d0-conc(il),(rbb(i,il),i=1,3),sxb(il)
        endif
        if(ill.eq.ninprc(ipl)) then
          ipl = ipl+1
          ill = 0
        endif
      enddo
c
      write(6,'(/''   Right region              '')')
      write(6,'('' ================ '')')
      ipl=1
      do il=1,ninprc(nprc+1)
        ill=il+(nextra+1)*ninprc(0)+nintfc
        write(6,'(i3,i4,f7.3,2x,3f10.5,2x,2f7.4,2x,3f7.3,f8.3,2x,i4)')
     &   ipl,il,concr(il),(cvec(ill,i),i=1,3),swsr(il),smtr(il),
     &   (rbr(j,il),j=1,3),sxr(il),iesubl(ill)
      enddo
c
c***********************************************************************
c                            FINAL CHECKS
c***********************************************************************
c------------
ccc   Checks for geometry (further graphics).
c
      open(unit=11,file='readgeom.out',status='unknown')
c
      if(igraph.ne.0) then
        do il=1,ntotal
          write(11,'(i3,3f16.9)') il,(cvec(il,i),i=1,3)
        enddo
        write(11,'(/,''Please, check manually that left-most I-layers'',
     &             /,''and right-most R-layers are bulk-like,'',
     &             /,''(automated version not yet ready)'')')
c
c     In case igraph=0, generate input_geo file and stop.
      else
        write(11,'(2i3,20x,''    igraph, iunit'')') -1,1
        write(11,'(2i3,20x,''    nprc, nextra'')') nprc,nextra
        write(11,'(a4,21x,''     idgroup '')') idgroup
        write(11,'(i4,f11.7,16x,''  i2dlat, a2d '')') i2dlat,a2d
        write(11,'(2f7.4,14x,''  b2d,d2d '')') b2d,d2d
        if(linbw) write(11,
     &    '(''.true.'',2f7.2,2i3,10x,''  linbw, phi1-2, lay1-2'')') 
     &    phi1*180.d0/pi,phi2*180.d0/pi,lay1,lay2
        if(.not.linbw) write(11,
     &    '(''.false.  0.  0.  1  1         linbw, phi1-2, lay1-2'')') 
        write(11,*) (ninprc(i),i=1,nprc),'             ninprc'
        write(11,'(2i3,20x,''    ninprcl, nbulkl'')') ninprc(0),nbulkl
        write(11,'(2i3,20x,''    ninprcr, nbulkr'')') ninprc(nprc+1),
     &  nbulkr
        write(11,'(e7.1,20x,''   ctol'')') ctol
        write(11,'(f7.4,20x,''   runit'')') 1.d0
c
        write(11,'(''      =================  '')')
        write(11,
     & '(''PL lay conc cvec(x,y,z) sws smt B(th,phi) SOC idpot'')')
        write(11,'(''      =================  '')')
c
        write(11,'(''          * L *          '')')
        ninprcl = ninprc(0)*(nextra+1)
        ipl=1
        il=1
        do ill=1,ninprcl
          write(11,'(i3,i4,f7.3,2x,3f13.8,2x,2f8.5,2x,4f7.2,2x,
     &      2f7.2,2a4)')
     &      ipl,il,concl(il),(cvec(ill,i),i=1,3),swsl(il),smtl(il),
     &      btha,bpha,btha,bpha,sxl(il),sxl(il),idpotla(il)(1:3),
     &      idpotlb(il)(1:3)
          il=il+1
          if(il.gt.ninprc(0)) then
            il=1
            ipl=ipl+1
          endif
        enddo
c
        write(11,'(''          * I *          '')')
        ninprcl = ninprc(0)*(nextra+1)
        ipl=1
        lpl=1
        do il=1,nintfc
          ill = il+ninprcl
          lpl = lpl+1
          write(11,'(i3,i4,f7.3,2x,3f13.8,2x,2f8.5,2x,4f7.2,2x,
     &      2f7.2,2a4)')
     &      ipl,il,conc(il),(cvec(ill,i),i=1,3),sws(il),smt(il),
     &      btha,bpha,btha,bpha,sxa(il),sxb(il),idpota(il)(1:3),
     &      idpotb(il)(1:3)
          if(lpl.gt.ninprc(ipl)) then
            lpl=1
            ipl=ipl+1
          endif
        enddo
c
        write(11,'(''          * R *          '')')
        ninprcr = ninprcl+nintfc+1
        ipl=1
        il=1
        do ill=ninprcr,ntotal
          write(11,'(i3,i4,f7.3,2x,3f13.8,2x,2f8.5,2x,4f7.2,2x,
     &      2f7.2,2a4)')
     &      ipl,il,concr(il),(cvec(ill,i),i=1,3),swsr(il),smtr(il),
     &      btha,bpha,btha,bpha,sxr(il),sxr(il),idpotra(il)(1:3),
     &      idpotrb(il)(1:3)
          il=il+1
          if(il.gt.ninprc(nprc+1)) then
            il=1
            ipl=ipl+1
          endif
        enddo
        stop
c
      endif
      close(11)
c
c
c------------
ccc   Checks for ISBZ integration
c----------------------------------------------------------------------
c     The a-matrices for point-group operations are stored for further
c     usage in kmesh2d.
c----------------------------------------------------------------------
c
c     Point-group symmetry operations corresponding to idgroup.
      call rgroup(i2dlat,idgroup,ng0,a0mat,melem)
c
c     Point-group symmetry operations corresponding to i2dlat.
c     (this means the minimum possible IBZ of the crystal group).
      imgroup = igroup(i2dlat)
      call rgroup(i2dlat,imgroup,ng,amat,melem)
c
      if(intbz.eq.1) then
c
        nmin=nextra*ninprc(0)+1
        nmax=ntotal-nextra*ninprc(nprc+1)
        do 20 il=nmin,nmax
c Laszlo 15/07/2003
          do i=1,2
            rlay0(i) = cvec(il,i)
          enddo
c         write(6,'(//''il='',i3,''  subl='',i3,''  r='',2f14.5)') 
c    &    il,iesubl(il),rlay0(1),rlay0(2)
          do ig=1,ng0
            do i=1,2
              rd0(i) =
     &        a0mat(i,1,ig)*rlay0(1) +a0mat(i,2,ig)*rlay0(2)
            enddo
c           write(6,'(/'' ig='',i4,''  rr='',2f14.5)') ig,rd0(1),rd0(2)
c Laszlo 15/07/2003
            eqsubl=.false.
            do 30 jl=nmin,nmax
              if(iesubl(il).ne.iesubl(jl)) goto 30
              do i=1,2
                rd(i)=rd0(i)-cvec(jl,i)
              enddo
              rdn=(a2(2)*rd(1)-a2(1)*rd(2))/
     &            (a2(2)*a1(1)-a2(1)*a1(2))
              rdm=(a1(2)*rd(1)-a1(1)*rd(2))/
     &            (a2(1)*a1(2)-a2(2)*a1(1))
              nrd=nint(rdn)
              mrd=nint(rdm)
              dnrd=dabs(float(nrd)-rdn)
              dmrd=dabs(float(mrd)-rdm)
c             write(6,'(/''jl='',i3,''  subl='',i3,''  r='',2f14.5)') 
c    &        jl,iesubl(jl),cvec(jl,1),cvec(jl,2)
c             write(6,'(2f14.5)') rd
c             write(6,'(2d12.5)') dnrd,dmrd
              if(dnrd.lt.ctol.and.dmrd.lt.ctol) then
                eqsubl=.true.
                goto 35
              end if
 30         continue
 35         if(.not.eqsubl) then
              write(6,'(/'' WARNING!!'',i4,''  ig'',i3)') il,ig
              write(6,'(
     &        ''   Irreducible representation of cvec not valid'',/,
     &        ''   Change to full BZ integration'')')
              intbz=2
            endif
          enddo
 20     continue
      endif
c
c ===============================
c Check about ideal bulk geometry
c ===============================
c
c     in total there are nprc+2+2*nextra PL's
c     the size of the first nextra+1 PL's is ninprc(0)
c     the size of the last  nextra+1 PL's is ninprc(nprc+1)
c     the size of the PL no. iprc (iprc=1,nprc) is ninprc(iprc)
c
      do iprc=1,nprc+1
        if(ninprc(iprc).ne.ninprc(iprc-1)) goto 100
      end do
      ninprc0=ninprc(0)
      cshift(1)=cvec(ninprc0+1,1)-cvec(1,1)
      cshift(2)=cvec(ninprc0+1,2)-cvec(1,2)
      cshift(3)=cvec(ninprc0+1,3)-cvec(1,3)
      il=0
      do iprc=1,nprc+1+2*nextra
        do inprc=1,ninprc0
          il=il+1
          do i=1,3
            if(i.eq.2) then
c             write(6,*) 'bulkgeo test',il+2,i
c             write(6,*) dabs(cvec(il+3,i)-cvec(il+2,i))
            end if
            if(dabs(cvec(il+ninprc0,i)-cvec(il,i)-cshift(i)).gt.ctol) 
     >        then 
c             write(6,*) 'bulkgeo fails in il=',il
c             write(6,*) 'bulkgeo fails in i=',i
              goto 100
            end if
          end do
        end do
      end do
      goto 105
c
  100 bulkgeo=.false.
      goto 110
  105 bulkgeo=.true.
  110 continue
      if(bulkgeo) then
        write(6,
     >  '(/'' *************************************************'')')
        write(6,
     >  '('' Ideal bulk geometry applies over the whole system'')')
        write(6,
     >  '('' *************************************************'')')
      else
        write(6,'(/'' **********************************'')')
        write(6,'('' Ideal bulk geometry does not apply'')')
        write(6,'('' **********************************'')')
      end if
c
c------------
ccc   Checks for potentials
c
c     if(layshift2.gt.nintfc) stop
c    >  'Cannot deal with layshift bigger than nintfc!'
c     if(layb0f.gt.nintfc) stop
c    >  'Cannot deal with layb0 bigger than nintfc!'
c
      return
      end
c
c================================================================
c
      subroutine bvec(bth0,bph0,fact,rb0)
c
      implicit real*8(a-h,o-z)
      dimension rb0(3)
c
      bth=fact*bth0
      bph=fact*bph0
      rb0(1)=dsin(bth)*dcos(bph)
      rb0(2)=dsin(bth)*dsin(bph)
      rb0(3)=dcos(bth)
c
      return
      end
c=======================================================================
c
      subroutine sphcheck(arot,rbl,rbr,rba,rbb,nintfc)
c
      implicit real*8(a-h,o-z)
      include '../param.h'
c
      character*4 idgroup
      dimension rbl(3,minprc),rbr(3,minprc)
      dimension rba(3,mintfc),rbb(3,mintfc)
      dimension ctry(2)
c
      common/latt2d/idgroup,i2dlat,a2d,b2d,d2d
      common/lay2d/cvec(mtotal,3),nextra,nbulkl,nbulkr,
     &             nprc,ninprc(0:mprc+1)
c
c       Some checks
        if(i2dlat.eq.1) stop 
     &  'Oblique lattice should not need any change of idgroup'
        if(i2dlat.eq.2.or.i2dlat.eq.3) then
          if(idgroup.ne.'Cs') stop
     &    'Rectangular lattice may need to change idgroup just for Cs'
        endif
        if(i2dlat.eq.4) then
          if(idgroup.ne.'Cs'.and.idgroup.ne.'C2v') stop
     &    'Square lattice may need to change idgroup just for Cs, C2v'
        endif
        if(i2dlat.eq.5) then
          if(idgroup.ne.'Cs'.and.idgroup.ne.'C2v') stop
     &  'Hexagonal lattice may need to change idgroup just for Cs,C2v'
        endif
c
        ntotal=nintfc+ninprc(0)*(nextra+1)+ninprc(nprc+1)*(nextra+1)
        do il=1,ntotal
          ctry(1)=cvec(il,1)
          ctry(2)=cvec(il,2)
          call vecrot(ctry,arot)
          cvec(il,1)=ctry(1)
          cvec(il,2)=ctry(2)
        enddo
c
        do il=1,ninprc(0)
          ctry(1)=rbl(1,il)
          ctry(2)=rbl(2,il)
          call vecrot(ctry,arot)
          rbl(1,il)=ctry(1)
          rbl(2,il)=ctry(2)
        end do
        do il=1,nintfc
          ctry(1)=rba(1,il)
          ctry(2)=rba(2,il)
          call vecrot(ctry,arot)
          rba(1,il)=ctry(1)
          rba(2,il)=ctry(2)
          ctry(1)=rbb(1,il)
          ctry(2)=rbb(2,il)
          call vecrot(ctry,arot)
          rbb(1,il)=ctry(1)
          rbb(2,il)=ctry(2)
        end do
        do il=1,ninprc(nprc+1)
          ctry(1)=rbr(1,il)
          ctry(2)=rbr(2,il)
          call vecrot(ctry,arot)
          rbr(1,il)=ctry(1)
          rbr(2,il)=ctry(2)
        end do 
c
      call flush(6)
      return
      end
