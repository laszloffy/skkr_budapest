c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine struc3d(
     &          lat3d,lface,a3d,nprc,nextra,              !input
     &          mintfc,minprc,mprc,mtotal,
     &          i2dlat,idgroup,a2d,b2d,d2d,               !output
     &          nbulkl,nbulkr,ninprc,nintfc,ntotal,
     &          concl,concr,conc,
     &          idpotla,idpotlb,idpotra,idpotrb,idpota,idpotb,
     &          cvec,ctol,swsl,swsr,sws,smtl,smtr,smt)
c
c====================================================================
c Generate 2D structure and parameters from semiinfinite 3D-bulk.
c Default values:
c - conc = 1 for all layers (later changed in pothandle)
c - common Wigner-Seitz and muffin-tin radius for all atoms
c - idpot only accounts for chemical stacking sequence, but does
c   NOT correspond to the actual chemical elements (set in pothandle)
c
c Only ideal HCP lattices considered (see cratio).
c====================================================================
c
      implicit real*8 (a-h,o-z)
      parameter(mbulk=4)
c
      character*3 idchem0(mbulk)
      dimension rlay0(mbulk,3)
c
      character*4 idgroup
      character*5 lat3d
      character*10 idpotla(minprc),idpotra(minprc),idpota(mintfc)
      character*10 idpotlb(minprc),idpotrb(minprc),idpotb(mintfc)
      dimension conc(mintfc),concl(minprc),concr(minprc)
      dimension swsl(minprc),sws(mintfc),swsr(minprc)
      dimension smtl(minprc),smt(mintfc),smtr(minprc)
      dimension cvec(mtotal,3),ninprc(0:mprc+1)
c
      s2 = dsqrt(2.d0)
      s3 = dsqrt(3.d0)
      s6 = dsqrt(6.d0)
      pi = dacos(-1.d0)
c
c     Some default values.
      conc0 = 1.d0
      ctol = 1.0d-4
      b2d = 1.d0
      d2d = 0.d0
c
c ===========
c  SC
c ===========
c
      if(lat3d.eq.'SC') then
        sws3d = a3d*(3.d0/(4.d0*pi))**(1.d0/3.d0)
        smt3d = a3d/2.d0
        nbulk = 1
        idchem0(1) = ' A '
c
        if(lface.eq.100) then
          i2dlat = 4
          idgroup = 'C4v'
          a2d = a3d
          ninprc3d = 2
          rlay0(1,1) = 0.d0
          rlay0(1,2) = 0.d0
          rlay0(1,3) = a3d
c
        elseif(lface.eq.110) then
          i2dlat = 3
          idgroup = 'C2v'
          a2d = a3d*s2
          b2d = 1/s2
          ninprc3d = 2
          rlay0(1,1) = a3d/s2
          rlay0(1,2) = 0.d0
          rlay0(1,3) = a3d/s2
c
        elseif(lface.eq.111) then
          i2dlat = 5
          idgroup = 'C3vB'
          a2d = a3d*s2
          ninprc3d = 3
          rlay0(1,1) = 0.d0
          rlay0(1,2) = a3d*s2/s3
          rlay0(1,3) = a3d/s3
c
        else
          goto 10
        endif
c
c ===========
c  FCC
c ===========
c
      elseif(lat3d.eq.'FCC') then
        sws3d = a3d*0.5d0*(3.d0/(2.d0*pi))**(1.d0/3.d0)
        smt3d = a3d*0.5d0/s2
        nbulk = 1
        idchem0(1) = ' A '
c
        if(lface.eq.100) then
          i2dlat = 4
          idgroup = 'C4v'
          a2d = a3d/s2
          ninprc3d = 3
          rlay0(1,1) = a3d*0.5d0/s2
          rlay0(1,2) = a3d*0.5d0/s2
          rlay0(1,3) = a3d*0.5d0
c
        elseif(lface.eq.110) then
          i2dlat = 3
          idgroup = 'C2v'
          a2d = a3d
          b2d = 1.d0/s2
          ninprc3d = 3
          rlay0(1,1) = a3d*0.5d0
          rlay0(1,2) = a3d*0.5d0/s2
          rlay0(1,3) = a3d*0.5d0/s2
c
        elseif(lface.eq.111) then
          i2dlat = 5
          idgroup = 'C3vB'
          a2d = a3d/s2
          ninprc3d = 3
          rlay0(1,1) = 0.d0
          rlay0(1,2) = a3d/s6
          rlay0(1,3) = a3d/s3
c
        else
          goto 10
        endif
c
c ===========
c  BCC
c ===========
c
      elseif(lat3d.eq.'BCC') then
        sws3d = a3d*0.5d0*(3.d0/pi)**(1.d0/3.d0)
        smt3d = a3d*s3/4.d0
        nbulk = 1
        idchem0(1) = ' A '
c
        if(lface.eq.100) then
          i2dlat = 4
          idgroup = 'C4v'
          a2d = a3d
          ninprc3d = 3
          rlay0(1,1) = a3d/2.d0
          rlay0(1,2) = a3d/2.d0
          rlay0(1,3) = a3d/2.d0
c
        elseif(lface.eq.110) then
          i2dlat = 2
          idgroup = 'C2v'
          a2d = a3d*s2
          b2d = 1.d0/s2
          d2d = 0.5d0
          ninprc3d = 2
          rlay0(1,1) = a3d/s2
          rlay0(1,2) = 0.d0
          rlay0(1,3) = a3d/s2
c
        elseif(lface.eq.111) then
          i2dlat = 5
          idgroup = 'C3vB'
          a2d = a3d*s2
          ninprc3d = 4
          rlay0(1,1) = 0.d0
          rlay0(1,2) = a3d*s2/s3
          rlay0(1,3) = a3d*0.5d0/s3
c
        else
          goto 10
        endif
c
c ===========
c  HCP
c ===========
c
      elseif(lat3d.eq.'HCP') then
        cratio = dsqrt (8.d0/3.d0)                ! c/a
c
        sws3d = a3d*0.5d0*( (cratio*3.d0*s3)/(2.d0*pi) )**(1.d0/3.d0)
        smt3d = a3d/2.d0
        fac = 0.5d0*dsqrt( 4.d0/3.d0 + (cratio)**2 )
        if(fac.le.1.d0) smt3d = fac*smt3d
        nbulk = 2
        idchem0(1) = ' A '
        idchem0(2) = ' A '
c
        if(lface.eq.0001) then
          i2dlat = 5
          idgroup = 'C3vB'
          a2d = a3d
          ninprc3d = 2
          rlay0(1,1) = a3d*0.5d0
          rlay0(1,2) = a3d*0.5d0/s3
          rlay0(1,3) = a3d*cratio*0.5d0
          rlay0(2,1) =-a3d*0.5d0
          rlay0(2,2) =-a3d*0.5d0/s3
          rlay0(2,3) = a3d*cratio*0.5d0
c
        elseif(lface.eq.110) then
          i2dlat = 3
          idgroup = 'C2v'
          a2d = a3d*cratio
          b2d = 1.d0/cratio
          ninprc3d = 4
          rlay0(1,1) = a3d*cratio*0.5d0
          rlay0(1,2) = a3d*0.5d0
          rlay0(1,3) = a3d*0.5d0/s3
          rlay0(2,1) = a3d*cratio*0.5d0
          rlay0(2,2) = 0.d0
          rlay0(2,3) = a3d/s3
c
        else
          goto 10
        endif
c
c ===========
c  NaCl
c ===========
c
      elseif(lat3d.eq.'NaCl') then
        sws3d = a3d*0.5d0*(3.d0/(4.d0*pi))**(1.d0/3.d0)
        smt3d = a3d/4.d0
        nbulk = 2
        idchem0(1) = ' A '
        idchem0(2) = ' B '
c
        if(lface.eq.100) then
          i2dlat = 4
          idgroup = 'C4v'
          a2d = a3d/s2
          ninprc3d = 4
          rlay0(1,1) = a3d*0.5d0/s2
          rlay0(1,2) = a3d*0.5d0/s2
          rlay0(1,3) = 0.0d0
          rlay0(2,1) = 0.d0
          rlay0(2,2) = 0.d0
          rlay0(2,3) = a3d/2.d0
c
        elseif(lface.eq.110) then
          i2dlat = 3
          idgroup = 'C2v'
          a2d = a3d
          b2d = 1.d0/s2
          ninprc3d = 4
          rlay0(1,1) = a3d*0.5d0
          rlay0(1,2) = 0.0d0
          rlay0(1,3) = 0.0d0
          rlay0(2,1) = 0.0d0
          rlay0(2,2) = a3d*0.5d0/s2
          rlay0(2,3) = a3d*0.5d0/s2
c
        elseif(lface.eq.111) then
          i2dlat = 5
          idgroup = 'C3vB'
          a2d = a3d*0.5d0/s2
          ninprc3d = 4
          do ii=1,nbulk
            rlay0(ii,1) = 0.0d0
            rlay0(ii,2) = a3d*0.5d0/s6
            rlay0(ii,3) = a3d*0.5d0/s3
          enddo
c
        else
          goto 10
        endif
c
c ===========
c  Cu3Au
c ===========
c
      elseif(lat3d.eq.'Cu3Au') then
        sws3d = a3d*0.5d0*(3.d0/(2.d0*pi))**(1.d0/3.d0)
        smt3d = a3d*0.5d0/s2
        nbulk = 4
        idchem0(1) = ' A '
        idchem0(2) = ' B '
        idchem0(3) = ' B '
        idchem0(4) = ' B '
c
        if(lface.eq.100) then
          i2dlat = 4
          idgroup = 'C2vx'
          a2d = a3d
          ninprc3d = 8
          rlay0(1,1) = a3d*0.5d0
          rlay0(1,2) = a3d*0.5d0
          rlay0(1,3) = 0.0d0
          rlay0(2,1) = a3d*0.5d0
          rlay0(2,2) = 0.0d0
          rlay0(2,3) = a3d*0.5d0
          rlay0(3,1) =-a3d*0.5d0
          rlay0(3,2) =-a3d*0.5d0
          rlay0(3,3) = 0.0d0
          rlay0(4,1) =-a3d*0.5d0
          rlay0(4,2) = 0.0d0
          rlay0(4,3) = a3d*0.5d0
c
        elseif(lface.eq.110) then
          i2dlat = 3
          idgroup = 'Csx'
          a2d = a3d*s2
          b2d = 1.d0/s2
c         ninprc3d = 8
          ninprc3d = 4
          rlay0(1,1) = a3d/s2
          rlay0(1,2) = 0.0d0
          rlay0(1,3) = 0.0d0
          rlay0(2,1) = a3d*s2/4.d0
          rlay0(2,2) = a3d*0.5d0
          rlay0(2,3) = a3d*s2/4.d0
          rlay0(3,1) =-a3d/s2
          rlay0(3,2) = 0.0d0
          rlay0(3,3) = 0.0d0
          rlay0(4,1) =-a3d*s2/4.d0
          rlay0(4,2) =-a3d*0.5d0
          rlay0(4,3) = a3d*s2/4.d0
c
        elseif(lface.eq.111) then
          i2dlat = 5
          idgroup = 'C3vB'
          a2d = a3d*s2
c         ninprc3d = 12
          ninprc3d = 4 
          rlay0(1,1) = a3d/s2
          rlay0(1,2) = 0.0d0
          rlay0(1,3) = 0.0d0
          rlay0(2,1) = a3d*0.5d0/s2
          rlay0(2,2) = a3d*s3*0.5d0/s2
          rlay0(2,3) = 0.0d0
          rlay0(3,1) =-a3d/s2
          rlay0(3,2) = 0.0d0
          rlay0(3,3) = 0.0d0
          rlay0(4,1) = a3d*0.5d0/s2
          rlay0(4,2) =-a3d*0.5d0/s6
          rlay0(4,3) = a3d/s3
c
        else
          goto 10
        endif
c
c ===========
c  CsCl
c ===========
c
      elseif(lat3d.eq.'CsCl') then
        sws3d = a3d*0.5d0*(3.d0/pi)**(1.d0/3.d0)
        smt3d = a3d*s3/4.d0
        nbulk = 2
        idchem0(1) = ' A '
        idchem0(2) = ' B '
c
        if(lface.eq.100) then
          i2dlat = 4
          idgroup = 'C4v'
          a2d = a3d
          ninprc3d = 4
          do ii=1,nbulk
            rlay0(ii,1) = a3d*0.5d0
            rlay0(ii,2) = a3d*0.5d0
            rlay0(ii,3) = a3d*0.5d0
          enddo
c
        elseif(lface.eq.110) then
          i2dlat = 3
          idgroup = 'C2v'
          a2d = a3d*s2
          b2d = 1.d0/s2
          ninprc3d = 4
          rlay0(1,1) = a3d/s2
          rlay0(1,2) = a3d*0.5d0
          rlay0(1,3) = 0.0d0
          rlay0(2,1) = 0.0d0
          rlay0(2,2) = a3d*0.5d0
          rlay0(2,3) = a3d/s2
c
        elseif(lface.eq.111) then
          i2dlat = 5
          idgroup = 'C3vB'
          a2d = a3d*s2
          ninprc3d = 4
          do ii=1,nbulk
            rlay0(ii,1) = 0.0d0
            rlay0(ii,2) = a3d*s2/s3
            rlay0(ii,3) = a3d*0.5d0/s3
          enddo
c
        else
          goto 10
        endif
c
c ================
c  ZnS or Diamond
c ================
c
      elseif(lat3d.eq.'ZnS'.or.lat3d.eq.'Diam') then
        sws3d = a3d*0.5d0*(3.d0/(4.d0*pi))**(1.d0/3.d0)
        smt3d = a3d*s3/8.d0
        nbulk = 2
        idchem0(1) = ' A '
        idchem0(2) = ' B '
        if(lat3d.eq.'Diam') idchem0(2) = ' A '
c
        if(lface.eq.100) then
          i2dlat = 4
          idgroup = 'C2v'
          a2d = a3d*s2
          ninprc3d = 4
          rlay0(1,1) = 0.0d0
          rlay0(1,2) = a3d/s2
          rlay0(1,3) = a3d/4.0d0
          rlay0(2,1) = a3d/s2
          rlay0(2,2) = 0.0d0
          rlay0(2,3) = a3d/4.0d0
c
        elseif(lface.eq.110) then
          i2dlat = 3
          idgroup = 'Csx'
          a2d = a3d
          b2d = 1.d0/s2
          ninprc3d = 4
          rlay0(1,1) =-a3d/4.d0
          rlay0(1,2) = a3d*0.5d0/s2
          rlay0(1,3) = 0.0d0
          rlay0(2,1) =-a3d/4.d0
          rlay0(2,2) = 0.0d0
          rlay0(2,3) = a3d*0.5d0/s2
c
        elseif(lface.eq.111) then
          i2dlat = 5
          idgroup = 'C3vB'
          a2d = a3d/s2
          ninprc3d = 6
          rlay0(1,1) = 0.0d0
          rlay0(1,2) = a3d/s6
          rlay0(1,3) = a3d/(4.d0*s3)
          rlay0(2,1) = 0.0d0
          rlay0(2,2) = 0.0d0
          rlay0(2,3) = a3d*s3/4.d0
c
        else
          goto 10
        endif
c
c ================
c
      else
        write(6,
     &  '(''Revise 3D-lattice code. Only implemented:''/,
     &  ''SC, FCC, BCC, HCP, NaCl, Cu3Au, CsCl, ZnS, Diam'')')
        stop
      endif
c
c Now, generate layer by layer structure.
c
c     Check
      if(nbulk.gt.mbulk) stop 'Increase mbulk in struc3d'
c
      nintfc = 0
      do i=0,nprc+1
        ninprc(i) = ninprc3d
        if(i.ge.1.and.i.le.nprc) nintfc = nintfc + ninprc(i)
      enddo
      ntotal=nintfc+ninprc(0)*(nextra+1)+ninprc(nprc+1)*(nextra+1)
      nbulkl = nbulk
      nbulkr = nbulk
c
      il = 0
      inum = 1
      do j=1,nprc
        do i=1,ninprc3d
          il=il+1
          sws(il) = sws3d
          smt(il) = smt3d
          conc(il) = conc0
          idpota(il) = idchem0(inum)
          idpotb(il) = idpota(il)
c
          if(j.eq.1) then
            swsl(i) = sws3d
            smtl(i) = smt3d
            concl(i) = conc0
            idpotla(i) = idchem0(inum)
            idpotlb(i) = idpotla(i)
c
            swsr(i) = sws3d
            smtr(i) = smt3d
            concr(i) = conc0
            idpotra(i) = idchem0(inum)
            idpotrb(i) = idpotra(i)
          endif
c 
          inum = inum+1
          if(inum.gt.nbulk) inum = 1
        enddo
      enddo
c
      do i=1,3
        cvec(1,i) = 0.d0
      enddo
      inum = 1
      do j=2,ntotal
        do i=1,3
          cvec(j,i) = cvec(j-1,i) + rlay0(inum,i)
        enddo
        inum = inum+1
        if(inum.gt.nbulk) inum = 1
      enddo
c
      return
c
  10  write(6,
     & '(''Revise 3D-lattice orientation. Only implemented:''/,
     & ''100, 110, 111  and for HCP: 0001, 110'')')
      stop
      end
