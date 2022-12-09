c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
c*****************************************************************
c Assumed that the left-(right-)most nbulk layers are bulk-like.
c*****************************************************************
      subroutine madelungb(qmom,nintfc,lmax,sigma0,bulk,side,
     >                     vmad,vleft,vright)
c=========================
c
c Calculate layer dependent Madelung potential for a bulk with complex lattice
c
c input:  
c         qmom - moments of the charge density
c         nintfc - number of layers
c         lmax - l-cutoff
c         sigma0 - dimensionless Ewald parameter
c         bulk - if .true., scf bulk calculation
c         side - if not bulk, 'L'eft or 'R'ight
c output: vmad - layer resolved Madelung potentials
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
      logical bulk
c
      character*1 side
c
      dimension cbulk(3,mbulk),vbulk(mbulk),vmad(mintfc)
c
      complex*16 qmom(lmsup,minprc),qmombulk(lmsup,mbulk)
      complex*16 sum,conmad(lmsup,lmsup,mbulk,mbulk)
c
      common/brav3dl/a1l(3),a2l(3),a3l(3)
      common/dirvec3dl/rvecl(3,mdir),volwsl,distrl,maxrl
      common/recvec3dl/veckl(3,mrec),volbzl,distkl,maxkl
c
      common/brav3dr/a1r(3),a2r(3),a3r(3)
      common/dirvec3dr/rvecr(3,mdir),volwsr,distrr,maxrr
      common/recvec3dr/veckr(3,mrec),volbzr,distkr,maxkr
c
      common/lay2d/clay(mtotal,3),nextra,nbulkl,nbulkr,
     &             nprc,ninprc(0:mprc+1) 
      common/crystl2d/vol2d,volbz2d
      common/madbulk/madvers
c
      common/test/itest
      common/madelung_energy/emad
      data tol/1.0d-12/
c----------------------------------------------------------------
c
      scale=dsqrt(vol2d)
      sigma=sigma0*scale
      lmaxs=2*lmax
      lmmaxs=(lmaxs+1)*(lmaxs+1)
c
      if(bulk.or.(side.eq.'L')) then
        ninprcb=ninprc(0)
        nbulk=nbulkl
c index of topmost layer in the left bulk
        nshift=(nextra+1)*ninprc(0)
      else if(side.eq.'R') then
        ninprcb=ninprc(nprc+1)
        nbulk=nbulkr
c index of layer lying by 'nbulk' layers below the topmost layer
c in the intermediate region 
        nshift=(nextra+1)*ninprc(0)+nintfc-nbulk
      end if
c
      nunits=ninprcb/nbulk
c
c Assign charges, multipole moments and non-primitive translation vectors 
c for a bulk unit
c
      do ibulk=1,nbulk
c
        do l=1,lmmaxs
          qmombulk(l,ibulk)=(0.0d0,0.0d0)
        end do
        do iunit=1,nunits
          ilay=(iunit-1)*nbulk+ibulk
          do l=1,lmmaxs
            qmombulk(l,ibulk)=qmombulk(l,ibulk)+qmom(l,ilay)
          end do
        end do
        do l=1,lmmaxs
          qmombulk(l,ibulk)=qmombulk(l,ibulk)/nunits
        end do
c
        do i=1,3
          cbulk(i,ibulk)=clay(nshift+ibulk,i)-clay(nshift+1,i)
        end do
c
      end do
c
c
c Calculate Madelung constants
c
      if(madvers.eq.1.or.madvers.eq.3) then
      call czero(conmad,lmsup*lmsup*mbulk*mbulk)
      if(bulk.or.side.eq.'L') then
         call bulkmad1(nbulk,mbulk,lmax,sigma,volwsl,cbulk,
     >                 maxkl,veckl,distkl,maxrl,rvecl,distrl,
     >                 conmad)
      else if(side.eq.'R') then
         call bulkmad1(nbulk,mbulk,lmax,sigma,volwsr,cbulk,
     >                 maxkr,veckr,distkr,maxrr,rvecr,distrr,
     >                 conmad)
      end if
      end if
c
      if(madvers.eq.3) then
      write(6,'(/'' Madelung constants 1'')')
      do i=1,nbulk
      do j=1,nbulk
        write(6,'('' Atoms'',2i4)') i,j
        call outmat1(conmad(1,1,i,j),lmmaxs,lmmaxs,lmsup,tol,6)
      end do
      end do
      end if
c
      if(madvers.eq.2.or.madvers.eq.3) then
      call czero(conmad,lmsup*lmsup*mbulk*mbulk)
      if(bulk.or.side.eq.'L') then
         call bulkmad2(nbulk,lmax,sigma,vol2d,cbulk,a3l,conmad)
      else if(side.eq.'R') then
         call bulkmad2(nbulk,lmax,sigma,vol2d,cbulk,a3r,conmad)
      end if
      end if
c
      if(madvers.ge.3) then
      write(6,'(/'' Madelung constants 2'')')
      do i=1,nbulk
      do j=1,nbulk
        write(6,'('' Atoms'',2i4)') i,j
        call outmat1(conmad(1,1,i,j),lmmaxs,lmmaxs,lmsup,tol,6)
      end do
      end do
      end if
c
c Calculate Madelung potentials end energy
c
      write(6,'(/'' Madelung potentials for bulk'')')
      emad_bulk=0.0d0
      do i=1,nbulk
        sum=(0.d0,0d0)
        do j=1,nbulk
        do l=1,lmmaxs
          sum=sum+conmad(1,l,i,j)*qmombulk(l,j)
          do lp=1,lmmaxs
            emad_bulk=emad_bulk+
     >      dconjg(qmombulk(lp,i))*conmad(lp,l,i,j)*qmombulk(l,j)
          end do
        end do
        end do
        write(6,'('' Sublattice'',i2,'' :'', 2d17.7)') i,sum
        vbulk(i)=sum
      end do
      emad_bulk=emad_bulk/2.0d0
c
c Assign Madelung potentials to all the bulk-like layers involved in
c the calculation
c
      if(bulk) then
        do ibulk=1,nbulk
        do iunit=1,nunits
          ilay=(iunit-1)*nbulk+ibulk
          vmad(ilay)=vbulk(ibulk)
        end do
        end do
        emad=nunits*emad_bulk
      else if(side.eq.'L') then
          vleft=vbulk(nbulk)
      else if(side.eq.'R') then
          vright=vbulk(1)
      end if
c       
      return
      end
