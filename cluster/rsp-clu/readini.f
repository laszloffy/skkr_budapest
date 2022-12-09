c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine readini(
     >imesh,ne,npanel,nepanel,ne1,ne2,ne3,ebottom,etop,eps,lmax,
     >eta,sigma,!park,kset,ksetmax,itcpam,cpatol,cpatest,
     >for006,for007,for008,for009,!for010,iwrite,
!    >leftpot,rightpot,laypot,
     >laycore,!leftmat,rightmat,laymat,
     >wrel,lms,!bulk,
     >dos,!rightm,newvac,ivacpot,v0,v00,
     > vscreen,iscreen,
     >itscfmax,tolerr,tolef,tolen,!intbz,kunning,
     >orbpol,!opotl,opotr,
     >opot,lliter,
     >isscf,rotind,spinorb,dmmin)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
c MPI parameters
      integer root, myrank, nprocs
      common/mpi0/root,myrank,nprocs
c
      logical bulk,dos,wrel,lms,cpatest
      logical orbpol,opotl,opotr,opot
      character*25 job
      character*25 jobclu
      character*30 leftpot,rightpot,laypot,laycore
      character*30 leftmat,rightmat,laymat,leftmom,rightmom
      character*30 for006,for007,for008,for009,for010
      character*1 rightm
      dimension nepanel(5),ebottom(5),etop(5),eps(5)
      dimension ksettab(20),intab(20),kset(me)
      dimension park(2),r(3,minprc)
      dimension isigba(mintfc),isigbb(mintfc)
c Laszloffy 06/04/17
      logical isscf,rotind,spinorb
      real*8  dmmin
c
      common/test/itest
      common/iterpar/errmax,itermaxl,itermaxr,ichk
      common/broypot/cmix,wbr,nbrmax(100)
      common/broydir/dmix,wbrdir,nbrdir
      common/xch/ixch
      common/extf/qc,cpc
c     common/pot/potshift,b0,layshift1,layshift2,layb0,layb0f,
c    &           ib0,isigba,isigbb
      common/madbulk/madvers
c
c***************************************************************
c INPUT/OUTPUT files
c***************************************************************
c -initialize jobname
c
      read(5,'(a25)') jobclu
      read(5,'(a25)') job
c
c ******* open files according to job 
c
c     ---------------------------------------------------
      call jobnam(job,for006,for007,for008,for009,for010)
c     ---------------------------------------------------
      if (myrank == root) then
      open(6,file=for006,status='unknown')
      else
      open(6,file='/dev/null',status='unknown')
      end if
c
      write(6,*)
c
      write(6,'(2x,''routine READINI>''/)')
      call flush(6)
      write(6,'('' controll and test output is in file     : '',
     >          a)') for006
      write(6,'('' new potentials are in file              : '',
     >          a)') for007
      write(6,'('' densities of states are in file         : '',
     >          a)') for008
c     write(6,'('' effective t-matrices are in file        : '',
c    >          a)') for010
c
c -read in filenames for potentials and charge densities
c
      read(5,'(a)') leftpot
      read(5,'(a)') laypot
      read(5,'(a)') rightpot
      read(5,'(a)') laycore
      read(5,'(a)') leftmat
      read(5,'(a)') laymat
      read(5,'(a)') rightmat
      read(5,'(a)') leftmom
      read(5,'(a)') rightmom
c
      write(6,'('' potential for lhs bulk from             : '',
     >          a)') leftpot
      write(6,'('' potential for interfacial layers from   : '',
     >          a)') laypot
      write(6,'('' potential for rhs bulk from             : '',
     >          a)') rightpot
      write(6,'('' core charge densities for layers from   : '',
     >          a)') laycore
      write(6,'('' effective t-matrices for lhs bulk       : '',
     >          a)') leftmat
      write(6,'('' effective t-matrices for interface      : '',
     >          a)') laymat
      write(6,'('' effective t-matrices for rhs bulk       : '',
     >          a)') rightmat
c
c -test parameters
c
      read(5,*) itest,iwrite
c     (temporal)
      IF(IWRITE.NE.0) IWRITE=0
c
c
c***************************************************************
c Conditions of calculation
c***************************************************************
c
c -weak relativistic?
c -(lms) representation on output?
c
c     read(5,*) wrel,lms
      read(5,*) lms
      wrel=.false. ! 01/03/2021 Laszloffy
      write(6,'('' No weak relativistic approach (deleted)!!!'')')
      if(lms) then
        write(6,'('' Dos and magnetization in (lms)'',
     >'' representation on output'')')
      else
        write(6,'('' Dos and magnetization in (kappa,my)'',
     >'' representation on output'')')
      endif
c
c -bulk or surface
c
      read(5,*) bulk,dos
      bulk=.false.
c     if(bulk) write(6,'('' Bulk calculation!!!'')')
      if(dos) write(6,'('' Only DOS will be calculated !!!'')')
c
c -enter lmax 
c
      read(5,*) lmax
      write(6,'(/'' lmax='',i2)') lmax
      if(lmax.gt.lmaxp) then
         write(6,*) '<readini>: STOP: Increase lmaxp in param.h to',lmax
         stop 
      endif
c
c - exchange
c
      read(5,*) ixch
      write(6,'(/'' ixch='',i3)') ixch  
      if(ixch.le.9) call setxcp(ixch,2)
c
c -ewald parameter for 2D structure constants and Madelung constants
c
      read(5,*) eta,sigma!,madvers 01/03/2021 Laszloffy
      madvers=2
      write(6,'('' eta: '',f8.4,5x,''sigma: '',f8.4/)') eta,sigma
c     if(madvers.eq.1) write(6,'(''Calculate bulk Madelung '',
c    >''constants by means of 3D Ewald method'')')
c     if(madvers.eq.2) write(6,'(''Calculate bulk Madelung '',
c    >''constants by means of 2D Ewald method'')')
c     if(madvers.eq.3) write(6,'(''Calculate bulk Madelung '',
c    >''constants by means of 3D and 2D Ewald methods'')')
c
c -input for screening and alpha's
c
      read(5,*) iscreen,vscreen
      iscreen=0
      if(iscreen.le.0) write(6,'(/'' No screening at all !'')')
c     if(iscreen.ge.1) then
c       write(6,'(/'' Screening with square well potentials !'')')
c       write(6,'('' Height of potential: '',f5.2,'' ry'')') vscreen
c     end if
c
c  with or without orbital polarization
c
      read(5,*) orbpol,opotl,opotr,opot
      if(orbpol) write(6,'(/'' Orbital polarization included !!!'')')
c
c Laszloffy 06/04/17
c -rotate the induced moments or not
c
c     write(6,'('' Pre successed ! '')')
      read(5,*) isscf,rotind,spinorb,dmmin
c     write(6,'('' After successed ! '')')
      if(.not.isscf) write(6,'(/'' No SCF calculations !!!'')')
      if(rotind) write(6,'(/'' Direction of induces spins will be
     >     updated !'')')
      if(spinorb) write(6,'(/'' Total magnetic moment used for new dir!'
     >')')
c
c -vacuum or bulk on the right
c
      read(5,*) rightm,v0,newvac
c     if(bulk) rightm='B'
c     write(6,'('' from the right is: '',a)') rightm
c     if(rightm.eq.'V') then
c       write(6,'('' Constant potential on the vacuum side :'',
c    >            f10.5,'' ryd'')') v0
c       if(newvac.eq.0) then
c         write(6,'('' Vacuum potential level is kept constant'')')
c       else
c         write(6,'('' Vacuum potential level will be updated'')')
c       end if
c     end if 
c
      read(5,*) qc,cpc
      qc=0.d0
      cpc=0.d0
c     write(6,'(/'' Charge and position(a.u.) of external capacitor'',
c    >          2f8.4)') qc,cpc
c
c -starting with 'square well' potentials for vacuum layers
c
      read(5,*) ivacpot,v00
      read(5,*) layshift1,layshift2,potshift
c     if(layshift1.gt.layshift2) then
c       layshift = layshift2
c       layshift2=layshift1
c       layshift1=layshift
c     endif
c     if(ivacpot.eq.1) 
c    > write(6,'('' Starting potential for vacuum layers : '',
c    >           f10.5,'' ryd'')') v00
c     if(layshift1.ne.0)
c    > write(6,'('' Shift potential on layers'',2i3,'' by'',
c    > f6.2,'' ry'')') layshift1,layshift2,potshift
c
c -starting with auxilliary exchange splitting
c
      read(5,*) ib0,b0
      read(5,*) layb0,layb0f
c     if(layb0.gt.layb0f) then
c       layb0i = layb0f
c       layb0f=layb0
c       layb0=layb0i
c     endif
c     if(layb0.eq.0) then
c       if(ib0.eq.0) then
c         write(6,'(/'' Add magnetic field: '',f10.5,'' ryd '',
c    >    ''at all layers'')') b0
c       else
c         write(6,'(/'' Start with magnetic field: '',f10.5,'' ryd '',
c    >    ''at all layers'')') b0
c       end if
c     else
c       if(ib0.eq.0) then
c         write(6,'(/'' Add magnetic field: '',f10.5,'' ryd '',
c    >    ''at layers'',2i3)') b0,layb0,layb0f
c       else
c       write(6,'(/'' Start with magnetic field: '',f10.5,'' ryd '',
c    >    ''at layers'',2i3)') b0,layb0,layb0f
c       end if
c     end if
c
c -iteration parameters for surface Green's function
c
      read(5,*) itermaxl,itermaxr,errmax,ichk
c     write(6,'(/'' Calculation of surface greens function:''/
c    >          ''  itermaxl= '',i3,''  itermaxr= '',i3,
c    >          ''  errmax='',1pd10.2)') itermaxl,itermaxr,errmax
c
c - CPA parameters
c
      read(5,*) itcpam,cpatol,cpatest
c     read(5,*) itcpam,cpatol,cpatest,lliter
c
c     write(6,'(/'' Maximum number of CPA iterations: '',i3,
c    >          /'' CPA tolerance: '',d12.4)') itcpam,cpatol
      lliter=1
c     write(6,'('' Number of iterations for Lloyds formula: '',i3)')
c    > lliter
c
c -parameters for Broyden's method
c
      read(5,*) cmix,wbr,nbrcase
      if(nbrcase.gt.100) 
     >  write(6,*) 
     >  '<readini>: STOP: Effective dimension of nbrmax too large'
      if(nbrcase.gt.100) stop 'Effective dimension of nbrmax too large'
      read(5,*) (nbrmax(i),i=1,nbrcase)
c
      itscfmax=0
      do i=1,nbrcase
        itscfmax=itscfmax+nbrmax(i)
      end do
      if(dos) itscfmax=1
c
      if(.not.dos) then
      write(6,'(/'' Mixing factors:'',2f10.5)') cmix,wbr
      write(6,'(/'' Broyden moduls:'',20i3)') (nbrmax(i),i=1,nbrcase)
      write(6,'(/'' Total number of iterations:'',i5)') itscfmax
      end if
c
c -mixing factor for orientation of spin fields
c
      read(5,*) dmix,wbrdir,nbrdir
      if(.not.dos) then
      write(6,'(/'' Mixing parameters for scf-spin reorientation: '',
     >2f10.5,i5)') dmix,wbrdir,nbrdir
      end if
c
c do scf iterations until precision as follows reached
c
      read(5,*) tolerr,tolef,tolen
      write(6,'(/'' Relative accuracy for scf potentials'',d10.3)')
     >tolerr
      if(bulk) write(6,'('' Absolute accuracy for Fermi energy'',
     >d10.3)') tolef
      if(bulk) write(6,'('' Absolute accuracy for total energy'',
     >d10.3)') tolen
c
c**********************************************************************
c  E- and k-mesh
c**********************************************************************
c -enter parameters for complex energy contour
c
      read(5,*) imesh,npanel
c
      if(dos.and.(imesh.ne.0)) then
       imesh=0
       write(6,'(/''<readini>: DOS calculation: imesh set to zero!!'')')
      endif
c
      if(imesh.ge.2.and.npanel.ne.1) then
        write(6,*) '<readini>: STOP:  Number of panels should be 1'
        stop ' Number of panels should be 1'
      end if
c
      if(npanel.gt.5) then
         write(6,*) '<readini>: STOP: Max. npanel allowed is 5'
         stop 'Max. npanel allowed is 5'
      end if
c
      ne=0
      do n=1,npanel
       read(5,*) nepanel(n),ne1,ne2,ne3,ebottom(n),etop(n),eps(n)
       write(6,'(/'' Panel:'',i2)') n
       if(imesh.eq.0)  then
         write(6,'('' Paralell to real axis'')')
         write(6,'('' energy range:'',t35,2f15.10)') ebottom(n),etop(n)
         write(6,'('' imaginary part:'',t35,f12.8)') eps(n)
         write(6,'('' No. of energy points:'',t35,i3)') nepanel(n)
       else if(imesh.eq.1)  then
         write(6,'('' Semicircle contour'')')
         write(6,'('' energy range:'',t35,2f15.10)') ebottom(n),etop(n)
         write(6,'('' epsilon     :'',t35,f15.10)') eps(n)
         write(6,'('' No. of energy points:'',t35,i2)') nepanel(n)
       else if(imesh.ge.2)  then
         write(6,'('' Matsubara poles'')')
         write(6,'('' Temperature (K):'',t35,f15.10)') eps(1)
         write(6,'('' Fermi level:'',t35,2f15.10)') etop(1)
         write(6,'('' No. of poles:'',t35,i2)') nepanel(1)
         if(imesh.eq.3)  then
         write(6,'('' Points on the quarter-circle:'',t35,i3)') ne1
         write(6,'('' Points parallel to real axis :'',t35,i3)') ne2
         write(6,'('' Points around the Fermi level :'',t35,i3)') ne3
         end if
       else
         write(6,*) '<readini>: STOP: Unknown imesh'
         stop 'Unknown imesh'
       end if
       if(imesh.eq.3) then
               ne=ne1+ne2+ne3+nepanel(1)
       else   
               ne=ne+nepanel(n)
       end if
      end do
c
      if(ne.gt.me) then
         write(6,*) '<readini>: STOP: Increase me in param.h to',me
         stop 
      endif
c
c -enter k-parallel
c  Follow Cunnigham convention for reference axes
c
      read(5,*) intbz,kunning,park
      read(5,*) ntab
      if(ntab.gt.20) then
        write(6,*) '<readini>: STOP: Max. ntab allowed is 20'
        stop 'Max. ntab allowed is 20'
      end if
c
      do itab=1,ntab
        read(5,*) ksettab(itab),intab(itab)
      end do
c
c     if(intbz.eq.0) then
c       write(6,'(/'' IBZ integration'')')
c     elseif(intbz.eq.1) then
c       write(6,'(/'' Point group operations applied'')')
c     else
c       write(6,'(/'' Brute force full BZ integration'')')
c     end if
c
      ie=0
      do itab=1,ntab
        do k=1,intab(itab)
           ie=ie+1
        end do
      end do
      if(ie.ne.ne) then 
         write(6,*) '<readini>: STOP: check table for nr. of e-points'
         write(6,*) '<readini>: STOP: actual number ie=',ie
         write(6,*) '<readini>: STOP: set number ne=',ne
         stop 'check table for nr. of e-points'
      end if
c
      call flush(6)
c
      return
      end
