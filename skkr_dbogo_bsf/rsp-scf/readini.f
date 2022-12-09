c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine readini(
     >imesh,ne,npanel,nepanel,ne1,ne2,ne3,ebottom,etop,eps,lmax,
     >eta,sigma,park,kset,ksetmax,itcpam,cpatol,cpatest,cpain,
     >printout,potout,dosout,pdosout,dosoutimp,pdosoutimp,tripletout,
     >tmatout,momout,
     >leftpot,rightpot,laypot,laycore,
     >leftmat,rightmat,laymat,leftmom,rightmom,
     >wrel,lms,bulk,dos,rightm,newvac,ivacpot,v0,v00,vrsh,
     >vscreen,iscreen,itscfmax,tolerr,tolef,tolen,intbz,kunning,
     >orbpol,opotl,opotr,opot,lliter,
c     >E_Fermi,leftdelta,rightdelta,laydelta,singlesite,reg,c_light,
     >E_Fermi,singlesite,reg,c_light,
     >clucalc,hostprint,nimp,npair,npair0,npair1)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
      logical bulk,dos,wrel,lms,cpatest,cpain,singlesite,reg
      logical orbpol,opotl,opotr,opot
      character*25 job
      character*30 leftpot,rightpot,laypot,laycore
c      character*30 leftdelta,rightdelta,laydelta
      character*30 leftmat,rightmat,laymat
      character*30 leftmom,rightmom           
      character*30 printout,potout,dosout,pdosout,tmatout,momout
      character*30 dosoutimp
      character*34 pdosoutimp
      character*34 tripletout
      character*1 rightm
      dimension nepanel(5),ebottom(5),etop(5),eps(5)
      dimension ksettab(20),intab(20),kset(me)
      dimension park(2),r(3,minprc)
      dimension isigba(mintfc),isigbb(mintfc)
c
      common/test/itest
      common/iterpar/errmax,itermaxl,itermaxr,ichk
      common/broypot/cmix,wbr,nbrmax(100)
      common/broydir/dmix,wbrdir,nbrdir
      common/xch/ixch
      common/extf/qc,cpc
      common/pot/potshift,b0,layshift1,layshift2,layb0,layb0f,
     &           ib0,isigba,isigbb,nomag
      common/madbulk/madvers
c MPI
      integer root, myrank, nprocs
      common/mpi/root,myrank,nprocs
c Version setting
      integer IOstatus
      logical clucalc,hostprint
      integer nimp,npair,npair1
c
c***************************************************************
c INPUT/OUTPUT files
c***************************************************************
c -initialize jobname
c
      read(5,'(a25)') job
c
c ******* open files according to job 
c
c     --------------------------------------------------------------
      call jobnam(job,printout,potout,dosout,pdosout,dosoutimp,tmatout,
     >            momout,pdosoutimp,tripletout)
c     --------------------------------------------------------------
      if (myrank == root) then
      open(6,file=printout,status='unknown')
      else
      open(6,file='/dev/null',status='unknown')
      end if
c
      write(6,*)
c
      write(6,'(2x,''routine READINI>''/)')
      write(6,'('' control and test output      : '',
     >          a)') printout
      write(6,'('' new potentials                : '',
     >          a)') potout
      write(6,'('' densities of states           : '',
     >          a)') dosout
      write(6,'('' partial densities of states   : '',
     >          a)') pdosout
      write(6,'('' effective t-matrices          : '',
     >          a)') tmatout
      write(6,'('' moments of charge density     : '',
     >          a)') momout
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
      write(6,'('' moments of charge density for lhs bulk  : '',
     >          a)') leftmom
      write(6,'('' moments of charge density for rhs bulk  : '',
     >          a)') rightmom
c
c -test parameter
      read(5,*) itest
c
c***************************************************************
c Conditions of calculation
c***************************************************************
!
! local coordinate systems or global?  --> Localmode is still not working
      ! localmode is defined in param.h!
      if(localmode) then
        write(6,'(/" Local coordinate systems in use. NOT IMPLENTED!")')
        stop 'Change localmode false in param.h'
      else
        write(6,'(/" Global coordinate systems in use.")')
      end if

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
      if(bulk) write(6,'('' Bulk calculation!!!'')')
      if(dos) write(6,'('' Only DOS will be calculated !!!'')')
c
c -enter lmax 
c
      read(5,*) lmax
      write(6,'(/'' lmax='',i2)') lmax
      if(lmax.gt.lmaxp) then
         write(6,*) 'Increase lmaxp in param.h to',lmax
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
      read(5,*) eta,sigma!,madvers 03/03/2021 Laszloffy
      madvers=2
      write(6,'('' eta: '',f8.4,5x,''sigma: '',f8.4/)') eta,sigma
      if(madvers.eq.1) write(6,'('' Calculate bulk Madelung '',
     >''constants by means of 3D Ewald method'')')
      if(madvers.eq.2) write(6,'('' Calculate bulk Madelung '',
     >''constants by means of 2D Ewald method'')')
      if(madvers.eq.3) write(6,'('' Calculate bulk Madelung '',
     >''constants by means of 3D and 2D Ewald methods'')')
c
c -input for screening and alpha's
c
      read(5,*) iscreen,vscreen
      if(iscreen.le.0) write(6,'(/'' No screening at all !'')')
      if(iscreen.eq.1) then
        write(6,'(/'' Screening with square well potentials !'')')
        write(6,'('' Height of potential: '',f5.2,'' ry'')') vscreen
      end if
      if(iscreen.gt.1)
     >  write(6,'(/'' Screening with hard-sphere potentials !'')')
c
c  with or without orbital polarization
c
      read(5,*) orbpol,opotl,opotr,opot
      if(orbpol) write(6,'(/'' Orbital polarization included !!!'')')
c
c -vacuum or bulk on the right
c
      read(5,*) rightm,v0,newvac
      if(bulk) rightm='B'
      write(6,'('' from the right is: '',a)') rightm
      if(rightm.eq.'V') then
        write(6,'('' Constant potential on the vacuum side :'',
     >            f10.5,'' ryd'')') v0
        if(newvac.eq.0) then
          write(6,'('' Vacuum potential level is kept constant'')')
        else
          write(6,'('' Vacuum potential level will be updated'')')
        end if
      end if 
c
      read(5,*) qc,cpc
      write(6,'(/'' Charge and position(a.u.) of external capacitor'',
     >          2f8.4)') qc,cpc
c
c -starting with 'square well' potentials for vacuum layers
c
      read(5,*) ivacpot,v00,vrsh
      if(ivacpot.eq.1) 
     > write(6,'(/'' Starting potential for vacuum layers : '',
     >           f10.5,'' ryd'')') v00
      if(bulk) write(6,'(/'' Average potential at WS radii:'',
     >f8.3,'' Ryd'')') vrsh
      call flush(6)
c
      read(5,*) layshift1,layshift2,potshift
      if(layshift1.gt.layshift2) then
        layshift = layshift2
        layshift2=layshift1
        layshift1=layshift
      endif
      if(layshift1.ne.0)
     > write(6,'(/'' Shift potential on layers'',2i3,'' by'',
     > f6.2,'' ry'')') layshift1,layshift2,potshift
      call flush(6)
c
c -starting with auxilliary exchange splitting
c
      read(5,*) ib0,b0,nomag
      if(nomag.eq.1) 
     >write(6,'(/'' No exchange field in the input potential file!'')')
      read(5,*) layb0,layb0f
      if(layb0.gt.layb0f) then
        layb0i = layb0f
        layb0f=layb0
        layb0=layb0i
      endif
      if(layb0.eq.0) then
        if(ib0.eq.0) then
          write(6,'(/'' Add magnetic field: '',f10.5,'' ryd '',
     >    ''at all layers'')') b0
        else
          write(6,'(/'' Start with magnetic field: '',f10.5,'' ryd '',
     >    ''at all layers'')') b0
        end if
      else
        if(ib0.eq.0) then
          write(6,'(/'' Add magnetic field: '',f10.5,'' ryd '',
     >    ''at layers'',2i3)') b0,layb0,layb0f
        else
        write(6,'(/'' Start with magnetic field: '',f10.5,'' ryd '',
     >    ''at layers'',2i3)') b0,layb0,layb0f
        end if
      end if
      call flush(6)
c
c -iteration parameters for surface Green's function
c
      read(5,*) itermaxl,itermaxr,errmax,ichk
      write(6,'(/'' Calculation of surface greens function:''/
     >          ''  itermaxl= '',i3,''  itermaxr= '',i3,
     >          ''  errmax='',1pd10.2)') itermaxl,itermaxr,errmax
      call flush(6)
c
c - CPA parameters
c
      read(5,*) itcpam,cpatol,cpatest,cpain
c     read(5,*) itcpam,cpatol,cpatest,lliter
c
      write(6,'(/'' Maximum number of CPA iterations: '',i3,
     >          /'' CPA tolerance: '',d12.4)') itcpam,cpatol
      lliter=1
c     write(6,'('' Number of iterations for Lloyds formula: '',i3)')
c    > lliter
      if(cpain) write(6,'(/'' Force to read in effective t-matrices'',
     >'' in case of CPA & DOS calculation'')')
      call flush(6)
c
c -parameters for Broyden's method
c
      read(5,*) cmix,wbr,nbrcase
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
      call flush(6)
c
c -mixing factor for orientation of spin fields
c
      read(5,*) dmix,wbrdir,nbrdir
      if(.not.dos) then
      write(6,'(/'' Mixing parameters for scf-spin reorientation: '',
     >2f10.5,i5)') dmix,wbrdir,nbrdir
      end if
      call flush(6)
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
      call flush(6)
c
c**********************************************************************
c  E- and k-mesh
c**********************************************************************
c -enter parameters for complex energy contour
c
      read(5,*) imesh,npanel
c     if(dos.and.(imesh.ne.0)) then
c       imesh=0
c       write(6,'(/'' DOS calculation: imesh set to zero!!'')')
c     endif
      if(imesh.ge.2.and.npanel.ne.1) stop 
     >' Number of panels should be 1'
      if(npanel.gt.5) stop 'Max. npanel allowed is 5'
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
         stop 'Unknown imesh'
       end if
       if(imesh.eq.3) then
               ne=ne1+ne2+ne3+nepanel(1)
       else   
               ne=ne+nepanel(n)
       end if
      end do
      if(ne.gt.me) then
         write(6,*) 'Increase me in param.h to',me
         stop 
      endif
      call flush(6)
c
c -enter k-parallel
c  Follow Cunnigham convention for reference axes
c
      read(5,*) intbz,kunning,park
      read(5,*) ntab
      if(ntab.gt.20) stop 'Max. ntab allowed is 20'
      do itab=1,ntab
        read(5,*) ksettab(itab),intab(itab)
      end do
c
      if(intbz.eq.0) then
        write(6,'(/'' IBZ integration'')')
      elseif(intbz.eq.1) then
        write(6,'(/'' Point group operations applied'')')
      else
        write(6,'(/'' Brute force full BZ integration'')')
      end if
      call flush(6)
c
      ie=0
      ksetmax=0
      do itab=1,ntab
        do k=1,intab(itab)
           ie=ie+1
           kset(ie)=ksettab(itab)
           if(kset(ie).gt.ksetmax) ksetmax=kset(ie)
        end do
      end do
      if(ie.ne.ne) then
         write(6,*) ie,ne,'check table for nr. of k-points'
         stop 'check table for nr. of k-points'
      end if 
c
      call flush(6)
c
      read(5,*) E_Fermi, c_light
      write(6,'(/'' Fermi Energy ='',f8.4)') E_Fermi
      write(6,'(/'' Speed of light ='',f8.4)') c_light
      read(5,*) singlesite,reg
      if(singlesite) write(6,'(/'' Single-site calculation!
     > Tau is replaced by the t-matrix!'')')
      if(reg) write(6,'(/''OnlyReg is .true. --> 
     > Only the regular wavefunctions are calculated!'')')
      if(.not.reg) 
     > write(6,'(/'' Irregular wavefunctions are calculated!'')')
c Version setting parameters
      read(5,*,IOSTAT=IOstatus) clucalc,hostprint
      if(IOstatus.ne.0) then
        clucalc=.false.
        hostprint=.false.
      end if
      if(clucalc.or.hostprint) then
        close(5)
        open(unit=5,file='cluster_geo.in',status='old')
        read(5,*) nimp
        npair0=nimp*(nimp-1)
        if(npair0.eq.0) npair0=1
      else
        nimp=0
        npair=0
        npair0=1
        npair1=1
      end if
      if(clucalc) then
        write(6,'(/'' Embedded cluster calculations
     > <-- clucalc is .true. For  '', i5, '' impurities '')') nimp
        if(hostprint) then
          write(6,'(/''   host t and tau matrices will be read from
     > files <-- hostprint is .true.'')')
        else
          write(6,'(/''   host t and tau matrices will be calculated
     > on the fly <-- hostprint is .false.'')')
        end if
      else
        write(6,'(/'' Bulk/layer calculations
     > <-- clucalc is .false.'')')
        if(hostprint) then
          write(6,'(/''   host t and tau matrices will be written to
     > files <-- hostprint is .true.'',i5,'' impurities '')')nimp
        end if
      end if
c
      call flush(6)
c
      return
      end