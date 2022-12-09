c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine readini(
     >imesh,ne,npanel,nepanel,ne1,ne2,ne3,ebottom,etop,
     >eps,lmax,eta,sigma,park,kset,itcpam,cpatol,cpatest,cpain,
     >printout,potout,dosout,pdosout,tmatout,momout,
     >leftpot,rightpot,laypot,laycore,
     >leftmat,rightmat,laymat,leftmom,rightmom,
     >bulk,nspheres,rel,dos,rightm,newvac,ivacpot,
     >v0,v00,vrsh,iscreen,vscreen,itscfmax,tolerr,tolef,intbz,kunning,
     >dlm,dlml,dlmr)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
c MPI parameters (esimon)
      integer root, myrank, nprocs
      common/mpi/root,myrank,nprocs
c      
      logical rel,bulk,nspheres,dos,cpatest,dlm,dlml,dlmr,cpain
      character*25 job
      character*30 leftpot,rightpot,laypot,laycore
      character*30 leftmat,rightmat,laymat
      character*30 leftmom,rightmom
      character*30 printout,potout,dosout,pdosout,tmatout,momout
      character*1 rightm
      dimension nepanel(5),ebottom(5),etop(5),eps(5)
      dimension ksettab(20),intab(20),kset(me)
      dimension park(2)
c
      common/test/itest
      common/iterpar/errmax,itermaxl,itermaxr,ichk
      common/broy/cmix,wbr,nbrmax(100),nbrcase
      common/xch/ixch,ixchout
      common/extf/qc,cpc
      common/pot/potshift,b0,layshift1,layshift2,layb0,layb0f,ib0
      common/madbulk/madvers
c
c***************************************************************
c INPUT/OUTPUT files
c*************************************************************** 
c
c -initialize jobname
c
      read(5,'(a25)') job
c
c ******* open files according to job
c
c     ------------------------------------------------------------------
      call jobnam(job,printout,potout,dosout,pdosout,tmatout,momout)
c     ------------------------------------------------------------------
      if (myrank == root) then
      open(6,file=printout,status='unknown')
      else
      open(6,file='/dev/null',status='unknown')
      end if
c
      write(6,*)
      write(6,'(2x,''routine READINI>''/)')
      write(6,'('' controll and test output      : '',
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
c
      read(5,*) itest
c
c***************************************************************
c Conditions of calculation
c***************************************************************
c
c -non/scalar relativistic
c
      read(5,*) rel
      if(rel) then
        write(6,'(/'' Scalar relativistic calculation'')')
       else
        write(6,'(/'' Nonrelativistic calculation'')')
      end if
c
c -bulk or surface
c
      read(5,*) bulk,dos,nspheres
      if(bulk) then
         write(6,'(/'' Bulk calculation!!!'')')
         if(nspheres) then
           write(6,'('' Force neutral spheres!!!'')')
         end if
      else
         nspheres=.false.
      end if
      if(dos) write(6,'(/'' Only DOS will be calculated !!!'')')
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
c -exchange
c
      read(5,*) ixch,ixchout
      write(6,'(/'' ixch='',i3)') ixch  
      if(ixch.le.9) call setxcp(ixch,2)
      if(ixch.eq.10) 
     >write(6,'('' Exchange-correlation type: PW92'')') 
      if(ixch.eq.11) 
     >write(6,'('' Exchange-correlation type: GGA PW91'')') 
      if(ixch.eq.12) 
     >write(6,'('' Exchange-correlation type: GGA PBE '')') 
      if(ixch.ge.13) 
     >write(6,'('' Exchange-correlation type: GGA PBEsol '')') 
      write(6,'('' ixchout='',i3/)') ixchout
c
c -ewald parameter for 2D structure constants and Madelung constants
c
      read(5,*) eta,sigma,madvers
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
c
c -vacuum or bulk on the right
c
      read(5,*) rightm,v0,newvac
      if(bulk) rightm='B'
      write(6,'('' from the right is: '',a)') rightm
      if(rightm.eq.'V') then
        write(6,'('' Constant potential on the vacuum side :'',
     >            f10.5,'' ryd'')') v0
        if(.not.dos) then
        if(newvac.eq.0) then
          write(6,'('' Vacuum potential level is kept constant'')')
        else
          write(6,'('' Vacuum potential level will be updated'')')
        end if
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
     >write(6,'(/'' Starting potential for vacuum layers : '',
     >          f10.5,'' ryd'')') v00
      if(bulk) write(6,'(/'' Average potential at WS radii:'',
     >f8.3,'' Ryd'')') vrsh
c
      read(5,*) layshift1,layshift2,potshift
      if(layshift1.gt.layshift2) then
        layshift = layshift2
        layshift2=layshift1
        layshift1=layshift
      endif
      if(layshift1.ne.0) 
     >  write(6,'(/'' Shift potential on layers'',2i3,'' by'',
     >  f6.2,'' ry'')') layshift1,layshift2,potshift
c
c -starting with auxilliary exchange splitting
c
      read(5,*) ib0,b0
      read(5,*) layb0,layb0f
      if(layb0.gt.layb0f) then
        layb0i = layb0f
        layb0f=layb0
        layb0=layb0i
      endif
      if(layb0.eq.0) then
        if(ib0.eq.0) then
          write(6,'(/'' Add magnetic field: '',f10.5,'' ryd
     >     at all layers'')') b0
        else
          write(6,'(/'' Start with magnetic field: '',f10.5,'' ryd
     >     at all layers'')') b0
        end if
      else
        if(ib0.eq.0) then
          write(6,'(/'' Add magnetic field: '',f10.5,'' ryd
     >     at layers'',2i3)') b0,layb0,layb0f
        else
        write(6,'(/'' Start with magnetic field: '',f10.5,'' ryd
     >     at layers'',2i3)') b0,layb0,layb0f
        end if
      end if 
c
c - option for DLM
c
      read(5,*) dlm,dlml,dlmr
      if(dlm) write(6,'(/'' Disordered Local Moment calculation'')')
      if(.not.bulk) then
        if(dlml) write(6,'('' Left substrate DLM '')')
        if(rightm.eq.'B'.and.dlmr) 
     >  write(6,'('' Right substrate DLM '')')
      end if
c
c -iteration parameters for surface Green's function
c
      read(5,*) itermaxl,itermaxr,errmax,ichk
      write(6,'(/'' Calculation of surface greens function:''/
     >          ''   itermaxl= '',i3,''  itermaxr= '',i3,
     >          ''  errmax='',1pd10.2)') itermaxl,itermaxr,errmax
c
c - CPA parameters
c
      read(5,*) itcpam,cpatol,cpatest,cpain
c
      write(6,'(/'' Maximum number of CPA iterations: '',i3,
     >          /'' CPA tolerance: '',d12.4)') itcpam,cpatol
      if(cpain) write(6,'(/'' Force to read in effective t-matrices'',
     >'' in case of CPA & DOS calculation'')')
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
      write(6,'(/'' Mixing factors:'',3f10.5)') cmix,wbr
      write(6,'('' Broyden moduls:'',20i3)') (nbrmax(i),i=1,nbrcase)
      write(6,'('' Total number of iterations:'',i5)') itscfmax
      end if
c
c do scf iterations until precision as follows reached
c
      read(5,*) tolerr,tolef
      write(6,'(/'' Relative accuracy for scf potentials'',d10.3)')
     >tolerr
      if(bulk) write(6,'('' Absolute accuracy for Fermi energy'',
     >d10.3)') tolef
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
c
c     if(.not.dos.and.ne.gt.me1) then
c       write(6,'('' ne='',i2,'' greater then me1='',i2)') ne,me1
c       stop
c     end if
c
c -enter k-parallel
c  Follow Cunnigham convention for reference axe.
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
      if(ie.ne.ne) stop 'check table for nr. of k-points' 
c
      return
      end
