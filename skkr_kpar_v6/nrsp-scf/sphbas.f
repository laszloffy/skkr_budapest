c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine sphbas(dg,dg1,ng,invg)
c======================
c
c Calculates the matrix representations of the groups:
c      C1,C2,C3,C4,C6,Cs,C2v,C3vA,C3vB,C4v,C6v
c Source:  S. Altman & P. Herzig: Point-group theory tables
c
c List of parameters:
c Input : idgroup  : group identifier. 
c
c Output: dg(lmmaxp,lmmaxp,melem) : the matrix reprs. of the group elem.
c         dg1(lmmaxp,lmmaxp,melem): the inverse matrices
c         ng : number of group elements
c         invg : labels the inverse transformations in the group
c Intermediate : o : transforms spherical harmonics in Condon-Shortly
c                    convention to the symmetry adopted bases.
c
c Subroutine matrep: generates the irreducible matrix representations 
c
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      parameter (ndim=25)
c
      character*4 idgroup
c
      dimension invg(melem),idimrep(6),irep(ndim)
c
      complex*16 dg(lmmaxp,lmmaxp,melem)
      complex*16 dg1(lmmaxp,lmmaxp,melem)
      complex*16 o(ndim,ndim),u(ndim,ndim)
      complex*16 d(ndim,ndim)
      complex*16 rep(2,2,6,melem)
      complex*16 one,sqrm1,facn,faci
c
      common/latt2d/idgroup,ilat,a2d,b2d,d2d
      data tol/1.0d-14/
c
c     write(6,'(/,'' idgroup : '',a4,/)')idgroup
      if(ndim.lt.lmmaxp) then
        write(6,*) 'SPHBAS: Increase parameter ndim to',lmmaxp
        stop
      endif
c
c First set up matrix for basis transformation and get information
c on the sequence of irreducible sub-spaces
c
c  numrep: number of irred. sub-spaces up to a matrix dimension of ndim 
c  irep(ir): index of irred. representation corresponding to sub-space
c            labelled by ir (1 <= ir <= numrep)
c
      one  =dcmplx(1.0d0,0.0d0)
      sqrm1=dcmplx(0.0d0,1.0d0)
      facn =dcmplx(1.d0/dsqrt(2.0d0),0.0d0)
      faci =dcmplx(0.0d0,1.d0/dsqrt(2.0d0))
      call czero(o,ndim*ndim)
c
      if(idgroup.eq.'C1') then
c
        do i=1,ndim
          o(i,i)=1.0d0
        end do
        numrep=ndim
        do ir=1,numrep
          irep(ir)=1
        end do
c
      else if(idgroup.eq.'C2') then
c
        do i=1,ndim
          o(i,i)=1.0d0
        end do
        numrep=ndim
        irep(1)= 1
        irep(2)= 2
        irep(3)= 1
        irep(4)= 2
        irep(5)= 1
        irep(6)= 2
        irep(7)= 1
        irep(8)= 2
        irep(9)= 1
        irep(10)=2
        irep(11)=1
        irep(12)=2
        irep(13)=1
        irep(14)=2
        irep(15)=1
        irep(16)=2
        irep(17)=1
        irep(18)=2
        irep(19)=1
        irep(20)=2
        irep(21)=1
        irep(22)=2
        irep(23)=1
        irep(24)=2
        irep(25)=1
c
      else if(idgroup.eq.'C3') then
c
        do i=1,ndim
          o(i,i)=1.0d0
        end do
        numrep=ndim
        irep(1)= 1
        irep(2)= 3
        irep(3)= 1
        irep(4)= 2
        irep(5)= 2
        irep(6)= 3
        irep(7)= 1
        irep(8)= 2
        irep(9)= 3
        irep(10)=1
        irep(11)=2
        irep(12)=3
        irep(13)=1
        irep(14)=2
        irep(15)=3
        irep(16)=1
        irep(17)=3
        irep(18)=1
        irep(19)=2
        irep(20)=3
        irep(21)=1
        irep(22)=2
        irep(23)=3
        irep(24)=1
        irep(25)=2
c
      else if(idgroup.eq.'C4') then
c
        do i=1,ndim
          o(i,i)=1.0d0
        end do
        numrep=ndim
        irep(1)= 1
        irep(2)= 4
        irep(3)= 1
        irep(4)= 3
        irep(5)= 2
        irep(6)= 4
        irep(7)= 1
        irep(8)= 3
        irep(9)= 2
        irep(10)=3
        irep(11)=2
        irep(12)=4
        irep(13)=1
        irep(14)=3
        irep(15)=2
        irep(16)=4
        irep(17)=1
        irep(18)=3
        irep(19)=2
        irep(20)=4
        irep(21)=1
        irep(22)=3
        irep(23)=2
        irep(24)=4
        irep(25)=1
c
      else if(idgroup.eq.'C6') then
c
        do i=1,ndim
          o(i,i)=1.0d0
        end do
        numrep=ndim
        irep(1)= 1
        irep(2)= 4
        irep(3)= 1
        irep(4)= 3
        irep(5)= 5
        irep(6)= 4
        irep(7)= 1
        irep(8)= 3
        irep(9)= 6
        irep(10)=2
        irep(11)=5
        irep(12)=4
        irep(13)=1
        irep(14)=3
        irep(15)=6
        irep(16)=2
        irep(17)=6
        irep(18)=2
        irep(19)=5
        irep(20)=4
        irep(21)=1
        irep(22)=3
        irep(23)=6
        irep(24)=2
        irep(25)=5
c
      else if(idgroup.eq.'Cs') then
c
        irep(1)= 1
          o(1,1)=one
        irep(2)= 1
          o(2,2)=-facn
          o(4,2)=facn
        irep(3)= 2
          o(2,3)=facn
          o(4,3)=facn
        irep(4)= 1
          o(3,4)=one 
        irep(5)= 1
          o(5,5)=facn
          o(9,5)=facn
        irep(6)= 2
          o(5,6)=-facn
          o(9,6)=facn
        irep(7)= 1
          o(6,7)=-facn
          o(8,7)=facn
        irep(8)= 2
          o(6,8)=facn
          o(8,8)=facn
        irep(9)= 1
          o(7,9)=one 
        irep(10)=1
          o(10,10)=-facn
          o(16,10)=facn
        irep(11)=2
          o(10,11)=facn
          o(16,11)=facn
        irep(12)=1
          o(11,12)=facn
          o(15,12)=facn
        irep(13)=2
          o(11,13)=-facn
          o(15,13)=facn
        irep(14)=1
          o(12,14)=-facn
          o(14,14)=facn
        irep(15)=2
          o(12,15)=facn
          o(14,15)=facn
        irep(16)=1
          o(13,16)=one 
        irep(17)=1
          o(17,17)=facn
          o(25,17)=facn
        irep(18)=2
          o(17,18)=-facn
          o(25,18)=facn
        irep(19)=1
          o(18,19)=-facn
          o(24,19)=facn
        irep(20)=2
          o(18,20)=facn
          o(24,20)=facn
        irep(21)=1
          o(19,21)=facn
          o(23,21)=facn
        irep(22)=2
          o(19,22)=-facn
          o(23,22)=facn
        irep(23)=1
          o(20,23)=-facn
          o(22,23)=facn
        irep(24)=2
          o(20,24)=facn
          o(22,24)=facn
        irep(25)=1
          o(21,25)=one
c
        numrep=25
c
      else if(idgroup.eq.'C2v') then
c
        irep(1)=1
          o(1,1)=one
        irep(2)=3
          o(2,2)=-facn
          o(4,2)=facn 
        irep(3)=4
          o(2,3)=facn
          o(4,3)=facn
        irep(4)=1
          o(3,4)=one
        irep(5)=1
          o(5,5)=facn
          o(9,5)=facn
        irep(6)=2
          o(5,6)=-facn
          o(9,6)=facn
        irep(7)=3
          o(6,7)=-facn
          o(8,7)=facn
        irep(8)=4
          o(6,8)=facn 
          o(8,8)=facn
        irep(9)=1
          o(7,9)=one
        irep(10)=3
          o(10,10)=-facn
          o(16,10)=facn
        irep(11)=4
          o(10,11)=facn
          o(16,11)=facn
        irep(12)=1
          o(11,12)=facn
          o(15,12)=facn
        irep(13)=2
          o(11,13)=-facn
          o(15,13)=facn
        irep(14)=3
          o(12,14)=-facn
          o(14,14)=facn
        irep(15)=4
          o(12,15)=facn
          o(14,15)=facn
        irep(16)=1
          o(13,16)=one
        irep(17)=1
          o(17,17)=facn
          o(25,17)=facn
        irep(18)=2
          o(17,18)=-facn
          o(25,18)=facn
        irep(19)=3
          o(18,19)=-facn
          o(24,19)=facn
        irep(20)=4
          o(18,20)=facn
          o(24,20)=facn
        irep(21)=1
          o(19,21)=facn
          o(23,21)=facn
        irep(22)=2
          o(19,22)=-facn
          o(23,22)=facn
        irep(23)=3
          o(20,23)=-facn
          o(22,23)=facn
        irep(24)=4
          o(20,24)=facn
          o(22,24)=facn
        irep(25)=1
          o(21,25)=one
c
        numrep=25
c
      else if(idgroup.eq.'C3vA') then
c
        irep(1)=1
          o(1,1)=one
        irep(2)=3
          o(4,2)=one
          o(2,3)=one
        irep(3)=1
          o(3,4)=one
        irep(4)=3
          o(5,5)=one 
          o(9,6)=-one 
        irep(5)=3
          o(8,7)=one
          o(6,8)=one  
        irep(6)=1
          o(7,9)=one
        irep(7)=1
          o(10,10)=-facn
          o(16,10)=facn
        irep(8)=2
          o(10,11)=facn
          o(16,11)=facn
        irep(9)=3
          o(11,12)=one
          o(15,13)=-one
        irep(10)=3
          o(14,14)=one
          o(12,15)=one 
        irep(11)=1
          o(13,16)=one
        irep(12)=3
          o(25,17)=one
          o(17,18)=-one
        irep(13)=1
          o(18,19)=-facn
          o(24,19)=facn
        irep(14)=2
          o(18,20)=facn
          o(24,20)=facn
        irep(15)=3
          o(19,21)=one
          o(23,22)=-one
        irep(16)=3
          o(22,23)=one
          o(20,24)=one
        irep(17)=1
          o(21,25)=one
c
        numrep=17
c
      else if(idgroup.eq.'C3vB') then
c
        irep(1)=1
          o(1,1)=one
        irep(2)=3
          o(4,2)=one
          o(2,3)=-one
        irep(3)=1
          o(3,4)=one
        irep(4)=3
          o(5,5)=one 
          o(9,6)=-one 
        irep(5)=3
          o(8,7)=one
          o(6,8)=-one  
        irep(6)=1
          o(7,9)=one
        irep(7)=1
          o(10,10)=facn
          o(16,10)=facn
        irep(8)=2
          o(10,11)=-facn
          o(16,11)=facn
        irep(9)=3
          o(11,12)=one
          o(15,13)=-one
        irep(10)=3
          o(14,14)=one
          o(12,15)=-one 
        irep(11)=1
          o(13,16)=one
        irep(12)=3
          o(25,17)=one
          o(17,18)=-one
        irep(13)=1
          o(18,19)=facn
          o(24,19)=facn
        irep(14)=2
          o(18,20)=-facn
          o(24,20)=facn
        irep(15)=3
          o(19,21)=one
          o(23,22)=-one
        irep(16)=3
          o(22,23)=one
          o(20,24)=-one
        irep(17)=1
          o(21,25)=one
c
        numrep=17
c
      else if(idgroup.eq.'C4v') then
c 
        irep(1)=1
          o(1,1)=one
        irep(2)=1
          o(3,2)=one
        irep(3)=5
          o(4,3)=one
          o(2,4)=-one
        irep(4)=3
          o(5,5)=facn
          o(9,5)=facn 
        irep(5)=4
          o(5,6)=-facn
          o(9,6)=facn  
        irep(6)=5
          o(8,7)=one
          o(6,8)=-one
        irep(7)=1
          o(7,9)=one
        irep(8)=5
          o(10,10)=one
          o(16,11)=-one 
        irep(9)=3
          o(11,12)=facn
          o(15,12)=facn
        irep(10)=4
          o(11,13)=-facn
          o(15,13)=facn 
        irep(11)=5
          o(14,14)=one
          o(12,15)=-one
        irep(12)=1
          o(13,16)=one
        irep(13)=1
          o(17,17)=facn
          o(25,17)=facn
        irep(14)=2
          o(17,18)=-facn
          o(25,18)=facn
        irep(15)=5
          o(18,19)=one
          o(24,20)=-one
        irep(16)=3
          o(19,21)=facn
          o(23,21)=facn
        irep(17)=4
          o(19,22)=-facn
          o(23,22)=facn
        irep(18)=5
          o(22,23)=one
          o(20,24)=-one
        irep(19)=1
          o(21,25)=one
c
        numrep=19
c
      else if(idgroup.eq.'C6v') then
c
        irep(1)=1
          o(1,1)=one
        irep(2)=5
          o(4,2)=one
          o(2,3)=-one
        irep(3)=1
          o(3,4)=one
        irep(4)=6
          o(5,5)=one 
          o(9,6)=-one 
        irep(5)=5
          o(8,7)=one
          o(6,8)=-one  
        irep(6)=1
          o(7,9)=one
        irep(7)=3
          o(10,10)=-facn
          o(16,10)=facn
        irep(8)=4
          o(10,11)=facn
          o(16,11)=facn
        irep(9)=6
          o(11,12)=one
          o(15,13)=-one
        irep(10)=5
          o(14,14)=one
          o(12,15)=-one 
        irep(11)=1
          o(13,16)=one
        irep(12)=6
          o(25,17)=one
          o(17,18)=-one
        irep(13)=3
          o(18,19)=-facn
          o(24,19)=facn
        irep(14)=4
          o(18,20)=facn
          o(24,20)=facn
        irep(15)=6
          o(23,21)=one
          o(19,22)=-one
        irep(16)=5
          o(22,23)=one
          o(20,24)=-one
        irep(17)=1
          o(21,25)=one
c
        numrep=17
c
      else
c
        write(6,'(/'' SPHBAS: unknown point group '',a4)') idgroup
        stop
c
      end if
c
      do i=1,ndim
      do j=1,ndim
        u(i,j)=dconjg(o(j,i))
      end do
      end do
c
c     write(6,'('' Matrix of basis transformation:'')')
c     call outmat1(o,lmmaxp,lmmaxp,ndim,tol,6)
c     write(6,'(/'' Inverse matrix'')')
c     call outmat1(u,lmmaxp,lmmaxp,ndim,tol,6)
c
c Now get irreducible matrix representations
      call matrep(idgroup,ng,melem,invg,idimrep,rep)
c
      do ig=1,ng
c
        call czero(d,ndim*ndim)
        idim=0
        do ir=1,numrep
          do i=1,idimrep(irep(ir))
          do j=1,idimrep(irep(ir))
            d(idim+i,idim+j)=rep(i,j,irep(ir),ig)
          end do
          end do
          idim=idim+idimrep(irep(ir))
        end do
        if(idim.ne.ndim) then
          write(6,'(/'' SPHBAS: dimension mismatch'',2i3)') idim,ndim
          stop
        end if
c
c       write(6,'(/'' Group element'',i3)') ig
c       call outmat1(d,lmmaxp,lmmaxp,ndim,tol,6)
c
c transform to (l,m) basis
        call tripmt(o,d,u,lmmaxp,lmmaxp,ndim)
c
        do i=1,lmmaxp
        do j=1,lmmaxp
          dg(i,j,ig)=d(i,j)
          dg1(i,j,ig)=dconjg(d(j,i))
        end do 
        end do 
c
      end do
c
      return
      end
