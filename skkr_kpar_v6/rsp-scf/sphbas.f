c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine sphbas(dg,dg1,ng,invg)
c======================
c
c Calculates the matrix representations of the double groups:
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
      parameter (ndim=50)
c
      character*4 idgroup
c
      dimension invg(melem),idimrep(6),irep(ndim)
c
      complex*16 dg(kmymaxp,kmymaxp,melem)
      complex*16 dg1(kmymaxp,kmymaxp,melem)
      complex*16 o(ndim,ndim),u(ndim,ndim)
      complex*16 d(ndim,ndim)
      complex*16 rep(2,2,6,melem)
      complex*16 one,sqrm1,facn,faci
c
      common/latt2d/idgroup,ilat,a2d,b2d,d2d
      data tol/1.0d-14/
c
c     write(6,*) idgroup,ilat,a2d,b2d,d2d
c     call flush(6)
      if(ndim.lt.kmymaxp) then
        write(6,*) 'SPHBAS: Increase parameter ndim to',kmymaxp
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
c l=0 j=1/2
        irep(1)= 1
        irep(2)= 2
c l=1 j=1/2
        irep(3)= 1
        irep(4)= 2
c l=1 j=3/2
        irep(5)= 2
        irep(6)= 1
        irep(7)= 2
        irep(8)= 1
c l=2 j=3/2
        irep(9)= 2
        irep(10)=1
        irep(11)=2
        irep(12)=1
c l=2 j=5/2
        irep(13)=1
        irep(14)=2
        irep(15)=1
        irep(16)=2
        irep(17)=1
        irep(18)=2
c l=3 j=5/2
        irep(19)=1
        irep(20)=2
        irep(21)=1
        irep(22)=2
        irep(23)=1
        irep(24)=2
c l=3 j=7/2
        irep(25)=2
        irep(26)=1
        irep(27)=2
        irep(28)=1
        irep(29)=2
        irep(30)=1
        irep(31)=2
        irep(32)=1
c l=4 j=7/2
        irep(33)=2
        irep(34)=1
        irep(35)=2
        irep(36)=1
        irep(37)=2
        irep(38)=1
        irep(39)=2
        irep(40)=1
c l=4 j=9/2
        irep(41)=1
        irep(42)=2
        irep(43)=1
        irep(44)=2
        irep(45)=1
        irep(46)=2
        irep(47)=1
        irep(48)=2
        irep(49)=1
        irep(50)=2
c
      else if(idgroup.eq.'C3') then
c
        do i=1,ndim
          o(i,i)=1.0d0
        end do
        numrep=ndim
c l=0 j=1/2
        irep(1)= 1
        irep(2)= 2
c l=1 j=1/2
        irep(3)= 1
        irep(4)= 2
c l=1 j=3/2
        irep(5)= 3
        irep(6)= 1
        irep(7)= 2
        irep(8)= 3
c l=2 j=3/2
        irep(9)= 3
        irep(10)=1
        irep(11)=2
        irep(12)=3
c l=2 j=5/2
        irep(13)=2
        irep(14)=3
        irep(15)=1
        irep(16)=2
        irep(17)=3
        irep(18)=1
c l=3 j=5/2
        irep(19)=2
        irep(20)=3
        irep(21)=1
        irep(22)=2
        irep(23)=3
        irep(24)=1
c l=3 j=7/2
        irep(25)=1
        irep(26)=2
        irep(27)=3
        irep(28)=1
        irep(29)=2
        irep(30)=3
        irep(31)=1
        irep(32)=2
c l=4 j=7/2
        irep(33)=1
        irep(34)=2
        irep(35)=3
        irep(36)=1
        irep(37)=2
        irep(38)=3
        irep(39)=1
        irep(40)=2
c l=4 j=9/2
        irep(41)=3
        irep(42)=1
        irep(43)=2
        irep(44)=3
        irep(45)=1
        irep(46)=2
        irep(47)=3
        irep(48)=1
        irep(49)=2
        irep(50)=3
c
      else if(idgroup.eq.'C4') then
c
        do i=1,ndim
          o(i,i)=1.0d0
        end do
        numrep=ndim
c l=0 j=1/2
        irep(1)= 2
        irep(2)= 1
c l=1 j=1/2
        irep(3)= 2
        irep(4)= 1
c l=1 j=3/2
        irep(5)= 3
        irep(6)= 2
        irep(7)= 1
        irep(8)= 4
c l=2 j=3/2
        irep(9)= 3
        irep(10)=2
        irep(11)=1
        irep(12)=4
c l=2 j=5/2
        irep(13)=4
        irep(14)=3
        irep(15)=2
        irep(16)=1
        irep(17)=4
        irep(18)=3
c l=3 j=5/2
        irep(19)=4
        irep(20)=3
        irep(21)=2
        irep(22)=1
        irep(23)=4
        irep(24)=3
c l=3 j=7/2
        irep(25)=1
        irep(26)=4
        irep(27)=3
        irep(28)=2
        irep(29)=1
        irep(30)=4
        irep(31)=3
        irep(32)=2
c l=4 j=7/2
        irep(33)=1
        irep(34)=4
        irep(35)=3
        irep(36)=2
        irep(37)=1
        irep(38)=4
        irep(39)=3
        irep(40)=2
c l=4 j=9/2
        irep(41)=2
        irep(42)=1
        irep(43)=4
        irep(44)=3
        irep(45)=2
        irep(46)=1
        irep(47)=4
        irep(48)=3
        irep(49)=2
        irep(50)=1
c
      else if(idgroup.eq.'C6') then
c
        do i=1,ndim
          o(i,i)=1.0d0
        end do
        numrep=ndim
c l=0 j=1/2
        irep(1)= 1
        irep(2)= 2
c l=1 j=1/2
        irep(3)= 1
        irep(4)= 2
c l=1 j=3/2
        irep(5)= 4
        irep(6)= 1
        irep(7)= 2
        irep(8)= 3
c l=2 j=3/2
        irep(9)= 4
        irep(10)=1
        irep(11)=2
        irep(12)=3
c l=2 j=5/2
        irep(13)=5
        irep(14)=4
        irep(15)=1
        irep(16)=2
        irep(17)=3
        irep(18)=6
c l=3 j=5/2
        irep(19)=5
        irep(20)=4
        irep(21)=1
        irep(22)=2
        irep(23)=3
        irep(24)=6
c l=3 j=7/2
        irep(25)=6
        irep(26)=5
        irep(27)=4
        irep(28)=1
        irep(29)=2
        irep(30)=3
        irep(31)=6
        irep(32)=5
c l=4 j=7/2
        irep(33)=6
        irep(34)=5
        irep(35)=4
        irep(36)=1
        irep(37)=2
        irep(38)=3
        irep(39)=6
        irep(40)=5
c l=4 j=9/2
        irep(41)=3
        irep(42)=6
        irep(43)=5
        irep(44)=4
        irep(45)=1
        irep(46)=2
        irep(47)=3
        irep(48)=6
        irep(49)=5
        irep(50)=4
c
      else if(idgroup.eq.'Cs') then
c
c l=0 j=1/2
        irep(1)= 1
          o(1,1)=facn
          o(2,1)=facn
        irep(2)= 2
          o(1,2)=-facn
          o(2,2)=facn
c l=1 j=1/2
        irep(3)= 2
          o(3,3)=facn
          o(4,3)=facn
        irep(4)= 1
          o(3,4)=-facn
          o(4,4)=facn
c l=1 j=3/2
        irep(5)= 2
          o(5,5)=facn
          o(8,5)=facn
        irep(6)= 1
          o(5,6)=-facn
          o(8,6)=facn
        irep(7)= 1
          o(6,7)=facn
          o(7,7)=facn
        irep(8)= 2
          o(6,8)=-facn
          o(7,8)=facn
c l=2 j=3/2
        irep(9)= 1
          o(9,9)=facn
          o(12,9)=facn
        irep(10)=2
          o(9,10)=-facn
          o(12,10)=facn
        irep(11)=2
          o(10,11)=facn
          o(11,11)=facn
        irep(12)=1
          o(10,12)=-facn
          o(11,12)=facn
c l=2 j=5/2
        irep(13)=1
          o(13,13)=facn
          o(18,13)=facn
        irep(14)=2
          o(13,14)=-facn
          o(18,14)=facn
        irep(15)=2
          o(14,15)=facn
          o(17,15)=facn
        irep(16)=1
          o(14,16)=-facn
          o(17,16)=facn
        irep(17)=1
          o(15,17)=facn
          o(16,17)=facn
        irep(18)=2
          o(15,18)=-facn
          o(16,18)=facn
c l=3 j=5/2
        irep(19)=2
          o(19,19)=facn
          o(24,19)=facn
        irep(20)=1
          o(19,20)=-facn
          o(24,20)=facn
        irep(21)=1
          o(20,21)=facn
          o(23,21)=facn
        irep(22)=2
          o(20,22)=-facn
          o(23,22)=facn
        irep(23)=2
          o(21,23)=facn
          o(22,23)=facn
        irep(24)=1
          o(21,24)=-facn
          o(22,24)=facn
c l=3 j=7/2
        irep(25)=2
          o(25,25)=facn
          o(32,25)=facn
        irep(26)=1
          o(25,26)=-facn
          o(32,26)=facn
        irep(27)=1
          o(26,27)=facn
          o(31,27)=facn
        irep(28)=2
          o(26,28)=-facn
          o(31,28)=facn
        irep(29)=2
          o(27,29)=facn
          o(30,29)=facn
        irep(30)=1
          o(27,30)=-facn
          o(30,30)=facn
        irep(31)=1
          o(28,31)=facn
          o(29,31)=facn
        irep(32)=2
          o(28,32)=-facn
          o(29,32)=facn
c l=4 j=7/2
        irep(33)=1
          o(33,33)=facn
          o(40,33)=facn
        irep(34)=2
          o(33,34)=-facn
          o(40,34)=facn
        irep(35)=2
          o(34,35)=facn
          o(39,35)=facn
        irep(36)=1
          o(34,36)=-facn
          o(39,36)=facn
        irep(37)=1
          o(35,37)=facn
          o(38,37)=facn
        irep(38)=2
          o(35,38)=-facn
          o(38,38)=facn
        irep(39)=2
          o(36,39)=facn
          o(37,39)=facn
        irep(40)=1
          o(36,40)=-facn
          o(37,40)=facn
c l=4 j=9/2
        irep(41)=1
          o(41,41)=facn
          o(50,41)=facn
        irep(42)=2
          o(41,42)=-facn
          o(50,42)=facn
        irep(43)=2
          o(42,43)=facn
          o(49,43)=facn
        irep(44)=1
          o(42,44)=-facn
          o(49,44)=facn
        irep(45)=1
          o(43,45)=facn
          o(48,45)=facn
        irep(46)=2
          o(43,46)=-facn
          o(48,46)=facn
        irep(47)=2
          o(44,47)=facn
          o(47,47)=facn
        irep(48)=1
          o(44,48)=-facn
          o(47,48)=facn
        irep(49)=1
          o(45,49)=facn
          o(46,49)=facn
        irep(50)=2
          o(45,50)=-facn
          o(46,50)=facn
c
        numrep=ndim
c
      else if(idgroup.eq.'C2v') then
c         
c l=0 j=1/2
        irep(1)=1
          o(2,1)=one
          o(1,2)=one
c l=1 j=1/2
        irep(2)=1
          o(4,3)=one
          o(3,4)=-one
c l=1 j=3/2
        irep(3)=1
          o(7,5)=one
          o(6,6)=one
        irep(4)=1
          o(5,7)=one
          o(8,8)=one
c l=2 j=3/2
        irep(5)=1
          o(11,9)=one
          o(10,10)=-one
        irep(6)=1
          o(9,11)=one
          o(12,12)=-one
c l=2 j=5/2
        irep(7)=1
          o(16,13)=one
          o(15,14)=one
        irep(8)=1
          o(14,15)=one
          o(17,16)=one
        irep(9)=1
          o(18,17)=one
          o(13,18)=one
c l=3 j=5/2
        irep(10)=1
          o(22,19)=one
          o(21,20)=-one
        irep(11)=1
          o(20,21)=one
          o(23,22)=-one
        irep(12)=1
          o(24,23)=one
          o(19,24)=-one
c l=3 j=7/2
        irep(13)=1
          o(29,25)=one
          o(28,26)=one
        irep(14)=1
          o(27,27)=one
          o(30,28)=one
        irep(15)=1
          o(31,29)=one
          o(26,30)=one
        irep(16)=1
          o(25,31)=one
          o(32,32)=one
c l=4 j=7/2
        irep(17)=1
          o(37,33)=one
          o(36,34)=-one
        irep(18)=1
          o(35,35)=one
          o(38,36)=-one
        irep(19)=1
          o(39,37)=one
          o(34,38)=-one
        irep(20)=1
          o(33,39)=one
          o(40,40)=-one
c l=4 j=9/2
        irep(21)=1
          o(46,41)=one
          o(47,42)=one
        irep(22)=1
          o(44,43)=one
          o(47,44)=one
        irep(23)=1
          o(48,45)=one
          o(43,46)=one
        irep(24)=1
          o(42,47)=one
          o(49,48)=one
        irep(25)=1
          o(50,49)=one
          o(41,50)=one
c
        numrep=25
c
      else if(idgroup.eq.'C3vA') then
c
c l=0 j=1/2
        irep(1)=1
          o(2,1)=one
          o(1,2)=one
c l=1 j=1/2
        irep(2)=1
          o(4,3)=one
          o(3,4)=-one
c l=1 j=3/2
        irep(3)=1
          o(7,5)=one
          o(6,6)=one
        irep(4)=2
          o(5,7)=faci
          o(8,7)=facn
        irep(5)=3
          o(5,8)=-faci
          o(8,8)=facn
c l=2 j=3/2
        irep(6)=1
          o(11,9)=one
          o(10,10)=-one
        irep(7)=2
          o(9,11)=-faci
          o(12,11)=facn
        irep(8)=3
          o(9,12)=faci
          o(12,12)=facn
c l=2 j=5/2
        irep(9)=1
          o(16,13)=one
          o(15,14)=one
        irep(10)=2
          o(14,15)=faci
          o(17,15)=facn
        irep(11)=3
          o(14,16)=-faci
          o(17,16)=facn
        irep(12)=1
          o(13,17)=one
          o(18,18)=-one
c l=3 j=5/2
        irep(13)=1
          o(22,19)=one
          o(21,20)=-one
        irep(14)=2
          o(20,21)=-faci
          o(23,21)=facn
        irep(15)=3
          o(20,22)=faci
          o(23,22)=facn
        irep(16)=1
          o(19,23)=one
          o(24,24)=one
c l=3 j=7/2
        irep(17)=1
          o(29,25)=one
          o(28,26)=one
        irep(18)=2
          o(27,27)=faci
          o(30,27)=facn
        irep(19)=3
          o(27,28)=-faci
          o(30,28)=facn
        irep(20)=1
          o(26,29)=one
          o(31,30)=-one
        irep(21)=1
          o(32,31)=one
          o(25,32)=-one
c l=4 j=7/2
        irep(22)=1
          o(37,33)=one
          o(36,34)=-one
        irep(23)=2
          o(35,35)=-faci
          o(38,35)=facn
        irep(24)=3
          o(35,36)=faci
          o(38,36)=facn
        irep(25)=1
          o(34,37)=one
          o(39,38)=one
        irep(26)=1
          o(40,39)=one
          o(33,40)=one
c l=4 j=9/2
        irep(27)=1
          o(46,41)=one
          o(45,42)=one
        irep(28)=2
          o(44,43)=faci
          o(47,43)=facn
        irep(29)=3
          o(44,44)=-faci
          o(47,44)=facn
        irep(30)=1
          o(43,45)=one
          o(48,46)=-one
        irep(31)=1
          o(49,47)=one
          o(42,48)=-one
        irep(32)=2
          o(41,49)=facn
          o(50,49)=faci
        irep(33)=3
          o(41,50)=facn
          o(50,50)=-faci
c
        numrep=33
c
      else if(idgroup.eq.'C3vB') then
c
c l=0 j=1/2
        irep(1)=1
          o(2,1)=one
          o(1,2)=one
c l=1 j=1/2
        irep(2)=1
          o(4,3)=one
          o(3,4)=-one
c l=1 j=3/2
        irep(3)=1
          o(7,5)=one
          o(6,6)=one
        irep(4)=2
          o(5,7)=-facn
          o(8,7)=facn
        irep(5)=3
          o(5,8)=facn
          o(8,8)=facn
c l=2 j=3/2
        irep(6)=1
          o(11,9)=one
          o(10,10)=-one
        irep(7)=2
          o(9,11)=facn
          o(12,11)=facn
        irep(8)=3
          o(9,12)=-facn
          o(12,12)=facn
c l=2 j=5/2
        irep(9)=1
          o(16,13)=one
          o(15,14)=one
        irep(10)=2
          o(14,15)=-facn
          o(17,15)=facn
        irep(11)=3
          o(14,16)=facn
          o(17,16)=facn
        irep(12)=1
          o(13,17)=one
          o(18,18)=one
c l=3 j=5/2
        irep(13)=1
          o(22,19)=one
          o(21,20)=-one
        irep(14)=2
          o(20,21)=facn
          o(23,21)=facn
        irep(15)=3
          o(20,22)=-facn
          o(23,22)=facn
        irep(16)=1
          o(19,23)=one
          o(24,24)=-one
c l=3 j=7/2
        irep(17)=1
          o(29,25)=one
          o(28,26)=one
        irep(18)=2
          o(27,27)=-facn
          o(30,27)=facn
        irep(19)=3
          o(27,28)=facn
          o(30,28)=facn
        irep(20)=1
          o(26,29)=one
          o(31,30)=one
        irep(21)=1
          o(32,31)=one
          o(25,32)=one
c l=4 j=7/2
        irep(22)=1
          o(37,33)=one
          o(36,34)=-one
        irep(23)=2
          o(35,35)=facn
          o(38,35)=facn
        irep(24)=3
          o(35,36)=-facn
          o(38,36)=facn
        irep(25)=1
          o(34,37)=one
          o(39,38)=-one
        irep(26)=1
          o(40,39)=one
          o(33,40)=-one
c l=4 j=9/2
        irep(27)=1
          o(46,41)=one
          o(45,42)=one
        irep(28)=2
          o(44,43)=-facn
          o(47,43)=facn
        irep(29)=3
          o(44,44)=facn
          o(47,44)=facn
        irep(30)=1
          o(43,45)=one
          o(48,46)=one
        irep(31)=1
          o(49,47)=one
          o(42,48)=one
        irep(32)=2
          o(41,49)=-facn
          o(50,49)=facn
        irep(33)=3
          o(41,50)=facn
          o(50,50)=facn
c
        numrep=33
c
      else if(idgroup.eq.'C4v') then
c 
c l=0 j=1/2
        irep(1)=1
          o(2,1)=one
          o(1,2)=one
c l=1 j=1/2
        irep(2)=1
          o(4,3)=one
          o(3,4)=-one
c l=1 j=3/2
        irep(3)=1
          o(7,5)=one
          o(6,6)=one
        irep(4)=2
          o(5,7)=one
          o(8,8)=one
c l=2 j=3/2
        irep(5)=1
          o(11,9)=one
          o(10,10)=-one
        irep(6)=2
          o(9,11)=one
          o(12,12)=-one
c l=2 j=5/2
        irep(7)=1
          o(16,13)=one
          o(15,14)=one
        irep(8)=2
          o(14,15)=one
          o(17,16)=one
        irep(9)=2
          o(18,17)=one
          o(13,18)=one
c l=3 j=5/2
        irep(10)=1
          o(22,19)=one
          o(21,20)=-one
        irep(11)=2
          o(20,21)=one
          o(23,22)=-one
        irep(12)=2
          o(24,23)=one
          o(19,24)=-one
c l=3 j=7/2
        irep(13)=1
          o(29,25)=one
          o(28,26)=one
        irep(14)=2
          o(27,27)=one
          o(30,28)=one
        irep(15)=2
          o(31,29)=one
          o(26,30)=one
        irep(16)=1
          o(25,31)=one
          o(32,32)=one
c l=4 j=7/2
        irep(17)=1
          o(37,33)=one
          o(36,34)=-one
        irep(18)=2
          o(35,35)=one
          o(38,36)=-one
        irep(19)=2
          o(39,37)=one
          o(34,38)=-one
        irep(20)=1
          o(33,39)=one
          o(40,40)=-one
c l=4 j=9/2
        irep(21)=1
          o(46,41)=one
          o(45,42)=one
        irep(22)=2
          o(44,43)=one
          o(47,44)=one
        irep(23)=2
          o(48,45)=one
          o(43,46)=one
        irep(24)=1
          o(42,47)=one
          o(49,48)=one
        irep(25)=1
          o(50,49)=one
          o(41,50)=one
c
        numrep=25
c
      else if(idgroup.eq.'C6v') then
c
c l=0 j=1/2
        irep(1)=1
          o(2,1)=one
          o(1,2)=one
c l=1 j=1/2
        irep(2)=1
          o(4,3)=one
          o(3,4)=-one
c l=1 j=3/2
        irep(3)=1
          o(7,5)=one
          o(6,6)=one
        irep(4)=2
          o(8,7)=one
          o(5,8)=-one 
c l=2 j=3/2
        irep(5)=1
          o(11,9)=one
          o(10,10)=-one
        irep(6)=2
          o(12,11)=one 
          o(9,12)=one 
c l=2 j=5/2
        irep(7)=1
          o(16,13)=one
          o(15,14)=one
        irep(8)=2
          o(17,15)=one
          o(14,16)=-one
        irep(9)=3
          o(13,17)=one
          o(18,18)=one
c l=3 j=5/2
        irep(10)=1
          o(22,19)=one
          o(21,20)=-one
        irep(11)=2
          o(23,21)=one 
          o(20,22)=one 
        irep(12)=3
          o(19,23)=one
          o(24,24)=-one
c l=3 j=7/2
        irep(13)=1
          o(29,25)=one
          o(28,26)=one
        irep(14)=2
          o(30,27)=one
          o(27,28)=-one
        irep(15)=3
          o(26,29)=one
          o(31,30)=one
        irep(16)=3
          o(32,31)=one
          o(25,32)=one
c l=4 j=7/2
        irep(17)=1
          o(37,33)=one
          o(36,34)=-one
        irep(18)=2
          o(38,35)=one
          o(35,36)=one
        irep(19)=3
          o(34,37)=one
          o(39,38)=-one
        irep(20)=3
          o(40,39)=one
          o(33,40)=-one
c l=4 j=9/2
        irep(21)=1
          o(46,41)=one
          o(45,42)=one
        irep(22)=2
          o(47,43)=one
          o(44,44)=-one
        irep(23)=3
          o(43,45)=one
          o(48,46)=one
        irep(24)=3
          o(49,47)=one
          o(42,48)=one
        irep(25)=2
          o(41,49)=one
          o(50,50)=-one
c
        numrep=25
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
c     call outmat1(o,kmymaxp,kmymaxp,ndim,tol,6)
c     write(6,'(/'' Inverse matrix'')')
c     call outmat1(u,kmymaxp,kmymaxp,ndim,tol,6)
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
c       call outmat1(d,kmymaxp,kmymaxp,ndim,tol,6)
c
c transform to (kap,my) basis
        call tripmt(o,d,u,kmymaxp,kmymaxp,ndim)
c
        do i=1,kmymaxp
        do j=1,kmymaxp
          dg(i,j,ig)=d(i,j)
          dg1(i,j,ig)=dconjg(d(j,i))
        end do 
        end do 
c
      end do
c
      return
      end
