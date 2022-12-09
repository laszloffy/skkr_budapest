c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine matrep(idgroup,ng,melem,invg,idimrep,rep)
c======================
c
c fill up matrix representations of the group identified 
c by idgroup: C1,C2,C3,C4,C6,Cs,C2v,C3vA,C3vB,C4v,C6v
c 
c
c  idgroup: point group identifier
c  ng: order of the point group
c  melem : maximum of order of the point group
c  invg: index for inverse of point group operations
c  idimrep: dimension of irreducible representations
c  rep(i,j,indirr,indelem) :  matrix representations
c      | |   |      |
c      | |   |    index of
c      | |   |    point group
c      | |   |    operation
c      | |   |
c      | |  index of
c      | |  irreducible
c      | |  representation
c      | |
c    indices
c    for symmetry
c    adapted basis
c
      implicit real*8 (a-h,o-z)
c
      character*4 idgroup
c
      dimension invg(melem),idimrep(6)
c
      complex*16 rep(2,2,6,melem)
      complex*16 v,vstar,one,sqrm1
c
      pi=4.0d0*datan(1.0d0)
      one  =dcmplx(1.0d0,0.0d0)
      sqrm1=dcmplx(0.0d0,1.0d0)
      call czero(rep,24*melem) 
c
      if(idgroup.eq.'C1') then
c
c the order of the group is 1
c there is one 1D irreducible representation: A
c
        ng=1
        idimrep(1)=1
c
c group element: E
        invg(1)=1
        rep(1,1,1,1)=1
c
      else if(idgroup.eq.'C2') then
c
c the order of the group is 2
c there are two 1D irreducible representations: A, B
c
        ng=2
        idimrep(1)=1
        idimrep(2)=1
c
c group element: E
        invg(1)=1
        rep(1,1,1,1)=one
        rep(1,1,2,1)=one
c
c group element: C2
        invg(2)=2
        rep(1,1,1,2)=one
        rep(1,1,2,2)=-one
c
      else if(idgroup.eq.'C3') then
c
c the order of the group is 3
c there are three 1D irreducible representations: A, 1E, 2E
c
        ng=3
        idimrep(1)=1
        idimrep(2)=1
        idimrep(3)=1
        v=cdexp(2.d0*pi*sqrm1/3.d0)
        vstar=dconjg(v)
c
c group element: E
        invg(1)=1
        rep(1,1,1,1)=one
        rep(1,1,2,1)=one
        rep(1,1,3,1)=one
c
c group element: C3+
        invg(2)=3
        rep(1,1,1,2)=one
        rep(1,1,2,2)=vstar       
        rep(1,1,3,2)=v
c
c group element: C3-
        invg(3)=2
        rep(1,1,1,3)=one
        rep(1,1,2,3)=v
        rep(1,1,3,3)=vstar    
c
      else if(idgroup.eq.'C4') then
c
c the order of the group is 4
c there are four 1D irreducible representations: A, B, 1E, 2E
c
        ng=4
        idimrep(1)=1
        idimrep(2)=1
        idimrep(3)=1
        idimrep(4)=1
        v=cdexp(2.d0*pi*sqrm1/8.d0)
        vstar=dconjg(v)
c
c group element: E
        invg(1)=1
        rep(1,1,1,1)=one
        rep(1,1,2,1)=one
        rep(1,1,3,1)=one
        rep(1,1,4,1)=one
c
c group element: C4+
        invg(2)=3
        rep(1,1,1,2)=one
        rep(1,1,2,2)=-one    
        rep(1,1,3,2)=-sqrm1
        rep(1,1,4,2)=sqrm1
c
c group element: C4-
        invg(3)=2
        rep(1,1,1,3)=one
        rep(1,1,2,3)=-one  
        rep(1,1,3,3)=sqrm1
        rep(1,1,4,3)=-sqrm1
c
c group element: C2
        invg(4)=4
        rep(1,1,1,4)=one
        rep(1,1,2,4)=one
        rep(1,1,3,4)=-one
        rep(1,1,4,4)=-one
c
      else if(idgroup.eq.'C6') then
c
c the order of the group is 6
c there are six 1D irred. representations: A, B, 1E1, 2E1, 1E2, 2E2
c
        ng=6
        idimrep(1)=1
        idimrep(2)=1
        idimrep(3)=1
        idimrep(4)=1
        idimrep(5)=1
        idimrep(6)=1
        v=cdexp(2.d0*pi*sqrm1/3.d0)
        vstar=dconjg(v)
c
c group element: E
        invg(1)=1
        rep(1,1,1,1)=one
        rep(1,1,2,1)=one
        rep(1,1,3,1)=one
        rep(1,1,4,1)=one
        rep(1,1,5,1)=one
        rep(1,1,6,1)=one
c
c group element: C6+
        invg(2)=3
        rep(1,1,1,2)=one
        rep(1,1,2,2)=-one
        rep(1,1,3,2)=-v
        rep(1,1,4,2)=-vstar
        rep(1,1,5,2)=v
        rep(1,1,6,2)=vstar
c
c group element: C6-
        invg(3)=2
        rep(1,1,1,3)=one
        rep(1,1,2,3)=-one
        rep(1,1,3,3)=-vstar
        rep(1,1,4,3)=-v
        rep(1,1,5,3)=vstar
        rep(1,1,6,3)=v
c
c group element: C3+
        invg(4)=5
        rep(1,1,1,4)=one
        rep(1,1,2,4)=one
        rep(1,1,3,4)=vstar   
        rep(1,1,4,4)=v
        rep(1,1,5,4)=vstar
        rep(1,1,6,4)=v
c
c group element: C3-
        invg(5)=4
        rep(1,1,1,5)=one
        rep(1,1,2,5)=one
        rep(1,1,3,5)=v
        rep(1,1,4,5)=vstar
        rep(1,1,5,5)=v
        rep(1,1,6,5)=vstar
c
c group element: C2
        invg(6)=6
        rep(1,1,1,6)=one
        rep(1,1,2,6)=-one
        rep(1,1,3,6)=-one 
        rep(1,1,4,6)=-one
        rep(1,1,5,6)=one
        rep(1,1,6,6)=one
c
      else if(idgroup.eq.'Cs') then
c
c the order of the group is 2
c there are two 1D irreducible representations: A', A"
c
        ng=2
        idimrep(1)=1
        idimrep(2)=1
c
c group element: E
        invg(1)=1
        rep(1,1,1,1)=1
        rep(1,1,2,1)=1
c
c group element: Sx
        invg(2)=2
        rep(1,1,1,2)=one
        rep(1,1,2,2)=-one
c
      else if(idgroup.eq.'C2v') then
c
c the order of the group is 4
c there are four 1D irreducible representations:  A1, A2, B1, B2
c
        ng=4
        idimrep(1)=1
        idimrep(2)=1
        idimrep(3)=1
        idimrep(4)=1
c
c group element: E
        invg(1)=1
        rep(1,1,1,1)=one
        rep(1,1,2,1)=one
        rep(1,1,3,1)=one
        rep(1,1,4,1)=one
c
c group element: C2
        invg(2)=2
        rep(1,1,1,2)=one
        rep(1,1,2,2)=one
        rep(1,1,3,2)=-one
        rep(1,1,4,2)=-one
c
c group element: Sx
        invg(3)=3
        rep(1,1,1,3)=one
        rep(1,1,2,3)=-one
        rep(1,1,3,3)=-one
        rep(1,1,4,3)=one
c
c group element: Sy
        invg(4)=4
        rep(1,1,1,4)=one
        rep(1,1,2,4)=-one
        rep(1,1,3,4)=one
        rep(1,1,4,4)=-one
c
      else if(idgroup.eq.'C3vA'.or.idgroup.eq.'C3vB') then
c
c the order of the group is 6
c there are two 1D irreducible representations:  A1, A2
c       and one 2D irreducible representation:   E
c
        ng=6
        idimrep(1)=1
        idimrep(2)=1
        idimrep(3)=2
        v=cdexp(2.0d0*sqrm1*pi/3.0d0)
        vstar=dconjg(v)
c
c group element: E
          invg(1)=1
          rep(1,1,1,1)=one
          rep(1,1,2,1)=one
          rep(1,1,3,1)=one
          rep(2,2,3,1)=one
c
c group element: C3+
          invg(2)=3
          rep(1,1,1,2)=one
          rep(1,1,2,2)=one
          rep(1,1,3,2)=vstar
          rep(2,2,3,2)=v
c
c group element: C3-
          invg(3)=2
          rep(1,1,1,3)=one
          rep(1,1,2,3)=one
          rep(1,1,3,3)=v
          rep(2,2,3,3)=vstar
c
c group element: Sv1
          invg(4)=4
          rep(1,1,1,4)=one
          rep(1,1,2,4)=-one
          rep(1,2,3,4)=-one
          rep(2,1,3,4)=-one
c
c group element: Sv2
          invg(5)=5
          rep(1,1,1,5)=one
          rep(1,1,2,5)=-one
          rep(1,2,3,5)=-v
          rep(2,1,3,5)=-vstar
c
c group element: Sv3
          invg(6)=6
          rep(1,1,1,6)=one
          rep(1,1,2,6)=-one
          rep(1,2,3,6)=-vstar
          rep(2,1,3,6)=-v
c
      else if(idgroup.eq.'C4v ') then
c
c the order of the group is 8
c there are four 1D irreducible representations: A1, A2, B1, B2
c       and one 2D irreducible representation:   E  
c
          ng=8
          idimrep(1)=1
          idimrep(2)=1
          idimrep(3)=1
          idimrep(4)=1
          idimrep(5)=2
          v=cdexp(2.0d0*sqrm1*pi/8.0d0)
          vstar=dconjg(v)
c
c group element: E
          invg(1)=1
          rep(1,1,1,1)=one
          rep(1,1,2,1)=one
          rep(1,1,3,1)=one
          rep(1,1,4,1)=one
          rep(1,1,5,1)=one
          rep(2,2,5,1)=one
c
c group element: C4+
          invg(2)=3
          rep(1,1,1,2)=one
          rep(1,1,2,2)=one
          rep(1,1,3,2)=-one
          rep(1,1,4,2)=-one
          rep(1,1,5,2)=-sqrm1
          rep(2,2,5,2)=sqrm1
c
c group element: C4-
          invg(3)=2
          rep(1,1,1,3)=one
          rep(1,1,2,3)=one
          rep(1,1,3,3)=-one
          rep(1,1,4,3)=-one
          rep(1,1,5,3)=sqrm1
          rep(2,2,5,3)=-sqrm1
c
c group element: C2
          invg(4)=4
          rep(1,1,1,4)=one
          rep(1,1,2,4)=one
          rep(1,1,3,4)=one
          rep(1,1,4,4)=one
          rep(1,1,5,4)=-one
          rep(2,2,5,4)=-one
c
c group element: Sv1
          invg(5)=5
          rep(1,1,1,5)=one
          rep(1,1,2,5)=-one
          rep(1,1,3,5)=one
          rep(1,1,4,5)=-one
c         rep(1,2,5,5)=-one
c         rep(1,2,5,5)=-one
          rep(1,2,5,5)=-one
          rep(2,1,5,5)=-one
c
c group element: Sv2
          invg(6)=6
          rep(1,1,1,6)=one
          rep(1,1,2,6)=-one
          rep(1,1,3,6)=one
          rep(1,1,4,6)=-one
c         rep(1,2,5,6)=one
c         rep(1,2,5,6)=one
          rep(1,2,5,6)=one
          rep(2,1,5,6)=one
c
c group element: Sd1
          invg(7)=7
          rep(1,1,1,7)=one
          rep(1,1,2,7)=-one
          rep(1,1,3,7)=-one
          rep(1,1,4,7)=one
c         rep(1,2,5,7)=sqrm1
c         rep(1,2,5,7)=-sqrm1
          rep(1,2,5,7)=sqrm1
          rep(2,1,5,7)=-sqrm1
c
c group element: Sd2
          invg(8)=8
          rep(1,1,1,8)=one
          rep(1,1,2,8)=-one
          rep(1,1,3,8)=-one
          rep(1,1,4,8)=one
c         rep(1,2,5,8)=-sqrm1
c         rep(1,2,5,8)=sqrm1
          rep(1,2,5,8)=-sqrm1
          rep(2,1,5,8)=sqrm1
c
      else if(idgroup.eq.'C6v ') then
c
c the order of the group is 12
c there are four 1D irreducible representations: A1, A2, B1, B2
c       and two 2D irreducible representation:   E1, E2  
c
          ng=12
          idimrep(1)=1
          idimrep(2)=1
          idimrep(3)=1
          idimrep(4)=1
          idimrep(5)=2
          idimrep(6)=2
          v=cdexp(2.0d0*sqrm1*pi/3.0d0)
          vstar=dconjg(v)
c
c group element: E
          invg(1)=1
          rep(1,1,1,1)=one
          rep(1,1,2,1)=one
          rep(1,1,3,1)=one
          rep(1,1,4,1)=one
          rep(1,1,5,1)=one
          rep(2,2,5,1)=one
          rep(1,1,6,1)=one
          rep(2,2,6,1)=one
c
c group element: C6+
          invg(2)=3
          rep(1,1,1,2)=one
          rep(1,1,2,2)=one
          rep(1,1,3,2)=-one
          rep(1,1,4,2)=-one
          rep(1,1,5,2)=-v 
          rep(2,2,5,2)=-vstar
          rep(1,1,6,2)=v
          rep(2,2,6,2)=vstar
c
c group element: C6-
          invg(3)=2
          rep(1,1,1,3)=one
          rep(1,1,2,3)=one
          rep(1,1,3,3)=-one
          rep(1,1,4,3)=-one
          rep(1,1,5,3)=-vstar 
          rep(2,2,5,3)=-v
          rep(1,1,6,3)=vstar
          rep(2,2,6,3)=v
c
c group element: C3+
          invg(4)=5
          rep(1,1,1,4)=one
          rep(1,1,2,4)=one
          rep(1,1,3,4)=one
          rep(1,1,4,4)=one
          rep(1,1,5,4)=vstar 
          rep(2,2,5,4)=v
          rep(1,1,6,4)=vstar
          rep(2,2,6,4)=v
c
c group element: C3-
          invg(5)=4
          rep(1,1,1,5)=one
          rep(1,1,2,5)=one
          rep(1,1,3,5)=one
          rep(1,1,4,5)=one
          rep(1,1,5,5)=v
          rep(2,2,5,5)=vstar
          rep(1,1,6,5)=v
          rep(2,2,6,5)=vstar
c
c group element: C2
          invg(6)=6
          rep(1,1,1,6)=one
          rep(1,1,2,6)=one
          rep(1,1,3,6)=-one
          rep(1,1,4,6)=-one
          rep(1,1,5,6)=-one
          rep(2,2,5,6)=-one 
          rep(1,1,6,6)=one
          rep(2,2,6,6)=one
c
c group element: Sd1
          invg(7)=7
          rep(1,1,1,7)=one
          rep(1,1,2,7)=-one
          rep(1,1,3,7)=-one
          rep(1,1,4,7)=one
          rep(1,2,5,7)=-one
          rep(2,1,5,7)=-one 
          rep(1,2,6,7)=-one
          rep(2,1,6,7)=-one
c
c group element: Sd2
          invg(8)=8
          rep(1,1,1,8)=one
          rep(1,1,2,8)=-one
          rep(1,1,3,8)=-one
          rep(1,1,4,8)=one
          rep(1,2,5,8)=-v
          rep(2,1,5,8)=-vstar
          rep(1,2,6,8)=-v
          rep(2,1,6,8)=-vstar
c
c group element: Sd3
          invg(9)=9
          rep(1,1,1,9)=one
          rep(1,1,2,9)=-one
          rep(1,1,3,9)=-one
          rep(1,1,4,9)=one
          rep(1,2,5,9)=-vstar
          rep(2,1,5,9)=-v
          rep(1,2,6,9)=-vstar
          rep(2,1,6,9)=-v
c
c group element: Sv1
          invg(10)=10
          rep(1,1,1,10)=one
          rep(1,1,2,10)=-one
          rep(1,1,3,10)=one
          rep(1,1,4,10)=-one
          rep(1,2,5,10)=one
          rep(2,1,5,10)=one 
          rep(1,2,6,10)=-one
          rep(2,1,6,10)=-one
c
c group element: Sv2
          invg(11)=11
          rep(1,1,1,11)=one
          rep(1,1,2,11)=-one
          rep(1,1,3,11)=one
          rep(1,1,4,11)=-one
          rep(1,2,5,11)=v
          rep(2,1,5,11)=vstar
          rep(1,2,6,11)=-v
          rep(2,1,6,11)=-vstar
c
c group element: Sv3
          invg(12)=12
          rep(1,1,1,12)=one
          rep(1,1,2,12)=-one
          rep(1,1,3,12)=one
          rep(1,1,4,12)=-one
          rep(1,2,5,12)=vstar
          rep(2,1,5,12)=v
          rep(1,2,6,12)=-vstar
          rep(2,1,6,12)=-v
c
      else 
c
          write(6,'(/'' GROUP: unknown point group '',a4)') idgroup
          stop
c
      end if
c
      return
      end
