c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine matr2(tlm,lmax,dmat,dmat1,ddph,ddph1,d2dph,d2dph1,
     >                 ddth,ddth1,d2dth,d2dth1,d2dthph,d2dthph1,d)
c ===================
c
c R=(tvec,Phi)
c dmat = D(R) and dmat1 = D(R)+ 
c
      implicit real*8 (a-h,o-z)
      include '../param.h'
c     parameter(jdim=2*lmaxp+2)
c
      dimension tlm(4)
      complex*16 dmat(kmymaxp,kmymaxp),dmat1(kmymaxp,kmymaxp)
      complex*16 ddph(kmymaxp,kmymaxp),ddph1(kmymaxp,kmymaxp)
      complex*16 d2dph(kmymaxp,kmymaxp),d2dph1(kmymaxp,kmymaxp)
      complex*16 ddth(kmymaxp,kmymaxp),ddth1(kmymaxp,kmymaxp)
      complex*16 d2dth(kmymaxp,kmymaxp),d2dth1(kmymaxp,kmymaxp)
      complex*16 d2dthph(kmymaxp,kmymaxp),d2dthph1(kmymaxp,kmymaxp)
      complex*16 ddthe(kmaxp+1,kmaxp+1),d2dthe(kmaxp+1,kmaxp+1)
      complex*16 cmat(kmymaxp,kmymaxp)
      complex*16 d(kmaxp+1,kmaxp+1),zi
      parameter (zi=(0.d0,1.d0))
      common/test/itest
c
      kmax=2*lmax+1
      kmymax=2*(lmax+1)*(lmax+1)
c
c Set up matrix of rotation
c
      do i=1,kmymax
      do j=1,kmymax
        dmat(i,j)=(0.d0,0.d0)
      end do
      end do
      ist=0
      do j2=1,2*lmax-1,2
        il=j2+1
	call rotmat2(d,trd,ddthe,d2dthe,il,tlm,1)
	do icase=1,2
          do m1=1,il
          do m2=1,il
            dmat(ist+m1,ist+m2)=d(il+1-m1,il+1-m2)
            ddph(ist+m1,ist+m2)=zi*(m2-m1)*d(il+1-m1,il+1-m2)
            d2dph(ist+m1,ist+m2)=-(m2-m1)**2*d(il+1-m1,il+1-m2)
            ddth(ist+m1,ist+m2)=ddthe(il+1-m1,il+1-m2)
            d2dthph(ist+m1,ist+m2)=zi*(m2-m1)*ddthe(il+1-m1,il+1-m2)
            d2dth(ist+m1,ist+m2)=d2dthe(il+1-m1,il+1-m2)
          end do
          end do
          ist=ist+il
        end do
      end do
      j2=2*lmax+1
      il=j2+1
      call rotmat2(d,trd,ddthe,d2dthe,il,tlm,1)
      do m1=1,il
      do m2=1,il
c fill up matrix according to ascending m index
        dmat(ist+m1,ist+m2)=d(il+1-m1,il+1-m2)
        ddph(ist+m1,ist+m2)=zi*(m2-m1)*d(il+1-m1,il+1-m2)
        d2dph(ist+m1,ist+m2)=-(m2-m1)**2*d(il+1-m1,il+1-m2)
        ddth(ist+m1,ist+m2)=ddthe(il+1-m1,il+1-m2)
        d2dthph(ist+m1,ist+m2)=zi*(m2-m1)*ddthe(il+1-m1,il+1-m2)
        d2dth(ist+m1,ist+m2)=d2dthe(il+1-m1,il+1-m2)
      end do
      end do
c
      do i=1,kmymax
      do j=1,kmymax
        dmat1(i,j)=dconjg(dmat(j,i))
        ddph1(i,j)=dconjg(ddph(j,i))
        d2dph1(i,j)=dconjg(d2dph(j,i))
        ddth1(i,j)=dconjg(ddth(j,i))
        d2dthph1(i,j)=dconjg(d2dthph(j,i))
        d2dth1(i,j)=dconjg(d2dth(j,i))
      end do
      end do
c
      if(itest.lt.4) return
      write(6,'(/'' Matrix of rotation'')') 
      call outmat(dmat,kmymax,kmymax,kmymaxp,6)
      write(6,'(/'' Inverse'')') 
      call outmat(dmat1,kmymax,kmymax,kmymaxp,6)
      call repl(cmat,dmat,kmymax,kmymaxp)
      call doubmt(cmat,dmat1,kmymax,kmymaxp)
      write(6,'(/'' D(R) * D(R**-1) '')')
      call outmat(cmat,kmymax,kmymax,kmymaxp,6)
      write(6,'(/'' dR/dphi'')') 
      call outmat(ddph,kmymax,kmymax,kmymaxp,6)
      write(6,'(/'' d2R/dphi2'')') 
      call outmat(d2dph,kmymax,kmymax,kmymaxp,6)
      write(6,'(/'' dR/dtheta'')') 
      call outmat(ddth,kmymax,kmymax,kmymaxp,6)
      write(6,'(/'' d2R/dtheta2'')') 
      call outmat(d2dth,kmymax,kmymax,kmymaxp,6)
      write(6,'(/'' d2R/dthetadphi'')') 
      call outmat(d2dthph,kmymax,kmymax,kmymaxp,6)
c
      return
      end
