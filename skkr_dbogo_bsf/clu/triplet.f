c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine tripletcalc(lmax,lms,coefg,coeff,sig,ggrmat,fgrmat,rb,
     >                       tripletmat)
c======================
c
c         Temporary function to print out the whole G matrix        
c input:  lmax  - maximum of angular momentum index
c         lms - if .true. output in (l,m,s) representation
c         coefg - angular momentum matrix for the large components
c         coeff - angular momentum matrix for the small components
c         sig - sign of the contribution from the small components
c         ggrmat, fgrmat - matrices of the radial integrals of the Green function
c         rb - orientation of exchange field, z in local mode, general in global mode
c output: tripletmat
      implicit none
c
      include '../param.h'
c
      logical lms
      integer lmax
      real*8 sig,rb(3)
      complex*16 coefg(kmymaxp,kmymaxp)
      complex*16 coeff(kmymaxp,kmymaxp)
      complex*16 ggrmat(kmymaxp,kmymaxp)
      complex*16 fgrmat(kmymaxp,kmymaxp)
      complex*16 zmom(kmymaxp)
      complex*16 utr(lmmaxp,lmmaxp)
      complex*16 utrp(lmmaxp,lmmaxp)
c
      integer kmymax,i,nl,nl2
      integer orb1, orb2
      real*8 pi,fac
      complex*16 zgmat(kmymaxp,kmymaxp)
      complex*16 zfmat(kmymaxp,kmymaxp)
      complex*16 zmat(kmymaxp,kmymaxp)
      complex*16 zlms(kmymaxp,kmymaxp)
      real*8 tripletmat(lmmaxp,lmmaxp,3)
c
      real*8 tol
      data tol/1.d-8/
c
      nl=lmax+1
      nl2=nl*nl
      pi=4.d0*datan(1.d0)
      fac=-1.d0/pi
      kmymax=2*(lmax+1)*(lmax+1)
c
      zgmat(1:kmymax,1:kmymax)=
     >matmul(coefg(1:kmymax,1:kmymax),ggrmat(1:kmymax,1:kmymax))
      zfmat(1:kmymax,1:kmymax)=
     >matmul(coeff(1:kmymax,1:kmymax),fgrmat(1:kmymax,1:kmymax))
      zmat(1:kmymax,1:kmymax)=fac*(zgmat(1:kmymax,1:kmymax)
     >                        +sig*zfmat(1:kmymax,1:kmymax))
c
c find diagonals in the required representation
c
      if (lms) then 
        if (localmode .or. abs(rb(3)-1.d0)<tol) then
           ! do usual lms transformation
           call replmsf(zlms,zmat,lmax)
        else
           ! change spinor basis after lms transformation
           call replmsf_global(zlms,zmat,lmax,rb)
        end if
        call ytrafo(utr,utrp,lmax)
        call tripmt(utr,zlms(1:lmmaxp,1:lmmaxp),utrp,nl2,nl2,lmmaxp)
        call tripmt(utr,zlms(lmmaxp+1:kmymaxp,1:lmmaxp),
     >       utrp,nl2,nl2,lmmaxp)
        call tripmt(utr,zlms(1:lmmaxp,lmmaxp+1:kmymaxp),
     >       utrp,nl2,nl2,lmmaxp)
        call tripmt(utr,zlms(lmmaxp+1:kmymaxp,lmmaxp+1:kmymaxp),
     >              utrp,nl2,nl2,lmmaxp)
        do orb1=1,lmmaxp
          do orb2=1,lmmaxp
            tripletmat(orb1,orb2,1) = dimag(0.25*(
     >           zlms(orb1,orb2)-zlms(orb2,orb1)
     >     -zlms(lmmaxp+orb1,lmmaxp+orb2)+zlms(lmmaxp+orb2,lmmaxp+orb1)
     >           ))
            tripletmat(orb1,orb2,2) = dreal(-0.25*(
     >           zlms(orb1,orb2)-zlms(orb2,orb1)
     >     +zlms(lmmaxp+orb1,lmmaxp+orb2)-zlms(lmmaxp+orb2,lmmaxp+orb1)
     >           ))
            tripletmat(orb1,orb2,3) = dimag(0.25*(
     >           zlms(orb1,lmmaxp+orb2)-zlms(lmmaxp+orb2,orb1)
     >          +zlms(lmmaxp+orb1,orb2)-zlms(orb2,lmmaxp+orb1)
     >           ))
            end do
        end do
      end if
c
      return
      end
      subroutine tripletprint(tripletmat,nimp,ne,cear,tripletout)
      implicit none
      include '../param.h'
      character*34 tripletout
      integer*4 ie,iimp,k,nimp,ne
      complex*16 cear(me)
      real*8 er
      real*8 tripletmat(lmmaxp,lmmaxp,3,nimp,ne)
      open(31,file=tripletout,status='unknown')
      write(31,*) nimp, 'energy,  xy  yz  zz  xz  xx-yy antisymm. mat.
     > rowwise'
      do iimp=1,nimp
      write(31,*) 'IIMP',iimp
      do ie=1,ne
        er=dreal(cear(ie))
        do k=1,3
          write(31,'(e18.8,10e18.8)') er,tripletmat(5,6,k,iimp,ie),
     >                                   tripletmat(5,7,k,iimp,ie),
     >                                   tripletmat(5,8,k,iimp,ie),
     >                                   tripletmat(5,9,k,iimp,ie),
     >                                   tripletmat(6,7,k,iimp,ie),
     >                                   tripletmat(6,8,k,iimp,ie),
     >                                   tripletmat(6,9,k,iimp,ie),
     >                                   tripletmat(7,8,k,iimp,ie),
     >                                   tripletmat(7,9,k,iimp,ie),
     >                                   tripletmat(8,9,k,iimp,ie)
        end do
      end do
      end do
      close(31)
      return
      end
