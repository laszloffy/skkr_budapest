c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine moment(lmax,lms,coefg,coeff,sig,ggrmat,fgrmat,zmom,rb)
c======================
c
c input:  lmax  - maximum of angular momentum index
c         lms - if .true. output in (l,m,s) representation
c         coefg - angular momentum matrix for the large components
c         coeff - angular momentum matrix for the small components
c         sig - sign of the contribution from the small components
c         ggrmat, fgrmat - matrices of the radial integrals of the Green function
c         rb - orientation of exchange field, z in local mode, general in global mode
c output: zmom - (kappa,my)-like or (l,m,s)-like resolution of the moment 
c
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
      real*8 pi,fac
      complex*16 zgmat(kmymaxp,kmymaxp)
      complex*16 zfmat(kmymaxp,kmymaxp)
      complex*16 zmat(kmymaxp,kmymaxp)
      complex*16 zlms(kmymaxp,kmymaxp)
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
c          call replms(zmom,zmat,lmax)
           call replmsf(zlms,zmat,lmax)
        else
           ! change spinor basis after lms transformation
c          call replms_global(zmom,zmat,lmax,rb)
           call replmsf_global(zlms,zmat,lmax,rb)
        end if
        call ytrafo(utr,utrp,lmax)
        call tripmt(utr,zlms(1:lmmaxp,1:lmmaxp),utrp,nl2,nl2,lmmaxp)
        call tripmt(utr,zlms(lmmaxp+1:kmymaxp,lmmaxp+1:kmymaxp),
     >              utrp,nl2,nl2,lmmaxp)
        do i=1,kmymax
          zmom(i)=zlms(i,i)
        end do
      else
        do i=1,kmymax
          zmom(i)=zmat(i,i)
        end do
      end if
c
      return
      end
