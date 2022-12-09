c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine bulkmad2(nbulk,lmax,sigma,vol,cbulk,cbulk0,conmad)
c========================
c
c Calculate Madelung constants for a bulk with complex lattice
c
c input:  nbulk - number of sublattices
c         lmax - l-cutoff
c         sigma - Ewald parameter
c         vol - volume of unit cell
c         cbulk - non-equivalent positions
c         cbulk0 - 'third' vector from Bravais matrix
c output: conmad - Madelung constants
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
      dimension cbulk(3,mbulk),cbulk0(3),cvec(3)
c
      complex*16 conmad(lmsup,lmsup,mbulk,mbulk)
c
      common/test/itest
c
      data tol/1.0d-8/,kmax/100/
c
      pi=4.d0*datan(1.d0)
      fac=4.d0*pi/vol
c
c Calculate Madelung potential for layers in a bulk unit 
c
      do i=1,nbulk
      do j=1,nbulk
c
c         write(6,'('' Atoms:'',2i3)') i,j
c 1. G.ne.0 contributions
c
c
          sum1=0.d0
          sum2=0.d0
          sum3=0.d0
          do k=0,kmax
c Loop over bulk units upwards
c
c Separation vector between layers (i,0) and (j,k)
c
            cvec(1)=cbulk(1,i)-cbulk(1,j)-k*cbulk0(1) 
            cvec(2)=cbulk(2,i)-cbulk(2,j)-k*cbulk0(2) 
            cvec(3)=cbulk(3,i)-cbulk(3,j)-k*cbulk0(3)
c
c           ------------------------------
            call phi(cvec,sigma,fi,psi,xi)
c           ------------------------------
c
            sum1=sum1+fi
            sum2=sum2+psi
            sum3=sum3+xi  
c           write(6,'(i3,2x,2d15.7,2x,2d15.7,5x,2d15.7,2x,2d15.7,
c    >      5x,2d15.7,2x,2d15.7)') k,fi,sum1,psi,sum2,xi,sum3
            if(dabs(fi).lt.tol) goto 1 
c
          end do   
          write(6,'(/'' WARNING!!! Madelungb: ''/'' layer '',
     >    ''sum did not converge up to bulk unit no.'',i4/)') k
    1     continue
c
          do k=-1,-kmax,-1
c Loop over bulk units downwards
c
c Separation vector between layers (i,0) and (j,k)
c
            cvec(1)=cbulk(1,i)-cbulk(1,j)-k*cbulk0(1) 
            cvec(2)=cbulk(2,i)-cbulk(2,j)-k*cbulk0(2) 
            cvec(3)=cbulk(3,i)-cbulk(3,j)-k*cbulk0(3)
c
c           ------------------------------
            call phi(cvec,sigma,fi,psi,xi)
c           ------------------------------
c
            sum1=sum1+fi 
            sum2=sum2+psi
            sum3=sum3+xi 
c           write(6,'(i3,2x,2d15.7,2x,2d15.7,5x,2d15.7,2x,2d15.7,
c    >      5x,2d15.7,2x,2d15.7)') k,fi,sum1,psi,sum2,xi,sum3
            if(dabs(fi).lt.tol) goto 2 
c
          end do
          write(6,'(/'' WARNING!!! Madelungb: ''/'' layer '',
     >    ''sum did not converge up to bulk unit no.'',i4/)') k
   2      continue
c
c 2. G=0 contribution
c        
          cperp=cbulk(3,i)-cbulk(3,j)
c         write(6,*) cperp,cbulk0(3)
          sum1=sum1+fac*(cperp*cperp/cbulk0(3)
     >                   -dabs(cperp)+cbulk0(3)/6.0d0)
          if(cperp.lt.-tol) sum2=sum2-fac*(2.0d0*cperp/cbulk0(3)+1.0d0)
          if(cperp.gt.tol) sum2=sum2-fac*(2.0d0*cperp/cbulk0(3)-1.0d0)
          sum3=-sum3*dsqrt(5.0d0/4.0d0)+fac*dsqrt(5.d0)/(3.d0*cbulk0(3))
c
          conmad(1,1,i,j)=sum1
          conmad(1,3,i,j)=sum2*dsqrt(3.0d0)
          conmad(3,1,i,j)=sum2*dsqrt(3.0d0)
c         conmad(1,7,i,j)=sum3
c         conmad(7,1,i,j)=sum3
c
      end do
      end do
c
      return
      end
