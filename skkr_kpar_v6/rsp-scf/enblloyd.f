c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine enblloyd(enbll,we,tau,ptinv,tinv,alphakkr,lmax,lliter)
c=======================
c
c Band energy of an impurity via Lloyd formula treated iteratively
c
c input:  we           ... weighting factor for energy integration 
c         tau          ... unscreened tau matrix 
c         ptinv        ... inverse of unscreened impurity t-matrix 
c         tinv         ... inverse of screened host t-matrix
c         alphakkr     ... screening parameters 
c         lmax         ... actual max. value of angular momentum index l
c output: enbll        ... integrand for the band energy
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      complex*16 det,we
      complex*16 alphakkr(0:lmaxp)
      complex*16 tau(kmymaxp,kmymaxp)
      complex*16 ptinv(kmymaxp,kmymaxp)
      complex*16 tinv(kmymaxp,kmymaxp)
      complex*16 tinvc(kmymaxp,kmymaxp)
      complex*16 cmat(kmymaxp,kmymaxp)
      complex*16 xmat(kmymaxp,kmymaxp)
c
      data tol/1.0d-10/,tiny/1.0d-3/
      integer ldex(50)
      data ldex/0,0,
     *          1,1,1,1,1,1,
     *          2,2,2,2,2,2,2,2,2,2,
     *          3,3,3,3,3,3,3,3,3,3,3,3,3,3,
     *          4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4/
c
      pi=4.0d0*datan(1.0d0)
      lmmax=(lmax+1)**2
      kmymax=2*lmmax
c
c unscreening of inverse t-matrix
      call repl(tinvc,tinv,kmymax,kmymaxp)
      call gjinv(tinvc,kmymax,kmymaxp,det)
      do kmy=1,kmymax
        tinvc(kmy,kmy)=tinvc(kmy,kmy)+alphakkr(ldex(kmy))
      end do
      call gjinv(tinvc,kmymax,kmymaxp,det)
c
c     write(6,'('' phys. host m-matrix'')')
c     call outmat1(tinvc,kmymax,kmymax,kmymaxp,tol,6)
c     write(6,'('' phys. imp. m-matrix'')')
c     call outmat1(ptinv,kmymax,kmymax,kmymaxp,tol,6)
c
      call repl(cmat,tinvc,kmymax,kmymaxp)
      call submat(cmat,ptinv,kmymax,kmymaxp)
c     write(6,'('' difference m-matrix'')')
c     call outmat1(cmat,kmymax,kmymax,kmymaxp,tol,6)
      call doubmt(cmat,tau,kmymax,kmymaxp)
c     write(6,'('' m*tau'')')
c     call outmat1(cmat,kmymax,kmymax,kmymaxp,tol,6)
c
      call czero(xmat,kmymaxp*kmymaxp)
      do kmy=1,kmymax
        xmat(kmy,kmy)=(1.0d0,0.0d0)
      end do
c
      enbll=0.0d0
      do l=1,lliter 
        call doubmt(xmat,cmat,kmymax,kmymaxp)
        do kmy=1,kmymax
          enbll=enbll+dimag(we*xmat(kmy,kmy))/(pi*float(l))
        end do
        if(l.ge.2.and.dabs(enbll0-enbll).lt.tiny) goto 20
        write(6,'(i3,d15.6)') l,enbll
        enbll0=enbll
      end do
      write(6,'(''Iteration of Lloyd formula not converged'')')
  20  continue
c
      return
      end
