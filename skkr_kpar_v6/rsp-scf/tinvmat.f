c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine tinvmat(tau1,tau2,lmax,iopt)
c
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
      integer kapdex(50),mydex(50),ldex(50),sdex(50)
      complex*16 tau1(kmymaxp,kmymaxp),tau2(kmymaxp,kmymaxp)
c
      data kapdex/-1,-1, 
     *             1, 1,-2,-2,-2,-2, 
     *             2, 2, 2, 2,-3,-3,-3,-3,-3,-3,
     *             3, 3, 3, 3, 3, 3,-4,-4,-4,-4,-4,-4,-4,-4,
     *             4, 4, 4, 4, 4, 4, 4, 4,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5/
      data  mydex/-1, 1,
     *            -1, 1,-3,-1, 1, 3,
     *            -3,-1, 1, 3,-5,-3,-1, 1, 3, 5,
     *            -5,-3,-1, 1, 3, 5,-7,-5,-3,-1, 1, 3, 5, 7,
     *            -7,-5,-3,-1, 1, 3, 5, 7,-9,-7,-5,-3,-1, 1, 3, 5, 7, 9/
      data   sdex/-1,-1, 
     *             1, 1,-1,-1,-1,-1,
     *             1, 1, 1, 1,-1,-1,-1,-1,-1,-1,         
     *             1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1,       
     *             1, 1, 1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1/
      data   ldex/ 0, 0, 
     *             1, 1, 1, 1, 1, 1,
     *             2, 2, 2, 2, 2, 2, 2, 2, 2, 2,         
     *             3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,       
     *             4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4/
c
c
      call czero(tau2,kmymaxp*kmymaxp)
      kmymax=2*(lmax+1)**2
c
      do i2=1,kmymax
      do j2=1,kmymax
        s1=sdex(i2)
        s2=sdex(j2)
        if(iopt.eq.1) then
          i1=i2-mydex(i2)
          j1=j2-mydex(j2)
          k=(mydex(i2)+mydex(j2))/2
          if(mod(k,2).eq.0) then
            sm=-1.0d0
          else
            sm=1.0d0
          end if
        else
          i1=j2-mydex(j2)
          j1=i2-mydex(i2)
          sm=1.0d0
        end if
        tau2(i2,j2)=sm*s1*s2*tau1(i1,j1)
      end do
      end do
c
      return
      end
