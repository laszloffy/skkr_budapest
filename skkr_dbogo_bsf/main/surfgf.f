c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine surfgfd(ml,m0,mmr,x,kdim,itermax,errmax,ichck,sgfconv)
c====================================================================
c
c solve surface green's function: f(x)=ml*(m0-x)**(-1)*mmr
c method: decimation technique
c
c input:  ml,m0,mr - complex square matrices
c         kdim     - actual dimension of matrices
c output: x        - result, matrix of same type as before
c        
c
      implicit real*8(a-h,o-z)
c
      include '../param.h'
      parameter(mdim=mdimdbogo)
c
      logical sgfconv
      complex*16 ml,m0,mmr,x
      complex*16 alfa,beta,eps,y1,y2,y3,y4,detl
      dimension ml(mdim,mdim),m0(mdim,mdim),mmr(mdim,mdim),
     &   x(mdim,mdim),alfa(mdim,mdim),beta(mdim,mdim),eps(mdim,mdim),
     &   y1(mdim,mdim),y2(mdim,mdim),y3(mdim,mdim),y4(mdim,mdim)
      data errmin/1.0d+10/
c
      if(ichck.gt.0) write(6,'('' Surface tau-matrix'')')
      sgfconv=.true.
c
      call czero(x,mdim*mdim)
      call czero(alfa,mdim*mdim)
      call czero(beta,mdim*mdim)
      call czero(eps,mdim*mdim)
      call czero(y1,mdim*mdim)
      call czero(y2,mdim*mdim)
      call czero(y3,mdim*mdim)
      call czero(y4,mdim*mdim)
c
      do i=1,kdim
      do j=1,kdim
        eps(i,j)=m0(i,j)
        alfa(i,j)=-ml(i,j)
        beta(i,j)=-mmr(i,j)
        x(i,j)=m0(i,j)
      end do
      end do
c
      iter=1
    1 continue
      do i=1,kdim
      do j=1,kdim
        y1(i,j)=eps(i,j)
      end do
      end do
      call gjinv(y1,kdim,mdim,detl)
      do i=1,kdim
      do j=1,kdim
        y2(i,j)=y1(i,j)
        y3(i,j)=y1(i,j)
        y4(i,j)=y1(i,j)
      end do
      end do
      call tripmt(alfa,y1,alfa,kdim,kdim,mdim)
      call tripmt(beta,y2,beta,kdim,kdim,mdim)
      call tripmt(alfa,y3,beta,kdim,kdim,mdim)
      call tripmt(beta,y4,alfa,kdim,kdim,mdim)
c
      sum=0.d0
      do i=1,kdim
      do j=1,kdim
c
        alfa(i,j)=y1(i,j)
        beta(i,j)=y2(i,j)
        eps(i,j)=eps(i,j)-y3(i,j)-y4(i,j)
        x(i,j)=x(i,j)-y3(i,j)
c
        xre=dreal(alfa(i,j))
        xim=dimag(alfa(i,j))
        sum=sum+xre*xre+xim*xim
      end do
      end do
      err=dsqrt(sum)
c
      if(ichck.gt.0) then
        write(6,'('' iter='',i4,''  error='',d12.7)') iter,err
        call flush(6)
      end if
c
      if(err.lt.errmax.or.iter.gt.itermax) goto 2
      if(err.gt.errmin) then
        write(6,'('' Surfgf: enormous error '',d12.7,
     &  '' on iter'',i3)') err,iter
        sgfconv=.false.
        goto 2
      endif
      iter=iter+1
      goto 1
c
    2 continue
c
      if(iter.gt.itermax) then 
        write(6,'('' itermax too small !  iter='',i3)') iter
        stop
      end if
c
      call gjinv(x,kdim,mdim,detl)
c      
      return
      end
