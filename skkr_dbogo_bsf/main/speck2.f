c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine speck2(nset,xk,yk,wk,nk,ilat,alat,beta,delta)
c ======================
c
c k points and (unnormalized) weights for the minimum IBZ of each 2D 
c crystal system (according to Cunningham, PRB 10 (1974), fig.1), 
c using approximately an equidistant point sampling.
c Caution: k values on output are multiplied by pi/a !
c
c input:  ilat - type of lattice
c         alat - lattice parameter
c         beta - lattice asymmetry ratio 
c         nset - determines nr. of k points
c output: xk,yk - array for k points
c         wk   - array for weights (without normalization)
c         nk   - number of points
c
      implicit real*8 (a-h,o-z)
c
      dimension xk(*),yk(*),wk(*)
      dimension v(2),dk1(2),dk2(2)
c
      data t/1.732050807568877d0/
      data pi/3.1415926535897932d0/
c
      if(ilat.le.0.or.ilat.gt.5) 
     > stop ' *** non-existing 2D lattice type ***'
c
c    --------
      n=nset
c    --------
c
c ********************
c oblique lattice (C2)
c ********************
c
      if(ilat.eq.1) then
c
        if(delta.lt.0.d0) then
          del=-delta
        else
          del=delta
        end if
        m=int(n/sqrt(beta*beta+del*del))+1
        dk1(1)=2.0d0/n
        dk1(2)=-2.0d0*del/(n*beta)
        dk2(1)=0.0d0
        dk2(2)=2.0d0/(m*beta)
        bod=beta/del
        bo1md=beta/(1.0d0-del)
        tiny=1.0d-6*dk2(2)
c
        xmax=1.0d0+(1.0d0/bod+1.0d0/bo1md)/(bod+bo1md)
        nk=0
        do 5 i=1,n
        do 5 j=-m+1,m
c
          x=(i-0.5d0)*dk1(1)+(j-0.5d0)*dk2(1)
          y=(i-0.5d0)*dk1(2)+(j-0.5d0)*dk2(2)
          if(x.gt.xmax) goto 5
c
          ymax1=1.0d0/beta
          ymax2=-bo1md*(x-1.0d0)+1.0d0/bo1md
          ymax=dmin1(ymax1,ymax2)
          ymin1=-1.0d0/beta
          ymin2=bod*(x-1.0d0)-1.0d0/bod
          ymin=dmax1(ymin1,ymin2)
          dymax1=ymax1-y
          dymax2=ymax2-y
          dymin1=y-ymin1
          dymin2=y-ymin2
c
          if(y.le.ymax.and.y.ge.ymin) then
            nk=nk+1
            xk(nk)=x
            yk(nk)=y
            if(dymax1.lt.tiny) then
              if(dymax2.lt.tiny) then
                wk(nk)=0.25d0
              else
                wk(nk)=0.5d0
              end if
            else if(dymax2.lt.tiny) then
              if(dymin1.lt.tiny) then
                wk(nk)=0.25d0
              else
                wk(nk)=0.5d0
              end if
            else if(dymin1.lt.tiny) then
              if(dymin2.lt.tiny) then
                wk(nk)=0.25d0
              else
                wk(nk)=0.5d0
              end if
            else if(dymin2.lt.tiny) then
              wk(nk)=0.5d0
            else
              wk(nk)=1.0d0
            end if
          else
            goto 5
          end if
c
    5   continue
c
        if(delta.lt.0.d0) then
          do ik=1,nk
            yk(ik)=-yk(ik)
          end do
        end if
c
      end if
c
c **********************************
c centered rectangular lattice (C2v)
c **********************************
c nset values of x
c
      if(ilat.eq.2) then
c
        dkx = 2.d0/dfloat(n+1)
        dky = 2.d0/dfloat(n+2)/beta
        tiny = dkx*1.d-6
        nk=0
        x=dkx/2.d0
        dowhile(x.le.2.d0)
          difx=dabs(x-2.d0)
          y=dky/2.d0
          ymax = -x*beta+(beta+1.d0/beta)
          dowhile(y.le.ymax)
            dify=dabs(y-ymax)
            nk=nk+1
            xk(nk)=x
            yk(nk)=y
            if(difx.lt.tiny) then
              wk(nk)=0.5d0
              if(dify.lt.tiny) wk(nk)=0.25d0
            else
              wk(nk)=1.d0
              if(dify.lt.tiny) wk(nk)=0.5d0
            endif
            y=y+dky
          enddo
          x=x+dkx
        enddo
c
      end if
c
c ***********************************
c primitive rectangular lattice (C2v)
c ***********************************
c nset values of x 
c
      if(ilat.eq.3) then
c  
        dkx = 1.d0/dfloat(n)
        dky = (1.d0/beta)/dfloat(n)
c       tinx = dkx*1.d-6
c       tiny = dky*1.d-6
        nk=0
        x=dkx/2.0d0
        dowhile(x.le.1.d0)
c         difx=dabs(x-1.d0)
          y=dky/2.0d0
          dowhile(y.le.1.d0/beta)
c           dify=dabs(y-1.d0/beta)
            nk=nk+1
            xk(nk)=x
            yk(nk)=y
            wk(nk)=1.d0
c           if(dify.lt.tiny) wk(nk)=0.5d0
            y=y+dky
          enddo
          x=x+dkx
        enddo
c
      end if
c
c ********************
c square lattice (C4v)
c ********************
c nset values of x 
c
      if(ilat.eq.4) then
c
        w=1.d0
        ik=0
        do i=1,n
          x=2*i-1
          x=0.5d0*x/n
          do j=1,i
            ik=ik+1
            y=2*j-1
            y=0.5d0*y/n
            xk(ik)=x
            yk(ik)=y
            wk(ik)=w
            if(i.eq.j) wk(ik)=0.5d0*w
          end do
        end do
        nk=ik
c
      end if
c
c ************************
c hexagonal lattice: (C6v)
c ************************
c nset values of y 
c
      if(ilat.eq.5) then
c
        ik=0
        dxk=2.0d0/dfloat(n)
        dyk=dxk/dsqrt(3.d0)
        do nu=1,2
          anu3=dfloat(nu)/3.d0
          do i=1,n
          do j=i,n
            iy=i+2*j+nu-3
            if(iy.gt.n) go to 10
            ik=ik+1
            ajm=1.d0
            if(i.eq.j.or.iy.eq.n) ajm=2.d0
            if(i.eq.j.and.iy.eq.n) ajm=6.d0
            wk(ik)=1.0/ajm
            xk(ik)=(dfloat(i-1)+anu3)*dxk
            yk(ik)=iy*dyk
   10       continue
          end do
          end do
        end do
        nk=ik
c
c
c       Rotate by -pi/3 to place one edge of IBZ along X
c
        phi=-pi/3.0d0
        do ik=1,nk
          v(1)=xk(ik)
          v(2)=yk(ik)
          call vecrot(v,phi)
          xk(ik)=v(1)
          yk(ik)=v(2)
        end do
c
      end if
c
c ****************************************************
c
      const=pi/alat
      do ik=1,nk
        xk(ik)=const*xk(ik)
        yk(ik)=const*yk(ik)
      end do
c
      return
      end
