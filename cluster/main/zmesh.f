c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine zmesh(imesh,npanel,nepanel,ne1,ne2,ne3,ebottom,etop,
     >                 epsilon,z,ww)
c===================== 
c
c set up contour and grid on complex energy plane 
c
c input:  imesh - if 0 paralell to real axis
c                    1 semicircle
c                    2 Matsubara poles for finite T
c                    3 Matsubara poles for finite T + semi-circle
c         npanel - number of energy panels for imesh=1 or 2
c         ne1,ne2,ne3 - number of energy points on different contour pathes
c                    for imesh=3
c         ebottom - lower boundary on real energy axis
c         etop  - upper boundary on real energy axis (Fermi level)
c         epsilon - parameter for log. gaussian mesh on semicircle
c output: z  - array containing complex energies
c         ww - weights
c
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
      dimension nepanel(5),ebottom(5),etop(5),epsilon(5)
      dimension z(me),ww(me),x(me),w(me),z1(me),w1(me)
      complex*16 z,ww,z1,w1,zbottom,zedge,ztop,z0,sqrtm1
c
      common/test/itest
      data sqrtm1/(0.d0,1.d0)/,eps0/1.d-5/
      data ekelvin/6.333382d-06/,temp0/0.1d0/,akt/8.d0/
c
      if(itest.ge.2) write(6,'(/'' < ZMESH'')')
      pi=4.d0*datan(1.d0)
c
      if(imesh.eq.0) then
        i=0
        do n=1,npanel
         if(nepanel(n).eq.1) then
            i=i+1
            z(i)=dcmplx(ebottom(n),epsilon(n))
            ww(i)=dcmplx(1.0d0,0.d0)
         else 
            der=(etop(n)-ebottom(n))/(nepanel(n)-1)
            do ii=1,nepanel(n)
              i=i+1
              z(i)=dcmplx(ebottom(n)+(ii-1)*der,epsilon(n))
              ww(i)=dcmplx(der,0.d0)
            end do
            ww(i-nepanel(n)+1)=dcmplx(0.5d0*der,0.d0)
            ww(i)=dcmplx(0.5d0*der,0.d0)
         end if
        end do
      end if
c
      if(imesh.eq.1) then
       i=0
       do n=1,npanel
        if(epsilon(n).gt.eps0) eps=1.d0/epsilon(n)
        r=(etop(n)-ebottom(n))/2.d0
        phi1=pi
        phi2=0.d0
        if(epsilon(n).gt.eps0) then
          y1=-dlog(1.d0+phi1/eps)
          y2=-dlog(1.d0+phi2/eps)
        end if
        z0=dcmplx(ebottom(n)+r,0.d0)
        call legzero(nepanel(n),x,w)
        do ii=1,nepanel(n)
          i=i+1
          if(epsilon(n).gt.eps0) then
            y=0.5d0*(y2-y1)*x(ii)+0.5d0*(y2+y1)
            phi=eps*(dexp(-y)-1.d0)
            z(i)=r*cdexp(sqrtm1*phi)+z0
            ww(i)=-0.5d0*(y2-y1)*eps*dexp(-y)*sqrtm1*(z(i)-z0)*w(ii)
          else
            phi=0.5d0*(phi2-phi1)*x(ii)+0.5d0*(phi2+phi1)
            z(i)=r*cdexp(sqrtm1*phi)+z0
            ww(i)=0.5d0*(phi2-phi1)*sqrtm1*(z(i)-z0)*w(ii)
          end if
          if(itest.ge.2) write(6,'(4d20.8)') z(i),ww(i)
        end do
       end do
      end if
c
      if(imesh.ge.2) then
c  Matsubara poles
      if(npanel.ne.1) stop 'Only one energy panel accessible!'
      ne=nepanel(1)
      efermi=etop(1)
      temp=epsilon(1)
      if(ne.le.0) stop ' ne.le.0: How many Matsubara-poles ???'
      if(temp.lt.temp0) stop ' Temperature is too small !!!'
      ekbt=ekelvin*temp
      delta=2.0d0*pi*ekbt*ne
      do j=1,ne
        j3=ne+1-j
        z(j3)=dcmplx(efermi,pi*(2*j-1)*ekbt)
        ww(j3)=-2.d0*pi*sqrtm1*ekbt
      end do
      end if
c
      if(imesh.eq.3) then
c  Matsubara poles and contour integration
        ne4=ne
        ne=ne1+ne2+ne3+ne4
        do i=1,ne4
          z1(i)=z(i)
          w1(i)=ww(i)
        end do
        do i=1,ne4
          z(ne1+ne2+ne3+i)=z1(i)
          ww(ne1+ne2+ne3+i)=w1(i)
        end do
c   semi circle contour
        call legzero(ne1,x,w)
        do i = 1,ne1
         y = pi/4.d0*(3.d0 - x(i))
         z(i) = delta*(1.d0 + cdexp(sqrtm1*y)) + ebottom(1)
         z0 = 1.d0/(1.d0 + cdexp((z(i) - etop(1))/ekbt))
         ww(i) = -sqrtm1*pi/4.d0*delta*cdexp(sqrtm1*y)*z0*w(i)
        end do
        ie = ne1
c   line 2piNkT above the real axis
        a = (etop(1) - ebottom(1) - delta - akt*ekbt)/2.d0
        b = (etop(1) + ebottom(1) + delta - akt*ekbt)/2.d0  
        call legzero(ne2,x,w)
        do i = 1,ne2
         ie = ie + 1
         z(ie) = dcmplx(a*x(i) + b, delta)
         y = 1.d0/(1.d0 + dexp((a*x(i)+b-etop(1))/ekbt))
         ww(ie) = dcmplx(w(i)*a*y,0.d0)
        end do
c   line 2piNkT above the real axis around eF
        call legzero(ne3,x,w)
        do i = 1,ne3
         ie = ie + 1
         z(ie)  = dcmplx(akt*ekbt*x(i)+etop(1),delta)
         y = 1.d0/(1.d0 + dexp(akt*x(i)))
         ww(ie) = dcmplx(akt*ekbt*y*w(i),0.d0)
        end do
      end if 
c
      return
      end
