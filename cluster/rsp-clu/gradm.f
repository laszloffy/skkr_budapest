c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine gradm(rba,nimp,gradth,gradph,theta0,phi0,
     >                 dedmx,dedmy,dedmz)
c
c ==================================================================      
c The routines calculates the dE/dM derivatives from dE/dtheta and
c dE/dphi 
c    rba(3,mimp) : unit vector parallel to the magntization
c    nimp        : number of atoms in the cluster
c    gradth(mimp): dE/dtheta
c    gradph(mimp): dE/dphi
c    theta0, phi0: asimutal and polar angles for the orientation of 
c                  the reference system 
c    dedmx       : dE/dM_x 
c    dedmy       : dE/dM_y 
c    dedmz       : dE/dM_z 
c ==================================================================      
      implicit none
      include '../param.h'
      real*8 rba(3,mimp)
      real*8 gradth(mimp), gradph(mimp)
      real*8 theta0,phi0
      real*8 dedmx(mimp),dedmy(mimp),dedmz(mimp)
c

      real*8 GRADPHI,DMDZ

      integer*4 I,NIMP,J,K
c
      real*8 rn(3),eth(3),eph(3),s,en
c
      rn(1) = dsin(theta0)*dcos(phi0)      
      rn(2) = dsin(theta0)*dsin(phi0)      
      rn(3) = dcos(theta0)
c
      do i = 1,nimp
       en = 0.d0
       do j = 1,3
        en = en + rba(j,i)*rn(j)
       end do 
       s = 1.d0 - en*en
       do j = 1,3
        eth(j) = (rba(j,i)*en - rn(j))/dsqrt(s)
       end do
       eph(1) = (rn(2)*rba(3,i)-rn(3)*rba(2,i))/s
       eph(2) = (rn(3)*rba(1,i)-rn(1)*rba(3,i))/s
       eph(3) = (rn(1)*rba(2,i)-rn(2)*rba(1,i))/s
       write(6,'(3f10.4,4x,3f10.4)')(eth(j),j=1,3),
     >                              (eph(j)*dsqrt(s),j=1,3) 
       dedmx(i) = gradth(i)*eth(1) + gradph(i)*eph(1)
       dedmy(i) = gradth(i)*eth(2) + gradph(i)*eph(2)
       dedmz(i) = gradth(i)*eth(3) + gradph(i)*eph(3)
      end do      
      return
      end
      subroutine printjij(vij,gradth,gradph,dedmx,dedmy,dedmz,nimp)
      implicit none
      include '../param.h'
      integer*4 I,NIMP,J,K
c
      real*8 vij(4,mimp,mimp)
      real*8 gradth(mimp),gradph(mimp)
      real*8 dedmx(mimp), dedmy(mimp), dedmz(mimp)
c      
      write(6,*) ' Gradients'
      write(6,'(''   gradth'',5x,''   gradph'',7x,
     >      ''    dE/dM_x'',5x,''  dE/dM_y'', 5x,''  dE/dM_z'')')
      do i = 1,nimp
       write(6,'(2e14.6,2x,3e14.6)') gradth(i),gradph(i),
     >       dedmx(i),dedmy(i),dedmz(i)
      end do      
      write(6,*)
      write(6,*) 'Second derivatives'
      write(6,*) 
     > '       phi/phi      theta/theta   phi/theta     theta/phi' 
      do i = 1,nimp
       do j = 1,nimp
        write(6,'(2i3,4e14.6)') i,j,(vij(k,i,j),k=1,4)
       end do
      end do
      return
      end
