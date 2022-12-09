c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine jij(nintfc,kmymax,we0,tauij,
     >               dtmph,dtmth,d2tmph,d2tmth,d2tmthph, 
     >               vij,gradph,gradth)
c
c Input:
c  nintfc      : number of atoms in the cluster
c  we0         : weight of the energy integral
c  tauij()     : physical tau matrix of the cluster
c  dtmph           : dma90p(dt/dphi)dma90
c  dtmth           : dma90p(dt/dtheta)dma90
c  d2tmph           : dma90p(dt^2/dphi^2)dma90
c  d2tmth           : dma90p(dt^2/dtheta^2)dma90
c  d2tmthph         : dma90p(dt^2/dtheta dphi)dma90
c Input, Output:
c  diage()          : diagonal terms of the second derivatives
c   diage(1, )      :  phi,phi
c   diage(2, )      :  theta,theta
c   diage(3, )      :  phi,theta
c  vij()            : second derivatives 
c   vij(1, , )      :  phi,phi
c   vij(2, , )      :  theta,theta
c   vij(3, , )      :  phi,theta
c   vij(4, , )      :  theta,phi
c
      include '../param.h'
c     parameter (mimp=100)
      complex*16 we,we0
      complex*16 tauij(kmymaxp,kmymaxp,mimp,mimp)
c
      complex*16 dtmph(kmymaxp,kmymaxp,mimp)
      complex*16 dtmth(kmymaxp,kmymaxp,mimp)
      complex*16 d2tmph(kmymaxp,kmymaxp,mimp)
      complex*16 d2tmth(kmymaxp,kmymaxp,mimp)
      complex*16 d2tmthph(kmymaxp,kmymaxp,mimp)
c
      real*8 vij(4,mimp,mimp)
      real*8 tol,pi
c
      complex*16 tmwrk(kmymaxp,kmymaxp),wrk1(kmymaxp,kmymaxp)
      complex*16 wrk2(kmymaxp,kmymaxp,2)
      complex*16 wrk3(kmymaxp,kmymaxp,2)
      real*8     diag(4)
      complex*16 trdbmt
      integer kmymax
      real*8  gradph(mimp),gradth(mimp)
c
c Compute diagonal terms
c
! TEST
      tol=1.d-8
      pi = 4.d0*datan(1.d0)
      we = we0/pi
!     write(6,*) 'JIJ'
!     write(6,*) 'nintfc',nintfc
!     write(6,*) 'kmymax',kmymax
!     write(6,*) 'we',we
!     call outmat1(tauij(1,1,1,1),kmymax,kmymax,kmymaxp,tol,6)
!     call outmat1(tauij(1,1,1,2),kmymax,kmymax,kmymaxp,tol,6)
!     call outmat1(tauij(1,1,2,2),kmymax,kmymax,kmymaxp,tol,6)
!     call outmat1(tauij(1,1,2,1),kmymax,kmymax,kmymaxp,tol,6)
!
!     write(6,*) 'dtmph 1'
!     call outmat1(dtmph(1,1,1),kmymax,kmymax,kmymaxp,tol,6)
!     write(6,*) 'dtmph 2'
!     call outmat1(dtmph(1,1,2),kmymax,kmymax,kmymaxp,tol,6)
!     write(6,*) 'dtmth 1'
!     call outmat1(dtmth(1,1,1),kmymax,kmymax,kmymaxp,tol,6)
!     write(6,*) 'dtmth 2'
!     call outmat1(dtmth(1,1,2),kmymax,kmymax,kmymaxp,tol,6)
!     write(6,*) 'd2tmph 1'
!     call outmat1(d2tmph(1,1,1),kmymax,kmymax,kmymaxp,tol,6)
!     write(6,*) 'd2tmth 1'
!     call outmat1(d2tmth(1,1,1),kmymax,kmymax,kmymaxp,tol,6)
!     write(6,*) 'd2tmthph 1'
!     call outmat1(d2tmthph(1,1,1),kmymax,kmymax,kmymaxp,tol,6)
!     call flush(6)
! TEST
       do li = 1,nintfc
        ip = li
         do i = 1,4
          diag(i) = 0.d0
         end do
        do i=1,kmymax
          do j=1,kmymax
           diag(1) = diag(1) +
     >       dimag(we*tauij(i,j,li,li)*d2tmph(j,i,li))
           diag(2) = diag(2) +
     >       dimag(we*tauij(i,j,li,li)*d2tmth(j,i,li))
           diag(3) = diag(3) +
     >       dimag(we*tauij(i,j,li,li)*d2tmthph(j,i,li))
c Gradients
           gradph(li) = gradph(li) +
     >       dimag(we*tauij(i,j,li,li)*dtmph(j,i,li))
           gradth(li) = gradth(li) +
     >       dimag(we*tauij(i,j,li,li)*dtmth(j,i,li))
          end do
         end do
c TEST
c     write(6,*) 'JIJ 2h'
c     call flush(6)
c     call outmat1(tauij(1,1,li,li),kmymax,kmymax,kmymaxp,tol,6)
c     write(6,*) 'JIJ 2h1'
c     call flush(6)
c     call outmat1(dtmph(1,1,li),kmymax,kmymax,kmymaxp,tol,6)
c     write(6,*) 'JIJ 2h2'
c     call flush(6)
c     call outmat1(wrk2(1,1,1),kmymax,kmymax,kmymaxp,tol,6)
c     write(6,*) 'JIJ 2h3'
c     call flush(6)
c TEST
!        call doubmt1(tauij(1,1,li,li),dtmph(1,1,li),wrk2(1,1,1),
!    >                   kmymax,kmymaxp)
c TEST
c     write(6,*) 'JIJ 2i'
c     call flush(6)
c TEST
!        call doubmt1(tauij(1,1,li,li),dtmth(1,1,li),wrk2(1,1,2),
!    >                   kmymax,kmymaxp)
c TEST
c     write(6,*) 'JIJ 2j'
c     call flush(6)
c TEST
c TEST
c     write(6,*) 'JIJ 2k'
c     call flush(6)
c TEST
!        do i = 1,kmymax
c TEST
c     write(6,*) 'JIJ 2l'
c     call flush(6)
c TEST
!         do j = 1,kmymax
c TEST
c     write(6,*) 'JIJ 2m'
c     call flush(6)
c TEST
!          diag(1) = diag(1) +
!    >        dimag(we*wrk2(i,j,1)*wrk2(j,i,1))
c TEST
c     write(6,*) 'JIJ 2n'
c     call flush(6)
c TEST
!          diag(2) = diag(2) +
!    >        dimag(we*wrk2(i,j,2)*wrk2(j,i,2))
c TEST
c     write(6,*) 'JIJ 2o'
c     call flush(6)
c TEST
!          diag(3) = diag(3) +
!    >        dimag(we*wrk2(i,j,1)*wrk2(j,i,2))
c TEST
c     write(6,*) 'JIJ 2p'
c     call flush(6)
c TEST
!         end do
c TEST
c     write(6,*) 'JIJ 2q'
c     call flush(6)
c TEST
!        end do
!        diage(1,li) = diage(1,li) + diag(1)
!        diage(2,li) = diage(2,li) + diag(2)
!        diage(3,li) = diage(3,li) + diag(3)
!       end do
c TEST
c     write(6,*) 'JIJ 3'
c     call flush(6)
c TEST
c 
c Off-diagonal terms
c 
c       do ip = 1,nintfc
         do iq = 1,nintfc
          call doubmt1(tauij(1,1,ip,iq),dtmph(1,1,iq),wrk2(1,1,1),
     >                kmymax,kmymaxp)
          call doubmt1(tauij(1,1,ip,iq),dtmth(1,1,iq),wrk2(1,1,2),
     >                kmymax,kmymaxp)
          call doubmt1(tauij(1,1,iq,ip),dtmph(1,1,ip),wrk3(1,1,1),
     >                kmymax,kmymaxp)
          call doubmt1(tauij(1,1,iq,ip),dtmth(1,1,ip),wrk3(1,1,2),
     >                kmymax,kmymaxp)
          vij(1,ip,iq) = vij(1,ip,iq) - dimag(we*trdbmt(wrk2(1,1,1),
     >                wrk3(1,1,1), kmymax,kmymaxp))
          vij(2,ip,iq) = vij(2,ip,iq) - dimag(we*trdbmt(wrk2(1,1,2),
     >                wrk3(1,1,2), kmymax,kmymaxp))
          vij(3,ip,iq) = vij(3,ip,iq) - dimag(we*trdbmt(wrk2(1,1,1),
     >                wrk3(1,1,2), kmymax,kmymaxp))
          vij(4,ip,iq) = vij(4,ip,iq) - dimag(we*trdbmt(wrk2(1,1,2),
     >                wrk3(1,1,1), kmymax,kmymaxp))
          if(ip.eq.iq) then
           vij(1,ip,iq) = vij(1,ip,iq) + diag(1)
           vij(2,ip,iq) = vij(2,ip,iq) + diag(2)
           vij(3,ip,iq) = vij(3,ip,iq) + diag(3)
           vij(4,ip,iq) = vij(4,ip,iq) + diag(3)
          endif   
         end do          
!     write(6,*) 'diag1=',diag(1)
!     write(6,*) 'diag2=',diag(2)
!     write(6,*) 'diag3=',diag(3)
!     write(6,*) 'vij1=',vij(1,1,1)
!     write(6,*) 'vij2=',vij(2,1,1)
!     write(6,*) 'vij3=',vij(3,1,1)
!     write(6,*) 'vij4=',vij(4,1,1)
        end do          
c TEST
c     write(6,*) 'JIJ 4'
c     call flush(6)
c TEST
        return
        end
