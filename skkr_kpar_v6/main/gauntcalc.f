c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine gauntcalc(lmax,gau)
      implicit none
!i/o
      integer lmax
      double precision gau(*),pi,blm
      parameter (pi=3.1415926535897932384626d0) 
!local       
      integer o,i,j,k,l,l1,l2,m,m1,m2
!=================      
      o=0
      do l1=0,lmax
         do m1=-l1,l1
            do l2=0,lmax
               do m2=-l2,l2
                  do l=0,lmax*2
                     do m=-l,l
                        o=o+1
                        gau(o)=blm(l2,m2,l1,-m1,l,m)*4.0*pi* 
     &                       (-1.d0)**((-l2-l+l1+2*m1)/2)
                     enddo ! m
                  enddo ! l
               enddo ! m2
            enddo ! l2
         enddo ! m1
      enddo ! l1
      return
      end
