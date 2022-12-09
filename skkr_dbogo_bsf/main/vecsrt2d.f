c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
C***vecsrt.spg  processed by SPAG 5.11R  at 12:38 on 17 Feb 2000
      subroutine vecsrt2d(maxx,maxind,r,rho,dist,ind,nsh,maxsh)
c --------------------------------------------------------------------**
c program to index vectors (r+rho) according to length
c and find vectors of equal length
c
c input:
c   maxx   nr. of lattice vectors
c   maxind size of index array
c   r      set of lattice vectors
c   rho    non-primitive vector in the unit cell
c   dist   maximum length of vector to be indexed
c
c output:
c   ind    index array: r(ind(n+1))+rho >= r(ind(n))+rho
c   nsh    number of vectors of equal length
c   maxsh  number of shells of vectors of equal length
c          (=number of elements in nsh)
c ------------------------- dummy arguments --------------------------**
c                                                                     **
c  input  - maxx                                                      **
c  input  - maxind                                                    **
c  input  - r                                                         **
c  input  - rho                                                       **
c  input  - dist                                                      **
c  output - ind                                                       **
c  output - nsh                                                       **
c  output - maxsh                                                     **
c                                                                     **
c ------------------------- common variables -------------------------**
c                                                                     **
c  modifies    ** nothing **                                          **
c  uses value  ** nothing **                                          **
c                                                                     **
c ----------------------- external subprograms -----------------------**
c                                                                     **
c  calls       indekx                                                 **
c  called by                                                          **
c                                                                     **
c --------------------------------------------------------------------**
c
      implicit real*8 (a-h,o-z)
c     implicit none
C*** Start of declarations inserted by SPAG
c     real*8 dist,dist2,r,r1,r2,rho,rmod,tol
c     integer i,ind,indexx,k,maxx,maxind,maxsh,nsh
C*** End of declarations inserted by SPAG
c
      include '../param.h'
c
      parameter (tol=1.D-06)
      dimension r(2,maxx),ind(maxx),nsh(maxind),rho(2)
      dimension rmod(mdr),indexx(mdr)
c
      if (mdr.lt.maxx) stop 'mdr too small'
      dist2=dist*dist
c
      do i=1,maxx
        ind(i)=0
        rmod(i)=(r(1,i)+rho(1))**2+(r(2,i)+rho(2))**2
      enddo
c
      call indekx(maxx,rmod,indexx)
c
      do i=1,maxind
        nsh(i)=0
      enddo
c
      r1=-1.D0
      k=0
      do i=1,maxx
        r2=rmod(indexx(i))
        if (r2.gt.dist2) goto 100
        if (i.eq.maxx) then
          if (dabs(r2-r1).le.tol) k=k-1
          goto 100
        endif
        ind(i)=indexx(i)
        if (dabs(r2-r1).gt.tol) k=k+1
        if (k.gt.maxind) then
          write(6,'(a55,i6,a3)') 
     &      "ERROR <vecsrt>: increase either mr or mk in param.h (k=",k,
     &      ") !"
          stop
        endif
        nsh(k)=nsh(k)+1
        r1=r2
      enddo
c
  100 continue
      maxsh=k
c     write(6,'(" no. of shells = ",i5)') maxsh 
      return
      end
