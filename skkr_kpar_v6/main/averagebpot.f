c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine averagebpot(nsl,nintfc,nbulk,vr,vrl,vrr)
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      dimension vrl(nrad,minprc),vrr(nrad,minprc),vr(nrad,mintfc)
      dimension nsl(minprc)
c
      nunits=nintfc/nbulk
c
      do 1 iu=2,nunits
        li = nbulk*(iu-1)
        do 1 il=1,nbulk
          if(nsl(il+li).ne.nsl(il)) stop 
     >      'radial grid mismatch in averagbpot!!'
 1    continue
c                
      do 10 il=1,nbulk
        do 10 i=1,nsl(il)
          vrl(i,il)=0.d0
 10   continue
c
      do 11 iu=1,nunits
        li = nbulk*(iu-1)
        do 11 il=1,nbulk
          do 11 i=1,nsl(il)
            vrl(i,il) = vrl(i,il) + vr(i,il+li)
 11   continue
c
      do il=1,nbulk
        do i=1,nsl(il)
          vrl(i,il) = vrl(i,il)/nunits
        enddo
      enddo
c
      do 12 iu=1,nunits
        li = nbulk*(iu-1)
        do 12 il=1,nbulk
          do 12 i=1,nsl(il)
            vrl(i,il+li) = vrl(i,il)
            vrr(i,il+li) = vrl(i,il+li)
            vr(i,il+li)  = vrl(i,il+li)
 12   continue
c
      return
      end
c
c======================
      subroutine averagebpotsp(nsl,nintfc,nbulk,vr,vrl,vrr)
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      dimension vrl(nrad,2,minprc),vrr(nrad,2,minprc),vr(nrad,2,mintfc)
      dimension nsl(minprc)
c
      nunits=nintfc/nbulk
c                
      do 10 il=1,nbulk
        do 10 is=1,2
        do 10 i=1,nsl(il)
          vrl(i,is,il)=0.d0
 10   continue
c
      do 11 iu=1,nunits
        li = nbulk*(iu-1)
        do 11 il=1,nbulk
          do 11 is=1,2
          do 11 i=1,nsl(il)
            vrl(i,is,il) = vrl(i,is,il) + vr(i,is,il+li)
 11   continue
c
      do il=1,nbulk
        do i=1,nsl(il)
          vrl(i,1,il) = vrl(i,1,il)/nunits
          vrl(i,2,il) = vrl(i,2,il)/nunits
        enddo
      enddo
c
      do 12 iu=1,nunits
        li = nbulk*(iu-1)
        do 12 il=1,nbulk
          do 12 is=1,2
          do 12 i=1,nsl(il)
            vrl(i,is,il+li) = vrl(i,is,il)
            vrr(i,is,il+li) = vrl(i,is,il+li)
            vr(i,is,il+li)  = vrl(i,is,il+li)
 12   continue
c
      return
      end
