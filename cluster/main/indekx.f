c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine indekx(n,a,ind)
c======================
c
c find index array ind such that a(ind(n+1)) >= a(ind(n))
c see w. press et al, numerical recipes, cambridge 1986
c
      implicit real*8 (a-h,o-z)
      dimension a(n),ind(n)
c
      do 1 i=1,n
    1 ind(i)=i
c
      l=n/2+1
      ir=n
    2 continue
      if(l.gt.1) then
      l=l-1
      indxt=ind(l)
      q=a(indxt)
         else
      indxt=ind(ir)
      q=a(indxt)
      ind(ir)=ind(1)
      ir=ir-1
      if(ir.eq.1) then
      ind(1)=indxt
      return
      endif
      endif
c
      i=l
      j=l+l
c
    3 if(j.le.ir) then
      if(j.lt.ir) then
      if(a(ind(j)).lt.a(ind(j+1))) j=j+1
      endif
      if(q.lt.a(ind(j))) then
      ind(i)=ind(j)
      i=j
      j=j+j
        else
      j=ir+1
      endif
      goto 3
      endif
      ind(i)=indxt
      goto 2
      end
