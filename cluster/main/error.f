c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
c
      subroutine error(n,m,msg)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      character*17 msg
c
      write(6,
     >'(2x,''***'',a,'' equal to '',i5,'' which is less than '',i5)') 
     > msg,m,n
c
      stop
      end
