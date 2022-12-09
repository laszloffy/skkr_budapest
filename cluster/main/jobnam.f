c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine jobnam(poti,aa1,aa2,aa3,aa4,aa5)
      character*25 poti
      character*30 aa1,aa2,aa3,aa4,aa5
c
      ij=index(poti,' ')-1
      aa1=poti(1:ij)//'.prn'
      aa2=poti(1:ij)//'.pot'
      aa3=poti(1:ij)//'.dos'
      aa4=poti(1:ij)//'.pdos'
      aa5=poti(1:ij)//'.tcpa'
      return
      end
