c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine jobnam(poti,aa1,aa2,aa3,aa4,aa5,aa6,aa7,aa8,aa9)
      character*25 poti
      character*30 aa1,aa2,aa3,aa4,aa5,aa6,aa7
      character*34 aa8,aa9
      ij=index(poti,' ')-1
      aa1=poti(1:ij)//'.prn'
      aa2=poti(1:ij)//'.pot'
      aa3=poti(1:ij)//'.dos'
      aa4=poti(1:ij)//'.pdos'
      aa5=poti(1:ij)//'.dos_imp'
      aa6=poti(1:ij)//'.tcpa'
      aa7=poti(1:ij)//'.mom'
      aa8=poti(1:ij)//'.pdos_imp'
      aa9=poti(1:ij)//'.triplet'
      return
      end
