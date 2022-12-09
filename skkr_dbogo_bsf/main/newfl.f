c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine newfl(itscf,kmy0,kmymax,kmymmax,ne,nbulk,conc,
     &           dosa,dosb,qca,qcb,qva,qvb,za,zb,efermi,defermi)
c=====================
c
      implicit real*8 (a-h,o-z)
      include '../param.h'
c
      logical cpalay
      dimension dosa(kmy0:kmymmax,mintfc,me),
     &          dosb(kmy0:kmymmax,mintfc,me),
     &          qca(mintfc),qcb(mintfc),qva(mintfc),qvb(mintfc),
     &          conc(mintfc),za(mintfc),zb(mintfc)
      common/test/itest
      data tiny/1.0d-6/
c
      efold=efermi
c
      dos = 0.d0
      q = 0.d0
      z  = 0.d0
      do 30 il=1,nbulk
        cpalay=1.d0-conc(il).gt.tiny
c
        q = q + (qca(il)+qva(il))*conc(il)
        z = z + za(il)*conc(il)
        if(cpalay) then
          q = q + (qcb(il)+qvb(il))*(1.d0-conc(il))
          z = z + zb(il)*(1.d0-conc(il))
        endif
c
        do 30 k=kmy0,kmymax
          dos = dos + dosa(k,il,ne)*conc(il)
          if(cpalay) dos = dos +dosb(k,il,ne)*(1.d0-conc(il))
 30   continue
c
      defermi=(z-q)/dos
      efermi=efold+defermi
c
      if(itest.ge.1) then
        write(6,'(/4e15.6)') z,q,dos,defermi
        write(6,
     > '('' I'',i4,''  EFermi'',3x,f15.10,5x,'' DEFermi'',3x,f15.10)')
     > itscf,efermi,defermi
      end if
c 
      return
      end
