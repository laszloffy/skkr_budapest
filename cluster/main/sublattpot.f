c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine sublattpot(nintfc,ns,vr,iesublatt)
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
c  isublatt:  number of sublattices
c  nsublatt(is): number of sites in sublattice # is
c  iesubl(i,is): label of site # i in sublattice # is
c
      dimension ns(mintfc),iesublatt(mintfc)
      dimension nsublatt(mintfc),iesubl(mintfc,mintfc)
      dimension vr(nrad,mintfc),vrsubl(nrad)
c
      isublatt=1
      nsublatt(1)=1
      iesubl(1,1)=1
      do li=2,nintfc
        do is=1,isublatt
        do i=1,nsublatt(is)
          if(iesublatt(li).eq.iesublatt(iesubl(i,is))) then
            nsublatt(is)=nsublatt(is)+1
            iesubl(nsublatt(is),is)=li
            goto 10
          end if
        end do
        end do
        isublatt=isublatt+1
        nsublatt(isublatt)=1
        iesubl(1,isublatt)=li
  10    continue
      end do
c     write(6,*) isublatt
c     write(6,*) (nsublatt(is),is=1,isublatt)
c     do is=1,isublatt
c       write(6,*) (iesubl(i,is),i=1,nsublatt(is))
c     end do
        
      do is=1,isublatt
        do j=1,nrad
          vrsubl(j)=0.0d0
        end do
        do i=1,nsublatt(is)
          li=iesubl(i,is)
          do j=1,ns(li)
            vrsubl(j)=vrsubl(j)+vr(j,li)
          end do
        end do
        do j=1,nrad
          vrsubl(j)=vrsubl(j)/nsublatt(is)
        end do
        do i=1,nsublatt(is)
          li=iesubl(i,is)
          do j=1,ns(li)
            vr(j,li)=vrsubl(j)
          end do
        end do
      end do
c
      return
      end
