c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine dmag(pp,qp,p,q,ce,r,csqinv,soc,kmymax,urp,urm)
c
c     ***********************************************************
c     * Derivatives of the wavefunctions from                   *
c     * the radial dirac equation                               *
c     ***********************************************************
c
      implicit none
c
      include '../param.h'
c
c I/O variables
      integer kmymax
      real*8 r,csqinv,soc
      complex*16 ce
      complex*16 p(kmymaxp),pp(kmymaxp)
      complex*16 q(kmymaxp),qp(kmymaxp)
      complex*16 urp(kmymaxp,kmymaxp),urm(kmymaxp,kmymaxp)
c
c internal variables
      integer i,j
      real*8 xkap
      complex*16 re,rep
c
      re=r*ce
      rep=csqinv*re+r
c
      do i=1,kmymax
        xkap=soc*(1.0d0+kapdex(i))-1.0d0
        qp(i)=xkap*q(i)-re*p(i)
        pp(i)=-xkap*p(i)+rep*q(i)
        do j=1,kmymax
          qp(i)=qp(i)+urp(i,j)*p(j)
          pp(i)=pp(i)-urm(i,j)*q(j)
        end do
      end do
c
      return
      end
c
c
c
c
      subroutine dmagb(pp,qp,p,q,pph,qph,ph,qh,ce,ceh,r,csqinv,
     >                 soc,kmymax,urp,urm,curp,curm,deltar)
c
c     ***********************************************************
c     * Derivatives of the wavefunctions from                   *
c     * the radial dirac equation                               *
c     ***********************************************************
c
      implicit none
c
      include '../param.h'
c
c I/O variables
      integer kmymax
      real*8 r,csqinv,soc
      complex*16 ce,ceh
      complex*16 p(kmymaxp),pp(kmymaxp)
      complex*16 q(kmymaxp),qp(kmymaxp)
      complex*16 urp(kmymaxp,kmymaxp),urm(kmymaxp,kmymaxp)
      complex*16 ph(kmymaxp),pph(kmymaxp)
      complex*16 qh(kmymaxp),qph(kmymaxp)
      complex*16 curp(kmymaxp,kmymaxp),curm(kmymaxp,kmymaxp)
      complex*16 deltar(kmymaxp,kmymaxp),cdeltar(kmymaxp,kmymaxp)
c
c internal variables
      integer i,j
      real*8 xkap
      complex*16 re,rep,reh,reph
c
      re=r*ce
      rep=csqinv*re+r
      reh=r*ceh
      reph=csqinv*reh+r
      cdeltar=conjg(transpose(deltar(:,:)))
c
      do i=1,kmymax
        xkap=soc*(1.0d0+kapdex(i))-1.0d0
        qp(i)=xkap*q(i)-re*p(i)
        pp(i)=-xkap*p(i)+rep*q(i)
        qph(i)=xkap*qh(i)-reh*ph(i)
        pph(i)=-xkap*ph(i)+reph*qh(i)
        do j=1,kmymax
          qp(i)=qp(i)+urp(i,j)*p(j)+deltar(i,j)*ph(j)
          pp(i)=pp(i)-urm(i,j)*q(j)-deltar(i,j)*qh(j)*csqinv
          qph(i)=qph(i)+curp(i,j)*ph(j)-cdeltar(i,j)*p(j)
          pph(i)=pph(i)-curm(i,j)*qh(j)+cdeltar(i,j)*q(j)*csqinv 
        end do
      end do
c
      return
      end
