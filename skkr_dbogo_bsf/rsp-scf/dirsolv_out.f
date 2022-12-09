c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine dirsolv_out(ce,ceh,E_Fermi,delta,
     >                       singrat,urat,drat,cl,soc,kmymax,kmy0,
     >                       tz,urp,urm,curp,curm,dx,rs,ns,g,f,iswitch)
c
c     *****************************************************************
c     * Integration of the radial Dirac-BdG equation by Adams method  *
c     * Integrate outwards!                                           *
c     *****************************************************************
c
      implicit none
c
      include '../param.h'
c
c I/O variables
      logical iswitch
      integer kmymax,kmy0,ns
      real*8 dx,rs,cl,soc,tz,E_Fermi
      real*8 singrat,urat,drat
      complex*16 urp(kmymaxp,kmymaxp,nrad),urm(kmymaxp,kmymaxp,nrad)
      complex*16 curp(kmymaxp,kmymaxp,nrad),curm(kmymaxp,kmymaxp,nrad)
      complex*16 ce,ceh,delta(nrad)
      complex*16 g(dbogomaxp,nrad),f(dbogomaxp,nrad)
c
c internal variables
      integer i,j,n,nit,nitmax,imode,lm,lmp,m,mp,ll,llp
      integer l0,kap0,kappai,kappaj,sk,myi,myj
      real*8 x,x0,dx2,c,csq,csqinv,r,tzoc2,xk0,ac,test,tiny
      real*8 rn,rmid,rnp1
c
      complex*16 urpn(kmymaxp,kmymaxp),urmn(kmymaxp,kmymaxp)
      complex*16 urpnp1(kmymaxp,kmymaxp),urmnp1(kmymaxp,kmymaxp)
      complex*16 urpmid(kmymaxp,kmymaxp),urmmid(kmymaxp,kmymaxp)
      complex*16 curpn(kmymaxp,kmymaxp),curmn(kmymaxp,kmymaxp)
      complex*16 curpnp1(kmymaxp,kmymaxp),curmnp1(kmymaxp,kmymaxp)
      complex*16 curpmid(kmymaxp,kmymaxp),curmmid(kmymaxp,kmymaxp)
      complex*16 deltan(kmymaxp,kmymaxp)
      complex*16 deltamid(kmymaxp,kmymaxp)
      complex*16 deltanp1(kmymaxp,kmymaxp)
      complex*16 deltamat(kmymaxp,kmymaxp,nrad)
      complex*16 deltarloc(nrad)
c
      complex*16 p(kmymaxp,nrad),pp(kmymaxp,nrad)
      complex*16 q(kmymaxp,nrad),qp(kmymaxp,nrad)
      complex*16 ph(kmymaxp,nrad),pph(kmymaxp,nrad)
      complex*16 qh(kmymaxp,nrad),qph(kmymaxp,nrad)
c
      complex*16 pn(kmymaxp),qn(kmymaxp)
      complex*16 pnh(kmymaxp),qnh(kmymaxp)
      complex*16 pmid(kmymaxp),qmid(kmymaxp)
      complex*16 pmidh(kmymaxp),qmidh(kmymaxp)
      complex*16 k1(kmymaxp),m1(kmymaxp)
      complex*16 k2(kmymaxp),m2(kmymaxp)
      complex*16 k3(kmymaxp),m3(kmymaxp)
      complex*16 k4(kmymaxp),m4(kmymaxp)
      complex*16 k1h(kmymaxp),m1h(kmymaxp)
      complex*16 k2h(kmymaxp),m2h(kmymaxp)
      complex*16 k3h(kmymaxp),m3h(kmymaxp)
      complex*16 k4h(kmymaxp),m4h(kmymaxp)
c
      complex*16 psum(kmymaxp),qsum(kmymaxp)
      complex*16 psumh(kmymaxp),qsumh(kmymaxp)
      complex*16 p0(kmymaxp),q0(kmymaxp)
      complex*16 p0h(kmymaxp),q0h(kmymaxp)
      complex*16 p1(kmymaxp),q1(kmymaxp)
      complex*16 p1h(kmymaxp),q1h(kmymaxp)
      complex*16 ppnp1(kmymaxp),qpnp1(kmymaxp)
      complex*16 ppnp1h(kmymaxp),qpnp1h(kmymaxp)
      complex*16 xlms(lmmaxp,lmmaxp,2,2)
      complex*16 deltakmy(kmymaxp,kmymaxp)
c
c
      data nitmax/50/,test/1.0d8/,tiny/1.0d-10/
c
!
      write(6,*) 'Ask for the original version of this file from the 
     >owners of the Copyright, e.g. from laszloffy.andras(at)wigner.hu'
!
      return
      end
