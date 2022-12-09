c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      function gaunt(l,lp,lpp,m,mp,mpp)
c===================
c
c calculate integral over three spherical harmonics y*lm yl'm' yl"m"
c by numerical integration
c subroutine init has to be called first to calculate array yg
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      common/cgaunt/yg(lshape,lmshape)
c
      gaunt=0.d0
      if(m-mp-mpp.ne.0) goto 2
      if(l+lp.lt.lpp) goto 2
      if(l+lpp.lt.lp) goto 2
      if(lp+lpp.lt.l) goto 2
      if(mod(l+lp+lpp,2).eq.1) goto 2
      k1=l*l+l+m+1
      k2=lp*lp+lp+mp+1
      k3=lpp*lpp+lpp+mpp+1
c
      do 1 i=1,lshape
    1 gaunt=gaunt+yg(i,k1)*yg(i,k2)*yg(i,k3)
      return
c
    2 write(6,*) ' gaunt ?',l,lp,lpp,m,mp,mpp
      return
      end
