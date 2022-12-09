c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine exchange(rho1,rho2,rho,iexch,v1,v2,exc)
c========================
c
c rho1, rho2, rho=rho1+rho2 are the true charge densities!
c 
      implicit real*8 (a-h,o-z)
c
      real*8 rhop(2), rhopp(2)
c
      data alph0/0.6108870577108572d0/
      data small/1.d-15/,alpha/1.0d0/
c
      v1=0.d0
      v2=0.d0
      exc=0.d0
c
      if(rho1.lt.small) return
      if(rho2.lt.small) return
      if(rho.lt.small) return
c
      rs=(0.2387324146378430d0/rho)**(1.d0/3.d0)
      v0=-2.d0*alph0/rs
c
c no gradient corrections !!!
c
      rhop(1)=0.0d0
      rhop(1)=0.0d0
      rhopp(2)=0.0d0
      rhopp(2)=0.0d0
      rr=0.0d0
c
      if(iexch.eq.11) then
c Xalpha
        v1=1.5d0*alpha*v0
        v2=v1
      else if(iexch.eq.12) then
c Gunarsson-Lundquist 
        beta=1.d0+0.0545d0*rs*dlog((rs+11.4d0)/rs)
        delta=1.d0-0.036d0*rs-1.36d0*rs/(1.d0+10.d0*rs)
        dkhi=(rho1-rho2)/rho
        v1=v0*(beta+delta*dkhi/(1.d0+0.297d0*dkhi)/3.d0)
        v2=v0*(beta-delta*dkhi/(1.d0-0.297d0*dkhi)/3.d0)
      else if(iexch.eq.1) then
        call xcpot(iexch,RHO1,RHO2,RHO,RHOP,RHOPP,RR,V1,V2,EXC)
      else if(iexch.eq.2) then
        call xcpot(iexch,RHO1,RHO2,RHO,RHOP,RHOPP,RR,V1,V2,EXC)
      else if(iexch.eq.3) then
        call xcpot(iexch,RHO1,RHO2,RHO,RHOP,RHOPP,RR,V1,V2,EXC)
      else if(iexch.eq.4) then
        call xcpot(iexch,RHO1,RHO2,RHO,RHOP,RHOPP,RR,V1,V2,EXC)
      else if(iexch.eq.6) then
        call xcpot(iexch,RHO1,RHO2,RHO,RHOP,RHOPP,RR,V1,V2,EXC)
      else if(iexch.eq.7) then
        call xcpot(iexch,RHO1,RHO2,RHO,RHOP,RHOPP,RR,V1,V2,EXC)
      else if(iexch.eq.8) then
        call xcpot(iexch,RHO1,RHO2,RHO,RHOP,RHOPP,RR,V1,V2,EXC)
      else if(iexch.eq.9) then
        call xcpot(iexch,RHO1,RHO2,RHO,RHOP,RHOPP,RR,V1,V2,EXC)
      else
        stop '*** Un-implemented Exchange ***'
      end if
c
      return
      end
      subroutine xcpot(ixc,rho1,rho2,rho,rhop,rhopp,rr,v1,v2,exc)
C   ******************************************************************
C   *                                                                *
C   *    Calculates exchange-correlation potential according to the  *
C   *    value of IXC = IXCH+1. The constants used in the various    *
C   *    expressions are set in SETXCP and in the data statement     *
C   *    below.                                                      *
C   *                                                                *
C   *    INPUT:RHO totale charge densiti,                            *
C   *          RHO1,RHO2 spin densities,                             *
C   *          RHOP(1),RHOP(2) first derivates,                      *
C   *          RHOPP(1),RHOPP(2) second derivates,                   *
C   *          R radius,                                             *
C   *          IXC see 'setxcp.f'                                    *
C   *                                                                *
C   *    HLS: 15-08-93                                               *
C   ******************************************************************
      implicit real*8(a-h,o-z)
      parameter(told=1.d-10,toldd=1.d-10)
      common/sxcp/xccp,xccf,xcrp,xcrf,xalpha,oth,fth,aa,bb,alpm,blpm,
     1            fourpi,aw,bw,cw,betab,fif,pp1,pp2,nss
      common/cap/acau,acap,bcau,bcap,ccau,ccap,dcau,dcap,fca,ocau,
     1           ocap,pcau,pcap,qcau,qcap,rcau,rcap,scau,scap,tcau,
     2           tcap,ucau,ucap
      dimension rhop(2),rhopp(2)
      data ap,xp0,bp,cp,qp,cp1,cp2,cp3/0.0621814,-0.10498,3.72744,
     1     12.9352,6.1519908,1.2117833,1.1435257,-0.031167608/
      data af,xf0,bf,cf,qf,cf1,cf2,cf3/0.0310907,-0.32500,7.060428,
     2     18.0578,4.7309269,2.9847935,2.7100059,-0.1446006/
  100 format(/,' XCPOT:**  RR less than 1. D-10 for IXC =',I3)
C
      if(rho.lt.1.d-20) then
         v1=0.D0
         v2=0.D0
         exc=0.D0
         return
      endif
      pi=fourpi/4.D0
      rs1=(fourpi*rho/3.)**oth
      rs=1./rs1
      go to (1,2,1,3,4,5,6,6,3),ixc
C
C     Barth-Hedin  J. PHYS. C5,1629(1972)
C
    1 rsf=rs/xcrf
      rsf2=rsf*rsf
      rsf3=rsf2*rsf
      rsp=rs/xcrp
      rsp2=rsp*rsp
      rsp3=rsp2*rsp
      fcf=(1.d0+rsf3)*log(1.d0+1.d0/rsf)+.5d0*rsf-rsf2-oth
      fcp=(1.d0+rsp3)*log(1.d0+1.d0/rsp)+.5d0*rsp-rsp2-oth
      epscp=-xccp*fcp
      epscf=-xccf*fcf
      epsxp=-.91633059d0/rs
      cny=5.1297628d0*(epscf-epscp)
      x=rho1/rho
      fx=(x**fth+(1.d0-x)**fth-aa)/bb
      exc=epsxp+epscp+fx*(cny+fth*epsxp)/5.1297628d0
      ars=-1.22177412d0/rs+cny
      brs=-xccp*log(1.d0+xcrp/rs)-cny
      trx1=(2.d0*x)**oth
      v1=ars*trx1+brs
      trx2=(2.d0*rho2/rho)**oth
      v2=ars*trx2+brs
      return
C
C     Slater X-Alpha
C
    2 exc=-0.75d0*xalpha*(0.5d0*rho)**oth
      v1=-xalpha*(rho1)**oth
      v2=-xalpha*(rho2)**oth
      return
C
C     Vosko-Wilk-Nusair  Can. J. Phys. 58,1200(1980)
C
    3 x=sqrt(rs)
      xpx=x*x+bp*x+cp
      xfx=x*x+bf*x+cf
      s=(rho2-rho1)/rho
      sp=1.d0+s
      sm=1.d0-s
      s4=s**4-1.d0
      fs=(sp**fth+sm**fth-2.d0)/aa
      beta=1.d0/(2.74208d0+3.182d0*x+0.09873d0*x*x+0.18268d0*x**3)
      dfs=fth*(sp**oth-sm**oth)/aa
      dbeta=-(0.27402d0*x+0.09873d0+1.591d0/x)*beta**2
      atnp=atan(qp/(2.d0*x+bp))
      atnf=atan(qf/(2.d0*x+bf))
      ecp=ap*(log(x*x/xpx)+cp1*atnp-cp3*(log((x-xp0)**2/xpx)+cp2*atnp
     1))
      ecf=af*(log(x*x/xfx)+cf1*atnf-cf3*(log((x-xf0)**2/xfx)+cf2*atnf
     1))
      ec=ecp+fs*(ecf-ecp)*(1.d0+s4*beta)
      tp1=(x*x+bp*x)/xpx
      tf1=(x*x+bf*x)/xfx
      ucp=ecp-ap/3.d0*(1.d0-tp1-cp3*(x/(x-xp0)-tp1-xp0*x/xpx))
      ucf=ecf-af/3.d0*(1.d0-tf1-cf3*(x/(x-xf0)-tf1-xf0*x/xfx))
      uc0=ucp+(ucf-ucp)*fs
      uc20=uc0+(ecf-ecp)*sm*dfs
      uc10=uc0-(ecf-ecp)*sp*dfs
      duc=(ucf-ucp)*beta*s4*fs+(ecf-ecp)*(-rs/3.d0)*dbeta*s4*fs
      duc2=duc+(ecf-ecp)*beta*sm*(4.d0*s**3*fs+s4*dfs)
      duc1=duc-(ecf-ecp)*beta*sp*(4.d0*s**3*fs+s4*dfs)
      uc1=uc10+duc1
      uc2=uc20+duc2
      epx=-0.91633059d0/rs*(1.d0+fth*fs/5.1297628d0)
      amyx2=-1.22177412d0/rs*sp**oth
      amyx1=-1.22177412d0/rs*sm**oth
      exc=ec+epx
      v1=uc1+amyx1
      v2=uc2+amyx2
      if(ixc.eq.4) return   
      excl=exc
      vxcl1=v1
      vxcl2=v2
      Go to 8
C
C     Langreth-Perdew-Mehl non-local corrections
C     Phys. Rev. B28,1809(1983)
C     The local potential is Barth-Hedin
C
    4 if(rr.lt.1.d-10) then
         write(6,100) ixc
         stop
      endif
      rsf=rs/xcrf
      rsf2=rsf*rsf
      rsf3=rsf2*rsf
      rsp=rs/xcrp
      rsp2=rsp*rsp
      rsp3=rsp2*rsp
      fcf=(1.d0+rsf3)*log(1.d0+1.d0/rsf)+.5d0*rsf-rsf2-oth
      fcp=(1.d0+rsp3)*log(1.d0+1.d0/rsp)+.5d0*rsp-rsp2-oth
      epscp=-xccp*fcp
      epscf=-xccf*fcf
      epsxp=-.91633059d0/rs
      cny=5.1297628d0*(epscf-epscp)
      x=rho1/rho
      fx=(x**fth+(1.d0-x)**fth-aa)/bb
      excloc=epsxp+epscp+fx*(cny+fth*epsxp)/5.1297628d0
      ars=-1.22177412d0/rs+cny
      brs=-xccp*log(1.d0+xcrp/rs)-cny
      trx1=(2.d0*x)**oth
      v1loc=ars*trx1+brs
      trx2=(2.d0*rho2/rho)**oth
      v2loc=ars*trx2+brs
c
C     Non-local corrections
C
      rhod=rhop(1)+rhop(2)
      rhodd=rhopp(1)+rhopp(2)
      if(abs(rhod).gt.told) then
         flpm=blpm*abs(rhod)/rho**(7.d0/6.d0)
         elpm=2.d0*exp(-flpm)
         excnl=alpm*(elpm-7.d0/9.d0)*rhod**2/rho**(7.d0/3.d0)
         xn=rhod/rho
         yn=(2.d0*rhod/rr+rhodd)/rhO
         p1=(7.d0/9.d0)*(yn-2.d0*xn*xn/3.d0)
         p2=(1.d0-flpm/2.d0)*yn
         p3=(2.d0/3.d0-11.d0*flpm/6.d0+7.d0*flpm*flpm/12.d0)*xn*xn
         p4=-0.5d0*flpm*(flpm-3.d0)*rhod*rhodd/rho/abs(rhod)
         p5=elpm*(p2-p3+p4)
         vnl=2.d0*alpm*(p1-p5)/rho**oth
         exc=excloc+excnl
         v1=v1loc+vnl
         v2=v2loc+vnl
      else
         exc=excloc
         v1=v1loc
         v2=v2loc
      endif
      return
C
C     Wigner expression
C
    5 rs78=1.d0/(rs+7.8d0)
      exc=-0.916d0*rs1-0.88d0*rs78
      v1=cw*rs78*rs78-aw*rs1-bw*rs78
      v2=v1
      return
C
C     Ceperley-Alder. Parametrization by Perdew and Zunger.
C
    6 if(rs.ge.1.d0) then
         sqrtrs=sqrt(rs)
         denom1=1.d0/(1.d0+acau*sqrtrs+bcau*rs)
         ecu=-0.2846d0*denom1
         vcu=ecu*(1.d0+ccau*sqrtrs+dcau*rs)*denom1
      else
         rslog=log(rs)
         rsln=rs*rslog
         ecu=-ocau+pcau*rslog-qcau*rs+rcau*rsln
         vcu=-scau+pcau*rslog-tcau*rs+ucau*rsln
      endif
      ec=ecu
      vc1=vcu
      vc2=vc1
      exu=-0.9164d0*rs1
      ex=exu
      vx1=fca*exu
      vx2=fca*exu
      xi=(rho1-rho2)/rho
      if(nss.eq.2) then
       if(rs.ge.1.d0) then
          sqrtrs=sqrt(rs)
          denom2=1.d0/(1.d0+acap*sqrtrs+bcap*rs)
          ecp=-0.1686d0*denom2
          vcp=ecp*(1.d0+ccap*sqrtrs+dcap*rs)*denom2
       else           
          rslog=log(rs)
          rsln=rs*rslog          
          ecp=-ocap+pcap*rslog-qcap*rs+rcap*rsln
          vcp=-scap+pcap*rslog-tcap*rs+ucap*rsln
       endif
       ec=ecu+fxz(xi)*(ecp-ecu)
       vc1=vcu+fxz(xi)*(vcp-vcu)+(ecp-ecu)*(1.d0-xi)*fxzd(xi)
       vc2=vcu+fxz(xi)*(vcp-vcu)-(ecp-ecu)*(1.d0+xi)*fxzd(xi)
       exxp=exu*2.d0**oth
       ex=exu+(exxp-exu)*fxz(xi)
       vx1=-1.22177412d0/rs*((1.d0+xi)**oth)
       vx2=-1.22177412d0/rs*((1.d0-xi)**oth)
      end if
      exc=ex+ec
      v1=vx1+vc1
      v2=vx2+vc2
      if(ixc.eq.7) return
C
C     Perdew-Wang generalized gradient approximation,
C                 Phys. Rev. B33,8800(1986),
C                 Phys. Rev. B33,8822(1986).
C
      ecl=ec
      vcl1=vc1
      vcl2=vc2
C
C     Total exchange by Perdew-Wang
C
      if(nss.eq.1) then
         rhx=rho
         rhod=rhop(1)+rhop(2)
         rhodd=rhopp(1)+rhopp(2)
         call pewax(rhx,rhod,rhodd,rr,ex,vx)
         vx1=vx
         vx2=vx
      else
         rhx=2.d0*rho1
         rhod=2.d0*rhop(1)
         rhodd=2.d0*rhopp(1)
         call pewax(rhx,rhod,rhodd,rr,ex1,vx)
         vx1=vx
         rhx=2.d0*rho2
         rhod=2.d0*rhop(2)
         rhodd=2.d0*rhopp(2)
         call pewax(rhx,rhod,rhodd,rr,ex2,vx)
         vx2=vx
         ex=(rho1*ex1+rho2*ex2)/rho
      endif
C
C     Local correlation by Ceperley-Alder 
C
C     Non local correlation by Perdew
C
      call perdew(rho1,rho2,rhop,rhopp,rr,vcnl1,vcnl2,ecnl)
C
      ec=ecl+ecnl
      exc=ex+ec
      v1=vx1+vcl1+vcnl1
      v2=vx2+vcl2+vcnl2
      return
C
C     Becke-Perdew non-local approximation,
C                 Phys. Rev. A38,3098(1988),
C                 Phys. Rev. B33,8822(1986).
C
C     Local exchange-correlation by Vosko-Wilk-Nusair
C     Non local exchange by Becke
C
    8 rxb=rho1
      rxbd=rhop(1)
      rxbdd=rhopp(1)
      call becke(rho,rxb,rxbd,rxbdd,rr,vbnl,ebnl)
      exnl1=ebnl
      vxnl1=vbnl
      rxb=rho2
      rxbd=rhop(2)
      rxbdd=rhopp(2)
      call becke(rho,rxb,rxbd,rxbdd,rr,vbnl,ebnl)
      exnl2=ebnl
      vxnl2=vbnl
      exnl=exnl1+exnl2
C
C     Non local correlation by Perdew
C
      call perdew(rho1,rho2,rhop,rhopp,rr,vcnl1,vcnl2,ecnl)
C
      exc=exnl+ecnl+excl
      v1=vxcl1+vxnl1+vcnl1
      v2=vxcl2+vxnl2+vcnl2
      return
      end
C********************************************************************
C   Calculation of the Perdew non local correlation potential
C********************************************************************
      subroutine perdew(rho1,rho2,rhop,rhopp,rr,vcnl1,vcnl2,ecnl)
      implicit real*8(a-h,o-z)
      common/sxcp/xccp,xccf,xcrp,xcrf,xalpha,oth,fth,aa,bb,alpm,blpm,
     1            fourpi,aw,bw,cw,betab,fif,pp1,pp2,nss
      dimension rhop(2),rhopp(2)
C
      rho=rho1+rho2
      xi=(rho1-rho2)/rho
      fca=4.d0/3.d0
      rs1=(fourpi*rho/3.)**oth
      rs=1./rs1
      rhod=rhop(1)+rhop(2)
      rhodd=rhopp(1)+rhopp(2)
      ddn=rhodd+2.d0*rhod/rr
      dn=rhod
      dna=dabs(dn)
      rho43=rho**fca
      fcps=fcpp(rs)
      dp=fcz(xi)
      dpm=2.d0/dp
      if(dna.gt.1.d-20) then
       dn2=dna*dna
       dnddn=dna*rhodd
       fccr=-fourpi/9.d0*(rs**4.d0)
       rho76=rho**pp2
       rho73=rho76*rho76
       fcpds=fccr*fcppd(rs)
       fi=fif*dna/(fcps*rho76)
       fi2=fi*fi
       eafi=exp(-fi)
       tn1=(2.d0-fi)*ddn
       tn2=(fca-pp1*fi+pp2*fi2)*dn2/rho
       tn3=fi*(fi-3.d0)*dnddn/dnA
       tn4=(fi2-fi-1.d0)*dn2
       if(nss.eq.2) then
        rho13=rho**oth
        oth2=2.d0*oth
        rhoh1=rho1**oth2
        rhoh2=rho2**oth2
        tn5=(1.d0-fi)*rho2*dn2-(2.d0-fi)*rho*rhop(2)*dn
C
        tn6=1.32283421d0*rho13*(rhoh1-rhoh2)/(rho**3.d0)/(dp**2)
c
        vcnl=dpm*eafi*(fcps*(tn1-tn2+tn3-tn5*tn6)-fcpds*tn4)
        vcnl1=-vcnl/rho43
        tn5=(1.d0-fi)*rho1*dn2-(2.d0-fi)*rho*rhop(1)*dn
        vcnl=dpm*eafi*(fcps*(tn1-tn2+tn3+tn5*tn6)-fcpds*tn4)
        vcnl2=-vcnl/rho43
       else
        vcnl=dpm*eafi*(fcps*(tn1-tn2+tn3)-fcpds*tn4)
        vcnl1=-vcnl/rho43
        vcnl2=vcnl1
       endif
       ecnl=dpm*eafi*fcps*dn2/rho73
      else
       ecnl=0.d0
C       VCNL1=-2.D0*DPM*FCPS*DDN/RHO43
        vcnl1=0.d0
        vcnl2=0.d0
C       VCNL2=-2.D0*DPM*FCPS*DDN/RHO43
      endif
      return
      end
C********************************************************************
C   Calculation of the Perdew-Wang exchange potential
C********************************************************************
      subroutine pewax(rho,rhod,rhodd,rr,ex,vx)
      implicit real*8(a-h,o-z)
      parameter(one=1)
C
      oth=1.d0/3.d0
      fca=4.d0/3.d0
      pi=acos(-one)
      ckf=(3.d0*rho*pi*pi)**oth
      ckf2=2.d0*ckf
      tk=rho*(ckf2**2.d0)
      ddn=rhodd+2.d0*rhod/rr
      t=ddn/tk
      dn=rhod
      dna=dabs(dn)
      if(dna.gt.1.d-20) then
       dnddn=dna*rhodd
       sk=ckf2*rho
       s=dna/sk
       uk=tk*sk
       u=dnddn/uk
       s2=s*s
       fs=fpw(s)
       fds=fpwd(s)
       fdds=fpwdd(s)
       ex=-1.47711753d0*(rho**oth)*fs
       g1=u/s2/fca-s-t/s/fca
       g2=s2-u/s/fca
       VX=-1.96949004D0*(RHO**OTH)*(FS+G1*FDS+G2*FDDS)
      ELSE
       EX=-1.47711753D0*(RHO**OTH)
       VX=-1.47711753D0*FCA*(RHO**OTH)
c       VX=-1.47711753D0*(FCA-0.1728D0*T)*(RHO**OTH)
      END IF
      RETURN
      END
C********************************************************************
C   Calculation of the Becke non local exchange potential
C********************************************************************
      subroutine becke(rho,rx,rxd,rxdd,rr,vxnl,exnl)
      implicit real*8(a-h,o-z)
      common/sxcp/xccp,xccf,xcrp,xcrf,xalpha,oth,fth,aa,bb,alpm,blpm,
     1            fourpi,aw,bw,cw,betab,fif,pp1,pp2,nss
C
      fca=4.d0/3.d0
      rh43=rx**fca
      dn=rxd
      dna=dabs(dn)
      ddn=rxdd+2.d0*rxd/rr
      if(dna.gt.0.5d-20) then
       dnddn=dna*rxdd
       rh3=rx**3.d0
       rh13=rx**oth
       rhr=rh43/rho
       xh=dna/rh43
       xh2=xh*xh
       xh3=xh*xh2
       ffh=sqrt(1.d0+xh*xh)
       fh=1.d0/ffh
       yh=log(xh+ffh)
       gh=1.d0+6.d0*betab*xh*yh
       exnl=-betab*rhr*xh2/gh
       exnl=2.d0*exnl
C
       th1=fca*xh2
       th2=(1.d0+3.d0*betab*xh*(yh-xh*fh))*ddn/rh43/gh
       th3=6.d0*betab*(dnddn/rh3-fca*xh3)
       th4=3.d0*(yh/gh)*(1.d0+2.d0*betab*xh*yh)
       th5=4.d0*xh*fh*(1.d0-3.d0*betab*xh2*fh)/gh
       th6=xh*(fh**3.d0)/gh+(th4+th5)/gh
       th7=(th1+th3*th6)*rh13-2.d0*th2
       vxnl=-2.d0*betab*th7/gh
      else
       exnl=0.d0
C       vxnl=4.d0*BETAB*DDN/RH43
        vxnl=0.d0
      endif
      return
      end
C********************************************************************
C********************************************************************
C     Function Perdew-Wang F(s)
C********************************************************************
      function fpw(x)
      implicit real*8(a-h,o-z)
      common/pwb/apw,bpw,cpw,pmpw,alfa,beta,gamma,delta,beta1,
     1           ccd0,ccd1,ccd2,ccd3,ccd4
C
      f=(1.d0+apw*(x**2.d0)+bpw*(x**4.d0)+cpw*(x**6.d0))**pmpw
      fpw=f
      return
      end
C********************************************************************
C     Function derivate Perdew- Wang F'(s)
C********************************************************************
      function fpwd(x)
      implicit real*8(a-h,o-z)
      common/pwb/apw,bpw,cpw,pmpw,alfa,beta,gamma,delta,beta1,
     1           ccd0,ccd1,ccd2,ccd3,ccd4
C
      w=(1.d0+apw*(x**2.d0)+bpw*(x**4.d0)+cpw*(x**6.d0))**(pmpw-1.d0)
      f=pmpw*w*(2.d0*apw*x+4.d0*bpw*(x**3.d0)+6.d0*cpw*(x**5.d0))
      fpwd=f
      return
      end
C********************************************************************
C     Function second-order derivate Perdew-Wang F''(s)
C********************************************************************
      function fpwdd(x)
      implicit real*8(a-h,o-z)
      common/pwb/apw,bpw,cpw,pmpw,alfa,beta,gamma,delta,beta1,
     1           ccd0,ccd1,ccd2,ccd3,ccd4
C
      q=(2.d0*apw+12.d0*bpw*(x**2.d0)+30.d0*cpw*(x**4.d0))
      v=(1.d0+apw*(x**2.d0)+bpw*(x**4.d0)+cpw*(x**6.d0))**(pmpw-1.d0)
      p=pmpw*v*q
      r=((2.d0*apw*x+4.d0*bpw*(x**3.d0)+6.d0*cpw*(x**5.d0))**2.d0)
      t=(1.d0+apw*(x**2.d0)+bpw*(x**4.d0)+cpw*(x**6.d0))**(pmpw-2.d0)
      f=pmpw*(pmpw-1.d0)*t*r+p
      fpwdd=f
      return
      end
C********************************************************************
C     Function Perdew C(n)
C********************************************************************
      function fcpp(x)
      implicit real*8(a-h,o-z)
      common/pwb/apw,bpw,cpw,pmpw,alfa,beta,gamma,delta,beta1,
     1           ccd0,ccd1,ccd2,ccd3,ccd4
C
      f1=0.002568d0+alfa*x+beta*(x**2.d0)
      f2=1.d0+gamma*x+delta*(x**2.d0)+beta1*(x**3.d0)
      f=0.001667d0+f1/f2
      fcpp=f
      return
      end
C********************************************************************
C     Function derivate Perdew C'(n)
C********************************************************************
      function fcppd(x)
      implicit real*8(a-h,o-z)
      common/pwb/apw,bpw,cpw,pmpw,alfa,beta,gamma,delta,beta1,
     1           ccd0,ccd1,ccd2,ccd3,ccd4
C
      f1=ccd0+ccd1*x+ccd2*(x**2.d0)+ccd3*(x**3.d0)+ccd4*(x**4.d0)
      f2=1.d0+gamma*x+delta*(x**2.d0)+beta1*(x**3.d0)
      f22=f2*f2
      f=f1/f22
      fcppd=f
      return
      end
C********************************************************************
C     Function Perdew-Zunger F(xi)
C********************************************************************
      function fxz(x)
      implicit real*8(a-h,o-z)
C
      fca=4.d0/3.d0
      f1=((1.d0+x)**fca)+((1.d0-x)**fca)-2.d0
      f2=0.5198421d0
      fxz=f1/f2
      return
      end
C********************************************************************
C     Function derivate Perdew-Zunger F'(xi)
C********************************************************************
      function fxzd(x)
      implicit real*8(a-h,o-z)
C
      oth=1.d0/3.d0
      f1=((1.d0+x)**oth)-((1.d0-x)**oth)
      f2=2.5648814d0
      fxzd=f1*f2
      return
      end
C********************************************************************
C     Function Perdew D(xi)
C********************************************************************
      function fcz(x)
      implicit real*8(a-h,o-z)
C
      cpc=5.d0/3.d0
      f1=((1.d0+x)**cpc)+((1.d0-x)**cpc)
      f2=f1/2.d0
      fcz=sqrt(f2)
      return
      end
      subroutine setxcp(ixc,ns)
C   ******************************************************************
C   *                                                                *
C   *    Initialize constants and text for XCPOT.                    *
C   *                                                                *
C   *   *On entry:                                                   *
C   *                                                                *
C   *    EXCHF : Slater exchange factor equal to 1 for full exchange *
C   *            currently set to permanent 1.0 /bu 29.07.93/        *
C   *    IXC   :=1        Barth-Hedin                                *
C   *            2        Slater X-Alpha                             *
C   *            3        Barth-Hedin-Janak                          *
C   *            4        Vosko-Wilk-Nusair                          *
C   *            5        Langreth-Perdew-Mehl                       *
C   *            6        Wigner exchange                            *
C   *            7        Ceperley-Alder                             *
C   *            8        Perdew-Wang                                *
C   *            9        Becke-Perdew                               *
C   *    NS    : Number of spins                                     *
C   *                                                                *
C   *    HLS:  14-11-92                                              *
C   ******************************************************************
      implicit real*8(a-h,o-z)
      parameter(one=1)
      character txch*3
      common/sxcp/xccp,xccf,xcrp,xcrf,xalpha,oth,fth,aa,bb,alpm,blpm,
     1            fourpi,aw,bw,cw,betab,fif,pp1,pp2,nss
      common/cap/acau,acap,bcau,bcap,ccau,ccap,dcau,dcap,fca,ocau,
     1           ocap,pcau,pcap,qcau,qcap,rcau,rcap,scau,scap,tcau,
     2           tcap,ucau,ucap
      common/pwb/apw,bpw,cpw,pmpw,alfa,beta,gamma,delta,beta1,
     1           ccd0,ccd1,ccd2,ccd3,ccd4
  101 format(/,' setxcp:** Spin polarization not implemented for IXC',
     1       ' =',I3,' Xcpot =',A)
  102 format(/,' SETXCP:',3X,'Slater exchange alpha =',F10.6)
  103 format(/,' SETXCP:** IXC =',I3,' not implemented')
C
      pi=acos(-one)
      fourpi=4.d0*pi
      oth=1.d0/3.d0
      nss=ns
c
      exchf=1.0d0
c
C
      if(ixc.eq.1) then
C
C        Barth-Hedin J. Phys. C5,1629(1972)
C
         txch='B-H'
         xccp=0.0504D0
         xccf=0.0254D0
         xcrp=30.D0
         xcrf=75.D0
         fth=4.D0/3.D0
         aa=0.5d0**oth
         bb=1.d0-aa
      elseif(ixc.eq.2) then
C
C        Slater X-Alpha
C
         txch='X-A'
         xalpha=6.d0*exchf*(3.d0/fourpi)**oth
         write(6,102) exchf
      elseif(ixc.eq.3) then
C
C        Barth-Hedin-Janak Phys. Rev. B12,1257(1975)
C
         txch='BHJ'
         xccp=0.045d0
         xccf=0.0225d0
         xcrp=21.d0
         xcrf=53.d0
         fth=4.d0/3.d0
         aa=0.5d0**oth
         bb=1.d0-aa
      elseif(ixc.eq.4) then
C
C        Vosko-Wilk-Nusair Can. J. Phys. 58,1200(1980)
C
         txch='VWN'
         fth=4.d0/3.d0
         aa=2.d0**fth-2.d0
      elseif(ixc.eq.5) then
C
C        Langreth-Perdew-Mehl Phys. Rev. B28,1809(1983)
C
         txch='LPM'
         if(ns.eq.2) then
            write(6,101) ixc,txch
            stop
         endif
         xccp=0.0504d0
         xccf=0.0254d0
         xcrp=30.d0
         xcrf=75.d0
         fth=4.d0/3.
         aa=0.5d0**oth
         bb=1.d0-aa
         smallf=0.15d0
         alpm=4.287115843d-03
         blpm=1.745415107d0*smallf
      elseif(ixc.eq.6) then
C
C        Wigner exchange
C
         txch='WXC'
         if(ns.eq.2) then
            write(6,101) ixc,txch
            stop
         endif
         aw=0.916d0*4.d0/3.d0
         bw=0.88d0*4.d0/3.d0
         cw=0.88d0*7.8d0/3.d0
      elseif(ixc.eq.7) then
C
C        Ceperley-Alder. Parametrization by Perdew and Zunger
C
         txch='C-A'
c
         acau=1.0529d0
         bcau=0.3334d0
         ccau=7.d0*acau/6.d0
         dcau=4.d0*bcau/3.d0
         fca=4.d0/3.d0
         ocau=0.096d0
         pcau=0.0622d0
         qcau=0.0232d0
         rcau=0.004d0
         scau=ocau+pcau/3.d0
         tcau=(2.d0*qcau+rcau)/3.d0
         ucau=2.d0*rcau/3.d0
         acap=1.3981d0
         bcap=0.2611d0
         ccap=7.d0*acap/6.d0
         dcap=4.d0*bcap/3.d0
         ocap=0.0538d0
         pcap=0.0311d0
         qcap=0.0096d0
         rcap=0.0014d0
         scap=ocap+pcap/3.d0
         tcap=(2.d0*qcap+rcap)/3.d0
         ucap=2.d0*rcap/3.d0
      elseif(ixc.eq.8) then
C
C        Perdew-Wang.Local correlation by CA. 
C
         txch='P-W'
         acau=1.0529d0
         bcau=0.3334d0
         ccau=7.d0*acau/6.d0
         dcau=4.d0*bcau/3.d0
         fca=4.d0/3.d0
         ocau=0.096d0
         pcau=0.0622d0
         qcau=0.0232d0
         rcau=0.004d0
         scau=ocau+pcau/3.d0
         tcau=(2.d0*qcau+rcau)/3.d0
         ucau=2.d0*rcau/3.d0
         acap=1.3981d0
         bcap=0.2611d0
         ccap=7.d0*acap/6.d0
         dcap=4.d0*bcap/3.d0
         ocap=0.0538d0
         pcap=0.0311d0
         qcap=0.0096d0
         rcap=0.0014d0
         scap=ocap+pcap/3.d0
         tcap=(2.d0*qcap+rcap)/3.d0
         ucap=2.d0*rcap/3.d0
C
         apw=1.296d0
         bpw=14.d0
         cpw=0.2d0
         pmpw=1.d0/15.d0
C
         cp1=0.001667d0
         cp2=0.002568d0
         alfa=0.023266d0
         beta=7.389d-6
         gamma=8.723d0
         delta=0.472d0
         beta1=7.389d-2
         ftild=0.11d0
         cinf=cp1+cp2
         fif=1.745d0*ftild*cinf
         ccd0=alfa-cp2*gamma
         ccd1=2.d0*(beta-cp2*delta)
         ccd2=beta*gamma-alfa*delta-3.d0*cp2*beta1
         ccd3=-2.d0*alfa*beta1
         ccd4=-beta*beta1
         pp1=11.d0/3.d0
         pp2=7.d0/6.d0
      elseif(ixc.eq.9) then
C
C        Becke-Perdew. Local correlation by VWN.
C
         txch='B-P'
         fca=4.d0/3.d0
         fth=fcA
         aa=2.d0**fth-2.d0
C
         betab=0.0042d0
c
         cp1=0.001667d0
         cp2=0.002568d0
         alfa=0.023266d0
         beta=7.389d-6
         gamma=8.723d0
         delta=0.472d0
         beta1=7.389d-2
         ftild=0.11d0
         cinf=cp1+cp2
         fif=1.745d0*ftild*cinf
         ccd0=alfa-cp2*gamma
         ccd1=2.d0*(beta-cp2*delta)
         ccd2=beta*gamma-alfa*delta-3.d0*cp2*beta1
         ccd3=-2.d0*alfa*beta1
         ccd4=-beta*beta1
         pp1=11.d0/3.d0
         pp2=7.d0/6.d0
C
c      ELSE
c         WRITE(6,103) IXC
c         STOP
      endif
c
      write(6,'('' Exchange-correlation type:'',a3/)') txch
c
      return
      end
