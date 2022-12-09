c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine initzero(qvpa,qvpb,
     & qva,qvha,qveha,qvhea,qvteha,qvthea,
     & qvb,qvhb,qvehb,qvheb,qvtehb,qvtheb,
     & spin_magvpa,spin_magvpb,spin_magva,spin_magvb,
     & orb_magvpa,orb_magvpb,orb_magva,orb_magvb,
     & enba,enbb,enorba,enorbb,qmoma,qmomb,
     & ddphenba,ddthenba,d2dphenba,d2dthenba,d2dthphenba,
     & ddphenbb,ddthenbb,d2dphenbb,d2dthenbb,d2dthphenbb,
     & rhova,rhoveha,rhovhea,rhovb,rhospa,rhospb,rhodspa,rhodspb,
     & rhomaga,rhomagb,enbifc,qvifc,omifc)
c
      implicit real*8 (a-h,o-z) 
      include '../param.h'
c
      dimension qva(mintfc),qvha(mintfc),qvhb(mintfc),qvb(mintfc)
      dimension qvpa(kmymaxp,mintfc),qvpb(kmymaxp,mintfc)     
      dimension spin_magvpa(kmymaxp,mintfc,3)
      dimension spin_magvpb(kmymaxp,mintfc,3)
      dimension spin_magva(mintfc,3),spin_magvb(mintfc,3)
      dimension orb_magvpa(kmymaxp,mintfc,3)
      dimension orb_magvpb(kmymaxp,mintfc,3)
      dimension orb_magva(mintfc,3),orb_magvb(mintfc,3)  
      dimension enba(mintfc),enbb(mintfc) 
      dimension ddphenba(mintfc),ddphenbb(mintfc)
      dimension ddthenba(mintfc),ddthenbb(mintfc)
      dimension d2dphenba(mintfc),d2dphenbb(mintfc)
      dimension d2dthenba(mintfc),d2dthenbb(mintfc)
      dimension d2dthphenba(mintfc),d2dthphenbb(mintfc)
      dimension enorba(mintfc),enorbb(mintfc)
      dimension rhova(nrad,mintfc),rhovb(nrad,mintfc)
      dimension rhospa(nrad,2,mintfc),rhospb(nrad,2,mintfc)
      dimension rhodspa(nrad,2,mintfc),rhodspb(nrad,2,mintfc)
      dimension rhomaga(nrad,mintfc),rhomagb(nrad,mintfc) 
      complex*16 qmoma(lmsup,mintfc),qmomb(lmsup,mintfc)
      complex*16 qveha(mintfc),qvhea(mintfc)
      complex*16 qvehb(mintfc),qvheb(mintfc)
      complex*16 qvteha(mintfc),qvthea(mintfc)
      complex*16 qvtehb(mintfc),qvtheb(mintfc)
      complex*16 rhoveha(nrad,mintfc)
      complex*16 rhovhea(nrad,mintfc)
c
      qvehb=(0.0d0,0.0d0)
      qveha=(0.0d0,0.0d0)
      qvheb=(0.0d0,0.0d0)
      qvhea=(0.0d0,0.0d0)
      qvtehb=(0.0d0,0.0d0)
      qvteha=(0.0d0,0.0d0)
      qvtheb=(0.0d0,0.0d0)
      qvthea=(0.0d0,0.0d0)
      rhoveha=(0.0d0,0.0d0)
      rhovhea=(0.0d0,0.0d0)
c
      call rzero(qvpa,kmymaxp*mintfc)
      call rzero(qvpb,kmymaxp*mintfc)
      call rzero(qvha,mintfc)
      call rzero(qvhb,mintfc)
      call rzero(qva,mintfc)
      call rzero(qvb,mintfc)
c
      call rzero(spin_magvpa,3*kmymaxp*mintfc)
      call rzero(spin_magvpb,3*kmymaxp*mintfc)
      call rzero(spin_magva,3*mintfc)
      call rzero(spin_magvb,3*mintfc)
c
      call rzero(orb_magvpa,3*kmymaxp*mintfc)
      call rzero(orb_magvpb,3*kmymaxp*mintfc)
      call rzero(orb_magva,3*mintfc)
      call rzero(orb_magvb,3*mintfc)
c
      call rzero(enba,mintfc)
      call rzero(enbb,mintfc)
      call rzero(ddphenba,mintfc)
      call rzero(ddthenba,mintfc)
      call rzero(d2dphenba,mintfc)
      call rzero(d2dthenba,mintfc)
      call rzero(d2dthphenba,mintfc)
      call rzero(ddphenbb,mintfc)
      call rzero(ddthenbb,mintfc)
      call rzero(d2dphenbb,mintfc)
      call rzero(d2dthenbb,mintfc)
      call rzero(d2dthphenbb,mintfc)
      call rzero(enorba,mintfc)
      call rzero(enorbb,mintfc)
c
      call czero(qmoma,lmsup*mintfc)
      call czero(qmomb,lmsup*mintfc)
c
      call rzero(rhova,nrad*mintfc)
      call rzero(rhovb,nrad*mintfc)
      call rzero(rhospa,2*nrad*mintfc)
      call rzero(rhospb,2*nrad*mintfc)
      call rzero(rhodspa,2*nrad*mintfc)
      call rzero(rhodspb,2*nrad*mintfc)
      call rzero(rhomaga,nrad*mintfc)
      call rzero(rhomagb,nrad*mintfc)         
c
      enbifc=0.d0
      qvifc=0.d0
      omifc=0.d0
c
      return
      end
