c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      function cerf (z) 
!     =============                                                     
!                                                                       
!----- evaluate the complex error function using either                 
!      1 power series                                                   
!      2 continued fraction                                             
!      3 hermite integration                                            
! ----------------------------------------------------------            
      implicit real*8 (a-h,o-z) 
      parameter (p=0.3275911, a1=0.254829592, a2=-0.284496736) 
      parameter (a3=1.421413741, a4=-1.453152027, a5=1.061405429) 
!                                                                       
      complex *16 z,zz,sum,zzs,xzzs,h1,h2,h3,u1,u2,u3,term1,   
     &            term2,cerf,cian,phase                                 
      dimension expn2(24) 
!                                                                       
      data ne /0/ 
!                                                                       
      save expn2,ne 
!                                                                       
      emach=1.0d-10 
      pi=4.0d0*atan(1.0d0) 
!                                                                       
      eps=5.0*emach 
      api=1.0/pi 
      absz=abs(z) 
c     write(6,*)'z,absz',z,absz
      cerf=1.0 
      if (absz.lt.emach) return 
!                                                                       
!----- the argument is translated to the first quadrant from            
!      the nn'th quadrant, before the method for the function           
!      evaluation is chosen                                             
!                                                                       
      x=dreal(z) 
      y=dimag(z) 
      ax=abs(x) 
      ay=abs(y) 
      zz=dcmplx(ax,ay) 
      zzs=zz*zz 
      nn=1 
      if (x.ne.ax) nn=2 
      if (y.ne.ay) nn=5-nn 
      if (nn.eq.3 .or. nn.eq.4) then
             if (abs(zzs).lt.300.d0) then
                xzzs=exp(-zzs)
c               write(6,*)'IN CERF: zzs,xzzs',zzs,xzzs
             else
                xzzs=(0.d0,0.d0)
c               write(6,*)'IN CERF:zzs,xzzs set 0.d0',zzs
             end if
      end if 
!     if (ax.gt.3.9.or.ay.gt.3.0) go to 80                             
      if (ay.ge.1.0.or.absz.ge.4.0) go to 40 
!                                                                       
!----- power series: see abramowitz and stegun's handbook of            
!      mathematical functns, p297                                       
!                                                                       
   10 q=1.0 
      xzzs=exp(-zzs) 
      md=1 
      factn=-1.0 
      factd=1.0 
      term1=zz 
      sum=zz 
!                                                                       
!----- sum up power series in groups of 5 terms                         
!                                                                       
   20 do 30 n=1,5 
      factn=factn+2.0 
      factd=factd+2.0 
      fact=factn/(q*factd) 
      term1=fact*zzs*term1 
      sum=sum+term1 
   30 q=q+1.0 
!                                                                       
!----- test convergence of power series                                 
!                                                                       
      abterm=abs(term1) 
      if (abterm.ge.eps) then 
      if (q.lt.1000.0) go to 20 
      write (6,*) ' error: in cerf - series not converged' 
      stop 
      endif 
      fact=2.0*sqrt(api) 
      sum=fact*(0.0,1.0)*sum 
      cerf=xzzs+xzzs*sum 
      go to 90 
!                                                                       
   40 continue 
!                                                                       
!----- approximation to w(z) based on abramowitz and stegun equation    
!      7.1.29                                                           
!      summing 24 terms, if convergence not achieved then try series    
!                                                                       
      md=2 
      if (ne.eq.0) then 
!                                                                       
!----- generate the exponential prefactors 2*exp(-n*n/4)/pi             
!                                                                       
      do 50 n=1,24 
      an=dfloat(n)/2.0 
      arg=-an*an 
      expn2(n)=2.0*api*exp(arg) 
   50 continue 
      ne=1 
!                                                                       
      endif 
!                                                                       
      emx2=exp(-ax*ax) 
      epx=exp(ax) 
      ex=epx 
      twoy=2.0*ay 
      phase=cos(ax*twoy)-(0.0,1.0)*sin(ax*twoy) 
      t=1.0/(1.0+p*ay) 
!                                                                       
!----- abramowitz and stegun equation 7.1.26 for the error functn of    
!      a real variable                                                  
!                                                                       
!      cerf=phase*t*(a1+t*(a2+t*(a3+t*(a4+t*a5))))-api*(phase-1.0)/twoy 
      cerf=phase*funerr(t)-api*(phase-1.0)/twoy 
!                                                                       
      sum=cerf 
      do 60 n=1,24 
      an=dfloat(n) 
      cian=(0.0,1.0)*an 
      denom=an*an+twoy*twoy 
      term1=expn2(n)*(twoy*phase+0.5*((cian-twoy)/ex-(cian+twoy)*ex))   
     & /denom                                                           
      sum=sum-term1 
      if (abs(term1/sum).lt.eps) go to 70 
      ex=ex*epx 
   60 continue 
!                                                                       
!----- not converged try series                                         
!                                                                       
      write (6,*) ' warning: w(z) approximate expansion not converged' 
      go to 10 
!                                                                       
   70 continue 
      cerf=sum*emx2 
!                                                                       
      go to 90 
!                                                                       
!----- asymptotic series: see abramowitz and stegun, p328               
!                                                                       
   80 if (ax.gt.6.or.ay.gt.6) then 
      cerf=0.5124242/(zzs-0.2752551)+0.05176536/(zzs-2.724745) 
      else 
         cerf=(0.4613135/(zzs-0.1901635)+0.09999216/(zzs-1.7844927)+ 
     &                0.002883894/(zzs-5.5253437)) 
      endif 
      cerf=(0.0,1.0)*zz*cerf 
      md=3 
!                                                                       
!----- symmetry relations are now used to transform the functn          
!      back to quadrant nn                                              
!                                                                       
   90 go to (130,110,120,100), nn 
  100 cerf=2.0*xzzs-cerf 
  110 cerf=dconjg(cerf) 
      return 
  120 cerf=2.0*xzzs-cerf 
  130 return 
      END                                           
                                                                        
                                                                        
                                                                        
      function csevl(x,cs,n) 
      implicit real*8 (a-h,o-z) 
!***begin prologue  csevl                                               
!***date written   770401   (yymmdd)                                    
!***revision date  861211   (yymmdd)                                    
!***category no.  c3a2                                                  
!***keywords  library=slatec(fnlib),                                    
!             type=single precision(csevl-s dcsevl-d),chebyshev,        
!             special functions                                         
!***author  fullerton, w., (lanl)                                       
!***purpose  evaluate the n-term chebyshev series cs at x.              
!***description                                                         
!                                                                       
! evaluate the n-term chebyshev series cs at x.  adapted from           
! r. broucke, algorithm 446, c.a.c.m., 16, 254 (1973). also see fox     
! and parker, chebyshev polynomials in numerical analysis, oxford press,
! page 56.                                                              
!                                                                       
!       input arguments --                                              
! x    value at which the series is to be evaluated.                    
! cs   array of n terms of a chebyshev series.  in eval-                
!      uating cs, only half the first coefficient is summed.            
! n    number of terms in array cs.                                     
!***references  (none)                                                  
!***routines called  xerror                                             
!***end prologue  csevl                                                 
!                                                                       
       dimension cs(*) 
       save 
!***first executable statement  csevl                                   
!1234567                                                                
       if(n.lt.1) stop 'csevl   number of terms le 0' 
       if(n.gt.1000) stop 'csevl   number of terms gt 1000' 
       if (x.lt. -1.0 .or. x.gt. 1.0)                          
     & stop 'csevl   x outside (-1,+1)'                                 
!                                                                       
       b1=0. 
       b0=0. 
       twox=2.*x 
       do 10 i=1,n 
       b2=b1 
       b1=b0 
       ni=n+1-i 
       b0=twox*b1-b2+cs(ni) 
   10  continue 
!                                                                       
       csevl = 0.5 * (b0-b2) 
!                                                                       
       return 
       END                                           
      function erfce (x,emach) 
!     ==============                                                    
!                                                                       
!----- calculates exp(x*x)*erfc(x), modified from lanl clams routine    
!      erfc(x) by jmm june 90                                           
!                                                                       
      implicit real*8(a-h,o-z) 
      dimension erfcs(13), erfccs(24), erc2cs(23) 
      save 
      data erfcs( 1) /   -.049046121234691808e0 / 
      data erfcs( 2) /   -.14226120510371364e0 / 
      data erfcs( 3) /    .010035582187599796e0 / 
      data erfcs( 4) /   -.000576876469976748e0 / 
      data erfcs( 5) /    .000027419931252196e0 / 
      data erfcs( 6) /   -.000001104317550734e0 / 
      data erfcs( 7) /    .000000038488755420e0 / 
      data erfcs( 8) /   -.000000001180858253e0 / 
      data erfcs( 9) /    .000000000032334215e0 / 
      data erfcs(10) /   -.000000000000799101e0 / 
      data erfcs(11) /    .000000000000017990e0 / 
      data erfcs(12) /   -.000000000000000371e0 / 
      data erfcs(13) /    .000000000000000007e0 / 
      data erc2cs( 1) /   -.069601346602309501e0 / 
      data erc2cs( 2) /   -.041101339362620893e0 / 
      data erc2cs( 3) /    .003914495866689626e0 / 
      data erc2cs( 4) /   -.000490639565054897e0 / 
      data erc2cs( 5) /    .000071574790013770e0 / 
      data erc2cs( 6) /   -.000011530716341312e0 / 
      data erc2cs( 7) /    .000001994670590201e0 / 
      data erc2cs( 8) /   -.000000364266647159e0 / 
      data erc2cs( 9) /    .000000069443726100e0 / 
      data erc2cs(10) /   -.000000013712209021e0 / 
      data erc2cs(11) /    .000000002788389661e0 / 
      data erc2cs(12) /   -.000000000581416472e0 / 
      data erc2cs(13) /    .000000000123892049e0 / 
      data erc2cs(14) /   -.000000000026906391e0 / 
      data erc2cs(15) /    .000000000005942614e0 / 
      data erc2cs(16) /   -.000000000001332386e0 / 
      data erc2cs(17) /    .000000000000302804e0 / 
      data erc2cs(18) /   -.000000000000069666e0 / 
      data erc2cs(19) /    .000000000000016208e0 / 
      data erc2cs(20) /   -.000000000000003809e0 / 
      data erc2cs(21) /    .000000000000000904e0 / 
      data erc2cs(22) /   -.000000000000000216e0 / 
      data erc2cs(23) /    .000000000000000052e0 / 
      data erfccs( 1) /    .0715179310202925e0 / 
      data erfccs( 2) /   -.026532434337606719e0 / 
      data erfccs( 3) /    .001711153977920853e0 / 
      data erfccs( 4) /   -.000163751663458512e0 / 
      data erfccs( 5) /    .000019871293500549e0 / 
      data erfccs( 6) /   -.000002843712412769e0 / 
      data erfccs( 7) /    .000000460616130901e0 / 
      data erfccs( 8) /   -.000000082277530261e0 / 
      data erfccs( 9) /    .000000015921418724e0 / 
      data erfccs(10) /   -.000000003295071356e0 / 
      data erfccs(11) /    .000000000722343973e0 / 
      data erfccs(12) /   -.000000000166485584e0 / 
      data erfccs(13) /    .000000000040103931e0 / 
      data erfccs(14) /   -.000000000010048164e0 / 
      data erfccs(15) /    .000000000002608272e0 / 
      data erfccs(16) /   -.000000000000699105e0 / 
      data erfccs(17) /    .000000000000192946e0 / 
      data erfccs(18) /   -.000000000000054704e0 / 
      data erfccs(19) /    .000000000000015901e0 / 
      data erfccs(20) /   -.000000000000004729e0 / 
      data erfccs(21) /    .000000000000001432e0 / 
      data erfccs(22) /   -.000000000000000439e0 / 
      data erfccs(23) /    .000000000000000138e0 / 
      data erfccs(24) /   -.000000000000000048e0 / 
      data sqrtpi /1.7724538509055160e0/ 
      data nterf, nterfc, nterc2, xsml, xmax, sqeps /3*0, 3*0./ 
!                                                                       
      if (x.lt.0) stop ' error in erce x<0' 
      if (nterf.ne.0) go to 10 
      eta = 0.1*emach 
      xmax=1.0/emach 
      nterf = inits (erfcs, 13, eta) 
      nterfc = inits (erfccs, 24, eta) 
      nterc2 = inits (erc2cs, 23, eta) 
!                                                                       
      sqeps = sqrt (2.0*emach) 
!                                                                       
   10 if (x.gt.emach) go to 20 
!                                                                       
!----- erfce(x) = 1.0 for x .lt. emach                                  
!                                                                       
      erfce = 1.0 
      return 
!                                                                       
   20 if (x.gt.xmax) go to 40 
      y = abs(x) 
      if (y.gt.1.0) go to 30 
!                                                                       
!----- erfce(x) = exp(x*x)*(1.0 - erf(x)) for  x .le. 1.                
!                                                                       
      if (y.lt.sqeps) erfce =   1.0 - 2.0*x/sqrtpi 
      if (y.ge.sqeps) erfce = ( 1.0 -                            
     &  x*(1.0 + csevl (2.*x*x-1., erfcs, nterf) ))*exp(x*x)            
      return 
!                                                                       
!----- erfc(x) = 1.0 - erf(x) for 1. .lt. abs(x) .le. xmax              
!                                                                       
   30 y = y*y 
      if (y.le.4.) erfce = 1.0/abs(x) * (0.5 + csevl ((8./y-5.)/3.,
     &  erc2cs, nterc2) )                                               
      if (y.gt.4.) erfce = 1.0/abs(x) * (0.5 + csevl (8./y-1.,    
     &  erfccs, nterfc) )                                               
      return 
!                                                                       
   40 erfce = 0.0 
      return 
!                                                                       
      END                                           
!                                                                       
      function funerr(t) 
!     ==========                                                        
!                                                                       
!------ calculate f(t)                                                  
!                                                                       
      implicit real*8 (a-h,o-z) 
      parameter (p=0.3275911d+00) 
!                                                                       
      x=(1.0/t-1.0)/p 
      emach=1.0e-12 
      funerr=erfce(x,emach) 
      return 
      END                                           
!                                                                       
      function inits(os,nos,eta) 
      implicit real*8 (a-h,o-z) 
!***begin prologue  inits                                               
!***date written   770401   (yymmdd)                                    
!***revision date  861211   (yymmdd)                                    
!***category no.  c3a2                                                  
!***keywords  library=slatec(fnlib),                                    
!             type=single precision(inits-s initds-d),chebyshev,        
!             initialize,orthogonal polynomial,orthogonal series,series,
!             special functions                                         
!***author  fullerton, w., (lanl)                                       
!***purpose  initializes an orthogonal series so that it defines the    
!            number of terms to carry in the series to meet a specified 
!            error.                                                     
!***description                                                         
!                                                                       
! initialize the orthogonal series so that inits is the number of terms 
! needed to insure the error is no larger than eta.  ordinarily, eta    
! will be chosen to be one-tenth machine precision.                     
!                                                                       
!             input arguments --                                        
! os     array of nos coefficients in an orthogonal series.             
! nos    number of coefficients in os.                                  
! eta    requested accuracy of series.                                  
!***references  (none)                                                  
!***routines called  xerror                                             
!***end prologue  inits                                                 
      dimension os(nos) 
      save 
!***first executable statement  inits                                   
      if (nos.lt.1) stop 'inits   number of coefficients  < 1' 
!                                                                       
      err = 0. 
      do 10 ii=1,nos 
        i = nos + 1 - ii 
        err = err + abs(os(i)) 
        if (err.gt.eta) go to 20 
   10 continue 
!                                                                       
   20 if (i.eq.nos) write(6,*) 'inits   eta may be too small' 
      inits = i 
!                                                                       
      return 
      END                                           
