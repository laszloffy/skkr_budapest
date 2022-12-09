c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine derivt(nintfc,kmymax,ptminv,
     > dmat,dmatp,
     > ddph,ddphp,d2dph,d2dphp,
     > ddth,ddthp,d2dth,d2dthp,
     > d2dthph,d2dthphp,
     > dtmph, dtmth, d2tmph, d2tmth, d2tmthph)
c
c Input: 
c  ptminv           : physical t^{-1} matrix in local frame
c  ddph,ddphp       : dR/dphi , dR'/dphi
c  d2dph, d2dphp    : d^2R/dphi^2  , d^2R'/dphi^2
c  ddth, ddthp      : dR/dtheta  , dR'/dtheta
c  d2dth, d2dthp    : d^2R/dtheta^2  , d^2R'/dtheta^2
c  d2dthph, d2dthphp: d^2R/(dtheta dphi)  , d^2R'/(dtheta dphi)
c
c  dmat, dmatp     : R(theta,phi) in the new frame
c ---------------------------------------------------------------
c On output:
c  dtmph           : dt/dphi
c  dtmth           : dt/dtheta
c  d2tmph          : dt^2/dphi^2
c  d2tmth          : dt^2/dtheta^2
c  d2tmthph        : dt^2/(dtheta dphi)
c
c  dR/dthetha = dR(pi/2,phi)/dtheta R(theta-pi/2,0)      
c  dR/dphi    = dR(pi/2,phi)/dphi R(theta-pi/2,0)      
c
c     implicit none
      include '../param.h'
c      
      INTEGER*4 NINTFC ,KMYMAX ,LI 
c      
      complex*16 ddph(kmymaxp,kmymaxp,mimp)
      complex*16 ddphp(kmymaxp,kmymaxp,mimp)
      complex*16 ddth(kmymaxp,kmymaxp,mimp)
      complex*16 ddthp(kmymaxp,kmymaxp,mimp)
      complex*16 d2dph(kmymaxp,kmymaxp,mimp)
      complex*16 d2dphp(kmymaxp,kmymaxp,mimp)
      complex*16 d2dth(kmymaxp,kmymaxp,mimp)
      complex*16 d2dthp(kmymaxp,kmymaxp,mimp)
      complex*16 d2dthph(kmymaxp,kmymaxp,mimp)
      complex*16 d2dthphp(kmymaxp,kmymaxp,mimp)
      complex*16 dmat(kmymaxp,kmymaxp,mimp)
      complex*16 dmatp(kmymaxp,kmymaxp,mimp)
c
      complex*16 ptminv(kmymaxp,kmymaxp,mimp)
      complex*16 tmwrk(kmymaxp,kmymaxp),wrk1(kmymaxp,kmymaxp)
c      
      complex*16 dtmph(kmymaxp,kmymaxp,mimp)
      complex*16 dtmth(kmymaxp,kmymaxp,mimp)
      complex*16 d2tmph(kmymaxp,kmymaxp,mimp)
      complex*16 d2tmth(kmymaxp,kmymaxp,mimp)
      complex*16 d2tmthph(kmymaxp,kmymaxp,mimp)
c
c
      do li=1,nintfc
c       
       call repl(tmwrk,ptminv(1,1,li),kmymax,kmymaxp)
c
c Derivatives of the inverse of the t matrix with respect of phi and theta
c  dtmph contains R'ph(phys t^{-1})Rp + R(phys t^{-1})Rp'ph
c  dtmth contains R'th(phys t^{-1})Rp + R(phys t^{-1})Rp'th
c
        call repl(wrk1,tmwrk,kmymax,kmymaxp)
        call repl(dtmth(1,1,li),wrk1,kmymax,kmymaxp)
        call tripmt(ddth(1,1,li),wrk1,dmatp(1,1,li),
     >              kmymax,kmymax,kmymaxp)
        call tripmt(dmat(1,1,li),dtmth(1,1,li),ddthp(1,1,li),
     >              kmymax,kmymax,kmymaxp)
        call addmat(dtmth(1,1,li),wrk1,kmymax,kmymaxp)
c        
        call repl(wrk1,tmwrk,kmymax,kmymaxp)
c
        call repl(dtmph(1,1,li),wrk1,kmymax,kmymaxp)
        call tripmt(ddph(1,1,li),wrk1,dmatp(1,1,li),
     >              kmymax,kmymax,kmymaxp)
        call tripmt(dmat(1,1,li),dtmph(1,1,li),ddphp(1,1,li),
     >              kmymax,kmymax,kmymaxp)
        call addmat(dtmph(1,1,li),wrk1,kmymax,kmymaxp)
c        
        call repl(wrk1,tmwrk,kmymax,kmymaxp)
c
c Second derivatives of the inverse of the t matrix
c with respect of phi and theta
c  d2tmph contains R''ph(phys t^{-1})Rp + R(phys t^{-1})Rp''ph +
c  + 2*R'ph(phys t^{-1})Rp'ph
c  d2tmth contains R''th(phys t^{-1})Rp + R(phys t^{-1})Rp''th +
c  + 2*R'th(phys t^{-1})Rp'th
c  d2tmthph contains R''thph(phys t^{-1})Rp + R(phys t^{-1})Rp''thph +
c  + R'th(phys t^{-1})Rp'ph + R'ph(phys t^{-1})Rp'th
c
        call repl(d2tmph(1,1,li),wrk1,kmymax,kmymaxp)
        call tripmt(d2dph(1,1,li),wrk1,dmatp(1,1,li),
     >              kmymax,kmymax,kmymaxp)
        call tripmt(dmat(1,1,li),d2tmph(1,1,li),d2dphp(1,1,li),
     >              kmymax,kmymax,kmymaxp)
        call addmat(d2tmph(1,1,li),wrk1,kmymax,kmymaxp)
        call repl(wrk1,tmwrk,kmymax,kmymaxp)
        call tripmt(ddph(1,1,li),wrk1,ddphp(1,1,li),
     >              kmymax,kmymax,kmymaxp)
        call addmata(d2tmph(1,1,li),wrk1,2.d0,kmymax,kmymaxp)
c
c
        call repl(wrk1,tmwrk,kmymax,kmymaxp)
c
        call repl(d2tmth(1,1,li),wrk1,kmymax,kmymaxp)
        call tripmt(d2dth(1,1,li),wrk1,dmatp(1,1,li),
     >              kmymax,kmymax,kmymaxp)
        call tripmt(dmat(1,1,li),d2tmth(1,1,li),d2dthp(1,1,li),
     >              kmymax,kmymax,kmymaxp)
        call addmat(d2tmth(1,1,li),wrk1,kmymax,kmymaxp)
        call repl(wrk1,tmwrk,kmymax,kmymaxp)
        call tripmt(ddth(1,1,li),wrk1,ddthp(1,1,li),
     >              kmymax,kmymax,kmymaxp)
        call addmata(d2tmth(1,1,li),wrk1,2.d0,kmymax,kmymaxp)
c
        call repl(wrk1,tmwrk,kmymax,kmymaxp)
c
        call repl(d2tmthph(1,1,li),wrk1,kmymax,kmymaxp)
        call tripmt(d2dthph(1,1,li),wrk1,dmatp(1,1,li),
     >              kmymax,kmymax,kmymaxp)
        call tripmt(dmat(1,1,li),d2tmthph(1,1,li),d2dthphp(1,1,li),
     >              kmymax,kmymax,kmymaxp)
        call addmat(d2tmthph(1,1,li),wrk1,kmymax,kmymaxp)
        call repl(wrk1,tmwrk,kmymax,kmymaxp)
        call tripmt(ddth(1,1,li),wrk1,ddphp(1,1,li),
     >              kmymax,kmymax,kmymaxp)
        call addmat(d2tmthph(1,1,li),wrk1,kmymax,kmymaxp)
        call repl(wrk1,tmwrk,kmymax,kmymaxp)
        call tripmt(ddph(1,1,li),wrk1,ddthp(1,1,li),
     >              kmymax,kmymax,kmymaxp)
        call addmat(d2tmthph(1,1,li),wrk1,kmymax,kmymaxp)
c
      end do
      return
      end
