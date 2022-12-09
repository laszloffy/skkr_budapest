c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine initzero1(kmymax,nimp,ne,
     > dosimp,doshimp,dosehimp,
     > dosheimp,dostehimp,dostheimp,
     > qvpimp,qvhpimp,qvimp,qvhimp,qvehimp,qvheimp,
     > qvtehimp,qvtheimp,
     > qmomimpha,qmomimpeha,qmomimphea,qmomimpteha,qmomimpthea)
c
      include '../param.h'
c
      real*8 dosimp(kmymaxp,nimp,me),doshimp(kmymaxp,nimp,me),
     > dosehimp(kmymaxp,nimp,me),dosheimp(kmymaxp,nimp,me),
     > dostehimp(kmymaxp,nimp,me),dostheimp(kmymaxp,nimp,me),
     > qvpimp(kmymaxp,nimp),qvhpimp(kmymaxp,nimp),qvimp(nimp),
     > qvhimp(nimp),qvehimp(nimp),qvheimp(nimp),qvtehimp(nimp),
     > qvtheimp(nimp)
c
      complex*16 qmomimpha(lmsup,nimp),qmomimpeha(nimp),
     > qmomimphea(nimp),qmomimpteha(nimp),qmomimpthea(nimp)
c
      dosimp=(0.0d0)
      doshimp=(0.0d0)
      dosehimp=(0.0d0)
      dosheimp=(0.0d0)
      dostehimp=(0.0d0)
      dostheimp=(0.0d0)
      qvpimp=(0.0d0)
      qvhpimp=(0.0d0)
      qvimp=(0.0d0)
      qvhimp=(0.0d0)
      qvehimp=(0.0d0)
      qvheimp=(0.0d0)
      qvtehimp=(0.0d0)
      qvtheimp=(0.0d0)
      qmomimpha=(0.0d0,0.0d0)
      qmomimpeha=(0.0d0,0.0d0)
      qmomimphea=(0.0d0,0.0d0)
      qmomimpteha=(0.0d0,0.0d0)
      qmomimpthea=(0.0d0,0.0d0)
c
      return
      end
