c*********************************************************
      integer lmaxp
      integer l1maxp
      integer lmmaxp
      integer kmaxp
      integer kmymaxp
      integer lsup
      integer lmsup
      integer lshape
      integer lmshape
      integer minprc
      integer mprc
      integer mprc1
      integer mintfc
      integer mextra
      integer mtotal
      integer mtotal1
      integer mnum
      integer mbulk
      integer mdimnr
      integer mdimr
      integer mkpar
      integer me
      integer mekr
      integer mdir
      integer mrec
      integer mr
      integer mk
      integer mdr
      integer nrad
      integer iorb
      integer melem
      integer melem1
      integer nuzp
c
      integer mimp
      integer mpair
      integer net
c
      integer mll1, mll2
c related to angular momentum space
      parameter (lmaxp=2)
      parameter (l1maxp=lmaxp+1)
      parameter (lmmaxp=(lmaxp+1)**2 )
      parameter (kmaxp=2*lmaxp+1 )
      parameter (kmymaxp=2*lmmaxp)
      parameter (lsup=2*lmaxp)
      parameter (lmsup=(lsup+1)**2)
      parameter (lshape=4*lmaxp)
      parameter (lmshape=(lshape+1)**2)
c         to geometry
      parameter (minprc=3)
      parameter (mprc=8,mprc1=8) 
      parameter (mintfc=mprc*minprc)
      parameter (mextra=1)
      parameter (mtotal=(mprc+2+mextra*2)*minprc)
      parameter (mtotal1=(mprc1+2+mextra*2)*minprc)
      parameter (mnum=mintfc+2*minprc)            ! gstore => ipl
      parameter (mbulk=4)                         ! madelung
c         to both
      parameter (mdimnr=minprc*lmmaxp)
      parameter (mdimr=minprc*kmymaxp)
c         to BZ-integral
      parameter (mkpar=50000)
c         to energy mesh
      parameter (me=16)
c         to a combination of e- and k-mesh 
c         (see in cpacoord)
      parameter (mekr=2200)      
c     parameter (mekr=1000)      
c         to Ewald summations
      parameter (mdir=20000,mrec=20000) ! adviced if too small
      parameter (mr=mdir,mk=mrec)       ! effective values for genlatt2d
      parameter (mdr=mrec)              ! vecsrt: must be max(mdir,mrec)
c         to others
      parameter (nrad=1100)              ! effective value for rcore
      parameter (iorb=30)
      parameter (melem=12)       
      parameter (melem1=8)       
      parameter (nuzp=2)
c
c - cluster
      parameter (net=16)
      parameter (mimp=100)
      parameter (mpair=mimp*mimp)
c
      parameter (mll1=50, mll2 = 100)
