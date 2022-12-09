c*********************************************************
c Declared parameters
c-------------------------------------
      integer iorb
      integer kmaxp
      integer kmymaxp
      integer bogomaxp
      integer dbogomaxp
      integer l1maxp
      integer lmaxp
      integer lmmaxp
      integer lmshape
      integer lmsup
      integer lshape
      integer lsup
      integer mbulk
      integer mdimnr
      integer mdimr
      integer mdimbogo
      integer mdimdbogo
      integer mdir
      integer mdr
      integer me
      integer mekr0
      integer me1
      integer melem
      integer melem1
      integer mextra
      integer minprc
      integer mintfc
      integer mk
      integer mkpar
      integer mkpar1
      integer mnum
      integer mprc
      integer mprc1
      integer mr
      integer mrec
      integer mtotal
      integer mtotal1
      integer nrad
      integer nuzp
c
      integer mimp
      integer mpair
      integer net
c whether to use local coordinate system
      logical localmode
      parameter (localmode=.false.) !leave it as it is, only globalmode is implemented
c related to angular momentum space
      parameter (lmaxp=2)  ! must be the same as lmax in input_rsp.in
      parameter (l1maxp=lmaxp+1)
      parameter (lmmaxp=(lmaxp+1)**2 )
      parameter (kmaxp=2*lmaxp+1 )
      parameter (kmymaxp=2*lmmaxp)
      parameter (bogomaxp=2*lmmaxp)   
      parameter (dbogomaxp=2*kmymaxp) 
      parameter (lsup=2*lmaxp)
      parameter (lmsup=(lsup+1)**2)
      parameter (lshape=4*lmaxp)
      parameter (lmshape=(lshape+1)**2)
c         to geometry
      parameter (minprc=4)
      parameter (mprc=4,mprc1=mprc) 
      parameter (mintfc=mprc*minprc)
      parameter (mextra=1)
      parameter (mtotal=(mprc+2+mextra*2)*minprc)
      parameter (mtotal1=(mprc1+2+mextra*2)*minprc)
      parameter (mnum=mintfc+2*minprc)            ! gstore => ipl
      parameter (mbulk=minprc)                         ! madelung
c         to both
      parameter (mdimnr=minprc*lmmaxp)
      parameter (mdimr=minprc*kmymaxp)
      parameter (mdimbogo=minprc*bogomaxp)  ! bogo
      parameter (mdimdbogo=minprc*dbogomaxp)  ! dirac-bogo
c         to BZ-integral
      parameter (mkpar=50000)
      parameter (mkpar1=136)
c         to energy mesh
      parameter (me=501)
      parameter (me1=20)
c         both to BZ integral and energy mesh
      parameter (mekr0=1000)
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
c     parameter (net=16)
c     parameter (mimp=100)
c     parameter (mpair=mimp*mimp)
c
      integer kapdex(72),ldex(72),jdex(72),mydex(72),lbdex(72),sdex(72)
      data kapdex/
     * -1,-1,
     *  1, 1,
     * -2,-2,-2,-2,
     *  2, 2, 2, 2,
     * -3,-3,-3,-3,-3,-3,
     *  3, 3, 3, 3, 3, 3,
     * -4,-4,-4,-4,-4,-4,-4,-4,
     *  4, 4, 4, 4, 4, 4, 4, 4,
     * -5,-5,-5,-5,-5,-5,-5,-5,-5,-5,
     *  5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
     * -6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6/
      data ldex/
     *  0, 0,
     *  1, 1,
     *  1, 1, 1, 1,
     *  2, 2, 2, 2,
     *  2, 2, 2, 2, 2, 2,
     *  3, 3, 3, 3, 3, 3,
     *  3, 3, 3, 3, 3, 3, 3, 3,
     *  4, 4, 4, 4, 4, 4, 4, 4,
     *  4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
     *  5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
     *  5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5/
      data lbdex/
     *  1, 1,
     *  0, 0,
     *  2, 2, 2, 2,
     *  1, 1, 1, 1,
     *  3, 3, 3, 3, 3, 3,
     *  2, 2, 2, 2, 2, 2,
     *  4, 4, 4, 4, 4, 4, 4, 4,
     *  3, 3, 3, 3, 3, 3, 3, 3,
     *  5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
     *  4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
     *  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6/
      data jdex/
     *  1, 1,
     *  1, 1,
     *  3, 3, 3, 3,
     *  3, 3, 3, 3,
     *  5, 5, 5, 5, 5, 5,
     *  5, 5, 5, 5, 5, 5,
     *  7, 7, 7, 7, 7, 7, 7, 7,
     *  7, 7, 7, 7, 7, 7, 7, 7,
     *  9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
     *  9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
     * 11,11,11,11,11,11,11,11,11,11,11,11/
      data mydex/
     * -1, 1,
     * -1, 1,
     * -3,-1, 1, 3,
     * -3,-1, 1, 3,
     * -5,-3,-1, 1, 3, 5,
     * -5,-3,-1, 1, 3, 5,
     * -7,-5,-3,-1, 1, 3, 5, 7,
     * -7,-5,-3,-1, 1, 3, 5, 7,
     * -9,-7,-5,-3,-1, 1, 3, 5, 7, 9,
     * -9,-7,-5,-3,-1, 1, 3, 5, 7, 9,
     *-11,-9,-7,-5,-3,-1, 1, 3, 5, 7, 9,11/
      data sdex/
     * -1,-1,
     *  1, 1,
     * -1,-1,-1,-1,
     *  1, 1, 1, 1,
     * -1,-1,-1,-1,-1,-1,
     *  1, 1, 1, 1, 1, 1,
     * -1,-1,-1,-1,-1,-1,-1,-1,
     *  1, 1, 1, 1, 1, 1, 1, 1,
     * -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     *  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     * -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1/
