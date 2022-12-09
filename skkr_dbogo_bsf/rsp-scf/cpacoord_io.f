c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
c Store E- and k- for magnetic bulk
      subroutine cpacoord(
c     ====================
     > itscf,ie,ce,lmax,nintfc,eta,rightm,bulk,bulkgeo,wrel,
     > kset,xk,wk,nk,intbz,iek,
     > conc,itcpam,cpatol,cpatest,
     > dmata,dmatb,dmatpa,dmatpb,
     > tminvl,tminvr,tminv,tminva,tminvb,taua,taub,gtaua,gtaub)
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
      parameter(mdim=mdimr) 
c -----------------------------
      character*1  temp1
      character*2  temp2
      character*50 sgfl_file
      character*50 sgfr_file
      character*50 g2d_file
c     character*15 sgfpath
c
      integer ksgfl
      integer ksgfr
      integer kg2d
      parameter(ksgfl=81)
      parameter(ksgfr=82)
      parameter(kg2d=83)
c -----------------------------------------------
c
      logical wrel,bulk,bulkgeo,sgfx,sgfy,strc,sgfconv
      logical cpacon,concpa,cpatest,cpalay
c
      character*1 rightm
c
      dimension park(2),xk(mkpar,2),wk(mkpar)
c
      dimension invg(melem),ieqg(melem),ieqgl(melem),ieqgr(melem)
      dimension igordl(melem),igordr(melem)
c
      dimension conc(mintfc)
c
      complex*16 ce,psq,fac
c
      complex*16 tminv(kmymaxp,kmymaxp,mintfc)
      complex*16 tminva(kmymaxp,kmymaxp,mintfc)
      complex*16 tminvb(kmymaxp,kmymaxp,mintfc)
      complex*16 tminvl(kmymaxp,kmymaxp,minprc)
      complex*16 tminvr(kmymaxp,kmymaxp,minprc)
      complex*16 ttmp(kmymaxp,kmymaxp,mintfc,melem)
      complex*16 ttmpl(kmymaxp,kmymaxp,minprc,melem)
      complex*16 ttmpr(kmymaxp,kmymaxp,minprc,melem)
c
      complex*16 xx(mdim,mdim),yy(mdim,mdim),
     > gg2d(lmmaxp,lmmaxp,-minprc:minprc,minprc,0:mprc1+1),
     > gg2dbulk(lmmaxp,lmmaxp,-minprc:minprc,minprc,0:2)
c
      complex*16 xio(mdim,mdim,mkpar1,melem1),
     >           yio(mdim,mdim,mkpar1,melem1),
     > g2dio(lmmaxp,lmmaxp,-minprc:minprc,minprc,0:mprc1+1,mkpar1)
c
      complex*16 tau(kmymaxp,kmymaxp,mintfc)
      complex*16 taua(kmymaxp,kmymaxp,mintfc)
      complex*16 taub(kmymaxp,kmymaxp,mintfc)
      complex*16 gtaua(kmymaxp,kmymaxp,mintfc)
      complex*16 gtaub(kmymaxp,kmymaxp,mintfc)
      complex*16 tau1(kmymaxp,kmymaxp,mintfc,melem)
      complex*16 tau2(kmymaxp,kmymaxp)
c
      complex*16 dmata(kmymaxp,kmymaxp,mintfc)
      complex*16 dmatpa(kmymaxp,kmymaxp,mintfc)
      complex*16 dmatb(kmymaxp,kmymaxp,mintfc)
      complex*16 dmatpb(kmymaxp,kmymaxp,mintfc)
c
      complex*16 rmat(kmymaxp,kmymaxp,melem)
      complex*16 rmatp(kmymaxp,kmymaxp,melem)
c
      complex*16 alphalkkr(0:lmaxp,minprc)
      complex*16 alpharkkr(0:lmaxp,minprc)
      complex*16 alphaintkkr(0:lmaxp,mintfc)
      common/scrpar/alphalkkr,alpharkkr,alphaintkkr
c
      common/relfac/fac
c
      common/test/itest
      common/lay2d/cvec(mtotal,3),nextra,nbulkl,nbulkr,
     &             nprc,ninprc(0:mprc+1)
c
      data tol/1.0d-8/ 
      data tiny/1.0d-6/
c
      save xio,yio,g2dio
c
c ********************
c initialize constants
c ********************
c
c---> c in rydberg units:
      c=274.072d0
      psq=ce+ce*ce/(c*c)
      fac=psq/ce
c
      irel=1
c
      nl=lmax+1
      kmymax=2*nl*nl
c
      ninprcl = ninprc(0)
      ninprcr = ninprc(nprc+1)
      ndiml=ninprcl*kmymax
      ndimr=ninprcr*kmymax
c
      nprc1=nprc
      if(bulkgeo) nprc1=1
      if(nprc1.gt.mprc1) stop ' increase mprc1 !!!'
c
c ****************************************************************
c initialize irreducible representations of point group operations
c ****************************************************************
c
c     --------------------------------------
      call sphbas(rmat,rmatp,ng,invg)
c     --------------------------------------
      if(intbz.eq.0.or.intbz.eq.2) ngeff=1
      if(intbz.eq.1) ngeff=ng
      if(ngeff.gt.melem1) then
        write(6,'(/'' Parameter melem1 is too small!!!'')')
        stop
      end if
c
c
c                        ***************************
c                        * starting CPA iterations *
c                        ***************************
c
      itcpa=1
      concpa=.true.        ! controls CPA tolerance
      cpacon=.true.        ! controls no. of CPA iterations
c
      iekstart=iek
  100 continue
      iek=iekstart
c
c     ----------------------------------------
      call czero(tau,kmymaxp*kmymaxp*mintfc)
c     ----------------------------------------
c
c *******************************
c adjust t-matrices for bulk case
c *******************************
c
      if(bulk) then
c The I region is supposed to comprise one PL only !
        if(kset.ge.1) then
c Force L and R t-matrices to be identical with the
c corresponding t-matrices in I
          do li=1,ninprcl
            call repl(tminvl(1,1,li),tminv(1,1,li),kmymax,kmymaxp)
            call repl(tminvr(1,1,li),tminv(1,1,li),kmymax,kmymaxp)
          enddo
        else
c This is for only spectral-DOS calculation: force I t-matrices
c to be identical with the corresponding self-consistent L t-matrices
          do li=1,nintfc
            call repl(tminv(1,1,li),tminvl(1,1,li),kmymax,kmymaxp)
          end do
        end if
      end if
c
c *******************************************************************
c transform inverse t-matrices with respect to point group operations
c and find degeneracies
c *******************************************************************
c
c     ---------------------------------------------------
      call sorttmat(
     > lmax,nintfc,ninprcl,ninprcr,rmat,rmatp,ngeff,invg,
     > tminvl,tminvr,tminv,ttmpl,ttmpr,ttmp,
     > ieqgl,ieqgr,ieqg,igordl,igordr,0)
c     ---------------------------------------------------
c
c *************************************************
c initialize saving of structure constants and SSPO
c *************************************************
c
      if(nk.gt.mkpar1) then
        write(6,'(''nk='',i3,'' greater than'',
     >  '' mkpar1='',i3)') nk,mkpar1
        stop
      end if
c
      if(bulk) then
c    calculate g2d only for the first CPA iteration
        strc=itcpa.eq.1
      else
c    calculate g2d only for the first scf and CPA iteration
        strc=itscf.eq.1.and.itcpa.eq.1
      end if
c
      if(bulk) then
c    always calculate SSPO
        sgfx=.true.
        sgfy=.true.
      else
c    calculate SSPO only for the first scf and CPA iteration
        sgfx=itcpa.eq.1.and.itscf.eq.1
        sgfy=itcpa.eq.1.and.itscf.eq.1
      end if
c    if vacuum, recalculate R SSPO after every tenth step
      if(rightm.eq.'V')
     > sgfx=(mod(itscf,10).eq.1).and.(itcpa.eq.1)
c
      if(ie.le.9) then
         write(temp1,'(i1)') ie
         sgfl_file='sgfl.'//temp1
         sgfr_file='sgfr.'//temp1
         g2d_file='g2d.'//temp1
      else
         write(temp2,'(i2)') ie
         sgfl_file='sgfl.'//temp2
         sgfr_file='sgfr.'//temp2
         g2d_file='g2d.'//temp2
      end if
      open(UNIT=ksgfl,FILE=sgfl_file,
     >     FORM='unformatted',STATUS='unknown')
      open(UNIT=ksgfr,FILE=sgfr_file,
     >     FORM='unformatted',STATUS='unknown')
      open(UNIT=kg2d,FILE=g2d_file,
     >     FORM='unformatted',STATUS='unknown')
c     write(6,*) '<cpacoord>: left sgf file is open  ', sgfl_file
c     write(6,*) '<cpacoord>: right sgf file is open ', sgfr_file
c     write(6,*) '<cpacoord>: g2d file is openi      ', g2d_file
c
c ****************************
c loop over k points in 2D IBZ
c ****************************
c
      do ik=1,nk
        park(1)=xk(ik,1)
        park(2)=xk(ik,2)
c       write(6,'('' K='',2f10.5)') park
c       call flush(6)
c
c ********************************************
c k-resolved layer-indexed structure constants
c ********************************************
c
        if(strc) then
c         -------------------------------------------------
          call gstore(park,psq,eta,bulkgeo,'L',lmax,nintfc,gg2d)   
c         -------------------------------------------------
          io=1
          call getg2d_io(gg2d,g2dio(1,1,-minprc,1,0,ik),
     >                   bulkgeo,lmax,nintfc,io,kg2d)
        else 
          io=-1
          if(itcpa.gt.1) io=0
          call getg2d_io(g2dio(1,1,-minprc,1,0,ik),gg2d,
     >                   bulkgeo,lmax,nintfc,io,kg2d)
        end if
c
c ********************************
c loop over point group operations
c ********************************
c
        do ig=1,ngeff
          ige=ieqg(ig)
          igel=ieqgl(ig)
          iger=ieqgr(ig)
          ig1=invg(ig)
c
c **********************
c surface Green function
c **********************
c
          if(igel.eq.ig) then
            if(sgfy) then
              if(.not.bulkgeo) then
              call gstore(park,psq,eta,.true.,'L',lmax,nintfc,gg2dbulk)
c             ----------------------------------------------------
              call dugo(lmax,kmymax,ttmpl(1,1,1,igel),kmymaxp,
     >                  gg2dbulk(1,1,-minprc,1,0),
     >                  gg2dbulk(1,1,-minprc,1,1),
     >                  yy,mdim,'L',.true.,irel,sgfconv)
c             ----------------------------------------------------
              else
c             ----------------------------------------------------
              call dugo(lmax,kmymax,ttmpl(1,1,1,igel),kmymaxp,
     >                  gg2d(1,1,-minprc,1,0),
     >                  gg2d(1,1,-minprc,1,1),
     >                  yy,mdim,'L',.true.,irel,sgfconv)
c             ----------------------------------------------------
              end if
              if(.not.sgfconv) then
                write(6,'('' Left SGF'')')
                write(6,'('' k='',2f15.6)') park
                stop
              end if
c
              io=1
              call getsfg_io(yy,yio(1,1,ik,igordl(igel)),
     >                       ndiml,mdim,io,ksgfl)
            else 
              io=-1
              if(itcpa.gt.1) io=0
              call getsfg_io(yio(1,1,ik,igordl(igel)),yy,
     >                       ndiml,mdim,io,ksgfl)
            end if
          else
              call repl(yy,yio(1,1,ik,igordl(igel)),ndiml,mdim)
          end if 
c
          if(iger.eq.ig) then
            if(sgfx) then
              if(.not.bulkgeo) then
              call gstore(park,psq,eta,.true.,'R',lmax,nintfc,gg2dbulk)
c             ----------------------------------------------------
              call dugo(lmax,kmymax,ttmpr(1,1,1,iger),kmymaxp,
     >                  gg2dbulk(1,1,-minprc,1,2),
     >                  gg2dbulk(1,1,-minprc,1,1),
     >                  xx,mdim,'R',.true.,irel,sgfconv)
c             ----------------------------------------------------
              else
c             ----------------------------------------------------
              call dugo(lmax,kmymax,ttmpr(1,1,1,iger),kmymaxp,
     >                  gg2d(1,1,-minprc,1,2),
     >                  gg2d(1,1,-minprc,1,1),
     >                  xx,mdim,'R',.true.,irel,sgfconv)
c             ----------------------------------------------------
              endif
              if(.not.sgfconv) then
                write(6,'('' Right SGF'')')
                write(6,'('' k='',2f15.6)') park
                stop
              end if
c
              io=1
              call getsfg_io(xx,xio(1,1,ik,igordr(iger)),
     >                       ndimr,mdim,io,ksgfr)
            else 
              io=-1
              if(itcpa.gt.1) io=0
              call getsfg_io(xio(1,1,ik,igordr(iger)),xx,
     >                       ndimr,mdim,io,ksgfr)
            end if 
          else
              call repl(xx,xio(1,1,ik,igordr(iger)),ndimr,mdim)
          end if 
c
          if(ige.eq.ig) then
c
c ***********************************
c k-resolved tau matrix for interface
c ***********************************
c
            wgeff = wk(ik)/ngeff
c           ----------------------------------------------------
            call tau2d(irel,wgeff,bulkgeo,lmax,kmymax,kmymaxp,nintfc,
     >                 ttmp(1,1,1,ig),gg2d,xx,yy,tau1(1,1,1,ig),mdim)
c           ----------------------------------------------------
c
          end if
c
c **************************************
c BZ sum for layer diagonal tau matrices
c **************************************
c
          do li=1,nintfc
c           -----------------------------------------------
            call repl(tau2,tau1(1,1,li,ige),kmymax,kmymaxp)
c           -----------------------------------------------
            call tripmt(rmatp(1,1,ig1),tau2,rmat(1,1,ig1),
     >                  kmymax,kmymax,kmymaxp)        
c           -----------------------------------------------
            call addmat(tau(1,1,li),tau2,kmymax,kmymaxp)
c           -----------------------------------------------
          end do
c
        end do
c *********************************
c end of loop over point operations
c *********************************
c
      end do
c *************************
c end of loop over k points
c *************************
c
c **********
c CPA solver
c **********
c     -----------------------------------------------------------------
      call cpacor(conc,tminv,tminva,tminvb,tau,taua,taub,kmymax,nintfc,
     >            itcpa,itcpam,cpatol,cpacon,concpa,cpatest,cpaerr)
c     -----------------------------------------------------------------
c
      if(concpa.and.cpacon) goto 100
      if(itest.ge.1) write(6,'('' CPA iterations:'',i3,2x,
     >'' ERROR:'',d15.6)') itcpa-1,cpaerr
c
c                        *************************
c                        * ending CPA iterations *
c                        *************************
c
c *************************************************************
c * loop over layers to transform layer diagonal tau matrices *
c * into physical representation and rotate to the local      *
c * frame of reference                                        *
c *************************************************************
c 
      do li=1,nintfc
         cpalay=(1.d0-conc(li)).gt.tiny
c
c        --------------------------------------------------------
         call phystau(taua(1,1,li),tminva(1,1,li),tminva(1,1,li),
     >                alphaintkkr(0,li),lmax,1)
c        --------------------------------------------------------
c
c rotate tau-matrix to local frame of reference
c
         call repl(gtaua(1,1,li),taua(1,1,li),kmymax,kmymaxp)
c        ------------------------------------------------------
         call tripmt(dmatpa(1,1,li),taua(1,1,li),dmata(1,1,li),
     >               kmymax,kmymax,kmymaxp)
c        ------------------------------------------------------
c
         if(itest.gt.2) then
            write(6,'(/'' tau -A '',i2)') li
            call outmat1(taua(1,1,li),kmymax,kmymax,kmymaxp,tol,6)
         end if
c
         if(cpalay) then
c+------------+
c+ BIG CPA IF +
c+------------+
c
c        --------------------------------------------------------
         call phystau(taub(1,1,li),tminvb(1,1,li),tminvb(1,1,li),
     >                alphaintkkr(0,li),lmax,1)
c        --------------------------------------------------------
c
c rotate tau-matrix to local frame of reference
c
         call repl(gtaub(1,1,li),taub(1,1,li),kmymax,kmymaxp)
c        ---------------------------------------------------
         call tripmt(dmatpb(1,1,li),taub(1,1,li),dmatb(1,1,li),
     >               kmymax,kmymax,kmymaxp)
c        ---------------------------------------------------
c
         if(itest.gt.2) then
            write(6,'(/'' tau -B '',i2)') li
            call outmat1(taub(1,1,li),kmymax,kmymax,kmymaxp,tol,6)
         end if
c
         end if
c+----------------+
c+ END BIG CPA IF +
c+----------------+
      end do
c *** end loop over layers ***
c
      return
      end
