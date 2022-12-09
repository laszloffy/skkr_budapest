c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
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
      parameter(mek=mekr) 
c
      logical wrel,bulk,bulkgeo,sgf,sgfx,sgfy,sgfconv
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
      complex*16 xsave(mdim,mdim,mek+1,melem1),
     >           ysave(mdim,mdim,mek+1,melem1),
     > g2d(lmmaxp,lmmaxp,-minprc:minprc,minprc,0:mprc1+1,mek+1)
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
      common/test/itest
      common/lay2d/cvec(mtotal,3),nextra,nbulkl,nbulkr,
     &             nprc,ninprc(0:mprc+1)
c
      common/relfac/fac
c
      data tol/1.0d-8/ 
      data tiny/1.0d-6/
c
      save xsave,ysave,g2d
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
c ****************************
c loop over k points in 2D IBZ
c ****************************
c
      do ik=1,nk
        park(1)=xk(ik,1)
        park(2)=xk(ik,2)
c
c *************************************************
c initialize saving of structure constants and SSPO
c *************************************************
c
        iek=iek+1
        if(iek.le.mek) then
          iekuse=iek
        else
          iekuse=mek+1
        end if
        if(iek.eq.mek+1) 
     >  write(6,'(/''WARNING: iek='',i3,'' greater than mek='',i3/
     >             ''         causes dramatic slow-down !'')') iek,mek
c
        if(bulk) then
c    calculate g2d only for the first CPA iteration
          nocalcstr=itcpa-1
        else
c    calculate g2d only for the first scf and CPA iteration
          nocalcstr=itscf+itcpa-2
        end if
c    do not store g2d for iek.gt.mek
        if(iek.gt.mek) nocalcstr=0
c
        if(bulk) then
c    always calculate SSPO
          sgf=.true.
        else
c    calculate SSPO only for the first scf and CPA iteration
c    or for iek.gt.mek
          sgf=(itcpa.eq.1.and.itscf.eq.1).or.(iek.gt.mek)
        end if
        sgfy=sgf
        sgfx=sgf
c    if vacuum, recalculate R SSPO after every tenth step
        if(rightm.eq.'V')
     >   sgfx=((mod(itscf,10).eq.1).and.(itcpa.eq.1)).or.(iek.gt.mek)
c
c ********************************************
c k-resolved layer-indexed structure constants
c ********************************************
c
        if(nocalcstr.ne.0) then
          call getg2d(g2d(1,1,-minprc,1,0,iekuse),gg2d,
     >                bulkgeo,lmax,nintfc)
        else
c         -------------------------------------------------
          call gstore(park,psq,eta,bulkgeo,'L',lmax,nintfc,gg2d)   
c         -------------------------------------------------
          call getg2d(gg2d,g2d(1,1,-minprc,1,0,iekuse),
     >                bulkgeo,lmax,nintfc)
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
          if((igel.eq.ig).and.sgfy) then
            if(.not.bulkgeo) then
            call gstore(park,psq,eta,.true.,'L',lmax,nintfc,gg2dbulk)   
c           ----------------------------------------------------
            call dugo(lmax,kmymax,ttmpl(1,1,1,igel),kmymaxp,
     >                gg2dbulk(1,1,-minprc,1,0),
     >                gg2dbulk(1,1,-minprc,1,1),
     >                yy,mdim,'L',.true.,irel,sgfconv)
c           ----------------------------------------------------
            else
c           ----------------------------------------------------
            call dugo(lmax,kmymax,ttmpl(1,1,1,igel),kmymaxp,
     >                gg2d(1,1,-minprc,1,0),
     >                gg2d(1,1,-minprc,1,1),
     >                yy,mdim,'L',.true.,irel,sgfconv)
c           ----------------------------------------------------
            end if
            if(.not.sgfconv) then
             write(6,'('' Left SGF'')')
             write(6,'('' k='',2f15.6)') park
             stop
            end if
c
            call repl(ysave(1,1,iekuse,igordl(igel)),yy,ndiml,mdim)
          else
            call repl(yy,ysave(1,1,iekuse,igordl(igel)),ndiml,mdim)
          end if 
c
          if((iger.eq.ig).and.sgfx) then
            if(.not.bulkgeo) then
            call gstore(park,psq,eta,.true.,'R',lmax,nintfc,gg2dbulk)   
c           ----------------------------------------------------
            call dugo(lmax,kmymax,ttmpr(1,1,1,iger),kmymaxp,
     >                gg2dbulk(1,1,-minprc,1,2),
     >                gg2dbulk(1,1,-minprc,1,1),
     >                xx,mdim,'R',.true.,irel,sgfconv)
c           ----------------------------------------------------
            else 
c           ----------------------------------------------------
            call dugo(lmax,kmymax,ttmpr(1,1,1,iger),kmymaxp,
     >                gg2d(1,1,-minprc,1,2),
     >                gg2d(1,1,-minprc,1,1),
     >                xx,mdim,'R',.true.,irel,sgfconv)
c           ----------------------------------------------------
            end if
            if(.not.sgfconv) then
             write(6,'('' Right SGF'')')
             write(6,'('' k='',2f15.6)') park
             stop
            end if
c
            call repl(xsave(1,1,iekuse,igordr(iger)),xx,ndimr,mdim)
          else
            call repl(xx,xsave(1,1,iekuse,igordr(iger)),ndimr,mdim)
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
