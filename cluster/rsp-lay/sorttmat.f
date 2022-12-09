c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine sorttmat(
c     ===================
     > lmax,nintfc,ninprcl,ninprcr,rmat,rmatp,ngeff,invg,
     > tminvl,tminvr,tminv,ttmpl,ttmpr,ttmp,
     > ieqgl,ieqgr,ieqg,igordl,igordr,isave0)
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      dimension invg(melem),ieqg(melem),ieqgl(melem),ieqgr(melem)
      dimension igordl(melem),igordr(melem)
c
      complex*16 tminv(kmymaxp,kmymaxp,mintfc)
      complex*16 tminvl(kmymaxp,kmymaxp,minprc)
      complex*16 tminvr(kmymaxp,kmymaxp,minprc)
      complex*16 ttmp(kmymaxp,kmymaxp,mintfc,melem)
      complex*16 ttmpl(kmymaxp,kmymaxp,minprc,melem)
      complex*16 ttmpr(kmymaxp,kmymaxp,minprc,melem)
c
      complex*16 rmat(kmymaxp,kmymaxp,melem)
      complex*16 rmatp(kmymaxp,kmymaxp,melem)
c
      common/test/itest
c
c     data tol/1.0d-8/  Default value
      data tol/1.0d-8 / 
      data tiny/1.0d-6/
c
c ********************
c initialize constants
c ********************
c
      nl=lmax+1
      nl2=nl*nl
      kmax=2*lmax+1
      kmymax=2*nl2
c
c *******************************************************************
c transform inverse t-matrices with respect to point group operations
c *******************************************************************
c
      do ig=1,ngeff
        ig1=invg(ig)
c
        do li=1,ninprcl
c         ----------------------------------------------------------
          call repl(ttmpl(1,1,li,ig),tminvl(1,1,li),kmymax,kmymaxp)
c         ----------------------------------------------------------
          call tripmt(rmat(1,1,ig1),ttmpl(1,1,li,ig),rmatp(1,1,ig1),
     >                kmymax,kmymax,kmymaxp)        
c         ----------------------------------------------------------
          if(itest.gt.2) then
            write(6,'(/'' tminv left'',i2)') li
            call outmat1(ttmpl(1,1,li,ig),kmymax,kmymax,kmymaxp,tol,6)
          end if
        end do
c
        do li=1,ninprcr
c         ----------------------------------------------------------
          call repl(ttmpr(1,1,li,ig),tminvr(1,1,li),kmymax,kmymaxp)
c         ----------------------------------------------------------
          call tripmt(rmat(1,1,ig1),ttmpr(1,1,li,ig),rmatp(1,1,ig1),
     >                kmymax,kmymax,kmymaxp)        
c         ----------------------------------------------------------
          if(itest.gt.2) then
            write(6,'(/'' tminv right'',i2)') li
            call outmat1(ttmpr(1,1,li,ig),kmymax,kmymax,kmymaxp,tol,6)
          end if
        end do
c
        do li=1,nintfc
c         ----------------------------------------------------------
          call repl(ttmp(1,1,li,ig),tminv(1,1,li),kmymax,kmymaxp)
c         ----------------------------------------------------------
          call tripmt(rmat(1,1,ig1),ttmp(1,1,li,ig),rmatp(1,1,ig1),
     >                kmymax,kmymax,kmymaxp)        
c         ----------------------------------------------------------
          if(itest.gt.2) then
            write(6,'(/'' tminv '',i2)') li
            call outmat1(ttmp(1,1,li,ig),kmymax,kmymax,kmymaxp,tol,6)
          end if
        end do
c
      end do
c
c ******************
c check degeneracies
c ******************
c
c   First check the WHOLE system:
c   if all "ic" within a "iig" loop are different from zero, that is, 
c   all symmetry rotated t-matrices are equivalent in each layer, 
c   then ieqg(ig)=1 for all ig.  
c                           ("ig" runs over ngeff symmetry operations)
c   else, ieqg(ig)=x, x being the smallest iig for which all "ic.ne.0".
c
      ieqg(1)=1
      do 7 ig=2,ngeff
        do 5 iig=1,ig-1
c
          do li=1,ninprcl
c           ------------------------------------------------
            call compmat(ttmpl(1,1,li,iig),ttmpl(1,1,li,ig),
     >                   kmymax,kmymaxp,tol,ic)
c           ------------------------------------------------
            if(ic.eq.0) goto 5
          enddo
c
          do li=1,ninprcr
c           ------------------------------------------------
            call compmat(ttmpr(1,1,li,iig),ttmpr(1,1,li,ig),
     >                   kmymax,kmymaxp,tol,ic)
c           ------------------------------------------------
            if(ic.eq.0) goto 5
          enddo
c
          do li=1,nintfc
c           ------------------------------------------------
            call compmat(ttmp(1,1,li,iig),ttmp(1,1,li,ig),
     >                   kmymax,kmymaxp,tol,ic)
c           ------------------------------------------------
            if(ic.eq.0) goto 5
          enddo
c
          goto 6
  5     continue
        ieqg(ig)=ig
        goto 7
  6     continue
        ieqg(ig)=iig
  7   continue
c
      if(itest.ge.2) write(6,'(/12i3)') (ieqg(ig),ig=1,ngeff)
c
c Non-magnetic bulk (ecore_save0)
c
      if(isave0.eq.1) return
c
c   LEFT side:
c   if all "ic" within a "iig" loop are different from zero, then
c   ieqgl(ig)=1 for all ig and igordl(1)=1.  
c   else, ieqgl(ig)=x, x being the smallest iig for which all "ic.ne.0",
c   and igordl(x)=1+ieqgl(x-1).
c
      ieqgl(1)=1
      igordl(1)=1
      do 17 ig=2,ngeff
        do 15 iig=1,ig-1
c
          do li=1,ninprcl
c           ------------------------------------------------
            call compmat(ttmpl(1,1,li,iig),ttmpl(1,1,li,ig),
     >                   kmymax,kmymaxp,tol,ic)
c           ------------------------------------------------
c
            if(ic.eq.0) goto 15
          enddo
c         
          goto 16
 15     continue
        ieqgl(ig)=ig
        igordl(ig)=igordl(ieqgl(ig-1))+1
        goto 17
 16     continue
        ieqgl(ig)=iig
 17   continue
c
c   RIGHT side:
c                    (same as LEFT side, now for ieqgr,igordr)
      ieqgr(1)=1
      igordr(1)=1
      do 27 ig=2,ngeff
        do 25 iig=1,ig-1
c
          do li=1,ninprcr
c           ------------------------------------------------
            call compmat(ttmpr(1,1,li,iig),ttmpr(1,1,li,ig),
     >                 kmymax,kmymaxp,tol,ic)
c           ------------------------------------------------
            if(ic.eq.0) goto 25
          enddo
c
          goto 26
 25     continue
        ieqgr(ig)=ig
        igordr(ig)=igordr(ieqgr(ig-1))+1
        goto 27
 26     continue
        ieqgr(ig)=iig
 27   continue
c
      if(itest.ge.2) then
       write(6,'(/12i3)') (ieqgl(ig),ig=1,ngeff)
       write(6,'(12i3)') (igordl(ieqgl(ig)),ig=1,ngeff)
       write(6,'(/12i3)') (ieqgr(ig),ig=1,ngeff)
       write(6,'(12i3)') (igordr(ieqgr(ig)),ig=1,ngeff)
      end if
c
      return
      end
