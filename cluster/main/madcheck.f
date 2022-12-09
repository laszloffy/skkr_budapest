c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
c=======================================================================
      subroutine madcheck(ctol,ntotal,igraph)
c
      implicit real*8(a-h,o-z)
      include '../param.h'
c
      dimension ctry(2)
c
      common/brav2d/a1(2),a2(2)
      common/lay2d/cvec(mtotal,3),nextra,nbulkl,nbulkr,
     &             nprc,ninprc(0:mprc+1)
c
      if(igraph.gt.-2.and.igraph.lt.2) then
c        write(6,'(/,'' WARNING!!'',/,
c    &   ''Check that all layers share the same 2-D cell '',
c    &   ''around the origin'')')
         return
      endif
c
      do il=1,ntotal
c
        cmod0=dsqrt(cvec(il,1)**2+cvec(il,2)**2)
        jok=-1
c
c       First, choose the smallest circle around Z.
        dowhile(jok.ne.0)
          do 20 n=-1,1
          do 20 m=-1,1
            if(n.eq.0.and.m.eq.0) goto 20
            do i=1,2
              ctry(i)=cvec(il,i)+dfloat(n)*a1(i)+dfloat(m)*a2(i)
            enddo
            cmod=dsqrt(ctry(1)**2+ctry(2)**2)
            if(cmod.lt.dabs(cmod0-ctol)) then
              do i=1,2
                cvec(il,i)=ctry(i)
              enddo
              cmod0=cmod
              jok=1
            endif
 20       continue
          if(jok.lt.0) then
            jok=0
          else
            jok=-1
          endif
        enddo
c
c       Now, select the closest Wigner-Seitz cells among all 
c       physical planes.
          cx = cvec(il,1)+ctol*1.0d-2
          cy = cvec(il,2)+ctol*1.0d-2
          cmod=dsqrt(cvec(il,1)**2+cvec(il,2)**2)
          if((cx.lt.0.d0).or.(cy.lt.0.d0)) then
            jok=-1
            if(cy.lt.0.d0) then
               jok=jok-1
               if(cx.lt.0.d0) jok=jok-1
            endif
            ctry(1)=cvec(il,1)
            ctry(2)=cvec(il,2)
            do n=-1,1
            do m=-1,1
              do i=1,2
                ctry(i)=cvec(il,i)+dfloat(n)*a1(i)+dfloat(m)*a2(i)
              enddo
              cmod0=dsqrt(ctry(1)**2+ctry(2)**2)
              if(dabs(cmod-cmod0).le.ctol) then
                cx=ctry(1)+ctol*1.0d-2
                cy=ctry(2)+ctol*1.0d-2
                if((jok.eq.-3).and.((cx.ge.0.d0).or.(cy.ge.0.d0))) then
                  cvec(il,1)=ctry(1)
                  cvec(il,2)=ctry(2)
                  if(cy.lt.0.d0) then
                    jok=-2
                  elseif(cx.lt.0.d0) then
                    jok=-1
                  else
                    jok=0
                  endif
                elseif((jok.eq.-2).and.(cy.ge.0.d0)) then
                  cvec(il,1)=ctry(1)
                  cvec(il,2)=ctry(2)
                  if(cx.lt.0.d0) then
                    jok=-1
                  else
                    jok=0
                  endif
                elseif((jok.eq.-1).and.(cx.ge.0.d0).and.(cy.ge.0.d0)) 
     &          then
                  cvec(il,1)=ctry(1)
                  cvec(il,2)=ctry(2)
                  jok=0
                endif
                if(jok.eq.0) goto 100
              endif
            enddo
            enddo
 100        continue
          endif
c
      enddo
c
      return
      end
