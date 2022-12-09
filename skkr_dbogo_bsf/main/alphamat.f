c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine alphamat(celoc,iscreen,vscreen,lmax,E_Fermi,
     >                    nintfc,ninprcl,ninprcr,c_light)
c=========================
c
c set up matrix alpha as the non-relativistic t-matrix of a
c square well
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      complex*16 alphalkkr,alpharkkr,alphaintkkr,
     >           ce,cs,ck,csk,tand,react,sqrm1,celoc
      complex*16 alphalkkrh,alpharkkrh,alphaintkkrh
      complex*16 fmb(0:l1maxp),fmn(0:l1maxp),fmh(0:l1maxp)
      complex*16 fmodb(0:l1maxp),fmodn(0:l1maxp),fmodh(0:l1maxp)
c
      common/muftin/smtl(minprc),smt(mintfc),smtr(minprc)
      common/scrpar/alphalkkr(0:lmaxp,minprc),alpharkkr(0:lmaxp,minprc),
     >              alphaintkkr(0:lmaxp,mintfc)
      common/scrparh/alphalkkrh(0:lmaxp,minprc),
     >              alpharkkrh(0:lmaxp,minprc),
     >              alphaintkkrh(0:lmaxp,mintfc)
      common/test/itest
c
      data sqrm1/(0.0d0,1.0d0)/
c
c      c=274.072d0 -- from input
      c=c_light
c
c
      ce=celoc+E_Fermi+
     >  (celoc+E_Fermi)*(celoc+E_Fermi)/(c*c)
      cs=ce-vscreen
      csk=cdsqrt(cs)
c     if(dimag(csk).lt.0.d0) csk=-csk
      ck=cdsqrt(ce)
c     if(dimag(ck).lt.0.d0) ck=-ck

c -- screening parameter for left and right

      do n=1,ninprcl
         call csbf(lmax+1,l1maxp,ck,smtl(n),fmb,fmn,fmh)
         call csbf(lmax+1,l1maxp,csk,smtl(n),fmodb,fmodn,fmodh)
         do l=0,lmax
           if(iscreen.le.1) then
             tand=(ck*fmb(l+1)*fmodb(l)-csk*fmb(l)*fmodb(l+1))/
     >            (ck*fmn(l+1)*fmodb(l)-csk*fmn(l)*fmodb(l+1))
           else
             tand=fmb(l)/fmn(l)
           end if
           react=-tand/ck
           alphalkkr(l,n)=react/(1.0d0+sqrm1*ck*react)
         end do
      end do
c
      do n=1,ninprcr
         call csbf(lmax+1,l1maxp,ck,smtr(n),fmb,fmn,fmh)
         call csbf(lmax+1,l1maxp,csk,smtr(n),fmodb,fmodn,fmodh)
         do l=0,lmax
           if(iscreen.le.1) then
             tand=(ck*fmb(l+1)*fmodb(l)-csk*fmb(l)*fmodb(l+1))/
     >            (ck*fmn(l+1)*fmodb(l)-csk*fmn(l)*fmodb(l+1))
           else
             tand=fmb(l)/fmn(l)
           end if
           react=-tand/ck
           alpharkkr(l,n)=react/(1.0d0+sqrm1*ck*react)
         end do
      end do
c
c -- now calculate the screening par. for the interface --
c
      do n=1,nintfc
         call csbf(lmax+1,l1maxp,ck,smt(n),fmb,fmn,fmh)
         call csbf(lmax+1,l1maxp,csk,smt(n),fmodb,fmodn,fmodh)
         do l=0,lmax
           if(iscreen.le.1) then
             tand=(ck*fmb(l+1)*fmodb(l)-csk*fmb(l)*fmodb(l+1))/
     >            (ck*fmn(l+1)*fmodb(l)-csk*fmn(l)*fmodb(l+1))
           else
             tand=fmb(l)/fmn(l)
           end if
           react=-tand/ck
           alphaintkkr(l,n)=react/(1.0d0+sqrm1*ck*react)
         end do
      end do
c
c       ------------ holes start -----------
c
      ce=-(celoc)+E_Fermi+
     >  (-(celoc)+E_Fermi)*(-(celoc)+E_Fermi)/(c*c)
      cs=ce-vscreen
      csk=-cdsqrt(cs)
c     if(dimag(csk).gt.0.d0) csk=-csk
      ck=-cdsqrt(ce)
c     if(dimag(ck).gt.0.d0) ck=-ck
c
c -- screening parameter for left and right
c
      do n=1,ninprcl
         call csbf(lmax+1,l1maxp,ck,smtl(n),fmb,fmn,fmh)
         call csbf(lmax+1,l1maxp,csk,smtl(n),fmodb,fmodn,fmodh)
         do l=0,lmax
           if(iscreen.le.1) then
             tand=(ck*fmb(l+1)*fmodb(l)-csk*fmb(l)*fmodb(l+1))/
     >            (ck*fmn(l+1)*fmodb(l)-csk*fmn(l)*fmodb(l+1))
           else
             tand=fmb(l)/fmn(l)
           end if
           react=-tand/ck
           alphalkkrh(l,n)=-react/(1.0d0+sqrm1*ck*react)
         end do
      end do
c
      do n=1,ninprcr
         call csbf(lmax+1,l1maxp,ck,smtr(n),fmb,fmn,fmh)
         call csbf(lmax+1,l1maxp,csk,smtr(n),fmodb,fmodn,fmodh)
         do l=0,lmax
           if(iscreen.le.1) then
             tand=(ck*fmb(l+1)*fmodb(l)-csk*fmb(l)*fmodb(l+1))/
     >            (ck*fmn(l+1)*fmodb(l)-csk*fmn(l)*fmodb(l+1))
           else
             tand=fmb(l)/fmn(l)
           end if
           react=-tand/ck
           alpharkkrh(l,n)=-react/(1.0d0+sqrm1*ck*react)
         end do
      end do
c
c -- now calculate the screening par. for the interface --
c
      do n=1,nintfc
         call csbf(lmax+1,l1maxp,ck,smt(n),fmb,fmn,fmh)
         call csbf(lmax+1,l1maxp,csk,smt(n),fmodb,fmodn,fmodh)
         do l=0,lmax
           if(iscreen.le.1) then
             tand=(ck*fmb(l+1)*fmodb(l)-csk*fmb(l)*fmodb(l+1))/
     >            (ck*fmn(l+1)*fmodb(l)-csk*fmn(l)*fmodb(l+1))
           else
             tand=fmb(l)/fmn(l)
           end if
           react=-tand/ck
           alphaintkkrh(l,n)=-react/(1.0d0+sqrm1*ck*react)
         end do
      end do
c
      if(itest.lt.3) return
      write(6,'(/2x,''routine ALPHAMAT>'')')
      do n=1,ninprcl
         write(6,'(''Electron part smtl(n):  n='',i2,4x,f14.12)') 
     >         n,smtl(n)
         write(6,'(10f12.8)') (alphalkkr(l,n),l=0,lmax)
         write(6,'(''Hole part smtl(n):  n='',i2,4x,f14.12)') n,smtl(n)
         write(6,'(10f12.8)') (alphalkkrh(l,n),l=0,lmax)
      enddo
      do n=1,ninprcr
         write(6,'(''Electron part smtr(n):  n='',i2,4x,f14.12)') 
     >         n,smtr(n)
         write(6,'(10f12.8)') (alpharkkr(l,n),l=0,lmax)
         write(6,'(''Hole part smtr(n):  n='',i2,4x,f14.12)') n,smtr(n)
         write(6,'(10f12.8)') (alpharkkrh(l,n),l=0,lmax)
      enddo
      do n=1,nintfc
         write(6,'(''Electron part smt(n):  n='',i2,4x,f14.12)') 
     >         n,smt(n)
         write(6,'(10f12.8)') (alphaintkkr(l,n),l=0,lmax)
         write(6,'(''Hole part smt(n):  n='',i2,4x,f14.12)') n,smt(n)
         write(6,'(10f12.8)') (alphaintkkrh(l,n),l=0,lmax)
      end do
c
      return
      end
