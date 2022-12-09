c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine pothandle(
     > bulk,linbw,rightm,lmax,nintfc,ninprcl,ninprcr,
     > vrsbulk,v00,ivacpot,opotl,opotr,opot,
     > leftpot,leftmom,concl,qmomla,qmomlb,dxl,nsl,rsl,
     > idpotla,vrla,brla,boprla,zla,idpotlb,vrlb,brlb,boprlb,zlb,
     > rightpot,rightmom,concr,qmomra,qmomrb,dxr,nsr,rsr,
     > idpotra,vrra,brra,boprra,zra,idpotrb,vrrb,boprrb,brrb,zrb,
     > laypot,conc,dx,ns,rs,
     > idpota,vra,bra,bopra,za,idpotb,vrb,brb,boprb,zb,igraph)
c
c -read Left,Right and Layer potentials
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      logical bulk,cpalay,opotl,opotr,opot,linbw,rspkkrin
      character*10 idpota(mintfc),idpotla(minprc),idpotra(minprc)
      character*10 idpotb(mintfc),idpotlb(minprc),idpotrb(minprc)
      character*30 leftpot,rightpot,laypot
      character*30 leftmom,rightmom
      character*1 rightm
c
      dimension vra(nrad,mintfc),vrb(nrad,mintfc)
      dimension bra(nrad,mintfc),brb(nrad,mintfc)
      dimension bopra(nrad,2,mintfc),boprb(nrad,2,mintfc)
      dimension dx(mintfc),ns(mintfc),rs(mintfc)
      dimension vrla(nrad,minprc),vrlb(nrad,minprc)
      dimension brla(nrad,minprc),brlb(nrad,minprc)
      dimension boprla(nrad,2,minprc),boprlb(nrad,2,minprc)
      dimension rsl(minprc),dxl(minprc),nsl(minprc)
      dimension vrra(nrad,minprc),vrrb(nrad,minprc)
      dimension brra(nrad,minprc),brrb(nrad,minprc)
      dimension boprra(nrad,2,minprc),boprrb(nrad,2,minprc)
      dimension rsr(minprc),dxr(minprc),nsr(minprc)
      dimension za(mintfc),zb(mintfc),zla(minprc),zlb(minprc),
     &          zra(minprc),zrb(minprc) 
      dimension conc(mintfc),concl(minprc),concr(minprc)
      dimension isigba(mintfc),isigbb(mintfc)
      dimension qq(2*lmsup)
c
      complex*16 qmomla(lmsup,minprc),qmomlb(lmsup,minprc)
      complex*16 qmomra(lmsup,minprc),qmomrb(lmsup,minprc)
c
      common/pot/potshift,b0,layshift1,layshift2,layb0,layb0f,
     &           ib0,isigba,isigbb,nomag
c
      data tiny/1.0d-6/
c
c -- rsl/~r is the Wigner-Seitz radius for the left and right bulk and
c    rs(i),i=1,nintfc is the WS-radius for the interface !!
c
      if(bulk) goto 100
c
      write(6,*) 'leftpot=',leftpot
      write(6,*) 'leftmom=',leftmom
      open(20,file=leftpot,status='old')
      open(21,file=leftmom,status='old')
      read(20,*) ns0,dx0,rspkkrin
      read(21,*)
      read(21,*) lmaxin
      lmaxs=2*lmaxin
      lmmaxs=(lmaxs+1)*(lmaxs+1)
c
      do li=1,ninprcl
c
        nsl(li)=ns0
        dxl(li)=dx0
        call readpot(
     >  idpotla(li),conc(li),vrla(1,li),brla(1,li),boprla(1,1,li),
     >  nomag,opotl,dxl(li),nsl(li),rsl(li),zla(li),potzero,rspkkrin)
c
        if(li.eq.1) vrsbulk=potzero
c
        read(21,*)
        read(21,*)
        read(21,*) (qq(i),i=1,2*lmmaxs)
        do i=1,lmmaxs
          qmomla(i,li)=dcmplx(qq(2*i-1),qq(2*i))
        end do
c
        cpalay=1.d0-concl(li).gt.tiny
        if(cpalay) then
c
          call readpot(
     >    idpotlb(li),conc0,vrlb(1,li),brlb(1,li),boprlb(1,1,li),
     >    nomag,opotl,dxl(li),nsl(li),rsl(li),zlb(li),potzero,rspkkrin)
c
          read(21,*)
          read(21,*)
          read(21,*) (qq(i),i=1,2*lmmaxs)
          do i=1,lmmaxs
            qmomlb(i,li)=dcmplx(qq(2*i-1),qq(2*i))
          end do
c
        else
c
          idpotlb(li)=idpotla(li)
          zlb(li)=zla(li)
          do i=1,lmmaxs
             qmomlb(i,li)=qmomla(i,li)
          end do
          do i=1,nsl(li)
            vrlb(i,li)=vrla(i,li)
            brlb(i,li)=brla(i,li)
            boprlb(i,1,li)=boprla(i,1,li)
            boprlb(i,2,li)=boprla(i,2,li)
          end do
        end if
        write(6,'(i4,2x,a3,''('',f5.3,'') - '',a3,''('',f5.3,'')'')')
     >  li,idpotla(li)(1:3),concl(li),idpotlb(li)(1:3),1.d0-concl(li)
        call flush(6)
      enddo
      close(20)
      close(21)
c
      if(linbw) then
        write(6,*) 'rightpot'
        do li=1,ninprcr
          idpotra(li)=idpotla(li)
          idpotrb(li)=idpotlb(li)
          zra(li)=zla(li)
          zrb(li)=zlb(li)
          nsr(li)=nsl(li)
          dxr(li)=dxl(li)
          rsr(li)=rsl(li)
          concr(li)=concl(li)
          do i=1,nsl(li)
            vrra(i,li)=vrla(i,li)
            vrrb(i,li)=vrlb(i,li)
            brra(i,li)=brla(i,li)
            brrb(i,li)=brlb(i,li)
            boprra(i,1,li)=boprla(i,1,li)
            boprra(i,2,li)=boprla(i,2,li)
            boprrb(i,1,li)=boprlb(i,1,li)
            boprrb(i,2,li)=boprlb(i,2,li)
          end do
          write(6,'(i4,2x,a3,''('',f5.3,'') - '',a3,''('',f5.3,'')'')')
     >    li,idpotra(li)(1:3),concr(li),idpotrb(li)(1:3),1.d0-concr(li)
          call flush(6)
        end do
        write(6,*) 'laypot'
        do li=1,nintfc
          lli=li-((li-1)/ninprcl)*ninprcl
          idpota(li)=idpotla(lli)
          idpotb(li)=idpotlb(lli)
          za(li)=zla(lli)
          zb(li)=zlb(lli)
          ns(li)=nsl(lli)
          dx(li)=dxl(lli)
          rs(li)=rsl(lli)
          conc(li)=concl(lli)
          do i=1,ns(li)
            vra(i,li)=vrla(i,lli)
            vrb(i,li)=vrlb(i,lli)
            bra(i,li)=brla(i,lli)
            brb(i,li)=brlb(i,lli)
            bopra(i,1,li)=boprla(i,1,lli)
            bopra(i,2,li)=boprla(i,2,lli)
            boprb(i,1,li)=boprlb(i,1,lli)
            boprb(i,2,li)=boprlb(i,2,lli)
          end do
          write(6,'(i4,2x,a3,''('',f5.3,'') - '',a3,''('',f5.3,'')'')')
     >    li,idpota(li)(1:3),conc(li),idpotb(li)(1:3),1.d0-conc(li)
          call flush(6)
        end do
        return
      end if

      if(rightm.eq.'V') then
c
        write(6,*) 'rightpot=  Vacuum' 
        do li=1,ninprcr
          concr(li)=1.d0
          rsr(li)=rsl(li)
          nsr(li)=nsl(li)
          dxr(li)=dxl(li)
          idpotra(li)='Vacuum    '
          zra(li)=0.d0
          idpotrb(li)='Vacuum    '
          zrb(li)=0.d0
          call czero(qmomra(1,li),lmsup)
          call czero(qmomrb(1,li),lmsup)
          write(6,'(2x,a3,''('',f5.3,'') - '',a3,''('',f5.3,'')'')')
     >    idpotra(li)(1:3),concr(li),idpotrb(li)(1:3),1.d0-concr(li)
        enddo
c
      else
c
        write(6,*) 'rightpot=',rightpot
        write(6,*) 'rightmom=',rightmom
        open(20,file=rightpot,status='old')
        open(21,file=rightmom,status='old')
        read(20,*) ns0,dx0,rspkkrin
        read(21,*)
        read(21,*) lmaxin
        lmaxs=2*lmaxin
        lmmaxs=(lmaxs+1)*(lmaxs+1)
c
        do li=1,ninprcr
          nsr(li)=ns0
          dxr(li)=dx0
          call readpot(
     >    idpotra(li),conc0,vrra(1,li),brra(1,li),boprra(1,1,li),
     >    nomag,opotr,dxr(li),nsr(li),rsr(li),zra(li),potzero,rspkkrin)
c
          read(21,*)
          read(21,*)
          read(21,*) (qq(i),i=1,2*lmmaxs)
          do i=1,lmmaxs
            qmomra(i,li)=dcmplx(qq(2*i-1),qq(2*i))
          end do
c
          cpalay=1.d0-concr(li).gt.tiny
          if(cpalay) then
c
          call readpot(
     >    idpotrb(li),conc0,vrrb(1,li),brrb(1,li),boprrb(1,1,li),
     >    nomag,opotr,dxr(li),nsr(li),rsr(li),zrb(li),potzero,rspkkrin)
c
          read(21,*)
          read(21,*)
          read(21,*) (qq(i),i=1,2*lmmaxs)
          do i=1,lmmaxs
            qmomrb(i,li)=dcmplx(qq(2*i-1),qq(2*i))
          end do
c
          else
c
            idpotrb(li)=idpotra(li)
            zrb(li)=zra(li)
            do i=1,lmmaxs
              qmomrb(i,li)=qmomra(i,li)
            end do
            do i=1,nsr(li)
              vrrb(i,li)=vrra(i,li)
              brrb(i,li)=brra(i,li)
              boprrb(i,1,li)=boprra(i,1,li)
              boprrb(i,2,li)=boprra(i,2,li)
            end do
c
          end if
          write(6,'(2x,a3,''('',f5.3,'') - '',a3,''('',f5.3,'')'')')
     >    idpotra(li)(1:3),concr(li),idpotrb(li)(1:3),1.d0-concr(li)
        enddo
        close(20)
        close(21)
      end if
c
  100 continue
      open(20,file=laypot,status='old')
      read(20,*) ns0,dx0,rspkkrin
      write(6,*) 'laypot=',laypot
      do li=1,nintfc
        ns(li)=ns0
        dx(li)=dx0
        call readpot(
     >  idpota(li),conc(li),vra(1,li),bra(1,li),bopra(1,1,li),
     >  nomag,opot,dx(li),ns(li),rs(li),za(li),potzero,rspkkrin)
c
        if(bulk.and.li.eq.1) vrsbulk=potzero
c
        if((li.ge.layshift1).and.(li.le.layshift2)) then
          x=dlog(rs(li))-(ns(li)-1)*dx(li)
          do i=1,ns(li)
            r=dexp(x)
            vra(i,li)=vra(i,li)+potshift*r
            x=x+dx(li)
          end do
        end if
        if(idpota(li)(1:3).eq.'Vac'.and.ivacpot.eq.1) then
          x=dlog(rs(li))-(ns(li)-1)*dx(li)
          do i=1,ns(li)
            r=dexp(x)
            vra(i,li)=v00*r
            bra(i,li)=0.d0 
            bopra(i,1,li)=0.d0 
            bopra(i,2,li)=0.d0 
            x=x+dx(li)
          end do
        end if
c
        cpalay=1.d0-conc(li).gt.tiny
        if(cpalay) then
c
          call readpot(
     >    idpotb(li),conc0,vrb(1,li),brb(1,li),boprb(1,1,li),
     >    nomag,opot,dx(li),ns(li),rs(li),zb(li),potzero,rspkkrin)
c
          if((li.ge.layshift1).and.(li.le.layshift2)) then
            x=dlog(rs(li))-(ns(li)-1)*dx(li)
            do i=1,ns(li)
              r=dexp(x)
              vrb(i,li)=vrb(i,li)+potshift*r
              x=x+dx(li)
            end do
          end if
          if(idpotb(li)(1:3).eq.'Vac'.and.ivacpot.eq.1) then
            x=dlog(rs(li))-(ns(li)-1)*dx(li)
            do i=1,ns(li)
              r=dexp(x)
              vrb(i,li)=v00*r
              brb(i,li)=0.d0 
              boprb(i,1,li)=0.d0 
              boprb(i,2,li)=0.d0 
              x=x+dx(li)
            end do
          end if
c
        else
c
          idpotb(li)=idpota(li)
          zb(li)=za(li)
          do i=1,ns(li)
            vrb(i,li)=vra(i,li)
            brb(i,li)=bra(i,li)
            boprb(i,1,li)=bopra(i,1,li)
            boprb(i,2,li)=bopra(i,2,li)
          end do
c
        end if
c
        if((layb0.eq.0).or.((li.ge.layb0).and.(li.le.layb0f))) then
          if(ib0.eq.0) then
            x=dlog(rs(li))-(ns(li)-1)*dx(li)
            do i=1,ns(li)
              r=dexp(x)
              bra(i,li)=bra(i,li)+r*b0
              brb(i,li)=brb(i,li)+r*b0
              x=x+dx(li)
            end do
          else
            x=dlog(rs(li))-(ns(li)-1)*dx(li)
            do i=1,ns(li)
              r=dexp(x)
              bra(i,li)=r*b0
              brb(i,li)=r*b0
              x=x+dx(li)
            end do
          end if
        end if
c
        if(isigba(li).ne.1) then
          do i=1,ns(li)
            bra(i,li)=-bra(i,li)
          end do
          write(6,'(''  Type A: sign of effective field changed'')')
        end if
        if(isigbb(li).ne.1) then
          do i=1,ns(li)
            brb(i,li)=-brb(i,li)
          end do
          write(6,'(''  Type B: sign of effective field changed'')')
        end if
c
        write(6,'(2x,a3,''('',f5.3,'') - '',a3,''('',f5.3,'')'')')
     >  idpota(li)(1:3),conc(li),idpotb(li)(1:3),1.d0-conc(li)
      end do
      close(20)
c
c
      if(.not.bulk) return
c
c limitation to 3D periodicity:
c
      do li=1,nintfc
        idpotla(li)=idpota(li)
        idpotra(li)=idpota(li)
        idpotlb(li)=idpotb(li)
        idpotrb(li)=idpotb(li)
        dxl(li)=dx(li)
        dxr(li)=dx(li)
        rsl(li)=rs(li)
        rsr(li)=rs(li)
        nsl(li)=ns(li)
        nsr(li)=ns(li)
        zla(li)=za(li)
        zra(li)=za(li)
        zlb(li)=zb(li)
        zrb(li)=zb(li)
        concl(li)=conc(li)
        concr(li)=conc(li)
      enddo
c
      return
      end
