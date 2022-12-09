c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine wafu(ce,psq,lmax,idpot,v0,vr,br,bopr,dx,ns,rs,
     >                tminv,gz,fz,gj,fj,nuz,indz,iflag,socsc)        
c====================
c
c input:  e - any complex energy
c         lmax  - maximum of angular momentum index
c         v0 - vacuum potential
c         idpot,vr,br,bopr,dx,ns,rs - as in 'readpot'
c output: tminv - inverse of single-site scattering matrix
c         gz,fz   - big and small component of regular radial solution*r
c         nuz - no. of (kap',my') components for (kap,my)
c         indz - selects (kap',my') for (kap,my)
c
      implicit real*8 (a-h,o-z)
c
      include '../param.h'
c
      character*10 idpot
      dimension vr(nrad),br(nrad),bopr(nrad,2)
      complex*16 tminv2(2,2)
      complex*16 tminv(kmymaxp,kmymaxp)
      complex*16 gz2(nrad,2,2),fz2(nrad,2,2)
      complex*16 gj2(nrad,2,2),fj2(nrad,2,2)
      complex*16 gz(nrad,nuzp,kmymaxp),fz(nrad,nuzp,kmymaxp)
      complex*16 gj(nrad,nuzp,kmymaxp),fj(nrad,nuzp,kmymaxp)
      complex*16 fb(0:l1maxp),fn(0:l1maxp),fh(0:l1maxp)
      complex*16 fb1(0:l1maxp),fn1(0:l1maxp),fh1(0:l1maxp)
      complex*16 ce,psq,p,cevac,psqvac,pvac,sk,dk,xk,react,sqrtm1
      dimension nuz(kmymaxp),indz(nuzp,kmymaxp)
      data sqrtm1/(0.d0,1.d0)/
c
c     write(6,'(a10)') idpot
c     write(6,'(2f10.6,i5)') dx,rs,ns
c     write(6,'(4d20.10)') (vr(j),j=1,ns)
c     write(6,'(4d20.10)') (br(j),j=1,ns)
c
      kmax=2*lmax+1
      kmymax=2*(lmax+1)*(lmax+1)
      call czero(tminv,kmymaxp*kmymaxp)
c
      if(idpot.eq.'Vacuum    ') then
c
        p=cdsqrt(psq)
        if(dimag(p).lt.0.d0) p=-p
        cevac=ce-v0
        pvac=cdsqrt(cevac)
        if(dimag(pvac).lt.0.d0) pvac=-pvac
        call csbf(lmax+1,l1maxp,p,rs,fb,fn,fh)
        call csbf(lmax+1,l1maxp,pvac,rs,fb1,fn1,fh1)
        do k=1,kmax
          l=k/2
          if(k.eq.2*l) then
            kap=l
            lb=l-1
            j=2*l-1
          else
            kap=-l-1
            lb=l+1
            j=2*l+1
          end if
          sk=dcmplx(dfloat(l-lb),0.d0)
          xk=sk*cevac*fb1(lb)/pvac
          xk=xk/fb1(l)
          sk=sk*ce/p
          dk=(xk*fb(l)-sk*fb(lb))/(xk*fn(l)-sk*fn(lb))
          react=-dk/p
          do my=-j,j,2
            kapmy=2*kap*kap+kap+(my+1)/2
            tminv(kapmy,kapmy)=1.d0/react+sqrtm1*p
          end do
        end do
c
      else  
c
        do l=0,lmax
          kap1=l
          kap2=-l-1
          do my=-2*l-1,2*l+1,2
c
            call spzwafu(socsc,ce,psq,l,my,vr,br,bopr,dx,ns,rs,
     >                   tminv2,gz2,fz2,gj2,fj2,iflag)
c
            if(iabs(my).eq.2*l+1) then
c
              kapmy=2*kap2*kap2+kap2+(my+1)/2
              tminv(kapmy,kapmy)=tminv2(2,2)
c
              if(iflag.eq.1) then
                nuz(kapmy)=1
                indz(1,kapmy)=kapmy
                do i=1,ns
                  gz(i,1,kapmy)=gz2(i,2,2)
                  fz(i,1,kapmy)=fz2(i,2,2)
                  gj(i,1,kapmy)=gj2(i,2,2)
                  fj(i,1,kapmy)=fj2(i,2,2)
                end do
              end if
c
            else
c
              kapmy1=2*kap1*kap1+kap1+(my+1)/2
              kapmy2=2*kap2*kap2+kap2+(my+1)/2
              tminv(kapmy1,kapmy1)=tminv2(1,1)
              tminv(kapmy1,kapmy2)=tminv2(1,2)
              tminv(kapmy2,kapmy1)=tminv2(2,1)
              tminv(kapmy2,kapmy2)=tminv2(2,2)
c
              if(iflag.eq.1) then
                nuz(kapmy1)=2
                nuz(kapmy2)=2
                indz(1,kapmy1)=kapmy1
                indz(2,kapmy1)=kapmy2
                indz(1,kapmy2)=kapmy2
                indz(2,kapmy2)=kapmy1
                do i=1,ns
                  gz(i,1,kapmy1)=gz2(i,1,1)
                  fz(i,1,kapmy1)=fz2(i,1,1)
                  gz(i,2,kapmy1)=gz2(i,2,1)
                  fz(i,2,kapmy1)=fz2(i,2,1)
                  gz(i,1,kapmy2)=gz2(i,2,2)
                  fz(i,1,kapmy2)=fz2(i,2,2)
                  gz(i,2,kapmy2)=gz2(i,1,2)
                  fz(i,2,kapmy2)=fz2(i,1,2)
                  gj(i,1,kapmy1)=gj2(i,1,1)
                  fj(i,1,kapmy1)=fj2(i,1,1)
                  gj(i,2,kapmy1)=gj2(i,2,1)
                  fj(i,2,kapmy1)=fj2(i,2,1)
                  gj(i,1,kapmy2)=gj2(i,2,2)
                  fj(i,1,kapmy2)=fj2(i,2,2)
                  gj(i,2,kapmy2)=gj2(i,1,2)
                  fj(i,2,kapmy2)=fj2(i,1,2)
                end do
              end if
c
            end if
c
          end do
        end do
c
      end if
c
      return
      end
