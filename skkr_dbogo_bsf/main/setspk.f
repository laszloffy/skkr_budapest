c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
       subroutine setspk(lat,iset,nspk)
c
       if(lat.eq.1) then         
         iset=min(iset,8)
         if(iset.eq.1) nspk=1
         if(iset.eq.2) nspk=2
         if(iset.eq.3) nspk=4
         if(iset.eq.4) nspk=8
         if(iset.eq.5) nspk=16
         if(iset.eq.6) nspk=32
         if(iset.eq.7) nspk=64
         if(iset.eq.8) nspk=128
         if(iset.eq.9) nspk=256
       elseif(lat.eq.3.or.lat.eq.2) then         
         iset=min(iset,4)
         if(iset.eq.1) nspk=1
         if(iset.eq.2) nspk=4
         if(iset.eq.3) nspk=16
         if(iset.eq.4) nspk=64
         if(iset.eq.5) nspk=256
       elseif(lat.eq.4) then
         iset=min(iset,5)
         if(iset.eq.1) nspk=1
         if(iset.eq.2) nspk=3
         if(iset.eq.3) nspk=10
         if(iset.eq.4) nspk=36
         if(iset.eq.5) nspk=136
         if(iset.eq.6) nspk=528
       elseif(lat.eq.5) then
         iset=min(iset,6)
         if(iset.eq.1) nspk=1
         if(iset.eq.2) nspk=3
         if(iset.eq.3) nspk=6
         if(iset.eq.4) nspk=18
         if(iset.eq.5) nspk=45
         if(iset.eq.6) nspk=135
         if(iset.eq.7) nspk=378
       end if
c
       return
       end
