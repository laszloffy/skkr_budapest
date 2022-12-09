c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
       subroutine rotvec(r,rp,vecn,phi,sign)
c=======================
c
       implicit real*8(a-h,o-z)
       integer sign
       dimension r(3),rp(3),vecn(3),cross(3)
       data tiny/1.0d-8/
c
       if(dabs(phi).lt.tiny) then
         rp(1)=r(1)
         rp(2)=r(2)
         rp(3)=r(3)
         return
       end if
c
       cosphi=dcos(phi)
       sinphi=dsin(phi)
       if(sign.lt.0) sinphi=-sinphi
c
       cross(1)=vecn(2)*r(3)-vecn(3)*r(2)
       cross(2)=vecn(3)*r(1)-vecn(1)*r(3)
       cross(3)=vecn(1)*r(2)-vecn(2)*r(1)
       dot=vecn(1)*r(1)+vecn(2)*r(2)+vecn(3)*r(3)
c
       rp(1)=cosphi*r(1)+sinphi*cross(1)+(1.0d0-cosphi)*dot*vecn(1)
       rp(2)=cosphi*r(2)+sinphi*cross(2)+(1.0d0-cosphi)*dot*vecn(2)
       rp(3)=cosphi*r(3)+sinphi*cross(3)+(1.0d0-cosphi)*dot*vecn(3)
c
       return
       end
