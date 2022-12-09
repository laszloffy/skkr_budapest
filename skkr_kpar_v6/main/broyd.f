c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
C***********************************************************************
c***************************<  BROYD  >*********************************
c***********************************************************************
      subroutine broyd(nm,fm,beta,m,l,nmme,fmme,df,deltan,deltaf,
     >                 alfa,win,n1,n2,n3,file1,file2,file3)
c
cx   called from: mixio                                           
cx   calls      : trns, dott, fo1aaf(NAGlib), scl                
cx   function : modified broyden's second method  mixing scheme
cx   reference: D.D. JOHNSON PRB 38, 12 807, 1988                 
cx              and also chapter 6 thesis                        
c
      implicit real*8 (a-h,o-z)
c     parameter( n1=21, n2=22, n3=23, itm=20 )                  
c laszlo
      parameter( itm=20 )                  
c
      character*5 file1,file2,file3
      real*8     nm(*),fm(*),nmme(*),fmme(*),df(*)             
     +,          deltan(*),deltaf(*)                           
      real*8     bkn(itm,itm),w(itm),fkm(itm),gm(itm)          
     +,          bkni(itm,itm)                                
      dimension  isc(itm)
c
      if(m.gt.itm) stop ' exceeded dimension in broyd '
c
C******                                                         
C     FILE N1 : /DELTAF(I)>                                     
C     FILE N2 : /U(I)>;                             VGL(13B)    
C     FILE N3 : /N(M)>,/F(M); <DELTAF(I)/DELTAF(J)> VGL(13A)   
C******                                                       
c     open(n1,file=file1,status='unknown')
c     open(n2,file=file2,status='unknown')
c     open(n3,file=file3,status='unknown')
c
      open(n1,file=file1,form='unformatted',status='unknown')
      open(n2,file=file2,form='unformatted',status='unknown')
      open(n3,file=file3,form='unformatted',status='unknown')
      rewind n1                                              
      rewind n2                                             
      rewind n3                                            
C                                                         
      w0 = 0.005D0                                         
c     W(M)= BETA                                        
c     IF ( BETA .LE. 1.0d0 ) THEN                        
c       PRINT*, ' FIRST ITERATION OF THIS SERIES '    
c       ALFA = BETA                                  
c       IF ( M .GT. 1 ) READ*, W(M)                 
c     ENDIF                                       
c simon
c      if(m.gt.1)then
c      open(1,file='w_in')
c      read(1,*)alfa,w(m)
c      close(1)
c      endif
c laszlo
       w(m)=win
C                                                
c     PRINT1, ALFA, BETA, W(M) ,m
c 1   FORMAT( ' ALFA= ',F5.3,' BETA= ',F8.3,' W(M)= ',F8.3,' m ',i5)
C                                                          
      if ( m .eq. 1 ) then                                
c     WRITE(N3,'(1pd25.16)') (NM(i),i=1,l),(FM(i),i=1,l)
      write(n3) (nm(i),i=1,l),(fm(i),i=1,l)
        call trns( nm, fm, alfa, l )                    
        return                                         
      endif                                           
C                                                    
c     read(N3,'(1pd25.16)') (NMME(i),i=1,l),(FMME(i),i=1,l)
      read(n3) (nmme(i),i=1,l),(fmme(i),i=1,l)
C                                                  
      do 10 i=1,l                                 
  10  df(i) = fm(i) - fmme(i)                    
      call dott( df, df, fnorm, l )             
      fnorm = 1.0D0/sqrt( fnorm )                
      do 20 i=1,l                             
      deltan(i) = fnorm*( nm(i) - nmme(i) )  
  20  deltaf(i) = fnorm*df(i)               
C                                          
C     STEL INVERSE BETA(K,N) OP    VGL. (13A)
C                                           
      bkni(m-1,m-1) = w0*w0 + w(m)*w(m)   
C                                        
      do 30 j=1, m-2                    
C       read(n3,'(1pd25.16)') ( BKNI(I,J), I= 1,J ) 
        read(n3) ( bkni(i,j), i= 1,j ) 
        w(j+1) = sqrt( bkni(j,j) - w0*w0 )
  30  continue                           
C                                       
      call dott( deltaf, fm, fkm(m-1), l )
      do 40 i= 1,m-2                     
C      read(N1,'(1pd25.16)') (DF(J), J=1,L)          
       read(n1) (df(j), j=1,l)          
       call dott( df, fm, fkm(i), l )  
       call dott( df, deltaf, bkni(i,m-1), l )
       bkni(i,m-1) = bkni(i,m-1)*w(m)*w(i+1) 
  40  continue                       
c     WRITE(N1,'(1pd25.16)') (DELTAF(J), J=1,L)  
      write(n1) (deltaf(j), j=1,l)  
C                                  
      rewind n3                   
c     WRITE(N3,'(1pd25.16)') (NM(i),i=1,l),(FM(i),i=1,l)
      write(n3) (nm(i),i=1,l),(fm(i),i=1,l)
      do 50 j=1,m-1             
c       WRITE(N3,'(1pd25.16)') (BKNI(I,J), I=1,J)
        write(n3) (bkni(i,j), i=1,j)
  50  continue                     
C                                 
      do 60 i= 1, m-2            
       do 60 j = i+1, m-1       
        bkni(j,i) = bkni(i,j)  
  60  continue                
C                            
C     BEPAAL DE WERKELIJKE BKN
C                            
      if ( m .eq. 2 ) then  
        bkn(1,1) = 1.0D0/bkni(1,1)
      else                     
c simon
       do 72 i=1,m-1
       do 73 j=1,m-1
   73  bkn(j,i)=0.0d0
   72  bkn(i,i)=1.0d0
       call gelim(bkni,isc,itm,m-1,1.0d-14)
       do 71 i=1,m-1
   71  call subs(bkni,isc,bkn(1,i),itm,m-1)
c       NERR = 0              
c       CALL F01AAF( BKNI, ITM, M-1, BKN, ITM, GM, NERR )
c       IF ( NERR .NE. 0 ) STOP ' INVERSION ERR BROYD ' 
      endif                                            
C                                                     
      do 70 i = 1,m-2                                
        do 70 j = i+1,m-1                           
         bkn(i,j) = 0.5D0 * ( bkn(i,j) + bkn(j,i) ) 
  70     bkn(j,i) =  bkn(i,j)                    
C                                               
C     STEL DE GAMMA(M,L)'S OP VGL. (15B)       
      do 80 n= 1,m-1                         
       gm(n) = 0.0D0                          
       do 80 k= 1,m-1                      
        gm(n) = gm(n) + w(k+1) * bkn(n,k) * fkm(k)
  80  continue                                   
C                                               
C     STEL DE SOM IN VLG. (15A) OP             
      call trns( deltan, deltaf, alfa, l )    
      call scl( deltaf, 0.0D0, l )            
      do 90 n=1,m-2                        
c      READ(N2,'(1pd25.16)') (DF(J), J=1,L)            
       read(n2) (df(j), j=1,l)            
       call scl( df, w(n+1)*gm(n), l )   
       call trns( deltaf, df, 1.0D0, l )  
  90  continue                         
c     WRITE(N2,'(1pd25.16)') (DELTAN(J), J=1,L)    
      write(n2) (deltan(j), j=1,l)    
C                                    
      call trns( deltaf, deltan, w(m)*gm(m-1), l )
      call trns( nm, fm, alfa, l )               
      call trns( nm, deltaf, -1.0D0, l )          
c     write(6,*)' exiting broyd'
      return                                   
      end
c***********************************************************************
c***************************<  TRNS  >**********************************
c***********************************************************************
      subroutine trns(a,b,c,n)                  
      implicit real*8 (a-h,o-z)
cx   called from: dmixp, broyd                 
cx   calls      : --                          
cx   function : add the vector B scaled by c to A. vector lenghts are N
cx   reference: --                                                    
      dimension a(n),b(n)
      do 1 i=1,n                                                    
 1    a(i)=a(i)+c*b(i)                                             
      return                                                      
      end                                                        
c***********************************************************************
c***************************<  GELIM  >*********************************
c***********************************************************************
      subroutine gelim(ar,nt,np,n,emach)
      implicit real*8(a-h,o-z)
      dimension ar(np,n),nt(n)
      if(n.lt.2) go to 15
      do 12 ii=2,n
      i=ii-1
      yrr=ar(i,i)
      in=i
      do 4 j=ii,n
      if(abs(yrr)-abs(ar(j,i)))3,4,4
3     yrr=ar(j,i)
      in=j
4     continue
      nt(i)=in
      if(in-i)5,7,5
5     do 6 j=i,n
      dum=ar(i,j)
      ar(i,j)=ar(in,j)
6     ar(in,j)=dum
7     if(abs(yrr)-emach)1,1,8
1     ar(i,i)=emach*emach
      go to 12
8     do 11 j=ii,n
      if(abs(ar(j,i))-emach)11,11,9
9     ar(j,i)=ar(j,i)/yrr
      do 10 k=ii,n
10    ar(j,k)=ar(j,k)-ar(i,k)*ar(j,i)
11    continue
12    continue
15    if(abs(ar(n,n))-emach)13,13,14
13    ar(n,n)=emach*emach
14    continue
      return
      end
c***********************************************************************
c***************************<  SUBS  >**********************************
c***********************************************************************
      subroutine subs(ar,nt,xr,np,n)
      implicit real*8(a-h,o-z)
      dimension ar(np,n),xr(n),nt(n)
      if(n.lt.2) go to 18
      do 20 ii=2,n
      i=ii-1
      if(nt(i)-i)16,17,16
16    in=nt(i)
      dum=xr(in)
      xr(in)=xr(i)
      xr(i)=dum
17    do 19 j=ii,n
      xr(j)=xr(j)-ar(j,i)*xr(i)
19    continue
20    continue
c        back substitution
18    do 25 ii=1,n
      i=n-ii+1
      ij=i+1
      if(i-n)21,25,21
21    do 22 j=ij,n
22    xr(i)=xr(i)-ar(i,j)*xr(j)
25    xr(i)=xr(i)/ar(i,i)
      return
      end
c***********************************************************************
c***************************<  SCL  >***********************************
c***********************************************************************
      subroutine scl(a,c,n)                               
      implicit real*8 (a-h,o-z)
cx   called from: dmixp, broyd                           
cx   calls      : --                                    
cx   function : SCALES THE N-DIMENSIONAL VECTOR A BY C 
cx   reference: --                                    
      real*8 a(n),c                                  
      do 1 i=1,n                                    
 1    a(i)=c*a(i)                                  
      return                                      
      end                                        
c***********************************************************************
c***************************<  DOTT  >**********************************
c***********************************************************************
      subroutine dott(a,b,c,n)
cx   called from: mixio, dmixp, broyd              
cx   calls      : --                              
cx   function : inproduct c between the vectors A and B of length N
cx   reference: --                                                
      implicit real*8 (a-h,o-z)
      dimension a(n),b(n)
      c=0.                                                     
      do 1 i=1,n                                              
 1    c = c+a(i)*b(i)                                        
      return                                                
      end                                                  
