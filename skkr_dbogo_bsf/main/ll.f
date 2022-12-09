c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      integer function ll(ilm) 
!- Returns l, given lm index                                            
! ----------------------------------------------------------------------
!i Inputs:                                                              
!i   ilm   :lm-value                                                    
!o Outputs:                                                             
!o   ll    :l-value                                                     
! ----------------------------------------------------------------------
      implicit none 
! Passed parameters:                                                    
      integer ilm 
! Local parameters:                                                     
      integer nlmax 
      parameter (nlmax=16) 
      integer lla(0:(2*nlmax-1)**2) 
      data lla/01*-1,01*00,03*01,05*02,07*03,09*04,     
     &               11*05,13*06,15*07,17*08,19*09,    
     &               21*10,23*11,25*12,27*13,29*14,   
     &               31*15,33*16,35*17,37*18,39*19,  
     &               41*20,43*21,45*22,47*23,49*24, 
     &               51*25,53*26,55*27,57*28,59*29,61*30/               
      ll = lla(ilm) 
      END                                           
