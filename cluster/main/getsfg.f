c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine getsfg(sfgin,sfgout,ndim,mdim,iwrite)
c
c---------------------------------------
c     iwrite = 0 -> store sfg
c     iwrite = 1 -> write sfg in file
c     iwrite =-1 -> read  sfg from file
c---------------------------------------  
c
      implicit real*8 (a-h,o-z)
      complex*16 sfgin(mdim,mdim),sfgout(mdim,mdim)
      data tol/1.0d-8/
c
ccc   if(iwrite.eq.-1) read(*,*) sfgin 
c
      call repl(sfgout,sfgin,ndim,mdim)
c
ccc   if(iwrite.eq.1) write(*,*) sfgout
c
ccccccccccc
ccc   WRITE(6,*) 'getsfg       ',ndim,mdim,tol
ccc   call outmat1(sfgin,ndim,ndim,mdim,tol,66)
ccc   WRITE(6,*) 'getsfg       ',ndim
ccc   call outmat1(sfgout,ndim,ndim,mdim,tol,67)
ccc   STOP
ccccccccccc
      return
      end
