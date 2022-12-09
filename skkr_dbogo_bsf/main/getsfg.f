c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine getsfg(sfgin,sfgout,ndim,mdim)
c
      implicit real*8 (a-h,o-z)
      complex*16 sfgin(mdim,mdim),sfgout(mdim,mdim)
      data tol/1.0d-8/
c
      call repl(sfgout,sfgin,ndim,mdim)
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
