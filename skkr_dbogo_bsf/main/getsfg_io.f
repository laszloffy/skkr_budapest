c     Name of Copyright owners: Gabor Csire, Andras Laszloffy, Bendeguz Nyari, Laszlo Szunyogh, and Balazs Ujfalussy
      subroutine getsfg_io(sfgin,sfgout,ndim,mdim,iwrite,ksgf)
c
c---------------------------------------
c     sfgout=sfgin
c     iwrite =-1 -> read  sfgin from file
c     iwrite = 1 -> write sfgout in file
c     else no i/o operation occurs
c---------------------------------------  
c
      implicit real*8 (a-h,o-z)
      complex*16 sfgin(mdim,mdim),sfgout(mdim,mdim)
      integer nn,k1,k2,ksgf
      data tol/1.0d-8/
c
      nn=mdim
c
      if(iwrite.eq.-1) then
        read(ksgf) ((sfgin(k1,k2),k1=1,nn),k2=1,nn)
      end if
c
      call repl(sfgout,sfgin,ndim,mdim)
c
      if(iwrite.eq.1) then
        write(ksgf) ((sfgout(k1,k2),k1=1,nn),k2=1,nn)
      end if
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
