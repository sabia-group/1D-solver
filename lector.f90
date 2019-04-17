      SUBROUTINE LECTOR(xi,xf,dx,nsize,V)
      use constants
      IMPLICIT NONE
        
      integer(kind=8) :: nsize,n
      real(kind=8)    :: xi ,xf, dx,aux
      !real(kind=8), dimension(:), ALLOCATABLE :: V
      real(kind=8), dimension(NDIM) :: V
      save

      V(:) =0.d0
      OPEN (5,FILE='grid.dat')

      READ (5,*) nsize,xi,xf,dx
      READ (5,*)
      !ALLOCATE( V(nsize) )

      
      do n=1,nsize 
      READ (5,*) aux,V(n)
      enddo 
      CLOSE(5)

      V(:)=V(:)*eV2au
      xi = xi *A2au
      xf = xf *A2au
      dx = dx *A2au
      
      RETURN
      END

