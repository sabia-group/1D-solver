module math 
IMPLICIT NONE
contains
!----------------------------------------------------------------
!----------------------------------------------------------------
!https://web.stanford.edu/class/me200c/tutorial_90/08_subprograms.html


FUNCTION Factorial(n)  RESULT(fact)

REAL(kind=8) :: Fact
INTEGER(kind=8)            :: i
INTEGER(kind=8),INTENT(IN) :: n

   fact = 1.0
   DO  i= 1,n
      fact = fact * real(i)
   END DO
END FUNCTION 

!----------------------------------------------------------------
!RECURSIVE FUNCTION Factorial2(n)  RESULT(Fact)
!
!INTEGER(kind=8) :: Fact2
!INTEGER(kind=8),INTENT(IN) :: n
!
!IF (n == 0) THEN
!   Fact2 = 1.0
!ELSE
!   Fact2 = real(n) * Factorial2(n-1)
!END IF
!END FUNCTION 
!----------------------------------------------------------------
pure function delta(n,m)
INTEGER(KIND=8), INTENT(IN) :: n,m
REAL(KIND=8)      :: delta


IF (n == m) THEN
    delta = 1
ELSE
   delta = 0
END IF
end function
!----------------------------------------------------------------
pure function mysqrt(m)
REAL(KIND=8), INTENT(IN) :: m
REAL(KIND=8)      :: mysqrt

IF (m .GT.0) THEN
    mysqrt = SQRT(m)
ELSE
   mysqrt = 0
END IF
end function
!----------------------------------------------------------------
subroutine diag(a,lambda,n,verbose)

!Interface to lapack dsyev, only for symmetric matrix
! Compute the eigenvalues and eigenvectors.

  implicit none
  integer ( kind = 4 )   n,lwork,verbose 
  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) lambda(n)
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: work 

  integer ( kind = 4 ) info
  character jobz
  character uplo


  IF (verbose .GE. 1) THEN
     write(6,*) '---------------- Diagonalizing with dsyev (LAPACK) ---------------------' 
  ENDIF

  lwork = 3 * n - 1
  ALLOCATE( work(lwork))  
  jobz = 'V'
  uplo = 'U'
 
  IF (verbose .EQ. 2) THEN
  write(6,*) 'a'
  call r8mat_print ( n, n, a, '  The matrix is:' )
  ENDIF

  call dsyev ( jobz, uplo, n, a, n, lambda, work, lwork, info )

  if ( info /= 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  DSYEV returned nonzero INFO = ', info

  end if

  IF (verbose .EQ. 2) THEN 
    call r8vec_print ( n, lambda, '  The eigenvalues:' )
    call r8mat_print ( n, n, a, '  The eigenvector matrix:' )
  ENDIF


  return
end
!----------------------------------------------------------------
!----------------------------------------------------------------
subroutine diag2(a,eigv,n,verbose)
implicit none
REAL(KIND=8),DIMENSION(n,n)  :: a,aux_a 
!REAL(KIND=8), DIMENSION(np) :: eigv, dp
INTEGER ( kind =4 ) :: n,verbose ,i
INTEGER,DIMENSION(n) :: jndx
REAL(KIND=8), DIMENSION(n) :: DP
REAL(KIND=8),DIMENSION(n)     ::  eigv,aux_eigv
!REAL(KIND = 8 ) lambda(n)

IF (verbose .GE. 1) THEN
write(6,*) '--------------- Diagonalizing using numeral recipies --------------------' 
END IF

IF (verbose .EQ. 2) THEN 
call r8mat_print ( n, n, a, '  The matrix A:' )
ENDIF



CALL TRED2(a,n  ,n,EIGV,DP)
CALL TQLI (EIGV,DP,n,n,a)

aux_eigv = eigv
aux_a = a
!CALL SORTING (n,aux_eigv,JNDX)

DO i=1,n
  eigv(i) = aux_eigv(jndx(i))
  a(1:n,i) = aux_a(1:n,jndx(i))
ENDDO
write(6,*) 'JNDX',(JNDX(I),I=1,n)

IF (verbose .EQ. 2) THEN 
  call r8vec_print ( n, eigv, '  The eigenvalues:' )
  call r8mat_print ( n, n, a, '  The eigenvector matrix:' )
ENDIF

end subroutine
!----------------------------------------------------------------
!----------------------------------------------------------------
subroutine r8mat_print ( m, n, a, title )

!    Input, integer M, the number of rows in A.
!    Input, integer N, the number of columns in A.
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!    Input, character ( len = * ) TITLE, a title to be printed.

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows in A.
!
!    Input, integer N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
!----------------------------------------------------------------
!----------------------------------------------------------------
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

    end do

  end do

  write ( *, '(a)' ) ' '

  return
end
!------------------------------------------------------------------
!------------------------------------------------------------------
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC. 
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
! 
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
! 
!    Input, character ( len = * ) TITLE, an optional title.
! 
  implicit none
  
  integer ( kind = 4 ) n
  
  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title
  
  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if
  
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,2x,g16.8)' ) i, a(i)
  end do
  
  return
end
!------------------------------------------------------------------
!------------------------------------------------------------------
subroutine Fourier (n,nres,func, ft,ndim) 
INTEGER(kind=8) :: n,nres
INTEGER(kind=8) :: NDIM
COMPLEX(kind=8)    :: func(n*nres),ft(NDIM),i,zero
INTEGER(kind=8) :: j,jj
real(KIND=8)    :: t 
  zero = (0.0,0.0)
  i    = (0.0,1.0)

  DO j =1,NDIM
  ft(j)=zero
  ENDDO

  DO  j = 1,NDIM
  DO jj = 1,n*nres
       t = real(jj-1)/real(nres)     
       ft(j) = ft(j)+func(jj)*ZEXP(-i*real(j)*t )
  END DO
  END DO
  
  return
end
!------------------------------------------------------------------
!------------------------------------------------------------------
!pure function checkH(H)
!CHECK if H is symmetric
!What else
!end

!------------------------------------------------------------------
!------------------------------------------------------------------

end module math
