! =================================================================================================
!                                    FACTORIZATION MODULE
! =================================================================================================
!
! This modules contains various subroutines for computing matrix factorizations.

module factorization_module

    implicit none
    private
    
    integer,    parameter   :: dp = kind(1.0d0) ! double precision kind specification
    
    public  :: matrixdata, numrank, numrankQR
    
    
    
    
! =================================================================================================
! TYPE MATRIXDATA
! Contains all the relevant data for a matrix, including all relevant factorization quantities as
! as well as subroutines for inverting or applying the matrix.
!
!   Attributes:
!
!   factored            --- logical; determines if factorization data has been computed
!   mode                --- type of factorization, such as 'lu, 'svd', 'hss'
!   m,n                 --- size of original matrix is m x n
!   A                   --- m x n matrix; a matrix with the same size as the original matrix
!   rk                  --- integer; numerical rank of matrix (subject to some tolerance)
!   ipiv                --- pivot data; used in LU factorization
!
!
!   Type-bound Procedures:
!
!   factor              --- constructs factorization according to mode
!   matvec              --- applies the factored matrix to a vector
!   invvec              --- applies the inverse matrix to a vector (i.e. solves linear system)
    
    type matrixdata
    
        logical                                         :: factored = .false.
        character(len=3)                                :: mode
        integer                                         :: m,n = 0
        complex(dp),    dimension(:,:), allocatable     :: A
        
        integer                                         :: rk = 0

        integer,    dimension(:),   allocatable         :: ipiv
        
        procedure       :: factor   => mtrx_factor
        procedure       :: matvec   => mtrx_matvec
        procedure       :: invvec   => mtrx_invvec
        

    end type matrixdata
! =================================================================================================    
    
    
    
    
    
    contains
    
    
! =================================================================================================
! SUBROUTINE NUMRANK
! Computes the numerical rank of a dense, complex-double-valued matrix subject to some relative
! tolerance value. This is NOT intended to be used for large matrices or matrices of very low rank,
! as the numerical rank is computed by computed the SVD of the entire matrix. For efficient
! numerical rank calculations (but less exact), use pivoted QR techniques.
    subroutine numrank(rank, A, tol)
    
        implicit none
        
        integer,                        intent(out)     :: rank
        complex(dp),    dimension(:,:), intent(in)      :: A
        real(dp),                       intent(in)      :: tol
        
        real(dp),       dimension(:),   allocatable :: sigma
        complex(dp),    dimension(:),   allocatable :: work
        real(dp),       dimension(:),   allocatable :: rwork
        integer                                     :: n, lwork, info
        complex(dp),    dimension(:),   allocatable :: U, VT
        
        if (tol <= 0.0d0) then
            write(*, '(A)') 'Error in NUMRANK: Requested tolerance must be positive.'
        end if
        
        integer :: j             
        
        n = size(A,1)
        
        allocate(sigma(1:n))
        
        lwork = 3*n
        
        allocate(work(1:lwork)
        allocate(rwork(1:5*n))

        ! SUBROUTINE ZGESVD(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK, INFO)
       
        call zgesvd('n', 'n', n, n, A, n, sigma, U, 0, VT, 0, work, lwork, rwork, info)

        j = 1
        
        do while ( (j <= n) && (sigma(j)/sigma(1) > tol) )
            j = j + 1
        end do
        
        rank = j - 1
    
    

    end subroutine numrank
! =================================================================================================




! =================================================================================================
! SUBROUTINE NUMRANKQR
! Computes numerical rank of a dense matrix A using pivoted QR factorization.
! Dev note: For sparse A use SPQR.
    subroutine numrankqr(rank, A, tol)
    
        implicit none
        
        integer,                        intent(out)     :: rank
        complex(dp),    dimension(:,:), intent(in)      :: A
        real(dp),                       intent(in)      :: tol
        

    end subroutine numrankqr
! =================================================================================================


    
    

end module factorization_module