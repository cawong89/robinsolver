! =================================================================================================
!                                    FACTORIZATION MODULE
! =================================================================================================
!
! This is a submodule of the Derived Type Module and contains various subroutines bound to the
! class MATRIXDATA, used for computing matrix factorizations.
!
! Dev note: Procedure interfaces are contained in the Derived Type Module, so procedure argument
! list does not need to be repeated in the submodule (F2008 standard).

submodule (derived_type_module) factorization_module

    ! implicit none
    ! private
    
    ! integer,    parameter   :: dp = kind(1.0d0) ! double precision kind specification
    
    ! public  :: numrank, numrankQR

    
    
    
    contains
    
    
    
    
! =================================================================================================    
! SUBROUTINE MTRX_SET
! Initializes a variable of type MATRIXDATA with a dense matrix.
    module procedure mtrx_set
    
        implicit none
        ! class(matrixdata),                  intent(out) :: mtrx
        ! complex(dp),        dimension(:,:), intent(in)  :: A
        
        mtrx%factored = .false. ! Should be redundant
        
        mtrx%m = size(A,1)
        mtrx%n = size(A,2)
        
        allocate(mtrx%A(1:mtrx%m,1:mtrx%n))
        
        mtrx%A = A
    
    end procedure mtrx_set
! =================================================================================================    



    
    
! =================================================================================================    
! SUBROUTINE MTRX_FACTOR
! Produces factorization data according to the chosen factorization mode.
    module procedure mtrx_factor
    
        implicit none
        integer :: info
        
        mtrx%mode = mode
        
        ! Mode = LU factorization
        if (mtrx%mode == 'LU') then
        
            if (.not. allocated(mtrx%ipiv)) then
            
                allocate(mtrx%ipiv(1:min(mtrx%m,mtrx%n)))
                
            else
            
                deallocate(mtrx%ipiv)
                allocate(mtrx%ipiv(1:min(mtrx%m,mtrx%n)))
            
            end if
        
            ! ZGETRF( M, N, A, LDA, IPIV, INFO )
            call zgetrf(mtrx%m, mtrx%n, mtrx%A, mtrx%m, mtrx%ipiv, info)
            
        else
        
            write(*, '(A)') 'Error in MTRX_FACTOR: Factorization mode not supported.'
            stop
        
        end if
        
        mtrx%factored = .true.
        
        
    
    end procedure mtrx_factor
! =================================================================================================    





! =================================================================================================    
! SUBROUTINE MTRX_GETINV
! Explicitly produces the matrix inverse.
    module procedure mtrx_getinv
    
        implicit none
        integer                                         :: info, lwork
        integer,        dimension(:),   allocatable     :: ipiv
        complex(dp),    dimension(:),   allocatable     :: work
        
        
        ! Check dimensions
        if (mtrx%m /= mtrx%n) then
            write(*, '(A)') 'Error in MTRX_GETINV: Matrix must be square.'
            stop
        else if ((size(InvA,1) /= size(InvA,2)) .or. (size(InvA,1) /= mtrx%m)) then
            write(*, '(A)') 'Error in MTRX_GETINV: Dimension mismatch.'
            stop
        end if
        
        if (.not. mtrx%factored) then
            allocate(ipiv(1:mtrx%m))
            lwork = max(64, mtrx%m)
            allocate(work(1:lwork))
            
            InvA = mtrx%A
            
            call zgetrf(mtrx%m, mtrx%n, InvA, mtrx%m, ipiv, info)
            call zgetri(mtrx%m, InvA, mtrx%m, ipiv, work, lwork, info)
            
        else
            if (mtrx%mode == 'LU') then
            
                lwork = max(64, mtrx%m)
                allocate(work(1:lwork))
                
                InvA = mtrx%A
                
                call zgetri(mtrx%m, InvA, mtrx%m, mtrx%ipiv, work, lwork, info)
            
            else
            
                write(*, '(A)') 'Error in MTRX_GETINV: Factorization mode not supported.'
                stop
            
            end if
            
            
        end if
    
    end procedure mtrx_getinv
! =================================================================================================    






! =================================================================================================    
! SUBROUTINE MTRX_MATVEC
! Applies matrix-vector product y = alpha*A*x + beta*y.
    module procedure mtrx_matvec
    
        implicit none
        complex(dp),    dimension(:),   allocatable :: vec_temp
        integer :: j
        
        ! Check dimensions match
        if (.not. ( (size(y) == mtrx%m) .and. (size(x) == mtrx%n) ) ) then
        
            write(*, '(A)') 'Error in MTRX_MATVEC: Vector size mismatch.'
            stop
        
        end if
        
        ! Direct matvec if matrix is unfactored
        if (.not. mtrx%factored) then
        
            !SUBROUTINE ZGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
            
            call zgemv('n', mtrx%m, mtrx%n, alpha, mtrx%A, mtrx%m, x, 1, beta, y, 1)
        
        else
        
            if (mtrx%mode == 'LU') then
            
                ! Check matrix is square
                if (mtrx%m /= mtrx% n) then
                    write(*, '(A)') 'Error in MTRX_MATVEC: Matvec cannot be applied for nonsquare LU'
                    stop
                end if
                
                allocate(vec_temp(1:size(x)))
                
                vec_temp = x
            
                ! SUBROUTINE ZTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
                
                ! multiply by U factor
                call ztrmv('U', 'n', 'n', mtrx%m, mtrx%A, mtrx%m, vec_temp, 1)

                
                ! multiply by L factor
                call ztrmv('L', 'n', 'u', mtrx%m, mtrx%A, mtrx%m, vec_temp, 1)
                
                ! pivot
                call zlaswp(1, vec_temp, mtrx%m, 1, mtrx%m, mtrx%ipiv, 1)

                ! zaxpy
                y = beta*y
                call zaxpy(mtrx%m, alpha, vec_temp, 1, y, 1)
                
            else

                write(*, '(A)') 'Error in MTRX_MATVEC: Factorization mode not supported.'
                stop
            
            end if
            
        end if
    
    end procedure mtrx_matvec
! =================================================================================================





! =================================================================================================
! SUBROUTINE MTRX_MATMAT
! Computes the matrix-matrix product C = alpha * A * B + beta * C.
    module procedure mtrx_matmat
            
        implicit none
        
        complex(dp),    dimension(:,:),     allocatable :: mat_temp
        
        integer :: k
        
        k = size(B,2)
        
        ! Check dimensions
        if (.not. (size(B,1) == mtrx%n) ) then
            write(*, '(A)') 'Error in MTRX_MATMAT: Vector size mismatch.'
            stop
        else if (.not. (size(C,1) == mtrx%m) ) then
            write(*, '(A)') 'Error in MTRX_MATMAT: Vector size mismatch.'
            stop
        else if (.not. (size(C,2) == k) ) then
            write(*, '(A)') 'Error in MTRX_MATMAT: Vector size mismatch.'
            stop
        end if
        
        if (.not. mtrx%factored) then
        
            !SUBROUTINE ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
            
            call zgemm('n','n', mtrx%m, mtrx%n, k, alpha, mtrx%A, mtrx%m, &
                        B, mtrx%n, beta, C, mtrx%m)
        
        else
        
            if (mtrx%mode == 'LU') then
            
                ! allocate temporary array for U*B
                allocate(mat_temp(1:mtrx%m, 1:k))
                mat_temp = B
                
                !SUBROUTINE ZTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
                
                call ztrmm('L', 'U', 'n', 'n', mtrx%m, k, cmplx(1.0,0.0,dp), mtrx%A, &
                            mtrx%m, mat_temp, mtrx%m)
                            
                call ztrmm('L', 'L', 'n', 'u', mtrx%m, k, alpha, mtrx%A, &
                            mtrx%m, mat_temp, mtrx%m)
                            
                ! pivot
                call zlaswp(k, mat_temp, mtrx%m, 1, mtrx%m, mtrx%ipiv, 1)
                
                C = mat_temp + beta*C
                
            
            else
                write(*, '(A)') 'Error in MTRX_MATMAT: Factorization mode not supported.'
                stop 
            end if
        
        end if
    
    end procedure mtrx_matmat
! =================================================================================================






! =================================================================================================    
! SUBROUTINE MTRX_INVVEC
! If X is not present, then this routine computes B <- beta*inv(A)*B.
! If X is present, then computes X <- alpha*X + beta*inv(A)*B
    module procedure mtrx_invvec
    
        implicit none
        complex(dp),    dimension(:,:), allocatable     :: mat_temp, x_temp    
        integer                                         :: nrhs, info
        integer,        dimension(:),   allocatable     :: ipiv
        character(len=1)                                :: LR
        
        ! Check dimensions match
        if (.not. (size(b,1) == mtrx%m )) then
            write(*, '(A)') 'Error in MTRX_INVVEC: Vector size mismatch.'
            stop
        else if (mtrx%m /= mtrx%n) then
            write(*, '(A)') 'Error in MTRX_INVVEC: Matrix must be square.'
        else if (present(x)) then
            if (.not. (size(x,2)==size(b,2)) ) then
                write(*, '(A)') 'Error in MTRX_INVVEC: Number of RHS mismatch.'
                stop
            else if (.not. (size(x,1) == mtrx%n)) then
                write(*, '(A)') 'Error in MTRX_INVVEC: Vector size mismatch.'
                stop
            end if
        end if
        
        nrhs = size(b,2)
        
        if (.not. mtrx%factored) then
        
            if (present(x)) then
                allocate(x_temp(1:size(x,1), 1:size(x,2)))
                x_temp = b
            end if
            
            ! Copy matrix to temporary array for LU factorization
            allocate(mat_temp(1:mtrx%m, 1:mtrx%n))
            mat_temp = mtrx%A
            
            allocate(ipiv(1:mtrx%m))
            
            if (present(x)) then
                call zgesv(mtrx%m, nrhs, mat_temp, mtrx%m, ipiv, x_temp, mtrx%m, info)
                x = alpha*x + beta*x_temp
            else
                call zgesv(mtrx%m, nrhs, mat_temp, mtrx%m, ipiv, b, mtrx%m, info)
                b = beta*b
            end if

        
        else
        
            if (mtrx%mode == 'LU') then
            
                if (present(x)) then
                    allocate(X_temp(1:size(x,1), 1:size(x,2)))
                    X_temp = b
                    
                    call zgetrs('n', mtrx%m, nrhs, mtrx%A, mtrx%m, mtrx%ipiv, x_temp, mtrx%m, info)
                    x = alpha*x + beta*x_temp
                else
                    call zgetrs('n', mtrx%m, nrhs, mtrx%A, mtrx%m, mtrx%ipiv, b, mtrx%m, info)
                    b = beta*b
                end if
                
                
            
            else

                write(*, '(A)') 'Error in MTRX_INVVEC: Factorization mode not supported.'
                stop
            
            end if
        
        end if
        
        
    
    end procedure mtrx_invvec
! =================================================================================================    

    


    
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
        
        integer :: j   
        
        if (tol <= 0.0d0) then
            write(*, '(A)') 'Error in NUMRANK: Requested tolerance must be positive.'
        end if
        
          
        
        n = size(A,1)
        
        allocate(sigma(1:n))
        
        lwork = 3*n
        
        allocate(work(1:lwork))
        allocate(rwork(1:5*n))

        ! SUBROUTINE ZGESVD(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK, INFO)
       
        call zgesvd('n', 'n', n, n, A, n, sigma, U, 0, VT, 0, work, lwork, rwork, info)

        j = 1
        
        do while ( (j <= n) .and. (sigma(j)/sigma(1) > tol) )
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


    
    

end submodule factorization_module