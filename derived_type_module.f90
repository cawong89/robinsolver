! =================================================================================================
!                                     DERIVED TYPE MODULE
! =================================================================================================
!
! This module contains all the commonly used derived types for the elliptic solver using
! hierarchical merging of Robin operators. This includes the data structures used for representing
! a 2D or 3D elliptic problem on a box.

module derived_type_module

    implicit none
    private
    
    integer,    parameter   :: dp = kind(1.0d0) ! double precision kind specification

    public  :: coeff_proc, scalar_proc, sol_proc, elliptic_operator, robin_tree_node, &
                elliptic_rhs, bndry_rhs, bndry_operator, bndry_vector, solve_opts
                
    public  :: matrixdata
                
                
                
                


! =================================================================================================
! Functions interface
! 1. COEFF_PROC: matrix-valued functions that serve as the highest-order coefficients of an
! elliptic operator written in divergence form. Only the upper 2x2 portion of the matrix is read if
! the problem is set to d = 2.
! 2. SCALAR_PROC: complex-valued functions on the domain.
! 3. SOL_PROC: Similar to SCALAR_PROC, but when the function must have a 3-dimensional argument.

    abstract interface
        pure function coeff_proc(x) result(c)
            import dp
            real(dp),       dimension(:),   intent(in)      :: x
            complex(dp),    dimension(1:3,1:3)              :: c
        end function

        pure function scalar_proc(x) result(r)
            import dp
            real(dp),       dimension(:),   intent(in)      :: x
            complex(dp)                                     :: r
        end function
        
        pure function sol_proc(x) result(r)
            import dp
            real(dp),       dimension(1:3), intent(in)      :: x
            complex(dp)                                     :: r
        end function
    end interface

! =================================================================================================





! =================================================================================================
! This type encodes a discretization-free representation of an elliptic operator on a 2D or 3D box
!
!   d           --- Euclidean dimension; must be d = 2 or 3
!   domain      --- Format [x1 y1 z1 ; x2 y2 z2] for coordinates of box vertices
!   coeff       --- Procedure pointer to matrix-valued divergence-form elliptic coefficients
!   k           --- Procedure pointer to complex-valued potential function
!   robin_coeff --- Format [a1 a2; ... ]^T for the Robin condition* on each face: (a1*u + a2*u_n)
!
! 	Notes:
!	* Robin condition must be ordered as X1,X2,Y1,Y2,Z1,Z2, where e.g. X2 is the face furthest in
!		the x direction. Thus only the first four entries will be read if d = 2.
    type elliptic_operator
    
        integer                         ::  d

        real(dp),   dimension(1:2,1:3)  ::  domain

        procedure(coeff_proc),  pointer,    nopass  ::  coeff => NULL()
        procedure(scalar_proc), pointer,    nopass  ::  k     => NULL()

        complex(dp),    dimension(1:6,1:2)  :: robin_coeff


    end type elliptic_operator
! =================================================================================================





! =================================================================================================
    type bndry_rhs
    
        procedure(scalar_proc), pointer,    nopass  :: g_face   => NULL()

    end type bndry_rhs
! =================================================================================================





! =================================================================================================
! This type encodes the right-hand-sides for a linear elliptic system.
!
!
!	d			---	Euclidean dimension; must be d = 2 or 3
!
    type elliptic_rhs
    
        integer                         :: d            ! Euclidean dimension; must be d = 2 or 3
        procedure(scalar_proc), pointer,    nopass  :: f! RHS Scalar function for interior operator
        type(bndry_rhs),    dimension(1:6)  :: g        ! 

    end type elliptic_rhs



! =================================================================================================




! =================================================================================================
! Derived type for a boundary operator between two (possibly different) faces of a box

    type bndry_operator
    
        complex(dp),    allocatable,    dimension(:,:)  :: mat
    
    end type bndry_operator

! =================================================================================================




! =================================================================================================
! Derived type for a vector supported on a face of a box.
    
    type bndry_vector
    
        complex(dp),    allocatable,    dimension(:)    :: vec
    
    end type bndry_vector
! =================================================================================================





! =================================================================================================
! TYPE MATRIXDATA
! Contains all the relevant data for a matrix, including all relevant factorization quantities as
! as well as subroutines for inverting or applying the matrix.
!
!   Attributes:
!
!   factored            --- logical; determines if factorization data has been computed
!   mode                --- type of factorization, such as 'LU', 'SVD', 'HSS'
!   m,n                 --- size of original matrix is m x n
!   A                   --- m x n matrix; a matrix with the same size as the original matrix
!   rk                  --- integer; numerical rank of matrix (subject to some tolerance)
!   ipiv                --- pivot data; used in LU factorization
!
!
!   Type-bound Procedures:
!
!   set                 --- Takes an uninitialized matrixdata variable and sets its original matrix
!   factor              --- constructs factorization according to mode
!   matvec              --- applies the factored matrix to a vector
!   invvec              --- applies the inverse matrix to a vector (i.e. solves linear system)
!
!   Dev note: May eventually want to give some components the PRIVATE attribute
    
    type matrixdata
    
        logical                                         :: factored = .false.
        character(len=3)                                :: mode
        integer                                         :: m,n = 0
        complex(dp),    dimension(:,:), allocatable     :: A
        
        integer                                         :: rk = 0

        integer,    dimension(:),   allocatable         :: ipiv
        
        contains
        
        procedure       :: set      => mtrx_set
        procedure       :: factor   => mtrx_factor
        procedure       :: getinv   => mtrx_getinv
        procedure       :: matvec   => mtrx_matvec
        procedure       :: matmat   => mtrx_matmat
        procedure       :: invvec   => mtrx_invvec
        

    end type matrixdata
! =================================================================================================





! =================================================================================================
! TYPE ROBIN_TREE_NODE
! Derived type for every node of the Robin tree structure
!
!   node_id     --- A unique non-negative ID. Used for debugging purposes.
!   lvl         --- level = 0 is the root!
!   isleaf      --- Set to .true. if leaf
!   parent      --- pointer to parent node
!   child1,2    --- pointer to child nodes
!   robin-cond  --- Robin-to-Robin parameters
!   box         --- Format [x1 y1 z1 ; x2 y2 z2] for coordinates of box vertices
!   ptbox       --- [n_x, n_y, n_z] number of points allocated on each side of the box
!                   also referred to as 'grid_shape' in other routines
!   RtR         --- 6x6 array of arrays for the blocks of Robin-to-Robin operator. 
!                   Format [x1 x2 y1 y2 z1 z2]. Last two rows/cols unallocated for 2D.
!   iface       --- Identify which face is the interface. An x-face would be defined by x = const,
!                   so x would be the normal direction. x1 = 1, z2 = 6, etc.
!   isbndry     --- Array of logicals that says whether or not the face of the box associated with
!                   the node is part of the boundary of the full domain. 
!   D           --- The matrix representing inv(zeta - T^(2)_ii), as used in the inversion of M.*
!   S           --- The matrix representing inv(zeta - T^(1)_ii - eta^2 * D). This is the inverse
!                   of the Schur complement of inv(D) in the matrix M.*
!   nu          --- Parameter related to the Robin condition parameters, used in inversion of M.
!   g           --- Incoming boundary values. Formated as vectors supported on each of the six
!                   faces of the box. Only allocated during solution construction.
!   h           --- Outgoing boundary data. Formated as vectors supported on each of the six
!                   faces of the box. Only allocated during solution construction.
!   sol         --- Local solution to the fully inhomogeneous elliptic problem, supported on the
!                   interior of the box. Only allocated during solution construction.
!
! * The inverse of M is given by
!
!   inv(M) = [S, -nu*S*D; -nu*D*S, D + nu^2 *D*S*D] (this is written in row-major format)
!

    type robin_tree_node

        integer                         ::  node_id
        integer                         ::  lvl
        logical                         ::  isleaf
        type(robin_tree_node),  pointer ::  parent
        type(robin_tree_node),  pointer ::  child1, child2

        complex(dp),    dimension(1:2,1:2)  :: robin_cond ! <--- Still needs to be incorporated!

        real(dp),   dimension(1:2,1:3)  ::  box = 0.0d0
        integer,    dimension(1:3)      ::  ptbox = -1


        type(bndry_operator),   dimension(1:6,1:6)  :: RtR
        integer                                     :: iface
        logical,                dimension(1:6)      :: isbndry = .false.
        
        type(matrixdata)                            :: D,S
        complex(dp)                                 ::  nu  ! parameter related to inverse of M
        
        type(matrixdata)                            :: LocalOp ! Leaf-local soln operator
        
        type(bndry_vector),             dimension(1:6)  :: g
        type(bndry_vector),             dimension(1:6)  :: h
        complex(dp),    allocatable,    dimension(:)    :: sol




    end type robin_tree_node
! =================================================================================================





! =================================================================================================
! TYPE DEBUG_OPTS
! Subtype of SOLVE_OPTS, contains options related to debugging. This is by no means a full
! replacement for GDB, is used for quick debugging.
!
!   txt_output      --- Toggles verbose text output to command line
!   no_delete       --- Toggles whether to delete intermediate vector and matrix quantities in the
!                       direct solver computation

    type debug_opts
    
        logical                 :: txt_output = .false.
        logical                 :: no_delete = .false.

    end type debug_opts

! =================================================================================================





! =================================================================================================
! TYPE SOLVE_OPTS
! This derived type contains options for the local solver used.
!   h_tgt       --- Maximum spatial discretization distance.
!   h           --- The true discretization used. Differs in each direction.
!   ordering    --- Specify ordering scheme.
!   disc        --- Specify discretization type. 'fd' = finite diff
!   local_mode  --- Factorization mode for leaf-level matrices
!   R_mode      --- Factorization mode for blocks of RtR matrices
!   kb          --- k parameter for boundary condition in Robin elliptic BVP.
!   debug       --- Debug options derived type.
    type solve_opts
    
        real(dp)                    :: h_tgt
        real(dp),   dimension(1:3)  :: h
        integer                     :: ordering = 1
        character(len=2)            :: disc = 'fd'
        character(len=3)            :: local_mode = 'LU'
        character(len=3)            :: R_mode = 'LU'
        real(dp)                    :: kb = 1.0d0
        type(debug_opts)            :: debug
    
    end type solve_opts

! =================================================================================================
 




! =================================================================================================
! Factorization module interface
    interface
    
        module subroutine mtrx_set(mtrx, A)
        
            class(matrixdata),                  intent(out) :: mtrx
            complex(dp),        dimension(:,:), intent(in)  :: A
        
        end subroutine mtrx_set
              
        
        module subroutine mtrx_factor(mtrx, mode)
        
            class(matrixdata),  intent(inout)   :: mtrx
            character(len=*),   intent(in)      :: mode
        
        end subroutine mtrx_factor
        
        
        module subroutine mtrx_getinv(mtrx,InvA)
        
            class(matrixdata),                  intent(in)  :: mtrx
            complex(dp),        dimension(:,:), intent(out) :: InvA
            
        end subroutine mtrx_getinv
        
        
        module subroutine mtrx_matvec(mtrx, x, y, alpha, beta)
        
            class(matrixdata),                      intent(in)      :: mtrx
            complex(dp),        dimension(:),       intent(in)      :: x
            complex(dp),        dimension(:),       intent(inout)   :: y
            complex(dp),                            intent(in)      :: alpha, beta
        
        end subroutine mtrx_matvec
        
        
        module subroutine mtrx_matmat(mtrx, B, C, alpha, beta)
        
            class(matrixdata),                      intent(in)      :: mtrx
            complex(dp),        dimension(:,:),     intent(in)      :: B
            complex(dp),        dimension(:,:),     intent(inout)   :: C
            complex(dp),                            intent(in)      :: alpha, beta
        
        end subroutine mtrx_matmat
        
        
        module subroutine mtrx_invvec(mtrx, b, beta, x, alpha)
        
            class(matrixdata),                      intent(in)      :: mtrx
            complex(dp),        dimension(:,:),     intent(inout)   :: b
            complex(dp),                            intent(in)      :: beta
            complex(dp),        dimension(:,:),     intent(inout),  optional    :: x
            complex(dp),                            intent(in),     optional    :: alpha
        
        end subroutine mtrx_invvec
        
    
    end interface
! =================================================================================================


    contains
    
    

!
!
!
!! =================================================================================================
!! Builds
!    subroutine
!
!! =================================================================================================



end module derived_type_module
