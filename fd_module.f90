! =================================================================================================
!                                   FINITE DIFFERENCES MODULE
! =================================================================================================

! This module contains specific routines for indexing and discretizing elliptic operators using
! the method of finite differences as the underlying discretization. To be called by
! local_solver_module.
!
! List of subroutines:
!   fd_matrix       --- Construct sparse matrix for elliptic system.
!   fd_matrix_row   --- Used by fd_matrix to construct matrix elements row-wise.
!   fd_solve        --- Solve the homogeneous interior, inhomogeneous elliptic BVP.
!   fd_rhs_vec
!   fd_leaf_solve
!   fd_outgoing_bc
!   grid2real
!   fgrid2lin
!   igrid2lin
!   lin2fgrid
!   lin2igrid
!   grid_classify
!   get_pt_face
!   get_h


module fd_module

    use derived_type_module
    
    implicit none
    private
    
    integer,    parameter, private   :: dp = kind(1.0d0) ! double precision kind specification
    
    public  :: fd_matrix, fd_solve, fd_rhs_vec, fd_leaf_solve, fd_outgoing_bc, fd_leaf_update, &
                fd_bndry_vec_malloc, fd_make_bndry_vec, &
                grid2real, fgrid2lin, lin2fgrid, lin2igrid, get_pt_face, get_h
    
    contains
    
! =================================================================================================
! SUBROUTINE FD_MATRIX
!
! This routine produces a dense matrix representing the elliptic differential operator on the
! interior of a rectangular domain. Because this routine is designed for dense inversion, the
! desired matrix should be very small, ideally no more than 1024 x 1024. 

       
    subroutine FD_matrix(domain,grid_shape,ell_op,opts,A)
    
        implicit none
        
        real(dp),   dimension(1:2,1:3), intent(in)      :: domain
        integer,    dimension(1:3),     intent(in)      :: grid_shape
        type(elliptic_operator),        intent(in)      :: ell_op
        type(solve_opts),               intent(in)      :: opts
        
        complex(dp),dimension(:,:), allocatable,    intent(out) :: A
        
        integer                                 :: numpt
        integer                                 :: j
        
        ! Check that A is unallocated
        if (allocated(A)) then
            stop 'Error: Array A passed to FD_MATRIX must be unallocated.'
        end if
        
        ! Determine size of domain grid to allocate matrix
        
        numpt = product(grid_shape(1:ell_op%d))
        
        do j = 1,ell_op%d
        
            numpt = numpt + 2*product(grid_shape(1:ell_op%d)) / grid_shape(j)
        
        end do
        
        allocate(A(1:numpt,1:numpt))
        
        A = 0.0d0
        
        ! Construct matrix; stencil is dimension-dependent
        ! If grid_shape = (a,b,c)
        ! Loop between (0:a+1,0:b+1,0:c+1)
        
        do j = 1, numpt
        
            call FD_matrix_row(j, grid_shape, domain, ell_op, opts, A)
        
        end do     
        
        
    
    end subroutine FD_matrix
! =================================================================================================






! =================================================================================================    
    
    subroutine FD_matrix_row(eval_pt, grid_shape, domain, ell_op, opts, A)
    
        implicit none
        
        integer,                        intent(in)      :: eval_pt
        integer,    dimension(1:3),     intent(in)      :: grid_shape
        real(dp),   dimension(1:2,1:3), intent(in)      :: domain
        type(elliptic_operator),        intent(in)      :: ell_op
        type(solve_opts),               intent(in)      :: opts
        complex(dp),dimension(:,:),     intent(inout)   :: A
        
        complex(dp)                 :: alpha, beta
        real(dp),   dimension(1:3)  :: h
        integer                     :: d
        integer                     :: n_direction
        integer,    dimension(1:3)  :: grid_pt, sten_pt, f_nbr, shift
        integer                     :: eval_pt2, dist
        real(dp),   dimension(1:3)  :: real_pt
        real(dp),   dimension(1:3)  :: sten_wt_arr
        logical                     :: io_stat
        
        integer     :: j 

        d = ell_op%d
        h = opts%h
        
        ! Set boundary condition parameters alpha, beta
        alpha = cmplx(0.0,1.0,dp) * opts%kb
        beta = cmplx(1.0,0.0,dp)

        call lin2fgrid(eval_pt, grid_pt, grid_shape, ell_op%d, opts)
        
        call grid2real(real_pt, grid_pt, grid_shape, d, opts, domain)
        
        call grid_classify(dist, f_nbr, grid_pt, grid_shape, d)
        
        !print *, eval_pt, 'coord :', grid_pt, 'f_nbr :', f_nbr
        
        ! Interior stencil
        ! Enforce differential condition Au = f
        if (dist > 0) then
        
            ! center pt
            A(eval_pt, eval_pt) = A(eval_pt,eval_pt) + ell_op%k(real_pt)
            
            ! 1D stencil in each direction
            do j = 1,d
            
                if (f_nbr(j) == 2*j - 1) then
                    sten_wt_arr = [8.0d0/3.0d0, -4.0d0, 4.0d0/3.0d0]
                else if (f_nbr(j) == 2*j) then
                    sten_wt_arr = [4.0d0/3.0d0, -4.0d0, 8.0d0/3.0d0]
                else
                    sten_wt_arr = [1.0d0, -2.0d0, 1.0d0]
                end if
           
                ! + 0
                A(eval_pt,eval_pt) = A(eval_pt, eval_pt) + (sten_wt_arr(2))/ (h(j)**2)
                
                ! + 1
                shift = [0,0,0]
                shift(j) = 1
                
                sten_pt = grid_pt + shift
                
                call fgrid2lin(eval_pt2, sten_pt, grid_shape, d, opts, io_stat)
                
                A(eval_pt, eval_pt2) = A(eval_pt, eval_pt2) + sten_wt_arr(3)/(h(j)**2)
                
                ! -1
                shift = [0,0,0]
                shift(j) = -1
                
                sten_pt = grid_pt + shift
                
                call fgrid2lin(eval_pt2, sten_pt, grid_shape, d, opts, io_stat)
                
                A(eval_pt, eval_pt2) = A(eval_pt, eval_pt2) + sten_wt_arr(1)/(h(j)**2)
            
            end do
            
        ! boundary stencil
        ! use 'ghost point' method to enforce second-order accurate Robin BC
        ! alternatively: use one-sided 3-pt normal derivative with weights (8/3, -3, 1/3) for
        ! derivative pointing to the left
        else if (dist == 0) then
        
            ! determine normal direction
            sten_pt = grid_pt
            shift = [0,0,0]
            
            do j = 1,d
            
                if (grid_pt(j) == 0) then
                    sten_pt(j) = sten_pt(j) + 1
                    n_direction = j
                    shift(j) = 1
                else if (grid_pt(j) == grid_shape(j) + 1) then
                    sten_pt(j) = sten_pt(j) - 1
                    n_direction = j
                    shift(j) = -1
                end if
            
            end do
            
            A(eval_pt, eval_pt) = A(eval_pt,eval_pt) + alpha + (8.0d0/3.0d0)/h(n_direction)
            
            call fgrid2lin(eval_pt2, grid_pt + shift, grid_shape, d, opts, io_stat)
            
            A(eval_pt, eval_pt2) = A(eval_pt, eval_pt2) - 3.0d0/h(n_direction)
            
            call fgrid2lin(eval_pt2, grid_pt + 2*shift, grid_shape, d, opts, io_stat)
            
            A(eval_pt, eval_pt2) = A(eval_pt,eval_pt2) + (1.0d0/3.0d0)/h(n_direction)
            
        
        end if
        
        if (.not. io_stat) then
            stop 'Error in FGRID2LIN.'
        end if
        
    end subroutine FD_matrix_row
! =================================================================================================



! =================================================================================================
! SUBROUTINE FD_SOLVE
! Inverts the operator for an elliptic robin BVP on the rectangle using finite differences. Output
! is an allocated matrix X whose width is the size of the grid for the boundary of the rectangle.

    subroutine FD_solve(domain, ell_op, opts, X)
    
        implicit none
        
        real(dp),   dimension(1:2,1:3),     intent(in)      :: domain
        type(elliptic_operator),            intent(in)      :: ell_op
        type(solve_opts),                   intent(inout)   :: opts
        complex(dp),    dimension(:,:), allocatable,    intent(out)     :: X
		
        integer     :: j
		integer,	dimension(1:3)		:: grid_shape
        integer     :: num_bndry
		complex(dp),dimension(:,:),	allocatable	:: A
        
        ! lapack variables        
        integer     :: m, info
        integer,    dimension(:),   allocatable :: ipiv
        

        
		
		! determine grid size, set in opts		
		call get_h(opts%h, grid_shape, opts%h_tgt, ell_op%d, domain)
		
		! construct matrix
		call FD_matrix(domain,grid_shape,ell_op,opts,A)
        
        
       
        ! ! DEBUG: WRITE A TO FILE
        ! open(unit = 1, file = "A_out_real.txt", action = "write", status = "replace")
        ! open(unit = 2, file = "A_out_imag.txt", action = "write", status = "replace")
        ! do j = 1, size(A,1)
            ! write(1,*) real(A(j,:))
            ! write(2,*) imag(A(j,:))
        ! end do
        
        ! allocate and set X        
        num_bndry = 0
        do j = 1,ell_op%d      
            num_bndry = num_bndry + 2*product(grid_shape(1:ell_op%d)) / grid_shape(j)        
        end do
		
        allocate(X(1:size(A,1), 1:num_bndry))
        
        X = cmplx(0.0,0.0,dp)
        do j = 1,num_bndry
            X(j,j) = cmplx(1.0,0.0,dp)
        end do
        
		
		! construct inverse to A; invoke LAPACK routines
        m = size(A,1)
        allocate(ipiv(1:m))
               
        !print *, 'Size of A:', size(A,1), size(A,2)
        call zgesv(m, num_bndry, A, m, ipiv, X, m, info )
        
        ! display error if inversion fails
        if (info /= 0) then
            stop 'ZGETRI: INFO is nonzero.'
        end if
       
    
    end subroutine FD_solve

! =================================================================================================





! =================================================================================================
! SUBROUTINE FD_RHS_VEC
    subroutine fd_rhs_vec(rhs, node, opts, mode, vec_out)
    
        implicit none
        
        type(elliptic_rhs),                             intent(in)  :: rhs
        type(robin_tree_node),                          intent(in)  :: node    
        type(solve_opts),                               intent(in)  :: opts
        integer,                                        intent(in)  :: mode
        complex(dp),    dimension(:),   allocatable,    intent(inout) :: vec_out
        
        real(dp),   dimension(1:3)      :: real_pt, face_real_pt
        integer,    dimension(1:3)      :: grid_idx
        integer,    dimension(1:6)      :: pt_face
        integer                         :: bndry_size
        
        integer :: j
        
        
        ! Mode 1: Produce vector supported on the boundary of the domain
        if (mode == 1) then
        
            call get_pt_face(pt_face,node%ptbox,rhs%d)
            
            bndry_size = sum(pt_face(1:2*rhs%d))
            
            ! Allocate vector
            allocate(vec_out(1:bndry_size))
            
            
            do j = 1, bndry_size
            
                face_real_pt = 0.0d0
            
                call lin2fgrid(j, grid_idx, node%ptbox, rhs%d, opts)
                call grid2real(real_pt, grid_idx, node%ptbox, rhs%d, opts, node%box)
                
                
                ! x-face 1
                if (j <= pt_face(1)) then
                
                    face_real_pt(1:2) = [real_pt(2), real_pt(3)]
                    vec_out(j) = rhs%g(1)%g_face(face_real_pt)
                
                ! x-face 2
                else if ((j > pt_face(1)) .and. (j <= sum(pt_face(1:2)))) then
                
                    face_real_pt(1:2) = [real_pt(2), real_pt(3)]
                    vec_out(j) = rhs%g(2)%g_face(face_real_pt)
                ! y-face 1
                else if ((j > sum(pt_face(1:2))) .and. (j <= sum(pt_face(1:3)))) then
                
                    face_real_pt(1:2) = [real_pt(1), real_pt(3)]
                    vec_out(j) = rhs%g(3)%g_face(face_real_pt)
                
                ! y-face2
                else if ((j > sum(pt_face(1:3))) .and. (j <= sum(pt_face(1:4)))) then
                
                    face_real_pt(1:2) = [real_pt(1), real_pt(3)]
                    vec_out(j) = rhs%g(4)%g_face(face_real_pt)
                
                ! z-face 1
                else if ((j > sum(pt_face(1:4))) .and. (j <= sum(pt_face(1:5)))) then
                
                    face_real_pt(1:2) = [real_pt(1), real_pt(2)]
                    vec_out(j) = rhs%g(5)%g_face(face_real_pt)
                
                ! z-face 2
                else if ((j > sum(pt_face(1:5))) .and. (j <= sum(pt_face(1:6)))) then
                
                    face_real_pt(1:2) = [real_pt(1), real_pt(2)]
                    vec_out(j) = rhs%g(6)%g_face(face_real_pt)
                    
                end if
                
                
            
            end do
            
            
        
        else if (mode == 2) then
        
            ! Dev note: to be completed later
        
        else if (mode == 3) then
        
            stop 'Error: Mode 3 is not currently supported.'
        
        
        end if
        
    end subroutine fd_rhs_vec
! =================================================================================================




! =================================================================================================
! SUBROUTINE FD_LEAF_SOLVE
! Solves a specific elliptic system on a rectangular domain with homogeneous boundary condition and
! inhomogeneous interior (differential) condition.

    subroutine fd_leaf_solve(node, ell_op, opts, f_rhs)
    
        implicit none
        
        type(robin_tree_node),              intent(inout)   :: node    
        type(elliptic_operator),            intent(in)      :: ell_op
        type(solve_opts),                   intent(in)      :: opts
        procedure(scalar_proc), pointer,    intent(in)      :: f_rhs
        
        complex(dp),    dimension(:,:), allocatable     :: A, b
        integer                                         :: bndry_size
        integer                                         :: m, info
        integer,        dimension(:),   allocatable     :: ipiv
        integer,        dimension(1:3)                  :: grid_pt
        integer,        dimension(1:6)                  :: pt_face
        real(dp),       dimension(1:3)                  :: real_pt
        
        integer :: j

		
		! build matrix
		call FD_matrix(node%box,node%ptbox,ell_op,opts,A)
        		
        m = size(A,1)
        
        allocate(b(1:m,1:1))        
        
        ! compute # of points per face
        call get_pt_face(pt_face, node%ptbox, ell_op%d)

        
        ! allocate solution vector u = node%sol
        bndry_size = sum(pt_face(1:2*ell_op%d))
        allocate(node%sol(1:(bndry_size + product(node%ptbox(1:ell_op%d)))))
        
        ! check size calculation is compatible (unnecessary except as debug)
        if (size(node%sol,1) /= m) then
            stop 'Error: Size of node%sol does not match size of B in LEAF_SOLVE.'
        end if
        
        ! Set RHS vector
        b(1:bndry_size,1) = cmplx(0.0,0.0,dp)
        do j = 1, (m - bndry_size)
            call lin2igrid(j, grid_pt, node%ptbox, ell_op%d, opts)
            call grid2real(real_pt, grid_pt, node%ptbox, ell_op%d, opts, node%box)
            b(bndry_size+j,1) = f_rhs(real_pt)
        end do
        
        ! Linear solve by LAPACK
        allocate(ipiv(1:m))
        call zgesv(m,1,A,m,ipiv,b,m,info)
        node%sol(:) = b(:,1)
        
        if (info /= 0) then
            print *, 'Error in ZGESV. INFO:', info
            stop
        end if
    
    end subroutine fd_leaf_solve
    
! =================================================================================================





! =================================================================================================
! SUBROUTINE FD_OUTGOING_BC
! Given a solution u, computes the outgoing boundary data h = sigma*u + tau*u_n.
! Called by OUTGOING_BC.
    subroutine fd_outgoing_bc(node, sigma, tau, d, opts)
    
        implicit none
        
        type(robin_tree_node),                  intent(inout)   :: node    
        complex(dp),                            intent(in)      :: sigma, tau
        integer,                                intent(in)      :: d
        type(solve_opts),                       intent(in)      :: opts
        
        integer,    dimension(1:6)              :: pt_face
        integer,    dimension(1:3)              :: shift, eval_pt
        integer                                 :: lin_idx
        logical                                 :: io_stat

        integer :: face, j
        
        call get_pt_face(pt_face, node%ptbox, d)
        
        
        ! Compute outgoing data
        do face = 1, 2*d
        
            do j = 1, pt_face(face)
                        
                shift = [0,0,0]
                
                shift(1 + (face-1)/2) = -1 + 2*modulo(face,2)
                
                ! use the 3-pt one-sided normal derivative with weights (8/3, -3, 1/3) * h^(-1)
                
                lin_idx = sum(pt_face(1:face-1)) + j
                
                
                call lin2fgrid(lin_idx, eval_pt, node%ptbox, d, opts)
                
                
                ! + 0 term
                node%h(face)%vec(j) = sigma * node%sol(lin_idx) + &
                    tau * (8.0d0/3.0d0)/opts%h(1 + (face-1)/2) * node%sol(lin_idx)
                
                
                ! -1 term
                
                call fgrid2lin(lin_idx, eval_pt + shift, node%ptbox, d, opts, io_stat)
                
                node%h(face)%vec(j) = node%h(face)%vec(j) + &
                    tau * (-3.0d0)/opts%h(1 + (face-1)/2) * node%sol(lin_idx)
                    
                ! -2 term
                
                call fgrid2lin(lin_idx, eval_pt + 2*shift, node%ptbox, d, opts, io_stat)
                
                
                node%h(face)%vec(j) = node%h(face)%vec(j) + &
                    tau * (1.0d0/3.0d0)/opts%h(1 + (face-1)/2) * node%sol(lin_idx)
                    

                
            
            end do
        
        end do       
        
    
    end subroutine fd_outgoing_bc
! =================================================================================================





! =================================================================================================
! SUBROUTINE FD_LEAF_UPDATE
! Constructs the solution to an elliptic problem on a rectangular domain with homogeneous interior
! (differential) condition and inhomogeneous Robin boundary condition. Adds to the existing
! solution vector associated to the input node. Called by LEAF_UPDATE.
    subroutine fd_leaf_update(node, ell_op, opts)
    
        implicit none
        type(robin_tree_node),              intent(inout)   :: node    
        type(elliptic_operator),            intent(in)      :: ell_op
        type(solve_opts),                   intent(in)      :: opts
        
        complex(dp),    dimension(:,:),     allocatable :: A, b  
        integer                                         :: m, bndry_size, info
        integer,        dimension(:),       allocatable :: ipiv    
        integer,        dimension(1:6)                  :: pt_face

        integer :: face
        
        ! build matrix
		call FD_matrix(node%box,node%ptbox,ell_op,opts,A)
        		
        m = size(A,1)
        
        allocate(b(1:m,1:1))        
        
        ! compute # of points per face
        call get_pt_face(pt_face, node%ptbox, ell_op%d)
        bndry_size = sum(pt_face(1:2*ell_op%d))
        
        ! check size calculation is compatible (unnecessary except as debug)
        if (size(node%sol,1) /= m) then
            stop 'Error: Size of node%sol does not match size of B in LEAF_UPDATE.'
        end if
        
        ! Set RHS vector
        
        do face = 1, 2*ell_op%d
        
            b(1+sum(pt_face(1:face-1)):sum(pt_face(1:face)),1) = node%g(face)%vec
            
        end do
        
        b(bndry_size+1:,1) = cmplx(0.0,0.0,dp)
        
        ! Perform linear solve
        allocate(ipiv(1:m))
        call zgesv(m,1,A,m,ipiv,b,m,info)
        node%sol(:) = node%sol + b(:,1)
        
        if (info /= 0) then
            print *, 'Error in ZGESV. INFO:', info
            stop
        end if
    
    
    end subroutine fd_leaf_update
! =================================================================================================





! =================================================================================================
! SUBROUTINE FD_BNDRY_VEC_MALLOC
! Allocates a derived type BNDRY_VECTOR for the correct size when associated to h or g data for a 
! Robin tree node, using finite differences.
    subroutine fd_bndry_vec_malloc(bndry_vec, node, d, opts)
    
        implicit none
        type(bndry_vector),     dimension(1:6),     intent(inout)   :: bndry_vec
        type(robin_tree_node),                      intent(inout)   :: node    
        integer,                                    intent(in)      :: d
        type(solve_opts),                           intent(in)      :: opts
    
        integer,    dimension(1:6)  :: pt_face
        integer :: j

        
        call get_pt_face(pt_face, node%ptbox, d)
        
        !print *, pt_face
        
        do j = 1,2*d
            if (.not. allocated(bndry_vec(j)%vec)) then
                !print *, 'Allocating at face :', j
                !print *, 'Vector size :', pt_face(j)
                allocate(bndry_vec(j)%vec(1:pt_face(j)))
            end if
        end do
        
    end subroutine fd_bndry_vec_malloc
! =================================================================================================





! =================================================================================================    
! SUBROUTINE FD_MAKE_BNDRY_VEC
! Constructs a derived type BNDRY_VECTOR whose values are determined by function evaluations given
! by input of type BNDRY_RHS. Called by MAKE_BNDRY_VEC.

    subroutine fd_make_bndry_vec(node, g_rhs, d, opts)
    
        implicit none
        
        type(robin_tree_node),              intent(inout)   :: node
        type(bndry_rhs),    dimension(1:6), intent(in)      :: g_rhs        
        integer,                            intent(in)      :: d
        type(solve_opts),                   intent(in)      :: opts
        
        real(dp),   dimension(1:3)              :: real_pt, face_real_pt
        integer                                 :: lin_idx
        integer,    dimension(1:3)              :: grid_idx
        integer,    dimension(1:6)              :: pt_face
        
        integer :: face, j
        
        call get_pt_face(pt_face,node%ptbox,d)

        
        
        ! Assign vector values by looping over boundary faces
        do face = 1, 2*d    

            if (node%isbndry(face)) then
        
                do j = 1, pt_face(face)

                    face_real_pt = 0.0d0
                    
                    lin_idx = sum(pt_face(1:face-1)) + j
                    
                    call lin2fgrid(lin_idx, grid_idx, node%ptbox, d, opts)
                    call grid2real(real_pt, grid_idx, node%ptbox, d, opts, node%box)
                    
                    ! Project coordinates onto boundary local coordinates
                
                    
                    face_real_pt(1:(face-1)/2) = real_pt(1:(face-1)/2)
                    face_real_pt((face+1)/2:d-1) = real_pt((face+3)/2:)

                    node%g(face)%vec(j) = g_rhs(face)%g_face(face_real_pt)
            
                end do
                
            end if
        
        end do
    
    end subroutine fd_make_bndry_vec
! =================================================================================================





! =================================================================================================
! SUBROUTINE GRID2REAL
! Converts a grid point on the rectangular mesh into its real-valued coordinates.
    subroutine grid2real(real_pt, grid_pt, grid_shape, d, opts, domain)

        implicit none
        
        real(dp),   dimension(1:3),     intent(out) :: real_pt
        integer,    dimension(1:3),     intent(in)  :: grid_pt, grid_shape
        integer,                        intent(in)  :: d
        type(solve_opts),               intent(in)  :: opts
        real(dp),   dimension(1:2,1:3), intent(in)  :: domain
        
        integer     :: j
        
        real_pt = 0.0d0
        
        ! the -h/2 shift is used to align grid so that the closest interior point to the boundary
        ! is at a distance h/2 from the boundary.
        
        real_pt = domain(1,:) + grid_pt * opts%h - opts%h/2
            
        ! correct coordinates for the boundary case      
        do j = 1,d
        
            if (grid_pt(j) == 0) then            
                real_pt(j) = domain(1,j)                
            else if (grid_pt(j) == grid_shape(j) + 1) then            
                real_pt(j) = domain(2,j)            
            end if
            
        end do


    end subroutine grid2real
! =================================================================================================




! =================================================================================================
! SUBROUTINE FGRID2LIN
! This subroutine computes the linear index based on the integer d-dimensional grid index for the
! combined boundary and interior of a rectangular domain.
!
! If the grid shape is (a,b,c), then the first 2ab + 2bc + 2ac index values correspond to
! the boundary points. The remaining abc index values correspond to interior points.
! 
! Corner/edge points are not allowed in this indexing scheme. If one is requested, io_stat will 
! output the value .FALSE.

    subroutine fgrid2lin(lin_idx, grid_idx, grid_shape, d, opts, io_stat)
    
        implicit none
        
        integer,                        intent(out) :: lin_idx
        integer,    dimension(1:3),     intent(in)  :: grid_idx
        integer,    dimension(1:3),     intent(in)  :: grid_shape
        integer,                        intent(in)  :: d
        type(solve_opts),               intent(in)  :: opts
        logical,                        intent(out) :: io_stat
        
        integer                     :: num_bndry, b_idx
        integer,    dimension(1:6)  :: num_face
        integer                     :: edge_count, bndry_face, bndry_coord, dist
        integer,    dimension(1:3)  :: f_neighbor
        integer                     :: j
        
        io_stat = .true.
        
        call get_pt_face(num_face,grid_shape,d)
        num_bndry = sum(num_face(1:2*d))
                     

        ! Check if grid point is a corner
        ! Faces are ordered by [x1 x2 y1 y2 z1 z2], where x1 x2 correspond to the faces given by
        ! the level sets x = x_0, x = x_1, respectively.
        
        
        
        call grid_classify(dist, f_neighbor, grid_idx, grid_shape, d)
        
        edge_count = 0 
        
        if (dist == 0) then
            ! counts number of coordinates that are on an edge (>1 means corner)
            
            do j = 1,d
                if (grid_idx(j) == 0) then
                    bndry_face = 2*j-1
                    edge_count = edge_count + 1
                else if (grid_idx(j) == grid_shape(j)+1) then
                    bndry_face = 2*j
                    edge_count = edge_count + 1
                end if
            end do
            
            if (edge_count > 1) then
                print *, 'FGRID2LIN: Corner point requested.'
                io_stat = .false.
                return
            end if
        end if
        
        ! Compute linear index; divide into boundary and interior cases
                
        ! Boundary case (edge_count == 1)
        if (edge_count > 0) then
        
            ! Compute which is the orth. axis 
            bndry_coord = (bndry_face+1)/2
            
            ! Assign local face index (b_idx)
            
            call igrid2lin(b_idx, [grid_idx(:bndry_coord-1), grid_idx(bndry_coord+1:),0], &
                [grid_shape(:bndry_coord-1), grid_shape(bndry_coord+1:),0], d-1, opts)
                
            ! Compute global linear index
            
            lin_idx = sum(num_face(1:bndry_face-1)) + b_idx
                
                
        ! Interior case         
        else if (edge_count == 0) then
        
            call igrid2lin(lin_idx, grid_idx, grid_shape, d, opts)
            
            lin_idx = lin_idx + sum(num_face(1:2*d))           
        
        end if
        
        
    
    
    end subroutine fgrid2lin


! =================================================================================================





! =================================================================================================
! SUBROUTINE IGRID2LIN
! This subroutine computes the linear index based on an integer d-dimensional grid index for the
! interior of the rectangular domain. Several indexing schemes are available.
    subroutine igrid2lin(lin_idx, grid_idx, grid_shape, d, opts)
    
        implicit none
        
        integer,                        intent(out) :: lin_idx
        integer,    dimension(1:3),     intent(in)  :: grid_idx
        integer,    dimension(1:3),     intent(in)  :: grid_shape
        integer,                        intent(in)  :: d
        type(solve_opts),               intent(in)  :: opts
        
        integer             :: num_int
        integer             :: j
        
        lin_idx = 0
        
        ! Compute number of pts on the interior
        num_int = product(grid_shape(1:d))       
        
        
        ! Naive wrapping indexing scheme, first index incremented first
        ! E.g. (1,1,1) (2,1,1) (3,1,1) ...
        
        if (opts%ordering == 1) then
        
            do j = d,2,-1           
                lin_idx = lin_idx + product(grid_shape(1:j-1)) * (grid_idx(j) - 1)            
            end do
            
                lin_idx = lin_idx + grid_idx(1)
        
        else if (opts%ordering == 2) then
        
            print *, 'ordering=2 parameter chosen; stop' ! Dev note: later, alternative indexing scheme
            stop
        
        else
        
            stop 'Invalid ordering parameter in igrid2lin'
        
        end if
    
    end subroutine igrid2lin


! =================================================================================================




! =================================================================================================
! SUBROUTINE LIN2FGRID
! This subroutine computes the integer d-dimensional grid index given the linear index, which
! takes into account a grid that includes a cornerless boundary.

    subroutine lin2fgrid(lin_idx, grid_idx, grid_shape, d, opts)
    
        implicit none

        integer,                        intent(in)  :: lin_idx
        integer,    dimension(1:3),     intent(out) :: grid_idx
        integer,    dimension(1:3),     intent(in)  :: grid_shape
        integer,                        intent(in)  :: d
        type(solve_opts),               intent(in)  :: opts
        
        integer :: j
        integer,    dimension(1:6)  :: pt_face
        integer :: bndry_size
        integer,    dimension(1:3)  :: bndry_shape
        
        ! Error check
        if (lin_idx <= 0) then
            print *, 'LIN_IDX value: ', lin_idx
            stop 'Error in LIN2FGRID: LIN_IDX is not a positive integer.'
        end if

        ! determine size of boundary
        call get_pt_face(pt_face,grid_shape,d)
        
        bndry_size = sum(pt_face(1:2*d))
        
        ! linear index for boundary points
        if (lin_idx <= bndry_size) then
            ! face 1
            if (lin_idx <= pt_face(1)) then
                bndry_shape = [grid_shape(2), grid_shape(3),-1]
                call lin2igrid(lin_idx, grid_idx, bndry_shape, d-1, opts)
                grid_idx = [0, grid_idx(1), grid_idx(2)]
            
            ! face 2
            else if ((lin_idx > pt_face(1)) .and. (lin_idx <= sum(pt_face(1:2)))) then
                bndry_shape = [grid_shape(2), grid_shape(3),-1]
                call lin2igrid(lin_idx - pt_face(1), grid_idx, bndry_shape, d-1, opts)
                grid_idx = [grid_shape(1)+1,grid_idx(1), grid_idx(2)]
            
            ! face 3
            else if ((lin_idx > sum(pt_face(1:2))) .and.   (lin_idx <= sum(pt_face(1:3)))) then
                bndry_shape = [grid_shape(1), grid_shape(3),-1]
                call lin2igrid(lin_idx - sum(pt_face(1:2)), grid_idx, bndry_shape, d-1, opts)
                grid_idx = [grid_idx(1), 0, grid_idx(2)]
            
            ! face 4
            else if ((lin_idx > sum(pt_face(1:3))) .and.   (lin_idx <= sum(pt_face(1:4)))) then
                bndry_shape = [grid_shape(1), grid_shape(3),-1]
                call lin2igrid(lin_idx - sum(pt_face(1:3)), grid_idx, bndry_shape, d-1, opts)
                grid_idx = [grid_idx(1), grid_shape(2)+1, grid_idx(2)]
            
            ! face 5
            else if ((lin_idx > sum(pt_face(1:4))) .and.   (lin_idx <= sum(pt_face(1:5)))) then
                bndry_shape = [grid_shape(1), grid_shape(2),-1]
                call lin2igrid(lin_idx - sum(pt_face(1:4)), grid_idx, bndry_shape, d-1, opts)
                grid_idx = [grid_idx(1), grid_idx(2),0]
            
            ! face 6
            else if (lin_idx > sum(pt_face(1:5))) then
                bndry_shape = [grid_shape(1), grid_shape(2),-1]
                call lin2igrid(lin_idx - sum(pt_face(1:5)), grid_idx, bndry_shape, d-1, opts)
                grid_idx = [grid_idx(1), grid_idx(2),grid_shape(3)+1]
            
            end if
        
        ! linear index for interior points
        else if (lin_idx > bndry_size) then
            call lin2igrid(lin_idx - bndry_size, grid_idx, grid_shape, d, opts)
        end if
        
        
    
    end subroutine lin2fgrid
! =================================================================================================




! =================================================================================================
! SUBROUTINE LIN2IGRID
! This subroutine computes the integer d-dimensional grid index given the linear index. Several
! indexing schemes are available.

    subroutine lin2igrid(lin_idx, grid_idx, grid_shape, d, opts)
    
        implicit none
        
        integer,                        intent(in)  :: lin_idx
        integer,    dimension(1:3),     intent(out) :: grid_idx
        integer,    dimension(1:3),     intent(in)  :: grid_shape
        integer,                        intent(in)  :: d
        type(solve_opts),               intent(in)  :: opts
        
        integer :: lin
        integer :: j
        
        if (lin_idx <= 0) then
            stop 'Error: LIN_IDX must be a positive integer.'
        end if
        
        lin = lin_idx
        grid_idx = 1
        
        do j = d, 1, -1
            grid_idx(j) = 1 + (lin-1) / product(grid_shape(1:(j-1)))
            lin = lin - (grid_idx(j) - 1) * product(grid_shape(1:(j-1)))
        end do
    
    end subroutine lin2igrid


! =================================================================================================





! =================================================================================================
! SUBROUTINE GRID_CLASSIFY
! This subroutine determines various properties of a grid point, in particular distance to the
! boundary and which face it is closest to. Used to determine what stencil to use about that point.
! Note: subroutine will give nonsensical values if grid_pt and grid_shape are incompatible.
    subroutine grid_classify(dist, f_neighbor, grid_pt, grid_shape, d)
    
        implicit none
        
        integer,                        intent(out) :: dist
        integer,    dimension(1:3),     intent(out) :: f_neighbor
        integer,    dimension(1:3),     intent(in)  :: grid_pt, grid_shape
        integer,                        intent(in)  :: d
        
        integer :: j
        
        f_neighbor = 0
        
        dist = minval( [grid_pt(1:d),grid_shape(1:d) - grid_pt(1:d) + 1] )
        
        ! determine face neighbors (up to d of them) if applicable
        
        if (dist <= 1) then
        
            do j = 1,d
            
                if ( (grid_pt(j) == 0) .or. (grid_pt(j) == 1) ) then
                    f_neighbor(j) = 2*j - 1
                else if ( (grid_pt(j) == grid_shape(j) ) .or. &
                    (grid_pt(j) == grid_shape(j) +1 ) ) then
                    f_neighbor(j) = 2*j
                end if
            
            end do
        
        end if
        
        
        
        ! if (dist == 1) then
        
            ! do j = 1,d            
                ! if (grid_pt(j) == 1) then                
                    ! f_neighbor(j) = 2*j - 1              
                ! else if (grid_pt(j) == grid_shape(j)) then                
                    ! f_neighbor(j) = 2*j                
                ! end if            
            ! end do
            
        ! else if (dist == 0) then
        
            ! do j = 1,d
                ! if (grid_pt(j) == 0) then
                    ! f_neighbor(j) = 2*j - 1
                ! else if (grid_pt(j) == grid_shape(j) + 1) then
                    ! f_neighbor(j) = 2*j
                ! end if
            ! end do
        
        ! end if
    
    end subroutine grid_classify


! =================================================================================================




! =================================================================================================
! SUBROUTINE GET_PT_FACE
! Determine the number of points per face of a box, given the number of points per side of a box.
    subroutine get_pt_face(pt_face,pt_box,d)

        implicit none
        
        integer,    dimension(1:6), intent(out) :: pt_face
        integer,    dimension(1:3), intent(in)  :: pt_box
        integer,                    intent(in)  :: d
        
        integer :: j
        
        pt_face = 0
        do j = 1,d
            pt_face(2*j - 1) = product(pt_box(1:d)) / pt_box(j)
            pt_face(2*j) = pt_face(2*j-1)
        end do

    end subroutine get_pt_face
! =================================================================================================




! =================================================================================================
! SUBROUTINE GET_H
! This subroutine computes the spatial discretization step size given a target h value and a
! rectangular domain, such that the actual step size is an integer multiple of the domain width in
! each direction.
    subroutine get_h(h_out, grid_shape, h_target, d, domain)
    
        implicit none
        
        real(dp),   dimension(1:3),     intent(out)     :: h_out
        integer,    dimension(1:3),     intent(out)     :: grid_shape
        real(dp),                       intent(in)      :: h_target
        integer,                        intent(in)      :: d
        real(dp),   dimension(1:2,1:3), intent(in)      :: domain
        
        integer :: j
        
        h_out = 0.0d0
        grid_shape = -1
        
        do j = 1,d
        
            ! grid_shape value is at least 3 to avoid degenerate grids
            grid_shape(j) = max(ceiling( (domain(2,j) - domain(1,j))/h_target ), 3)            
            
            h_out(j) = ( domain(2,j) - domain(1,j) ) / grid_shape(j)
        
        end do
    
    end subroutine get_h
! =================================================================================================

end module fd_module