! =================================================================================================
!                                     LOCAL SOLVER MODULE
! =================================================================================================
!
! This module contains routines for computing local solution operators, as well as local Robin
! boundary operators, on a Robin tree node. This module can be viewed as containing all node-bound
! procedures.
!
! Almost all the procedures in this module are merely interfaces to specific instances depending on
! the chosen discretization.
!
! Currently supported discretizations:
!   - 2nd order finite differences
!
! Dev note: Possibly turn these routines into type-bound procedures for robin_tree_node. Would
! simplify things with OOP.

module local_solver_module

    use derived_type_module
    use fd_module

    implicit none
    private
    
    integer,    parameter, private   :: dp = kind(1.0d0) ! double precision kind specification
    
    public  :: elliptic_invert, elliptic_apply
    public  :: leaf_solve, outgoing_bc, leaf_update, make_bndry_vec
    
    ! public :: stencil
    
! ! =================================================================================================
! ! TYPE STENCIL
! ! This derived type encodes a particular stencil. Used in FD_Method type.
! ! dev note: feature is not complete. This type can only compile using IBM Fortran compiler.

    ! type stencil(n)
    
        ! integer,    len                             :: n
        ! integer                                     :: d
        ! real(dp),   dimension(-n:n, -n:n, -n:n)     :: weights   
    
    ! end type stencil

! ! =================================================================================================






! =================================================================================================
! TYPE FD_METHOD
! This derived type specifies a particular finite difference method.
! dev note: feature forthcoming

! =================================================================================================






    
    contains





! =================================================================================================
! SUBROUTINE ELLIPTIC_INVERT
! Solves a general elliptic robin BVP on a rectangular domain. This domain is a subdomain of the
! domain of definition for the elliptic operator. Robin boundary condition expression is
!               i*kb*u + u_n = g
! Interfaces with different elliptic solvers.

    subroutine elliptic_invert(node, ell_op, opts, X)
    
        implicit none
        
        type(robin_tree_node),              intent(inout)   :: node
        type(elliptic_operator),            intent(in)      :: ell_op
        type(solve_opts),                   intent(inout)   :: opts
        complex(dp),    dimension(:,:), allocatable,    intent(out) :: X
        
        if (opts%disc == 'fd') then
            call FD_solve(node, ell_op, opts, X)
        else
            stop 'Error in ELLIPTIC_INVERT: Other discretizations not currently supported.'
        end if
    
    end subroutine elliptic_invert
! =================================================================================================




! =================================================================================================
! SUBROUTINE ELLIPTIC_APPLY
! Applies the solution matrix X produces by subroutine ELLIPTIC_INVERT for a specific set 'g' of
! boundary values. The output is a vector representing the particular solution to the elliptic PDE.
!
! This subroutine is only used as a debug tool for ELLIPTIC_INVERT.

    subroutine elliptic_apply(domain, ell_op, opts, X, g, u)
    
        implicit none
        
        real(dp),   dimension(1:2,1:3),     intent(in)  :: domain
        type(elliptic_operator),            intent(in)  :: ell_op
        type(solve_opts),                   intent(in)  :: opts
        complex(dp),    dimension(:,:),     intent(in)  :: X
        type(bndry_rhs),    dimension(1:6), intent(in)  :: g
        
        complex(dp),    dimension(:),   allocatable,    intent(out) :: u
        
        integer,    dimension(1:3)  :: grid_shape, grid_idx
        real(dp),   dimension(1:3)  :: real_pt, face_real_pt
        integer,    dimension(1:6)  :: pt_face
        complex(dp),    dimension(:),   allocatable :: bndry_vec
        integer :: bndry_size
        integer :: j
        
        if (opts%disc == 'fd') then
            ! reproduce grid shape
            grid_shape = -1
            do j = 1,ell_op%d
                grid_shape(j) = ceiling( (domain(2,j) - domain(1,j))/opts%h_tgt )
            end do
            
            ! compute # of points per face
            call get_pt_face(pt_face,grid_shape,ell_op%d)
            
            ! allocate solution vector u
            bndry_size = sum(pt_face(1:2*ell_op%d))
            allocate(u(1:(bndry_size + product(grid_shape(1:ell_op%d)))))
            
            ! check boundary size is correct
            if (bndry_size /= size(X,2)) then
                print *, bndry_size
                print *, size(X,2)
                stop 'Error: Width of solution matrix X does not match specified grid size.'
            else if (size(u) /= size(X,1)) then
                print *, size(u)
                print *, size(X,1)
                stop 'Error: Length of solution vector u does not match height of solution matrix X.'
            end if
            
            ! construct vector of boundary values from boundary function g
            allocate(bndry_vec(1:bndry_size))
            
            do j = 1, bndry_size
            
                face_real_pt = 0.0d0
            
                call lin2fgrid(j, grid_idx, grid_shape, ell_op%d, opts)
                call grid2real(real_pt, grid_idx, grid_shape, ell_op%d, opts, domain)
                
                
                ! x-face 1
                if (j <= pt_face(1)) then
                
                    face_real_pt(1:2) = [real_pt(2), real_pt(3)]
                    bndry_vec(j) = g(1)%g_face(face_real_pt)
                
                else if ((j > pt_face(1)) .and. (j <= sum(pt_face(1:2)))) then
                
                    face_real_pt(1:2) = [real_pt(2), real_pt(3)]
                    bndry_vec(j) = g(2)%g_face(face_real_pt)
                    
                else if ((j > sum(pt_face(1:2))) .and. (j <= sum(pt_face(1:3)))) then
                
                    face_real_pt(1:2) = [real_pt(1), real_pt(3)]
                    bndry_vec(j) = g(3)%g_face(face_real_pt)
                    
                else if ((j > sum(pt_face(1:3))) .and. (j <= sum(pt_face(1:4)))) then
                
                    face_real_pt(1:2) = [real_pt(1), real_pt(3)]
                    bndry_vec(j) = g(4)%g_face(face_real_pt)
                    
                else if ((j > sum(pt_face(1:4))) .and. (j <= sum(pt_face(1:5)))) then
                
                    face_real_pt(1:2) = [real_pt(1), real_pt(2)]
                    bndry_vec(j) = g(5)%g_face(face_real_pt)
                    
                else if ((j > sum(pt_face(1:5))) .and. (j <= sum(pt_face(1:6)))) then
                
                    face_real_pt(1:2) = [real_pt(1), real_pt(2)]
                    bndry_vec(j) = g(6)%g_face(face_real_pt)
                    
                end if

            
            end do
            
            ! Solution vector by matvec
            !u = 0.0d0
            call zgemv('n',size(X,1),bndry_size,cmplx(1.0,0.0,dp),X,size(X,1),bndry_vec,1, cmplx(0.0,0.0,dp) ,u,1)
            
            print *, 'Solution produced. Vector size: ', size(u)
            
            
        else
            stop 'Other discretizations not yet supported.'        
        end if
    
    end subroutine elliptic_apply
! =================================================================================================




! =================================================================================================
! SUBROUTINE LEAF_SOLVE
! Subroutine finds the specific solution vector to an elliptic system on a rectangle with
! homogeneous boundary condition and inhomogeneous interior condition. Replaces the solution vector
! on the node with the new computed vector.
    subroutine leaf_solve(node, ell_op, opts, f_rhs)
        
        implicit none
        type(robin_tree_node),              intent(inout)   :: node
        type(elliptic_operator),            intent(in)      :: ell_op
        type(solve_opts),                   intent(in)      :: opts
        procedure(scalar_proc), pointer,    intent(in)      :: f_rhs

        
        if (allocated(node%sol)) then
            stop 'Error in LEAF_SOLVE: Passed local solution vector must be unallocated.'
        end if
        
        if (opts%disc == 'fd') then
        
            call fd_leaf_solve(node, ell_op, opts, f_rhs)
            
        else
        
            stop 'Error: Other discretizations not currently supported by LEAF_SOLVE'
        
        end if

        if (opts%debug%txt_output) then
            write(*, '(A, I0, A)') 'LEAF_SOLVE at node id ', node%node_id, ' successful'
        end if
        
    
    end subroutine leaf_solve
! =================================================================================================




! =================================================================================================
! SUBROUTINE OUTGOING_BC
! Given a solution u, computes the outgoing boundary data sigma*u + tau*u_n. Called by the
! subroutine FORWARD_SUB. 
    subroutine outgoing_bc(node, sigma, tau, d, opts)
    
        implicit none
        
        type(robin_tree_node),                  intent(inout)   :: node
        complex(dp),                            intent(in)      :: sigma, tau
        integer,                                intent(in)      :: d
        type(solve_opts),                       intent(in)      :: opts
        
        if (opts%disc == 'fd') then
        
            call fd_outgoing_bc(node, sigma, tau, d, opts)
            
        else
        
            stop 'Error: Other discretizations not supported by OUTGOING_BC.'
        
        end if
        
    
    end subroutine outgoing_bc
! =================================================================================================




! =================================================================================================
! SUBROUTINE LEAF_UPDATE
! Constructs the solution to an elliptic problem on a rectangular domain with homogeneous interior
! (differential) condition and inhomogeneous Robin boundary condition. Adds to the existing
! solution vector associated to the input node. Called by BACKWARD_SUB.
    subroutine leaf_update(node, ell_op, opts)
    
        implicit none
        type(robin_tree_node),              intent(inout)   :: node
        type(elliptic_operator),            intent(in)      :: ell_op
        type(solve_opts),                   intent(in)      :: opts
        
        ! Check the local node solution is already allocated
        if (.not. allocated(node%sol)) then
        
            stop 'Error: Local solution must already be allocated in LEAF_UPDATE.'
        
        end if
        
        if (opts%disc == 'fd') then
        
            call fd_leaf_update(node, ell_op, opts)
        
        else
        
            stop 'Error: Other discretizations not supported by LEAF_UPDATE.'
        
        end if
        
        if (opts%debug%txt_output) then
            write(*, '(A, I0, A)') 'LEAF_UPDATE at node id ', node%node_id, ' successful'
        end if
    
    
    end subroutine leaf_update
! =================================================================================================





! =================================================================================================
! SUBROUTINE MAKE_RHS_VEC
! Constructs a complex-valued vector given a discretization and the RHS functions of a linear
! elliptic system defined locally on a Robin tree node.
!
!   Mode 1: f = 0, g variable
!   Mode 2: f variable, g = 0
!   Mode 3: f and g variable (fully inhomogeneous)

    subroutine make_rhs_vec(rhs, node, opts, mode, vec_out)
    
        implicit none
        
        type(elliptic_rhs),                             intent(in)  :: rhs
        type(robin_tree_node),                          intent(in)  :: node
        type(solve_opts),                               intent(in)  :: opts
        integer,                                        intent(in)  :: mode
        
        complex(dp),    dimension(:),   allocatable,    intent(out) :: vec_out
        
               
        if (opts%disc == 'fd') then
        
            call fd_rhs_vec(rhs, node, opts, mode, vec_out)
        
        else
        
            stop 'Error in MAKE_RHS_VEC: Other discretizations not currently supported'
        end if
    
    end subroutine make_rhs_vec
! =================================================================================================






! =================================================================================================
! SUBROUTINE MAKE_BNDRY_VEC
! Constructs a derived type BNDRY_VECTOR whose values are determined by function evaluations given
! by input of type BNDRY_RHS. 

    subroutine make_bndry_vec(node, g_rhs, d, opts)
    
        implicit none
        
        type(robin_tree_node),              intent(inout)   :: node
        type(bndry_rhs),    dimension(1:6), intent(in)      :: g_rhs        
        integer,                            intent(in)      :: d
        type(solve_opts),                   intent(in)      :: opts
        
        if (opts%disc == 'fd') then
        
            call fd_make_bndry_vec(node, g_rhs, d, opts)
        
        else
            
            stop 'Error in MAKE_BNDRY_VEC: Other discretizations not currently supported.'
        
        end if
        
        if (opts%debug%txt_output) then
            write(*, '(A, I0, A)') 'MAKE_BNDRY_VEC at node id ', node%node_id, ' successful'
        end if
    
    end subroutine make_bndry_vec
! =================================================================================================





! =================================================================================================
! SUBROUTINE COEFF_GRID
! This subroutine evaluates the 3x3-matrix-valued coefficient in the elliptic operator given
! the interior grid index points on a specified domain which must be a subdomain of ell_op.

        ! dev note: routine will fail if boundary point is requested (until this subroutine is
        ! generalized to allow any ordering/grid pt scheme)

    subroutine coeff_grid(coeff_out,grid_pt,domain,opts,ell_op)
    
        implicit none
        
        complex(dp),    dimension(1:3,1:3),     intent(out)     :: coeff_out
        integer,        dimension(1:3),         intent(in)      :: grid_pt
        real(dp),       dimension(1:2,1:3),     intent(in)      :: domain
        type(solve_opts),                       intent(in)      :: opts
        type(elliptic_operator),                intent(in)      :: ell_op
        
        ! the -h/2 shift is used to align grid so that the closest interior point to the boundary
        ! is at a distance h/2 from the boundary.
        coeff_out = ell_op%coeff(domain(1,:) + grid_pt * opts%h - opts%h/2)
        
        

    end subroutine coeff_grid


! =================================================================================================










! =================================================================================================
! SUBROUTINE GET_KB
! This subroutine returns a constant value 'kb' for the purpose of the Robin boundary data used in
! the multi-level robin solver.

! Dev note: this is incomplete
    subroutine get_kb(opts, ell_op, method)
    
        implicit none
        
        type(solve_opts),           intent(inout)   :: opts
        type(elliptic_operator),    intent(in)      :: ell_op
        character(len=3),           intent(in)      :: method
        
        integer,    dimension(1:3)  :: n
        integer     :: jx, jy, jz
        real(dp),   dimension(1:3)  :: eval_pt
        
        real(dp)    :: cur_val
        
        cur_val = 0.0d0
        n = ceiling( (ell_op%domain(2,:) - ell_op%domain(1,:))/opts%h )        
        
        do jx = 0, n(1)           
        do jy = 0, n(2)                
        do jz = 0, n(3)
        
            eval_pt = ell_op%domain(1,:) + [jx, jy, jz] * opts%h
            
            if (method == 'max') then
                cur_val = max(cur_val, abs(ell_op%k(eval_pt)))
                
            ! default to 'avg' calculation for kb
            else
    
            end if
                
        end do                    
        end do  
        end do
        

        
    
    end subroutine get_kb
! =================================================================================================  






end module local_solver_module