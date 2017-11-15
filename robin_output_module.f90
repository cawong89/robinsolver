! =================================================================================================
!                                       ROBIN OUTPUT MODULE
! =================================================================================================
!
! This module contains procedures for reading the solution output using the Robin tree data
! structure.

module robin_output_module

    use derived_type_module
    use robin_tree_module
    use fd_module
    
    implicit none
    private
    
    integer,    parameter   :: dp = kind(1.0d0) ! double precision kind specification
    
    public  :: leaf_solution
    public  :: get_solution, check_solution, write_variable





! =================================================================================================
! TYPE LEAF_SOLUTION
! A solution vector living on a leaf of a Robin tree
!
! Dev note: In the future, possible support for solution vectors living on non-rectangular domains

    type leaf_solution
    
        real(dp),       dimension(1:2,1:3)                  :: box
        complex(dp),    dimension(:),       allocatable     :: u
    
    end type leaf_solution
    
! =================================================================================================




    contains
    
! =================================================================================================

    subroutine get_solution(T,sol)
    
        implicit none
        type(robin_tree),                                       intent(in)      :: T
        type(leaf_solution),    dimension(:),   allocatable,    intent(inout)   :: sol
        
        integer :: num_leaves, iter
        
        if (.not. T%factored) then
            stop 'Error in GET_SOLUTION: Robin Tree not factored.'
        end if
        
        ! Determine the number of leaves and allocate
        num_leaves = 0
        call count_leaves(T%root, num_leaves)
        
        allocate(sol(1:num_leaves))
        
        ! Build solution 
        iter = 1
        call get_node_sol(T%root, sol, iter)
        
        ! Output
        write(*, '(A, I0, A)') 'Solution vectors extracted on ', num_leaves, ' leaf nodes.'
        
        
    
    end subroutine get_solution

! =================================================================================================





! =================================================================================================
! SUBROUTINE CHECK_SOLUTION
! Compares the solution vectors living on a LEAF_SOLUTION array to a supplied complex-valued
! function.
    subroutine check_solution(d, leaf_sol, true_sol, opts, error_vec, norm, rel)
    
        implicit none
        integer,                                    intent(in)  :: d
        type(leaf_solution),    dimension(:),       intent(in)  :: leaf_sol
        procedure(sol_proc)                                     :: true_sol
        type(solve_opts),                           intent(in)  :: opts
        character(len=3),                           intent(in)  :: norm
        logical,                                    intent(in)  :: rel
        
        real(dp),   dimension(:),   allocatable,    intent(out) :: error_vec
        
        complex(dp),    dimension(:),   allocatable :: u_exact
        
        integer :: j
        
        ! Allocate error_vec
        allocate(error_vec(1:size(leaf_sol)))
        
        ! iterate through leaves
        do j = 1, size(leaf_sol)
            call proc2vec(d, leaf_sol(j)%box, true_sol, u_exact, opts)
            
            ! Compute norm of error on leaf node
            
            if (.not. rel) then
            
                if (norm == '1') then
                    error_vec(j) = sum(abs(leaf_sol(j)%u - u_exact))
                else if (norm == '2') then
                    error_vec(j) = norm2(abs(leaf_sol(j)%u - u_exact))
                else
                    error_vec(j) = maxval(abs(leaf_sol(j)%u - u_exact))
                end if
            
            else if (rel) then
            
                if (norm == '1') then
                    error_vec(j) = sum(abs((leaf_sol(j)%u - u_exact)/u_exact))
                else if (norm == '2') then
                    error_vec(j) = norm2(abs((leaf_sol(j)%u - u_exact)/u_exact))
                else
                    error_vec(j) = maxval(abs((leaf_sol(j)%u - u_exact)/u_exact))
                end if
            
            end if
        
        end do
        

    
    end subroutine check_solution
! =================================================================================================




! =================================================================================================
! SUBROUTINE WRITE_VARIABLE
! Writes one of the variables used in the elliptic solver calculation to file.
! Opens an output file with arguments:
!   unit = unitval, file = filename, action = "write", status = "replace"
! Possibly variables to be written are:
!   'g'     --- Incoming boundary condition.
!   'h'     --- Outgoing boundary data. May not exist if debugging is not on.
!   'R'     --- Robin-to-Robin operator.
!   'u'     --- Local elliptic solution. Only exists on leaf nodes.
! Optional index arguments 'idx1' and 'idx2' 


    subroutine write_variable(T, unitval, filename, nodeid, variable, idx1, idx2)
    
        implicit none
        type(robin_tree),       intent(in)  :: T
        integer,                intent(in)  :: unitval
        character(len=*),       intent(in)  :: filename
        integer,                intent(in)  :: nodeid
        character(len=3),       intent(in)  :: variable
        integer,                intent(in)  :: idx1, idx2
        
        type(robin_tree_node),  pointer :: node
        
        ! Target node
        if ((nodeid < 0) .or. (nodeid > (size(T%node_list) - 1) ) ) then
        
            write(*, '(A)') 'Error in WRITE_VARIABLE: Argument nodeid out of bounds'
            return
            
        end if
        
        if (.not. associated(T%node_list(nodeid)%ptr) ) then
        
            write(*, '(A)') 'Error in WRITE_VARIABLE: Requested node pointer not associated'
            return
            
        end if
        
        node => T%node_list(nodeid)%ptr
        
        ! Open output file
        open(unit = unitval, file = filename, action = "write", status = "replace")
        
        
        
        
        
        ! Write g data
        if (variable == 'g') then
        
            if ((idx1 <= 0) .or. (idx1 > 2*T%ell_op%d)) then
            
                write(*, '(A)') 'Error in WRITE_VARIABLE: Argument idx1 out of bounds'
                
            end if
            
        
        ! Write h data
        else if (variable == 'h') then
        
        ! Write RtR data
        else if (variable == 'R') then
        
        ! Write local solution data
        else if (variable == 'u') then
        
        else
        
            write(*, '(A)') 'Error in WRITE_VARIABLE: Invalid variable name specified.'
            return
            
        end if
        
        
    
    end subroutine write_variable
! =================================================================================================





! =================================================================================================
! SUBROUTINE PROC2VEC
! Similar to MAKE_RHS_VEC (in LOCAL_SOLVER_MODULE), this subroutine constructs a vector based upon
! a scalar procedure, evaluated on a rectangular domain. The discretization used depends on 
! SOLVE_OPTS.

    subroutine proc2vec(d, box, proc, vec_out, opts)
    
        implicit none
        integer,                                                intent(in)      :: d
        real(dp),               dimension(1:2,1:3),             intent(in)      :: box
        procedure(sol_proc)                                                     :: proc
        complex(dp),            dimension(:),   allocatable,    intent(out)     :: vec_out
        type(solve_opts),                                       intent(in)      :: opts
        
        integer,    dimension(1:3)  :: grid_shape, grid_idx
        real(dp),   dimension(1:3)  :: real_pt
        integer,    dimension(1:6)  :: pt_face
        real(dp),   dimension(1:3)  :: h_out
        
        integer :: lin_idx
        
        
        if (opts%disc == 'fd') then
            
            ! Build finite-difference grid parameters
            call get_h(h_out, grid_shape, opts%h_tgt, d, box)
            
            call get_pt_face(pt_face,grid_shape,d)
            
            if (.not. allocated(vec_out)) then
                allocate(vec_out(1:(sum(pt_face(1:2*d)) + product(grid_shape(1:d)))))
            end if
            
            ! Construct vector from procedure
            do lin_idx = 1, size(vec_out)
                call lin2fgrid(lin_idx, grid_idx, grid_shape, d, opts)
                call grid2real(real_pt, grid_idx, grid_shape, d, opts, box)
                vec_out(lin_idx) = proc(real_pt)
            end do
        
        else
            stop 'Error in PROC2VEC: Other discretizations not currently supported.'
        end if
    
    end subroutine proc2vec
! =================================================================================================




! =================================================================================================
! SUBROUTINE COUNT_LEAVES

    recursive subroutine count_leaves(node, num_leaves)
    
        implicit none
        type(robin_tree_node),      intent(in)      :: node
        integer,                    intent(inout)   :: num_leaves
        
        if (node%isleaf) then
            num_leaves = num_leaves + 1
        else
            call count_leaves(node%child1, num_leaves)
            call count_leaves(node%child2, num_leaves)
        end if
    
    end subroutine count_leaves

! =================================================================================================




! =================================================================================================
! SUBROUTINE GET_NODE_SOL
    
    recursive subroutine get_node_sol(node, sol, iter)
    
        implicit none
        type(robin_tree_node),                  intent(in)      :: node
        type(leaf_solution),    dimension(:),   intent(inout)   :: sol
        integer,                                intent(inout)   :: iter
        
        if (node%isleaf) then
            allocate(sol(iter)%u(1:size(node%sol)))
            sol(iter)%box = node%box
            sol(iter)%u = node%sol
            iter = iter + 1
        else
            call get_node_sol(node%child1, sol, iter)
            call get_node_sol(node%child2, sol, iter)
        end if
    
    end subroutine
! =================================================================================================


    


end module robin_output_module