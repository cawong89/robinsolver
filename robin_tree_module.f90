! =================================================================================================
!                                       ROBIN TREE MODULE
! =================================================================================================
!
!
! This module defines the tree structure and its core procedures used in the construction of a
! factorization and direct solution for a Helmholtz operator using Robin-to-Robin operators.
!
! References:
!
! -- Gillman, Barnett, Martinsson. "A spectrally accurate direct solution technique for frequency-
! domain scattering problems with variable media".
!
! -- Liu, de Hoop, Xia. "Interconnected hierarchical rank-structured method for direct elliptic
! solution and damped high-frequency Helmholtz preconditioning".
!
!
!
! Author: Christopher A. Wong
! Date: 26 October 2017
! Standard: Fortran 2008

module robin_tree_module

    use derived_type_module
    use fd_module
    use local_solver_module
    use time_module

    implicit none
    private

    integer,    parameter   :: dp = kind(1.0d0) ! double precision kind specification

    public                  :: robin_tree





! =================================================================================================
! TYPE NODE_PTR
    type node_ptr
        type(robin_tree_node),  pointer :: ptr
    end type node_ptr
! =================================================================================================





! =================================================================================================
! TYPE ROBIN_TREE
! An array containing the nodes of the Robin tree
    type robin_tree

        integer                                 :: nlvls    ! number of levels
        type(elliptic_operator)                 :: ell_op   ! associated elliptic operator
		type(solve_opts)						:: opts		! options related to linear solver
        type(robin_tree_node),  pointer         :: root     ! Root node (whole domain)
        complex(dp),    dimension(1:2,1:2)      :: bndry_params ! bndry values and data 
		logical									:: factored = .false. ! true if factoring complete
        real(dp)                                :: f_time   ! Time (in sec) for factorization
        real(dp)                                :: s_time   ! Time (in sec) for linear solve
        
        type(node_ptr), dimension(:),   allocatable :: node_list    ! array of pointers to nodes

        contains

        procedure       :: factor   => robin_tree_factor
        procedure       :: solve    => robin_tree_solve
        procedure       :: memsize  => robin_tree_memsize
        procedure       :: inspect  => robin_tree_inspect
        
        ! dev note: Per F2003 syntax, type-bound procedures must take the type as a class argument

    end type robin_tree
! =================================================================================================







    contains



! =================================================================================================
! SUBROUTINE ROBIN_TREE_FACTOR
!
! Main subroutine for constructing tree-based factorization of the elliptic BVP using
! Robin-to-Robin operators.

    subroutine robin_tree_factor(T, L, ell_op, opts)

        implicit none

        class(robin_tree),          intent(inout)               :: T
        integer,                    intent(in)                  :: L
        type(elliptic_operator),    intent(in)                  :: ell_op
        type(solve_opts),           intent(in)                  :: opts
        
        real(dp),   dimension(1:L)  :: lvl_mem  ! this array estimates bytes of memory per level
        real(dp),   parameter       :: min_det = (10.0)**(-3.0) ! min det for bndry_params matrix
        
        integer :: j
        logical :: io_stat
        character(len=6)    :: mem_char, lvl_char


        ! ====================================================
        
        call start_clock()
        
        ! Set T properties
        T%nlvls = L
        T%ell_op = ell_op
        T%opts = opts
        
        ! Check validity of ell_op
        
        if (minval(T%ell_op%domain(2,:) - T%ell_op%domain(1,:)) < 0) then
        
            stop 'Error in ROBIN_TREE_FACTOR: Invalid operator domain format.'
            
        else if ( .not. ( (T%ell_op%d == 2) .or. (T%ell_op%d == 3)  )  ) then
        
            stop 'Error in ROBIN_TREE_FACTOR: Dimension must be 2 or 3.'
        
        end if
        
        ! Check validity of bndry_params
        
        if ( abs(T%bndry_params(1,1) * T%bndry_params(2,2) - & 
            T%bndry_params(1,2) * T%bndry_params(2,1)) < min_det ) then
        
            stop 'Error: Boundary parameter matrix is nearly singular.'
        
        end if
        


        ! Generate nodes of the tree

        allocate(T%root)
        allocate(T%node_list(0:(2**(L+1) - 2)))
        
        T%root%node_id = 0
        T%root%lvl = 0
        T%root%isleaf = (0 == L)
		T%root%parent => NULL()
		
		T%root%box = T%ell_op%domain
        
        T%root%isbndry(1:2*T%ell_op%d) = .true.

		call bin_tree(T%root,L,T%ell_op%d)
        call node_id_ptr(T,T%root)
        
        write(*, '(A)') 'Domain partition completed.'
		
		! Begin RTR computation
        
        io_stat = .true.
        
        call compute_RTR(T, T%root)
        
        write(*, '(A)') 'RTR computation completed.'
        
        ! Dev note: variable io_stat tracking routine error currently not passed to compute_RTR
	

        if (io_stat) then
            T%factored = .true.
            write(*, '(A)') 'Factorization completed.'
        else
            T%factored = .false.
            write(*, '(A)') 'Factorization was not completed. Something went wrong.'
        end if
        
        call read_clock(T%f_time)


    end subroutine robin_tree_factor
! =================================================================================================




! =================================================================================================
! SUBROUTINE ROBIN_TREE_SOLVE
!
! Main subroutine for solving the boundary value elliptic problem using the factorization data
! contained with the Robin tree data structure. Can only be used once the Robin tree data structure
! has been fully populated via subroutine robin_tree_factor.

    subroutine robin_tree_solve(T,rhs)

        implicit none

        class(robin_tree),              intent(inout)      :: T
        type(elliptic_rhs),             intent(in)      :: rhs

        
        call start_clock()
        
        ! ==========================================================
        ! Check conditions
		if (.not. T%factored) then
			stop 'Error: Robin Tree must be factored first.'
        else if ( T%ell_op%d /= rhs%d ) then
            ! write(*,*) 'Error: Dimensions of Robin Tree and RHS do not match.'
            stop 'Error: Dimensions of Robin Tree and RHS do not match.'
        end if
		
		
		! Forward and backward tree reversal for solution construction
        
        write(*, '(A)') 'Commencing forward substitution...'
        
        call forward_sub(T, T%root, rhs)
        
        write(*, '(A)') 'Commencing backward substitution...'
        
        call backward_sub(T, T%root, rhs)

        write(*, '(A)') 'Elliptic solve completed.'
        
        call read_clock(T%s_time)
        


    end subroutine robin_tree_solve

! =================================================================================================




! =================================================================================================
! SUBROUTINE ROBIN_TREE_MEMSIZE
!   Checks the amount of memory required in the storage of the fully populated Robin operator tree.
    subroutine robin_tree_memsize(T)
    
        implicit none
        class(robin_tree),  intent(in)  :: T
        
        ! Dev note: (obviously) not complete
    
    end subroutine robin_tree_memsize

! =================================================================================================





! =================================================================================================
! SUBROUTINE ROBIN_TREE_INSPECT
!   Inspects a Robin Tree and outputs text reporting its status. Mainly used for debugging purposes.
!
!   Dev note: work in progress
    subroutine robin_tree_inspect(T)
    
        implicit none
        class(robin_tree),  intent(in)  :: T
        
        write(*, '(A,L1)') 'Factored: ', T%factored
        
        if (T%factored) then
            write(*, '(A,I0)') 'Number of levels: ', T%nlvls
        end if
    
    end subroutine robin_tree_inspect

! =================================================================================================





! =================================================================================================
! SUBROUTINE BIN_TREE
!
! Recursive subroutine for producing the binary tree partitioning of the domain, setting the
! non-factorization related data for each node of the Robin tree, up to level 'maxlvl'. 

    recursive subroutine bin_tree(node,maxlvl,d)

        implicit none

        type(robin_tree_node),  pointer,	intent(inout)	:: node
		integer,							intent(in)		:: maxlvl
		integer,							intent(in)		:: d
        
        integer,    save    :: id_counter
        
        if (node%lvl == 0) then
            id_counter = 0
        end if
		
		
		
		! Terminate at leaf node
		if ( node%isleaf .eqv. .true. ) then
		
			return
			
			
		! otherwise construct children
		else
		
			allocate(node%child1)
			allocate(node%child2)

			call node_divide(d,node)
		
            id_counter = id_counter + 1
            node%child1%node_id = id_counter
			node%child1%lvl = node%lvl + 1
			node%child1%isleaf = (node%child1%lvl == maxlvl)
			node%child1%parent => node
			
			
            id_counter = id_counter + 1
            node%child2%node_id = id_counter
			node%child2%lvl = node%lvl + 1
			node%child2%isleaf = (node%child2%lvl == maxlvl)
			node%child2%parent => node
            
            
            call bin_tree(node%child1,maxlvl,d)			
			call bin_tree(node%child2,maxlvl,d)
			
		end if

    end subroutine bin_tree


! =================================================================================================




! =================================================================================================
! SUBROUTINE NODE_ID_PTR
! Assigns each pointer in the node list to the node with corresponding node ID.
    recursive subroutine node_id_ptr(T,node)
    
        implicit none
        type(robin_tree),               intent(inout)   :: T
        type(robin_tree_node),  target, intent(in)      :: node
        
        T%node_list(node%node_id)%ptr => node
        
        if (.not. node%isleaf) then
            
            call node_id_ptr(T,node%child1)
            call node_id_ptr(T,node%child2)
        
        end if
    
    end subroutine node_id_ptr
! =================================================================================================





! =================================================================================================
! SUBROUTINE COMPUTE_RTR
!
! Subroutine that computes the Robin-to-Robin boundary operators on a given node of the Robin tree.

    recursive subroutine compute_RTR(T, node)
    
        implicit none
        
		type(robin_tree),						intent(inout)   :: T
        type(robin_tree_node),                  intent(inout)   :: node
        
        real(dp),   dimension(1:3)      :: diffbox
        real(dp),   dimension(1:2,1:3)  :: diffbox2
        integer,    dimension(1:6)      :: pt_face
        complex(dp) :: zeta, nu
        integer :: iface_size, rsize, csize
        integer :: j, row,col
        complex(dp),    dimension(:,:), allocatable :: C11,C12,C21,C22
        
        ! lapack variables        
        integer     :: info, lwork
        integer,    dimension(:),   allocatable ::  ipiv
        complex(dp),dimension(:),   allocatable ::  work
        
        ! ! Determine which bndry face of the node is the interface
        ! if ( associated( node%parent ) ) then
            ! diffbox2 = abs(node%parent%box - node%box)
                            
            ! node%iface = f_idx(maxloc(maxval(diffbox2(:,1:T%ell_op%d),2),1), &
                ! maxloc(maxval(diffbox2(:,1:T%ell_op%d),1),1))
        ! end if

        ! Leaf computation: Compute RTR directly from local elliptic solve
        if ( node%isleaf ) then
        
            call get_h(T%opts%h, node%ptbox, T%opts%h_tgt, T%ell_op%d, node%box)  
            
            call RTR_solve(node, T%ell_op, T%opts)
            
            
        else
        
            ! recurse
            call compute_RTR(T,node%child1)
            call compute_RTR(T,node%child2)

            ! determine width of merged box for current node
            diffbox = abs(box_width(node%box) - box_width(node%child1%box) )            
            node%ptbox = node%child1%ptbox            
            node%ptbox(maxloc(diffbox(1:T%ell_op%d),1)) = &
                2 * node%ptbox(maxloc(diffbox(1:T%ell_op%d),1))

            
            call get_pt_face(pt_face,node%ptbox,T%ell_op%d)

                        
			
			! M-matrix block allocation
            
            iface_size = size(node%child1%RtR(node%child1%iface,node%child1%iface)%mat,1)
            
            if ( .not. size(node%child1%RtR(node%child1%iface,node%child1%iface)%mat) == &
                size(node%child2%RtR(node%child2%iface,node%child2%iface)%mat)  ) then
                
                stop 'Error: Interface sizes of node children do not match'
                
            end if
            
            ! Lapack block workspace allocation
            lwork = max(64, iface_size)
            allocate(work(1:lwork))
            
            
            allocate(node%D(1:iface_size, 1:iface_size), node%S(1:iface_size, 1:iface_size)   )
            allocate(ipiv(1:iface_size))
            
            nu = (-T%bndry_params(1,1) * T%bndry_params(2,2) + & 
                T%bndry_params(1,2) * T%bndry_params(2,1)) / &
                (  2 * T%bndry_params(1,1) * T%bndry_params(1,2) )
                
            zeta = (T%bndry_params(1,1) * T%bndry_params(2,2) + & 
                T%bndry_params(1,2) * T%bndry_params(2,1)) / &
                (  2 * T%bndry_params(1,1) * T%bndry_params(1,2) )
                
            node%nu = nu
            
                
            ! D matrix block inversion, D = inv(zeta * I - T^(2)_ii)
            
            node%D = - node%child2%RtR(node%child2%iface,node%child2%iface)%mat
            
            do j = 1, iface_size            
                node%D(j,j) = node%D(j,j) + zeta            
            end do

            
            call zgetrf(iface_size, iface_size, node%D, iface_size, ipiv, info)


            call zgetri(iface_size, node%D, iface_size, ipiv, work, lwork, info)
            
            
            ! S matrix block inversion, S = inv(zeta*I - T^(1)_ii - nu^2 * D)
            
            node%S = - node%child1%RtR(node%child1%iface,node%child1%iface)%mat - &
                (nu*nu) * node%D
            
            do j = 1, iface_size
                node%S(j,j) = node%S(j,j) + zeta
            end do
            
            
            call zgetrf(iface_size, iface_size, node%S, iface_size, ipiv, info)
  
            call zgetri(iface_size, node%S, iface_size, ipiv, work, lwork, info)
            
            
            ! Deallocate T^(1)_ii and T^(2)_ii as they are no longer needed
            if (.not. T%opts%debug%no_delete) then
                deallocate(node%child1%RtR(node%child1%iface,node%child1%iface)%mat)
                deallocate(node%child2%RtR(node%child2%iface,node%child2%iface)%mat)
            end if
            
            ! Construct blocks of RtR matrix. This loop is broken down into subcases corresponding
            ! to faces of the boundary that are merged, and those are in the direction of the
            ! interface.

            ! Outer loop: partially form matrix blocks (from the right)
            do col = 1,2*T%ell_op%d            
                
                
                ! col case 1: form C11 and C21 if n2 != i1 (col = n2)
                if (col /= node%child1%iface) then
                
                    csize = size(node%child1%RtR(node%child1%iface,col)%mat,2)
                
                    ! Construct C11 = S * T^(1)_in2
                    allocate(C11(1:iface_size, 1:csize))
                        
                    call zgemm('n', 'n', iface_size, csize, iface_size, & 
                        cmplx(1.0,0.0,dp), node%S, &
                        iface_size, node%child1%RtR(node%child1%iface,col)%mat, iface_size, &
                        cmplx(0.0,0.0,dp), C11, iface_size)
                        
                    
                    ! Construct C21 = nu * D * S * T^(1)_in2
                    allocate(C21(1:iface_size, 1:csize))
                        
                    call zgemm('n', 'n', iface_size, csize, iface_size, nu, node%D, &
                        iface_size, C11, iface_size, cmplx(0.0,0.0,dp), C21, iface_size) 

                end if
                
                ! col case 2: form C12 and C22 if n2 != i2 (col = n2)
                
                if (col /= node%child2%iface) then
                
                    csize = size(node%child1%RtR(node%child2%iface,col)%mat,2)
                    
                    ! Partially construct C22 = D * T^(2)_in2
                    
                    allocate(C22(1:iface_size, 1:csize))
                    
                    call zgemm('n', 'n', iface_size, csize, iface_size, &
                        cmplx(1.0,0.0,dp), node%D, &
                        iface_size, node%child2%RtR(node%child2%iface,col)%mat, iface_size, &
                        cmplx(0.0,0.0,dp), C22, iface_size)
                    
                    ! Construct C12 = nu * S * D * T^(2)_in2 = nu * S * C22
                    
                    allocate(C12(1:iface_size, 1:csize))
                        
                    call zgemm('n', 'n', iface_size, csize, iface_size, &
                        nu, node%S, &
                        iface_size, C22, iface_size, cmplx(0.0,0.0,dp), C12, iface_size)

                    ! Construct C22 = (D + nu^2 * D*S*D) * T^(2)_in2 = C22 + nu*D*C12                    
                    
                    call zgemm('n', 'n', iface_size, csize, iface_size, &
                        nu, node%D, &
                        iface_size, C12, iface_size, cmplx(1.0,0.0,dp), C22, iface_size) 
                
                end if

                
                ! Inner loop: Allocate and calculate matrix block (row,col)
                do row = 1,2*T%ell_op%d              
                
                    
                    allocate(node%RtR(row,col)%mat(1:pt_face(row),1:pt_face(col)) )

                    ! (1,1) Block : T^(1)_n1n2 + T^(1)_n1i * S * T^(1)_in2
                    
                    if (allocated(C11) .and. (row /= node%child1%iface)) then
                    
                        rsize = size(node%child1%RtR(row,col)%mat,1)
                        csize = size(C11,2)
                        
                        call zgemm('n', 'n', rsize, csize, iface_size, cmplx(1.0,0.0,dp), &
                                node%child1%RtR(row, node%child1%iface)%mat, rsize, &
                                C11, iface_size, cmplx(0.0,0.0,dp), &
                                node%RtR(row,col)%mat(1:rsize,1:csize), &
                                rsize)              
                        
                        node%RtR(row,col)%mat(1:rsize,1:csize) = &
                                node%RtR(row,col)%mat(1:rsize,1:csize) + &
                                node%child1%RtR(row,col)%mat                        
                    
                    end if
                    
                    ! (1,2) Block : - T^(1)_n1i * S * D * T^(2)_in2
                    
                    if (allocated(C12) .and. (row /= node%child1%iface) ) then
                    
                        rsize = size(node%child1%RtR(row,col)%mat,1)
                        csize = size(C12,2)

                        call zgemm('n', 'n', rsize, csize, iface_size, cmplx(-1.0,0.0,dp), &
                                node%child1%RtR(row, node%child1%iface)%mat, rsize, &
                                C12, iface_size, cmplx(0.0,0.0,dp), &
                                node%RtR(row,col)%mat(1:rsize, &
                                (pt_face(col) - csize + 1):pt_face(col)), &
                                rsize)
                    
                    end if
                    
                    ! (2,1) Block : -T^(2)_n1i * D * S * T^(1)_in2
                    
                    if (allocated(C21) .and. (row /= node%child2%iface) ) then
                    
                        rsize = size(node%child2%RtR(row,col)%mat,1) 
                        csize = size(C21,2)
                        
                        call zgemm('n', 'n', rsize, csize, iface_size, cmplx(-1.0,0.0,dp), &
                                node%child2%RtR(row, node%child2%iface)%mat, rsize, &
                                C21, iface_size, cmplx(0.0,0.0,dp), &
                                node%RtR(row,col)%mat(pt_face(row) - rsize +1:pt_face(row), &
                                1:csize), rsize)
                    
                    end if
                    
                    ! (2,2) Block : T^(2)_n1n2 + T^(2)_n1i * (D + D * S * D) * T^(2)_in2
                    
                    if (allocated(C22) .and. (row /= node%child2%iface) ) then
                    
                        rsize = size(node%child2%RtR(row,col)%mat,1)
                        csize = size(C22,2)
                        
                        call zgemm('n', 'n', rsize, csize, iface_size, cmplx(1.0,0.0,dp), &
                                node%child2%RtR(row, node%child2%iface)%mat, rsize, &
                                C22, iface_size, cmplx(0.0,0.0,dp), &
                                node%RtR(row,col)%mat((pt_face(row) - rsize + 1):pt_face(row), &
                                (pt_face(col) - csize + 1):pt_face(col)), &
                                rsize)
                                
                        node%RtR(row,col)%mat((pt_face(row) - rsize + 1):pt_face(row), &
                                (pt_face(col) - csize + 1):pt_face(col)) = &
                                node%RtR(row,col)%mat((pt_face(row) - rsize + 1):pt_face(row), &
                                (pt_face(col) - csize + 1):pt_face(col)) + &
                                node%child2%RtR(row,col)%mat
                    
                    end if
                    
                end do
                
                ! Deallocate intermediate matrices before next step in outer loop                
                if (allocated(C11)) then
                    deallocate(C11, C21)
                end if
                
                if (allocated(C12)) then
                    deallocate(C12,C22)
                end if

            
            end do
        ! Deallocate T^(1)_nn and T^(2)_nn
        if (.not. T%opts%debug%no_delete) then
        
            do row = 1,T%ell_op%d    
                do col = 1,T%ell_op%d
                    ! Delete T^(1)_nn
                    if ((row /= node%child1%iface) .and. (col /= node%child1%iface)) then
                        deallocate(node%child1%RtR(row,col)%mat)
                    end if
                    ! Delete T^(2)_nn
                    if ((row /= node%child2%iface) .and. (col /= node%child2%iface)) then
                        deallocate(node%child2%RtR(row,col)%mat)
                    end if
                    
                end do        
            end do
        
        end if
              
        
        end if
    
    end subroutine compute_RTR



! =================================================================================================




! =================================================================================================
! SUBROUTINE RTR_UPDATE
! Computes one entry X(t,col) of RTR matrix X with outgoing boundary values sigma*u + tau*u_n, 
! where t is implicitly specified by point grid_idx
!
! Dev note: This routine needs to be reworked. Furthermore, possibly should be coupled to the
! related subroutine OUTGOING_BC.

    subroutine RTR_update(grid_idx,face,col,node,ell_op,opts,X)
    
        implicit none
        integer,    dimension(1:3),     intent(in)  :: grid_idx
        integer,                        intent(in)  :: face
        integer,                        intent(in)  :: col
        type(robin_tree_node),          intent(in)  :: node
        type(elliptic_operator),        intent(in)  :: ell_op
        type(solve_opts),               intent(in)  :: opts
        complex(dp),    dimension(:,:), intent(inout)  :: X
        
        logical :: io_stat
        integer,    dimension(1:3)  :: shift
        complex(dp) :: sigma, tau
        integer :: j
        integer :: row1, row2, n_direction
        
        sigma = cmplx(0.0,-1.0,dp) * opts%kb
        tau = cmplx(1.0,0.0,dp)
        
        n_direction = (face+1)/2
        
        if (opts%disc == 'fd') then
                    
            call fgrid2lin(row1, grid_idx, node%ptbox, ell_op%d, opts, io_stat)
            
            ! compute shift direction for 3-pt one-sided normal derivative calculation
            shift = [0,0,0]            
            shift(n_direction) = -1 + 2*modulo(face,2)
            

            
            ! + 0 term
            X(row1,col) = sigma * X(row1,col) + &
                            tau * (8.0d0/3.0d0)/opts%h(n_direction) * X(row1,col)
            
            
            ! -1 term
            
            call fgrid2lin(row2, grid_idx + shift, node%ptbox, ell_op%d, opts, io_stat)
            
            X(row1,col) = X(row1,col) + tau * (-3.0d0)/opts%h(n_direction) * X(row2,col)
            
            ! -2 term
            
            call fgrid2lin(row2, grid_idx + 2*shift, node%ptbox, ell_op%d, opts, io_stat)
            
            X(row1,col) = X(row1,col) + tau * (1.0d0/3.0d0)/opts%h(n_direction) * X(row2,col)
            
        
        end if
    
    end subroutine RTR_update
! =================================================================================================




! =================================================================================================
! SUBROUTINE RTR_SOLVE
! Construct RtR (impedance) operator for a leaf node of the Robin tree.

    subroutine RTR_solve(node, ell_op, opts)
    
        implicit none
        type(robin_tree_node),              intent(inout)   :: node
        type(elliptic_operator),            intent(in)      :: ell_op
        type(solve_opts),                   intent(inout)   :: opts
        
        integer :: row_start, col_start
        integer :: face, row, col, lin_idx
        integer :: j
        integer,    dimension(1:3)  :: grid_idx
        integer,    dimension(1:6)  :: pt_face
        complex(dp),    dimension(:,:), allocatable :: X
        logical, save :: size_display = .false.
        
        ! Dev note: following error check is only for debugging purposes
        if (.not. node%isleaf) then
            stop 'Error in RTR_SOLVE: Routine called on non-leaf node'
        end if
        
        call get_pt_face(pt_face, node%ptbox, ell_op%d)
        
        ! invert elliptic operator
        call elliptic_invert(node%box, ell_op, opts, X)
        
        if (.not. size_display) then
            write(*, '(A, I0, A, I0)') 'Size of matrix at leaf level: ', size(X,1), 'x', size(X,1)
            size_display = .true.
        end if
        
        
        
        ! construct outgoing Robin data -i*kb*u + u_n, operation is done IN PLACE
        ! then must be copied to array blocks in node%RtR
        
        ! dev note: This method only uses the finite differences discretization!!!
        ! dev note: Might be something wrong with this calculation here!
         
        do col = 1, size(X,2)
        do face = 1, 2*ell_op%d
        
            do j = 1, pt_face(face)
            
                lin_idx = j + sum(pt_face(1:face-1))
                call lin2fgrid(lin_idx, grid_idx, node%ptbox, ell_op%d, opts)
                call RTR_update(grid_idx, face, col, node, ell_op, opts, X)
            
            end do
                           
        
        end do 
        end do
        
        ! assign RTR blocks -- this doubles local storage of solution matrix
        
        do col = 1, 2*ell_op%d
            col_start = sum(pt_face(1:col-1)) + 1
            do row = 1, 2*ell_op%d
                row_start = sum(pt_face(1:row-1)) + 1
     
                allocate( node%RtR(row,col)%mat(1:pt_face(row),1:pt_face(col))  )
                
                node%RtR(row,col)%mat = X(row_start:row_start + pt_face(row) - 1, &
                    col_start:col_start + pt_face(col) - 1)
        
            end do
        end do
        
        deallocate(X) ! dev note: unnecessary in modern Fortran but just to be safe
        
        
    
    end subroutine RTR_solve

! =================================================================================================




! =================================================================================================
! SUBROUTINE FORWARD_SUB
! Constructs the outgoing data 'h' on the virtual interfaces contributed by the zero-BC elliptic
! system. Outgoing data takes the form
!       sigma*u + tau*u_n = 0
! Recurses from leaf down to root.

    recursive subroutine forward_sub(T,node,rhs)
    
        implicit none
        type(robin_tree),                   intent(in)      :: T
        type(robin_tree_node),              intent(inout)   :: node
        type(elliptic_rhs),                 intent(in)      :: rhs
        
        integer,        dimension(1:6)              :: pt_face
        complex(dp)                                 :: sigma, tau
        integer                                     :: iface_axis, vec_length
       
        integer :: j
        
        
        ! Dev note: These constants may eventually be replaced by input variables
        
        sigma = cmplx(0.0,-1.0,dp) * T%opts%kb
        tau = cmplx(1.0,0.0,dp)
        
        ! Allocate node%h  
        if ( (node%lvl > 0) .or. (node%isleaf) ) then        
            
            call bndry_vec_malloc(node%h, node, T%ell_op%d, T%opts) 
            
            
        end if
        
        ! allocate node%g
        call get_pt_face(pt_face, node%ptbox, T%ell_op%d)
        call bndry_vec_malloc(node%g, node, T%ell_op%d, T%opts)
        

        
        
        if (node%isleaf) then

        
            ! Perform local solve of the form Au = f, Bu = 0 on leaf
            call leaf_solve(node, T%ell_op, T%opts, rhs%f)
           
            
            call outgoing_bc(node, sigma, tau, T%ell_op%d, T%opts)
            
        
        
        else
        
            ! Recurse down to leaves
            call forward_sub(T, node%child1, rhs)
            call forward_sub(T, node%child2, rhs)
            
            ! Compute initial interface values of g based on update formula
            ! (g_i^(1), g_i^(2)) = inv(M) * (h_i^(1), h_i^(2))
            
            ! subroutine apply_invm(y1, y2, x1, x2, alpha, D, S)
            call apply_invm(node%child1%g(node%child1%iface)%vec, &
                            node%child2%g(node%child2%iface)%vec, &
                            node%child1%h(node%child1%iface)%vec, &
                            node%child2%h(node%child2%iface)%vec, &
                            cmplx(0.0,0.0,dp), &
                            node%D, node%S,node%nu)
            
         
            ! Compute final outgoing data h associated to current node
            ! h = [h_n^(1), h_n^(2)] + [T_ni^(1) g_i^(1), T_ni^(2) g_i^(2)]
            ! This step is ignored if node is root
            
            if (node%lvl > 0) then
            
                ! identify iface_axis (1, 2, or 3)
                
                iface_axis = s_idx(node%child1%iface)
                
                do j = 1,2*T%ell_op%d
                
                    ! along the axis normal to the interface, h only takes values from one child
                    if (s_idx(j) == iface_axis) then
                    
                        ! use child2 values
                        if (j == node%child1%iface) then
                        
                            ! compute h_j = h_j^(2) + T_ji^(2) g_i^(2)
                        
                            node%h(j)%vec = node%child2%h(j)%vec
                            
                            call zgemv('n', size(node%child2%RtR(j, node%child2%iface)%mat, 1), &
                                        size(node%child2%RtR(j, node%child2%iface)%mat, 2), &
                                        cmplx(1.0,0.0,dp), &
                                        node%child2%RtR(j, node%child2%iface)%mat, &
                                        size(node%child2%RtR(j, node%child2%iface)%mat, 1), &
                                        node%child2%g(node%child2%iface)%vec, 1, &
                                        cmplx(1.0,0.0,dp), &
                                        node%h(j)%vec, 1)
                                        
                            
                        ! use child1 values
                        else if (j == node%child2%iface) then
                        
                            ! compute h_j = h_j^(1) + T_ji^(1) g_i^(1)
                        
                            node%h(j)%vec = node%child1%h(j)%vec
                            
                            call zgemv('n', size(node%child1%RtR(j, node%child1%iface)%mat, 1), &
                                        size(node%child1%RtR(j, node%child1%iface)%mat, 2), &
                                        cmplx(1.0,0.0,dp), &
                                        node%child1%RtR(j, node%child1%iface)%mat, &
                                        size(node%child1%RtR(j, node%child1%iface)%mat, 1), &
                                        node%child1%g(node%child1%iface)%vec, 1, &
                                        cmplx(1.0,0.0,dp), &
                                        node%h(j)%vec, 1)
                        
                        else
                            stop 'Error: Could not find correct face values in FORWARD_SUB.'
                        end if
                    
                    ! along axes parallel to the interface, h takes half its values from each child
                    else
                    
                        ! check vector is divisible by 2; otherwise there's a partition problem
                        vec_length = size(node%h(j)%vec,1)
                        if (modulo(vec_length, 2) ==1 ) then
                            stop 'Error: Boundary vector H cannot be partitioned in two.'
                        end if
                        
                        ! initialize h_j = [h_j^(1), h_j^(2)]
                        node%h(j)%vec(1:vec_length/2) = node%child1%h(j)%vec
                        node%h(j)%vec(vec_length/2+1:) = node%child2%h(j)%vec
                        
                        ! compute [T_ji^(1) g_i^(1), T_ji^(2) g_i^(2)] + h_j --> h_j
                        
                        call zgemv('n', size(node%child1%RtR(j, node%child1%iface)%mat, 1), &
                                        size(node%child1%RtR(j, node%child1%iface)%mat, 2), &
                                        cmplx(1.0,0.0,dp), &
                                        node%child1%RtR(j, node%child1%iface)%mat, &
                                        size(node%child1%RtR(j, node%child1%iface)%mat, 1), &
                                        node%child1%g(node%child1%iface)%vec, 1, &
                                        cmplx(1.0,0.0,dp), &
                                        node%h(j)%vec(1:vec_length/2), 1)
                                        
                        call zgemv('n', size(node%child2%RtR(j, node%child2%iface)%mat, 1), &
                                        size(node%child2%RtR(j, node%child2%iface)%mat, 2), &
                                        cmplx(1.0,0.0,dp), &
                                        node%child2%RtR(j, node%child2%iface)%mat, &
                                        size(node%child2%RtR(j, node%child2%iface)%mat, 1), &
                                        node%child2%g(node%child2%iface)%vec, 1, &
                                        cmplx(1.0,0.0,dp), &
                                        node%h(j)%vec(vec_length/2+1:), 1)
                                    
                        
                    
                    end if
                
                end do
                
                
            
            end if
            
            
            
            ! deallocate h-data from children unless debug == .true.
            if (.not. T%opts%debug%no_delete) then

                do j = 1, 2*T%ell_op%d
                    deallocate(node%child1%h(j)%vec)
                    deallocate(node%child2%h(j)%vec)
                end do
            
            end if
            
            if (T%opts%debug%txt_output) then
                write(*, '(A, I0, A)') 'FORWARD_SUB at NODE ID ', node%node_id, ' successful'
            end if
            
        
        end if
        
    end subroutine forward_sub
! =================================================================================================




! =================================================================================================
! SUBROUTINE BACKWARD_SUB
! Constructs the incoming data 'g' on the virtual interfaces contributed by the true boundary
! condition on the root in order to compute the solution contributed by the boundary condition.
! This subroutine is only meant to be called after FORWARD_SUB.
! Recurses from root up to leaf.
    recursive subroutine backward_sub(T, node, rhs)
    
        implicit none
        type(robin_tree),                   intent(in)      :: T
        type(robin_tree_node),              intent(inout)   :: node
        type(elliptic_rhs),                 intent(in)      :: rhs
        
        complex(dp),    dimension(:),   allocatable         :: Tg1,Tg2
        integer                                             :: m,n            
        integer                                             :: iface_axis    
        integer,        dimension(1:6)                      :: ptface
        
        logical,    save    :: bndry_set = .false.    ! Flag for if boundary conditions have been applied
        
        integer :: face
        integer :: j   ! <--- only used in debug
        
        ! The first time routine is called, apply boundary conditions to populate g data for all
        ! faces on the full domain boundary (for nodes at all levels)
        
        if (.not. bndry_set) then
            
            write(*, '(A)') 'Applying boundary conditions...'
            
            ! ! DEV NOTE: THIS ROUTINE DOESN'T RESPECT THE PARTITIONED INDEXING!!!!
            ! call make_bndry_vec(node%g, node, rhs%g, T%ell_op%d, T%opts)
            
            call apply_bc(T, T%root, rhs)
            
            bndry_set = .true.
        
        end if
        
        
        ! For leaf computation, solve Au = 0, Bu = g
        if (node%isleaf) then
        
            call leaf_update(node, T%ell_op, T%opts)

        
        else
        
            ! assign g_n for children
            
            iface_axis = s_idx(node%child1%iface)
            call get_pt_face(ptface, node%ptbox, T%ell_op%d)
            
            do face = 1, 2*T%ell_op%d
            
                if (.not. node%isbndry(face)) then
            
                    ! orthogonal to interface: parent g_n is partitioned to make child data
                    if (s_idx(face) /= iface_axis) then


                        node%child1%g(face)%vec = node%g(face)%vec(1:ptface(face)/2)
                        node%child2%g(face)%vec = node%g(face)%vec(ptface(face)/2 + 1:)

                     
                    ! parallel to interface: parent g_n is assigned to only one child 
                    else

                                            
                        if (face == node%child1%iface) then                    
                            node%child2%g(face)%vec = node%g(face)%vec                        
                        else if (face == node%child2%iface) then                    
                            node%child1%g(face)%vec = node%g(face)%vec                    
                        end if

                    
                    
                    end if
                
                end if
            
            
            
            end do
            
            
            
            ! g_i update computation:
            
        
            allocate(Tg1(1:size(node%child1%g(node%child1%iface)%vec,1)))
            allocate(Tg2(1:size(node%child1%g(node%child1%iface)%vec,1)))
            
           
            ! Construct Tgj = T_in^(j) * g_n^(j) for j = 1,2
            Tg1 = cmplx(0.0,0.0,dp)
            Tg2 = cmplx(0.0,0.0,dp)
            do face = 1, 2*T%ell_op%d           
                if (face /= node%child1%iface) then
                    m = size(node%child1%RtR(node%child1%iface, face)%mat,1)
                    n = size(node%child1%RtR(node%child1%iface, face)%mat,2)
                    
                    call zgemv('n', m,n, cmplx(1.0,0.0, dp), &
                        node%child1%RtR(node%child1%iface, face)%mat, m, &
                        node%child1%g(face)%vec, 1, cmplx(1.0,0.0,dp), &
                        Tg1, 1)

                
                end if
                
                if (face /= node%child2%iface) then
                    m = size(node%child2%RtR(node%child2%iface, face)%mat,1)
                    n = size(node%child2%RtR(node%child2%iface, face)%mat,2)
                    
                    call zgemv('n', m, n, cmplx(1.0,0.0,dp), &
                        node%child2%RtR(node%child2%iface, face)%mat, m, &
                        node%child2%g(face)%vec, 1, cmplx(1.0,0.0,dp), &
                        Tg2, 1)

                
                end if               
            end do
            
            ! update g_i for children based on update formula
            
        
            call apply_invm(node%child1%g(node%child1%iface)%vec, &
                            node%child2%g(node%child2%iface)%vec, &
                            Tg1, Tg2, &
                            cmplx(1.0,0.0,dp), &
                            node%D, node%S, node%nu)
                            

            
            ! deallocate (avoid memory leak!)
            deallocate(Tg1, Tg2)
            
            ! deallocate g data for current (parent) node
            if (.not. T%opts%debug%no_delete) then
                do face = 1, 2*T%ell_op%d
                    deallocate(node%g(face)%vec)
                end do
            end if
            
        
            ! recurse
            call backward_sub(T,node%child1,rhs)
            call backward_sub(T,node%child2,rhs)       
        
        
        end if
        
        
    
    end subroutine backward_sub
! =================================================================================================





! =================================================================================================
! SUBROUTINE BNDRY_VEC_MALLOC
! Allocates the derived type array of type BNDRY_VECTOR associated with a particular Robin tree
! node. The allocation scheme is dependent on the solver opts.
    subroutine bndry_vec_malloc(bndry_vec, node, d, opts)
    
        implicit none
        type(bndry_vector),    dimension(1:6),      intent(inout)   :: bndry_vec
        type(robin_tree_node),                      intent(inout)   :: node
        integer,                                    intent(in)      :: d
        type(solve_opts),                           intent(in)      :: opts
              
        if (opts%disc == 'fd') then
        
            call fd_bndry_vec_malloc(bndry_vec, node, d, opts)
        
        else
        
            stop 'Other discretizations not supported by BNDRY_VEC_MALLOC.'
        
        end if
        
    
    end subroutine bndry_vec_malloc
! =================================================================================================




! =================================================================================================
! SUBROUTINE APPLY_BC
! Uses boundary condition data to populate the g data for all faces along the boundary of the full
! domain of the elliptic problem for nodes at all levels.
    recursive subroutine apply_bc(T, node, rhs)
    
        implicit none
        type(robin_tree),           intent(in)      :: T
        type(robin_tree_node),      intent(inout)   :: node
        type(elliptic_rhs),         intent(in)      :: rhs
        
        integer :: face, m
        
        ! For leaf, directly apply function evaluations from bndry_rhs to populate g data
        if (node%isleaf) then        
        
            call make_bndry_vec(node, rhs%g, T%ell_op%d, T%opts)
        
        ! Non-leaf, merge g data for faces on the boundary
        else if (.not. all(.not. node%isbndry) ) then
        
            call apply_bc(T, node%child1, rhs)
            call apply_bc(T, node%child2, rhs)
            
            do face = 1, 2*T%ell_op%d
            
                if (node%isbndry(face)) then
                
                    m = size(node%g(face)%vec)
                    
                    
                    ! Boundary face is only a part of child1
                    if (face == node%child2%iface) then
                        node%g(face)%vec = node%child1%g(face)%vec
                        
                    ! Boundary face is only a part of child2
                    else if (face == node%child1%iface) then
                        node%g(face)%vec = node%child2%g(face)%vec
                        
                    ! Boundary face is part of both children -- merge g data
                    else
                        node%g(face)%vec(1:m/2) = node%child1%g(face)%vec
                        node%g(face)%vec(m/2+1:) = node%child2%g(face)%vec
                        
                    end if
                
                end if
            
            end do
            
            
        
        end if
    
    
    end subroutine apply_bc
! =================================================================================================





! =================================================================================================
! SUBROUTINE APPLY_INVM
! Applies the inverse of the M matrix to a block vector x = [x_1, x_2], producing a matrix vector
! product y = alpha * y + inv(M) * x.
!
! The inverse of M must be precomputed, such that it is given by its block values D,S, so that
!   inv(M) = [S, -nu*SD; -nu*DS, D + nu*nu*DSD].
!
    subroutine apply_invm(y1, y2, x1, x2, alpha, D, S, nu)
    
        implicit none
        complex(dp),    dimension(:),   intent(inout)   :: y1, y2
        complex(dp),    dimension(:),   intent(in)      :: x1, x2
        complex(dp),                    intent(in)      :: alpha
        complex(dp),    dimension(:,:), intent(in)      :: D,S
        complex(dp),                    intent(in)      :: nu
        
        complex(dp),    dimension(:),   allocatable     :: interm, interm2 ! temp vector variables
        
        
        allocate(interm(1:size(x1,1)), interm2(1:size(x1,1)))
        
        !interm = cmplx(0.0,0.0,dp) ! <--- Should not be necessary
        
        ! Step 1: interm <- S*x1
        
        call zgemv('n', size(S,1), size(S,2), cmplx(1.0,0.0,dp), S, size(S,1), &
                x1, 1, cmplx(0.0,0.0,dp), interm, 1)

                
        ! Step 2: y1 <- alpha*y1 + interm = alpha*y1 + S*x1
        
        y1 = alpha*y1 + interm

        
        ! Step 3: y2 <- alpha*y2 - nu*D*interm = alpha*y2 - nu*D*S*x1
        
        call zgemv('n', size(D,1), size(D,2), -nu, D, size(D,1), &
                interm, 1, alpha, y2, 1)

        
        ! Step 4: interm <- D*x2
            
        call zgemv('n', size(D,1), size(D,2), cmplx(1.0,0.0,dp), D, size(D,1), &
                x2, 1, cmplx(0.0,0.0,dp), interm, 1)

                
        ! Step 5: y2 <- y2 + interm = alpha*y2 - nu*D*S*x1 + D*x2
        y2 = y2 + interm


        ! Step 6: interm <- S*interm = S*D*x2
        
        interm2 = interm
        call zgemv('n', size(S,1), size(S,2), cmplx(1.0,0.0,dp), S, size(S,1), &
                interm2, 1, cmplx(0.0,0.0,dp), interm, 1)

                
        ! Step 7: y1 <- y1 - nu*interm = alpha*x1 + S*x1 - nu*S*D*x2
        
        y1 = y1 - nu*interm


                
        ! Step 8: y2 <- y2 + nu^2*D*interm
        
        call zgemv('n', size(D,1), size(D,2), nu*nu, D, size(D,1), &
                interm, 1, cmplx(1.0,0.0,dp), y2, 1)

                
                
                
        deallocate(interm, interm2) ! <--- Unnecessary in modern Fortran!


    end subroutine apply_invm
! =================================================================================================




! =================================================================================================
! SUBROUTINE PBNDRY2REAL
! Given the linear index on a partitioned boundary, outputs the real-valued coordinate of the point.
! The linear indexing is divided into six contiguous sections, corresponding to the six planar
! faces on a cube. Within each face, indexing is recursively divided in half.
!
! 2D BOUNDARY EXAMPLE:
!
!  11 12 | 15 16
!  09 10 | 13 14
!  -------------
!  02 04 | 06 08
!  01 03 | 05 07
!
! In the above, the highest level interface is the horizontal one. The lower level interface is the
! vertical one. 
!
!    subroutine pbndry2real(T,domain,lin_idx,real_pt)
!

! =================================================================================================










! =================================================================================================
! SUBROUTINE NODE_DIVIDE
!
! Takes input node and determines what the boxes of its children must be by cutting along the
! longest direction of the node's box. Also assigns other geometric parameters such as
! the interface along each child box and whether or not any box faces are subsets of the full
! domain boundary.
!
! Notes: If two edges have equal length, the first one is chosen (in xyz order). d = 2 or 3

	subroutine node_divide(d, node)
	
		implicit none
		
		integer,							intent(in)	    :: d
        type(robin_tree_node),              intent(inout)   :: node
		
		real(dp),	dimension(1:3)						:: widths
		integer											:: j, edge
		
		widths = 0.0d0
        
        widths(1:d) = node%box(2,1:d) - node%box(1,1:d)
		
		! do j = 1,d
			! widths(j) = box_in(2,j) - box_in(1,j)
		! end do
		
        ! Find index of longest direction, along which to cut
		edge = maxloc(widths,1)
		
		node%child1%box = node%box
		node%child2%box = node%box
		
		node%child1%box(2,edge) = node%box(1,edge) + widths(edge) * 0.5d0
		node%child2%box(1,edge) = node%child1%box(2,edge)
        
        ! Identify interface for each child
        node%child1%iface = 2*edge
        node%child2%iface = 2*edge - 1
        
        ! Identify child faces that are on the true domain boundary
        do j = 1, 2*d
            if (node%isbndry(j)) then
            
                if (j /= node%child1%iface) then
                    node%child1%isbndry(j) = .true.
                end if
                
                if (j /= node%child2%iface) then
                    node%child2%isbndry(j) = .true.
                end if
            
            end if
        end do
        
		
	end subroutine node_divide

! =================================================================================================





! =================================================================================================
! FUNCTION BOX_WIDTH

    pure function box_width(box) result(width)
    
        implicit none
        real(dp),   dimension(1:2,1:3), intent(in)  :: box
        real(dp),   dimension(1:3)                  :: width
        
        width = box(2,:) - box(1,:)
    
    end function box_width

! =================================================================================================




! =================================================================================================
! These are just some boring re-indexing functions 
!
! FUNCTION F_IDX
! Maps from 2x3 index to linear index

    pure function f_idx(i,j) result(k)
    
        implicit none
    
        integer,    intent(in)  :: i,j
        integer                 :: k
        
        k = 2*(j-1) + i
    
    end function f_idx
    
! FUNCTION S_IDX
! Maps from linear face index 1,...,6 to the side index 1,...,3
! in other words decoding [x1 x2 y1 y2 z1 z2] into [x y z]

    pure function s_idx(j) result(k)
    
        implicit none
        
        integer,    intent(in)  :: j
        integer                 :: k
        
        k = (j-1)/2 + 1
        
    end function s_idx

! =================================================================================================




end module robin_tree_module
