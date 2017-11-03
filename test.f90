! I compile with:
! gfortran -march=native -O3 -c derived_type_module.f90 fd_module.f90 local_solver_module.f90 time_module.f90 robin_tree_module.f90 robin_output_module.f90
!
! gfortran -march=native -O3 -o output test.f90 derived_type_module.o fd_module.o local_solver.o time_module.o robin_tree_module.o robin_output_module.o -L C:\OpenBLAS-v0.2.15-Win64-int32\lib -lopenblas

program test

    use derived_type_module
    use local_solver_module
    use robin_tree_module
    use robin_output_module
    
    

    implicit none
    
    integer,    parameter   :: dp = kind(1.0d0) ! double precision kind specification
    
    type(elliptic_operator)             :: my_op
    type(solve_opts)                    :: my_opts
    type(bndry_rhs),    dimension(1:6)  :: g
    
    type(robin_tree)                    :: my_tree
    type(elliptic_rhs)                  :: rhs
    
    type(leaf_solution),    dimension(:),   allocatable :: u_out
    
    real(dp),   dimension(:),   allocatable :: errors    
    
    integer                             :: n, L, j, jj, outputnode, rowout
    
    
    ! Construct elliptic operator
    my_op%d = 2
    my_op%domain(:,1) = [0.0d0,1.0d0]
    my_op%domain(:,2) = [0.0d0,1.0d0]
    my_op%domain(:,3) = [0.0d0,1.0d0]
    my_op%coeff => coeff
    my_op%k => wavespeed
    
    ! Set boundary condition function
    
    print *, dp
    
    if (my_op%d == 2) then
        g(1)%g_face => zero2D
        g(2)%g_face => one2D
        g(3)%g_face => zero2D
        g(4)%g_face => one2D
    
    else if (my_op%d == 3) then
        
        g(1)%g_face => x0bc
        g(2)%g_face => x1bc
        g(3)%g_face => y0bc
        g(4)%g_face => y1bc
        g(5)%g_face => z0bc
        g(6)%g_face => z1bc
    end if
    
    ! Select opts
    n = 512
    my_opts%h_tgt = 1.0d0/real(n,dp)
    my_opts%kb = 1.0d0
    my_opts%debug%txt_output = .false.
    my_opts%debug%no_delete = .false.
    
    ! Build tree
    my_tree%bndry_params(1,:) = [cmplx(0.0,1.0,dp)*my_opts%kb, cmplx(1.0,0.0,dp)]
    my_tree%bndry_params(2,:) = [cmplx(0.0,-1.0,dp)*my_opts%kb, cmplx(1.0,0.0,dp)]
    
    ! Set RHS
    
    rhs%d = my_op%d
    rhs%f => zerointerior
    rhs%g = g
    
    
    
    ! Invoke elliptic solver
    L = 8
    call my_tree%factor(L,my_op, my_opts)
    call my_tree%solve(rhs)
    call my_tree%inspect()
    
    ! Check solution
    call get_solution(my_tree,u_out)

    
    if (my_op%d ==2) then
        call check_solution(my_tree%ell_op%d, u_out, true_sol2D, my_tree%opts, errors, 'Inf', .true.)
    else if (my_op%d ==3) then
        call check_solution(my_tree%ell_op%d, u_out, truesol_nonsym, my_tree%opts, errors, 'Inf', .true.)
    end if
    
    print *, 'Error: ', maxval(errors)
  
    ! ! Write impedance matrices to file
    
    ! outputnode = 2
    ! rowout = 1
    ! open(unit = 1, file = "T_comp_real.txt", action = "write", status = "replace")
    ! open(unit = 2, file = "T_comp_imag.txt", action = "write", status = "replace")
    ! do j = 1, size(my_tree%node_list(outputnode)%ptr%RTR(rowout,1)%mat,1)
        ! write(1, *) realpart(my_tree%node_list(outputnode)%ptr%RTR(rowout,1)%mat(j,:)), &
            ! realpart(my_tree%node_list(outputnode)%ptr%RTR(rowout,2)%mat(j,:)), &
            ! realpart(my_tree%node_list(outputnode)%ptr%RTR(rowout,3)%mat(j,:)), &
            ! realpart(my_tree%node_list(outputnode)%ptr%RTR(rowout,4)%mat(j,:))
        ! write(2, *) imagpart(my_tree%node_list(outputnode)%ptr%RTR(rowout,1)%mat(j,:)), &
            ! imagpart(my_tree%node_list(outputnode)%ptr%RTR(rowout,2)%mat(j,:)), &
            ! imagpart(my_tree%node_list(outputnode)%ptr%RTR(rowout,3)%mat(j,:)), &
            ! imagpart(my_tree%node_list(outputnode)%ptr%RTR(rowout,4)%mat(j,:))
    ! end do    
    
    
    contains
    
    pure function coeff(x) result(c)
    
        implicit none
        
        real(dp),   dimension(:),   intent(in)  :: x
        complex(dp),    dimension(1:3,1:3)      :: c
        
        c(:,1) = [1,0,0]
        c(:,2) = [0,1,0]
        c(:,3) = [0,0,1]
    
    end function coeff
    
    pure function wavespeed(x) result(r)
    
        implicit none
        
        real(dp),   dimension(:),   intent(in)  :: x
        complex(dp)                             :: r
        
        r = my_opts%kb
    
    end function wavespeed
    
! =================================================================================================

    pure function zerointerior(x) result(r)
    
        implicit none
        real(dp),   dimension(:),   intent(in)  :: x
        complex(dp)                             :: r
        
        r = cmplx(0.0,0.0,dp)
    
    end function zerointerior




! =================================================================================================

    
    pure function zero2D(x) result(r)
    
        implicit none
        real(dp),   dimension(:),   intent(in)  :: x
        complex(dp)                             :: r
        
        r = cmplx(0.0, 1.0, dp) * dcos(x(1) / sqrt(2.0d0))
        
    end function zero2D
    
    pure function zero3D(x) result(r)
    
        implicit none
        real(dp),   dimension(:),   intent(in)  :: x
        complex(dp)                             :: r
        
        r = cmplx(0.0, 1.0, dp) * dcos(x(1) / sqrt(3.0d0)) * dcos(x(2)/sqrt(3.0d0))
        
    end function zero3D
    
    pure function one2D(x) result(r)
    
        implicit none
        real(dp),   dimension(:),   intent(in)  :: x
        complex(dp)             :: r
        
        r = cmplx(0.0, 1.0, dp) * dcos(1.0d0/sqrt(2.0d0)) * dcos(x(1)/sqrt(2.0d0))  & 
            -1.0d0/sqrt(2.0d0) * dsin(1.0d0/sqrt(2.0d0)) * dcos( x(1) / sqrt(2.0d0)) 
        
    end function one2D
    
    pure function one3D(x) result(r)
    
        implicit none
        real(dp),   dimension(:),   intent(in)  :: x
        complex(dp)             :: r
        
        r = cmplx(0.0, 1.0, dp) * cos(1.0/sqrt(3.0)) * cos(x(1)/sqrt(3.0)) * cos(x(2)/sqrt(3.0)) & 
            -1.0/sqrt(3.0) * sin(1/sqrt(3.0)) * cos( x(1) / sqrt(3.0)) * cos( x(2) / sqrt(3.0))
        
    end function one3D
    
    pure function true_sol2D(x) result(r)
    
        implicit none
        real(dp),   dimension(1:3), intent(in)  :: x
        complex(dp)             :: r
        
        r = dcos(x(1)/sqrt(2.0)) * dcos(x(2)/sqrt(2.0))
    
    end function true_sol2D
    
    pure function true_sol3D(x) result(r)
        
        implicit none
        real(dp),   dimension(1:3), intent(in)  :: x
        complex(dp)             :: r
        
        r = cos(x(1)/sqrt(3.0d0)) * cos(x(2)/sqrt(3.0d0)) * cos(x(3)/sqrt(3.0d0))
    
    end function true_sol3D
    
! =================================================================================================
    
    
    pure function x0bc(x) result(r)
    
        implicit none
        real(dp),   dimension(:),   intent(in)  :: x
        complex(dp)                             :: r
        
        r = cmplx(0.0,1.0,dp) * sin(x(1)/sqrt(3.0d0)) * sin(x(2)/sqrt(3.0d0) + 1.0d0)
    
    end function x0bc
    
    
    
    pure function x1bc(x) result(r)
    
        implicit none
        real(dp),   dimension(:),   intent(in)  :: x
        complex(dp)                             :: r
        
        r = cmplx(0.0,1.0,dp) * cos(1.0d0/sqrt(3.0d0)) * sin(x(1)/sqrt(3.0d0)) * sin(x(2)/sqrt(3.0d0) + 1.0d0) &
            - 1.0d0/sqrt(3.0d0) * sin(1.0d0/sqrt(3.0d0)) * sin(x(1)/sqrt(3.0d0)) * sin(x(2)/sqrt(3.0d0) + 1.0d0)
            
    end function x1bc
    
    
    pure function y0bc(x) result(r)
    
        implicit none
        real(dp),   dimension(:),   intent(in)  :: x
        complex(dp)                             :: r
        
        r = -1.0d0/sqrt(3.0d0) * cos(x(1)/sqrt(3.0d0)) * sin(x(2)/sqrt(3.0d0) + 1.0d0)
    
    end function y0bc
    
    pure function y1bc(x) result(r)
    
        implicit none
        real(dp),   dimension(:),   intent(in)  :: x
        complex(dp)                             :: r
        
        r = cmplx(0.0,1.0,dp) * cos(x(1)/sqrt(3.0d0)) * sin(1.0d0/sqrt(3.0d0)) * sin(x(2)/sqrt(3.0d0) + 1.0d0) &
            + 1.0d0/sqrt(3.0d0) * cos(x(1)/sqrt(3.0d0)) * cos(1.0d0/sqrt(3.0d0)) * sin(x(2)/sqrt(3.0d0) + 1.0d0)
    
    end function y1bc
    
    pure function z0bc(x) result(r)
    
        implicit none
        real(dp),   dimension(:),   intent(in)  :: x
        complex(dp)                             :: r
        
        r = cmplx(0.0,1.0,dp) * cos(x(1)/sqrt(3.0d0)) * sin(x(2)/sqrt(3.0d0)) * sin(1.0d0) &
            - 1.0d0/sqrt(3.0d0) * cos(x(1)/sqrt(3.0d0)) * sin(x(2)/sqrt(3.0d0)) * cos(1.0d0)
    
    end function z0bc
    
    pure function z1bc(x) result(r)
    
        implicit none
        real(dp),   dimension(:),   intent(in)  :: x
        complex(dp)                             :: r
        
        r = cmplx(0.0,1.0,dp) * cos(x(1)/sqrt(3.0d0)) * sin(x(2)/sqrt(3.0d0)) * sin(1.0d0/sqrt(3.0d0)+1.0d0) &
            + 1.0d0/sqrt(3.0d0) * cos(x(1)/sqrt(3.0d0)) * sin(x(2)/sqrt(3.0d0)) * cos(1.0d0/sqrt(3.0d0) + 1.0d0)
    
    end function z1bc
    
    pure function truesol_nonsym(x) result(r)
    
        implicit none
        real(dp),   dimension(1:3), intent(in)  :: x
        complex(dp)             :: r
        
        r = cos(x(1)/sqrt(3.0d0)) * sin(x(2)/sqrt(3.0d0)) * sin(x(3)/sqrt(3.0d0) + 1.0d0)
    
    end function truesol_nonsym


end program test