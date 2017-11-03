# robinsolver

By Christopher A. Wong

This code solves the Helmholtz equation with variable potential on a rectangular domain based on
the hierarchical Poincare-Steklov method. It is written in Fortran 2008. This is purely a research
code and likely has lots of issues stemming from bad practices.

References:

-- Gillman, Barnett, Martinsson. "A spectrally accurate direct solution technique for
frequency-domain scattering problems with variable media".

-- Liu, de Hoop, Xia. "Interconnected hierarchical rank-structured method for direct elliptic
solution and damped high-frequency Helmholtz preconditioning".

Note: Everything is a work-in-progress. Lots of incomplete subroutines/features.

I. Compilation
==============

There is currently no makefile. Use your Fortran compiler to compile object files for the following,
in order:

    derived_type_module.f90
    time_module.f90
    fd_module.f90
    local_solver_module.f90
    robin_tree_module.f90
    robin_output_module.f90
    
Your source file that invokes the elliptic solver must be linked to the above object files, as well
as to BLAS.

II. Input Format
================

Any code that uses these routines must include the appropriate modules. Write
    use derived_type_module
    use local_solver_module
    use robin_tree_module
    use robin_output_module
The primary data structure in this code is the derived type ROBIN_TREE. At the beginning of your
code, declare:

    type(robin_tree) :: my_tree
    
Next you must declare an elliptic operator, which is done by declaring:

    type(elliptic_operator) :: my_operator
    
Solver options must also be declared by

    type(solve_opts)  :: my_opts
    
To begin, first specify what MY_OPERATOR is. The derived type ELLIPTIC_OPERATOR has the following
attributes to be set:

    d         --- integer; Euclidean dimension; must be d = 2 or 3
    domain    --- real double, dimension(1:2,1:3); The boundaries of the rectangular domain, written as [x1 y1 z1; x2 y2 z2]
    k         --- Pointer to complex-valued potential function
    robin_coeff --- complex double, dimension(1:6,1:2); Of the form [a1 a2; ... ]^T, for Robin condition* on each face: (a1*u + a2*u_n)

Next, we set the options. The main option to set is the discretizations. Currently not fully
supported, but for finite differences you can set the grid size with

    my_opts%h_tgt = 0.01
    
for example, in order to have a grid size of at most 1/100. The solver may adjust this value
in order to achieve an integer number of step sizes within the domain. Currently you will also need
to set

    my_opts%kb = 1.0d0
    
where this value is the value of 'k' used in the impedance boundary data. In the future this
will not be used.

III. Invoking Solver
====================

Now that the ROBIN_TREE is initialized and the ELLIPTIC_OPERATOR declared, we can perform the
factorization. This is done by

    call my_tree%factor(L, my_operator, my_opts)
    
where L specifies the number of levels in the hierarchical factorization. This subroutine
will populate the ROBIN_TREE data type with all the relevant intermediate matrices required to
build the factorization. Depending on the size of your problem, this step may take a long time and 
consume a large amount of memory. 

To solve a specific elliptic problem, we must specify a right-hand-side for the elliptic system.
This is done by declaring

    type(elliptic_rhs)    :: my_rhs
    
The derived type ELLIPTIC_RHS has the following attributes to be set:

    d       --- integer; Euclidean dimension of the problem
    f       --- pointer to a scalar function for the interior portion of the RHS
    g       --- TYPE(BNDRY_RHS), dimension(1:6)

Note that ELLIPTIC_RHS depends on a further subtype BNDRY_RHS(1:6). Each of the six elements of this
array is itself a pointer to a scalar function indicating the boundary condition on each of the six
faces of the rectangular domain. If 2D, only the first four elements are associated.

With the right-hand-side determined, we can invoke the elliptic solve with

    call my_tree%solve(my_rhs)
    
This subroutine will populate the ROBIN_TREE type with the solution vectors local to each leaf node.
To extract the solution, we must lastly declare yet another variable:

    type(leaf_solution),  dimension(:),   allocatable     :: my_u
    
Then we can obtain the numerical solution vectors by writing

    call get_solution(my_tree, my_u)
    
This allocates the array and populates it with the numerical solution vectors. Each element my_u(j)
has two attributes:

    box     --- The subdomain.
    u       --- The numerical solution vector supported on that subdomain.    
