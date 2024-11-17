program main
    use mod_utils
    use mod_io
    use mod_matrix
    use mod_valid

    implicit none
    type(matdef) :: mat
    type(meshdef) :: mesh
    logical :: debug_flag = .false.
    integer(kint) :: i
    real(kdouble), allocatable :: fluxes(:,:)

    call input_mesh(mesh)
    call debug('input', debug_flag)

    allocate(mat%A(mesh%nnode, mesh%nnode), source=0.0d0)
    allocate(mat%b(mesh%nnode), source=0.0d0)
    allocate(mat%x(mesh%nnode), source=0.0d0)

    call get_coef_mat(mesh, mat%A)
    call debug('get_coef', debug_flag)

    call get_RHS(mesh, mat%b)

    call apply_bc(mesh, mat)
    call debug('apply_bc', debug_flag)

    call solver(mat)
    call debug('solver', debug_flag)

    allocate(fluxes(mesh%nelem,2))
    call get_elem_fluxes(mesh, mat%x, fluxes)

    call output_vtk(mesh%node, mesh%elem, mat%x, mesh%nnode, mesh%nelem, fluxes)
    call output_data(mesh, mat)
    write(*,*)'***** Normal termination *****'

end program main