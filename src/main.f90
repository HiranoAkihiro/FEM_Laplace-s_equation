program main
    use mod_utils
    use mod_io
    use mod_matrix
    use mod_valid

    implicit none
    type(matdef) :: mat
    type(meshdef) :: mesh
    type(paramdef) :: param
    logical :: debug_flag = .true.
    integer(kint) :: i
    real(kdouble), allocatable :: fluxes(:,:)
    character(len=:), allocatable :: fname

    call get_flag(param)
    call input_mesh(mesh, param)
    call debug('input', debug_flag)

    allocate(mat%A(mesh%nnode, mesh%nnode), source=0.0d0)
    allocate(mat%b(mesh%nnode), source=0.0d0)
    allocate(mat%x(mesh%nnode), source=0.0d0)

    call get_coef_mat(mesh, param, mat%A)
    call debug('get_coef', debug_flag)

    call apply_bc(mesh, mat)
    call debug('apply_bc', debug_flag)

    call solver(mat)
    call debug('solver', debug_flag)

    allocate(fluxes(mesh%nelem,2))
    call get_elem_fluxes(mesh, mat%x, fluxes)

    if(param%flag == 10)then
        fname = 'laplace1'
    elseif(param%flag == 20)then
        fname = 'laplace2'
    endif

    call output_vtk(mesh%node, mesh%elem, mat%x, mesh%nnode, mesh%nelem, fluxes, fname)
    call debug('output', debug_flag)

end program main