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
    integer(kint) :: icel
    integer(kint), allocatable :: elem_id(:)
    real(kdouble), allocatable :: u(:)
    real(kdouble), allocatable :: strain(:,:)
    real(kdouble), allocatable :: stress(:,:)
    character(len=:), allocatable :: fname

    call input_param(param)
    call input_mesh(mesh, param)
    call debug('input', debug_flag)

    allocate(mat%A(ndof*mesh%nnode, ndof*mesh%nnode), source=0.0d0)
    allocate(mat%b(ndof*mesh%nnode), source=0.0d0)
    allocate(mat%x(ndof*mesh%nnode), source=0.0d0)

    call get_coef_mat(mesh, param, mat%A)
    call debug('get_coef', debug_flag)

    call get_load(mesh, mat%b)

    call apply_bc_elastic(mesh, mat)
    call debug('apply_bc', debug_flag)

    call solver(mat)
    call debug('solver', debug_flag)

    ! call get_elem_fluxes(mesh, mat%x, fluxes)

    allocate(u(ndof*mesh%nbase_func), source=0.0d0)
    allocate(strain(3, mesh%nelem), source=0.0d0)
    allocate(stress(3, mesh%nelem), source=0.0d0)
    allocate(elem_id(3))
    do icel=1, mesh%nelem
        call get_elem_id(mesh, elem_id, icel)
        u(1:2) = mat%x((elem_id(1)-1)*ndof+1:(elem_id(1)-1)*ndof+2)
        u(3:4) = mat%x((elem_id(2)-1)*ndof+1:(elem_id(2)-1)*ndof+2)
        u(5:6) = mat%x((elem_id(3)-1)*ndof+1:(elem_id(3)-1)*ndof+2)
        call get_strain_and_stress(mesh, param, elem_id, u, strain(:,icel), stress(:,icel))
    enddo

    if(param%flag == 1)then
        fname = 'square'
    elseif(param%flag == 2)then
        fname = 'spherehole'
    endif

    call output_vtk_elastic(mesh%node, mesh%elem, mat%x, mesh%nnode, mesh%nelem, strain, stress, fname)
    call output_data(mesh, mat)
    call debug('output', debug_flag)

    if(mesh%nelem == 2)then
        write(*,*)'<displacement>'
        write(*,*)
        do icel=1, mesh%nnode*ndof
            write(*,*)mat%x(icel)
        enddo
    endif
end program main