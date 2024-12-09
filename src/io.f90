module mod_io
    use mod_utils
    implicit none
    
contains

subroutine get_flag(param)
    implicit none
    type(paramdef) :: param

    open(20,file='input.dat',status='old')
        read(20,*)param%flag
    close(20)
end subroutine get_flag
subroutine input_param(param)
    implicit none
    type(paramdef) :: param
    
    open(20,file='input.dat',status='old')
        read(20,*)param%E
        read(20,*)param%nu
        read(20,*)param%flag
    close(20)
    param%elastic = .true.
end subroutine input_param

subroutine input_mesh(mesh, param)
    type(meshdef) :: mesh
    type(paramdef) :: param
    integer(kint) :: i, n

    open(20,file='node.dat',status='old')
    open(21,file='elem.dat',status='old')
    open(22,file='bc.dat',status='old')
    if(param%elastic)open(23,file='load.dat',status='old')

        read(20,*)mesh%nnode, ndof
        allocate(mesh%node(ndof, mesh%nnode), source=0.0d0)
        do i=1, mesh%nnode
            read(20,*)mesh%node(1,i), mesh%node(2,i)
        enddo

        read(21,*)mesh%nelem, mesh%nbase_func
        allocate(mesh%elem(mesh%nbase_func, mesh%nelem))
        do i=1, mesh%nelem
            read(21,*)mesh%elem(1,i), mesh%elem(2,i), mesh%elem(3,i)
        enddo

        read(22,*)mesh%nbc, n
        allocate(mesh%bc(mesh%nbc), source=0.0d0)
        allocate(mesh%bc_id(mesh%nbc))
        allocate(mesh%bc_dof(mesh%nbc))
        do i=1, mesh%nbc
            read(22,*)mesh%bc_id(i), mesh%bc_dof(i), mesh%bc(i)
        enddo

        if(param%elastic)then
            read(23,*)mesh%nload, n
            allocate(mesh%load(mesh%nload), source=0.0d0)
            allocate(mesh%load_id(mesh%nload))
            allocate(mesh%load_dof(mesh%nload))
            do i=1, mesh%nload
                read(23,*)mesh%load_id(i), mesh%load_dof(i), mesh%load(i)
            enddo
        endif
    close(20)
    close(21)
    close(22)
    if(param%elastic)close(23)
end subroutine input_mesh
    
subroutine output_vtk(Nm, elements, x, nn, ne, fluxes, fname)
    implicit none
    integer(4) :: i
    integer(4) :: nn, ne
    real(8), allocatable :: Nm(:,:)
    real(8), allocatable :: x(:)
    integer(4), allocatable :: elements(:,:)
    real(8), allocatable :: fluxes(:,:)
    character(len=*) :: fname
    character(len=:), allocatable :: dir_name

    dir_name = '../visual/'//fname//'.vtk'
    open(20,file=dir_name,status="replace")
        write(20,'(a)')'# vtk DataFile Version 4.1'
        write(20,'(a)')'Laplace Equation Solution with Fluxes'
        write(20,'(a)')'ASCII'
        write(20,'(a)')'DATASET UNSTRUCTURED_GRID'
        write(20,*)''
        
        write(20,'(a,a,i0,a,a)')'POINTS',' ',nn,' ','float'
        do i=1,nn
            write(20,'(3f12.6)')Nm(1,i),Nm(2,i),0.0d0
        enddo
        write(20,*)''
        
        write(20,'(a,a,i0,a,i0)')'CELLS',' ',ne,' ',ne*4
        do i=1,ne
            write(20,'(i0,3i6)')3,elements(1,i)-1,elements(2,i)-1,elements(3,i)-1
        enddo
        write(20,*)''
        
        write(20,'(a,a,i0)')'CELL_TYPES',' ',ne
        do i=1,ne
            write(20,'(i0)')5
        enddo
        write(20,*)''
        
        write(20,'(a,a,i0)')'POINT_DATA',' ',nn
        write(20,'(a)')'SCALARS Solution float 1'
        write(20,'(a)')'LOOKUP_TABLE default'
        do i=1,nn
            write(20,'(f12.6)')x(i)
        enddo
        write(20,*)''

        ! フラックスデータの出力（要素ごと）
        write(20,'(a,a,i0)')'CELL_DATA',' ',ne
        
        ! X方向フラックス
        write(20,'(a)')'SCALARS Flux_X float 1'
        write(20,'(a)')'LOOKUP_TABLE default'
        do i=1,ne
            write(20,'(f12.6)')fluxes(i,1)
        enddo
        write(20,*)''
        
        ! Y方向フラックス
        write(20,'(a)')'SCALARS Flux_Y float 1'
        write(20,'(a)')'LOOKUP_TABLE default'
        do i=1,ne
            write(20,'(f12.6)')fluxes(i,2)
        enddo

    close(20)
end subroutine output_vtk

subroutine output_vtk_elastic(Nm, elements, x, nn, ne, strain, stress, fname)
    implicit none
    integer(4) :: i
    integer(4) :: nn, ne
    real(8), allocatable :: Nm(:,:)
    real(8), allocatable :: x(:)
    integer(4), allocatable :: elements(:,:)
    real(kdouble), intent(in) :: strain(:,:)
    real(kdouble), intent(in) :: stress(:,:)
    character(len=*), intent(in) :: fname
    character(len=:), allocatable :: dir_name

    dir_name = '../visual/'//fname//'.vtk'
    open(20,file=dir_name,status="replace")
        write(20,'(a)')'# vtk DataFile Version 4.1'
        write(20,'(a)')'Laplace Equation Solution with Fluxes'
        write(20,'(a)')'ASCII'
        write(20,'(a)')'DATASET UNSTRUCTURED_GRID'
        write(20,*)''
        
        write(20,'(a,a,i0,a,a)')'POINTS',' ',nn,' ','float'
        do i=1,nn
            write(20,'(3f12.6)')Nm(1,i),Nm(2,i),0.0d0
        enddo
        write(20,*)''
        
        write(20,'(a,a,i0,a,i0)')'CELLS',' ',ne,' ',ne*4
        do i=1,ne
            write(20,'(i0,3i6)')3,elements(1,i)-1,elements(2,i)-1,elements(3,i)-1
        enddo
        write(20,*)''
        
        write(20,'(a,a,i0)')'CELL_TYPES',' ',ne
        do i=1,ne
            write(20,'(i0)')5
        enddo
        write(20,*)''
        
        write(20,'(a,a,i0)')'POINT_DATA',' ',nn
        write(20,'(a)')'VECTORS disp float'
        do i=1,nn
            write(20,'(f12.6)')x(ndof*(i-1)+1),x(ndof*(i-1)+2),0.0d0
        enddo
        write(20,*)''

        write(20,'(a,a,i0)')'CELL_DATA',' ',ne
        write(20,'(a)')'SCALARS stress_x float 1'
        write(20,'(a)')'LOOKUP_TABLE default'
        do i=1,ne
            write(20,'(f12.6)')stress(1,i)
        enddo

        write(20,'(a)')'SCALARS stress_y float 1'
        write(20,'(a)')'LOOKUP_TABLE default'
        do i=1,ne
            write(20,'(f12.6)')stress(2,i)
        enddo

        write(20,'(a)')'SCALARS stress_t float 1'
        write(20,'(a)')'LOOKUP_TABLE default'
        do i=1,ne
            write(20,'(f12.6)')stress(3,i)
        enddo

        write(20,'(a)')'SCALARS strain_x float 1'
        write(20,'(a)')'LOOKUP_TABLE default'
        do i=1,ne
            write(20,'(f12.6)')strain(1,i)
        enddo

        write(20,'(a)')'SCALARS strain_y float 1'
        write(20,'(a)')'LOOKUP_TABLE default'
        do i=1,ne
            write(20,'(f12.6)')strain(2,i)
        enddo

        write(20,'(a)')'SCALARS strain_t float 1'
        write(20,'(a)')'LOOKUP_TABLE default'
        do i=1,ne
            write(20,'(f12.6)')strain(3,i)
        enddo
        write(20,*)''

    close(20)
end subroutine output_vtk_elastic

subroutine output_data(mesh, mat)
    implicit none
    type(meshdef), intent(in) :: mesh
    type(matdef), intent(in) :: mat
    logical :: dir_exists
    character(len=:), allocatable :: dir_name
    character(len=100) :: cur_dir, str, command
    integer(kint) :: i, nx, ierr
    real(kdouble) :: var

    var = sqrt(dble(mesh%nnode))
    nx = int(var)
    call getcwd(cur_dir, ierr)
    write(str,'(i0)')(nx-1)**2
    cur_dir = 'poisson'
    dir_name = trim(cur_dir)//trim(str)

    call chdir("../data", ierr)
    inquire(file=dir_name, exist=dir_exists)
    if(.not.dir_exists)then
        command = 'mkdir '//dir_name
        call system(trim(command))
        if(ierr/=0)stop 'mkdir'
    endif

    open(20,file=dir_name//'/on_axis.dat',status='replace')
    open(21,file=dir_name//'/center.dat',status='replace')
        do i=nx,2,-1
            write(20,'(g0.17,a,g0.17)')-mesh%node(1,i),' ',mat%x(i)
        enddo
        do i=1,nx
            write(20,'(g0.17,a,g0.17)')mesh%node(1,i),' ',mat%x(i)
        enddo

        write(21,'(i0,a,g0.17)')(nx-1)**2,' ',mat%x(1)
    close(21)
    close(20)
end subroutine output_data

subroutine get_strain_and_stress(mesh, param, elem_id, u, strain, stress)
    implicit none
    type(meshdef) :: mesh
    type(paramdef) :: param
    real(kdouble), intent(in) :: u(:)
    real(kdouble), intent(inout) :: strain(:)
    real(kdouble), intent(inout) :: stress(:)
    integer(kint), intent(in) :: elem_id(:)
    real(kdouble), allocatable :: B_mat(:,:), D_mat(:,:)
    real(kdouble) :: det
    real(kdouble) :: x(3), y(3)
    integer(kint) :: i, j

    x = 0.0d0
    y = 0.0d0
    det = 0.0d0
    allocate(B_mat(3,6) ,source=0.0d0)
    allocate(D_mat(3,3) , source=0.0d0)

    do i=1,3
        x(i) = mesh%node(1, elem_id(i))
        y(i) = mesh%node(2, elem_id(i))
    enddo

    det = (x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1))

    B_mat(1,1) = (y(2) - y(3))/det
    B_mat(1,3) = (y(3) - y(1))/det
    B_mat(1,5) = (y(1) - y(2))/det
    B_mat(2,2) = (x(3) - x(2))/det
    B_mat(2,4) = (x(1) - x(3))/det
    B_mat(2,6) = (x(2) - x(1))/det
    B_mat(3,1) = (x(3) - x(2))/det
    B_mat(3,3) = (x(1) - x(3))/det
    B_mat(3,5) = (x(2) - x(1))/det
    B_mat(3,2) = (y(2) - y(3))/det
    B_mat(3,4) = (y(3) - y(1))/det
    B_mat(3,6) = (y(1) - y(2))/det

    strain = 0.0d0
    do i=1,3
        do j=1,6
            strain(i) = strain(i) + B_mat(i,j)*u(j)
        enddo
    enddo

    D_mat(1,1) = param%E/(1.0d0-param%nu**2.0d0)
    D_mat(1,2) = param%E*param%nu/(1.0d0-param%nu**2.0d0)
    D_mat(2,1) = param%E*param%nu/(1.0d0-param%nu**2.0d0)
    D_mat(2,2) = param%E/(1.0d0-param%nu**2.0d0)
    D_mat(3,3) = param%E/(2.0d0*(1.0d0+param%nu))

    stress = 0.0d0
    do i=1,3
        do j=1,3
            stress(i) = stress(i) + D_mat(i,j)*strain(j)
        enddo
    enddo
end subroutine get_strain_and_stress
end module mod_io