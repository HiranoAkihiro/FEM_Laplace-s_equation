program main
    implicit none
    integer(4), parameter :: kint = 4
    integer(4), parameter :: kdouble = 8
    integer(4) :: ndof

    type matdef
        real(8), allocatable :: A(:,:)
        real(8), allocatable :: b(:)
        real(8), allocatable :: x(:)
    endtype

    type meshdef
        integer(4) :: nnode
        integer(4) :: nelem
        integer(4) :: nbc
        integer(4) :: nbase_func
        real(8), allocatable :: node(:,:)
        integer(4), allocatable :: elem(:,:)
        integer(4), allocatable :: bc_id(:)
        real(8), allocatable :: bc(:)
    endtype
    type(matdef) :: mat
    type(meshdef) :: mesh
    integer(kint) :: i

    call input_mesh(mesh)
    call debug('input')

    allocate(mat%A(mesh%nnode, mesh%nnode), source=0.0d0)
    allocate(mat%b(mesh%nnode), source=0.0d0)
    allocate(mat%x(mesh%nnode), source=0.0d0)

    call get_coef_mat(mesh, mat%A)
    call debug('get_coef')
    do i=1, mesh%nnode
        write(*,*)mat%A(:,i)
    enddo
    call apply_bc(mesh, mat)
    call debug('apply_bc')

    call solver(mat)
    call debug('solver')

    call output_vtk(mesh%node, mesh%elem, mat%x, mesh%nnode, mesh%nelem)

    write(*,*)'Normal termination'

contains

subroutine input_mesh(mesh)
    type(meshdef) :: mesh
    integer(kint) :: i, n

    open(20,file='node.dat',status='old')
    open(21,file='elem.dat',status='old')
    open(22,file='bc.dat',status='old')

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

        read(22,*)mesh%nbc
        allocate(mesh%bc(mesh%nbc), source=0.0d0)
        allocate(mesh%bc_id(mesh%nbc))
        do i=1, mesh%nbc
            read(22,*)mesh%bc_id(i), n, mesh%bc(i)
        enddo
    close(20)
    close(21)
    close(22)
end subroutine input_mesh

subroutine apply_bc(mesh, mat)
    implicit none
    type(meshdef) :: mesh
    type(matdef) :: mat
    integer(kint) :: i, j, id

    do j=1, mesh%nbc
        id = mesh%bc_id(j)
        do i=1, mesh%nnode
            mat%A(i,id) = 0.0d0
        enddo

        do i=1, mesh%nnode
            mat%A(id,i) = 0.0d0
        enddo
        mat%A(id,id) = 1.0d0
        mat%b(id) = mesh%bc(j)
    enddo
end subroutine

subroutine solver(mat)
    type(matdef) :: mat
    real(kdouble), allocatable :: A(:,:)
    real(kdouble), allocatable :: b(:)
    real(kdouble), allocatable :: x(:)
    integer(kint) :: n, i

    n = size(mat%A, dim=1)
    allocate(A(n,n), source=0.0d0)
    allocate(b(n), source=0.0d0)
    allocate(x(n), source=0.0d0)

    A = mat%A
    b = mat%b

    do i=1, n
        write(*,*)A(:,i)
    enddo
    write(*,*)
    write(*,*)b

    call LU_decomposition(n, A, b, x)
    mat%x = x
end subroutine solver

subroutine LU_decomposition(n, A, b, x)
    implicit none
    integer(4), intent(in) :: n !> degrees of freedom
    real(8), intent(in) :: A(:,:) !> coefficient matrix
    real(8), intent(in) :: b(:) !> right-hand side vector
    real(8), intent(out) :: x(:) !> solution vector
    
    real(8), allocatable :: L(:,:) !> lower triangular matrix
    real(8), allocatable :: U(:,:) !> upper triangular matrix
    real(8), allocatable :: y(:) !> intermediate vector
    
    integer(4) :: i, j, k
    
    ! Allocate memory
    allocate(L(n,n), source = 0.0d0)
    allocate(U(n,n), source = 0.0d0)
    allocate(y(n), source = 0.0d0)
    
    ! LU decomposition
    do i = 1, n
        ! Calculate upper triangular part of U
        do j = i, n
            U(i, j) = A(i, j)
            if (i > 1) then
                do k = 1, i-1
                    U(i, j) = U(i, j) - L(i, k) * U(k, j)
                end do
            end if
        end do
        
        ! Error check: diagonal element should not be zero
        if (abs(U(i,i)) < 1.0d-10) then
            print *, "Error: Zero diagonal element in LU decomposition"
            stop
        end if
        
        ! Calculate lower triangular part of L
        L(i, i) = 1.0d0
        do j = i+1, n
            L(j, i) = A(j, i)
            if (i > 1) then
                do k = 1, i-1
                    L(j, i) = L(j, i) - L(j, k) * U(k, i)
                end do
            end if
            L(j, i) = L(j, i) / U(i, i)
        end do
    end do
    
    ! Forward substitution L * y = b
    do i = 1, n
        y(i) = b(i)
        if (i > 1) then
            do j = 1, i-1
                y(i) = y(i) - L(i, j) * y(j)
            end do
        end if
    end do
    
    ! Backward substitution U * x = y
    do i = n, 1, -1
        x(i) = y(i)
        if (i < n) then
            do j = i+1, n
                x(i) = x(i) - U(i, j) * x(j)
            end do
        end if
        x(i) = x(i) / U(i, i)
    end do
end subroutine LU_decomposition

subroutine debug(char)
    implicit none
    character(len=*), intent(in) :: char
    character(len=:), allocatable :: name
    
    name = '***** done '//trim(char)//' *****'
    write(*,*)name
end subroutine debug

subroutine get_coef_mat(mesh, mat)
    implicit none
    type(meshdef) :: mesh
    real(kdouble), intent(inout) :: mat(:,:)
    integer(kint) :: i
    real(kdouble), allocatable :: e_mat(:,:)
    integer(kint), allocatable :: elem_id(:)

    allocate(e_mat(mesh%nbase_func, mesh%nbase_func), source=0.0d0)
    allocate(elem_id(mesh%nbase_func), source=0)

    do i=1, mesh%nelem
        call get_elem_id(mesh, elem_id, i)
        call get_element_mat(mesh, elem_id, e_mat)
        write(*,*)i
        write(*,*)e_mat
        write(*,*)
        call add_global_mat(elem_id, e_mat, mat)
    enddo
end subroutine

subroutine get_elem_id(mesh, elem_id, icel)
    implicit none
    type(meshdef) :: mesh
    integer(kint), intent(in) :: icel
    integer(kint), intent(inout) :: elem_id(:)

    elem_id = 0
    elem_id(:) = mesh%elem(:,icel)
end subroutine get_elem_id

subroutine get_element_mat(mesh, elem_id, e_mat)
    implicit none
    type(meshdef) :: mesh
    integer(kint), intent(in) :: elem_id(:)
    real(kdouble), intent(inout) :: e_mat(:,:)
    real(kdouble), allocatable :: B_mat(:,:)
    real(kdouble) :: det, G
    real(kdouble) :: x(3), y(3)
    integer(kint) :: i, j

    x = 0.0d0
    y = 0.0d0
    det = 0.0d0
    allocate(B_mat(2,3) ,source=0.0d0)

    do i=1,3
        x(i) = mesh%node(1, elem_id(i))
        y(i) = mesh%node(2, elem_id(i))
    enddo

    det = (x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1))
    G = abs(det/2.0d0)

    B_mat(1,1) = (y(2) - y(3))/det
    B_mat(1,2) = (y(3) - y(1))/det
    B_mat(1,3) = (y(1) - y(2))/det
    B_mat(2,1) = (x(3) - x(2))/det
    B_mat(2,2) = (x(1) - x(3))/det
    B_mat(2,3) = (x(2) - x(1))/det

    do i=1,3
        do j=1,3
            e_mat(j,i) = G*(B_mat(1,i)*B_mat(1,j)+B_mat(2,i)*B_mat(2,j))
        enddo
    enddo
end subroutine

subroutine add_global_mat(elem_id, e_mat, A)
    integer(kint), intent(in) :: elem_id(:)
    real(kdouble), intent(in) :: e_mat(:,:)
    real(kdouble), intent(inout) :: A(:,:)
    integer(kint) :: i, j, ii, jj

    do j=1, 3
        do i=1, 3
            ii = elem_id(i)
            jj = elem_id(j)
            A(ii,jj) = A(ii,jj) + e_mat(i,j)
        enddo
    enddo
end subroutine add_global_mat

subroutine output_vtk(Nm, elements, x, nn, ne)
    implicit none
    integer(4) :: i
    integer(4) :: nn, ne
    real(8), allocatable :: Nm(:,:)
    real(8), allocatable :: x(:)
    integer(4), allocatable :: elements(:,:)

    open(20,file="../visual/laplace_solution.vtk",status="replace")
        write(20,'(a)')'# vtk DataFile Version 4.1'
        write(20,'(a)')'Laplace Equation Solution'
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
    close(20)
end subroutine output_vtk
    
end program main