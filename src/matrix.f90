module mod_matrix
    use mod_utils
    implicit none
    
contains
    subroutine get_coef_mat(mesh, param, mat)
        implicit none
        type(meshdef) :: mesh
        type(paramdef) :: param
        real(kdouble), intent(inout) :: mat(:,:)
        integer(kint) :: i
        real(kdouble), allocatable :: e_mat(:,:)
        integer(kint), allocatable :: elem_id(:)

        if(param%elastic)then
            allocate(e_mat(mesh%nbase_func*ndof, mesh%nbase_func*ndof), source=0.0d0)
        else
            allocate(e_mat(mesh%nbase_func, mesh%nbase_func), source=0.0d0)
        endif
        allocate(elem_id(mesh%nbase_func), source=0)

        do i=1, mesh%nelem
            call get_elem_id(mesh, elem_id, i)
            if(param%elastic)then
                call get_element_stiff_mat(mesh, param, elem_id, e_mat)
                call add_global_stiff_mat(elem_id, e_mat, mat)
            else
                call get_element_mat(mesh, elem_id, e_mat)
                call add_global_mat(elem_id, e_mat, mat)
            endif
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

    subroutine get_element_stiff_mat(mesh, param, elem_id, e_mat)
        implicit none
        type(meshdef) :: mesh
        type(paramdef) :: param
        integer(kint), intent(in) :: elem_id(:)
        real(kdouble), intent(inout) :: e_mat(:,:)
        real(kdouble), allocatable :: B_mat(:,:)
        real(kdouble), allocatable :: D_mat(:,:)
        real(kdouble), allocatable :: DB_mat(:,:)
        real(kdouble) :: det, G
        real(kdouble) :: x(3), y(3)
        integer(kint) :: i, j

        x = 0.0d0
        y = 0.0d0
        det = 0.0d0
        allocate(B_mat(3,6) ,source=0.0d0)
        allocate(D_mat(3,3) , source=0.0d0)
        allocate(DB_mat(3,6) , source=0.0d0)

        do i=1,3
            x(i) = mesh%node(1, elem_id(i))
            y(i) = mesh%node(2, elem_id(i))
        enddo

        det = (x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1))
        G = abs(det/2.0d0)

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

        D_mat(1,1) = param%E/(1.0d0-param%nu**2.0d0)
        D_mat(1,2) = param%E*param%nu/(1.0d0-param%nu**2.0d0)
        D_mat(2,1) = param%E*param%nu/(1.0d0-param%nu**2.0d0)
        D_mat(2,2) = param%E/(1.0d0-param%nu**2.0d0)
        D_mat(3,3) = param%E/(2.0d0*(1.0d0+param%nu))

        DB_mat = matmul(D_mat, B_mat)
        e_mat = matmul(transpose(B_mat), DB_mat)
        e_mat = e_mat*G
    end subroutine get_element_stiff_mat

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

    subroutine add_global_stiff_mat(elem_id, e_mat, A)
        integer(kint), intent(in) :: elem_id(:)
        real(kdouble), intent(in) :: e_mat(:,:)
        real(kdouble), intent(inout) :: A(:,:)
        integer(kint) :: i, j, k, l, ii, jj

        do j=1, 3
            do i=1, 3
                ii = elem_id(i)
                jj = elem_id(j)
                do k=1, ndof
                    do l=1, ndof
                        A(ndof*(ii-1)+k,ndof*(jj-1)+l) = A(ndof*(ii-1)+k,ndof*(jj-1)+l) + e_mat(ndof*(i-1)+k,ndof*(j-1)+l)
                    enddo
                enddo
            enddo
        enddo
    end subroutine add_global_stiff_mat

    subroutine get_RHS(mesh, b)
        implicit none
        type(meshdef) :: mesh
        real(kdouble), intent(inout) :: b(:)
        integer(kint) :: icel, in, i
        real(kdouble) :: S
        integer(kint), allocatable :: elem_id(:)

        allocate(elem_id(mesh%nbase_func), source=0)
        S = (0.5d0**2.0d0)/dble(mesh%nelem)
        do icel=1,mesh%nelem
            call get_elem_id(mesh, elem_id, icel)
            do i=1,mesh%nbase_func
                in = elem_id(i)
                b(in) = b(in) + S/3.0d0
            enddo
        enddo
    end subroutine get_RHS

    subroutine get_load(mesh, b)
        implicit none
        type(meshdef) :: mesh
        real(kdouble), intent(inout) :: b(:)
        integer(kint) :: i, in, dof

        do i=1, mesh%nload
            in = mesh%load_id(i)
            dof = mesh%load_dof(i)
            b((in-1)*ndof+dof) = mesh%load(i)
        enddo
    end subroutine get_load
end module mod_matrix