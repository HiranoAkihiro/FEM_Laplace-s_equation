module mod_matrix
    use mod_utils
    implicit none
    
contains
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
end module mod_matrix