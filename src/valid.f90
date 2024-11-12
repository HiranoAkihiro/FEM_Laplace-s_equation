module mod_valid
    use mod_utils
    use mod_matrix
    implicit none
    
contains

subroutine get_elem_fluxes(mesh, u, fluxes)
    type(meshdef) :: mesh
    real(kdouble), intent(in) :: u(:)  ! 形状: (num_nodes)
    integer(kint), allocatable :: elem_id(:)
    real(kdouble), intent(out) :: fluxes(:,:)  ! 形状: (num_elements, 2, 3)
    integer(kint) :: icel, j
    real(kdouble) :: B(2,3)
    real(kdouble) :: u_elem(3), x(3), y(3)
    real(kdouble) :: S

    allocate(elem_id(mesh%nbase_func))
    fluxes = 0.0d0

    do icel = 1, mesh%nelem
        call get_elem_id(mesh, elem_id, icel)
        u_elem(:) = u(elem_id(:))

        x(:) = mesh%node(1,elem_id(:))
        y(:) = mesh%node(2,elem_id(:))

        S = 0.5 * abs((x(2) - x(1)) * (y(3) - y(1)) - (x(3) - x(1)) * (y(2) - y(1)))

        B(1,1) = (y(2) - y(3))/(2.0d0*S)
        B(1,2) = (y(3) - y(1))/(2.0d0*S)
        B(1,3) = (y(1) - y(2))/(2.0d0*S)
        B(2,1) = (x(3) - x(2))/(2.0d0*S)
        B(2,2) = (x(1) - x(3))/(2.0d0*S)
        B(2,3) = (x(2) - x(1))/(2.0d0*S)

        do j=1,3
            fluxes(icel,1) = fluxes(icel,1) + B(1,j)*u_elem(j)
            fluxes(icel,2) = fluxes(icel,2) + B(2,j)*u_elem(j)
        enddo
    enddo
    fluxes = abs(fluxes)
end subroutine get_elem_fluxes
end module mod_valid