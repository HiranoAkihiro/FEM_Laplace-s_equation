module mod_valid
    use mod_utils
    implicit none
    
contains

subroutine get_element_fluxes(coords, elements, u, fluxes)
    real(kdouble), intent(in) :: coords(:,:)  ! 形状: (2, num_nodes)
    integer(kint), intent(in) :: elements(:,:)  ! 形状: (3, num_elements)
    real(kdouble), intent(in) :: u(:)  ! 形状: (num_nodes)
    real(kdouble), intent(out) :: fluxes(:,:)  ! 形状: (num_elements, 2, 3)

    integer(kint) :: num_elements, i
    real(kdouble) :: x(3), y(3), u_values(3), area, b(3), c(3)
    real(kdouble) :: grad_x, grad_y

    num_elements = size(elements, 2)

    do i = 1, num_elements
        ! 要素の節点座標を取得
        x = coords(1, elements(:,i))
        y = coords(2, elements(:,i))
        
        ! 要素の節点での解の値を取得
        u_values = u(elements(:,i))

        ! 要素の面積を計算
        area = 0.5 * abs((x(2) - x(1)) * (y(3) - y(1)) - (x(3) - x(1)) * (y(2) - y(1)))

        ! b と c の係数を計算
        b = [y(2) - y(3), y(3) - y(1), y(1) - y(2)]
        c = [x(3) - x(2), x(1) - x(3), x(2) - x(1)]

        ! 勾配を計算
        grad_x = sum(b * u_values) / (2.0 * area)
        grad_y = sum(c * u_values) / (2.0 * area)

        ! フラックスを計算して格納
        fluxes(i, 1) = -grad_x  ! x方向のフラックス
        fluxes(i, 2) = -grad_y  ! y方向のフラックス
    end do
end subroutine get_element_fluxes
end module mod_valid