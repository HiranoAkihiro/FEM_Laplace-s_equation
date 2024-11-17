module mod_utils
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
contains
    subroutine apply_bc(mesh, mat)
        implicit none
        type(meshdef) :: mesh
        type(matdef) :: mat
        integer(kint) :: i, j, id

        ! do j=1, mesh%nbc
        !     id = mesh%bc_id(j)
        !     do i=1, mesh%nnode
        !         mat%A(i,id) = 0.0d0
        !     enddo

        !     do i=1, mesh%nnode
        !         mat%A(id,i) = 0.0d0
        !     enddo
        !     mat%A(id,id) = 1.0d0
        !     mat%b(id) = mesh%bc(j)
        ! enddo

        do j=1, mesh%nbc
            id = mesh%bc_id(j)
            do i=1, mesh%nnode
                if (i /= id) then
                    mat%b(i) = mat%b(i) - mat%A(i,id) * mesh%bc(j)
                    mat%A(i,id) = 0.0d0
                end if
            enddo
            mat%A(id,:) = 0.0d0
            mat%A(id,id) = 1.0d0
            mat%b(id) = mat%b(id) + mesh%bc(j)
        enddo
    end subroutine

    subroutine solver(mat)
        type(matdef) :: mat
        real(kdouble), allocatable :: A(:,:)
        real(kdouble), allocatable :: b(:)
        real(kdouble), allocatable :: x(:)
        integer(kint) :: n

        n = size(mat%A, dim=1)
        allocate(A(n,n), source=0.0d0)
        allocate(b(n), source=0.0d0)
        allocate(x(n), source=0.0d0)

        A = mat%A
        b = mat%b
        ! call LU_decomposition(n, A, b, x)
        call gauss_elimination(A, b, x, n)
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

    subroutine gauss_elimination(a, b, x, n)
        implicit none
        integer, intent(in) :: n
        real(8), intent(inout) :: a(n,n), b(n)
        real(8), intent(out) :: x(n)
        real(8) :: factor
        integer :: i, j, k
    
        ! 前進消去
        do k = 1, n-1
            do i = k+1, n
                factor = a(i,k) / a(k,k)
                a(i,k+1:n) = a(i,k+1:n) - factor * a(k,k+1:n)
                b(i) = b(i) - factor * b(k)
            end do
        end do
    
        ! 後退代入
        x(n) = b(n) / a(n,n)
        do i = n-1, 1, -1
            x(i) = (b(i) - sum(a(i,i+1:n) * x(i+1:n))) / a(i,i)
        end do
    end subroutine gauss_elimination

    subroutine debug(char, flag)
        implicit none
        character(len=*), intent(in) :: char
        character(len=:), allocatable :: name
        logical, intent(in) :: flag
        
        if(flag)then
            name = '***** done '//trim(char)//' *****'
            write(*,*)name
            write(*,*)
        endif
    end subroutine debug
end module mod_utils