module mod_io
    use mod_utils
    implicit none
    
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

        read(22,*)mesh%nbc, n
        allocate(mesh%bc(mesh%nbc), source=0.0d0)
        allocate(mesh%bc_id(mesh%nbc))
        do i=1, mesh%nbc
            read(22,*)mesh%bc_id(i), n, mesh%bc(i)
        enddo
    close(20)
    close(21)
    close(22)
end subroutine input_mesh
    
subroutine output_vtk(Nm, elements, x, nn, ne, fluxes)
    implicit none
    integer(4) :: i
    integer(4) :: nn, ne
    real(8), allocatable :: Nm(:,:)
    real(8), allocatable :: x(:)
    integer(4), allocatable :: elements(:,:)
    real(8), allocatable :: fluxes(:,:)

    open(20,file="../visual/laplace_solution_with_fluxes.vtk",status="replace")
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
end module mod_io