subroutine READ_MESH_POLY_OPENFOAM(grid)

    use MOD_OPENFOAM_GRID
    use globalVariables, only : mpi_id, mpi_ierr

    implicit none

    type(mesh) :: grid
    character*256 :: outGrid

    if (mpi_id.eq.0) then
       write(*,*) " "
       write(*,*) "Mesh: assuming .foam format"
    endif

    ! maybe I should make those global?
    call READ_DIME_MESH_OPENFOAM( grid )

    allocate( grid%faceVector( grid%nb_faces ))
    call READ_DIME_MESH_FACES( grid )

    allocate( grid%faceConnVector(0: grid%nb_faces_conn))
    call READ_DIME_MESH_CONN_FACES( grid )

    allocate(grid%elemFacesDime(grid%nb_elems))
    call READ_DIME_MESH_ELEMENTS(grid)

    ! todo
    allocate(grid%elemFacesConn(0:grid%nb_elems_conn))
    allocate(grid%connElem(0:grid%nb_elems_conn))
    call READ_MESH_CONN_ELEMENTS( grid )


    allocate(grid%points(grid%nb_points,3))
    call READ_MESH_POINTS( grid )

    outGrid = "/home/alberto/phd/coding/fortran/intersection/test/Intersection/VTK/"//trim(adjustl(trim(grid%nameGrid)))//".vtk"
    call DUMP_MESH_VTK( grid, outGrid)


end subroutine
