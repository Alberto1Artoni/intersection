program MAIN

    use globalVariables, only : mpi_id, mpi_ierr
    use fast_speed_polyhedron ! intersection map and polyFluid
    use intersectionCGAL_new

    use MOD_OPENFOAM_GRID
    use MOD_INTERSECTION

    implicit none

    type(mesh) :: acoustic, fluid
    character*100 :: filename
    type(polyhedron_fast) :: polyAcu, polyFlu, inter

    ! sarebbe più elegante mettere tutto insieme
    type(polyhedron_fast), allocatable, dimension(:) :: acousticPoly
    type(polyhedron_fast), allocatable, dimension(:) :: fluidPoly
    integer*4 :: edge
    logical :: isIntersecting

    edge = 0

    call INITIALIZE_MPI()

    if (mpi_id .eq. 0) then
        write(*,*) "---------------------------------------------"
        write(*,*) "Wrapper to CGAL to perform mesh intersections"
        write(*,*) "---------------------------------------------"
    endif

    ! READ OPTIONS
    ! this has to be built
    ! system/controlDict

    call READ_OPTIONS( acoustic, fluid)

    if (mpi_id .eq. 0) then
        write(*,*) ""
        write(*,*) "Options read!"
        write(*,*) ""
    endif

    call READ_POLY_GRID(acoustic)
    allocate(acousticPoly(acoustic%nb_elems))     
    call MAKE_MESH_POLY(acoustic, acousticPoly, acoustic%nb_elems)


    call READ_POLY_GRID(fluid)
    allocate(fluidPoly(fluid%nb_elems))     
    call MAKE_MESH_POLY(fluid, fluidPoly, fluid%nb_elems)


    !@todo
    ! Intenzione: riciclare un lettore mio e metto griglia poliedrica formato OpenFOAM
    ! In un secondo momento, includo anche la griglia tipo SPEED

    ! PARTITION ACOUSTIC GRID

    ! partizionata la griglia devo ragionare element-wise. Questo non era fatto lato SPEED
    ! perché avevo le mappe e ricostruivo al momento gli elementi 

    ! READ FLUID GRID 

    ! PERFORM INTERSECTIONS
    if (mpi_id .eq. 0) then
        write(*,*) ""
        write(*,*) "Starting intersections..."
        write(*,*) ""
    endif
    call NAIVE_INTERSECTIONS(acousticPoly, acoustic%nb_elems, fluidPoly, fluid%nb_elems)


    ! DEALLOCATING

    ! RETURN
    call MPI_FINALIZE(mpi_ierr)

end program
