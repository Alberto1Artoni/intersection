program MAIN

    use globalVariables, only : mpi_id, mpi_ierr
    use fast_speed_polyhedron ! intersection map and polyFluid

    use MOD_OPENFOAM_GRID
    use MOD_INTERSECTION

    implicit none

    type(mesh) :: acoustic, fluid

    ! sarebbe più elegante mettere tutto insieme
    type(polyhedron_fast) , allocatable, dimension(:) :: acousticPoly, fluidPoly
    

    ! MPI INITIALIZATION
    call INITIALIZE_MPI()

    if (mpi_id .eq. 0) then
        write(*,*) "---------------------------------------------"
        write(*,*) "Wrapper to CGAL to perform mesh intersections"
        write(*,*) "---------------------------------------------"
    endif

    ! READ OPTIONS
    ! this has to be built
    ! system/controlDict

    call READ_OPTIONS()
    ! READ ACOUSTIC GRID
    !! => new format
    !! constant/polyMesh => mesh tipo OpenFOAM
    !! constant/old.mesh => mesh tipo SPEED

    acoustic%path = adjustl(trim("/home/alberto/phd/coding/fortran/intersection/test/OpenFOAM/Cube/"))
    acoustic%reader = "openfoam"
    acoustic%nameGrid = "acoustic"
    call READ_POLY_GRID(acoustic)
    allocate(acousticPoly(acoustic%nb_elems))     
    call MAKE_MESH_POLY(acoustic, acousticPoly, acoustic%nb_elems)


    fluid%path = adjustl(trim("/home/alberto/phd/coding/fortran/intersection/test/OpenFOAM/Voronoi2D/"))
    fluid%reader = "openfoam"
    fluid%nameGrid = "fluid"
    call READ_POLY_GRID(fluid)
    allocate(fluidPoly(fluid%nb_elems))     
    call MAKE_MESH_POLY(fluid, fluidPoly, fluid%nb_elems)


    call MPI_BARRIER( mpi_ierr )

    !@todo
    ! Intenzione: riciclare un lettore mio e metto griglia poliedrica formato OpenFOAM
    ! In un secondo momento, includo anche la griglia tipo SPEED

    ! PARTITION ACOUSTIC GRID

    ! partizionata la griglia devo ragionare element-wise. Questo non era fatto lato SPEED
    ! perché avevo le mappe e ricostruivo al momento gli elementi 

    ! READ FLUID GRID 

    ! PERFORM INTERSECTIONS
    call NAIVE_INTERSECTIONS(acousticPoly, acoustic%nb_elems, fluidPoly, fluid%nb_elems)
    

    ! DEALLOCATING

    ! RETURN
    call MPI_FINALIZE(mpi_ierr)

end program
