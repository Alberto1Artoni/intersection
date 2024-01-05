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

    call READ_OPTIONS()
!   ! READ ACOUSTIC GRID
!   !! => new format
!   !! constant/polyMesh => mesh tipo OpenFOAM
!   !! constant/old.mesh => mesh tipo SPEED

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
    write(*,*) fluid%nb_elems
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
    call NAIVE_INTERSECTIONS(acousticPoly, acoustic%nb_elems, fluidPoly, fluid%nb_elems)

!   ! intersect only two elements
!   filename = "./Solid/outCubo1.off"
!   call READ_OFF_FILE_DIME_FAST(filename, polyAcu%nV, polyAcu%nF, edge)
!   allocate(polyAcu%vX(polyAcu%nV), polyAcu%vY(polyAcu%nV), polyAcu%vZ(polyAcu%nV))
!   call READ_OFF_FILE_FAST(filename, polyAcu)

!   filename = "./Solid/outCubo2.off"
!   call READ_OFF_FILE_DIME_FAST(filename, polyFlu%nV, polyFlu%nF, edge)
!   allocate(polyFlu%vX(polyFlu%nV), polyFlu%vY(polyFlu%nV), polyFlu%vZ(polyFlu%nV))
!   call READ_OFF_FILE_FAST(filename, polyFlu)

!   write(*,*) "call to intersect"
!   call intersectPolyhedra(  polyAcu, polyFlu, inter, isIntersecting)

    ! DEALLOCATING

    ! RETURN
    call MPI_FINALIZE(mpi_ierr)

end program
