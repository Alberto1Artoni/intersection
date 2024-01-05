subroutine READ_OPTIONS()

    use globalVariables, only : mpi_id, mpi_ierr

    implicit none

    !! options will be read only on core one
    ! then broadcast them all

    if (mpi_id .eq. 0) then

        write(*,*) "Currently no options have been initialized"

    endif

    ! todo migliorare design: qui andrebbe un file con le opzioni griglia
    ! call READ_GRID_OPTIONS()
!   reader = "openfoam"
!   parallelOF = .false.
!   pathOF = "/home/alberto/phd/coding/fortran/intersection/test/readFluidGrid/"

end subroutine
