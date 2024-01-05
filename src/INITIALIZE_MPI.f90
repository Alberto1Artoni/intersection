subroutine INITIALIZE_MPI()
   
    use globalVariables, only : mpi_np, mpi_comm, mpi_id, mpi_ierr

    implicit none 


    ! Initialization of MPI ranks
    call MPI_INIT(mpi_ierr)
    call MPI_COMM_RANK(mpi_comm, mpi_id, mpi_ierr)
    call MPI_COMM_SIZE(mpi_comm, mpi_np, mpi_ierr)

    ! There is also the space for the initialization of Petsc:
    ! not done here atm

end subroutine
