subroutine READ_OPTIONS( acoustic, fluid)

    use globalVariables, only : mpi_id, mpi_ierr
    use MOD_OPENFOAM_GRID

    implicit none

    type(mesh) :: acoustic, fluid
    character*256 :: inline, acousticDict, fluidDict
    character*13  :: keyword

    integer*4 :: status
    integer*4 :: ileft,iright, order, stages
    integer*4 :: arglen,i
    integer*4 :: file_row = 0


    !! options currently are read from all the processors
    open(40,file="./system/controlDict")

    do 
        read(40, '(A)',IOSTAT = status) inline
        file_row = file_row + 1
        if (status.ne.0) exit

        !> Skip comment lines in the aeroOptions.input file
        if (inline(1:1) .eq. ' ') then
           cycle
        endif

        ileft = 1
        iright = len_trim(inline)

        !> Compute index to first non-keyword argument
        do i = 1,iright
            if (inline(i:i) .eq. ' ') exit
        enddo

        ileft = i + 1

        !> Get keyword out of whole line by cutting after the first ' ' afterwards
        keyword = inline(1:(ileft-2))

        !> Now ileft points after the first blank
        arglen = len_trim(inline(ileft:iright))

        !> Remove this for keywords without arguments!
        if (arglen .eq. 0) then
           write(*,'(A,I3,A,A,A)') 'FATAL in aeroOptions.input, row',file_row, &
                                   ': no argument given for command "', trim(keyword), '"'
           call EXIT()
        endif

        select case (keyword)
            case('meshAcoustic')
                read(inline(ileft:iright), '(A)') acousticDict
                call READ_GRID_OPTIONS(acoustic, acousticDict)
            case('meshFluid')
                read(inline(ileft:iright), '(A)') fluidDict
                call READ_GRID_OPTIONS(fluid, fluidDict)
        end select
    end do

    close(40)

end subroutine
