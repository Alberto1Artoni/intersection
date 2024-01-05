subroutine READ_GRID_OPTIONS( grid, filename)

    use MOD_OPENFOAM_GRID

    implicit none

    type(mesh) :: grid
    character*256 :: filename, inline, gridname, tmp
    character*13  :: keyword

    integer*4 :: status
    integer*4 :: ileft,iright, order, stages
    integer*4 :: arglen,i
    integer*4 :: file_row = 0

    gridname = adjustl(trim("./system/"//adjustl(trim(filename))))
    open(20,file=gridname)

    do 
        read(20, '(A)',IOSTAT = status) inline
        file_row = file_row + 1
        if (status.ne.0)  then
            exit
        endif

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
            case('path')
                read(inline(ileft:iright), '(A)') grid%path
            case('reader')
                read(inline(ileft:iright), '(A)') grid%reader
            case('nameGrid')
                read(inline(ileft:iright), '(A)') grid%nameGrid
        end select
    end do

    close(20)

end subroutine
