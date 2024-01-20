module MOD_INTERSECTION

    contains

    subroutine MAKE_PREPROCESSING_INTERSECTIONS( )
        ! this routine should accelerate intersection computations
        ! returns the maps of the intersecting elements

        use fast_speed_polyhedron
        use intersectionCGAL_new

        implicit none

    ! needs to be implemented


    end subroutine MAKE_PREPROCESSING_INTERSECTIONS

    subroutine NAIVE_INTERSECTIONS( acuPoly, nb_elemAcu, fluPoly, nb_elemFlu)

        use fast_speed_polyhedron
        use intersectionCGAL_new
        use globalVariables, only : mpi_id, mpi_ierr, mpi_np

        implicit none
include 'SPEED.MPI'

        integer*4 :: nb_elemAcu, nb_elemFlu
        type(polyhedron_fast), dimension(nb_elemAcu) :: acuPoly
        type(polyhedron_fast), dimension(nb_elemFlu) :: fluPoly
        type(polyhedron_fast) :: outIntersection

        integer*4 :: jAcu, jFlu, start, end, temp
        real*8 :: xAmin, xAmax, xFmin, xFmax, &
                  yAmin, yAmax, yFmax, yFmin, &
                  zAmin, zAmax, zFmin, zFmax, toll

        real*8 :: volume, vol

        logical   :: isIntersecting
        character*128  :: filename
        character*5000 :: outstring

        ! partitioning where everyone knows everything
        if ( mpi_id .lt. mod(nb_elemAcu, mpi_np)) then
            start = (nb_elemAcu / mpi_np)*  mpi_id    + mpi_id + 1
            end   = (nb_elemAcu / mpi_np)* (mpi_id+1) + mpi_id + 1
        else
            start = (nb_elemAcu / mpi_np)*  mpi_id    + mod(nb_elemAcu, mpi_np) + 1
            end   = (nb_elemAcu / mpi_np)* (mpi_id+1) + mod(nb_elemAcu, mpi_np)
        end if

        temp = 0
        volume = 0.0

        do jAcu= start, end

            call GET_POLY_BBOX(xAMin,yAMin,zAMin, xAMax, yAMax, zAMax, acuPoly(jAcu))

            toll = (xAMax - xAMin) / 1000.0;

            if (mod( jAcu - start, 25 ).eq.0) then
                write(*,'(A1,I0,A2,A24,F6.2,A1)') "[",mpi_id,"] ","Intersection progress   ", REAL(jAcu - start) / (end-start)*100, "%"
            endif

            do jFlu=1, nb_elemFlu

                call GET_POLY_BBOX(xFMin,yFMin,zFMin, xFMax, yFMax, zFMax, fluPoly(jFlu))

                if ( (.not.(( xAmin - xFmax .gt. toll) .or. (xFmin - xAmax .gt. toll) ) ) .and. &
                     (.not.(( yAmin - yFmax .gt. toll) .or. (yFmin - yAmax .gt. toll) ) ) .and. &
                     (.not.(( zAmin - zFmax .gt. toll) .or. (zFmin - zAmax .gt. toll) ) ) ) then


                    isIntersecting = .false.
                
                    call intersectPolyhedra(acuPoly(jAcu), fluPoly(jFlu), &
                                           outIntersection, isIntersecting)

                    if (isIntersecting) then
                        write(filename,'(A,I0,A,I0,A)') './DUMP/acu', jAcu ,"_", jFlu, ".vtk"
                        temp = temp + 1
                        call EXPORT_VTK_POLYDATA_DATA_FAST(filename, outIntersection, jAcu, temp)

                        call GET_POLY_VOLUME(outIntersection, vol)
                        volume = volume + vol

                        call POLY_DESTROY_FAST(outIntersection)
                    endif

                endif

            end do


        end do

        write(*,'(A1,I0,A2,A27)') "[",mpi_id,"] ","Intersection progress: done!"

        call MPI_ALLREDUCE(MPI_IN_PLACE, volume, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, mpi_ierr)

        write(*,*) "Volume: ", volume

    end subroutine


end module MOD_INTERSECTION
