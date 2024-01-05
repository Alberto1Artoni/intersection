module MOD_OPENFOAM_GRID
    ! module with OpenFOAM file format grid

    implicit none

        type :: mesh
            ! options
            ! potrebbe essere interessante salvare il percorso file della griglia qui
            
            character*256 :: path
            character*256 :: reader
            character*256 :: nameGrid

            ! grid 
            integer*4 :: nb_elems, nb_points, nb_faces
            integer*4 :: nb_internal_faces
            integer*4 :: nb_faces_conn
            integer*4 :: nb_elems_conn
            integer*4, dimension(:), allocatable :: faceVector        ! cambiare nome?
            integer*4, dimension(:), allocatable :: faceConnVector
            integer*4, dimension(:), allocatable :: elemFacesDime
            integer*4, dimension(:), allocatable :: elemFacesConn
            integer*4, dimension(:), allocatable :: faceOwner, faceNeig       ! these are needed for higher order
            integer*4, dimension(:), allocatable :: connElem
            real*8, dimension(:,:), allocatable  :: points
        end type mesh


    contains

    subroutine READ_DIME_MESH_OPENFOAM(grid)

        implicit none
      
        type(mesh) :: grid

        character*256 :: gridfile
        character*256 :: inline
        integer*4 :: jj,kk

        gridfile=trim(adjustl(trim(grid%path)//'constant/polyMesh/points'))

        !> Open mesh file
        open(40,file=gridfile)

        do jj=1,18
          read(40,'(A)') inline
        enddo
        read(40,'(A)') inline
        read(inline,*) grid%nb_points

        close(40)

        gridfile=trim(adjustl(trim(grid%path)//'constant/polyMesh/owner'))
        open(40,file=gridfile)

        do jj=1,12
          read(40,'(A)') inline
        enddo

        read(40,'(A)') inline
        jj = 1;
        do while (inline(jj:jj+6) .ne. 'nCells:')
        jj = jj+1;
        end do

        kk = jj+7;

        do while (inline(jj:jj+1) .ne. ' ')
        jj = jj+1;
        end do

        read(inline(kk:jj),*) grid%nb_elems

        close(40)

        gridfile=trim(adjustl(trim(grid%path)//'constant/polyMesh/faces'))
        open(40,file=gridfile)

        do jj=1,18
          read(40,'(A)') inline
        enddo

        read(40,'(A)') inline
        read(inline,*) grid%nb_faces

        close(40)


        gridfile=trim(adjustl(trim(grid%path)//'constant/polyMesh/neighbour'))
        open(40,file=gridfile)

        do jj=1,19
          read(40,'(A)') inline
        enddo
        read(40,'(A)') inline
        read(inline,*) grid%nb_internal_faces

        close(40)

    end subroutine READ_DIME_MESH_OPENFOAM

    subroutine READ_DIME_MESH_FACES(grid)

        implicit none
        type(mesh) :: grid

        character*256 :: gridfile, inline
        integer*4 :: jj,kk
        integer*4 :: ierror, intVal, val             ! this checks if the input is valid


        grid%nb_faces_conn = grid%nb_faces + grid%nb_faces

        gridfile=trim(adjustl(trim(grid%path)//'constant/polyMesh/faces'))
        open(40,file=gridfile)

        do jj=1,20
            read(40,'(A)') inline
        enddo
        ierror = 1
        do jj=1,grid%nb_faces

            read(40,'(A)') inline
            ! devo prendere la lunghezza fino alla parentesi aperta
            ! potrei avere due cifre
            val= INDEX(inline, "(")-1
            read(inline(1:val),*, iostat=ierror) intVal

            if ( ierror .ne. 0 ) then

                read(40,'(A)') inline
                read(inline(1:2),*) intVal
                ! write(*,*) inline, intVal

                grid%faceVector(jj) = intVal
                do kk = 1,(intVal+3)
                    read(40,'(A)') inline
                end do
                grid%nb_faces_conn = grid%nb_faces_conn + grid%faceVector(jj)

            else

                grid%faceVector(jj) = intVal
                grid%nb_faces_conn = grid%nb_faces_conn + grid%faceVector(jj)
            end if

        end do

        close(40)

    end subroutine READ_DIME_MESH_FACES

    subroutine READ_DIME_MESH_CONN_FACES(grid)

        implicit none

        type(mesh) :: grid
        character*256 :: gridfile, inline

        integer*4 :: jj,kk, iter
        integer*4 :: ierror, intVal, val             ! this checks if the input is valid
        integer*4 :: valIn, valOut, valTmp, valOld   ! iteratori sulla stringa

        grid%faceConnVector(0) = grid%nb_faces
        grid%faceConnVector(1) = grid%nb_faces + 1

        do jj=2,grid%nb_faces
            grid%faceConnVector(jj) = grid%faceConnVector(jj-1) + grid%faceVector(jj-1) + 1
            ! +1 needed for the face size
        end do

        iter = grid%nb_faces+1


        gridfile=trim(adjustl(trim(grid%path)//'constant/polyMesh/faces'))
        open(40,file=gridfile)

        do jj=1,20
            read(40,'(A)') inline
        enddo


        do jj=1,grid%nb_faces

       ! am also storing the number of caes in this vector
       !his is a redundant information, but avoids computations
            grid%faceConnVector(iter) =  grid%faceVector(jj)
            iter = iter+1

            read(40,'(A)') inline
            ! devo prendere la lunghezza fino alla parentesi aperta
            ! potrei avere due cifre
            valIn = INDEX(inline, "(")+1
            valOut= INDEX(inline, ")")-1

            if ( valOut .eq. -1 ) then

                read(40,'(A)') inline     ! read face dimension if gt then 11
                read(40,'(A)') inline     ! skip parenthesis

                do kk = 1,grid%faceVector(jj)
                    read(40,'(A)') inline
                    read(inline,*) intVal
                    grid%faceConnVector(iter) =  intVal
                    iter = iter+1
                end do
                read(40,'(A)') inline     ! skip parenthesis
                read(40,'(A)') inline     ! skip blank space

            else
                valOld = valIn
                do kk =1,grid%faceVector(jj)-1
                    valTmp = INDEX(inline(valOld:valOut), " ")

                    read(inline(valOld:valTmp+valOld),*) intVal
                    grid%faceConnVector(iter) =  intVal
                    iter = iter+1
                    valOld = valTmp+valOld
                end do

                read(inline(valOld:valOut),*) intVal

                grid%faceConnVector(iter) =  intVal
                iter = iter+1

            end if

        end do

        close(40)

    end subroutine READ_DIME_MESH_CONN_FACES


    subroutine READ_DIME_MESH_ELEMENTS(grid)

        implicit none

        type(mesh) :: grid

        character*256 :: gridfile, inline
        integer*4 :: jj,kk
        integer*4       :: ierror, intVal, val, neigVal

        grid%elemFacesDime = 0
        grid%nb_elems_conn = 0

        gridfile=trim(adjustl(trim(grid%path)//'constant/polyMesh/owner'))
        open(40,file=gridfile)

        do jj=1,21
            read(40,'(A)') inline
        enddo

        ! read the file
        ! access the element and sum one to the face dimensions

        do jj=1,grid%nb_faces

            read(40,'(A)') inline
            read(inline,*) intVal
            intVal = intVal + 1       ! from C++ to Fortran
            grid%elemFacesDime(intVal) = grid%elemFacesDime(intVal)+1  
        end do

        close(40)

        ! now read the neighbour elements
        gridfile=trim(adjustl(trim(grid%path)//'constant/polyMesh/neighbour'))
        open(40,file=gridfile)

        do jj=1,19
            read(40,'(A)') inline
        enddo

        read(40,'(A)') inline
        read(inline,*) neigVal
        read(40,'(A)') inline     ! skip parenthesis

        ! read the file
        ! access the element and sum one to the face dimensions

        do jj=1,neigVal

            read(40,'(A)') inline
            read(inline,*) intVal
            intVal = intVal + 1       ! from C++ to Fortran
            grid%elemFacesDime(intVal) = grid%elemFacesDime(intVal)+1  

        end do

        close(40)

        grid%nb_elems_conn = grid%nb_elems + grid%nb_elems;
        do jj=1,grid%nb_elems
            grid%nb_elems_conn = grid%elemFacesDime(jj) + grid%nb_elems_conn;
        end do

    end subroutine READ_DIME_MESH_ELEMENTS


    subroutine READ_MESH_CONN_ELEMENTS( grid )

        implicit none

        type(mesh) :: grid

        character*256 :: gridfile, inline
        integer*4 :: jj,kk
        integer*4, dimension(:), allocatable :: elemIter          ! saves the current shift

        integer*4 :: ierror, intVal, val
        integer*4 :: nF, j, jface, neigVal

        grid%elemFacesConn(0) = grid%nb_elems
        grid%elemFacesConn(1) = grid%nb_elems + 1
        allocate(elemIter(grid%nb_elems))
        elemIter = 0

        do jj=2,grid%nb_elems
            grid%elemFacesConn(jj) = grid%elemFacesConn(jj-1) + grid%elemFacesDime(jj-1) + 1
                             ! get offset from previous element index
                             ! add face dimension of the previous face
                             ! add one to the offset for the face dimension
        end do

        do jj=1,grid%nb_elems
            grid%elemFacesConn(grid%elemFacesConn(jj)) = grid%elemFacesDime(jj)
        end do

        gridfile=trim(adjustl(trim(grid%path)//'constant/polyMesh/owner'))
        open(40,file=gridfile)

        do jj=1,21
            read(40,'(A)') inline
        enddo

        ! read the file

        allocate(grid%faceOwner(grid%nb_faces))

        do jj=1,grid%nb_faces
            read(40,'(A)') inline
            read(inline,*) intVal
            intVal = intVal + 1       ! from C++ to Fortran
            elemIter(intVal) = elemIter(intVal) +1                       ! shift for the conn vector
            grid%elemFacesConn(grid%elemFacesConn(intVal) + elemIter(intVal) ) = jj  ! already Fortran
            grid%faceOwner(jj) = intVal
            ! shift to current element, skip the face number, 
        end do

        close(40)

        ! now read the neighbour elements
        gridfile=trim(adjustl(trim(grid%path)//'constant/polyMesh/neighbour'))
        open(40,file=gridfile)

        do jj=1,19
            read(40,'(A)') inline
        enddo

        read(40,'(A)') inline
        read(inline,*) neigVal

        read(40,'(A)') inline     ! skip parenthesis
        allocate(grid%faceNeig(neigVal))
       
        do jj=1,neigVal
            read(40,'(A)') inline
            read(inline,*) intVal
            intVal = intVal + 1       ! from C++ to Fortran
            elemIter(intVal) = elemIter(intVal) +1        ! shift for the conn vector
            grid%elemFacesConn(grid%elemFacesConn(intVal) + elemIter(intVal)) = jj
            ! shift to current element, skip the face number
            grid%faceNeig(jj) = intVal
        end do

        close(40)

        deallocate(elemIter)

        do j=0,grid%nb_elems
            grid%connElem(j) = grid%elemFacesConn(j)
        end do

        do j=1,grid%nb_elems

            nF =   grid%elemFacesConn(grid%elemFacesConn(j))
            grid%connElem( grid%connElem(j) ) = nF

            do jface = 1,nF

                if (grid%elemFacesConn(grid%elemFacesConn(j) + jface) .le. neigVal) then
                    ! internal face
                    if (grid%faceOwner(grid%elemFacesConn(grid%elemFacesConn(j) + jface)) .eq. j) then
                          grid%connElem( grid%connElem(j) + jface ) =  grid%faceNeig(grid%elemFacesConn(grid%elemFacesConn(j) + jface))
                    else
                          grid%connElem( grid%connElem(j) + jface ) = grid%faceOwner(grid%elemFacesConn(grid%elemFacesConn(j) + jface))
                    endif
                else
                    ! boundary face
                    grid%connElem( grid%connElem(j) + jface ) = -1

                endif
            end do
        end do


    end subroutine READ_MESH_CONN_ELEMENTS

    subroutine READ_MESH_POINTS( grid )

        implicit none

        type(mesh) :: grid

        character*256 :: gridfile, inline
        integer*4 :: jj
        real*8    :: xx,yy,zz

        gridfile=trim(adjustl(trim(grid%path)//'constant/polyMesh/points'))

        !> Open mesh file
        open(40,file=gridfile)
        do jj=1,20
            read(40,'(A)') inline
        enddo
        
        ! a more consistent implementation has 3 vectors for x,y,z

        do jj=1,grid%nb_points
            read(40,'(A)') inline
            read(inline(2:len_trim(inline)-1),*) xx,yy,zz
            grid%points(jj,1) = xx
            grid%points(jj,2) = yy
            grid%points(jj,3) = zz
        end do


        close(40)

    end subroutine


    subroutine DUMP_MESH_VTK( grid, filename)

        implicit none
        type(mesh) :: grid
        character(len=*), intent(in) :: filename

        integer*4 :: VTK_file_unit
        integer*4 :: nelem, nsize
        integer*4 :: ii, jj, kk, nF, nV, nTot, faceIndex, facePointer

        write(*,*) "Exporting to vtk..."

        open(newunit=VTK_file_unit, action='WRITE', file=filename, &
                      form='FORMATTED', status='replace')

        write(VTK_file_unit,'(A)')'# vtk DataFile Version 3.0'
        write(VTK_file_unit,'(A)')'VTKFile'
        write(VTK_file_unit,'(A)')'ASCII'
        write(VTK_file_unit,'(A)')
        write(VTK_file_unit,'(A)')'DATASET UNSTRUCTURED_GRID'

        write(VTK_file_unit,'(A7,I12,A7)')'POINTS ',grid%nb_points,' double';
            POINT_LOOP: do jj=1,grid%nb_points
                write(VTK_file_unit,'(3E16.8,1X)') grid%points(jj,1:3)
            end do POINT_LOOP
        write(VTK_file_unit,'(A)')


        ! @ old VTK:
        nelem = grid%elemFacesConn(0)
        nsize = 0

        do jj = 1,nelem

            nF = grid%elemFacesConn(grid%elemFacesConn(jj))

            nTot = 0

            do kk=1,nF ! loop over the faces. Sum up all the vertices
                faceIndex   = grid%elemFacesConn(grid%elemFacesConn(jj) + kk )
                facePointer = grid%faceConnVector(faceIndex)
                nTot = nTot + grid%faceConnVector(facePointer)
            end do

            nTot = nTot + nF  +1          
            nsize = nTot + nsize
        end do


        write(VTK_file_unit,'(A7,I12,I12)')'CELLS ',nelem, nsize+nelem


        do jj = 1,nelem
            nF = grid%elemFacesConn(grid%elemFacesConn(jj))

            nTot = 0

            do kk=1,nF ! loop over the faces. Sum up all the vertices
                faceIndex   = grid%elemFacesConn(grid%elemFacesConn(jj) + kk )
                facePointer = grid%faceConnVector(faceIndex)
                nTot = nTot + grid%faceConnVector(facePointer)
            end do

            nTot = nTot + nF + 1            ! devo considerare anche se stesso

            write(VTK_file_unit,'(I0)',advance='no') nTot
            write(VTK_file_unit,'(A)',advance='no') ' '
            write(VTK_file_unit,'(I0)',advance='no') nF
            write(VTK_file_unit,'(A)',advance='no') ' '

            do kk=1,nF ! loop over the faces. Dump all the connectivity
                faceIndex   = grid%elemFacesConn(grid%elemFacesConn(jj) + kk )
                facePointer = grid%faceConnVector(faceIndex)
                nV = grid%faceConnVector(facePointer)
                write(VTK_file_unit,'(I0)', advance='no') nV
                write(VTK_file_unit,'(A)',advance='no') ' '

                do ii=1,nV
                    write(VTK_file_unit,'(I0)', advance='no') grid%faceConnVector(facePointer+ii)
                    write(VTK_file_unit,'(A)',advance='no') ' '
                end do

            end do

            ! newline
            write(VTK_file_unit,'(A)')
        end do

        write(VTK_file_unit,'(A11,I12)')'CELL_TYPES ',nelem

        do jj = 1,nelem
            write(VTK_file_unit,'(I12)') 42     ! assume evryone to be a polyhedron
        end do

        close(unit=VTK_file_unit)

    end subroutine   DUMP_MESH_VTK

end module
