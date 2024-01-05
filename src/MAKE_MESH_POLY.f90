subroutine MAKE_MESH_POLY( grid, gridPoly, nelem)

    ! create a vector of the original polyhedrons. NO CUT

    use fast_speed_polyhedron
    use qsort_c_module
    use MOD_OPENFOAM_GRID

    implicit none
    ! input
    integer*4 :: nelem
    type(mesh) :: grid
    type(polyhedron_fast) ,  dimension(nelem) :: gridPoly

    ! aux variables
    integer*4 :: dimConnFace, dimConnElem, dimPoints
    real*8, dimension(3,3) :: local_face_coords
    real*8, dimension(3)   :: normal_face, bar

    integer*4 :: ne1v, ne2v, ne3v, nVF
    integer*4 :: nUniques , nVdup, nF, faceIndex, facePointer
    integer*4 :: iFlu, kFace, kUnique, jVert, it, jDup,i           ! iterators

    integer*4, dimension(:), allocatable :: globalFaceMap, globalVertexList, globalVertexMap,&
                                            localVertexMap, globalVertexListSorted, tmp

    ! interface with the current mesh
    dimConnFace = grid%nb_faces_conn
    dimConnElem = grid%nb_elems_conn
    dimPoints   = grid%nb_points

    ! first read dimension
    ! assemble then the connectivity
    ! global map: map with the indeces of the fluid mesh
    ! local  map: map with the indeces of the fluid element

    do iFlu = 1,nelem                               ! loop over all fluid element in the mesh

        ! retrieve global maps

        nF = grid%elemFacesConn(grid%elemFacesConn(iFlu))       ! get the number of faces of the current fluid element
        
        allocate(globalFaceMap(nF))

        nVdup = 0                        ! number of vertices (with duplicates)

        ! get faces (global indexing)
        do kFace=1,nF                       ! loop over the faces. Sum up all the vertices
            faceIndex   = grid%elemFacesConn(grid%elemFacesConn(iFlu) + kFace )
            globalFaceMap(kFace) = faceIndex
            facePointer = grid%faceConnVector(faceIndex)
            nVdup = nVdup + grid%faceConnVector(facePointer)
        end do

        ! get vertices (with duplicates)
        allocate(globalVertexList(nVdup))
        allocate(globalVertexMap(nVdup))
        it = 0
        do kFace=1,nF
            faceIndex = globalFaceMap(kFace)
            facePointer = grid%faceConnVector(faceIndex)
            do jVert=1,grid%faceConnVector(facePointer)
                it = it+1
                globalVertexList(it) = grid%faceConnVector(facePointer+jVert) 
            end do
        end do

        ! get vertices (without duplicates)
        ! SORT

        ! copy the original sorted list
        allocate(globalVertexListSorted(nVdup))
        do jVert=1,nVdup
            globalVertexListSorted(jVert)=globalVertexList(jVert)
        end do

        call QsortC(globalVertexListSorted)

        ! COUNT UNIQUE
        allocate(localVertexMap(nVdup))     ! just to spare CPU time
        localVertexMap = -1
        nUniques = 1
        localVertexMap(nUniques)= globalVertexListSorted(1) 
        do jVert=2,nVdup
            if (localVertexMap(nUniques) .eq. globalVertexListSorted(jVert)) then
            else
                nUniques = nUniques + 1
                localVertexMap(nUniques) = globalVertexListSorted(jVert)
            endif
        end do

        ! localVertexMap is such that after nUniques it has only -1
        ! if you want it exact, pay 2 loops

        ! invert localVertexMap 
        do jDup=1,nVdup
            do kUnique=1,nUniques
                if (localVertexMap(kUnique) .eq. globalVertexList(jDup)) then
                    globalVertexMap(jDup) = kUnique
                end if
            end do
        end do
        ! globalVertexMap: index i grid%points to index j of the local map

        gridPoly(iFlu)%nV = nUniques
        gridPoly(iFlu)%nF = nF
        gridPoly(iFlu)%nc = nF + nF + nVdup

        ! devo popolare conn
        allocate(gridPoly(iFlu)%conn(0:gridPoly(iFlu)%nc))

        gridPoly(iFlu)%conn(0) = nF
        gridPoly(iFlu)%conn(1) = nF+1
        it = 0                                                ! iterator for vertex indexing
        kFace = 1
        facePointer = grid%faceConnVector(globalFaceMap(kFace))                                      ! size of the cur face
        gridPoly(iFlu)%conn(gridPoly(iFlu)%conn(kFace)) = grid%faceConnVector(facePointer)         ! loc/glo nV
        do jVert=1,grid%faceConnVector(facePointer)                                                  ! loop over the vertices
            it = it+1
            gridPoly(iFlu)%conn(gridPoly(iFlu)%conn(kFace)+jVert) = globalVertexMap(it)      ! get the mapped
        end do

        do kFace=2,nF
            facePointer = grid%faceConnVector(globalFaceMap(kFace-1))                                    ! size of the prev face
            gridPoly(iFlu)%conn(kFace) = grid%faceConnVector(facePointer) + &
                                          gridPoly(iFlu)%conn(kFace-1) + 1                        ! local connectivity
            facePointer = grid%faceConnVector(globalFaceMap(kFace))                                      ! size of the cur face
            gridPoly(iFlu)%conn(gridPoly(iFlu)%conn(kFace)) = grid%faceConnVector(facePointer)         ! loc/glo nV
            do jVert=1,grid%faceConnVector(facePointer)               ! loop over the vertices
                it = it+1
                gridPoly(iFlu)%conn(gridPoly(iFlu)%conn(kFace)+jVert) = globalVertexMap(it)      ! get the mapped
            end do
        end do


        allocate(gridPoly(iFlu)%vX(nUniques))
        allocate(gridPoly(iFlu)%vY(nUniques))
        allocate(gridPoly(iFlu)%vZ(nUniques))

        do jVert=1,nUniques                                                 ! convert from C to fortran
            gridPoly(iFlu)%vX(jVert) = grid%points(localVertexMap(jVert)+1,1)
            gridPoly(iFlu)%vY(jVert) = grid%points(localVertexMap(jVert)+1,2)
            gridPoly(iFlu)%vZ(jVert) = grid%points(localVertexMap(jVert)+1,3)
        end do

        ! check if the normal are outward
        bar(1) = sum(gridPoly(iFlu)%vX(:))/gridPoly(iFlu)%nV
        bar(2) = sum(gridPoly(iFlu)%vY(:))/gridPoly(iFlu)%nV
        bar(3) = sum(gridPoly(iFlu)%vZ(:))/gridPoly(iFlu)%nV

        do kFace=1,nF

            ne1v = gridPoly(iFlu)%conn(gridPoly(iFlu)%conn(kFace)+1)
            ne2v = gridPoly(iFlu)%conn(gridPoly(iFlu)%conn(kFace)+2)
            ne3v = gridPoly(iFlu)%conn(gridPoly(iFlu)%conn(kFace)+3)

            local_face_coords(1,1) =  gridPoly(iFlu)%vX(ne1v)         
            local_face_coords(2,1) =  gridPoly(iFlu)%vY(ne1v)         
            local_face_coords(3,1) =  gridPoly(iFlu)%vZ(ne1v)         

            local_face_coords(1,2) =  gridPoly(iFlu)%vX(ne2v)         
            local_face_coords(2,2) =  gridPoly(iFlu)%vY(ne2v)         
            local_face_coords(3,2) =  gridPoly(iFlu)%vZ(ne2v)         

            local_face_coords(1,3) =  gridPoly(iFlu)%vX(ne3v)         
            local_face_coords(2,3) =  gridPoly(iFlu)%vY(ne3v)         
            local_face_coords(3,3) =  gridPoly(iFlu)%vZ(ne3v)         

            normal_face(1) = (local_face_coords(2,2)-local_face_coords(2,1)) * &
                             (local_face_coords(3,3)-local_face_coords(3,1)) - &
                             (local_face_coords(3,2)-local_face_coords(3,1)) * &
                             (local_face_coords(2,3)-local_face_coords(2,1)) 

            normal_face(2) = (local_face_coords(3,2)-local_face_coords(3,1)) * &
                             (local_face_coords(1,3)-local_face_coords(1,1)) - &
                             (local_face_coords(3,3)-local_face_coords(3,1)) * &
                             (local_face_coords(1,2)-local_face_coords(1,1)) 

            normal_face(3) = (local_face_coords(1,2)-local_face_coords(1,1)) * &
                             (local_face_coords(2,3)-local_face_coords(2,1)) - &
                             (local_face_coords(1,3)-local_face_coords(1,1)) * &
                             (local_face_coords(2,2)-local_face_coords(2,1)) 

        if (((bar(1)-local_face_coords(1,1)) * normal_face(1) + &
             (bar(2)-local_face_coords(2,1)) * normal_face(2) + &
             (bar(3)-local_face_coords(3,1)) * normal_face(3)) .lt. 0.0d0) then

                ! do nothing

            else
                ! make the switcheroo
                nVF = gridPoly(iFlu)%conn(gridPoly(iFlu)%conn(kFace))
                allocate(tmp(nVF))
                ! copy
                do i=1,nVF
                    tmp(i) = gridPoly(iFlu)%conn(gridPoly(iFlu)%conn(kFace)+i)
                end do
                ! switch order
                do i=1,nVF
                    gridPoly(iFlu)%conn(gridPoly(iFlu)%conn(kFace)+i) = tmp(nVF-i+1)
                end do

                deallocate(tmp)
        end if

        end do


        ! deallocate useless stuff
        deallocate(globalFaceMap,globalVertexMap, globalVertexList, localVertexMap, globalVertexListSorted)

    end do

end subroutine MAKE_MESH_POLY
