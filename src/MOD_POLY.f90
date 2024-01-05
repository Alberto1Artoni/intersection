module fast_speed_polyhedron

    !! CLEAN after porting

    implicit none

    type polyhedron_fast
        integer*4 :: nV
        integer*4 :: nF
        integer*4 :: nc
        integer*4 , allocatable, dimension(:) :: conn
        real*8 , allocatable, dimension(:)    :: vX
        real*8 , allocatable, dimension(:)    :: vY
        real*8 , allocatable, dimension(:)    :: vZ
    end type polyhedron_fast

    type polygon_fast
        integer*4 :: nV
        integer*4 :: nc
        integer*4 , allocatable, dimension(:) :: conn
        real*8 , allocatable, dimension(:)    :: vX
        real*8 , allocatable, dimension(:)    :: vY
        real*8 , allocatable, dimension(:)    :: vZ
    end type polygon_fast

    type gaussData           
        real*8 :: volume
        real*8, dimension(:), allocatable :: faceArea
        real*8, dimension(:,:), allocatable :: faceNormal
        real*8, dimension(:), allocatable :: distance
        real*8 :: field
        real*8, dimension(:), allocatable :: faceField
    end type

    type map_intersection_fast
        integer*4 :: nIntersection
        integer*4 :: locCounter
        integer*4 , dimension(:), allocatable :: getFlu
        type(polyhedron_fast) , dimension(:), allocatable :: getPoly
        type(gaussData),dimension(:), allocatable :: getGauss
    end type

    integer*4, allocatable, dimension(:) :: mapLocToGloFluid
    type(polygon_fast) , allocatable, dimension(:) :: polyFluidBoundary
    type(map_intersection_fast), dimension(:), allocatable :: intersectionMapFast

    contains

    subroutine DEALLOCATE_INTERSECTION_MAP(ne)
        implicit none
        integer*4 :: ne, i, iFlu

        do i=1,ne
            do iFlu=1,intersectionMapFast(i)%nIntersection
                call POLY_DESTROY_FAST(intersectionMapFast(i)%getPoly(iFlu))
                deallocate(intersectionMapFast(i)%getFlu)

         !      deallocate(intersectionMapFast(i)%getGauss(iFlu)%faceField)
         !      deallocate(intersectionMapFast(i)%getGauss(iFlu)%distance)
         !      deallocate(intersectionMapFast(i)%getGauss(iFlu)%faceNormal)
         !      deallocate(intersectionMapFast(i)%getGauss(iFlu)%faceArea)

            end do
        end do
        
    end subroutine
    

    subroutine READ_OFF_FILE_DIME_FAST(filename, nvertices, nfaces, nedge)

        implicit none

        character*100 :: filename
        character*80  :: inline

        integer*4 :: nvertices, nfaces, nedge

        filename=trim(filename)
        open(40,file=filename)

        read(40,'(A)') inline       ! header
        ! okay, maybe it can be more general...

        read(40,'(A)') inline
        read(inline(1:len_trim(inline)),*) nvertices, nfaces, nedge

        close(40)

    end subroutine READ_OFF_FILE_DIME_FAST

    subroutine READ_OFF_FILE_FAST(filename, polyElem)

        implicit none

        character*100 :: filename
        character*80 :: inline
        real*8    :: xx,yy,zz
        integer*4 :: nV, nF, pos
        integer*4 :: i, j, it, iter
        integer*4, allocatable, dimension(:) :: getFace
        integer*4, allocatable, dimension(:) :: getFaceDim
        type(polyhedron_fast) :: polyElem


        open(40,file=filename)
        read(40,'(A)') inline       ! header
        read(40,'(A)') inline       ! dimensions, already read
        read(40,'(A)') inline       ! empty


        nV = polyElem%nV
        nF = polyElem%nF

        do i=1,nV
            read(40,'(A)') inline
            read(inline(1:len_trim(inline)),*) xx,yy,zz
            polyElem%vX(i) = xx
            polyElem%vY(i) = yy
            polyElem%vZ(i) = zz
        end do

        ! read face dimensions:

        allocate(getFace(nF*30))    ! take it big enough
                                    ! if you want to it precisely you should use the same
                                    ! construct as with the polyhedron
                                    ! however the allocation and dellocation might cost too much
                                    ! this should be faster altough less robust
        allocate(getFaceDim(nF))

        it = 0
        do i=1,nF
            read(40,'(A)') inline
            inline = trim(inline)       ! remove trailing and leading space
            pos    = index(inline,' ')
            read(inline(1:pos),*) getFaceDim(i)

            do j=1,getFaceDim(i)
                inline = adjustl(trim(inline(pos:len_trim(inline))))
                pos    = index(inline,' ')
                it = it + 1
                read(inline(1:pos),*) getFace(it)
            end do
        end do

        allocate(polyElem%conn(0:nF+nF+it))
        polyElem%nc = nF+nF+it
        polyElem%conn(0) = nF
        polyElem%conn(1) = nF + 1
!       write(*,*) polyElem%conn(1)

        do i= 2,nF
            polyElem%conn(i) = getFaceDim(i-1) + polyElem%conn(i-1) + 1
                                     ! dim faccia     ! puntatore        ! +1 per il numero di facce 
        end do

        iter = 0

        do i=1,nF
            polyElem%conn(polyElem%conn(i)) = getFaceDim(i)
            do j=1,getFaceDim(i)
                iter = iter+1
                polyElem%conn(polyElem%conn(i)+j) = getFace(iter) + 1   ! c++ to fortran index conversion
            end do
        end do

        deallocate(getFace)
        deallocate(getFaceDim)
        close(40)

    end subroutine

    subroutine WRITE_OFF_POLY_FAST(poly, outString)
    ! this methods writes the polyhedron as a string
        implicit none
        type(polyhedron_fast) :: poly
        character(5000) :: outString        ! take a sufficiently big string
        character(100) :: tmp
        integer*4 :: i, j, dimFace

        ! the strategy is to build outString
        ! outString is a line containing the OFF polyhedron format

        write(outString,'(A)') "OFF"//NEW_LINE('a')
        write(outString,'(A,I0,A,I0,A,I0)') trim(adjustl(outString)),poly%nV," ", poly%nF," ", 0
        write(outString, *) trim(outString)//NEW_LINE('a')//NEW_LINE('a')

        do i=1,poly%nV 
            write(tmp, '(F15.8)') poly%vX(i)
            tmp = trim(adjustl(tmp))
            write(outString,*) trim(outString)//tmp

            write(tmp, '(F15.8)') poly%vY(i)
            tmp = trim(adjustl(tmp))
            write(outString,*) trim(outString)//" "//tmp

            write(tmp, '(F15.8)') poly%vZ(i)
            tmp = trim(adjustl(tmp))
            write(outString,*) trim(outString)//" "//tmp

            write(outString,*) trim(outString)//NEW_LINE('a')
        end do


        do i=1,poly%nF

            dimFace = poly%conn(poly%conn(i))
            write(tmp, '(I0)') dimFace          ! write face dimension

            do j=(poly%conn(i)+1) , (poly%conn(i) + dimFace)
                write(tmp, '(A,I0)') trim(tmp)//" ",poly%conn(j) -1
            end do

            write(outString, *) trim(outString)//trim(tmp)//NEW_LINE('a')

        end do

    end subroutine WRITE_OFF_POLY_FAST

    subroutine POLY_DESTROY_FAST(poly)
        implicit none
        type(polyhedron_fast) :: poly
        integer :: i

        deallocate(poly%vX)                              ! deallocate vertices
        deallocate(poly%vY)                              ! deallocate vertices
        deallocate(poly%vZ)                              ! deallocate vertices

        deallocate(poly%conn)     ! deallocate each face conn

    end subroutine

    subroutine EXPORT_VTK_POLYDATA_DATA_FAST(filename,poly,val, val2, val3)
        ! overload of EXPORT_VTK_POLYDATA in order to rapresent the fluid to acustic map
        ! export only ONE polyhedron

        implicit none

        type(polyhedron_fast) :: poly
        character*100 :: filename, tmp
        integer*4 :: i,j, ntot, val, val2
        real*8, optional :: val3
        
        open(50, file=filename)
        write(50,'(A)') "# vtk DataFile Version 2.0"
        write(50,'(A)') "VTK from Matlab"
        write(50,'(A)') "ASCII"
        write(50,'(A)') "DATASET POLYDATA"
        write(50,'(A,I0,A)') "POINTS ", poly%nV, " float"

        do i=1,poly%nV
            write(50,'(F22.14,1X,F22.14,1X,F22.14,1X)') poly%vX(i),poly%vY(i),poly%vZ(i)
        end do

        ntot = 0
        do i=1,poly%nF
            ntot = ntot + poly%conn(poly%conn(i))
        end do
        ntot = ntot + poly%nF

        write(50,'(A,I0,A,I0)') "POLYGONS ", poly%nF," ", ntot
        do i=1,poly%nF
            write(tmp, '(I0,A)') poly%conn(poly%conn(i))," "
            do j=1,poly%conn(poly%conn(i))
                write(tmp, '(A,I0)') trim(tmp)//" ", poly%conn(poly%conn(i)+j)-1
                ! C like convention
            end do
            write(50,'(A)') trim(adjustl(tmp))
        end do

        write(50,'(A,I0)') "CELL_DATA ", poly%nF
        write(50,'(A)') "SCALARS ACU_ID double 1"
        write(50,'(A)') "LOOKUP_TABLE default"
        do i=1,poly%nF
            write(50,'(I0)') val
        end do
        write(50,'(A)') "SCALARS FLU_ID double 1"
        write(50,'(A)') "LOOKUP_TABLE default"
        do i=1,poly%nF
            write(50,'(I0)') val2
        end do

        if (present(val3)) then
            write(50,'(A)') "SCALARS VAL double 1"
            write(50,'(A)') "LOOKUP_TABLE default"
            do i=1,poly%nF
                write(50,'(F12.5)') val3
            end do
        endif

        close(50)

    end subroutine EXPORT_VTK_POLYDATA_DATA_FAST

    subroutine GET_POLY_BBOX(xMin,yMin,zMin, xMax, yMax, zMax, poly)

        implicit none
        real*8 :: xMin,yMin,zMin, xMax, yMax, zMax
        type(polyhedron_fast) ::  poly

        ! sufficit to do it on the whole points stuff
        xMin = MINVAL(poly%vX(:))
        yMin = MINVAL(poly%vY(:))
        zMin = MINVAL(poly%vZ(:))
        xMax = MAXVAL(poly%vX(:))
        yMax = MAXVAL(poly%vY(:))
        zMax = MAXVAL(poly%vZ(:))

    end subroutine

end module

module intersectionCGAL_new
  implicit none

  private
  public :: intersectPolyhedra, INTEGRAL_MONOMIALS_POLYHEDRON_3D_FAST,  &
            GET_POLY_VOLUME, GET_POLY_FACES_AREA

  interface

    subroutine free(c_pointer) bind (C, name = "free")
    use, intrinsic :: iso_c_binding
        type(c_ptr), value              :: c_pointer
    end subroutine

    subroutine     nefInterface   (     &
                           nVAcu,  nFAcu,  ncAcu,               &
                           vxAcu,  vyAcu,  vzAcu,  connAcu,     &
                           nVFlu,  nFFlu,  ncFlu,               &
                           vxFlu,  vyFlu,  vzFlu,  connFlu,     &
                           nV,  nF,  nc,            &
                           vx,  vy,  vz,  conn) bind(C, name="nefInterface")
      use, intrinsic :: iso_c_binding
      integer*4 :: nVAcu, nFAcu, ncAcu
      integer*4 :: nVFlu, nFFlu, ncFlu
      real*8, dimension(nVAcu) :: vxAcu, vyAcu, vzAcu
      real*8, dimension(nVFlu) :: vxFlu, vyFlu, vzFlu
      integer*4, dimension(ncAcu) :: connAcu
      integer*4, dimension(ncFlu) :: connFlu
      integer*4 :: nV, nF, nc
      type(c_ptr) :: vx, vy, vz
      type(c_ptr) :: conn

    end subroutine nefInterface

  end interface

contains

    subroutine CONVERT_C_F_POINTER(nV, nF, nc, vx, vy, vz, conn, poly)
        use fast_speed_polyhedron
        use, intrinsic :: iso_c_binding
        implicit none
        type(polyhedron_fast) :: poly
        integer*4 :: nV, nF, nc, i, j
        integer(kind=c_int),pointer :: fp(:)
        real(kind=c_double),pointer :: fvx(:), fvy(:), fvz(:)
        type(c_ptr) :: vx, vy, vz
        type(c_ptr) :: conn

    !   write(*,*) nV, nF, nc                       
        
        call c_f_pointer(conn, fp, [nc+1])
        call c_f_pointer(vx, fvx, [nV])
        call c_f_pointer(vy, fvy, [nV])
        call c_f_pointer(vz, fvz, [nV])
        
        poly%nV = nV
        poly%nF = nF
        poly%nc = nc

        allocate(poly%conn(0:nc))

        poly%conn(0) = nF

!       write(*,*) fp(:)
!       write(*,*) fp(1)
!       write(*,*) fp(2:nF+1)
!       write(*,*) fp(nF+2:nc+1)

        do i=2,(nF+1)
!           write(*,*) i-1
            poly%conn(i-1) = fp(i)
            poly%conn(poly%conn(i-1)) = fp(fp(i)+1)

            do j= fp(i)+1, (fp(i)+fp(fp(i)+1))
!                write(*,*) j
                 poly%conn(j) = fp(j+1) + 1
            end do


        end do
   !    write(*,*) fp(:)
   !    write(*,*) poly%conn(:)

        allocate(poly%vX(nV), poly%vY(nV), poly%vZ(nV))
        do i=1,nV
            poly%vX(i) = fvx(i)
            poly%vY(i) = fvy(i)
            poly%vZ(i) = fvz(i)
        end do

    end subroutine CONVERT_C_F_POINTER

    subroutine intersectPolyhedra(polyAcu,polyFlu, poly, isIntersecting)

        use fast_speed_polyhedron
        use, intrinsic :: iso_c_binding
        implicit none
        type(polyhedron_fast) :: polyAcu, polyFlu, poly
        integer*4 :: nV, nF, nc, i, j
        type(c_ptr) :: vx, vy, vz
        type(c_ptr) :: conn
        integer(kind=c_int),pointer :: fp(:)
        real(kind=c_double),pointer :: fvx(:), fvy(:), fvz(:)
        logical :: isIntersecting

        
        call nefInterface(polyAcu%nV, polyAcu%nF, polyAcu%nc, &
                          polyAcu%vx, polyAcu%vy, polyAcu%vz, polyAcu%conn, &
                          polyFlu%nV, polyFlu%nF, polyFlu%nc, &
                          polyFlu%vx, polyFlu%vy, polyFlu%vz, polyFlu%conn, & 
                          nV, nF, nc, vx, vy, vz, conn)

        if ( nV .eq. -2) then
            write(*,*) "triangulate"
        endif                   

        if (nV .gt. 0) then
            isIntersecting = .true.

            call CONVERT_C_F_POINTER(nV, nF, nc, vx, vy, vz, conn, poly)

            call free(vx)
            call free(vy)
            call free(vz)
            call free(conn)

        else
            poly%nV = 0
            poly%nF = 0
            poly%nc = 0
            isIntersecting = .false.
        endif


    end subroutine intersectPolyhedra

    subroutine INTEGRAL_MONOMIALS_POLYHEDRON_3D_FAST(poly, integral_monomials, p)
        ! QUAD FREE RULE: Pennesi and Antonietti
        use fast_speed_polyhedron
        implicit none
      
        type(polyhedron_fast), intent(in) :: poly
        integer, intent(in) :: p !< Maximal order to be integrated
        real*8, dimension(0:p,0:p,0:p), intent(out) :: integral_monomials
        character(5000) :: outString
        !< Integral of monomials over the current (mapped) element

    ! Local variables

        integer :: i,j,k, no_edges,no_face_nodes,ierr, face_no, l_xi, l_eta, l_zeta,iface
        real*8 :: x1,x2,y1,y2,length, b_i, b_ij, b_ijk
        real*8, dimension(3) :: point, X0i, X0ij, v1, v2, v3, normal_face, normal_edge, normal_pt,bar
        real*8, dimension(0:p,0:p,0:p) :: face_integral_monomials, edge_integral_monomials
        real*8, allocatable, dimension(:,:) :: local_face_coords

    

        integral_monomials= 0.0d0
    
        bar(1) = sum(poly%vX(:))/poly%nV
        bar(2) = sum(poly%vY(:))/poly%nV
        bar(3) = sum(poly%vZ(:))/poly%nV

!       write(*,*) "Baricentro:"
!       write(*,*) bar
!       write(*,*)

        do iface = 1, poly%nF       ! loop over the number of faces
            
            face_integral_monomials = 0.0d0;
                            
            no_face_nodes = poly%conn(poly%conn(iface))
            
            no_edges = no_face_nodes ! now this is true because a face is a polygon
            
            ! Map the polyhedron to the reference element
            ! if d=3 then no_face_node is different from 2

            allocate(local_face_coords(3,no_face_nodes+1))
            ! local_face_coords:
            ! 1: no_face_nodes => nodi locali
            ! no_face_nodes+1  => poly%v(:,1)
            ! questo perché poi ci accede come se fossero edge, cioé
            
            do k = 1,no_face_nodes
                 local_face_coords(1,k) = poly%vX( poly%conn(poly%conn(iface)+k)  )
                 local_face_coords(2,k) = poly%vY( poly%conn(poly%conn(iface)+k)  )
                 local_face_coords(3,k) = poly%vZ( poly%conn(poly%conn(iface)+k)  )
            end do

            ! Choose X0i as the center of mass of the face
            X0i(1) = sum(local_face_coords(1,1:no_face_nodes))/no_face_nodes
            X0i(2) = sum(local_face_coords(2,1:no_face_nodes))/no_face_nodes
            X0i(3) = sum(local_face_coords(3,1:no_face_nodes))/no_face_nodes
            
           ! to check 
            normal_face(1) = (local_face_coords(2,2)-local_face_coords(2,1)) * &
                             (X0i(3)-local_face_coords(3,1)) - &
                             (local_face_coords(3,2)-local_face_coords(3,1)) * &
                             (X0i(2)-local_face_coords(2,1)) 

            normal_face(2) = (local_face_coords(3,2)-local_face_coords(3,1)) * &
                             (X0i(1)-local_face_coords(1,1)) - &
                             (X0i(3)-local_face_coords(3,1)) * &
                             (local_face_coords(1,2)-local_face_coords(1,1)) 

            normal_face(3) = (local_face_coords(1,2)-local_face_coords(1,1)) * &
                             (X0i(2)-local_face_coords(2,1)) - &
                             (X0i(1)-local_face_coords(1,1)) * &
                             (local_face_coords(2,2)-local_face_coords(2,1)) 

            if (((bar(1)-local_face_coords(1,1)) * normal_face(1) + &
                 (bar(2)-local_face_coords(2,1)) * normal_face(2) + &
                 (bar(3)-local_face_coords(3,1)) * normal_face(3)) .gt. 0.0d0) then
                    normal_face = -normal_face
            end if

            normal_face = normal_face / dsqrt(normal_face(1)**2 + normal_face(2)**2 + normal_face(3)**2)

            b_i = dot_product(normal_face,local_face_coords(:,1))
            
            ! Now loop over edges of this face
            local_face_coords(:,no_face_nodes+1) = local_face_coords(:,1)
            
            do j = 1,no_face_nodes
              
                  edge_integral_monomials = 0.0d0
                  
                  v1 = local_face_coords(:,j); v2 = local_face_coords(:,j+1)
                  
                  ! Choose X0ij, nij1, nij2, dij1, dij2
                  X0ij = v1 !and thanks to that, we don't need nij1 nor dij1
                
                  normal_pt = (v2-v1)


                  b_ijk = NORM2(normal_pt)
                  normal_pt = normal_pt/b_ijk

              
              ! Find nij, that is the normal to this edge lying on the hyperplane of the current face
                  normal_edge = dot_product(X0i-v1,normal_pt)*normal_pt - (X0i-v1)


                  normal_edge = normal_edge/NORM2(normal_edge)

                  ! Algebraic distance between X0i and this edge
                  b_ij = dot_product(v1-X0i,normal_edge)

                  ! Now start the homogeneous function integration
                  edge_integral_monomials(0,0,0) = b_ijk/(real(3-2,8))
                  face_integral_monomials(0,0,0) = face_integral_monomials(0,0,0) &
                    + b_ij*edge_integral_monomials(0,0,0)

                  do l_xi = 1,p

                    edge_integral_monomials(l_xi,0,0) = ( b_ijk*v2(1)**l_xi &
                      + X0ij(1)*real(l_xi,8)*edge_integral_monomials(l_xi-1,0,0) &
                      )/(real(3+l_xi-2,8))
                    
                    edge_integral_monomials(0,l_xi,0) = ( b_ijk*v2(2)**l_xi &
                      + X0ij(2)*real(l_xi,8)*edge_integral_monomials(0,l_xi-1,0) &
                      )/(real(3+l_xi-2,8))
                    
                    edge_integral_monomials(0,0,l_xi) = ( b_ijk*v2(3)**l_xi &
                      + X0ij(3)*real(l_xi,8)*edge_integral_monomials(0,0,l_xi-1) &
                      )/(real(3+l_xi-2,8))


                    face_integral_monomials(l_xi,0,0) = face_integral_monomials(l_xi,0,0) &
                      + b_ij*edge_integral_monomials(l_xi,0,0)
                    
                    face_integral_monomials(0,l_xi,0) = face_integral_monomials(0,l_xi,0) &
                      + b_ij*edge_integral_monomials(0,l_xi,0)
                    
                    face_integral_monomials(0,0,l_xi) = face_integral_monomials(0,0,l_xi) &
                      + b_ij*edge_integral_monomials(0,0,l_xi)

                  end do

                  do l_xi = 1,p
                    do l_eta = 1,p-l_xi

                      edge_integral_monomials(l_xi,l_eta,0) = ( b_ijk*(v2(1)**l_xi)*(v2(2)**l_eta) &
                        + X0ij(1)*real(l_xi,8)*edge_integral_monomials(l_xi-1,l_eta,0) &
                        + X0ij(2)*real(l_eta,8)*edge_integral_monomials(l_xi,l_eta-1,0) &
                        )/(real(3+l_xi+l_eta-2,8))

                      edge_integral_monomials(l_xi,0,l_eta) = ( b_ijk*(v2(1)**l_xi)*(v2(3)**l_eta) &
                        + X0ij(1)*real(l_xi,8)*edge_integral_monomials(l_xi-1,0,l_eta) &
                        + X0ij(3)*real(l_eta,8)*edge_integral_monomials(l_xi,0,l_eta-1) &
                        )/(real(3+l_xi+l_eta-2,8))

                      edge_integral_monomials(0,l_xi,l_eta) = ( b_ijk*(v2(2)**l_xi)*(v2(3)**l_eta) &
                        + X0ij(2)*real(l_xi,8)*edge_integral_monomials(0,l_xi-1,l_eta) &
                        + X0ij(3)*real(l_eta,8)*edge_integral_monomials(0,l_xi,l_eta-1) &
                        )/(real(3+l_xi+l_eta-2,8))


                      face_integral_monomials(l_xi,l_eta,0) = face_integral_monomials(l_xi,l_eta,0) &
                        + b_ij*edge_integral_monomials(l_xi,l_eta,0)

                      face_integral_monomials(l_xi,0,l_eta) = face_integral_monomials(l_xi,0,l_eta) &
                        + b_ij*edge_integral_monomials(l_xi,0,l_eta)
                      
                      face_integral_monomials(0,l_xi,l_eta) = face_integral_monomials(0,l_xi,l_eta) &
                        + b_ij*edge_integral_monomials(0,l_xi,l_eta)

                    end do
                  end do


                  do l_xi = 1,p
                    do l_eta = 1,p-l_xi
                      do l_zeta = 1,p-l_xi-l_eta

                        edge_integral_monomials(l_xi,l_eta,l_zeta) = ( &
                          b_ijk*(v2(1)**l_xi)*(v2(2)**l_eta)*(v2(3)**l_zeta) &
                          + X0ij(1)*real(l_xi,8)*edge_integral_monomials(l_xi-1,l_eta,l_zeta) &
                          + X0ij(2)*real(l_eta,8)*edge_integral_monomials(l_xi,l_eta-1,l_zeta) &
                          + X0ij(3)*real(l_zeta,8)*edge_integral_monomials(l_xi,l_eta,l_zeta-1) &
                          )/(real(3+l_xi+l_eta+l_zeta-2,8)) 

                        face_integral_monomials(l_xi,l_eta,l_zeta) = face_integral_monomials(l_xi,l_eta,l_zeta) &
                          + b_ij*edge_integral_monomials(l_xi,l_eta,l_zeta)

                      end do
                    end do
                  end do

            end do ! loop over edges

            face_integral_monomials(0,0,0) = face_integral_monomials(0,0,0)/(real(3-1,8))
            integral_monomials(0,0,0) = integral_monomials(0,0,0) + b_i*face_integral_monomials(0,0,0)/real(3,8)

            do l_xi = 1,p

                  face_integral_monomials(l_xi,0,0) = (face_integral_monomials(l_xi,0,0) &
                    + X0i(1)*real(l_xi,8)*face_integral_monomials(l_xi-1,0,0) )/(real(3+l_xi-1,8))

                  face_integral_monomials(0,l_xi,0) = (face_integral_monomials(0,l_xi,0) &
                    + X0i(2)*real(l_xi,8)*face_integral_monomials(0,l_xi-1,0) )/(real(3+l_xi-1,8))

                  face_integral_monomials(0,0,l_xi) = (face_integral_monomials(0,0,l_xi) &
                    + X0i(3)*real(l_xi,8)*face_integral_monomials(0,0,l_xi-1) )/(real(3+l_xi-1,8))


                  integral_monomials(l_xi,0,0) = integral_monomials(l_xi,0,0) &
                    + b_i*face_integral_monomials(l_xi,0,0)/(real(3+l_xi,8))
                  integral_monomials(0,l_xi,0) = integral_monomials(0,l_xi,0) &
                    + b_i*face_integral_monomials(0,l_xi,0)/(real(3+l_xi,8))
                  integral_monomials(0,0,l_xi) = integral_monomials(0,0,l_xi) &
                    + b_i*face_integral_monomials(0,0,l_xi)/(real(3+l_xi,8))

            end do


            do l_xi = 1,p
                  do l_eta = 1,p-l_xi

                    face_integral_monomials(l_xi,l_eta,0) = (face_integral_monomials(l_xi,l_eta,0) &
                      + X0i(1)*real(l_xi,8)*face_integral_monomials(l_xi-1,l_eta,0) &
                      + X0i(2)*real(l_eta,8)*face_integral_monomials(l_xi,l_eta-1,0) &
                      )/(real(3+l_xi+l_eta-1,8))

                    face_integral_monomials(l_xi,0,l_eta) = (face_integral_monomials(l_xi,0,l_eta) &
                      + X0i(1)*real(l_xi,8)*face_integral_monomials(l_xi-1,0,l_eta) &
                      + X0i(3)*real(l_eta,8)*face_integral_monomials(l_xi,0,l_eta-1) &
                      )/(real(3+l_xi+l_eta-1,8))

                    face_integral_monomials(0,l_xi,l_eta) = (face_integral_monomials(0,l_xi,l_eta) &
                      + X0i(2)*real(l_xi,8)*face_integral_monomials(0,l_xi-1,l_eta) &
                      + X0i(3)*real(l_eta,8)*face_integral_monomials(0,l_xi,l_eta-1) &
                      )/(real(3+l_xi+l_eta-1,8))


                    integral_monomials(l_xi,l_eta,0) = integral_monomials(l_xi,l_eta,0) &
                      + b_i*face_integral_monomials(l_xi,l_eta,0)/(real(3+l_xi+l_eta,8))
                    
                    integral_monomials(l_xi,0,l_eta) = integral_monomials(l_xi,0,l_eta) &
                      + b_i*face_integral_monomials(l_xi,0,l_eta)/(real(3+l_xi+l_eta,8))
                    
                    integral_monomials(0,l_xi,l_eta) = integral_monomials(0,l_xi,l_eta) &
                      + b_i*face_integral_monomials(0,l_xi,l_eta)/(real(3+l_xi+l_eta,8))

                  end do
            end do


            do l_xi = 1,p
              do l_eta = 1,p-l_xi
                do l_zeta = 1,p-l_xi-l_eta

                  face_integral_monomials(l_xi,l_eta,l_zeta) = ( &
                    face_integral_monomials(l_xi,l_eta,l_zeta) &
                    + X0i(1)*real(l_xi,8)*face_integral_monomials(l_xi-1,l_eta,l_zeta) &
                    + X0i(2)*real(l_eta,8)*face_integral_monomials(l_xi,l_eta-1,l_zeta) &
                    + X0i(3)*real(l_zeta,8)*face_integral_monomials(l_xi,l_eta,l_zeta-1) &
                    )/(real(3+l_xi+l_eta+l_zeta-1,8))

                  integral_monomials(l_xi,l_eta,l_zeta) = integral_monomials(l_xi,l_eta,l_zeta) &
                    + b_i*face_integral_monomials(l_xi,l_eta,l_zeta)/(real(3+l_xi+l_eta+l_zeta,8))

                end do
              end do
            end do

        
            deallocate(local_face_coords)


        end do ! loop over faces

    if ( isnan(integral_monomials(0,0,0)) ) then

        write(*,*) "Dumping off:"

        call WRITE_OFF_POLY_FAST(poly, outString)
        write(*,*) TRIM(outString)

        write(*,*) "Warning! During the integration a NaN was met. Monomials are set to zero."
        integral_monomials = 0.0d0
    endif

      
   end subroutine INTEGRAL_MONOMIALS_POLYHEDRON_3D_FAST

   subroutine GET_POLY_VOLUME(poly, vol)
        use fast_speed_polyhedron

        implicit none

        type(polyhedron_fast), intent(in) :: poly
        real*8, intent(out) :: vol

        real*8, dimension(:,:,:), allocatable :: integral_monomials

        ! compute poly volume
 
        allocate(integral_monomials(0:1, 0:1, 0:1))      ! forse 3*p?
        integral_monomials = 0.0d0
 
        call INTEGRAL_MONOMIALS_POLYHEDRON_3D_FAST(  poly, integral_monomials, 1)
 
        vol = integral_monomials(0,0,0)

   end subroutine GET_POLY_VOLUME

   subroutine TRIANGLE_AREA(poly, ne1,ne2, barFace, triangleArea)
        use fast_speed_polyhedron

        implicit none

        type(polyhedron_fast), intent(in) :: poly
        integer*4 :: ne1, ne2
        real*8, dimension(3) :: AB, AC, barFace
        real*8 :: triangleArea

        AB(1) = poly%vX(ne2) - poly%vX(ne1)
        AB(2) = poly%vY(ne2) - poly%vY(ne1)
        AB(3) = poly%vZ(ne2) - poly%vZ(ne1)

        AC(1) = barFace(1) - poly%vX(ne1)
        AC(2) = barFace(2) - poly%vY(ne1)
        AC(3) = barFace(3) - poly%vZ(ne1)

        triangleArea = 0.5 * sqrt(( AB(2)*AC(3) - AB(3)*AC(2) )**2 +    &
                                  ( AB(3)*AC(1) - AB(1)*AC(3) )**2 +    &
                                  ( AB(1)*AC(2) - AB(2)*AC(1) )**2 )

   end subroutine

   subroutine GET_POLY_FACES_AREA( poly, nF, faceArea)
        use fast_speed_polyhedron

        implicit none

        type(polyhedron_fast), intent(in) :: poly
        integer*4 , intent(in) :: nF
        real*8, dimension(nF), intent(out) :: faceArea
        integer*4 :: it, nFaceVerts, ne1, ne2, ne3, j
        real*8 :: triangleArea
        real*8, dimension(3) :: barFace
        
        faceArea = 0.0d0
        do j=1, nF

            nFaceVerts = poly%conn(poly%conn(j))
            barFace = 0.0d0
            do it=1,nFaceVerts
                barFace(1) = barFace(1) + poly%vX(poly%conn(poly%conn(j) + it))
                barFace(2) = barFace(2) + poly%vY(poly%conn(poly%conn(j) + it))
                barFace(3) = barFace(3) + poly%vZ(poly%conn(poly%conn(j) + it))
            end do

            barFace = barFace / nFaceVerts

            do it=1,nFaceVerts-1
                ne1 = poly%conn(poly%conn(j) + it  )
                ne2 = poly%conn(poly%conn(j) + it+1)
                ! compute area of the triangle
                call TRIANGLE_AREA(poly, ne1,ne2, barFace,triangleArea)
                faceArea(j) = faceArea(j) + triangleArea
            end do

            ne1 = poly%conn(poly%conn(j) + nFaceVerts)
            ne2 = poly%conn(poly%conn(j) + 1)

            call TRIANGLE_AREA(poly, ne1,ne2, barFace, triangleArea)
            faceArea(j) = faceArea(j) + triangleArea

        end do

   end subroutine 

   
end module  intersectionCGAL_new
